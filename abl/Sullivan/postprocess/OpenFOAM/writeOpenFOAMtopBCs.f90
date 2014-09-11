subroutine nearestPoints(x, y, ixl, ixu, iyl, iyu, wxlyl, wxuyl, wxlyu, wxuyu)
  double precision, intent(in) :: x, y
  integer, intent(inout) :: ixl, ixu, iyl, iyu
  double precision, intent(inout) :: wxlyl, wxuyl, wxlyu, wxuyu
  double precision :: wxl, wxu, wyl, wyu
  integer :: nx = 768, ny=768
  double precision :: xl=5120.0, yl=5120.0
  double precision :: dx, dy 
  
  dx = xl/dble(nx)
  dy = yl/dble(ny)
  ixl = mod(  (floor(x/dx)+nx), nx ) + 1
  ixu = mod(  (floor(x/dx)+nx+1), nx ) + 1
!  write(*,*) "Point ", x, " is located in between ", ixl, " and ", ixu
  wxu = dmod(x+xl,dx)/dx
  wxl = 1.0 - wxu

  iyl = mod(  (floor(y/dy)+ny), ny ) + 1
  iyu = mod(  (floor(y/dy)+ny+1), ny ) + 1
!  write(*,*) "Point ", y, " is located in between ", iyl, " and ", iyu
  wyu = dmod(y+yl,dy)/dy
  wyl = 1.0 - wyu

!  write(*,*) 'wxl = ', wxl, 'wxu = ', wxu, 'wyl = ', wyl, 'wyu = ', wyu
!  write(*,*) 'ixl = ', ixl, 'iyl = ', iyl, 'ixu = ', ixu, 'iyu = ', iyu 
  
  wxlyl = wxl*wyl
  wxuyl = wxu*wyl
  wxlyu = wxl*wyu
  wxuyu = wxu*wyu
  return
end subroutine nearestPoints

program writeOpenFOAMBCs
  !This is a program to read the wtbox files output by the LES code and write them in the format required by the timeVaryingMappedFvPatchField in OpenFOAM
  implicit none

  integer :: i,j,k,it,fileCounter,iz !Counter variables
  
  real :: xl, yl         ! x, y domain sizes
  real :: dx, dy, dz         ! x, y grid lengths
  integer :: nnx, nny, nnz      ! x, y grid dimensions
  integer :: nvar=5 !The number of variables stored in a plane at every time step
  integer :: nxy !nnx*nny
  integer :: wtExtentY, wtExtentZ !The spanwise and vertical extent of the wind turbine in number of grid points.
  integer :: wtx, wty !Number of wind turbines in each direction
  integer :: nt !Number of time steps
  real :: dt  !Time step
  real(kind=4), dimension(:,:,:), allocatable :: pA_xy !Plane Arrays 
  real(kind=4), dimension(:,:), allocatable :: u,v,t,e !Plane Arrays 
  real(kind=4), dimension(:,:,:), allocatable :: w !Plane Arrays 
  real, allocatable, dimension(:) :: umeanProfile, vmeanProfile, wmeanProfile, tmeanProfile !Profile of mean velocities and temperature from postprocessing of his.mp files
  integer :: zLevel !The z level at which this analysis is to be done.
  integer :: sizeOfReal=1 !The size of real in words on this system
  integer :: fileXY, fileDT  !The index of the files used to open the "xy.data" and the "dt" files
  integer :: fileUmeanProfile, fileTmeanProfile !The index of the file used to open the "uxym" file
  real(kind=4) :: timeStart, tLESnew, tLESold
  real(kind=4) :: dtLES
  integer :: itLES
  character :: dtFileReadLine
  integer :: meanOrTurb ! = 1 for writing turbulent velocity fields and = 0 for mean fields

  !Variables to interpolate to yawed co-ordinate system
  real(kind=4) :: yawAngle
  real(kind=4) :: xLocCur, yLocCur, zLocCur
  integer :: ixl, ixu, iyl, iyu !The integer locations of the 4 points used to interpolate the velocity from the LES of ABL field to the reqd. point in the FAST field
  double precision :: wxlyl, wxuyl, wxlyu, wxuyu !The weights of the 4 points used to interpolate the velocity from the LES of ABL field to the reqd. point in the FAST field
  double precision :: uRot, vRot, wRot, tRot, eRot !The interpolated value of velocity at the reqd. point in the FAST field


  character*25 :: command

  meanOrTurb = 1 !Write turbulent velocity fields
  if(meanOrTurb == 0) then
     write(*,*) 'Writing mean profiles!!!  '
     write(*,*) 'If you want turbulent flow field, set meanOrTurb = 1'
  else
     write(*,*) 'Writing turbulent velocity fields!!!  '
     write(*,*) 'If you want just the mean flow field, set meanOrTurb = 0'
  end if

  yawAngle = 21.5*3.14159265359/180.0
  xl = 5120.0
  yl = 5120.0
  nnx = 768
  nny = 768
  nnz = 50
  dx = 5120./768.0
  dy = 5120./768.0
  dz = 2048./256.0
  sizeOfReal = 1
  fileXY = 95
  fileDT = 91
  fileUMeanProfile = 83
  fileTMeanProfile = 84

  wtExtentZ = 35
  wtExtentY = 25
  wtx = 7
  wty = 7
  
  nxy = nnx*nny
  nt =  1500
  
  allocate(pA_xy(nvar,nnx,nny))
  allocate(u(nnx,nny))
  allocate(v(nnx,nny))
  allocate(w(nnx,nny,0:1))
  allocate(t(nnx,nny))
  allocate(e(nnx,nny))
  allocate(umeanProfile(nnz))
  allocate(vmeanProfile(nnz))
  allocate(wmeanProfile(nnz))
  allocate(tmeanProfile(nnz))
  w(:,:,0) = 0.0
  !Read in umean and uvar Profiles
  open(fileUMeanProfile,file="uxym")
  read(fileUMeanProfile,'(A1)') dtFileReadLine !Dummy don't bother
  do iz = 1, nnz
     read(fileUMeanProfile,*) umeanProfile(iz), vmeanProfile(iz), wmeanProfile(iz)
  end do
  close(fileUMeanProfile)

  open(fileTMeanProfile,file="txym")
  read(fileTMeanProfile,'(A1)') dtFileReadLine !Dummy don't bother
  do iz = 1, nnz
     read(fileTMeanProfile,*) tmeanProfile(iz), wmeanProfile(iz)
  end do
  close(fileTMeanProfile)

  !Writing list of points on inlet face
  open(unit=11,file="top/points")
  write(11,*)"/*--------------------------------*- C++ -*----------------------------------*\ "
  write(11,*)"| =========                 |                                                 | "
  write(11,*)"| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           | "
  write(11,*)"|  \\    /   O peration     | Version:  1.6                                   | "
  write(11,*)"|   \\  /    A nd           | Web:      http://www.OpenFOAM.org               | "
  write(11,*)"|    \\/     M anipulation  |                                                 | "
  write(11,*)"\*----------------------------------------------------------------------------*/"
  write(11,*)"FoamFile"
  write(11,*)"{"
  write(11,*)"    version     2.0;"
  write(11,*)"    format      ascii;"
  write(11,*)"    class       vectorField;"
  write(11,*)"    object      points;"
  write(11,*)"}"
  write(11,*)"// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //"
  write(11,*)"("
  do j=1,84
     do i= 1, 84
        write(11,*)"(", (i-1)*dx, " ", (j-1)*dy, " ", 250.0,")"
     end do
  end do
  write(11,*)")"
  write(11,*)"// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // "
  close(11)

  if(meanOrTurb == 1) then
     open(unit=fileXY, file="../viz.abl.094999_102000.xy.data", form="unformatted", access="direct", recl=nvar*nnx*nny*sizeOfReal)
  end if
  open(fileDT,file="dt",form="formatted")
  read(fileDT,'(A1)') dtFileReadLine
  read(fileDT,'(A1)') dtFileReadLine
  read(fileDT,'(A1)') dtFileReadLine
  read(fileDT,'(A1)') dtFileReadLine
  do k=1,1678    !Start at 38774.3s
     read(fileDT,'(e15.6,e15.6)') timeStart, dtLES
  end do
  tLESnew = timeStart
  do fileCounter=1678,nt+1678
     if(meanOrTurb == 1) then
        k = 32
        read(fileXY,rec=(fileCounter-1)*nnz + k) pA_xy        
        u(:,:) = pA_xy(1,:,:) + 7.5
        v(:,:) = pA_xy(2,:,:)
        w(:,:,1) = pA_xy(3,:,:)
        t(:,:) = pA_xy(4,:,:) + 300.0
        e(:,:) = pA_xy(5,:,:)
        read(fileXY,rec=(fileCounter-1)*nnz + k-1) pA_xy        
        w(:,:,0) = pA_xy(3,:,:)
     end if
     write(*,'(e15.6,e15.6)') tLESnew, dtLES
     write(command,'(F07.1)') tLESnew-timeStart
     call system("mkdir "//"top/"//trim(adjustl(command)))
     open(unit=11,file="top/"//trim(adjustl(command))//"/U")
     write(11,*)"/*--------------------------------*- C++ -*----------------------------------*\ "
     write(11,*)"| =========                 |                                                 | "
     write(11,*)"| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           | "
     write(11,*)"|  \\    /   O peration     | Version:  1.6                                   | "
     write(11,*)"|   \\  /    A nd           | Web:      http://www.OpenFOAM.org               | "
     write(11,*)"|    \\/     M anipulation  |                                                 | "
     write(11,*)"\*---------------------------------------------------------------------------*/ "
     write(11,*)"FoamFile"
     write(11,*)"{"
     write(11,*)"    version     2.0;"
     write(11,*)"    format      ascii;"
     write(11,*)"    class       vectorAverageField;"
     write(11,*)"    object      values;"
     write(11,*)"}"
     write(11,*)"// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // "
   
     write(11,*) "(0 0 0)"

     write(11,*) (84)*(84)
     write(11,*)"("
     do j=254-42+1,254+42
        yLocCur = (j-1)*dy     
        do i=251-42+1, 251+42 !Corresponds to turbine at 1636-1762m
           xLocCur = (i-1)*dx
           call nearestPoints(dble(xLocCur*cos(yawAngle) - yLocCur*sin(yawAngle)), dble(xLocCur*sin(yawAngle) + yLocCur*cos(yawAngle)), ixl, ixu, iyl, iyu, wxlyl, wxuyl, wxlyu, wxuyu)
           uRot = u(ixl,iyl)*wxlyl + u(ixu,iyl)*wxuyl + u(ixl,iyu)*wxlyu + u(ixu,iyu)*wxuyu
           vRot = v(ixl,iyl)*wxlyl + v(ixu,iyl)*wxuyl + v(ixl,iyu)*wxlyu + v(ixu,iyu)*wxuyu
           wRot = 0.5*(w(ixl,iyl,1)+w(ixl,iyl,0))*wxlyl + 0.5*(w(ixu,iyl,1)+w(ixu,iyl,0))*wxuyl + 0.5*(w(ixl,iyu,1)+w(ixl,iyu,0))*wxlyu + 0.5*(w(ixu,iyu,1)+w(ixu,iyu,0))*wxuyu
           if(meanOrTurb == 1) then
              write(11,*)"(", uRot*cos(yawAngle)+vRot*sin(yawAngle), " ", -uRot*sin(yawAngle)+vRot*cos(yawAngle), " ", wRot, ")"
           else
              write(11,*)"(",umeanProfile(k)*cos(yawAngle)+vmeanProfile(k)*sin(yawAngle), " ", -umeanProfile(k)*sin(yawAngle)+vmeanProfile(k)*cos(yawAngle), " ", 0.0, ")"
           end if
        end do
     end do
     write(11,*)")"
     close(11)

     open(unit=11,file="top/"//trim(adjustl(command))//"/T")
     write(11,*)"/*--------------------------------*- C++ -*----------------------------------*\ "
     write(11,*)"| =========                 |                                                 | "
     write(11,*)"| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           | "
     write(11,*)"|  \\    /   O peration     | Version:  1.6                                   | "
     write(11,*)"|   \\  /    A nd           | Web:      http://www.OpenFOAM.org               | "
     write(11,*)"|    \\/     M anipulation  |                                                 | "
     write(11,*)"\*---------------------------------------------------------------------------*/ "
     write(11,*)"FoamFile"
     write(11,*)"{"
     write(11,*)"    version     2.0;"
     write(11,*)"    format      ascii;"
     write(11,*)"    class       scalarField;"
     write(11,*)"    object      values;"
     write(11,*)"}"
     write(11,*)"// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //"
     write(11,*)"0"
     write(11,*) (84)*(84)
     write(11,*)"("
     do j=254-42+1,254+42
        yLocCur = (j-1)*dy     
        do i=251-42+1, 251+42 !Corresponds to turbine at 1636-1762m
           xLocCur = (i-1)*dx
           call nearestPoints(dble(xLocCur*cos(yawAngle) - yLocCur*sin(yawAngle)), dble(xLocCur*sin(yawAngle) + yLocCur*cos(yawAngle)), ixl, ixu, iyl, iyu, wxlyl, wxuyl, wxlyu, wxuyu)
           tRot = t(ixl,iyl)*wxlyl + t(ixu,iyl)*wxuyl + t(ixl,iyu)*wxlyu + t(ixu,iyu)*wxuyu 
           if (meanOrTurb == 1) then
              write(11,*) tRot
           else
              write(11,*) sum(t(:,:))/real(nxy)
           end if
        end do
     end do
     write(11,*)")"
     close(11)

     open(unit=11,file="top/"//trim(adjustl(command))//"/k")
     write(11,*)"/*--------------------------------*- C++ -*----------------------------------*\ "
     write(11,*)"| =========                 |                                                 | "
     write(11,*)"| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           | "
     write(11,*)"|  \\    /   O peration     | Version:  1.6                                   | "
     write(11,*)"|   \\  /    A nd           | Web:      http://www.OpenFOAM.org               | "
     write(11,*)"|    \\/     M anipulation  |                                                 | "
     write(11,*)"\*---------------------------------------------------------------------------*/ "
     write(11,*)"FoamFile"
     write(11,*)"{"
     write(11,*)"    version     2.0;"
     write(11,*)"    format      ascii;"
     write(11,*)"    class       scalarAverageField;"
     write(11,*)"    object      values;"
     write(11,*)"}"
     write(11,*)"// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // "
     write(11,*)"0"
     write(11,*) 84*wtExtentZ
     write(11,*)"("
     do j=254-42+1,254+42
        do i=251-42+1, 251+42 !Corresponds to turbine at 1636-1762m
           xLocCur = (i-1)*dx
           call nearestPoints(dble(xLocCur*cos(yawAngle) - yLocCur*sin(yawAngle)), dble(xLocCur*sin(yawAngle) + yLocCur*cos(yawAngle)), ixl, ixu, iyl, iyu, wxlyl, wxuyl, wxlyu, wxuyu)
           eRot = e(ixl,iyl)*wxlyl + e(ixu,iyl)*wxuyl + e(ixl,iyu)*wxlyu + e(ixu,iyu)*wxuyu 
           if (meanOrTurb == 1) then
              write(11,*) eRot
           else
              write(11,*) sum(e(:,:))/real(nxy)
           end if
        end do
     end do
     write(11,*)")"
     close(11)
     
     read(fileDT,'(e15.6,e15.6)') tLESnew, dtLES

  end do

end program writeOpenFOAMBCs
