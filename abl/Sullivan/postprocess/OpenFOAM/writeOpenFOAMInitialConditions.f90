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

program writeOpenFOAMInitialConditions
  !This is a program to read the wtbox files output by the LES code and write them in the format required by OpenFOAM for the initial conditions
  implicit none

  integer :: i,j,k,it,fileCounter,iz !Counter variables
  
  real :: xl, yl         ! x, y domain sizes
  real :: dx, dy, dz         ! x, y grid lengths
  integer :: nnx, nny, nnz      ! x, y grid dimensions
  integer :: nvar=5 !The number of variables stored in a plane at every time step
  integer :: nxy !nnx*nny
  integer :: wtExtentX, wtExtentY, wtExtentZ !The streamwise, spanwise and vertical extent of the domain
  integer :: nt !Number of time steps
  real :: dt  !Time step
  real(kind=4), dimension(:,:,:), allocatable :: pA_xy !Plane Arrays 
  real(kind=4), dimension(:,:,:), allocatable :: u,v,w,t,e !Plane Arrays 
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
  double precision :: tmpX, tmpY, tmpZ , tmpXM, tmpYM, tmpZM !Temporary variables


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

  wtExtentZ = 50
  wtExtentY = 416
  wtExtentX = 416
  
  nxy = nnx*nny
  nt =  1500
  
  allocate(pA_xy(nvar,nnx,nny))
  allocate(u(nnx,nny,nnz))
  allocate(v(nnx,nny,nnz))
  allocate(w(nnx,nny,0:nnz))
  allocate(t(nnx,nny,nnz))
  allocate(e(nnx,nny,nnz))
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

  if(meanOrTurb == 1) then
     open(unit=fileXY, file="../viz.abl.094999_102000.xy.data", form="unformatted", access="direct", recl=nvar*nnx*nny*sizeOfReal)
  end if
  open(fileDT,file="dt",form="formatted")
  read(fileDT,'(A1)') dtFileReadLine
  read(fileDT,'(A1)') dtFileReadLine
  read(fileDT,'(A1)') dtFileReadLine
  read(fileDT,'(A1)') dtFileReadLine
  read(fileDT,'(e15.6,e15.6)') timeStart, dtLES
  tLESnew = timeStart
  do fileCounter=1,nt+1
     if(meanOrTurb == 1) then
        do k = 1, nnz
           read(fileXY,rec=(fileCounter-1)*nnz + k) pA_xy        
           u(:,:,k) = pA_xy(1,:,:) + 7.5
           v(:,:,k) = pA_xy(2,:,:)
           w(:,:,k) = pA_xy(3,:,:)
           t(:,:,k) = pA_xy(4,:,:) + 300.0
           e(:,:,k) = pA_xy(5,:,:)
        end do
     end if
     write(*,'(e15.6,e15.6)') tLESnew, dtLES
     write(command,'(F07.1)') tLESnew-timeStart
     call system("mkdir "//"initialConditions/"//trim(adjustl(command)))
     open(unit=11,file="initialConditions/"//trim(adjustl(command))//"/U")
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
     write(11,*)"    class       volVectorField;"
     write(11,*)"    location    ",tLESnew-timeStart,"; "
     write(11,*)"    object      U;"
     write(11,*)"}"
     write(11,*)"// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // "
     write(11,*)" "
     write(11,*)"dimensions      [0 1 -1 0 0 0 0]; "
     write(11,*)" "
     write(11,*)"internalField   nonuniform List<vector> "
     write(11,*) wtExtentX*wtExtentY*wtExtentZ
     write(11,*)"("

     open(unit=12,file="initialConditions/"//trim(adjustl(command))//"/tke")
     write(12,*)"/*--------------------------------*- C++ -*----------------------------------*\ "
     write(12,*)"| =========                 |                                                 | "
     write(12,*)"| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           | "
     write(12,*)"|  \\    /   O peration     | Version:  1.6                                   | "
     write(12,*)"|   \\  /    A nd           | Web:      http://www.OpenFOAM.org               | "
     write(12,*)"|    \\/     M anipulation  |                                                 | "
     write(12,*)"\*---------------------------------------------------------------------------*/ "
     write(12,*)"FoamFile"
     write(12,*)"{"
     write(12,*)"    version     2.0;"
     write(12,*)"    format      ascii;"
     write(12,*)"    class       volScalarField;"
     write(12,*)"    location    ",tLESnew-timeStart,"; "
     write(12,*)"    object      tke;"
     write(12,*)"}"
     write(12,*)"// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // "
     write(12,*)" "
     write(12,*)"dimensions      [0 2 -2 0 0 0 0]; "
     write(12,*)" "
     write(12,*)"internalField   nonuniform List<scalar> "
     write(12,*) wtExtentX*wtExtentY*wtExtentZ
     write(12,*)"("


     do k=1,wtExtentZ
        do j=1,wtExtentY
           yLocCur = (j-1)*dy
           do i=1,wtExtentX
              xLocCur = (i-1)*dx              
              call nearestPoints(dble(xLocCur*cos(yawAngle) - yLocCur*sin(yawAngle)), dble(xLocCur*sin(yawAngle) + yLocCur*cos(yawAngle)), ixl, ixu, iyl, iyu, wxlyl, wxuyl, wxlyu, wxuyu)
              uRot = u(ixl,iyl,k)*wxlyl + u(ixu,iyl,k)*wxuyl + u(ixl,iyu,k)*wxlyu + u(ixu,iyu,k)*wxuyu
              vRot = v(ixl,iyl,k)*wxlyl + v(ixu,iyl,k)*wxuyl + v(ixl,iyu,k)*wxlyu + v(ixu,iyu,k)*wxuyu
              wRot = 0.5*(w(ixl,iyl,k)+w(ixl,iyl,k-1))*wxlyl + 0.5*(w(ixu,iyl,k)+w(ixu,iyl,k-1))*wxuyl + 0.5*(w(ixl,iyu,k)+w(ixl,iyu,k-1))*wxlyu + 0.5*(w(ixu,iyu,k)+w(ixu,iyu,k-1))*wxuyu
              tmpX = uRot*cos(yawAngle)+vRot*sin(yawAngle)
              tmpY = -uRot*sin(yawAngle)+vRot*cos(yawAngle)
              tmpZ = wRot
              tmpXM = umeanProfile(k)*cos(yawAngle)+vmeanProfile(k)*sin(yawAngle)
              tmpYM = -umeanProfile(k)*sin(yawAngle)+vmeanProfile(k)*cos(yawAngle)
              tmpZM = 0.0
              
              if(meanOrTurb == 1) then
                 write(11,*)"(", uRot*cos(yawAngle)+vRot*sin(yawAngle), " ", -uRot*sin(yawAngle)+vRot*cos(yawAngle), " ", wRot, ")"
                 write(12,*) 0.5 * ( (tmpX-tmpXM)*(tmpX-tmpXM) + (tmpY-tmpYM)*(tmpY-tmpYM) + (tmpZ-tmpZM)*(tmpZ-tmpZM) )
              else
                 write(11,*)"(",umeanProfile(k)*cos(yawAngle)+vmeanProfile(k)*sin(yawAngle), " ", -umeanProfile(k)*sin(yawAngle)+vmeanProfile(k)*cos(yawAngle), " ", 0.0, ")"
              end if
           end do
        end do
     end do
     write(11,*)");"
     write(11,*)"boundaryField "
     write(11,*)"{ "
     write(11,*)"bottomWall"
     write(11,*)"{"
     write(11,*)"type            fixedValue;"
     write(11,*)"value           uniform (0 0 0);"
     write(11,*)"}"
     write(11,*)"top"
     write(11,*)"{"
     write(11,*)"type            zeroGradient;"
     write(11,*)"}"
     write(11,*)"inlet "
     write(11,*)"{"
     write(11,*)"type              timeVaryingMappedFixedValue;"
     write(11,*)"setAverage        off;"
     write(11,*)"}"
     write(11,*)"outlet"
     write(11,*)"{"
     write(11,*)"type              zeroGradient;"
     write(11,*)"}"
     write(11,*)"sideInlet"
     write(11,*)"{"
     write(11,*)"type              zeroGradient;"
     write(11,*)"}"
     write(11,*)"sideOutlet"
     write(11,*)"{"
     write(11,*)"type              zeroGradient;"
     write(11,*)"}"
     write(11,*)"}"
     close(11)


     write(12,*)");"
     write(12,*)"boundaryField "
     write(12,*)"{ "
     write(12,*)"bottomWall"
     write(12,*)"{"
     write(12,*)"type            fixedValue;"
     write(12,*)"value           uniform 0;"
     write(12,*)"}"
     write(12,*)"top"
     write(12,*)"{"
     write(12,*)"type            zeroGradient;"
     write(12,*)"}"
     write(12,*)"inlet "
     write(12,*)"{"
     write(12,*)"type              zeroGradient;"
     write(12,*)"}"
     write(12,*)"outlet"
     write(12,*)"{"
     write(12,*)"type              zeroGradient;"
     write(12,*)"}"
     write(12,*)"sideInlet"
     write(12,*)"{"
     write(12,*)"type              zeroGradient;"
     write(12,*)"}"
     write(12,*)"sideOutlet"
     write(12,*)"{"
     write(12,*)"type              zeroGradient;"
     write(12,*)"}"
     write(12,*)"}"
     close(12)



     open(unit=11,file="initialConditions/"//trim(adjustl(command))//"/T")
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
     write(11,*)"    class       volScalarField;"
     write(11,*)"    location    ",tLESnew-timeStart,";"
     write(11,*)"    object      T;"
     write(11,*)"}"
     write(11,*)"// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // "
     write(11,*)" "
     write(11,*)"dimensions      [0 0 0 1 0 0 0]; "
     write(11,*)" "
     write(11,*)"internalField   nonuniform List<scalar> "
     write(11,*) wtExtentX*wtExtentY*wtExtentZ
     write(11,*)"("
     do k=1,wtExtentZ
        do j=1,wtExtentY
           yLocCur = (j-1)*dy
           do i=1,wtExtentX
              xLocCur = (i-1)*dx              
              call nearestPoints(dble(xLocCur*cos(yawAngle) - yLocCur*sin(yawAngle)), dble(xLocCur*sin(yawAngle) + yLocCur*cos(yawAngle)), ixl, ixu, iyl, iyu, wxlyl, wxuyl, wxlyu, wxuyu)
              tRot = t(ixl,iyl,k)*wxlyl + t(ixu,iyl,k)*wxuyl + t(ixl,iyu,k)*wxlyu + t(ixu,iyu,k)*wxuyu 
              if (meanOrTurb == 1) then
                 write(11,*) tRot
              else
                 write(11,*) sum(t(:,:,k))/real(nxy)
              end if
           end do
        end do
     end do
     write(11,*)");"
     write(11,*)"boundaryField "
     write(11,*)"{ "
     write(11,*)"bottomWall"
     write(11,*)"{"
     write(11,*)"type            fixedValue;"
     write(11,*)"value           uniform 303.79;"
     write(11,*)"}"
     write(11,*)"top"
     write(11,*)"{"
     write(11,*)"type            zeroGradient;"
     write(11,*)"}"
     write(11,*)"inlet "
     write(11,*)"{"
     write(11,*)"type              timeVaryingMappedFixedValue;"
     write(11,*)"setAverage        off;"
     write(11,*)"}"
     write(11,*)"outlet"
     write(11,*)"{"
     write(11,*)"type              zeroGradient;"
     write(11,*)"}"
     write(11,*)"sideInlet"
     write(11,*)"{"
     write(11,*)"type              zeroGradient;"
     write(11,*)"}"
     write(11,*)"sideOutlet"
     write(11,*)"{"
     write(11,*)"type              zeroGradient;"
     write(11,*)"}"
     write(11,*)"}"
     close(11)

     open(unit=11,file="initialConditions/"//trim(adjustl(command))//"/k")
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
     write(11,*)"    class       volScalarField;"
     write(11,*)"    location    ", tLESnew-timeStart,";"
     write(11,*)"    object      k;"
     write(11,*)"}"
     write(11,*)"// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // "
     write(11,*)" "
     write(11,*)"dimensions      [0 2 -2 0 0 0 0]; "
     write(11,*)" "
     write(11,*)"internalField   nonuniform List<scalar> "
     write(11,*) wtExtentX*wtExtentY*wtExtentZ
     write(11,*)"("
     do k=1,wtExtentZ
        do j=1,wtExtentY
           yLocCur = (j-1)*dy
           do i=1,wtExtentX
              xLocCur = (i-1)*dx              
              call nearestPoints(dble(xLocCur*cos(yawAngle) - yLocCur*sin(yawAngle)), dble(xLocCur*sin(yawAngle) + yLocCur*cos(yawAngle)), ixl, ixu, iyl, iyu, wxlyl, wxuyl, wxlyu, wxuyu)
              eRot = e(ixl,iyl,k)*wxlyl + e(ixu,iyl,k)*wxuyl + e(ixl,iyu,k)*wxlyu + e(ixu,iyu,k)*wxuyu 
              if (meanOrTurb == 1) then
                 write(11,*) eRot
              else
                 write(11,*) sum(e(:,:,k))/real(nxy)
              end if
           end do
        end do
     end do
     write(11,*)");"
     write(11,*)"boundaryField "
     write(11,*)"{ "
     write(11,*)"bottomWall"
     write(11,*)"{"
     write(11,*)"type            zeroGradient;"
     write(11,*)"}"
     write(11,*)"top"
     write(11,*)"{"
     write(11,*)"type            zeroGradient;"
     write(11,*)"}"
     write(11,*)"inlet "
     write(11,*)"{"
     write(11,*)"type              timeVaryingMappedFixedValue;"
     write(11,*)"setAverage        off;"
     write(11,*)"}"
     write(11,*)"outlet"
     write(11,*)"{"
     write(11,*)"type              zeroGradient;"
     write(11,*)"}"
     write(11,*)"sideInlet"
     write(11,*)"{"
     write(11,*)"type              zeroGradient;"
     write(11,*)"}"
     write(11,*)"sideOutlet"
     write(11,*)"{"
     write(11,*)"type              zeroGradient;"
     write(11,*)"}"
     write(11,*)"}"
     close(11)
              
     read(fileDT,'(e15.6,e15.6)') tLESnew, dtLES
  end do

end program writeOpenFOAMInitialConditions
