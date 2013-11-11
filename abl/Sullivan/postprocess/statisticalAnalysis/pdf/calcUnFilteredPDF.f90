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

program calcUnfilteredPDF

  !!!This is a program to compute the pdf at a plane.
  
  implicit none

  integer :: i,j,k,it,wtx,wty,fileCounter,iz !Counter variables
  
  real :: xl, yl         ! x, y domain sizes
  real :: dx, dy         ! x, y grid lengths
  integer :: nnx, nny, nnz      ! x, y grid dimensions
  real(kind=8), dimension(:), allocatable :: xArr, yArr !Array of x and y locations on the grid
  integer :: nvar=5 !The number of variables stored in a plane at every time step
  integer :: nxy !nnx*nny
  integer :: nt !Number of time steps
  real :: t !Current time
  real :: dt  !Time step
  real(kind=4), dimension(:,:,:), allocatable :: pA_xy !Plane Arrays 
  real(kind=4), dimension(:,:), allocatable :: u,v,w !Plane Arrays 
  real(kind=4), dimension(:,:), allocatable :: uRot,vRot,wRot !Plane Arrays 
  real(kind=4), dimension(:), allocatable :: uWrite,vWrite,wWrite !Plane Arrays in the format to be written into vtr files
  real, allocatable, dimension(:) :: umean, vmean, wmean
  real, allocatable, dimension(:) :: umeanProfile, vmeanProfile, wmeanProfile !Profile of mean velocities from postprocessing of his.mp files
  real, allocatable, dimension(:) :: uvarProfile, vvarProfile, wvarProfile !Profile of velocity variances from postprocessing of his.mp files
  integer :: zLevel !The z level at which the data is to be visualized.
  integer :: sizeOfReal=1 !The size of real in words on this system
  integer :: fileXY, fileDT  !The index of the files used to open the "xy.data" and the "dt" files

  !Variables to interpolate to yawed co-ordinate system
  real(kind=4) :: yawAngle
  real(kind=4) :: xLocCur, yLocCur, zLocCur
  integer :: ixl, ixu, iyl, iyu !The integer locations of the 4 points used to interpolate the velocity from the LES of ABL field to the reqd. point in the FAST field
  double precision :: wxlyl, wxuyl, wxlyu, wxuyu !The weights of the 4 points used to interpolate the velocity from the LES of ABL field to the reqd. point in the FAST field

  character *20  buffer  !Buffer to read in command line arguments
  character *40  writeFilename, readFileName

  !PDF
  real :: umin, umax, uvar, wmin, wmax, wvar
  real, allocatable, dimension(:) :: updf, wpdf
  real :: binSizeU, binSizeW
  integer :: binLoc, nbinsU, nbinsW

  
  yawAngle = 21.5*3.14159265359/180.0
  xl = 5120.0
  yl = 5120.0
  nnx = 768
  nny = 768
  dx = xl/nnx
  dy = yl/nny
  allocate(xArr(nnx))
  allocate(yArr(nny))
  do i=1,nnx
     xArr(i) = dble(i-1)*dx
     yArr(i) = dble(i-1)*dy
  end do
  nnz = 50
  sizeOfReal = 1
  fileXY = 95
  fileDT = 91
  
  nxy = nnx*nny
!  nt = 10
  nt =  5000

  allocate(pA_xy(nvar,nnx,nny))
  allocate(u(nnx,nny))
  allocate(v(nnx,nny))
  allocate(w(nnx,nny))
  allocate(uRot(nnx,nny))
  allocate(vRot(nnx,nny))
  allocate(wRot(nnx,nny))

  call getarg(1,buffer)
  read(buffer,*) zLevel
  write(*,*) zLevel

  !Initialize PDF variables
  umin = 0.0
  umax = 22.0
  wmin = -5.0
  wmax = 5.0
  binSizeW = 0.02
  binSizeU = 0.05
  nBinsU = int((umax-umin)/binSizeU)+1
  nBinsW = int((wmax-wmin)/binSizeW)+1
  allocate(uPDF(0:nBinsU))
  allocate(wPDF(0:nBinsW))
  uPDF = 0.0
  wPDF = 0.0

  open(unit=fileXY, file="../../viz.abl.094999_102000.xy.data", form="unformatted", access="direct", recl=nvar*nnx*nny*sizeOfReal)
  do fileCounter=1,nt
     write(*,*) 'Reading record number ', (fileCounter-1)*nnz + zLevel
     read(fileXY,rec=(fileCounter-1)*nnz + zLevel) pA_xy
     u = pA_xy(1,:,:) + 7.5
     v = pA_xy(2,:,:)
     w = pA_xy(3,:,:) !Temporary
     if(zLevel /= 1) then
        read(fileXY,rec=(fileCounter-1)*nnz + zLevel-1) pA_xy
        w = 0.5 * ( w + pA_xy(3,:,:))
     else
        w = 0.5 * ( w + 0.0)        
     end if
     read(fileXY,rec=(fileCounter-1)*nnz + zLevel) pA_xy

     umean = sum(u)/real(nxy)
     vmean = sum(v)/real(nxy)
     wmean = sum(w)/real(nxy)

     !Perform linear interpolation
     do j = 1,nny
        yLocCur = (j-1)*dy
        do i = 1,nnx
           xLocCur = (i-1)*dx
           call nearestPoints(dble(xLocCur*cos(yawAngle) - yLocCur*sin(yawAngle)), dble(xLocCur*sin(yawAngle) + yLocCur*cos(yawAngle)), ixl, ixu, iyl, iyu, wxlyl, wxuyl, wxlyu, wxuyu)
           uRot(i,j) = u(ixl,iyl)*wxlyl + u(ixu,iyl)*wxuyl + u(ixl,iyu)*wxlyu + u(ixu,iyu)*wxuyu 
           vRot(i,j) = v(ixl,iyl)*wxlyl + v(ixu,iyl)*wxuyl + v(ixl,iyu)*wxlyu + v(ixu,iyu)*wxuyu 
           wRot(i,j) = w(ixl,iyl)*wxlyl + w(ixu,iyl)*wxuyl + w(ixl,iyu)*wxlyu + w(ixu,iyu)*wxuyu
        end do
     end do
!      !Now rotate the velocity vector using 'u' as a temporary vector
     u = uRot*cos(yawAngle) + vRot*sin(yawAngle)
     vRot = -uRot*sin(yawAngle) + vRot*cos(yawAngle)
     uRot = u
     
     do j=1,nny
        do i=1,nnx
           binLoc = int((uRot(i,j)-umin)/binSizeU)
           updf(binLoc) = updf(binLoc) + 1
           binLoc = int((wRot(i,j)-wmin)/binSizeW)
           wpdf(binLoc) = wpdf(binLoc) + 1
        end do
     end do

  end do
  close(fileXY)

  updf = updf/(real(nt*nxy))
  wpdf = wpdf/(real(nt*nxy))

  write(writeFileName,'(I0.2)') zLevel
  call system("mkdir -p "//"zLevel"//trim(adjustl(writeFileName)))

  open(unit=23,file="zLevel"//trim(adjustl(writeFileName))//"/uUnfilteredPDF")
  write(23,"(A)",advance="no")'#  Bin , u unfiltered PDF'
  write(23,*) 
  do i=0,nbinsU
     write(23,"(F,F)",advance="no") umin+0.5*binSizeU + i*binSizeU, updf(i)/binSizeU
     write(23,*) 
  end do
  close(23)
  
  open(unit=23,file="zLevel"//trim(adjustl(writeFileName))//"/wUnfilteredPDF")
  write(23,"(A)",advance="no")'#  Bin , w unfiltered PDF'
  write(23,*) 
  do i=0,nbinsW
     write(23,"(F,F)",advance="no") wmin+0.5*binSizeW + i*binSizeW, wpdf(i)/binSizeW
     write(23,*) 
  end do
  close(23)

  
end program calcUnfilteredPDF
