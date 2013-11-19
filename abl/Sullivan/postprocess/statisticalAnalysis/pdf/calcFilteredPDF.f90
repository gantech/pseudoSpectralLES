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

program calcFilteredPDF

  !!!This is a program to compute the filtered pdf and the conditional means at a plane depending on thresholds in the pdf.
  
  implicit none
  include "/opt/fftw/3.3.0.0/x86_64/include/fftw3.f"

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
  double precision :: umean, vmean, wmean
  double precision :: uMeanFromProfile, vMeanFromProfile !The mean velocity from the profile.
  integer :: zLevel !The z level at which the data is to be visualized.
  integer :: sizeOfReal=1 !The size of real in words on this system
  integer :: fileXY  !The index of the files used to open the "xy.data" file
  integer :: fileUmeanProfile !The index of the file used to open the "uxym" file
  integer :: filewSpectrumPeakWavenumberProfile !The index of the file used to open the "wSpectrumPeakWavenumber" file
  character :: dtFileReadLine

  !Variables to interpolate to yawed co-ordinate system
  real(kind=4) :: yawAngle
  real(kind=4) :: xLocCur, yLocCur, zLocCur
  integer :: ixl, ixu, iyl, iyu !The integer locations of the 4 points used to interpolate the velocity from the LES of ABL field to the reqd. point in the FAST field
  double precision :: wxlyl, wxuyl, wxlyu, wxuyu !The weights of the 4 points used to interpolate the velocity from the LES of ABL field to the reqd. point in the FAST field

  !FFT variables
  double precision, allocatable, dimension(:,:) :: fft_in !Temporary array to store the input to the FFTW function
  double complex, allocatable, dimension(:,:) :: fft_out !Temporary array to store the output from the FFTW function
  real, allocatable, dimension(:) :: waveN ! Wavenumber
  real :: deltaWaveN !The difference between two successive wave numbers
  real :: pi = 3.14159265358979
  integer :: r ! Distance of 2d wave number from orgin
  integer*8 plan, plan_inv ; !Plan to create the FFT and it's inverse
  integer :: stat
  
  integer :: u_cutoff, v_cutoff, w_cutoff !Filter cutoff wavenumber for each velocity component

  character(10) ::  fileName !Output file name
  character *20  buffer  !Buffer to read in command line arguments
  character *40  writeFilename, readFileName

  !Thresholding to define downdraft, mean and updraft in horizontal and vertical velocity
  double precision :: u_ud_cdf, u_dd_cdf, w_ud_cdf, w_dd_cdf !Thresholds based on the Cumulative Density Function of filtered data
  double precision :: uDD, uUD !DD - DownDraft, UD - UpDraft
  double precision :: wDD, wUD !DD - DownDraft, UD - UpDraft
  integer*8 :: uDDcount, uUDcount !Integer count of number of points in these structures
  integer*8 :: wDDcount, wUDcount !Integer count of number of points in these structures
  double precision :: umeanCondUD, vmeanCondUD, wmeanCondUD !Conditional means in updraft
  double precision :: umeanCondDD, vmeanCondDD, wmeanCondDD !Conditional means in downdraft
  double precision :: umeanCondLS, vmeanCondLS, wmeanCondLS !Conditional means in low speed streaks
  double precision :: umeanCondHS, vmeanCondHS, wmeanCondHS !Conditional means in high speed streaks

  !PDF
  double precision :: umin, umax, uvar, wmin, wmax, wvar
  double precision, allocatable, dimension(:) :: updf, wpdf
  double precision :: binSizeU, binSizeW
  integer :: binLoc, nbinsU, nbinsW

  
  yawAngle = 21.5*3.14159265359/180.0
  xl = 5120.0
  yl = 5120.0
  nnx = 768
  nny = 768
  dx = xl/nnx
  dy = yl/nny
  nnz = 50
  sizeOfReal = 1
  fileXY = 95
  fileUMeanProfile = 83
  filewSpectrumPeakWavenumberProfile = 93

!  nt = 100
  nt =  5000

  allocate(pA_xy(nvar,nnx,nny))
  allocate(u(nnx,nny))
  allocate(v(nnx,nny))
  allocate(w(nnx,nny))
  allocate(uRot(nnx,nny))
  allocate(vRot(nnx,nny))
  allocate(wRot(nnx,nny))

  write(*,*) 'Usage: ExecutableName CutOffWaveNumber u_ud_cdf u_dd_cdf w_ud_cdf w_dd_cdf'
  call getarg(1,buffer)
  read(buffer,*) zLevel
  write(*,*) 'Working on zLevel = ', zLevel

  call getarg(2,buffer)
  read(buffer,*) u_ud_cdf

  call getarg(3,buffer)
  read(buffer,*) u_dd_cdf

  call getarg(4,buffer)
  read(buffer,*) w_ud_cdf

  call getarg(5,buffer)
  read(buffer,*) w_dd_cdf

  !Read in umean from the profile
  open(fileUMeanProfile,file="uxym")
  read(fileUMeanProfile,'(A1)') dtFileReadLine !Dummy don't bother
  do iz = 1, zLevel
     read(fileUMeanProfile,*) uMeanFromProfile, vMeanFromProfile, wmean
     write(*,*) iz, uMeanFromProfile, vMeanFromProfile, wmean
  end do
  close(fileUMeanProfile)
  umean = uMeanFromProfile*cos(yawAngle) + vMeanFromProfile*sin(yawAngle)
  vMeanFromProfile = -uMeanFromProfile*sin(yawAngle) + vMeanFromProfile*cos(yawAngle)
  uMeanFromProfile =  umean

  !Read in umean from the profile
  open(filewSpectrumPeakWavenumberProfile,file="wSpectrumPeakWavenumber")
  read(filewSpectrumPeakWavenumberProfile,'(A1)') dtFileReadLine !Dummy don't bother
  do iz = 1, zLevel
     read(filewSpectrumPeakWavenumberProfile,*) k, u_cutoff
  end do
  close(filewSpectrumPeakWavenumberProfile)
  u_cutoff = 3*u_cutoff !Cut off at 3 times the peak in the w spectrum
  v_cutoff = u_cutoff
  w_cutoff = u_cutoff
  write(*,*) 'Filtering at a wave number of ', u_cutoff

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

  if (u_cutoff .ne. 0) then
     !Initialise FFT variables
     write(*,*) 'Initializing and allocating FFTW variables'
     !!Allocate complex arrays
     allocate(waveN(nnx))
     allocate(fft_in(nnx,nny))
     allocate(fft_out(nnx/2+1,nny))
     
     !!Wave number
     deltaWaveN = 1/xl
     do i =1,nnx/2+1
        waveN(i) = (i-1)
     end do
     do i=nnx/2+2,nnx
        waveN(i) = i-nnx-1
     end do
     
     !!Create FFTW plan
     call dfftw_plan_dft_r2c_2d(plan, nnx, nny, fft_in, fft_out, FFTW_ESTIMATE)
     call dfftw_plan_dft_c2r_2d(plan_inv, nnx, nny, fft_out, fft_in, FFTW_ESTIMATE)
     write(*,*) 'Finished initializing and allocating FFTW variables. Beginning first time loop'  
  end if

  open(unit=fileXY, file="../../viz.abl.094999_102000.xy.data", form="unformatted", access="direct", recl=nvar*nnx*nny*sizeOfReal)
  do fileCounter=1,nt
     write(*,*) 'Reading record number ', (fileCounter-1)*nnz + zLevel
     read(fileXY,rec=(fileCounter-1)*nnz + zLevel) pA_xy
     u = pA_xy(1,:,:) + 7.5
     v = pA_xy(2,:,:)
     w = pA_xy(3,:,:) !Temporary
     if(zLevel /= 0) then
        read(fileXY,rec=(fileCounter-1)*nnz + zLevel-1) pA_xy
        w = 0.5 * ( w + pA_xy(3,:,:))
     else
        w = 0.5 * ( w + 0.0)        
     end if
     read(fileXY,rec=(fileCounter-1)*nnz + zLevel) pA_xy

     umean = sum(u)/dble(dble(nnx)*dble(nny))
     vmean = sum(v)/dble(dble(nnx)*dble(nny))
     wmean = sum(w)/dble(dble(nnx)*dble(nny))

     if(u_cutoff .ne. 0) then
        !U filter
        fft_in = u - umean
        call dfftw_execute_dft_r2c(plan, fft_in, fft_out)
        fft_out = fft_out/dble(dble(nnx)*dble(nny))
        do i = 1,nnx/2+1
           do j = 1,nny
              r = nint(sqrt(waveN(i)*waveN(i) + waveN(j)*waveN(j)))
              if( r .gt. u_cutoff) then
                 fft_out(i,j) = 0
              end if
           end do
        end do
        call dfftw_execute_dft_c2r(plan_inv, fft_out, fft_in)     
        u = fft_in + umean
        
        !     write(*,*) 'Finished u filtering'
        
        !V filter
        fft_in = v - vmean
        call dfftw_execute_dft_r2c(plan, fft_in, fft_out)
        fft_out = fft_out/dble(dble(nnx)*dble(nny))
        do i = 1,nnx/2+1
           do j = 1,nny
              r = nint(sqrt(waveN(i)*waveN(i) + waveN(j)*waveN(j)))
              if( r .gt. u_cutoff) then
                 fft_out(i,j) = 0
              end if
           end do
        end do
        call dfftw_execute_dft_c2r(plan_inv, fft_out, fft_in)     
        v = fft_in + vmean
        
        !     write(*,*) 'Finished v filtering'
        
        !W filter
        fft_in = w - wmean
        call dfftw_execute_dft_r2c(plan, fft_in, fft_out)
        fft_out = fft_out/dble(dble(nnx)*dble(nny))
        do i = 1,nnx/2+1
           do j = 1,nny
              r = nint(sqrt(waveN(i)*waveN(i) + waveN(j)*waveN(j)))
              if( r .gt. w_cutoff) then
                 fft_out(i,j) = 0
              end if
           end do
        end do
        call dfftw_execute_dft_c2r(plan_inv, fft_out, fft_in)     
        w = fft_in + wmean
     end if

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

  updf = updf/(dble(dble(nt)*dble(nnx)*dble(nny)))
  wpdf = wpdf/(dble(dble(nt)*dble(nnx)*dble(nny)))

  write(writeFileName,'(I0.2)') zLevel
  call system("mkdir "//"zLevel"//trim(adjustl(writeFileName)))

  open(unit=88,file="zLevel"//trim(adjustl(writeFileName))//"/uFilteredPDF")
  write(88,"(A)")'#  Bin , u filtered PDF'
  write(88,*) 
  do i=0,nbinsU
     write(88,"(F,F)") umin+0.5*binSizeU + i*binSizeU, updf(i)/binSizeU
     write(88,*) 
  end do
  close(88)

  do i=0,nbinsU
     !     write(*,*) wmin+0.5*binSize + i*binSize, ' ', i,  ' ', wpdf(i), ' ', sum(wpdf(:i)) 
     if(  (sum(updf(:i)) .gt. u_ud_cdf) .and. (uud .eq. 0)) then
        uud = umin+0.5*binSizeU + i*binSizeU
     end if
     if(  (sum(updf(:i)) .gt. u_dd_cdf) .and. (udd .eq. 0)) then
        udd = umin+0.5*binSizeU + i*binSizeU
     end if
  end do
  write(*,*) 'uud = ' , uud
  write(*,*) 'udd = ' , udd
  

  open(unit=88,file="zLevel"//trim(adjustl(writeFileName))//"/wFilteredPDF")
  write(88,"(A)")'#  Bin , w filtered PDF'
  do i=0,nbinsW
     write(88,"(F,F)") wmin+0.5*binSizeW + i*binSizeW, wpdf(i)/binSizeW
  end do
  close(88)

  do i=0,nbinsW
     !     write(*,*) wmin+0.5*binSize + i*binSize, ' ', i,  ' ', wpdf(i), ' ', sum(wpdf(:i)) 
     if(  (sum(wpdf(:i)) .gt. w_ud_cdf) .and. (wud .eq. 0)) then
        wud = wmin+0.5*binSizeW + i*binSizeW
     end if
     if(  (sum(wpdf(:i)) .gt. w_dd_cdf) .and. (wdd .eq. 0)) then
        wdd = wmin+0.5*binSizeW + i*binSizeW
     end if
  end do
  write(*,*) 'wud = ' , wud
  write(*,*) 'wdd = ' , wdd
  

  write(*,*) 'Finished computing and writing PDFs. Beginning second time loop'  
  
  !Initializing all conditional means to zero
  umeanCondLS = 0.0
  umeanCondHS = 0.0
  vmeanCondLS = 0.0
  vmeanCondHS = 0.0
  wmeanCondLS = 0.0
  wmeanCondHS = 0.0

  umeanCondUD = 0.0
  umeanCondDD = 0.0
  vmeanCondUD = 0.0
  vmeanCondDD = 0.0
  wmeanCondUD = 0.0
  wmeanCondDD = 0.0

  
  uDDcount = 0
  uUDcount = 0
  wDDcount = 0
  wUDcount = 0
  
  do fileCounter=1,nt
     write(*,*) 'Reading record number ', (fileCounter-1)*nnz + zLevel
     read(fileXY,rec=(fileCounter-1)*nnz + zLevel) pA_xy
     u = pA_xy(1,:,:) + 7.5
     v = pA_xy(2,:,:)
     w = pA_xy(3,:,:) !Temporary
     if(zLevel /= 0) then
        read(fileXY,rec=(fileCounter-1)*nnz + zLevel-1) pA_xy
        w = 0.5 * ( w + pA_xy(3,:,:))
     else
        w = 0.5 * ( w + 0.0)        
     end if
     read(fileXY,rec=(fileCounter-1)*nnz + zLevel) pA_xy

     umean = sum(u)/dble(dble(nnx)*dble(nny))
     vmean = sum(v)/dble(dble(nnx)*dble(nny))
     wmean = sum(w)/dble(dble(nnx)*dble(nny))

     if(u_cutoff .ne. 0) then
        !U filter
        fft_in = u - umean
        call dfftw_execute_dft_r2c(plan, fft_in, fft_out)
        fft_out = fft_out/dble(dble(nnx)*dble(nny))
        do i = 1,nnx/2+1
           do j = 1,nny
              r = nint(sqrt(waveN(i)*waveN(i) + waveN(j)*waveN(j)))
              if( r .gt. u_cutoff) then
                 fft_out(i,j) = 0
              end if
           end do
        end do
        call dfftw_execute_dft_c2r(plan_inv, fft_out, fft_in)     
        u = fft_in + umean
        
        !     write(*,*) 'Finished u filtering'
        
        !V filter
        fft_in = v - vmean
        call dfftw_execute_dft_r2c(plan, fft_in, fft_out)
        fft_out = fft_out/dble(dble(nnx)*dble(nny))
        do i = 1,nnx/2+1
           do j = 1,nny
              r = nint(sqrt(waveN(i)*waveN(i) + waveN(j)*waveN(j)))
              if( r .gt. u_cutoff) then
                 fft_out(i,j) = 0
              end if
           end do
        end do
        call dfftw_execute_dft_c2r(plan_inv, fft_out, fft_in)     
        v = fft_in + vmean
        
        !     write(*,*) 'Finished v filtering'
        
        !W filter
        fft_in = w - wmean
        call dfftw_execute_dft_r2c(plan, fft_in, fft_out)
        fft_out = fft_out/dble(dble(nnx)*dble(nny))
        do i = 1,nnx/2+1
           do j = 1,nny
              r = nint(sqrt(waveN(i)*waveN(i) + waveN(j)*waveN(j)))
              if( r .gt. w_cutoff) then
                 fft_out(i,j) = 0
              end if
           end do
        end do
        call dfftw_execute_dft_c2r(plan_inv, fft_out, fft_in)     
        w = fft_in + wmean
     end if

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

     umeanCondHS = umeanCondHS + sum(uRot(:,:),uRot(:,:) .gt. uDD)
     umeanCondLS = umeanCondLS + sum(uRot(:,:),uRot(:,:) .lt. uUD)
     umeanCondDD = umeanCondDD + sum(uRot(:,:),wRot(:,:) .lt. wDD)
     umeanCondUD = umeanCondUD + sum(uRot(:,:),wRot(:,:) .gt. wUD)
     
     vmeanCondHS = vmeanCondHS + sum(vRot(:,:),uRot(:,:) .gt. uDD)
     vmeanCondLS = vmeanCondLS + sum(vRot(:,:),uRot(:,:) .lt. uUD)
     vmeanCondDD = vmeanCondDD + sum(vRot(:,:),wRot(:,:) .lt. wDD)
     vmeanCondUD = vmeanCondUD + sum(vRot(:,:),wRot(:,:) .gt. wUD)
     
     uDDcount = uDDcount + count(uRot(:,:) .gt. uDD) 
     uUDcount = uUDcount + count(uRot(:,:) .lt. uUD) 
     
     wmeanCondHS = wmeanCondHS + sum(wRot(:,:),uRot(:,:) .gt. uDD)
     wmeanCondLS = wmeanCondLS + sum(wRot(:,:),uRot(:,:) .lt. uUD)
     wmeanCondDD = wmeanCondDD + sum(wRot(:,:),wRot(:,:) .lt. wDD)
     wmeanCondUD = wmeanCondUD + sum(wRot(:,:),wRot(:,:) .gt. wUD)
     wDDcount = wDDcount + count(wRot(:,:) .lt. wDD) 
     wUDcount = wUDcount + count(wRot(:,:) .gt. wUD) 

  end do

  umeanCondDD = umeanCondDD/dble(wDDcount)
  umeanCondUD = umeanCondUD/dble(wUDcount)
  umeanCondHS = umeanCondHS/dble(uDDcount)
  umeanCondLS = umeanCondLS/dble(uUDcount)

  vmeanCondDD = vmeanCondDD/dble(wDDcount)
  vmeanCondUD = vmeanCondUD/dble(wUDcount)
  vmeanCondHS = vmeanCondHS/dble(uDDcount)
  vmeanCondLS = vmeanCondLS/dble(uUDcount)

  wmeanCondHS = wmeanCondHS/dble(uDDcount)
  wmeanCondLS = wmeanCondLS/dble(uUDcount)
  wmeanCondDD = wmeanCondDD/dble(wDDcount)
  wmeanCondUD = wmeanCondUD/dble(wUDcount)


  open(unit=89,file="zLevel"//trim(adjustl(writeFileName))//"/conditionalMeanVelProfile")
  write(89,"(A)") 'zLevel umeanUD     umeanLS   umean      umeanDD  umeanHS    wmeanLS    wmeanUD    wmean       wmeanDD      wmeanHS'
  write(89,*)
  write(89,"(I,16F)") zLevel, umeanCondUD, umeanCondLS, uMeanFromProfile, umeanCondDD, umeanCondHS,&
       vmeanCondUD, vmeanCondLS, vMeanFromProfile, vmeanCondDD, vmeanCondHS, &
       wmeanCondUD, wmeanCondLS, 0.0, wmeanCondDD, wmeanCondHS
  write(89,*)
  close(89)

  close(fileXY)


end program calcFilteredPDF
