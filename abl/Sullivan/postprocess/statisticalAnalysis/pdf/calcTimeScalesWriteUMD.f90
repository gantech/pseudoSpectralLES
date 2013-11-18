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

program calcTimeScalesWriteUMD
  !!!This is a program to analyze the data at plane of hub height. It reads in the velocity data at the hub height at each instant in time. The program is capable of doing the following analyses
  
  !1. To filter the velocity field at each instant in time using a 2D spectral filter
  !2. To get the time scale of the updrafts and downdrafts using the pdf to determine the definition of each 
  !3. It also write the UMD data based on U and W at 210 locations.
  
  implicit none
  include "/opt/fftw/3.3.0.0/x86_64/include/fftw3.f"

  integer :: i,j,k,it,wtx,wty,fileCounter,iz !Counter variables
  
  real :: xl, yl         ! x, y domain sizes
  real :: dx, dy         ! x, y grid lengths
  integer :: nnx, nny, nnz      ! x, y grid dimensions
  real(kind=8), dimension(:), allocatable :: xArr, yArr !Array of x and y locations on the grid
  integer :: nvar=5 !The number of variables stored in a plane at every time step
  integer :: wtExtentX !The spacing between wind turbines in the streamwise direction
  integer :: wtExtentY !The spanwise extent of the wind turbine in number of grid points.
  integer :: nWtx, nWty !Number of wind turbines in each direction
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
  integer :: fileXY, fileDT  !The index of the files used to open the "xy.data" and the "dt" files
  integer :: fileUmeanProfile !The index of the file used to open the "uxym" file
  real(kind=4) :: tLESnew, tLESold
  real(kind=4) :: dtLES
  integer :: itLES
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

  integer, dimension(:,:), allocatable :: UMD !Array to determine whether each point is an Updraft, Mean or Downdraft point
  real, dimension(:,:,:), allocatable :: t_u_UMD !Time scale of DownDraft, UpDraft and Mean obtained from u
  real, dimension(:,:,:), allocatable :: t_w_UMD !Time scale of DownDraft, UpDraft and Mean obtained from w
  real, dimension(:,:,:,:), allocatable :: t_u_UMD_ind !Time scale of DownDraft, UpDraft and Mean obtained from u - Individual time scales
  real, dimension(:,:,:,:), allocatable :: t_w_UMD_ind !Time scale of DownDraft, UpDraft and Mean obtained from w - Individual time scales
  integer, dimension(:,:), allocatable :: dStatU, dStatUPrev !Draft status at current and previous time step
  integer, dimension(:,:,:), allocatable :: nUMD_U !number of updraft, mean or downdraft structures at each point
  integer, dimension(:,:), allocatable :: dStatW, dStatWPrev !Draft status at current and previous time step
  integer, dimension(:,:,:), allocatable :: nUMD_W !number of updraft, mean or downdraft structures at each point

  !PDF
  double precision :: umin, umax, uvar, wmin, wmax, wvar
  double precision, allocatable, dimension(:) :: updf, wpdf
  double precision :: binSizeU, binSizeW
  integer :: binLoc, nbinsU, nbinsW

  !PDF of time scales
  real :: t_u_ud_min, t_u_ud_max, t_u_dd_min, t_u_dd_max, t_u_m_min, t_u_m_max
  real :: t_w_ud_min, t_w_ud_max, t_w_dd_min, t_w_dd_max, t_w_m_min, t_w_m_max
  real, allocatable, dimension(:) :: t_u_ud_pdf, t_u_dd_pdf, t_u_m_pdf
  real, allocatable, dimension(:) :: t_w_ud_pdf, t_w_dd_pdf, t_w_m_pdf
  real :: t_u_ud_binSize, t_u_dd_binSize, t_u_m_binSize
  real :: t_w_ud_binSize, t_w_dd_binSize, t_w_m_binSize
  integer :: t_u_ud_nBins, t_u_dd_nBins, t_u_m_nBins
  integer :: t_w_ud_nBins, t_w_dd_nBins, t_w_m_nBins

  !Wind turbine locations
  integer :: wtXloc, wtYloc

  character(10) ::  fileName !Output file name
  character *20  buffer  !Buffer to read in command line arguments
  character *40  writeFilename, readFileName

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
  fileDT = 91
  fileUMeanProfile = 83

  wtExtentY = 25
  nWtx = 15
  nWty = 14
  
!  nt = 10
  nt =  5000

  allocate(pA_xy(nvar,nnx,nny))
  allocate(u(nnx,nny))
  allocate(v(nnx,nny))
  allocate(w(nnx,nny))
  allocate(uRot(nnx,nny))
  allocate(vRot(nnx,nny))
  allocate(wRot(nnx,nny))

  allocate(UMD(nnx,nny))
  allocate(dStatU(nnx,nny))
  allocate(dStatUPrev(nnx,nny))
  allocate(dStatW(nnx,nny))
  allocate(dStatWPrev(nnx,nny))

  allocate(t_u_UMD(-1:1,nnx,nny))
  allocate(t_w_UMD(-1:1,nnx,nny))
!   allocate(t_u_UMD_ind(-1:1,500,nnx,nny)) !An upper limit of 500 flow structures at each point
!   allocate(t_w_UMD_ind(-1:1,500,nnx,nny)) !An upper limit of 500 flow structures at each point
  allocate(nUMD_U(-1:1,nnx,nny))
  allocate(nUMD_W(-1:1,nnx,nny))

  if(allocated(dStatU)) then
     write(*,*) 'Allocated dStatU'
  end if

  write(*,*) 'Usage: ExecutableName CutOffWaveNumber u_ud_cdf u_dd_cdf w_ud_cdf w_dd_cdf'
  call getarg(1,buffer)
  read(buffer,*) zLevel
  write(*,*) 'Working on zLevel = ', zLevel
  u_cutoff = 24
  v_cutoff = 24
  w_cutoff = 24

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

  !Initializing thresholds for downdraft, mean and updraft in horizontal and vertical velocity
  write(*,*) 'Ideally.. I should read the thresholds for the structures from the corresponding zLevel directory...based on the computation of the filtered PDF.'
  write(*,*) 'However..I am hard coding them for now. Until the filteredPDF code is formatted to write them in an appropriate format'

  uUD =   13.025
  uDD =   14.525
  wUD =   0.3099
  wDD =  -0.3500

  !Initialize all time scales 
  t_u_UMD = 0.0
  t_w_UMD = 0.0
!   t_u_UMD_ind = 0.0
!   t_w_UMD_ind = 0.0

  !Initialize all draft structures
  nUMD_U(:,:,:) = 1
  nUMD_W(:,:,:) = 1
  dStatU(:,:) = 0
  dStatUPrev(:,:) = 10 !Some value other that -1,0,1 to initialize
  dStatW(:,:) = 0
  dStatWPrev(:,:) = 10 !Some value other that -1,0,1 to initialize

  open(fileDT,file="dt",form="formatted")
  read(fileDT,'(A1)') dtFileReadLine
  read(fileDT,'(A1)') dtFileReadLine
  read(fileDT,'(A1)') dtFileReadLine
  read(fileDT,'(A1)') dtFileReadLine
  open(unit=fileXY, file="../../viz.abl.094999_102000.xy.data", form="unformatted", access="direct", recl=nvar*nnx*nny*sizeOfReal)
  do fileCounter=1,nt

     read(fileDT,'(e15.6,e15.6)') tLESnew, dt
     write(*,*) t, ' ', fileCounter, ' '
!     write(*,*) 'Reading record number ', (fileCounter-1)*nnz + zLevel
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

     !Count number of updraft, mean and downdraft structures
     dStatU = 0
     where(uRot .gt. uDD) 
        dStatU = -1  !Downdrafts or high speed streaks
     end where
     where(uRot .lt. uUD) 
        dStatU = +1  !Updrafts or low speed streaks
     end where

     do wty=1,nWtx
        do wtx=1,nWty
           write(writeFilename,"(I3,A10)") (wtx-1)*nWty + wty, '_UMD_U.txt'
           open(unit=23,file=writeFilename,access="append")
           wtYloc = (wty-1)*54+wtExtentY/2 
           wtXloc = (wtx-1)*54+1 
           write(23,*) tLESnew, ' ',dStatU(wtXloc,wtYloc), ' ', uRot(wtXloc, wtYloc)
           close(23)
        end do
     end do

     do i=1,nnx
        do j=1,nny
           t_u_UMD(dStatU(i,j),i,j) = t_u_UMD(dStatU(i,j),i,j) + dt
           if((dStatU(i,j)-dStatUPrev(i,j)) .ne. 0) then
              nUMD_U(dStatU(i,j),i,j) = nUMD_U(dStatU(i,j),i,j) + 1
           end if
!           t_u_UMD_ind(dStatU(i,j),nUMD_U(dStatU(i,j),i,j),i,j) = t_u_UMD_ind(dStatU(i,j),nUMD_U(dStatU(i,j),i,j),i,j) + dt
        end do
     end do
     dStatUPrev = dStatU
     
     dStatW = 0
     where(wRot .lt. wDD) 
        dStatW = -1 !Downdrafts
     end where
     where(wRot .gt. wUD) 
        dStatW = +1 !Updrafts
     end where

     do wty=1,nWtx
        do wtx=1,nWty
           write(writeFilename,"(I3,A10)") (wtx-1)*nWty + wty, '_UMD_W.txt'
           open(unit=23,file=writeFilename,access="append")
           wtYloc = (wty-1)*54+wtExtentY/2 
           wtXloc = (wtx-1)*54+1 
           write(23,*) tLESnew, ' ',dStatW(wtXloc,wtYloc), ' ', wRot(wtXloc,wtYloc)
           close(23)
        end do
     end do

     do i=1,nnx
        do j=1,nny
           t_w_UMD(dStatW(i,j),i,j) = t_w_UMD(dStatW(i,j),i,j) + dt
           if((dStatW(i,j)-dStatWPrev(i,j)) .ne. 0) then
              nUMD_W(dStatW(i,j),i,j) = nUMD_W(dStatW(i,j),i,j) + 1
           end if
!           t_w_UMD_ind(dStatW(i,j),nUMD_W(dStatW(i,j),i,j),i,j) = t_w_UMD_ind(dStatW(i,j),nUMD_W(dStatW(i,j),i,j),i,j) + dt
        end do
     end do
     dStatWPrev = dStatW

     tLESold = tLESnew

  end do

  write(writeFileName,'(I0.2)') zLevel
  call system("mkdir -p "//"zLevel"//trim(adjustl(writeFileName)))

  open(unit=23,file="zLevel"//trim(adjustl(writeFileName))//"uTimeScalesUMD")
  write(23,'(I, 8F)',advance="no") u_cutoff, uud, sum(t_u_UMD(-1,:,:)/real(nUMD_U(-1,:,:)))/dble(dble(nnx)*dble(nny)), &
       udd,  sum(t_u_UMD(1,:,:)/real(nUMD_U(1,:,:)))/dble(dble(nnx)*dble(nny)), &
       uMeanFromProfile, sum(t_u_UMD(0,:,:)/real(nUMD_U(0,:,:)))/dble(dble(nnx)*dble(nny))
  close(23)
  
  open(unit=23,file="zLevel"//trim(adjustl(writeFileName))//"wTimeScalesUMD",access="APPEND")
  write(23,'(I, 8F)',advance="no") w_cutoff, wud, sum(t_w_UMD(1,:,:)/real(nUMD_W(1,:,:)))/dble(dble(nnx)*dble(nny)), &
       wdd,  sum(t_w_UMD(-1,:,:)/real(nUMD_W(-1,:,:)))/dble(dble(nnx)*dble(nny)), &
       0.0, sum(t_w_UMD(0,:,:)/real(nUMD_W(0,:,:)))/dble(dble(nnx)*dble(nny))
  close(23)


!   !Initialize PDF variables
!   t_u_ud_min = 1.0
!   t_u_ud_max = 1000.0
!   t_u_dd_min = 1.0
!   t_u_dd_max = 1000.0
!   t_u_m_min = 1.0
!   t_u_m_max = 1000.0
!   t_w_ud_min = 1.0
!   t_w_ud_max = 1000.0
!   t_w_dd_min = 1.0
!   t_w_dd_max = 1000.0
!   t_w_m_min = 1.0
!   t_w_m_max = 1000.0
!   t_u_ud_binSize = 2.0
!   t_u_dd_binSize = 2.0
!   t_u_m_binSize = 2.0
!   t_w_ud_binSize = 2.0
!   t_w_dd_binSize = 2.0
!   t_w_m_binSize = 2.0
!   t_u_ud_nbins = int((t_u_ud_max-t_u_ud_min)/t_u_ud_binSize)+1
!   t_u_dd_nbins = int((t_u_dd_max-t_u_ud_min)/t_u_ud_binSize)+1
!   t_u_m_nbins = int((t_u_m_max-t_u_m_min)/t_u_m_binSize)+1
!   t_w_ud_nbins = int((t_w_ud_max-t_w_ud_min)/t_w_ud_binSize)+1
!   t_w_dd_nbins = int((t_w_dd_max-t_w_dd_min)/t_w_dd_binSize)+1
!   t_w_m_nbins = int((t_w_m_max-t_w_m_min)/t_w_m_binSize)+1
!   allocate(t_u_ud_pdf(0:t_u_ud_nBins))
!   allocate(t_u_dd_pdf(0:t_u_dd_nBins))
!   allocate(t_u_m_pdf(0:t_u_m_nBins))
!   allocate(t_w_ud_pdf(0:t_w_ud_nBins))
!   allocate(t_w_dd_pdf(0:t_w_dd_nBins))
!   allocate(t_w_m_pdf(0:t_w_m_nBins))
!   t_u_ud_pdf = 0.0
!   t_u_dd_pdf = 0.0
!   t_u_m_pdf = 0.0
!   t_w_ud_pdf = 0.0
!   t_w_dd_pdf = 0.0
!   t_w_m_pdf = 0.0
  

!   do it=1,500
!      do i=1,nnx
!         do j=1,nny
!            if(t_u_UMD_ind(it,+1,i,j) .gt. dt) then
!               binLoc = int(( t_u_UMD_ind(it,+1,i,j) - t_u_ud_min)/t_u_ud_binSize)
!               t_u_ud_pdf(binLoc) = t_u_ud_pdf(binLoc) + 1
!            end if
!            if(t_u_UMD_ind(it,-1,i,j) .gt. dt) then
!               binLoc = int(( t_u_UMD_ind(it,-1,i,j) - t_u_dd_min)/t_u_dd_binSize)
!               t_u_dd_pdf(binLoc) = t_u_dd_pdf(binLoc) + 1
!            end if
!            if(t_u_UMD_ind(it,0,i,j) .gt. dt) then
!               binLoc = int(( t_u_UMD_ind(it,0,i,j) - t_u_m_min)/t_u_m_binSize)
!               t_u_m_pdf(binLoc) = t_u_m_pdf(binLoc) + 1
!            end if

!            if(t_w_UMD_ind(it,+1,i,j) .gt. dt) then
!               binLoc = int(( t_w_UMD_ind(it,+1,i,j) - t_w_ud_min)/t_w_ud_binSize)
!               t_w_ud_pdf(binLoc) = t_w_ud_pdf(binLoc) + 1
!            end if
!            if(t_w_UMD_ind(it,-1,i,j) .gt. dt) then
!               binLoc = int(( t_w_UMD_ind(it,-1,i,j) - t_w_dd_min)/t_w_dd_binSize)
!               t_w_dd_pdf(binLoc) = t_w_dd_pdf(binLoc) + 1
!            end if
!            if(t_w_UMD_ind(it,0,i,j) .gt. dt) then
!               binLoc = int(( t_w_UMD_ind(it,0,i,j) - t_w_m_min)/t_w_m_binSize)
!               t_w_m_pdf(binLoc) = t_w_m_pdf(binLoc) + 1
!            end if
!         end do
!      end do
!   end do
!   t_u_ud_pdf = t_u_ud_pdf/(sum(t_u_ud_pdf)*t_u_ud_binSize)
!   t_u_dd_pdf = t_u_dd_pdf/(sum(t_u_dd_pdf)*t_u_dd_binSize)
!   t_u_m_pdf = t_u_m_pdf/(sum(t_u_m_pdf)*t_u_m_binSize)
!   t_w_ud_pdf = t_w_ud_pdf/(sum(t_w_ud_pdf)*t_w_ud_binSize)
!   t_w_dd_pdf = t_w_dd_pdf/(sum(t_w_dd_pdf)*t_w_dd_binSize)
!   t_w_m_pdf = t_w_m_pdf/(sum(t_w_m_pdf)*t_w_m_binSize)


!   write(*,*) 'Writing PDFs to file '
!   open(unit=23,file="zLevel"//trim(adjustl(writeFileName))//"t_u_ud_PDF")
!   write(23,*)'#  Bin      t_u_ud_PDF'
!   do i=0,t_u_ud_nBins
!      write(23,*) t_u_ud_min+0.5*t_u_ud_binSize + i*t_u_ud_binSize, t_u_ud_pdf(i)
!   end do
!   close(23)
!   open(unit=23,file="zLevel"//trim(adjustl(writeFileName))//"t_u_dd_PDF")
!   write(23,*)'#  Bin      t_u_dd_PDF'
!   do i=0,t_u_dd_nBins
!      write(23,*) t_u_dd_min+0.5*t_u_dd_binSize + i*t_u_dd_binSize, t_u_dd_pdf(i)
!   end do
!   close(23)
!   open(unit=23,file="zLevel"//trim(adjustl(writeFileName))//"t_u_m_PDF")
!   write(23,*)'#  Bin      t_u_m_PDF'
!   do i=0,t_u_m_nBins
!      write(23,*) t_u_m_min+0.5*t_u_m_binSize + i*t_u_m_binSize, t_u_m_pdf(i)
!   end do
!   close(23)
!   open(unit=23,file="zLevel"//trim(adjustl(writeFileName))//"t_w_ud_PDF")
!   write(23,*)'#  Bin      t_w_ud_PDF'
!   do i=0,t_w_ud_nBins
!      write(23,*) t_w_ud_min+0.5*t_w_ud_binSize + i*t_w_ud_binSize, t_w_ud_pdf(i)
!   end do
!   close(23)
!   open(unit=23,file="zLevel"//trim(adjustl(writeFileName))//"t_w_dd_PDF")
!   write(23,*)'#  Bin      t_w_dd_PDF'
!   do i=0,t_w_dd_nBins
!      write(23,*) t_w_dd_min+0.5*t_w_dd_binSize + i*t_w_dd_binSize, t_w_dd_pdf(i)
!   end do
!   close(23)
!   open(unit=23,file="zLevel"//trim(adjustl(writeFileName))//"t_w_m_PDF")
!   write(23,*)'#  Bin      t_w_m_PDF'
!   do i=0,t_w_m_nBins
!      write(23,*) t_w_m_min+0.5*t_w_m_binSize + i*t_w_m_binSize, t_w_m_pdf(i)
!   end do
!   close(23)

! !Time scales

!   write(*,*) 'Time scale using u'
!   write(*,*) 'Updrafts = ', sum(t_u_UMD(1,:,:)/real(nUMD_U(1,:,:)))*dt/dble(dble(nnx)*dble(nny))
!   write(*,*) 'Downdrafts = ', sum(t_u_UMD(0,:,:)/real(nUMD_U(0,:,:)))*dt/dble(dble(nnx)*dble(nny))
!   write(*,*) 'Mean = ', sum(t_u_UMD(-1,:,:)/real(nUMD_U(-1,:,:)))*dt/dble(dble(nnx)*dble(nny))

!   write(*,*) 'Number of updraft structures U = ', nUMD_U(1,1,1) 
!   write(*,*) 'Number of mean structures U = ', nUMD_U(0,1,1) 
!   write(*,*) 'Number of downdraft structures U = ', nUMD_U(-1,1,1) 
  
!   write(*,*) 'Number of updraft structures W = ', nUMD_W(1,1,1) 
!   write(*,*) 'Number of mean structures W = ', nUMD_W(0,1,1) 
!   write(*,*) 'Number of downdraft structures W = ', nUMD_W(-1,1,1) 

!   write(*,*) 'Time scale using w'
!   write(*,*) 'Updrafts = ', sum(t_w_UMD(1,:,:)/real(nUMD_W(1,:,:)))*dt/dble(dble(nnx)*dble(nny))
!   write(*,*) 'Downdrafts = ', sum(t_w_UMD(0,:,:)/real(nUMD_W(0,:,:)))*dt/dble(dble(nnx)*dble(nny))
!   write(*,*) 'Mean = ', sum(t_w_UMD(-1,:,:)/real(nUMD_W(-1,:,:)))*dt/dble(dble(nnx)*dble(nny))

  close(fileXY)


end program calcTimeScalesWriteUMD
