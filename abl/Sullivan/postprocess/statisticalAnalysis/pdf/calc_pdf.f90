program conditionalMeanVelProfiles
  
  !!!This is a program to compute the conditional mean vel profiles 
    !based on thresholding the velocity field. The thresholding values 
    !are decided by the pdf of the filtered velocity field with the 
    !filter applied at 3 times the peak in the w velocity spectrum at that level.
  
  !1. Find overall mean vel profile
  !2. Filter vel field at desired cutoff wavenumber, find pdf and fix threshold level at each z level.
  !3. Find conditional mean velocity based on thresholds
  
  implicit none
  include "/usr/global/fftw/3.2.2/pgi/include/fftw3.f"

  integer :: i,j,k,iz,fcounter !Counter variables

  ! GRID PARAMETERS
  integer :: nnx, nny, nnz      ! x, y, and z grid dimensions 
  integer :: nxy                ! number of grid points in a horizontal plane
  real :: xl, yl, zl         ! x, y, and z domain sizes
  real :: dx, dy, dz         ! x, y, and z grid lengths
  integer :: nscalars !Number of scalars
  logical :: qscal_sav, cscal_sav
  integer :: itsfc, iqsfc, icsfc
  real :: wtsfc, wqsfc, wcsfc
  real :: tsfcc, qsfcc, csfcc
  real :: ugtop, ugbot, vgtop, vgbot
  real :: dtdzf, divgls, fcor, amonin, utau
  real :: time_start, dt, z0
  real :: ugal, vgal
  integer :: stat

  real(kind=4), allocatable, dimension(:,:,:) :: u,v,w,t,e,p,q,c
  real(kind=4), allocatable, dimension(:,:) :: velTemp

  real :: umean, vmean, wmean
  real, dimension(:,:), allocatable :: velmean
  real, allocatable, dimension(:) :: umeanCondUD, vmeanCondUD, wmeanCondUD !Conditional means in updraft
  real, allocatable, dimension(:) :: umeanCondDD, vmeanCondDD, wmeanCondDD !Conditional means in downdraft
  real, allocatable, dimension(:) :: umeanCondLS, vmeanCondLS, wmeanCondLS !Conditional means in low speed streaks
  real, allocatable, dimension(:) :: umeanCondHS, vmeanCondHS, wmeanCondHS !Conditional means in high speed streaks
  
  
  !FFT variables
  double precision, allocatable, dimension(:,:) :: fft_in !Temporary array to store the input to the FFTW function
  double complex, allocatable, dimension(:,:) :: fft_out !Temporary array to store the output from the FFTW function
  real, allocatable, dimension(:) :: waveN ! Wavenumber
  real :: deltaWaveN !The difference between two successive wave numbers
  real :: pi = 3.14159265358979
  integer :: r ! Distance of 2d wave number from orgin
  integer*8 :: plan, plan_inv !Plan to create the FFT and it's inverse
  integer, dimension(:), allocatable :: u_cutoff, v_cutoff, w_cutoff !Filter cutoff wavenumber for each velocity component at each z level

  !Thresholding to define downdraft, mean and updraft in horizontal and vertical velocity
  real :: u_ud_cdf, u_dd_cdf, w_ud_cdf, w_dd_cdf !Thresholds based on the Cumulative Density Function of filtered data
  real, dimension(:), allocatable :: uDD, uUD !DD - DownDraft, UD - UpDraft
  real, dimension(:), allocatable :: wDD, wUD !DD - DownDraft, UD - UpDraft
  integer*8, dimension(:), allocatable :: uDDcount, uUDcount !Integer count of number of points in these structures
  integer*8, dimension(:), allocatable :: wDDcount, wUDcount !Integer count of number of points in these structures


  !PDF
  real :: umin, umax, uvar, wmin, wmax, wvar
  real, allocatable, dimension(:,:) :: updf, wpdf
  real :: binSizeU, binSizeW
  integer :: binLoc, nbinsU, nbinsW

  !Wind turbine locations
  integer :: wtXloc1, wtXloc2, wtYloc1, wtYloc2

  character(10) ::  fileName !Output file name
  character*20 :: buffer  !Buffer to read in command line arguments


  open(unit=75,file="fieldfiles")
  read(75,*,iostat=stat)buffer
  !call getarg(1,buffer)
  open(unit=50,file=buffer,form="unformatted")
  read(50) time_start, nnx, nny, nnz, xl, yl, zl, qscal_sav, cscal_sav
  read(50) dt, z0, itsfc, wtsfc, tsfcc, amonin, utau
  read(50) iqsfc, wqsfc, qsfcc, icsfc, wcsfc, csfcc
  read(50) dtdzf, divgls, fcor, ugtop, ugbot, vgtop, vgbot 
  close(50)
  ugal = (max(0.,max(ugtop,ugbot))+min(0.,min(ugtop,ugbot)))*0.5
  vgal = (max(0.,max(vgtop,vgbot))+min(0.,min(vgtop,vgbot)))*0.5
  
  dx = xl/nnx
  dy = yl/nny
  dz = zl/nnz
  nxy = nnx*nny
!  nnz = int((160./dz)+0.5) !Read and write data only at the levels below the top of the wind turbine.
  
  allocate(u(nnx,nny,nnz))
  allocate(v(nnx,nny,nnz))
  allocate(w(nnx,nny,nnz))
  allocate(t(nnx,nny,nnz))
  allocate(e(nnx,nny,nnz))
  allocate(p(nnx,nny,nnz))
  allocate(velTemp(nnx,nny))
  
  allocate(velmean(3,nnz))

  !Initialize PDF variables
  umin = -10.0
  umax = 22.0
  wmin = -5.0
  wmax = 5.0
  binSizeW = 0.02
  binSizeU = 0.1
  nBinsU = int((umax-umin)/binSizeU)+1
  nBinsW = int((wmax-wmin)/binSizeW)+1
  allocate(uPDF(0:nBinsU,nnz))
  allocate(wPDF(0:nBinsW,nnz))
  velmean = 0.0
  uPDF = 0.0
  wPDF = 0.0
  fcounter = 0
  read(75,*,iostat=stat)buffer
  do while (stat .eq. 0) 
     open(unit=50,file=buffer,form="unformatted")
     read(50) time_start, nnx, nny, nnz, xl, yl, zl, qscal_sav, cscal_sav
     read(50) dt, z0, itsfc, wtsfc, tsfcc, amonin, utau
     read(50) iqsfc, wqsfc, qsfcc, icsfc, wcsfc, csfcc
     read(50) dtdzf, divgls, fcor, ugtop, ugbot, vgtop, vgbot
     write(*,*) buffer
!     nnz = int((160./dz)+0.5) !Read and write data only at the levels below the top of the wind turbine.
     do iz=1,nnz
        read(50) u(:,:,iz), v(:,:,iz), w(:,:,iz), t(:,:,iz), e(:,:,iz), &
             p(1:nnx,1:nny,iz)
        u(:,:,iz) = u(:,:,iz) + ugal
        v(:,:,iz) = v(:,:,iz) + vgal
        umean = sum(u(:,:,iz))/real(nxy) 
        vmean = sum(v(:,:,iz))/real(nxy)
        wmean = sum(w(:,:,iz))/real(nxy)
        velmean(1,iz) = velmean(1,iz) + umean
        velmean(2,iz) = velmean(2,iz) + vmean        
        velmean(3,iz) = velmean(3,iz) + wmean
     end do
     do iz=1,nnz-1
        velTemp = 0.5*(u(:,:,iz) + u(:,:,iz+1))
        do i=1,nnx
           do j=1,nny
              binLoc = int((velTemp(i,j)-umin)/binSizeU)
              updf(binLoc,iz) = updf(binLoc,iz) + 1
              binLoc = int((w(i,j,iz)-wmin)/binSizeW)
              wpdf(binLoc,iz) = wpdf(binLoc,iz) + 1
           end do
        end do
     end do

     close(50)
     read(75,*,iostat=stat)buffer
     fcounter = fcounter + 1
  end do
  close(75)
  velmean = velmean/fcounter

  updf = updf/(real(fcounter*nxy))
  wpdf = wpdf/(real(fcounter*nxy))

  open(unit=23,file="uPDF")
  write(23,"(A,175I)",advance="no")'#  Bin ', (iz,iz=1,nnz-1)
  write(23,*) 
  do i=0,nbinsU
     write(23,"(F,175F)",advance="no") umin+0.5*binSizeU + i*binSizeU, (updf(i,iz)/binSizeU ,iz=1,nnz-1)
     write(23,*) 
  end do
  close(23)
  
  open(unit=23,file="wPDF")
  write(23,"(A,175I)",advance="no")'#  Bin ', (iz,iz=1,nnz-1)
  write(23,*) 
  do i=0,nbinsW
     write(23,"(F,175F)",advance="no") wmin+0.5*binSizeW + i*binSizeW, (wpdf(i,iz)/binSizeW ,iz=1,nnz-1)
     write(23,*) 
  end do
  close(23)

end program conditionalMeanVelProfiles
