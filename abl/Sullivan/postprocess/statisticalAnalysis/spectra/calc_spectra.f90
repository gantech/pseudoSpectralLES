program calc_spectra
!This program reads the files that contain the 3d data (velocity, pressure, scalars) in the ".xydata" files and calculates the velocity spectra

implicit none
include "/opt/fftw/3.3.0.0/x86_64/include/fftw3.f"

integer :: i,j,k !Counter variables
integer :: iz
integer :: fileXY !File index for the xy.data file 
integer :: sizeOfReal=1 !The size of real in words
! GRID PARAMETERS
integer :: nnx, nny, nnz      ! x, y, and z grid dimensions 
integer :: nxy                ! number of grid points in a horizontal plane
real :: xl, yl, zl         ! x, y, and z domain sizes
real :: dx, dy, dz         ! x, y, and z grid lengths
integer :: nvar            ! Number of variables written out at each plane

real :: ugal, vgal

integer :: fcounter !Number of files used in the averaging process
real(kind=4), dimension(:,:,:), allocatable :: pA_xy !Plane Arrays                                                                                                     
real, allocatable, dimension(:) :: umean, vmean, wmean
real, allocatable, dimension(:) :: uvar, vvar, wvar
real, allocatable, dimension(:) :: l111, l112, l221, l222, l331, l332

!FFT variables
double precision, allocatable, dimension(:,:) :: fft_in !Temporary array to store the input to the FFTW function
double complex, allocatable, dimension(:,:) :: fft_out !Temporary array to store the output from the FFTW function
real, allocatable, dimension(:,:) :: uhat, vhat, what !Array to store the spectrum at pre-defined z levels
real, allocatable, dimension(:) :: waveN ! Wavenumber
real, allocatable, dimension(:) :: nCount !Number of points in each wave number ring.
integer :: nLevels = 50 !Number of levels at which the spectrum is to be measured.
!The levels at which the spectrum is to be measured
integer, dimension(50) :: zLevels = (/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50/) 
real :: deltaWaveN !The difference between two successive wave numbers
real :: pi = 3.14159265358979
integer :: r ! Distance of 2d wave number from orgin
integer*8 plan ; !Plan to create the FFT
integer :: stat


fileXY = 95
fcounter = 0
nnx = 768
nny = 768
nxy = nnx*nny
xl = 5120.0
yl = 5120.0
zl = 2048.0
nnz = 50
nvar = 5

allocate(pA_xy(nvar,nnx,nny))
allocate(umean(nnz))
allocate(vmean(nnz))
allocate(wmean(nnz))

write(*,*) 'nnx = ', nnx
write(*,*) 'nny = ', nny
write(*,*) 'nnz = ', nnz

dx = xl/nnx
dy = yl/nny
dz = zl/dble(256)

write(*,*) 'dx = ', dx
write(*,*) 'dy = ', dy
write(*,*) 'dz = ', dz

!Initialise FFT variables
!!Allocate complex arrays
allocate(waveN(nnx))
allocate(nCount(0:nnx/2))
allocate(fft_in(nnx,nny))
allocate(fft_out(nnx/2+1,nny))
allocate(uhat(0:nnx/2, nLevels))
allocate(vhat(0:nnx/2, nLevels))
allocate(what(0:nnx/2, nLevels))

allocate(uvar(nLevels))
allocate(vvar(nLevels))
allocate(wvar(nLevels))
allocate(l111(nLevels))
allocate(l112(nLevels))
allocate(l221(nLevels))
allocate(l222(nLevels))
allocate(l331(nLevels))
allocate(l332(nLevels))
l111 = 0
l112 = 0
l221 = 0
l222 = 0
l331 = 0
l332 = 0
uvar = 0
vvar = 0
wvar = 0

!!Wave number
deltaWaveN = 1/xl
do i =1,nnx/2+1
   waveN(i) = (i-1)
end do
do i=nnx/2+2,nnx
   waveN(i) = i-nnx-1
end do
nCount = 0
uhat = 0
vhat = 0
what = 0

!!Create FFTW plan
call dfftw_plan_dft_r2c_2d(plan, nnx, nny, fft_in, fft_out, FFTW_ESTIMATE)

open(fileXY,file="../../viz.abl.094999_102000.xy.data" ,form='unformatted',access='direct',recl=nvar*nnx*nny*sizeOfReal)
do fcounter = 1, 5000, 50
   do iz=1,nnz
      !Read data
      write(*,*) 'fcounter = ', fcounter, ' iz = ', iz
      read(fileXY,rec=(fcounter-1)*nnz+iz) pA_xy
      !Getting the means
      umean(iz) = sum(pA_xy(1,:,:))/real(nxy)
      vmean(iz) = sum(pA_xy(2,:,:))/real(nxy)
      wmean(iz) = sum(pA_xy(3,:,:))/real(nxy)
      
      !U spectra
      fft_in = pA_xy(1,:,:) - umean(zLevels(iz))
      call dfftw_execute_dft_r2c(plan, fft_in, fft_out)
      fft_out = fft_out/real(nxy)
      do i = 2,nnx/2
         do j = 1,nny
            r = nint(sqrt(waveN(i)*waveN(i) + waveN(j)*waveN(j)))
            if( r .le. nnx/2) then
               nCount(r) = nCount(r) + 2
               uhat(r,iz) = uhat(r,iz) + 2.0*fft_out(i,j)*conjg(fft_out(i,j))
            end if
         end do
      end do
      do i = 1,nnx/2+1,nnx/2
         do j = 1,nny
            r = nint(sqrt(waveN(i)*waveN(i) + waveN(j)*waveN(j)))
            if( r .le. nnx/2) then
               nCount(r) = nCount(r) + 1
               uhat(r,iz) = uhat(r,iz) + fft_out(i,j)*conjg(fft_out(i,j))
            end if
         end do
      end do
      
      l111(iz) = l111(iz) + sum(fft_out(1,:)*conjg(fft_out(1,:)))
      l112(iz) = l112(iz) + sum(fft_out(:,1)*conjg(fft_out(:,1)))
      uvar(iz) = uvar(iz) + sum(abs( fft_out(:,:) )**2 ) + sum(abs(fft_out(2:nnx/2,:))**2)
      
      !V spectra
      fft_in = pA_xy(2,:,:) - vmean(zLevels(iz))
      call dfftw_execute_dft_r2c(plan, fft_in, fft_out)
      fft_out = fft_out/real(nxy)
      do i = 2,nnx/2
         do j = 1,nny
            r = nint(sqrt(waveN(i)*waveN(i) + waveN(j)*waveN(j)))
            if( r .le. nnx/2) then
               nCount(r) = nCount(r) + 2
               vhat(r,iz) = vhat(r,iz) + 2.0*fft_out(i,j)*conjg(fft_out(i,j))
            end if
         end do
      end do
      do i = 1,nnx/2+1,nnx/2
         do j = 1,nny
            r = nint(sqrt(waveN(i)*waveN(i) + waveN(j)*waveN(j)))
            if( r .le. nnx/2) then
               nCount(r) = nCount(r) + 1
               vhat(r,iz) = vhat(r,iz) + fft_out(i,j)*conjg(fft_out(i,j))
            end if
         end do
      end do
      
      l221(iz) = l221(iz) + sum(fft_out(1,:)*conjg(fft_out(1,:)))
      l222(iz) = l222(iz) + sum(fft_out(:,1)*conjg(fft_out(:,1)))
      vvar(iz) = vvar(iz) + sum(abs( fft_out(:,:) )**2 ) + sum(abs(fft_out(2:nnx/2,:))**2)
      
      !W spectra
      fft_in = pA_xy(3,:,:) - wmean(zLevels(iz))
      call dfftw_execute_dft_r2c(plan, fft_in, fft_out)
      fft_out = fft_out/real(nxy)
      do i = 2,nnx/2
         do j = 1,nny
            r = nint(sqrt(waveN(i)*waveN(i) + waveN(j)*waveN(j)))
            if( r .le. nnx/2) then
               nCount(r) = nCount(r) + 2
               what(r,iz) = what(r,iz) + 2.0*fft_out(i,j)*conjg(fft_out(i,j))
            end if
         end do
      end do
      do i = 1,nnx/2+1,nnx/2
         do j = 1,nny
            r = nint(sqrt(waveN(i)*waveN(i) + waveN(j)*waveN(j)))
            if( r .le. nnx/2) then
               nCount(r) = nCount(r) + 1
               what(r,iz) = what(r,iz) + fft_out(i,j)*conjg(fft_out(i,j))
            end if
         end do
      end do
      
      l331(iz) = l331(iz) + sum(fft_out(1,:)*conjg(fft_out(1,:)))
      l332(iz) = l332(iz) + sum(fft_out(:,1)*conjg(fft_out(:,1)))
      wvar(iz) = wvar(iz) + sum(abs( fft_out(:,:) )**2 ) + sum(abs(fft_out(2:nnx/2,:))**2)
      
      !Verification of spectra by calculating variance in 3 different ways
!      write(*,*) 'Level ', zLevels(iz), sum((pA_xy(3,:,:) - wmean(zLevels(iz)))*(pA_xy(3,:,:) - wmean(zLevels(iz))))/real(nxy-1), sum(abs( fft_out(:,:) )**2 ) + sum(abs(fft_out(2:nnx/2,:))**2)
!           , sum(uhat(:,iz))
   end do
end do
close(75)

!Averaging over all the field files read so far
do i=1,nnx/2
   uhat(i,:) = uhat(i,:) * pi*((i+0.5)**2 - (i-0.5)**2) *3*nLevels*100/nCount(i)
   vhat(i,:) = vhat(i,:) * pi*((i+0.5)**2 - (i-0.5)**2) *3*nLevels*100/nCount(i)
   what(i,:) = what(i,:) * pi*((i+0.5)**2 - (i-0.5)**2) *3*nLevels*100/nCount(i)
end do

write(*,*) (nCount(j),j=0,nnx/2)
write(*,*) (waveN(j),j=1,nnx)
call dfftw_destroy_plan(plan)

open(unit=15, file="u_spectra.dat")
write(15,"(A, 50I)", advance="no") '#k ', (zLevels(iz),iz=1,nLevels)
write(15,*)
do j=0,nnx/2
   write(15,"(I,50F)", advance="no") j, (uhat(j,iz),iz=1,nLevels)
   write(15,*)
end do
close(15)


open(unit=15, file="v_spectra.dat")
write(15,"(A, 50I)", advance="no") '#k ', (zLevels(iz),iz=1,nLevels)
write(15,*)
do j=0,nnx/2
   write(15,"(I, 50F)", advance="no") j, (vhat(j,iz),iz=1,nLevels)
   write(15,*)
end do
close(15)


open(unit=15, file="w_spectra.dat")
write(15,"(A, 50I)", advance="no") '#k ', (zLevels(iz),iz=1,nLevels)
write(15,*)
do j=0,nnx/2
   write(15,"(I, 50F)", advance="no") j, (what(j,iz),iz=1,nLevels)
   write(15,*)
end do
close(15)

end program calc_spectra
