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

program generateHDFfiles
  !!!This is a program to write the velocity field data into HDF5 files
  
  use HDF5
  implicit none

  integer :: i,j,k,it,wtx,wty,fileCounter,iz !Counter variables
  
  real :: xl, yl         ! x, y domain sizes
  real :: dx, dy         ! x, y grid lengths
  integer :: nnx, nny, nnz      ! x, y grid dimensions
  real(kind=8), dimension(:), allocatable :: xArr, yArr !Array of x and y locations on the grid
  integer :: nvar=5 !The number of variables stored in a plane at every time step
  integer :: nt !Number of time steps
  real :: t !Current time
  real :: dt  !Time step
  real(kind=4), dimension(:,:,:), allocatable :: pA_xy !Plane Arrays 
  real(kind=4), dimension(:,:), allocatable :: u,v,w !Plane Arrays 
  real(kind=4), dimension(:,:), allocatable :: uRot,vRot,wRot !Plane Arrays 
  real(kind=4), dimension(:,:,:), allocatable :: uWrite,vWrite,wWrite !Plane Arrays in the format to be written into vtr files
  real :: umean, vmean, wmean
  real, allocatable, dimension(:) :: umeanProfile, vmeanProfile, wmeanProfile !Profile of mean velocities from postprocessing of his.mp files
  integer :: zLevel !The z level at which the data is to be visualized.
  integer :: sizeOfReal=1 !The size of real in words on this system
  integer :: fileXY, fileDT  !The index of the files used to open the "xy.data" and the "dt" files
  integer :: fileUmeanProfile, fileUVarProfile !The index of the files used to open the "uxym" and "uvar" files
  real(kind=4) :: tLESnew, tLESold
  real(kind=4) :: dtLES
  integer :: itLES
  character :: dtFileReadLine

  !Variables to interpolate to yawed co-ordinate system
  real(kind=4) :: yawAngle
  real(kind=4) :: xLocCur, yLocCur, zLocCur
  integer :: ixl, ixu, iyl, iyu !The integer locations of the 4 points used to interpolate the velocity from the LES of ABL field to the reqd. point in the FAST field
  double precision :: wxlyl, wxuyl, wxlyu, wxuyu !The weights of the 4 points used to interpolate the velocity from the LES of ABL field to the reqd. point in the FAST field

  character(10) ::  fileName !Output file name
  character *20  buffer  !Buffer to read in command line arguments
  character *40  writeFilename, readFileName

  ! HDF5 write variables

  CHARACTER(LEN=4), PARAMETER :: Uxdsetname = "Ux"     ! Dataset name
  CHARACTER(LEN=4), PARAMETER :: Uydsetname = "Uy"     ! Dataset name
  CHARACTER(LEN=4), PARAMETER :: Uzdsetname = "Uz"     ! Dataset name

  INTEGER(HID_T) :: file_id       ! File identifier
  INTEGER(HID_T) :: dset_id       ! Dataset identifier
  INTEGER(HID_T) :: dspace_id     ! Dataspace identifier

  INTEGER(HSIZE_T), DIMENSION(3) :: dims = (/1,1,1/) ! Dataset dimensions
  INTEGER     ::   rank = 3                        ! Dataset rank

  INTEGER, DIMENSION(10,11,12) :: dset_data, data_out ! Data buffers
  INTEGER(HSIZE_T), DIMENSION(3) :: data_dims

  integer(kind = 4) :: iErr = 0 !To store the output of each HDF5 command

  double precision :: temp !Temporary variable
  

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
  fileUMeanProfile = 83
  fileUVarProfile = 84
  
!  nt = 10
  nt =  400

  !
  ! Initialize HDF FORTRAN interface.
  !
  CALL h5open_f(iErr)
  data_dims(1) = nnx
  data_dims(2) = nny
  data_dims(3) = nnz

  allocate(pA_xy(nvar,nnx,nny))
  allocate(u(nnx,nny))
  allocate(v(nnx,nny))
  allocate(w(nnx,nny))
  allocate(uRot(nnx,nny))
  allocate(vRot(nnx,nny))
  allocate(wRot(nnx,nny))
  allocate(uWrite(nnx,nny,nnz))
  allocate(vWrite(nnx,nny,nnz))
  allocate(wWrite(nnx,nny,nnz))  
  allocate(umeanProfile(nnz))
  allocate(vmeanProfile(nnz))
  allocate(wmeanProfile(nnz))

  write(*,*) 'Writing unfiltered fields'
  
  !Read in umean and uvar Profiles
  open(fileUMeanProfile,file="uxym")
  read(fileUMeanProfile,'(A1)') dtFileReadLine !Dummy don't bother
  do iz = 1, nnz
     read(fileUMeanProfile,*) umeanProfile(iz), vmeanProfile(iz), wmeanProfile(iz)
     temp = umeanProfile(iz)*cos(yawAngle) + vmeanProfile(iz)*sin(yawAngle)
     vmeanProfile(iz) = -umeanProfile(iz)*sin(yawAngle) + vmeanProfile(iz)*cos(yawAngle)
     umeanProfile(iz) = temp
  end do
  close(fileUMeanProfile)

  open(unit=fileXY, file="../viz.abl.094999_102000.xy.data", form="unformatted", access="direct", recl=nvar*nnx*nny*sizeOfReal)
  open(fileDT,file="dt",form="formatted")
  read(fileDT,'(A1)') dtFileReadLine
  read(fileDT,'(A1)') dtFileReadLine
  read(fileDT,'(A1)') dtFileReadLine
  read(fileDT,'(A1)') dtFileReadLine
  do fileCounter=1,nt
     read(fileDT,'(e15.6,e15.6)') tLESnew, dtLES
     write(*,*) 'Writing time ', fileCounter, ' ', tLESnew
     do zLevel = 1,nnz
!        write(*,*) 'Reading record number ', (fileCounter-1)*nnz + zLevel
        read(fileXY,rec=(fileCounter-1)*nnz + zLevel) pA_xy
        u = pA_xy(1,:,:) + 7.5
        v = pA_xy(2,:,:)
        w = pA_xy(3,:,:)
        if(zLevel /= 1) then
           read(fileXY,rec=(fileCounter-1)*nnz + zLevel-1) pA_xy
           w = 0.5 * ( w + pA_xy(3,:,:))
        else
           w = 0.5 * ( w + 0.0)
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
        !Now rotate the velocity vector using 'u' as a temporary vector
        u = uRot*cos(yawAngle) + vRot*sin(yawAngle)
        vRot = -uRot*sin(yawAngle) + vRot*cos(yawAngle)
        uRot = u
   
        uWrite(:,:,zLevel) = uRot(:,:)
        vWrite(:,:,zLevel) = vRot(:,:)
        wWrite(:,:,zLevel) = wRot(:,:)
        
     end do

     
     !
     ! Create a new file using default properties.
     !
     write(writeFilename,'(F07.1)') tLESnew
     writeFilename = 'velData_'//trim(adjustl(writeFilename))//'.h5'
     CALL h5fcreate_f(writeFilename, H5F_ACC_TRUNC_F, file_id, iErr)

     !
     ! Create the dataspace.
     !
     CALL h5screate_simple_f(rank, data_dims, dspace_id, iErr)
     !
     ! Create the dataset with default properties.
     !
     CALL h5dcreate_f(file_id, Uxdsetname, H5T_NATIVE_REAL, dspace_id, &
          dset_id, iErr)
     !
     ! Write the dataset.
     !
     CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, uWrite, data_dims, iErr)
     !
     ! Close the dataset.
     !
     CALL h5dclose_f(dset_id, iErr)
     
     !
     ! Create the dataspace.
     !
     CALL h5screate_simple_f(rank, data_dims, dspace_id, iErr)
     !
     ! Create the dataset with default properties.
     !
     CALL h5dcreate_f(file_id, Uydsetname, H5T_NATIVE_REAL, dspace_id, &
          dset_id, iErr)
     !
     ! Write the dataset.
     !
     CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, vWrite, data_dims, iErr)
     !
     ! Close the dataset.
     !
     CALL h5dclose_f(dset_id, iErr)

     !
     ! Create the dataspace.
     !
     CALL h5screate_simple_f(rank, data_dims, dspace_id, iErr)
     !
     ! Create the dataset with default properties.
     !
     CALL h5dcreate_f(file_id, Uzdsetname, H5T_NATIVE_REAL, dspace_id, &
          dset_id, iErr)
     !
     ! Write the dataset.
     !
     CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, wWrite, data_dims, iErr)
     !
     ! Close the dataset.
     !
     CALL h5dclose_f(dset_id, iErr)


  end do

  close(fileXY)

  !
  ! Close the HDF5 file.
  !
  CALL h5fclose_f(file_id, iErr)
  !
  ! Close HDF5 FORTRAN interface.
  !
  CALL h5close_f(iErr)
  
end program generateHDFfiles
