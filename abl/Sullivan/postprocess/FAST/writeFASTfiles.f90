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


program	writeFASTFiles
  !!!Write input files for FAST computations with ABL inflow fields
  !  use LIB_VTK_IO
  integer nnx, nny, nnz !Dimenions of the domain
  integer nvar_o !Number of variables written out at each plane and time step
  real(kind=4), dimension(:,:,:), allocatable, target :: pA_xyOld, pA_xyNew !Plane Arrays 
  logical fileReadNow !A bit used to determine whether data is to be read in at the current time step or not 
  integer fileXY  !Index of files used in the program
  integer fileDT  !Index of file "dt"
  character dtFileReadLine !To skip first few lines in the dt file
  character *120 fileXYname !The name of the file to opened for reading
  integer sizeOfReal !Size of real in bytes = 1
  integer iTstart, ntFAST, ntLES !The starting and the total number of time steps.
  integer iReadVizFile !Flag to indicate whether the viz file has to be read or not
  character*6 strItStart !The starting time step as a string through command line argument
  integer(kind=4) itLES, itFAST, iVar !Iterator variables
  integer i,j,k !More iterator variables
  real(kind=4) tLESnew, tLESold, tFAST, tStart !Current time of LES and FAST files that are open
  real(kind=4) dtFAST !Time step for FAST simulation
  real(kind=4) :: xl, yl, zl
  real(kind=4) :: dx, dy, dz
  real(kind=4) :: yawAngle
  
  real(kind=4), allocatable, dimension(:) :: xLoc, yLoc, zLoc

  !Variables needed to write into Turbsim - .bts format
  integer fileWTbts !The file index for the bts file currently being written into
  character *120 writeFileBTSname !The name of the bts file to be written into
  real(kind=4), dimension(3) :: vslope, vintercept
  integer(kind=2)  ::id=7
  real(kind=4) :: uHub, hubHt, zBottom
  integer(kind=4) :: nChar
  integer(kind=1) :: ci
  integer :: wtx, wty !Counter for the wind turbine number
  integer :: nWtx, nWty !Number of wind turbines in each direction
  integer :: wtxStart, wtxEnd !Starting and ending turbine row in x direction
  integer(kind=4) :: wtExtentY, wtExtentZ  
  !wtExtentY - Number of the points in the horizontal direction across 1 wind turbine.
  !wtExtentZ - Number of the points in the vertical direction across 1 wind turbine.
  real(kind=4) :: dt 
  real(kind=4) :: xLocCur, yLocCur, zLocCur
  integer :: ixl, ixu, iyl, iyu !The integer locations of the 4 points used to interpolate the velocity from the LES of ABL field to the reqd. point in the FAST field
  double precision :: wxlyl, wxuyl, wxlyu, wxuyu !The weights of the 4 points used to interpolate the velocity from the LES of ABL field to the reqd. point in the FAST field
  double precision :: uVel, vVel, wVel !The interpolated value of velocity at the reqd. point in the FAST field
  double precision :: uVelNew, vVelNew, wVelNew !The interpolated value of velocity at the reqd. point in the FAST field
  double precision :: uVelOld, vVelOld, wVelOld !The interpolated value of velocity at the reqd. point in the FAST field
  double precision :: tmp !Temporary variable

  character *20  buffer  !Buffer to read in command line arguments

  nnx = 768
  nny = 768
  nnz = 256
  xl = 5120.0
  yl = 5120.0
  zl = 2048.0
  yawAngle = 21.5*3.14159265359/180.0 
  dx = xl/real(nnx)
  dy = yl/real(nny)
  dz = zl/real(nnz)
  allocate(xLoc(nnx))
  allocate(yLoc(nny))
  allocate(zLoc(nnz))
  do i = 1,nnx
     xLoc(i) = 0.5*dx + real(i-1)*dx
  end do
  do i = 1,nny
     yLoc(i) = 0.5*dy + real(i-1)*dy
  end do
  do i = 1,nnz
     zLoc(i) = 0.5*dz + real(i-1)*dz
  end do
  nvar_o = 5
  fileXY = 95
  fileDT = 93
  sizeOfReal = 1
  write(*,*) 'Taking sizeOfReal = 1. Specified in word counts rather than bytes'
!  write(*,*) 'Taking sizeOfReal = 4. Specified in bytes rather than word counts'
  nTFAST = 7500  
  nt = 7500
  nTLES = 7000

  call getarg(1,buffer)
  read(buffer,*) wtxStart

  call getarg(2,buffer)
  read(buffer,*) wtxEnd
  write(*,*) 'Writing turbines from row', wtxStart, ' to ', wtxEnd
  
!  nWtx = 15
!  nWtx = 1
  nWty = 14
!  nWty = 1
  dt = 0.2
  
  !Turbsim - .bts files format
  vslope(1) = 100.0
  vslope(2) = 100.0
  vslope(3) = 100.0
  vintercept(1) = 2000.0
  vintercept(2) = 2000.0
  vintercept(3) = 2000.0
  hubHt = 90.0
  uHub = 14.0
  zBottom = 0.5*dz
  wtExtentY = 25
  wtExtentZ = 25

  do wtx=wtxStart, wtxEnd
    do wty=1,nWty
       write(writeFileBTSname,"(I3.3)") (wtx-1)*nWty + wty
       writeFileBTSname = trim(writeFileBTSname)//'.bts'
!       write(*,*) writeFilename
       open(unit=40,file=writeFileBTSname,form="binary",access="sequential", RECL=2)
       write(40) id
       write(40) wtExtentZ
       write(40) wtExtentY
       write(40) 0
       write(40) nt
       write(40) dz
       write(40) dy
       write(40) dt
       write(40) uHub
       write(40) hubHt
       write(40) zBottom
       do i=1,3
          write(40) vslope(i)
          write(40) vintercept(i)
       end do
       write(40) int(17,kind=4)
       writeFileBTSname = "LES-MCBL 20 mins "
       do i = 1,17
          write(40) int(ICHAR(writeFileBTSname(i:i)),kind=1)
       end do
       close(40)
    end do
 end do
 

 allocate(pA_xyNew(nvar_o,nnx,nny)) !Allocate array
 allocate(pA_xyOld(nvar_o,nnx,nny)) !Allocate array
 
 open(fileXY,file="../viz.abl.094999_102000.xy.data",form='unformatted',access='direct',recl=nvar_o*nnx*nny*sizeOfReal)

 !Open file containing dt info about LES of ABL run
 open(fileDT,file="dt",form="formatted")
 read(fileDT,'(A1)') dtFileReadLine
 read(fileDT,'(A1)') dtFileReadLine
 read(fileDT,'(A1)') dtFileReadLine
 read(fileDT,'(A1)') dtFileReadLine
 read(fileDT,'(e15.6,e15.6,I11)') tLESnew, dtLES

 tStart = 38219.4
 itLES = 0
 itFAST = 0
 tLESnew = tStart
 tLESold = tStart
 tFAST = tStart
 dtFAST = 0.2

 !Write the starting time step data first. No interpolation required
 do k = 1,wtExtentZ
    read(fileXY,rec=k) pA_xyNew
    do wtx=wtxStart, wtxEnd
       i = (wtx-1)*54+1
       xLocCur = (i-1)*dx
       do wty=1,nWty
          write(writeFileBTSname,"(I3.3)") (wtx-1)*nWty + wty
          writeFileBTSname = trim(writeFileBTSname)//'.bts'
          open(unit=40,file=writeFileBTSname,form="binary",access="append", RECL=2)
          do j = (wty-1)*54+1, (wty-1)*54+wtExtenty
             yLocCur = (j-1)*dy
             call nearestPoints(dble(xLocCur*cos(yawAngle)-yLocCur*sin(yawAngle)), dble(xLocCur*sin(yawAngle)+yLocCur*cos(yawAngle)), ixl, ixu, iyl, iyu, wxlyl, wxuyl, wxlyu, wxuyu)
             uVel = pA_xyNew(1,ixl,iyl)*wxlyl + pA_xyNew(1,ixu,iyl)*wxuyl + pA_xyNew(1,ixl,iyu)*wxlyu + pA_xyNew(1,ixu,iyu)*wxuyu + 7.5
             vVel = pA_xyNew(2,ixl,iyl)*wxlyl + pA_xyNew(2,ixu,iyl)*wxuyl + pA_xyNew(2,ixl,iyu)*wxlyu + pA_xyNew(2,ixu,iyu)*wxuyu
             wVel = pA_xyNew(3,ixl,iyl)*wxlyl + pA_xyNew(3,ixu,iyl)*wxuyl + pA_xyNew(3,ixl,iyu)*wxlyu + pA_xyNew(3,ixu,iyu)*wxuyu
             tmp = uVel*cos(yawAngle) + vVel*sin(yawAngle)
             vVel = -uVel*sin(yawAngle) + vVel*cos(yawAngle)
             uVel = tmp
             write(40) INT(  uVel*vslope(1) + vintercept(1), kind=2)
             write(40) INT(  vVel*vslope(2) + vintercept(3), kind=2)
             write(40) INT(  wVel*vslope(3) + vintercept(3), kind=2)
          end do
       end do
       close(40)
    end do
 end do
 
 tFAST = tStart
 tLES = tStart
 fileReadNow = .false.
 do itFAST = 1, ntFAST     
    write(*,*) 'iTFAST = ', iTFAST, ' itLES = ', itLES
    write(*,*) 'tFAST = ', tFAST, ' tLESnew = ', tLESnew
    tFAST = tFAST + dtFAST
    if(tFAST > tLESnew)  then
       fileReadNow = .true.
       tLESold = tLESnew
       write(*,*) 'Read LES time information now'
       read(fileDT,'(e15.6,e15.6)') tLESnew, dtLES
       itLES = itLES + 1
       write(*,*) 'tLESnew = ', tLESnew, ' tLESold = ', tLESold, ' itLES = ', itLES
    end if
    do k=1,wtExtentZ
       read(fileXY,rec=iTLES*50+k) pA_xyNew
       read(fileXY,rec=(iTLES-1)*50+k) pA_xyOld
       do wtx=wtxStart, wtxEnd
          i = (wtx-1)*54+1
          xLocCur = (i-1)*dx
          do wty=1,nWty
             write(writeFileBTSname,"(I3.3)") (wtx-1)*nWty + wty
             writeFileBTSname = trim(writeFileBTSname)//'.bts'
             open(unit=40,file=writeFileBTSname,form="binary",access="append", RECL=2)
             do j = (wty-1)*54+1, (wty-1)*54+wtExtenty
                yLocCur = (j-1)*dy
                call nearestPoints(dble(xLocCur*cos(yawAngle)-yLocCur*sin(yawAngle)), dble(xLocCur*sin(yawAngle)+yLocCur*cos(yawAngle)), ixl, ixu, iyl, iyu, wxlyl, wxuyl, wxlyu, wxuyu)
                uVelNew = pA_xyNew(1,ixl,iyl)*wxlyl + pA_xyNew(1,ixu,iyl)*wxuyl + pA_xyNew(1,ixl,iyu)*wxlyu + pA_xyNew(1,ixu,iyu)*wxuyu + 7.5
                vVelNew = pA_xyNew(2,ixl,iyl)*wxlyl + pA_xyNew(2,ixu,iyl)*wxuyl + pA_xyNew(2,ixl,iyu)*wxlyu + pA_xyNew(2,ixu,iyu)*wxuyu
                wVelNew = pA_xyNew(3,ixl,iyl)*wxlyl + pA_xyNew(3,ixu,iyl)*wxuyl + pA_xyNew(3,ixl,iyu)*wxlyu + pA_xyNew(3,ixu,iyu)*wxuyu
                tmp = uVelNew*cos(yawAngle) + vVelNew*sin(yawAngle)
                vVelNew = -uVelNew*sin(yawAngle) + vVelNew*cos(yawAngle)
                uVelNew = tmp

                uVelOld = pA_xyOld(1,ixl,iyl)*wxlyl + pA_xyOld(1,ixu,iyl)*wxuyl + pA_xyOld(1,ixl,iyu)*wxlyu + pA_xyOld(1,ixu,iyu)*wxuyu + 7.5
                vVelOld = pA_xyOld(2,ixl,iyl)*wxlyl + pA_xyOld(2,ixu,iyl)*wxuyl + pA_xyOld(2,ixl,iyu)*wxlyu + pA_xyOld(2,ixu,iyu)*wxuyu
                wVelOld = pA_xyOld(3,ixl,iyl)*wxlyl + pA_xyOld(3,ixu,iyl)*wxuyl + pA_xyOld(3,ixl,iyu)*wxlyu + pA_xyOld(3,ixu,iyu)*wxuyu
                tmp = uVelOld*cos(yawAngle) + vVelOld*sin(yawAngle)
                vVelOld = -uVelOld*sin(yawAngle) + vVelOld*cos(yawAngle)
                uVelOld = tmp

!                 write(*,*) uVelNew, vVelNew, wVelNew
!                 write(40) INT(  (uVelNew + 7.5)*vslope(1) + vintercept(1), kind=2)
!                 write(40) INT(  vVelNew*vslope(2) + vintercept(3), kind=2)
!                 write(40) INT(  wVelNew*vslope(3) + vintercept(3), kind=2)
!                write(*,*) (tLESnew-tFAST), uVelOld, (tFAST-tLESold), uVelNew, (tLESnew - tLESold)
                write(40) INT( ( ((tLESnew-tFAST) * uVelOld + (tFAST-tLESold) * uVelNew)/(tLESnew - tLESold) )  * vslope(1) + vintercept(1), kind=2)
                write(40) INT( ( ((tLESnew-tFAST) * vVelOld + (tFAST-tLESold) * vVelNew)/(tLESnew - tLESold) )  * vslope(2) + vintercept(2), kind=2)
                write(40) INT( ( ((tLESnew-tFAST) * wVelOld + (tFAST-tLESold) * wVelNew)/(tLESnew - tLESold) )  * vslope(3) + vintercept(3), kind=2)
             end do
          end do
          close(40)
       end do
    end do
    fileReadNow = .false.
 end do
 
 close(fileXY)
 
end program writeFASTFiles
