      program rddns
c
c ---------------- read les history file from gabls runs
c
      parameter(nq=110000)
      real t(nq), dt(nq),utau(nq),
     +     ziavg(nq), amonin(nq), holtop(nq),
     +     tsfcc(nq), uwsfc(nq), vwsfc(nq),
     +     divgmax(nq), wt_min(nq), wt_le(nq),
     +     ucfl(nq), vcfl(nq), wcfl(nq), wtsfc(nq),
     +     usqr(nq), vsqr(nq), wsqr(nq),tsqr(nq),
     +     uwle(nq), uwsb(nq), uwtot(nq),
     +     vwle(nq), vwsb(nq), vwtot(nq),
     +     wtle(nq), wtsb(nq), wttot(nq),
     +     englez(nq), eavg(nq), wabs(nq)
c
      character path*120
      character done*4, dum*1, xyfile*16
c
      done='nope'
      knt = 0
      lu = 99
      luxy = 98
      xyfile='abl.his.xy'
      ibreak = 0
    1 continue
      call rdfile(lu,path,done,ibreak)
      if(done .eq. 'done') go to 1000
      write(6,6200) path
 6200 format(' working on file = ',a120)
c
    2 continue
c
      if(ibreak .eq. 1) then
         ibreak = 0
         read(lu,10,end=999) dum
         read(lu,10,end=999) dum
   10    format(a1)
         go to 2
      endif
      knt = knt + 1
      read(lu,6000,end=999) t(knt), dt(knt),utau(knt),
     +             ziavg(knt), amonin(knt), holtop(knt),
     +             tsfcc(knt), uwsfc(knt), vwsfc(knt),
     +             divgmax(knt), wt_min(knt), wt_le(knt),
     +             ucfl(knt), vcfl(knt), wcfl(knt), wtsfc(knt),
     +             usqr(knt), vsqr(knt), wsqr(knt),tsqr(knt),
     +             uwle(knt), uwsb(knt), uwtot(knt),
     +             vwle(knt), vwsb(knt), vwtot(knt),
     +             wtle(knt), wtsb(knt), wttot(knt),
     +             englez(knt), eavg(knt), wabs(knt)
 6000 format(5e17.8)
      go to 2
c
  999 continue
c
c ------------------- at end of file
c
      knt = knt - 1
      go to 1
c
c --------------- at end of all files
c
 1000 continue
      knt = knt - 1
      write(6,6100) knt
 6100 format(' Total number of all points = ',i6)
c
c -------------------- write time series now
c
      open(luxy,file='timeHistoryFiles/dt',form='formatted')
      write(luxy,2000) (t(i), dt(i),i=1,knt)
 2000 format('#k "dt variation"',/,
     +       '#lx "t "',/,
     +       '#ly "\gdt"',/,
     +       '#k "DT   "',/,
     +       (2e15.6))
      close(luxy)
c
      open(luxy,file='timeHistoryFiles/tsfcc',form='formatted')
      write(luxy,2001) (t(i), tsfcc(i),i=1,knt)
 2001 format('#n',/,
     +       '#k "tsfcc variation"',/,
     +       '#lx "t "',/,
     +       '#ly "\gh\sbsfc\eb"',/,
     +       '#k "TSFCC"',/,
     +       (2e15.6))
      close(luxy)
c
      open(luxy,file='timeHistoryFiles/utau',form='formatted')
      write(luxy,1001) (t(i), utau(i),i=1,knt)
 1001 format('#n',/,
     +       '#k "utau variation"',/,
     +       '#lx "t "',/,
     +       '#ly "u\sb\gt\eb"',/,
     +       '#k "UTAU"',/,
     +       (2e15.6))
      close(luxy)
c
      open(luxy,file='timeHistoryFiles/wtsfc',form='formatted')
      write(luxy,1002) (t(i), wtsfc(i),i=1,knt)
 1002 format('#n',/,
     +       '#k "q_* variation"',/,
     +       '#lx "t "',/,
     +       '#ly "q\sb*\eb"',/,
     +       '#k "QSTAR"',/,
     +       (2e15.6))
      close(luxy)
c
      open(luxy,file='timeHistoryFiles/ziavg',form='formatted')
      write(luxy,1003) (t(i), ziavg(i),i=1,knt)
 1003 format('#n',/,
     +       '#k "zi variation"',/,
     +       '#lx "t "',/,
     +       '#ly "zi"',/,
     +       '#k "ZI"',/,
     +       (2e15.6))
      close(luxy)
c
      open(luxy,file='timeHistoryFiles/amonin',form='formatted')
      write(luxy,1004) (t(i), amonin(i),i=1,knt)
 1004 format('#n',/,
     +       '#k "L variation"',/,
     +       '#lx "t "',/,
     +       '#ly "L"',/,
     +       '#k "L-MO"',/,
     +       (2e15.6))
      close(luxy)
c
      open(luxy,file='timeHistoryFiles/holtop',form='formatted')
      write(luxy,1005) (t(i), holtop(i),i=1,knt)
 1005 format('#n',/,
     +       '#k "zi/L variation"',/,
     +       '#lx "t "',/,
     +       '#ly "zi/L"',/,
     +       '#k "zi/L"',/,
     +       (2e15.6))
      close(luxy)
c
      open(luxy,file='timeHistoryFiles/divgmax',form='formatted')
      write(luxy,2002) (t(i), divgmax(i),i=1,knt)
 2002 format('#n',/,
     +       '#k "divergence variation"',/,
     +       '#lx "t "',/,
     +       '#y l',/,
     +       '#ly "\dl\sp2\epU"',/,
     +       '#k "DIVG "',/,
     +       (2e15.6))
      close(luxy)
c
      open(luxy,file='timeHistoryFiles/uw',form='formatted')
      write(luxy,2003) (t(i), uwle(i), uwsb(i), uwtot(i),i=1,knt)
 2003 format('#n',/,
     +       '#k "uw flux variation, first resolved part"',/,
     +       '#lx "t "',/,
     +       '#ly "uw"',/,
     +       '#k "UWLE "',/,
     +       (4e15.6))
      close(luxy)
      
c
      open(luxy,file='timeHistoryFiles/vw',form='formatted')
      write(luxy,3001) (t(i), vwle(i), vwsb(i), vwtot(i), i=1,knt)
 3001 format('#n',/,
     +       '#k "vw flux variation, first resolved part"',/,
     +       '#lx "t "',/,
     +       '#ly "vw"',/,
     +       '#k "VWLE "',/,
     +       (4e15.6))
      close(luxy)
c
      open(luxy,file='timeHistoryFiles/wt',form='formatted')
      write(luxy,4001) (t(i), wtle(i), wtsb(i), wttot(i), i=1,knt)
 4001 format('#n',/,
     +       '#k "wt flux variation, first resolved part"',/,
     +       '#lx "t "',/,
     +       '#ly "w\gh"',/,
     +       '#k "WTLE "',/,
     +       (4e15.6))
      close(luxy)
c
      open(luxy,file='timeHistoryFiles/uvar',form='formatted')
      write(luxy,5001) (t(i), usqr(i), vsqr(i), wsqr(i), i=1,knt)
 5001 format('#n',/,
     +       '#k "u variances"',/,
     +       '#lx "t "',/,
     +       '#ly "u\sp2\ep, v\sp2\ep, w\sp2\ep"',/,
     +       '#k "UPS  "',/,
     +       (4e15.6))
      close(luxy)

c
      open(luxy,file='timeHistoryFiles/cfl',form='formatted')
      write(luxy,6001) (t(i), ucfl(i), vcfl(i), wcfl(i), i=1,knt)
 6001 format('#n',/,
     +       '#k "ucfl"',/,
     +       '#lx "t "',/,
     +       '#ly "umx/dx, vmx/dy, wmx/dz"',/,
     +       '#k "UCFL "',/,
     +       (4e15.6))
      close(luxy)
c
      open(luxy,file='timeHistoryFiles/uwsfc',form='formatted')
      write(luxy,7001) (t(i), uwsfc(i),i=1,knt)
 7001 format('#n',/,
     +       '#k "uwsfc"',/,
     +       '#lx "t "',/,
     +       '#ly "uw\sbsfc\eb, vw\sbsfc\eb"',/,
     +       '#k "UWSFC"',/,
     +       (2e15.6))
      close(luxy)

      open(luxy,file='timeHistoryFiles/vwsfc',form='formatted')
      write(luxy,7002) (t(i), vwsfc(i),i=1,knt)
 7002 format('#k "vwsfc"',/,
     +       '#m 2',/,
     +       '#k "VWSFC"',/,
     +       (2e15.6))
      close(luxy)

c
      open(luxy,file='timeHistoryFiles/wt_min',form='formatted')
      write(luxy,7003) (t(i), wt_min(i),i=1,knt)
 7003 format('#n',/,
     +       '#k "wt_min"',/,
     +       '#lx "t "',/,
     +       '#ly "w\gh"',/,
     +       '#m 2',/,
     +       '#k "WTMIN"',/,
     +       (2e15.6))
      close(luxy)

      open(luxy,file='timeHistoryFiles/wtle',form='formatted')
      write(luxy,7004) (t(i), wt_le(i),i=1,knt)
 7004 format('#k "wt_le"',/,
     +       '#m 1',/,
     +       '#k "WTLEM"',/,
     +       (2e15.6))
      close(luxy)
c
      open(luxy,file='timeHistoryFiles/englez',form='formatted')
      write(luxy,7005) (t(i), englez(i),i=1,knt)
 7005 format('#n',/,
     +       '#k "englez"',/,
     +       '#lx "t "',/,
     +       '#ly "TKE"',/,
     +       '#k "RESTKE"',/,
     +       (2e15.6))
      close(luxy)

      open(luxy,file='timeHistoryFiles/eavg',form='formatted')
      write(luxy,7006) (t(i), eavg(i),i=1,knt)
 7006 format('#k "eavg"',/,
     +       '#m 2',/,
     +       '#k "TKESGS"',/,
     +       (2e15.6))
      close(luxy)
c
      open(luxy,file='timeHistoryFiles/wabs',form='formatted')
      write(luxy,7008) (t(i), wabs(i),i=1,knt)
 7008 format('#n',/,
     +       '#k "max w"',/,
     +       '#lx "t "',/,
     +       '#ly "|w|"',/,
     +       '#k "WABS"',/,
     +       (2e15.6))
      close(luxy)
c
      stop
      end
      subroutine rdfile(lu,path,done,ibreak)
c
      character path*120, done*4
      data kount /0/
      data ilast,inc, icurrent /102001, 0500, -499/
      save kount, icurrent, ilast, inc
c
      iskip = 0
      ibreak = 0
      kount = kount + 1
      icurrent = icurrent + inc
      iend     = icurrent + inc 
c     iend     = icurrent + inc -1
      path = '/lustre/scratch/gantech/les_abl/Sullivan/data/'//
     +    'sullivanTest/case037/his.mp.abl.095001_095501'
      indx   = index(path,'.mp.')
      write(path(indx+8:indx+13),'(i6.6)') icurrent
      write(path(indx+15:indx+20),'(i6.6)') iend
      if(icurrent .ge. ilast) then
        done='done'
      endif
      close(lu)
      open(lu,file=path,form='formatted')
c
      return
      end
