      program rdles
c
c ------------- read les history file with profile information
c               specific to 2D MPI code. Included w kurtosis
c               for data comparison
c               modified for 64 bit pgi addressing
c
      include 'par_prof.f'
c
      character path(2)*120
      character done*4, case*3
      character xyfile*40, title*120
      integer lu(2)
      integer*8  ist_1, ien_1, ip_old, ip_int  ! for 64 bit addressing
c      integer iZi, i0p1Zi, i0p5Zi, i0p8Zi !Initial location of zi and 0p1, 0p5 and 0p8 zi
c      real   fcor !The coriolis parameter
c
      done='nope'
      knt = 0
      lu(1) = 99
      lu(2) = 98
      luDqDt = 72
      lu0p1Zi = 93
      lu0p5Zi = 92
      lu0p8Zi = 91      
      luxy  = 88
      luProfiles = 87
      ig    = 1
      iscl  = 1
      num_tot = 0
c
c ------------ set i_mac = 1 for mac (32 bit history files)
c              set i_mac = 0 for sgi
c
      imac = 1
      if(imac .eq. 0) then
         i_byte = 1
      else
         i_byte = 4
      endif
c
c -------------- set iopt depending on how you want the data output

c  iopt = 0, output profiles from database nearest to time t1_start
c            with no interpolation between profiles
c       = 1, output profiles at selected frequency t_freq starting 
c            at t1_start and ending at t2_end using interpolation
c       = 2, output profiles by averaging between times 
c            t1_start and t2_end
c
c -------------- set case and averaging time option
c
c
      case = 'abl'   ! MCBL simulation for First Generation CWF simulations
c
c
c
      if(case .eq. 'abl') then
         xyfile   = 'abl.xy'
         u_gal    = 7.5
         t_ref    = 300.0
         t1_start = 38000.0
         t2_end   = 40000.0
         iopt     = 2
         nnx      = 756
         s_max_i  = 0.08
         fcor = 0.0001
         iZi = 6*maxnz/10
         i0p1Zi = iZi/10
         i0p5Zi = iZi/2
         i0p8Zi = 8*iZi/10
c        t_freq   = 200.0
      else
        write(6,8181)
      endif
      if(iopt .eq. 0) then
         write(6,9001) t1_start
      elseif(iopt .eq. 1) then
         write(6,9002) t1_start, t2_end, t_freq
      elseif(iopt .eq. 2) then
         write(6,9003) t1_start, t2_end
      else
         write(6,0010) iopt
         stop
      endif
c
      nnz = maxnz
c
c ---------- just for this analysis
c
      xl  = 5120.0
      zl  = 2048.0
      dx  = xl/float(nnx)
      dy  = dx
      dzt = zl/float(maxnz)
      delta_f = (dx*1.5*dy*1.5*dzt)**(1.0/3.0)

c
c ---------- Open dQdt file and files to write info at various levels 0p1Zi, 0p5Zi, 0p8Zi
c
      open(luDqDt, file='timeHistoryFiles/dQdt.dat', form='formatted')
      write(luDqDt, 1232)
 1232 format('#t  dQxDt   dQyDt')
      open(lu0p1Zi, file='timeHistoryFiles/0p1iZiTime.dat', 
     +     form='formatted')
      write(lu0p1Zi, 1233) 
 1233 format('#t  uxym  vxum  uule  vvle  wwle  uwle  vwle  wcube'// 
     +   '  wtle  uwsb  vwsb  ttle  wtsb  englez  engsbz')
      open(lu0p5Zi, file='timeHistoryFiles/0p5iZiTime.dat', 
     +     form='formatted')
      write(lu0p5Zi, 1234) 
 1234 format('#t  uxym  vxum  uule  vvle  wwle  uwle  vwle  wcube'// 
     +   '  wtle  uwsb  vwsb  ttle  wtsb  englez  engsbz')
      open(lu0p8Zi, file='timeHistoryFiles/0p8iZiTime.dat', 
     +     form='formatted')
      write(lu0p8Zi, 1235)
 1235 format('#t  uxym  vxum  uule  vvle  wwle  uwle  vwle  wcube'// 
     +   '  wtle  uwsb  vwsb  ttle  wtsb  englez  engsbz')
      
c
c -------- be careful computing length of history file
c
      ist_1   = loc(wwsb(1))
      ien_1   = loc(tke_sprod(maxnz))
      ip_len = (ien_1 - ist_1)/4 + 1
      write(6,8080) ist_1, ip_len
c
      ip_old = loc(wwsb_o(1))
      ip_int = loc(wwsb_i(1))
c
      open(luxy,file=xyfile,form='formatted')
      nnz     = maxnz
      num_avg = 0
      ifile   = 0
      more    = 0
    1 continue
         ifile = ifile + 1
         if(case .eq. 'abl') then
            call rdfile_abl(lu,path,done,iskip,ip_len,ist_1,i_byte)
         else
            write(6,7777)
            stop
         endif
         if(iskip .eq. 1) go to 1
         if(done .eq. 'done') go to 1000
         write(6,6200) path(1), path(2)
         irec = 0
    2    continue
            num_avg = num_avg + 1
            irec = irec + 1
            call swap(ip_len,ist_1,ip_old,more)
            call get_dat_p(lu(1),ist_1,ip_len,irec,more)
            call get_dat_a(lu(2),num_avg,more)
            call get_dz
c
            if(more .eq. 0) then
               call computeWriteDQdt(num_avg, luDqDt, dzt, zl)
               call writeVariousZiLevelData(lu0p1Zi, lu0p5Zi, lu0p8Zi, 
     +              num_avg)
               weit = float(num_avg)
               weit1 = 1.0 - 1.0/weit
               weit2 = 1.0/weit
               call get_zi(num_avg)
               if(iopt .eq. 0) then
                 call opt_0(num_avg,ip_len,ist_1,ip_old,ip_int)
               elseif(iopt .eq. 1) then
                 call opt_1(num_avg,ip_len,ist_1,ip_old,ip_int)
               elseif(iopt .eq. 2) then
                 call opt_2(num_avg,num_tot,
     +                      ip_len,ist_1,ip_old,ip_int)
               endif
c
c ------------- find max gradient
c
                 call save_grad_max(num_avg)
                 call save_theta(num_avg)
               go to 2
            else
               num_avg = num_avg - 1
               go to 1
            endif
 1000 continue
      if(iopt .eq. 2) then
         write(6,9100) num_tot
         call make_prof
         call theta_special(num_avg)
         call wrt_dat(t1_start)
c
         call wrt_theta(num_avg)
      endif
c
c ------------ write file with zi and other statistics
c
      call wrt_zi(num_avg)
c
      write(6,3000) num_avg 
c
      i_gabls = 0
      if(i_gabls .eq. 1) then
         call get_tke(case,ihour,num_tke)
         call wrt_gabls(title,case,ihour,num_avg,num_tke)
      endif
c
c     call xyplot(luxy)
c
      close(luxy)
      close(luDqDt)
      close(lu0p1Zi)
      close(lu0p5Zi)
      close(lu0p8Zi)
      
      stop
 0010 format(' 0010, Illegal iopt: iopt = ',i5)
 3000 format(' number of averages = ',i5)
 6200 format(' working on files = ',/,5x,a120,/,5x,a120)
 7777 format(' bad case, stop')
 8080 format(' 8080 ist_1 = ',i10,' ip_len = ',i6)
 8181 format('STOP, bad case')
 9001 format(' will find profiles closest to t = ',e15.6)
 9002 format(' will find profiles between t1 = ',e15.6,
     +          'and t2 = ',e15.6,/,
     +          ' t frequency = ',e15.6)
 9003 format(' will average profiles between t1 = ',e15.6,
     +          'and t2 = ',e15.6)
 9100 format(' total number of profiles averaged = ',i5)
      end
      subroutine save_grad_max(num_avg)
c
c ------------ save maximum gradient in theta
c
      include 'par_prof.f'
      real grad_th(maxnz)
c
c --------- get max gradient at "zi"
c
      do iz=1,maxnz-1
         if(t_ziavg(num_avg) .gt. zu(iz) .and. 
     +      t_ziavg(num_avg) .le. zu(iz+1)) then
               ipt = iz
               go to 100
         endif
      enddo 
c
  100 continue
      s_max(num_avg) = (txym(ipt+1,1) - txym(ipt,1))/
     +                 (zu(ipt+1) - zu(ipt))
c
c ---------- alternate computation of max gradient
c
c     do iz=1,maxnz-1
c        grad_th(iz) =  (txym(iz+1,1) - txym(iz,1))/
c    +                 (zu(iz+1) - zu(iz))
c     enddo
c
c     iz_min = 5
c     s_temp = grad_th(iz_min)
c     do iz=iz_min,maxnz-5
c        if(grad_th(iz) .gt. s_temp) then
c          s_temp = grad_th(iz)
c        endif
c     enddo
c
c     s_max(num_avg) = s_temp
c
      return
      end
      subroutine save_theta(num_avg)
c
c ------------ save theta profile
c
      include 'par_prof.f'
c
c --------- get max gradient at "zi"
c
      do iz=1,maxnz
         thp(iz,num_avg) = txym(iz,1)
      enddo
c
      do iz=1,maxnz
         wtp(iz,num_avg)   = wttot(iz,1)
         wtsbp(iz,num_avg) = wtsb(iz,1)
      enddo
c
      return
      end
      subroutine make_prof
c
c ----------- make additional profiles from database
c             including eddy viscosity profiles
c
c             edit g/theta_o to fit your problem
c
      include 'par_prof.f'
      data vk /0.4/
      save vk
c
      grav    = 9.81
      theta_o = 300.0
      buoy    = grav/theta_o
c
      do iz=1,nnz-1
         usqr_a = 0.5*(usqr_i(iz) + usqr_i(iz+1))
         vsqr_a = 0.5*(vsqr_i(iz) + vsqr_i(iz+1))
         tke_w(iz) = 0.5*(usqr_a + vsqr_a + wsqr_i(iz))
      enddo
c
c -------- tidy engz_i so that it conforms to the total from usqr.
c          compute at the u points
c
      do iz=2,nnz
         wsqr_a = 0.5*(wsqr_i(iz) + wsqr_i(iz-1))
         engz_i(iz) = 0.5*(usqr_i(iz) + vsqr_i(iz) + wsqr_a)
      enddo
      wsqr_a = 0.5*wsqr_i(1)
      engz_i(1) = 0.5*(usqr_i(1) + vsqr_i(1) + wsqr_a)
c
      do iz=1,nnz
         uttot(iz) = utle_i(iz,iscl) + utsb_i(iz,iscl)
         vttot(iz) = vtle_i(iz,iscl) + vtsb_i(iz,iscl)
         wskew(iz) = 0.0
         wkurt(iz) = 0.0
         tskew(iz) = 0.0
         if(wsqr_i(iz) .le. 0.0) then
            wskew(iz) = 0.0
            wkurt(iz) = 0.0
            tskew(iz) = 0.0
         else
            wskew(iz) = wcube_i(iz)/(wsqr_i(iz)**1.5)
            wkurt(iz) = wfour_i(iz)/(wsqr_i(iz)**2)
            tskew(iz) = tcube_i(iz,iscl)/(tsqr_i(iz,iscl)**1.5)
         endif
      enddo
c
c --------- sweep thru profiles looking for max w**3 and w**2
c
      wsqr_max = 0.0
      isqr_max = 0
      wcub_max = 0.0
      icub_max = 0
      do iz=1,nnz
         if(wsqr_i(iz) .gt. wsqr_max) then
            wsqr_max = wsqr_i(iz)
            isqr_max = iz
         endif
         if(wcube_i(iz) .gt. wcub_max) then
            wcub_max = wcube_i(iz)
            icub_max = iz
         endif
      enddo
c
c --------- get average zi, ustar, etc for this period
c
      n_avg = i2_end + 1 - i1_start
      fn    = 1.0/(n_avg)
      ustar_bar = 0.0
      zi_bar    = 0.0
      qstar_bar = 0.0
      amo_bar   = 0.0
      wmax_bar  = 0.0
      do i=i1_start,i2_end
         weit = 1.0/float(i+1-i1_start)
         weit1 = 1.0 - weit
         ustar_bar = ustar_bar*weit1 + t_utau(i)*weit
         qstar_bar = qstar_bar*weit1 + t_wtsfc(i)*weit
         zi_bar    = zi_bar*weit1 + t_ziavg(i)*weit
         amo_bar   = amo_bar*weit1 + t_amonin(i)*weit
         wmax_bar  = wmax_bar*weit1 + t_wabs(i)*weit
      enddo
c
c -------- use running average from above
c
c     qstar_bar = qstar_bar*fn
c     ustar_bar = ustar_bar*fn
c     amo_bar   = amo_bar*fn
c     zi_bar    = zi_bar*fn
c     wmax_bar  = wmax_bar*fn
      wstar_bar = (zi_bar*qstar_bar*buoy)**(1.0/3.0)
      write(6,6000) n_avg, i1_start, i2_end,
     +              ustar_bar, qstar_bar, zi_bar, amo_bar,
     +              wstar_bar, wmax_bar
c
c --------- get MO functions
c
      z_sl = 0.2*zi_bar
      z_sl = 1.0*zi_bar
      do iz=1,nnz
         if(zw(iz) .lt. z_sl) then
           n_sl = iz
         endif
      enddo
      if(qstar_bar .eq. 0.) then
         tstar_inv = 0.0 
      else
         tstar_inv = -ustar_bar/qstar_bar
      endif
      do iz=1,n_sl
         dudz = (uxym_i(iz+1) - uxym_i(iz))/dzu(iz+1)
         dvdz = (vxym_i(iz+1) - vxym_i(iz))/dzu(iz+1)
         dtdz = (txym_i(iz+1,1) - txym_i(iz,1))/dzu(iz+1)
         phi_m(iz) = vk*zw(iz)*sqrt(dudz**2 + dvdz**2)/ustar_bar
         phi_s(iz) = vk*zw(iz)*dtdz*tstar_inv
      enddo
c
c ------------ make businger mo functions
c
      n_mo = 100
      n_mo = 200
      dz_mo = z_sl/float(n_mo)
      do iz=1,n_mo
         ztemp = float(iz)*dz_mo
         if(qstar_bar .eq. 0.) then
            zeta = 0.0
         else
            zeta = ztemp/amo_bar
         endif
         call busngr(zeta,phi_m_b(iz),phi_s_b(iz),psim,psis)
         call fzol(zeta,phi_m_l(iz),phi_s_l(iz),psim,psis)
         z_mo(iz) = ztemp/zi_bar
      enddo
c
c ------------ get momentun and scalar viscosities
c
      do iz=1,nnz
         flux_tot = sqrt(uwtot_i(iz)**2 + vwtot_i(iz)**2)
         dudz = (uxym_i(iz+1) - uxym_i(iz))/dzu(iz+1)
         dvdz = (vxym_i(iz+1) - vxym_i(iz))/dzu(iz+1)
         dtdz = (txym_i(iz+1,1) - txym_i(iz,1))/dzu(iz+1)
         eddy_m(iz) = (flux_tot/sqrt(dudz**2 + dvdz**2))/
     +                (ustar_bar*zi_bar)
         eddy_s(iz) = (-wttot_i(iz,1)/dtdz)/
     +                (ustar_bar*zi_bar)
      enddo
c
c ---------- sweep thru the heat flux profile and look for the minimum
c            to find the usual z_i
c
      wt_min = qstar_bar
      iz_min = 1
      do iz=1,nnz
         if(wttot_i(iz,1) .lt. wt_min) then
            iz_min = iz
            wt_min = wttot_i(iz,1)
         endif
      enddo
      zi_min_bar = zw(iz_min)
c
c --------- sweep thru dissipation profile and output
c           normalized result
c
      z_mid = zi_bar*0.5
      do iz=1,nnz
         if(zw(iz) .lt. z_mid .and.
     +      zw(iz+1) .gt. z_mid) then
            i_mid = iz
         endif
      enddo
      weit  = (z_mid - zw(i_mid))/(zw(i_mid+1) - zw(i_mid))
      weit1 = 1.0 - weit
      avg_dissp(2) = weit1*tke_diss_i(i_mid) + weit*tke_diss_i(i_mid+1)
      avg_prod(2)  = weit1*shrz_i(i_mid) + weit*shrz_i(i_mid+1)
c
      z_sl = zi_bar*0.1
      do iz=1,nnz
         if(zw(iz) .lt. z_sl .and.
     +      zw(iz+1) .gt. z_sl) then
            i_sl = iz
         endif
      enddo
      weit  = (z_sl - zw(i_sl))/(zw(i_sl+1) - zw(i_sl))
      weit1 = 1.0 - weit
      avg_dissp(1) = weit1*tke_diss_i(i_sl) + weit*tke_diss_i(i_sl+1)
      avg_prod(1)  = weit1*shrz_i(i_sl) + weit*shrz_i(i_sl+1)
c
      z_zi = zi_bar*0.9
      do iz=1,nnz
         if(zw(iz) .lt. z_zi .and.
     +      zw(iz+1) .gt. z_zi) then
            i_zi = iz
         endif
      enddo
      weit  = (z_zi - zw(i_zi))/(zw(i_zi+1) - zw(i_zi))
      weit1 = 1.0 - weit
      avg_dissp(3) = weit1*tke_diss_i(i_zi) + weit*tke_diss_i(i_zi+1)
      avg_prod(3)  = weit1*shrz_i(i_zi) + weit*shrz_i(i_zi+1)
c
c --------- sweep thru total tke profile and output
c           result
c
      do iz=1,nnz
         if(zu(iz) .lt. z_mid .and.
     +      zu(iz+1) .gt. z_mid) then
            i_mid = iz
         endif
      enddo
      weit  = (z_mid - zu(i_mid))/(zu(i_mid+1) - zu(i_mid))
      weit1 = 1.0 - weit
      total_tke(2) = weit1*(engz_i(i_mid+1) + engsbz_i(i_mid+1)) +
     +            weit*(engz_i(i_mid) + engsbz_i(i_mid))
      avg_ebar(2)  = weit1*engsbz_i(i_mid+1) + weit*engsbz_i(i_mid)
c
c -------- at zi
c
      do iz=1,nnz
         if(zu(iz) .lt. z_zi .and.
     +      zu(iz+1) .gt. z_zi) then
            i_zi = iz
         endif
      enddo
c
      weit  = (z_zi - zu(i_zi))/(zu(i_zi+1) - zu(i_zi))
      weit1 = 1.0 - weit
      total_tke(3) = weit1*(engz_i(i_zi+1) + engsbz_i(i_zi+1)) +
     +            weit*(engz_i(i_zi) + engsbz_i(i_zi))
      avg_ebar(3)  = weit1*engsbz_i(i_zi+1) + weit*engsbz_i(i_zi)
c
c -------- at sl
c
      do iz=1,nnz
         if(zu(iz) .lt. z_sl .and.
     +      zu(iz+1) .gt. z_sl) then
            i_sl = iz
         endif
      enddo
c
      weit  = (z_sl - zu(i_sl))/(zu(i_sl+1) - zu(i_sl))
      weit1 = 1.0 - weit
      total_tke(1) = weit1*(engz_i(i_sl+1) + engsbz_i(i_sl+1)) +
     +            weit*(engz_i(i_sl) + engsbz_i(i_sl))
      avg_ebar(1)  = weit1*engsbz_i(i_sl+1) + weit*engsbz_i(i_sl)
c
      return
 6000 format(' Results of averaging',/,
     +       ' n_avg = ',i6,' i1_s = ',i6,' i2_end = ',i6,/,
     +       ' u_star = ',e15.6,' q_star = ',e15.6,/,
     +       ' zi_bar = ',e15.6,' L_mo   = ',e15.6,' w_star = ',e15.6,/,
     +       ' wmax_bar = ',e15.6)
      end
      subroutine theta_special(num_avg)
c
c ----------- get theta profile at a special time
c
      include 'par_prof.f'
c
      time_loc = 15.0*zi_bar/wstar_bar
      write(6,7702) zi_bar, wstar_bar, time_loc
 7702 format(' zi_bar = ',e15.6,/,
     +       ' wstar_bar = ',e15.6,/,
     +       ' time_loc  = ',e15.6)
      write(6,7703) num_avg
 7703 format(' num avg = ',i6)
c
      do k=1,num_avg-1
         if(time(k) .lt. time_loc .and. 
     +      time(k+1) .ge. time_loc) then
            ipt = k
         endif
      enddo
c
      write(6,7701) ipt, time(ipt)
 7701 format(' time location = ',i8,/,
     +       ' time = ',e15.6)
      do k=1,maxnz
         theta_sp(k) = thp(k,ipt)
      enddo
c
c     write(6,8701) (iz,zu(iz),theta_sp(iz),iz=1,maxnz)
 8701 format(' iz ',5x,' zu ',5x,' theta ',/,(i5,2e15.6))
c
      return
      end
      subroutine wrt_theta(num_avg)
c
c ----------- get theta profile at a special time
c
      include 'par_prof.f'
c
c     do ipt=1,num_avg,50
      do ipt=1,num_avg,200
         write(78,7801) ipt, (thp(k,ipt),zu(k),k=1,maxnz)
 7801    format('#n',/,
     +          '#k "ipt = ',i4.4,'"',/,
     +          '#y 0 1500 500',/,
     +          '#x 300 310',/,
     +          (2e15.6))
      enddo
c
      fac = 1.0/0.24
c     do ipt=1,num_avg,50
      do ipt=1,num_avg,200
         write(79,7901) ipt, (wtp(k,ipt)*fac,zw(k),k=1,maxnz)
 7901    format('#n',/,
     +          '#k "ipt = ',i4.4,'"',/,
     +          '#y 0 1500 500',/,
     +          '#x -0.4 1.0',/,
     +          '#m 1',/,
     +          (2e15.6))
         write(79,7902) (wtsbp(k,ipt)*fac,zw(k),k=1,maxnz)
 7902    format('#k',/,
     +          '#m 11',/,
     +          (2e15.6))
      enddo
c
      return
      end
      subroutine get_dz
c
      include 'par_prof.f'
c
      do iz=1,nnz
         dzw(iz) = zw(iz) - zw(iz-1)
         dzu(iz) = zu(iz) - zu(iz-1)
      enddo
      dzw(0) = dzw(1)
      dzw(nnz+1) = dzw(nnz)
      dzu(0) = dzu(1)
      dzu(nnz+1) = dzu(nnz)
c
      return
      end
      subroutine rdfile_abl(lu,path,done,iskip,ip_len,ist,i_byte)
c
c ------- RF7
c
      character path(2)*120, done*4
      integer lu(2)
      data kount /0/
      data ilast,inc, icurrent /102001, 0500, -000499/
      save kount, icurrent, ilast, inc
c
      iskip = 0
      kount = kount + 1
      icurrent = icurrent + inc
      iend     = icurrent + inc
c
      path(1)='/lustre/scratch/gantech/les_abl/Sullivan/data/'//
     +    'sullivanTest/case037/his.mp.abl.095001_095501.ieee'
      path(2)='/lustre/scratch/gantech/les_abl/Sullivan/data/'//
     +    'sullivanTest/case037/his.mp.abl.095001_095501'
c
      indx   = index(path(1),'.mp.')
      do j=1,2
         write(path(j)(indx+8:indx+13),'(i6.6)') icurrent
         write(path(j)(indx+15:indx+20),'(i6.6)') iend
      enddo
      if(icurrent .ge. ilast) then
        done='done'
      else
c
      close(lu(1))
      open(lu(1),file=path(1),form='unformatted',access='direct',
     +     recl=ip_len*i_byte)
c
      close(lu(2))
      open(lu(2),file=path(2),form='formatted')
c
c     write(6,8001) lu(1), path(1), lu(2), path(2)
 8001 format(' lu = ',i4,' path(1) = ',a120,/,
     +       ' lu = ',i4,' path(2) = ',a120)
      endif
c
      return
      end
      subroutine get_dat_p(lu,ist,ip_len,irec,more)
c
      include 'par_prof.f'
      pointer (ist,fprof(ip_len))
c
      more = 0
      read(lu,rec=irec,err=10) (fprof(i),i=1,ip_len)
      go to 99999
c
   10 continue
      more = 1
      write(6,6000) irec
 6000 format(' Hit end of file, irec = ',i5)
c
99999 continue
      return
      end
      subroutine get_dat_a(lu,knt,more)
c
      include 'par_prof.f'
c
      more = 0
c     read(lu,6000,err=10,end=10)
      read(lu,*,err=10,end=10)
     +         time(knt),t_dt(knt),t_utau(knt),
     +         t_ziavg(knt),t_amonin(knt),t_holtop(knt),
     +         t_tsfcc(knt),t_uwsfc(knt),t_vwsfc(knt),
     +         t_divgmax(knt), t_wt_min(knt), t_wt_le(knt),
     +         t_ucfl(knt), t_vcfl(knt), t_wcfl(knt),
     +         t_wtsfc(knt),
     +         t_ups(knt),t_vps(knt),t_wps(knt),t_tps(knt),
     +         t_uwle(knt),t_uwsb(knt),t_uw_tot(knt),
     +         t_vwle(knt),t_vwsb(knt),t_vw_tot(knt),
     +         t_wtle(knt),t_wtsb(knt),t_wt_tot(knt),
     +         t_englez(knt),t_eavg(knt), t_wabs(knt)
c
 6000 format(5e17.8)
c
c ----------- make a time series of the skewness
c
      iz_max = 1
      wsk_max = 0.0
      do iz=1,nnz
         wskew(iz) = 0.0
         if(wsqr(iz) .le. 0.0) then
            wskew(iz) = 0.0
         else
            wskew(iz) = wcube(iz)/(wsqr(iz)**1.5)
            if(wskew(iz) .gt. wsk_max) then
               wsk_max = wskew(iz)
               iz_max  = iz
            endif
         endif
      enddo
      t_wskew(knt) = wsk_max
c
      go to 99999
c
   10 continue
      more = 1
      write(6,6100) knt
 6100 format(' Hit end of file, knt = ',i5)
c
99999 continue
      return
      end
      subroutine computeWriteDQdt(knt, luDqDt, dzt, zl)
c
c ----------- compute dQdt which indicates whether the simulation has reached a quasi steady state or not.
c     
      include 'par_prof.f'

c      integer knt, luDqDt
      real dqxdt, dqydt
      real u1Mag
      u1Mag = sqrt((uxym(1)+u_gal)**2 + vxym(1)**2)
      dqxdt = - (uxym(1)+u_gal)/u1Mag
     + + (fcor * sum(vxym(1:maxnz))*dzt
     + /(t_utau(knt)*t_utau(knt)) )
      dqydt = - vxym(1)/u1Mag 
     + + (fcor * 2.0 * u_gal * zl
     + - fcor * (sum( uxym(1:maxnz)+u_gal ))*dzt)
     + /(t_utau(knt)*t_utau(knt))
c      write(6,2345) t, dqxdt, dqydt
      write(luDqDt,2345)time(knt),dqxdt,dqydt,dzt
 2345 format((4e17.8))
      return
      end
      subroutine writeVariousZiLevelData(lu0p1Zi, lu0p5Zi, lu0p8Zi, knt)

      include 'par_prof.f'

      write(lu0p1Zi,2346)time(knt),uxym(i0p1Zi)+u_gal, vxym(i0p1Zi),
     +     usqr(i0p1Zi), vsqr(i0p1Zi), wsqr(i0p1Zi), uwle(i0p1Zi), 
     +     vwle(i0p1Zi), wcube(i0p1Zi), wtle(i0p1Zi,1), uwsb(i0p1Zi), 
     +     vwsb(i0p1Zi), tsqr(i0p1Zi,1),  wtsb(i0p1Zi,1), 
     +     englez(i0p1Zi), engsbz(i0p1Zi) 
      write(lu0p5Zi,2346)time(knt),uxym(i0p5Zi)+u_gal, vxym(i0p5Zi),
     +     usqr(i0p5Zi), vsqr(i0p5Zi), wsqr(i0p5Zi), uwle(i0p5Zi), 
     +     vwle(i0p5Zi), wcube(i0p5Zi), wtle(i0p5Zi,1), uwsb(i0p5Zi), 
     +     vwsb(i0p5Zi), tsqr(i0p5Zi,1),  wtsb(i0p5Zi,1), 
     +     englez(i0p5Zi), engsbz(i0p5Zi) 
      write(lu0p8Zi,2346)time(knt),uxym(i0p8Zi)+u_gal, vxym(i0p8Zi),
     +     usqr(i0p8Zi), vsqr(i0p8Zi), wsqr(i0p8Zi), uwle(i0p8Zi), 
     +     vwle(i0p8Zi), wcube(i0p8Zi), wtle(i0p8Zi,1), uwsb(i0p8Zi), 
     +     vwsb(i0p8Zi), tsqr(i0p8Zi,1),  wtsb(i0p8Zi,1), 
     +     englez(i0p8Zi), engsbz(i0p8Zi)
 2346 format((16e17.8))
      end
      subroutine wrt_zi(num_avg)
c
c ----------- append xyplot data file with zi 
c             and other time dependent data
c
      include 'par_prof.f'
c
      write(luxy,6000) (time(i),t_ziavg(i),i=1,num_avg)
 6000 format('#n',/,
     +       '#k ',/,
     +       '#lx "time (s)"',/,
     +       '#ly "z\sbi\eb"',/,
     +       '#m 1',/,
     +       '#lw 0.3',/,
     +       '#k "ZI"',/,
     +       (2e17.8))
c
      fac = 1.0/s_max_i
      write(luxy,6001) (time(i),s_max(i)*fac,i=1,num_avg)
 6001 format('#n',/,
     +       '#k ',/,
     +       '#lx "time (s)"',/,
     +       '#ly "Max theta grad"',/,
     +       '#m 1',/,
     +       '#lw 0.3',/,
     +       (2e17.8))
c
      write(luxy,6010) (time(i),t_wtsfc(i),i=1,num_avg)
 6010 format('#n',/,
     +       '#k ',/,
     +       '#lx "time (s)"',/,
     +       '#ly "Q\sb*\eb"',/,
     +       '#m 1',/,
     +       '#lw 0.3',/,
     +       (2e17.8))
c
      write(luxy,6020) (time(i),t_utau(i),i=1,num_avg)
 6020 format('#n',/,
     +       '#k ',/,
     +       '#lx "time (s)"',/,
     +       '#ly "u\sb*\eb (m/s)"',/,
     +       '#m 1',/,
     +       '#y 0.15 0.23 .02',/,
     +       '#lw 0.3',/,
     +       (2e17.8))
c
      write(luxy,6030) (time(i),t_amonin(i),i=1,num_avg)
 6030 format('#n',/,
     +       '#k ',/,
     +       '#lx "time (s)"',/,
     +       '#ly "L (m)"',/,
     +       '#y l 100 1000',/,
     +       '#m 1',/,
     +       '#lw 0.3',/,
     +       (2e17.8))
c
      write(luxy,6040) (time(i),t_wabs(i),i=1,num_avg)
 6040 format('#n',/,
     +       '#k ',/,
     +       '#lx "time (s)"',/,
     +       '#ly "|w| (m/s)"',/,
     +       '#m 1',/,
     +       '#lw 0.3',/,
     +       (2e17.8))
c
      write(luxy,6050) (time(i),t_tsfcc(i),i=1,num_avg)
 6050 format('#n',/,
     +       '#k ',/,
     +       '#lx "time (s)"',/,
     +       '#ly "\gh\sbsfc\eb "',/,
     +       '#y 3 9 2',/,
     +       '#m 1',/,
     +       '#lw 0.3',/,
     +       (2e17.8))
c
      write(luxy,6060) (time(i),t_wskew(i),i=1,num_avg)
 6060 format('#n',/,
     +       '#k ',/,
     +       '#lx "time (s)"',/,
     +       '#ly "wskew"',/,
     +       '#y 0 3 1',/,
     +       '#m 1',/,
     +       '#lw 0.3',/,
     +       (2e17.8))
c
      return
      end
      subroutine wrt_dat(time_output)
c
c ----------- write xyplot data files with profiles at time_output
c
      include 'par_prof.f'
c
      data ionce /0/
      save ionce
      if(ionce .ne. 0) then
        write(luxy,9001)
 9001   format('#n')
      else
        ionce = 1
      endif
c
      iscl = 1
      write(6,5301) u_gal 
 5301 format(' u gal = ',e15.6)
c
      write(luxy,0100) time_output
 0100 format('#k "profiles output at time = ',e15.6)
      write(luxy,0200) zi_bar, qstar_bar, amo_bar,
     +                 wstar_bar, ustar_bar, zi_min_bar,
     +                 wmax_bar, delta_f, wsqr_max, zw(isqr_max),
     +                 wcub_max, zw(icub_max), grav/theta_o,
     +                 avg_dissp(1), avg_prod(1), avg_ebar(1),
     +                 avg_dissp(2), avg_prod(2), avg_ebar(2),
     +                 avg_dissp(3), avg_prod(3), avg_ebar(3),
     +                 total_tke(1),total_tke(2),
     +                 total_tke(3)
 0200 format('#k "Zi_bar = ',e15.6,/,
     +       '#k "Qstar_bar = ',e15.6,/,
     +       '#k "Amo_bar = ',e15.6,/,
     +       '#k "Wstar_bar = ',e15.6,/,
     +       '#k "Ustar_bar = ',e15.6,/,
     +       '#k "Zi_min_bar = ',e15.6,/,
     +       '#k "W_max_bar = ',e15.6,/,
     +       '#k "Delta_f  = ',e15.6,/,
     +       '#k "W_sqrmax_bar = ',e15.6,/,
     +       '#k "Z_w_sqrmax = ',e15.6,/,
     +       '#k "W_cubmax_bar = ',e15.6,/,
     +       '#k "Z_w_cubmax = ',e15.6,/,
     +       '#k "Buoyancy = ',e15.6,/,
     +       '#k "Avg_dissp_sl = ',e15.6,/,
     +       '#k "Avg_prod_sl = ',e15.6,/,
     +       '#k "Avg_ebar_sl = ',e15.6,/,
     +       '#k "Avg_dissp_mid = ',e15.6,/,
     +       '#k "Avg_prod_mid = ',e15.6,/,
     +       '#k "Avg_ebar_mid = ',e15.6,/,
     +       '#k "Avg_dissp_zi = ',e15.6,/,
     +       '#k "Avg_prod_zi = ',e15.6,/,
     +       '#k "Avg_ebar_zi = ',e15.6,/,
     +       '#k "Total_tke_sl   = ',e15.6,/,
     +       '#k "Total_tke_mid  = ',e15.6,/,
     +       '#k "Total_tke_zi  = ',e15.6,/,
     +       '#k "End_of_constants')

      open(luProfiles, file="profiles/uxym", form='formatted')
      write(luProfiles,1000) (uxym_i(iz)+u_gal,vxym_i(iz),
     +     zu(iz),iz=1,nnz)
 1000 format('#<u>, <v>, zu',/,
     +       (3e15.6))
      close(luProfiles)
c
      open(luProfiles, file="profiles/txym", form='formatted')
      write(luProfiles,1003) (txym_i(iz,iscl),zu(iz),iz=1,nnz)
 1003 format('#TXYM  zu',/,
     +       (2e15.6))
      close(luProfiles)
c
      open(luProfiles, file="profiles/uvar", form='formatted')
      write(luProfiles,2000) (usqr_i(iz),vsqr_i(iz),wsqr_i(iz),
     +   zu(iz),iz=1,nnz)
 2000 format('# <u2>, <v2>, <w2>, zu',/,
     +       (4e15.6))
      close(luProfiles)
c
      open(luProfiles, file="profiles/tsqr", form='formatted')
      write(luProfiles,2003) (tsqr_i(iz,iscl),zu(iz),iz=1,nnz)
 2003 format('#<T2> zu ',/,
     +       (2e15.6))
      close(luProfiles)
c
      open(luProfiles, file="profiles/wkurt", form='formatted')
      write(luProfiles,6005) (wkurt(iz), zw(iz), iz=1,nnz)
 6005 format('#WKURTNORM  zw',/,
     +       (2e15.6))
      close(luProfiles)
c
      open(luProfiles, file="profiles/wcube", form='formatted')
      write(luProfiles,2015) (wcube_i(iz), zw(iz), iz=1,nnz)
 2015 format('#WCUBE zw',/,
     +       (2e15.6))
      close(luProfiles)
c
      open(luProfiles, file="profiles/wSkew", form='formatted')
      write(luProfiles,2020) (wskew(iz), zw(iz), iz=1,nnz)
 2020 format('#WSKEW zw',/,
     +       (2e15.6))
      close(luProfiles)
c
      open(luProfiles, file="profiles/tskew", form='formatted')
      write(luProfiles,2005) (tskew(iz), zw(iz),iz=1,nnz)
 2005 format('#TSKEW zw',/,
     +       (2e15.6))
      close(luProfiles)
c
      open(luProfiles, file="profiles/uw", form='formatted')
      write(luProfiles,4000) (uwle_i(iz), uwsb_i(iz), 
     +     uwtot_i(iz), zw(iz),iz=1,nnz)
 4000 format('#UWLE UWSB UWTOT zw',/,
     +       (4e15.6))
      close(luProfiles)
c
      open(luProfiles, file="profiles/vw", form='formatted')
      write(luProfiles,4010) (vwle_i(iz), vwsb_i(iz), 
     +     vwtot_i(iz), zw(iz),iz=1,nnz)
 4010 format('#VWLE VWSB VWTOT zw',/,
     +       (4e15.6))
      close(luProfiles)
c
      open(luProfiles, file="profiles/wt", form='formatted')
      write(luProfiles,4100) (wtle_i(iz,iscl),
     +     wtsb_i(iz,iscl), wttot_i(iz,iscl), 
     +     zw(iz),iz=1,nnz)
 4100 format('#WTLE WTSB WTTOT zw',/,
     +       (4e15.6))
      close(luProfiles)
c
      open(luProfiles, file="profiles/ut", form='formatted')
      write(luProfiles,4200) (utle_i(iz,iscl),
     +     utsb_i(iz,iscl), uttot(iz), 
     +     zw(iz),iz=1,nnz)
 4200 format('#UTLE UTSB UTTOT zu',/,
     +       (4e15.6))
      close(luProfiles)
c
      open(luProfiles, file="profiles/vt", form='formatted')
      write(luProfiles,4300) (vtle_i(iz,iscl),
     +     vtsb_i(iz,iscl), vttot(iz), 
     +     zu(iz),iz=1,nnz)
 4300 format('#VTLE VTSB VTTOT zu',/,
     +       (4e15.6))
      close(luProfiles)
c
      open(luProfiles, file="profiles/englez", form='formatted')
      write(luProfiles,5100) (englez_i(iz),
     +     zw(iz),iz=1,nnz)
 5100 format('#ENGLEZ zw',/,
     +       (2e15.6))
      close(luProfiles)
c
      open(luProfiles, file="profiles/tke", form='formatted')
      write(luProfiles,5101) (engz_i(iz), engsbz_i(iz), 
     + engsbz_i(iz) + engz_i(iz), zu(iz),iz=1,nnz)
 5101 format('#TKELE TKESGS TKETOT zu',/,
     +       (4e15.6))
      close(luProfiles)
c
      open(luProfiles, file="profiles/tke_w", form='formatted')
      write(luProfiles,5104) (tke_w(iz),zw(iz),iz=1,nnz)
 5104 format('#TKEWLE zw',/,
     +       (2e15.6))
      close(luProfiles)
c
c ----------- write terms in tke budget
c
      fac = zi_bar/wstar_bar**3
      write(luxy,8000) (tke_rprod_i(iz)*fac,zw(iz),iz=1,nnz-1)
 8000 format('#n',/,
     +       '#k " tke budget analysis"',/,
     +       '#k " Shear production"',/,
     +       '#y 0 1500 500',/,
     +       '#ly "z (m)"',/,
     +       '#m 1',/,
     +       '#lw 1.0',/,
     +       '#k "KE-SR"',/,
     +       (2e15.6))
      write(luxy,8001) (tke_wp_i(iz)*fac,zw(iz),iz=1,nnz-1)
 8001 format('#k',/,
     +       '#k " pressure transport"',/,
     +       '#y 0 1500 500',/,
     +       '#m 2',/,
     +       '#lw 1.0',/,
     +       '#k "KE-WP"',/,
     +       (2e15.6))
      write(luxy,8002) (tke_wq_i(iz)*fac,zw(iz),iz=1,nnz-1)
 8002 format('#k',/,
     +       '#k " turbulent transport"',/,
     +       '#y 0 1500 500',/,
     +       '#m 3',/,
     +       '#lw 1.0',/,
     +       '#k "KE-WQ"',/,
     +       (2e15.6))
      write(luxy,8003) (tke_buoy_i(iz)*fac,zw(iz),iz=1,nnz-1)
 8003 format('#k',/,
     +       '#k " buoyancy production"',/,
     +       '#y 0 1500 500',/,
     +       '#m 4',/,
     +       '#lw 1.0',/,
     +       '#k "KE-BU"',/,
     +       (2e15.6))
      write(luxy,8004) (-tke_diss_i(iz)*fac,zw(iz),iz=1,nnz-1)
 8004 format('#k',/,
     +       '#k " viscous dissipation"',/,
     +       '#y 0 1500 500',/,
     +       '#m 1',/,
     +       '#lw 0.5',/,
     +       '#k "KE-DISSP"',/,
     +       (2e15.6))
c
c ------------ large value at first grid point?
c
      write(luxy,8005) (tke_sprod_i(iz)*fac,zw(iz),iz=2,nnz-1)
 8005 format('#k',/,
     +       '#k " sgs production"',/,
     +       '#y 0 1500 500',/,
     +       '#m 1',/,
     +       '#lw 0.2',/,
     +       '#k "KE-SS"',/,
     +       (2e15.6))
      write(luxy,8006) (tke_tau_i(iz)*fac,zw(iz),iz=1,nnz-1)
 8006 format('#k',/,
     +       '#k " sgs transport"',/,
     +       '#y 0 1500 500',/,
     +       '#m 2',/,
     +       '#lw 0.2',/,
     +       '#k "KE-STRANS"',/,
     +       (2e15.6))
c
c -------------- output mo profiles
c
      open(luProfiles, file="profiles/phim", form='formatted')
      write(luProfiles,9000) (phi_m(iz),zw(iz)/zi_bar,iz=1,n_sl)
 9000 format('#n',/,
     +       '#k " mo function for momentum"',/,
     +       '#lx "\gf\sbm\eb"',/,
     +       '#ly "z/z\sbi\eb"',/,
     +       '#y 0 0.2 0.1',/,
     +       '#x 0 1 1.0',/,
     +       '#cs 0.8',/,
     +       '#c "\bu"',/,
     +       '#lw 0.6',/,
     +       '#k "PHIM"',/,
     +       (2e15.6))
      close(luProfiles)
c
      open(luProfiles, file="profiles/phim_b", form='formatted')
      write(luProfiles,9010) (phi_m_b(iz),z_mo(iz),iz=1,n_mo)
 9010 format('#k',/,
     +       '#k " mo function for momentum businger"',/,
     +       '#c0',/,
     +       '#m 1',/, 
     +       '#lw 0.6',/,
     +       '#k "PHIM_B"',/,
     +       (2e15.6))
      close(luProfiles)
c
      open(luProfiles, file="profiles/phim_l", form='formatted')
      write(luProfiles,9020) (phi_m_l(iz),z_mo(iz),iz=1,n_mo)
 9020 format('#k',/,
     +       '#k " mo function for momentum large"',/,
     +       '#m 2',/, 
     +       '#lw 1.0',/,
     +       '#k "PHIM_L"',/,
     +       (2e15.6))
      close(luProfiles)
c
      open(luProfiles, file="profiles/phis", form='formatted')
      write(luProfiles,9100) (phi_s(iz),zw(iz)/zi_bar,iz=1,n_sl)
 9100 format('#n',/,
     +       '#k " mo function for heat"',/,
     +       '#lx "\gf\sbs\eb"',/,
     +       '#ly "z/z\sbi\eb"',/,
     +       '#y 0 0.2 0.1',/,
     +       '#x 0 1 0.5',/,
     +       '#cs 0.8',/,
     +       '#c "\bu"',/,
     +       '#lw 0.6',/,
     +       '#k "PHIS"',/,
     +       (2e15.6))
      close(luProfiles)
c
      open(luProfiles, file="profiles/phis_b", form='formatted')
      write(luProfiles,9110) (phi_s_b(iz),z_mo(iz),iz=1,n_mo)
 9110 format('#k',/,
     +       '#k " mo function for scalar businger"',/,
     +       '#c0',/,
     +       '#m 1',/, 
     +       '#lw 0.6',/,
     +       '#k "PHIS_B"',/,
     +       (2e15.6))
      close(luProfiles)
c
      open(luProfiles, file="profiles/phis_l", form='formatted')
      write(luProfiles,9120) (phi_s_l(iz),z_mo(iz),iz=1,n_mo)
 9120 format('#k',/,
     +       '#k " mo function for scalar large"',/,
     +       '#m 2',/, 
     +       '#lw 1.0',/,
     +       '#k "PHIS_L"',/,
     +       (2e15.6))
      close(luProfiles)
c
c -------------- output eddy viscosity profiles
c
      open(luProfiles, file="profiles/eddym", form='formatted')
      write(luProfiles,9300) (eddy_m(iz),zw(iz)/zi_bar,iz=1,nnz/2)
 9300 format('#n',/,
     +       '#k " eddy viscosity for momentum"',/,
     +       '#lx "K\sbm\eb/u\sb*\ebz\sbi\eb"',/,
     +       '#ly "z/z\sbi\eb"',/,
     +       '#y 0 1.0 0.2',/,
     +       '#x 0 1.0 0.2',/,
     +       '#c0',/,
     +       '#m 1',/,
     +       '#lw 1.0',/,
     +       '#k "EDDYM"',/,
     +       (2e15.6))
      close(luProfiles)
c
      open(luProfiles, file="profiles/eddys", form='formatted')
      write(luProfiles,9305) (eddy_s(iz),zw(iz)/zi_bar,iz=1,nnz/2)
 9305 format('#n',/,
     +       '#k " eddy viscosity for scalar"',/,
     +       '#lx "K\sbs\eb/u\sb*\ebz\sbi\eb"',/,
     +       '#ly "z/z\sbi\eb"',/,
     +       '#y 0 1.0 0.2',/,
     +       '#x 0 1.0 0.2',/,
     +       '#c0',/,
     +       '#m 1',/,
     +       '#lw 1.0',/,
     +       '#k "EDDYS"',/,
     +       (2e15.6))
      close(luProfiles)

      return
      end
      subroutine get_zi(num_avg)
c
c ------------ find zi from min uw and min wt flux rule
c
      include 'par_prof.f'
c
      real flux_mo(maxnz)
      data smal /1.e-08/
      data crit, zi_low, zi_max /0.05, 50.0, 300.0/
c
      iz = 1
      zimin = zw(iz)
      wtmin = wttot(iz,iscl)
      do iz=2,nnz
         if(wttot(iz,iscl) .lt. wtmin) then
           zimin = zw(iz)
           wtmin = wttot(iz,iscl)
         endif
      enddo
      t_zimin(num_avg) = zimin
c
c ---------------- uw flux rule
c
      surf_mo = t_utau(num_avg)**2
      if(surf_mo .lt. smal) surf_mo = 1.0
      s_inv = 1.0/surf_mo
      do iz=1,nnz
         flux_mo(iz) = s_inv*sqrt(uwle(iz)**2 + vwle(iz)**2)
      enddo
      izloc = 2
      do iz=2,nnz
         if(flux_mo(iz) .lt. crit) then
           izloc = iz
           go to 99
         endif
      enddo
   99 continue
      slope = (zu(izloc) - zu(izloc-1))/
     +        (flux_mo(izloc) - flux_mo(izloc-1))
      zimin = slope*(crit - flux_mo(izloc-1)) + zu(izloc-1)
      zimin = amax1(zi_low, zimin)
      zimin = amin1(zi_max, zimin)
      t_ziuw(num_avg) = zimin/0.95
c
      return
      end
      subroutine opt_0(num_avg,ip_len,ip_st,ip_old,ip_int)
c
c ------------ locate profiles nearest selected time
c
      include 'par_prof.f'
      pointer (ip_st,fprof(1))
      pointer (ip_old,fprof_old(1))
      pointer (ip_int,fprof_i(1))
c
      ifound = 0
c
c ------------ trap profiles around desired time
c
      if(time(num_avg) .gt. t1_start .and.
     +   time(num_avg-1) .le. t1_start) then
        i1_start = num_avg
        ifound = 1
        diff_2 = time(num_avg) - t1_start
        diff_1 = t1_start - time(num_avg-1) 
        if(diff_1 .lt. diff_2) then
           i1_start = num_avg - 1
           i2_end   = i1_start
           do i=1,ip_len
              fprof_i(i) = fprof_old(i)
           enddo
           write(6,6000) t1_start, time(num_avg-1)
 6000      format(' found profiles ',/,
     +            ' selected time           = ',e15.6,/,
     +            ' output profiles at time = ',e15.6)
c
c ---------- output data in xyplot format
c
           call make_prof
           call wrt_dat(time(num_avg-1))
        else
           i1_start = num_avg
           i2_end   = i1_start
           do i=1,ip_len
              fprof_i(i) = fprof(i)
           enddo
           write(6,6100) t1_start, time(num_avg)
 6100      format(' found profiles ',/,
     +            ' selected time           = ',e15.6,/,
     +            ' output profiles at time = ',e15.6)
c
c ---------- output data in xyplot format
c
           call make_prof
           call wrt_dat(time(num_avg))
        endif
      endif
c
      return
      end
      subroutine opt_1(num_avg,ip_len,ip_st,ip_old,ip_int)
c
c ------------ locate all profiles between t1 and t2
c
      include 'par_prof.f'
      pointer (ip_st,fprof(1))
      pointer (ip_old,fprof_old(1))
      pointer (ip_int,fprof_i(1))
c
      data ionce, time_run /0,-1.0/
      save ionce, time_run
c
      ifound = 0
      if(time(num_avg) .lt. t1_start .or.
     +   time(num_avg) .gt. t2_end) go to 999
      if(ionce .eq. 0) then
        time_run = t1_start
        ionce = 1
      endif
c
c ------------ trap profiles around desired time
c
      if(time(num_avg) .ge. time_run .and.
     +   time(num_avg-1) .le. time_run) then
        i1_start = num_avg
        i2_end   = i1_start
        ifound = 1
        weit = (time_run - time(num_avg - 1))/
     +         (time(num_avg) - time(num_avg - 1))
        weit1 = 1.0 - weit
        do i=1,ip_len
           fprof_i(i) = weit*fprof(i) + weit1*fprof_old(i)
        enddo
        write(6,6000) num_avg, time(num_avg), time(num_avg-1), weit
 6000   format(' found time',/,
     +         ' num_avg = ',i6,' time(i) = ',e15.6,
     +         ' time(i-1) = ',e15.6,' weit = ',e15.6)
c
c ---------- output profiles
c
        call make_prof
        call wrt_dat(time_run)
        time_run = time_run + t_freq
      endif
c
  999 continue
c
      return
      end
      subroutine opt_2(num_avg,num_tot,
     +                 ip_len,ip_st,ip_old,ip_int)
c
c ------------ average profiles between t1 and t2
c
      include 'par_prof.f'
      pointer (ip_st,fprof(1))
      pointer (ip_old,fprof_old(1))
      pointer (ip_int,fprof_i(1))
c
      if(time(num_avg) .lt. t1_start) go to 999
      if(time(num_avg) .gt. t2_end)   go to 999
c
c ----------- start averaging
c
      if(num_tot .eq. 0) then
        do i=1,ip_len
           fprof_i(i) = 0.0
        enddo
        num_tot = 1
        write(6,5900) time(num_avg)
 5900   format(' Started averaging at time = ',e15.6)
        i1_start = num_avg
      endif
      i2_end = num_avg
c
c ------------ average profiles between desired times
c
      weit = 1.0/float(num_tot)
      weit1 = 1.0 - weit
      do i=1,ip_len
         fprof_i(i) = weit*fprof(i) + weit1*fprof_i(i)
      enddo
      num_tot = num_tot + 1
c
  999 continue
c
      return
      end
      subroutine swap(ip_len, ip_st, ip_old, more)
c
c ------------ swap new and old profiles for possible interpolation
c              in time
c
      include 'par_prof.f'
      pointer (ip_st,fprof(1))
      pointer (ip_old,fprof_old(1))
c
      data ionce /0/
      save ionce
c
c ----------- special for first time thru
c
      if(ionce .eq. 0) then
         ionce = 1
         do i=1,ip_len
            fprof(i) = 0.0
         enddo
      endif
c
      if(more .eq. 0) then
         do i=1,ip_len
            fprof_old(i) = fprof(i)
         enddo
      endif
c
      return
      end
      subroutine busngr(zeta,phim,phis,psim,psis)
c
c **** Businger's version of similarity theory
c
      data pih /1.57079633/
      if(zeta .lt. 0.) then
         x=(1.-15.*zeta)**0.25
         phim = 1.0/x
         psim=2.*alog((1.+x)/2.)+alog((1.+x*x)/2.)-2.*atan(x)+pih
         if(psim.gt.2.0)psim=2.0
         y=sqrt(1.-9.*zeta)
         phis = 0.74/y
         psis=alog((1.+y)/2.)*2.
      else if(zeta .gt. 0) then
         phim = 1.0 + 4.7*zeta
         phis = 0.74 + 4.7*zeta
         psim = -4.7*zeta
         psis = -4.7*zeta
      else
         phim = 1.0
         phis = 0.74
         psim = 0.0
         psis = 0.0
      endif
      return
      end
      subroutine fzol(zeta,phim,phis,psim,psis)
c        estimate the stability functions for momentum, m
c                                         and scalars,  c
c        from input of the stability parameter zeta = z/L

      data c1/5./
      data a3,b3,a4,b4/1.258,8.382,-28.862,98.9545/
      data zetam,zetas/-0.2,-1.0/
c
      psimu(Y)  = 1.571+2.*(alog(.5*(1.+Y))-atan(Y))+alog(.5+.5*Y**2)
      psisu(Y)  = 2. * alog(0.5+0.5*Y)
      psicu(Y,G)= (1-g)*alog(abs(Y-1.))
     &          + 0.5*(g+2.)*alog(abs(Y**2+Y+1))
     &          - (2.*g+1.) / sqrt(3.) * atan((Y+0.5)*2./sqrt(3.))
      Xm(zol) = (1.-16. * zol)**0.25
      Xs(zol) = sqrt(1.-16. * zol)
      Xc(zol,f) =  abs(1. - f * zol)**(4./3.) / (1. - f * zol)
c
      if(zeta.ge.0.0)       then
c                                          STABLE
      if(zeta.le.1.0) then
        phim = 1 + c1 * zeta
        psim = - c1 * zeta
        phis = phim
        psis = psim
                      else
c                                   use limiting form
        phim = c1 + zeta
        psim = (1.-c1) * (1. + alog(zeta) ) - zeta
        phis = phim
        psis = psim
                      endif

                            else
c                                         UNSTABLE
c                                                  momentum
       if(zeta.ge.zetam) then
         phim = 1. / Xm(zeta)
         psim = psimu(Xm(zeta))
                         else
c                            use convective limit for momentum
         X = (1. - b3/a3 * zeta)**(1./3.)

         fm = a3**(-1./3.)
         phim = fm / Xc(zeta,b3/a3)
         psim = psimu(Xm(zetam))
     *        + psicu(Xc(zeta,b3/a3),fm)
     *        - psicu(Xc(zetam,b3/a3),fm)
                         endif

c                                         UNSTABLE scalars
       if(zeta.ge.zetas) then
         phis = 1. / Xs(zeta)
         psis = psisu(Xs(zeta))
                         else
c                              use convective limit for scalars
         fs =   abs(a4)**(-1./3.) * abs(a4)/a4
         phis = (a4 - b4 * zeta)**(-1./3.)
         psis = psisu(Xs(zetas))
     *        + psicu(Xc(zeta,b4/a4),fs)
     *        - psicu(Xc(zetas,b4/a4),fs)
                         endif

                            endif
       return
       end
       subroutine wrt_gabls(title,case,ihour,num_avg,num_tke)
c
c --------------- write intercomparison files in specific format
c
      include 'par_prof.f'
c
      character case*3
      character title*120, pathg*120
      real fint(1:maxnz+1)
      data bad_val /-0.9999999e07/
c
c ------------- write "A" dataset
c
      nnz = maxnz
      lug = 31
      pathg = './'//case//'/XX'
      ipos = index(pathg,'X')
      write(pathg(ipos+1:ipos+1),'(i1.1)') ihour
c
c ----------- write "A" data set
c
      pathg(ipos:ipos) = 'A'
      close(lug)
      open(lug,file=pathg,form='formatted')
      write(lug,1000) title
 1000 format(a120)
      write(lug,1001) nnz
 1001 format(i5)
      write(lug,1002) (zu(iz),iz=1,nnz)
 1002 format(10e15.7)
      write(lug,1002) (uxym_i(iz)+u_gal,iz=1,nnz)
      write(lug,1002) (vxym_i(iz),iz=1,nnz)
      write(lug,1002) (txym_i(iz,1),iz=1,nnz)
c
c ----------- write "B" data set
c
      pathg(ipos:ipos) = 'B'
      close(lug)
      open(lug,file=pathg,form='formatted')
      write(lug,1000) title
      write(lug,1001) nnz
      write(lug,1002) (zu(iz),iz=1,nnz)
      write(lug,1002) (usqr_i(iz),iz=1,nnz)
      write(lug,1002) (vsqr_i(iz),iz=1,nnz)
      call lin_upts(wsqr_i,fint,nnz)
      write(lug,1002) (fint(iz),iz=1,nnz)
      call lin_upts(wskew,fint,nnz)
      write(lug,1002) (fint(iz),iz=1,nnz)
      write(lug,1002) (engsbz_i(iz),iz=1,nnz)
      write(lug,1002) (tsqr_i(iz,1),iz=1,nnz)
c
c ----------- write "C" data set
c
      pathg(ipos:ipos) = 'C'
      close(lug)
      open(lug,file=pathg,form='formatted')
      write(lug,1000) title
      write(lug,1001) nnz-1
      write(lug,1002) (zw(iz),iz=1,nnz-1)
      write(lug,1002) (uwle_i(iz),iz=1,nnz-1)
      write(lug,1002) (uwsb_i(iz),iz=1,nnz-1)
      write(lug,1002) (vwle_i(iz),iz=1,nnz-1)
      write(lug,1002) (vwsb_i(iz),iz=1,nnz-1)
      write(lug,1002) (wtle_i(iz,1),iz=1,nnz-1)
      write(lug,1002) (wtsb_i(iz,1),iz=1,nnz-1)
      call lin_wpts(utle_i,fint,nnz)
      write(lug,1002) (fint(iz),iz=1,nnz-1)
      call lin_wpts(utsb_i,fint,nnz)
      write(lug,1002) (fint(iz),iz=1,nnz-1)
      call lin_wpts(vtle_i,fint,nnz)
      write(lug,1002) (fint(iz),iz=1,nnz-1)
      call lin_wpts(vtsb_i,fint,nnz)
      write(lug,1002) (fint(iz),iz=1,nnz-1)
c
c ----------- write "D" data set
c
      pathg(ipos:ipos) = 'D'
      close(lug)
      open(lug,file=pathg,form='formatted')
      write(lug,1000) title
      write(lug,1001) num_tke
      write(lug,1002) (z_tke(iz),iz=1,num_tke)
      write(lug,1002) (shrz_tke(iz),iz=1,num_tke)
      write(lug,1002) (bad_val,iz=1,num_tke)
      write(lug,1002) (bprod(iz),iz=1,num_tke)
      write(lug,1002) (tran(iz),iz=1,num_tke)
      write(lug,1002) (diss(iz),iz=1,num_tke)
      write(lug,1002) (tend(iz),iz=1,num_tke)
c
c ----------- write "E" data set
c
      pathg(ipos:ipos) = 'E'
      close(lug)
      open(lug,file=pathg,form='formatted')
      write(lug,1000) title
      write(lug,1001) num_avg
      write(lug,1002) (time(i),i=1,num_avg)
      write(lug,1002) (t_ziuw(i),i=1,num_avg)
      write(lug,1002) (t_wtsfc(i),i=1,num_avg)
      write(lug,1002) (t_utau(i),i=1,num_avg)
      write(lug,1002) (t_amonin(i),i=1,num_avg)
      write(lug,1002) (t_wabs(i),i=1,num_avg)
c
      return
      end
      subroutine lin_upts(f,fi,nnz)
c
c ----------- take data at w points and average
c             to u points
c
      real f(nnz),fi(nnz)
c
      do iz=2,nnz
         fi(iz) = (f(iz) + f(iz-1))*0.5
      enddo
      fi(1) = 1.5*f(1) -0.5*f(2)
c
      return
      end
      subroutine lin_wpts(f,fi,nnz)
c
c ----------- take data at u points and average
c             to w points
c
      real f(*),fi(*)
c
      do iz=1,nnz-1
         fi(iz) = (f(iz) + f(iz+1))*0.5
      enddo
      fi(nnz) = 1.5*f(nnz) -0.5*f(nnz-1)
c
      return
      end
      subroutine get_tke(case,ihour,num_tke)
c
c --------------- write intercomparison files in specific format
c
      include 'par_prof.f'
c
      character case*3, cvar*20
      character pathg*120
c
c ------------- get tke dataset
c              e.g., sb1.hr8.xy
c
      lu = 31
      pathg = './tke_budg/'//case//'.hrX'//'.xy'
      ipos = index(pathg,'X')
      write(pathg(ipos:ipos),'(i1.1)') ihour
      open(lu,file=pathg,form='formatted')
      write(6,9001) pathg
 9001 format(' pathg = ',a120)
c
      nmatch=5
      cvar = ' KE-S'
      call getdat(lu,cvar,shr_xy,z_tke,num_tke,nmatch)
      cvar = ' KE-P'
      call getdat(lu,cvar,pre_xy,z_tke,num_tke,nmatch)
      cvar = ' KE-T'
      call getdat(lu,cvar,tra_xy,z_tke,num_tke,nmatch)
      cvar = ' KE-D'
      call getdat(lu,cvar,dis_xy,z_tke,num_tke,nmatch)
      nmatch=8
      cvar = ' KE-BUOY'
      call getdat(lu,cvar,buo_xy,z_tke,num_tke,nmatch)
      cvar = ' KE-TEND'
      call getdat(lu,cvar,ten_xy,z_tke,num_tke,nmatch)
c
      do iz=1,num_tke
         diss(iz)      = -dis_xy(iz)
         shrz_tke(iz)  = shr_xy(iz)
         bprod(iz)     = buo_xy(iz)
         tran(iz)      = pre_xy(iz) + tra_xy(iz)
         tend(iz)      = ten_xy(iz)
      enddo
c
      return
      end
      subroutine getdat(lu,cvar,stuff,zstuff,n,nmatch)
      dimension stuff(*),zstuff(*)
      character cvar*20
      character temp*128
c
      rewind lu
      n = 0
    1 continue
         read(lu,100,end=999) temp
  100    format(a128)
         iloc = index(temp,cvar(1:nmatch))
      if(iloc .eq. 0) go to 1
c
c ----------------- found data token
c
    2 continue
        read(lu,100,end=900) temp
        if(temp(1:1) .eq. '#') then
          go to 900
        else
          n = n + 1
          read(temp,*) stuff(n),zstuff(n)
        endif
      go to 2
  900 continue
      return
  999 continue
      write(6,6000) cvar
 6000 format(' trouble, cannot find variable ',a20)
      stop
      end
