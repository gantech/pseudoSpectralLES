

c
c --------------- for case aa3
c



c
c --------------- for case rf4
c



c
c --------------- for case rf0
c



c
c --------------- for case rg0
c



c
c --------------- for case rf1
c



c
c --------------- for case rg1
c



c
c --------------- for case rf5
c



c
c --------------- for case rf2
c



c
c --------------- for case rf6
c



c
c --------------- for case rf7
c



c
c --------------- for case rf3
c



c
c --------------- for case sf3
c



c
c --------------- for case rg3
c



c
c --------------- for case rg4
c



c
c --------------- for case rg5
c



c
c --------------- for case rh1
c



c
c --------------- for case rh2
c



c
c --------------- for case ba1
c



c
c --------------- for case abl
c

      parameter (maxnz=256,maxnz1=maxnz+1)

c
c --------------- for case abl
c
c
      parameter (nscl=1)
      parameter (nq=110000)
      parameter (nmox = 200)
c ----------------------------------------------------------------------
      common /iia/ fcor, ig, iscl, nnz, iopt, luxy, i1_start, i2_end, n_sl,
     + n_mo, isqr_max, icub_max, iZi, i0p1Zi, i0p5Zi, i0p8Zi
c ----------------------------------------------------------------------
      common /aa/ dzw(0:maxnz1),dzu(0:maxnz1),
     + wwsb(maxnz),engz(0:maxnz1),
     + engsbz(0:maxnz1),
     + englez(maxnz),uxym(0:maxnz1),
     + vxym(0:maxnz1),wxym(0:maxnz1),
     + txym(0:maxnz1,nscl),divz(0:maxnz1),
     + utle(maxnz,nscl), utsb(maxnz,nscl),
     + vtle(maxnz,nscl), vtsb(maxnz,nscl),
     + wtle(maxnz,nscl), wtsb(maxnz,nscl),
     + wttot(maxnz,nscl),
     + zw(0:maxnz1),zu(0:maxnz1),
     + shrz(maxnz),buyz(maxnz),
     + triz(maxnz),
     + uwsb(maxnz),vwsb(maxnz),
     + uwle(maxnz),vwle(maxnz),
     + uwtot(maxnz),vwtot(maxnz),
     + wcube(maxnz), wfour(maxnz),
     + tcube(maxnz,nscl),
     + usqr(maxnz), vsqr(maxnz),
     + wsqr(maxnz), tsqr(maxnz,nscl),
     + tke_rprod(maxnz), tke_wq(maxnz),
     + tke_wp(maxnz), tke_tau(maxnz),
     + tke_tran(maxnz), tke_buoy(maxnz),
     + tke_diss(maxnz), tke_sprod(maxnz)
c ----------------------------------------------------------------------
      common /aa_o/ wwsb_o(maxnz),engz_o(0:maxnz1),
     + engsbz_o(0:maxnz1),
     + englez_o(maxnz),uxym_o(0:maxnz1),
     + vxym_o(0:maxnz1),wxym_o(0:maxnz1),
     + txym_o(0:maxnz1,nscl),divz_o(0:maxnz1),
     + utle_o(maxnz,nscl), utsb_o(maxnz,nscl),
     + vtle_o(maxnz,nscl), vtsb_o(maxnz,nscl),
     + wtle_o(maxnz,nscl), wtsb_o(maxnz,nscl),
     + wttot_o(maxnz,nscl),
     + zw_o(0:maxnz1),zu_o(0:maxnz1),
     + shrz_o(maxnz),buyz_o(maxnz),
     + triz_o(maxnz),
     + uwsb_o(maxnz),vwsb_o(maxnz),
     + uwle_o(maxnz),vwle_o(maxnz),
     + uwtot_o(maxnz),vwtot_o(maxnz),
     + wcube_o(maxnz), wfour_o(maxnz),
     + tcube_o(maxnz,nscl),
     + usqr_o(maxnz), vsqr_o(maxnz),
     + wsqr_o(maxnz), tsqr_o(maxnz,nscl),
     + tke_rprod_o(maxnz), tke_wq_o(maxnz),
     + tke_wp_o(maxnz), tke_tau_o(maxnz),
     + tke_tran_o(maxnz), tke_buoy_o(maxnz),
     + tke_diss_o(maxnz), tke_sprod_o(maxnz)
c ----------------------------------------------------------------------
      common /aa_i/ wwsb_i(maxnz),engz_i(0:maxnz1),
     + engsbz_i(0:maxnz1),
     + englez_i(maxnz),uxym_i(0:maxnz1),
     + vxym_i(0:maxnz1),wxym_i(0:maxnz1),
     + txym_i(0:maxnz1,nscl),divz_i(0:maxnz1),
     + utle_i(maxnz,nscl), utsb_i(maxnz,nscl),
     + vtle_i(maxnz,nscl), vtsb_i(maxnz,nscl),
     + wtle_i(maxnz,nscl), wtsb_i(maxnz,nscl),
     + wttot_i(maxnz,nscl),
     + zw_i(0:maxnz1),zu_i(0:maxnz1),
     + shrz_i(maxnz),buyz_i(maxnz),
     + triz_i(maxnz),
     + uwsb_i(maxnz),vwsb_i(maxnz),
     + uwle_i(maxnz),vwle_i(maxnz),
     + uwtot_i(maxnz),vwtot_i(maxnz),
     + wcube_i(maxnz), wfour_i(maxnz),
     + tcube_i(maxnz,nscl),
     + usqr_i(maxnz), vsqr_i(maxnz),
     + wsqr_i(maxnz), tsqr_i(maxnz,nscl),
     + tke_rprod_i(maxnz), tke_wq_i(maxnz),
     + tke_wp_i(maxnz), tke_tau_i(maxnz),
     + tke_tran_i(maxnz), tke_buoy_i(maxnz),
     + tke_diss_i(maxnz), tke_sprod_i(maxnz),
     + thp(maxnz,nq), wtp(maxnz,nq), wtsbp(maxnz,nq)
c ----------------------------------------------------------------------
      common /bb_i/ wskew(maxnz), uttot(maxnz), vttot(maxnz),
     + tskew(maxnz), wkurt(maxnz)
c ----------------------------------------------------------------------
      common /cc/ grav, theta_o, t_ref, u_gal, ustar_bar, qstar_bar,
     + wstar_bar, zi_bar, amo_bar, wmax_bar, s_max_i,
     + delta_f, wsqr_max, total_tke(3),
     + wcub_max, avg_dissp(3), avg_prod(3),
     + avg_ebar(3),
     + phi_m(maxnz), phi_s(maxnz),
     + phi_m_b(nmox), phi_s_b(nmox),
     + phi_m_l(nmox), phi_s_l(nmox),
     + z_mo(nmox), zi_min_bar,
     + t1_start, t2_end, t_freq,
     + eddy_m(maxnz), eddy_s(maxnz),
     + theta_sp(maxnz), tke_w(maxnz)
c ----------------------------------------------------------------------
      common /ee/ time(nq),t_dt(nq),t_utau(nq),
     + t_ziavg(nq),t_amonin(nq),t_holtop(nq),
     + t_tsfcc(nq),t_uwsfc(nq),t_vwsfc(nq),
     + t_divgmax(nq), t_wt_min(nq), t_wt_le(nq),
     + t_ucfl(nq), t_vcfl(nq), t_wcfl(nq),
     + t_wtsfc(nq),
     + t_ups(nq),t_vps(nq),t_wps(nq),t_tps(nq),
     + t_uwle(nq),t_uwsb(nq),t_uw_tot(nq),
     + t_vwle(nq),t_vwsb(nq),t_vw_tot(nq),
     + t_wtle(nq),t_wtsb(nq),t_wt_tot(nq),
     + t_englez(nq),t_eavg(nq), t_wabs(nq), t_zimin(nq),
     + t_ziuw(nq), t_wskew(nq), s_max(nq)
c ----------------------------------------------------------------------
      common /ff/ z_tke(maxnz), shrz_tke(maxnz), bprod(maxnz),
     + tran(maxnz), diss(maxnz), tend(maxnz),
     + shr_xy(maxnz), buo_xy(maxnz), tra_xy(maxnz),
     + pre_xy(maxnz), dis_xy(maxnz), ten_xy(maxnz)
c ----------------------------------------------------------------------
