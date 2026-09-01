!=======================================================================

! Diagnostic information output during run
!
! authors: Elizabeth C. Hunke, LANL

      module icedrv_diagnostics

      use icedrv_kinds
      use icedrv_constants, only: nu_diag, nu_diag_out
      use icedrv_domain_size, only: nx, nilyr, nslyr, ncat
      use icepack_intfc, only: c0
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_query_parameters
      use icepack_intfc, only: icepack_query_tracer_flags, icepack_query_tracer_indices
      use icedrv_system, only: icedrv_system_abort
      use icepack_therm_mushy, only: permeability
      implicit none
      private
      public :: runtime_diags, &
                init_mass_diags, &
                icedrv_diagnostics_debug, &
                print_state

      ! diagnostic output file
      character (len=char_len), public :: diag_file

      ! point print data

      logical (kind=log_kind), public :: &
         print_points         ! if true, print point data

      integer (kind=int_kind), parameter, public :: &
         npnt = 2             ! total number of points to be printed

      character (len=char_len), dimension(nx), public :: nx_names

      ! for water and heat budgets
      real (kind=dbl_kind), dimension(nx) :: &
         pdhi             , & ! change in mean ice thickness (m)
         pdhs             , & ! change in mean snow thickness (m)
         pde                  ! change in ice and snow energy (W m-2)

!=======================================================================

      contains

!=======================================================================

! Writes diagnostic info (max, min, global sums, etc) to standard out
!
! authors: Elizabeth C. Hunke, LANL
!          Bruce P. Briegleb, NCAR
!          Cecilia M. Bitz, UW

      subroutine runtime_diags (dt)

      use icedrv_flux, only: evap, fsnow, frazil
      use icedrv_flux, only: fswabs, flw, flwout, fsens, fsurf, flat
      use icedrv_flux, only: frain
      use icedrv_flux, only: Tair, Qa, fsw, fcondtop
      use icedrv_flux, only: fbot, fcondbot
      use icedrv_flux, only: meltt, meltb, meltl, snoice, phin, w_diag, dSdt_diag
      use icedrv_flux, only: dh0_cumul, da0_cumul  
      use icedrv_flux, only: meltt_cumul, meltb_cumul, melts_cumul, congel_cumul
      use icedrv_flux, only: snoice_cumul, frazil_cumul, meltl_cumul
      use icedrv_flux, only: melttn_cumul, meltbn_cumul, meltsn_cumul, congeln_cumul
      use icedrv_flux, only: dsnown_cumul, snoicen_cumul, meltln_cumul 
      use icedrv_flux, only: dsnow, congel, sst, sss, Tf, fhocn, Tref
      use icedrv_state, only: aice, vice, vsno, trcr, aicen, vicen, vsnon
      use icedrv_state, only: trcrn, g0n, g1n, hLn, hRn
      use icedrv_arrays_column, only: albicen, albsnon, albpndn    

      
      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      ! local variables

      integer (kind=int_kind) :: &
         n, k, nil

      logical (kind=log_kind) :: &
         calc_Tsfc

      ! fields at diagnostic points
      real (kind=dbl_kind) :: & 
         pTair, pfsnow, pfrain, &
         paice, hiavg, hsavg, hbravg, psalt, pTsfc, &
         pevap, pfhocn, perm

      real (kind=dbl_kind), dimension (nx) :: &
         work1, work2
      real (kind=dbl_kind), dimension (nilyr) :: &
         Tinterni   , &
         zqin_print
      real (kind=dbl_kind), dimension (nslyr) :: &
         Tinterns        

      real (kind=dbl_kind) :: &
         Tffresh, rhos, rhow, rhoi

      logical (kind=log_kind) :: tr_brine
      integer (kind=int_kind) :: nt_fbri, nt_Tsfc
      integer (kind=int_kind) :: nt_qice, nt_sice, nt_qsno
      
      character(len=*), parameter :: subname='(runtime_diags)'

      !-----------------------------------------------------------------
      ! query Icepack values
      !-----------------------------------------------------------------

      call icepack_query_parameters(calc_Tsfc_out=calc_Tsfc)
      call icepack_query_tracer_flags(tr_brine_out=tr_brine)
      call icepack_query_tracer_indices(nt_fbri_out=nt_fbri, nt_Tsfc_out=nt_Tsfc, &
				nt_qice_out=nt_qice, nt_qsno_out=nt_qsno, &
				nt_sice_out=nt_sice)
      call icepack_query_parameters(Tffresh_out=Tffresh, rhos_out=rhos, &
           rhow_out=rhow, rhoi_out=rhoi)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      !-----------------------------------------------------------------
      ! NOTE these are computed for the last timestep only (not avg)
      !-----------------------------------------------------------------

      call total_energy (work1)
      call total_salt   (work2)
      
      do n = 1, nx
        !call print_state('Printing state',n)
        pTair = Tair(n) - Tffresh ! air temperature
        pfsnow = fsnow(n)*dt/rhos ! snowfall
        pfrain = frain(n)*dt/rhow ! rainfall
        
        paice = aice(n)           ! ice area           
        hiavg = c0                ! avg snow/ice thickness
        hsavg = c0
        hbravg = c0               ! avg brine thickness
        psalt = c0 
        if (paice /= c0) then
          hiavg = vice(n)/paice
          hsavg = vsno(n)/paice
          if (tr_brine) hbravg = trcr(n,nt_fbri)* hiavg
        endif
        if (vice(n) /= c0) psalt = work2(n)/vice(n)
        pTair = Tref(n)- Tffresh 
        pTsfc = trcr(n,nt_Tsfc)   ! ice/snow sfc temperature
        pevap = evap(n)*dt/rhoi   ! sublimation/condensation
        pdhi(n) = vice(n) - pdhi(n)  ! ice thickness change
        pdhs(n) = vsno(n) - pdhs(n)  ! snow thickness change
        pde(n) =-(work1(n)- pde(n))/dt ! ice/snow energy change
        pfhocn = -fhocn(n)        ! ocean heat used by ice
        call calc_vertical_profile(nilyr,    nslyr,   &
                 aice(n),  vice(n),   vsno(n),  zqin_print, &
                 trcr (n,nt_qice:nt_qice+nilyr-1), Tinterni,    &
                 trcr (n,nt_qsno:nt_qsno+nslyr-1), Tinterns,    &
                 trcr (n,nt_sice:nt_sice+nilyr-1))
        

        !-----------------------------------------------------------------
        ! start spewing
        !-----------------------------------------------------------------
        
        write(nu_diag_out+n-1,899) nx_names(n)
        
        write(nu_diag_out+n-1,*) '                         '
        write(nu_diag_out+n-1,*) '----------atm----------'
        write(nu_diag_out+n-1,900) 'air temperature (C)    = ',pTair
        write(nu_diag_out+n-1,900) 'specific humidity      = ',Qa(n)
        write(nu_diag_out+n-1,900) 'snowfall (m)           = ',pfsnow
        write(nu_diag_out+n-1,900) 'rainfall (m)           = ',pfrain
        if (.not.calc_Tsfc) then
          write(nu_diag_out+n-1,900) 'total surface heat flux= ', fsurf(n)
          write(nu_diag_out+n-1,900) 'top sfc conductive flux= ',fcondtop(n)
          write(nu_diag_out+n-1,900) 'latent heat flux       = ',flat(n)
        else
          write(nu_diag_out+n-1,900) 'shortwave radiation sum= ',fsw(n)
          write(nu_diag_out+n-1,900) 'longwave radiation     = ',flw(n)
        endif
        write(nu_diag_out+n-1,*) '----------ice----------'
        write(nu_diag_out+n-1,900) 'area fraction          = ',aice(n)! ice area
        write(nu_diag_out+n-1,900) 'avg ice thickness (m)  = ',hiavg
        write(nu_diag_out+n-1,900) 'tot ice volume (m)     = ',vice(n) 
        write(nu_diag_out+n-1,900) 'avg snow depth (m)     = ',hsavg
        write(nu_diag_out+n-1,900) 'avg salinity (ppt)     = ',psalt
        write(nu_diag_out+n-1,900) 'avg brine thickness (m)= ',hbravg
        
        if (calc_Tsfc) then
          write(nu_diag_out+n-1,900) 'surface temperature(C) = ',pTsfc ! ice/snow !pTsfc
          write(nu_diag_out+n-1,900) 'absorbed shortwave flx = ',fswabs(n)
          write(nu_diag_out+n-1,900) 'outward longwave flx   = ',flwout(n)
          write(nu_diag_out+n-1,900) 'sensible heat flx      = ',fsens(n)
          write(nu_diag_out+n-1,900) 'latent heat flx        = ',flat(n)
        endif
        write(nu_diag_out+n-1,900) 'subl/cond (m ice)      = ',pevap   ! sublimation/condensation
        write(nu_diag_out+n-1,900) 'top melt (m)           = ',meltt(n)
        write(nu_diag_out+n-1,900) 'bottom melt (m)        = ',meltb(n)
        write(nu_diag_out+n-1,900) 'lateral melt (m)       = ',meltl(n)
        write(nu_diag_out+n-1,900) 'new ice (m)            = ',frazil(n) ! frazil
        write(nu_diag_out+n-1,900) 'congelation (m)        = ',congel(n)
        write(nu_diag_out+n-1,900) 'snow-ice (m)           = ',snoice(n)
        write(nu_diag_out+n-1,900) 'snow change (m)        = ',dsnow(n)
        write(nu_diag_out+n-1,900) 'effective dhi (m)      = ',pdhi(n)   ! ice thickness change
        write(nu_diag_out+n-1,900) 'effective dhs (m)      = ',pdhs(n)   ! snow thickness change
        write(nu_diag_out+n-1,900) 'intnl enrgy chng(W/m^2)= ',pde (n)   ! ice/snow energy change
        write(nu_diag_out+n-1,*) '----------ocn----------'
        write(nu_diag_out+n-1,900) 'sst (C)                = ',sst(n)    ! sea surface temperature
        write(nu_diag_out+n-1,900) 'sss (ppt)              = ',sss(n)    ! sea surface salinity
        write(nu_diag_out+n-1,900) 'freezing temp (C)      = ',Tf(n)     ! freezing temperature
        write(nu_diag_out+n-1,900) 'heat used (W/m^2)      = ',pfhocn    ! ocean heat used by ice (MP: rather, leftover??)
        write(nu_diag_out+n-1,900) 'bottom conductive flux = ',fcondbot(n)  ! bottom conductive heat flux
        write(nu_diag_out+n-1,900) 'bottom heat flux       = ',fbot(n)      ! heat flux at bottom surface of ice (excluding excess) (W/m^2)
        write(nu_diag_out+n-1,900) 'cumul snowmelt (m)          = ',melts_cumul(n)        
        write(nu_diag_out+n-1,900) 'cumul topmelt (m)           = ',meltt_cumul(n)
        write(nu_diag_out+n-1,900) 'cumul bottommelt (m)        = ',meltb_cumul(n)
        write(nu_diag_out+n-1,900) 'cumul lateralmelt (m)       = ',meltl_cumul(n)
        write(nu_diag_out+n-1,900) 'cumul cat vol loss (m)      = ', dh0_cumul(n)
        write(nu_diag_out+n-1,900) 'cumul newice (m)            = ',frazil_cumul(n) ! frazil
        write(nu_diag_out+n-1,900) 'cumul congel (m)            = ',congel_cumul(n)
        write(nu_diag_out+n-1,900) 'cumul snowice (m)           = ',snoice_cumul(n)
        write(nu_diag_out+n-1,900) 'cumul cat area loss (m)     = ', da0_cumul(n)
        do k = 1, nslyr	  
          write(nu_diag_out+n-1,900) 'int layer snow temp (ppt) = ',Tinterns(k)  !internal snow temperature in layer k 
	enddo       
	do k = 1, nilyr
          write(nu_diag_out+n-1,900) 'int layer salinity (ppt) = ',trcr(n,nt_sice+k-1)  ! salinit in each layer k
          write(nu_diag_out+n-1,900) 'int layer ice temp (ppt) = ',Tinterni(k)  ! internal ice temperature in layer k
          write(nu_diag_out+n-1,900) 'int lyr ice entlpy (ppt) = ',zqin_print(k)*1d-8  ! internal ice temperature in layer k
          !print *, 'internal layer enthalpy is: ', trcr(n,nt_qice+k-1), Tinterni(k) 
	enddo
        do k = 1, ncat
          write(nu_diag_out+n-1,900) 'ice cat. Tsfc 		= ',trcrn(n,nt_Tsfc,k)  ! Temp at top surface
          write(nu_diag_out+n-1,900) 'ice cat. areafrac 	= ',aicen(n,k)  ! ocean heat used by ice
          write(nu_diag_out+n-1,900) 'ice cat. volume (m) 	= ',vicen(n,k)  ! ice volume in category k    
          write(nu_diag_out+n-1,900) 'ice cat. snow vol. (m) 	= ',vsnon(n,k)  ! snow volume in category k
          write(nu_diag_out+n-1,900) 'cat. snowmelt (m)         = ',meltsn_cumul(n,k)        
          write(nu_diag_out+n-1,900) 'cat. topmelt (m)          = ',melttn_cumul(n,k)
          write(nu_diag_out+n-1,900) 'cat. bottommelt (m)       = ',meltbn_cumul(n,k)
          write(nu_diag_out+n-1,900) 'cat. lateralmelt (m)      = ',meltln_cumul(n,k)
          write(nu_diag_out+n-1,900) 'cat. congel (m)           = ',congeln_cumul(n,k)
          write(nu_diag_out+n-1,900) 'cat. snowice (m)          = ',snoicen_cumul(n,k) 
          write(nu_diag_out+n-1,900) 'cat. darcy vel. (m)       = ',w_diag(n,k) ! Darcy velocity in the ice
          write(nu_diag_out+n-1,900) 'alb. ice (m)              = ',albicen(n,k)
          write(nu_diag_out+n-1,900) 'alb. sno (m)              = ',albsnon(n,k)
          write(nu_diag_out+n-1,900) 'alb. pnd (m)              = ',albpndn(n,k)    
          write(nu_diag_out+n-1,900) 'g0 itd constant 		= ',g0n(n,k)     ! ITD constants
          write(nu_diag_out+n-1,900) 'g1 itd slope 		= ',g1n(n,k)        ! ITD category slope     
          write(nu_diag_out+n-1,900) 'hL itd left limit 	= ',hLn(n,k)   ! ITD category left boundary
          write(nu_diag_out+n-1,900) 'hR itd right limit 	= ',hRn(n,k)  ! ITD category right boundary                 
          do nil = 1, nilyr
            write(nu_diag_out+n-1,900) 'liquid fraction           = ', phin(n,nil,k)  ! liquid fraction
            write(nu_diag_out+n-1,900) 'salinity changes          = ', dSdt_diag(n,nil,k)*3.6d3*2.4d1  ! dSdt due to brine physics, PSU per day
            perm = permeability(phin(n,nil,k))
            write(nu_diag_out+n-1,900) 'layer permeability        = ', perm  ! ice layer permeability 
          enddo
	enddo 
        
      end do
899   format (43x,a24)
900   format (a25,2x,f24.17)

      end subroutine runtime_diags

!=======================================================================

! Computes global combined ice and snow mass sum
!
! author: Elizabeth C. Hunke, LANL

      subroutine init_mass_diags

      use icedrv_domain_size, only: nx
      use icedrv_state, only: vice, vsno

      integer (kind=int_kind) :: i

      real (kind=dbl_kind), dimension (nx) :: work1

      character(len=*), parameter :: subname='(init_mass_diags)'

      call total_energy (work1)
      do i = 1, nx
         pdhi(i) = vice (i)
         pdhs(i) = vsno (i)
         pde (i) = work1(i)
      enddo

      end subroutine init_mass_diags

!=======================================================================

! Computes total energy of ice and snow in a grid cell.
!
! authors: E. C. Hunke, LANL

      subroutine total_energy (work)

      use icedrv_domain_size, only: ncat, nilyr, nslyr, nx
      use icedrv_state, only: vicen, vsnon, trcrn

      real (kind=dbl_kind), dimension (nx), intent(out) :: &
         work      ! total energy

      ! local variables

      integer (kind=int_kind) :: &
        i, k, n

      integer (kind=int_kind) :: nt_qice, nt_qsno

      character(len=*), parameter :: subname='(total_energy)'

      !-----------------------------------------------------------------
      ! query Icepack values
      !-----------------------------------------------------------------

         call icepack_query_tracer_indices(nt_qice_out=nt_qice, nt_qsno_out=nt_qsno)
         call icepack_warnings_flush(nu_diag)
         if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
             file=__FILE__,line= __LINE__)

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

         work(:) = c0

      !-----------------------------------------------------------------
      ! Aggregate
      !-----------------------------------------------------------------

         do n = 1, ncat
            do k = 1, nilyr
               do i = 1, nx
                  work(i) = work(i) &
                          + trcrn(i,nt_qice+k-1,n) &
                          * vicen(i,n) / real(nilyr,kind=dbl_kind)
               enddo            ! i
            enddo               ! k

            do k = 1, nslyr
               do i = 1, nx
                  work(i) = work(i) &
                          + trcrn(i,nt_qsno+k-1,n) &
                          * vsnon(i,n) / real(nslyr,kind=dbl_kind)
               enddo            ! i
            enddo               ! k
         enddo                  ! n

      end subroutine total_energy

!=======================================================================

! Computes bulk salinity of ice and snow in a grid cell.
! author: E. C. Hunke, LANL

      subroutine total_salt (work)

      use icedrv_domain_size, only: ncat, nilyr, nx
      use icedrv_state, only: vicen, trcrn

      real (kind=dbl_kind), dimension (nx),  &
         intent(out) :: &
         work      ! total salt

      ! local variables

      integer (kind=int_kind) :: &
        i, k, n

      integer (kind=int_kind) :: nt_sice

      character(len=*), parameter :: subname='(total_salt)'

      !-----------------------------------------------------------------
      ! query Icepack values
      !-----------------------------------------------------------------

         call icepack_query_tracer_indices(nt_sice_out=nt_sice)
         call icepack_warnings_flush(nu_diag)
         if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
             file=__FILE__,line= __LINE__)

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

         work(:) = c0

      !-----------------------------------------------------------------
      ! Aggregate
      !-----------------------------------------------------------------

         do n = 1, ncat
            do k = 1, nilyr
               do i = 1, nx
                  work(i) = work(i) &
                          + trcrn(i,nt_sice+k-1,n) &
                          * vicen(i,n) / real(nilyr,kind=dbl_kind)
               enddo            ! i
            enddo               ! k
         enddo                  ! n

      end subroutine total_salt

!=======================================================================
!
! Wrapper for the print_state debugging routine.
! Useful for debugging in the main driver (see ice.F_debug)
!
! author Elizabeth C. Hunke, LANL
!
      subroutine icedrv_diagnostics_debug (plabeld)

      use icedrv_calendar, only: istep1

      character (*), intent(in) :: plabeld

      character(len=*), parameter :: subname='(icedrv_diagnostics_debug)'

      ! printing info for routine print_state

      integer (kind=int_kind), parameter :: &
         check_step = 1439, & ! begin printing at istep1=check_step
         ip = 3               ! i index

      if (istep1 >= check_step) then
         call print_state(plabeld,ip)
      endif

      end subroutine icedrv_diagnostics_debug

!=======================================================================

! This routine is useful for debugging.
! Calls to it should be inserted in the form (after thermo, for example)
!     plabel = 'post thermo'
!     if (istep1 >= check_step) call print_state(plabel,ip)
! 'use ice_diagnostics' may need to be inserted also
! author: Elizabeth C. Hunke, LANL

      subroutine print_state(plabel,i)

      use icedrv_calendar,  only: istep1, time
      use icedrv_domain_size, only: ncat, nilyr, nslyr
      use icedrv_state, only: aice0, aicen, vicen, vsnon, uvel, vvel, trcrn
      use icedrv_flux, only: uatm, vatm, potT, Tair, Qa, flw, frain, fsnow
      use icedrv_flux, only: fsens, flat, evap, flwout
      use icedrv_flux, only: swvdr, swvdf, swidr, swidf, rhoa
      use icedrv_flux, only: frzmlt, sst, sss, Tf, Tref, Qref, Uref
      use icedrv_flux, only: uocn, vocn
      use icedrv_flux, only: fsw, fswabs, fswint_ai, fswthru, scale_factor
      use icedrv_flux, only: alvdr_ai, alvdf_ai, alidf_ai, alidr_ai

      character (*), intent(in) :: plabel

      integer (kind=int_kind), intent(in) :: & 
          i              ! horizontal index

      ! local variables

      real (kind=dbl_kind) :: &
          eidebug, esdebug, &
          qi, qs, Tsnow, &
          puny, Lfresh, cp_ice, &
          rhoi, rhos

      integer (kind=int_kind) :: n, k

      integer (kind=int_kind) :: nt_Tsfc, nt_qice, nt_qsno

      character(len=*), parameter :: subname='(print_state)'

      !-----------------------------------------------------------------
      ! query Icepack values
      !-----------------------------------------------------------------

      call icepack_query_tracer_indices(nt_Tsfc_out=nt_Tsfc, nt_qice_out=nt_qice, &
           nt_qsno_out=nt_qsno)
      call icepack_query_parameters(puny_out=puny, Lfresh_out=Lfresh, cp_ice_out=cp_ice, &
           rhoi_out=rhoi, rhos_out=rhos)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      !-----------------------------------------------------------------
      ! write diagnostics
      !-----------------------------------------------------------------

      write(nu_diag,*) trim(plabel)
      write(nu_diag,*) 'istep1, i, time', &
                        istep1, i, time
      write(nu_diag,*) ' '
      write(nu_diag,*) 'aice0', aice0(i)
      do n = 1, ncat
         write(nu_diag,*) ' '
         write(nu_diag,*) 'n =',n
         write(nu_diag,*) 'aicen', aicen(i,n)
         write(nu_diag,*) 'vicen', vicen(i,n)
         write(nu_diag,*) 'vsnon', vsnon(i,n)
         if (aicen(i,n) > puny) then
            write(nu_diag,*) 'hin', vicen(i,n)/aicen(i,n)
            write(nu_diag,*) 'hsn', vsnon(i,n)/aicen(i,n)
         endif
         write(nu_diag,*) 'Tsfcn',trcrn(i,nt_Tsfc,n)
         write(nu_diag,*) ' '
      enddo                     ! n

      eidebug = c0
      do n = 1,ncat
         do k = 1,nilyr
            qi = trcrn(i,nt_qice+k-1,n)
            write(nu_diag,*) 'qice, cat ',n,' layer ',k, qi
            eidebug = eidebug + qi
            if (aicen(i,n) > puny) then
               write(nu_diag,*)  'qi/rhoi', qi/rhoi
            endif
         enddo
         write(nu_diag,*) ' '
      enddo
      write(nu_diag,*) 'qice(i)',eidebug
      write(nu_diag,*) ' '

      esdebug = c0
      do n = 1,ncat
         if (vsnon(i,n) > puny) then
            do k = 1,nslyr
               qs = trcrn(i,nt_qsno+k-1,n)
               write(nu_diag,*) 'qsnow, cat ',n,' layer ',k, qs
               esdebug = esdebug + qs
               Tsnow = (Lfresh + qs/rhos) / cp_ice
               write(nu_diag,*) 'qs/rhos', qs/rhos
               write(nu_diag,*) 'Tsnow', Tsnow
            enddo
            write(nu_diag,*) ' '
         endif
      enddo
      write(nu_diag,*) 'qsnow(i)',esdebug
      write(nu_diag,*) ' '

      write(nu_diag,*) 'uvel(i)',uvel(i)
      write(nu_diag,*) 'vvel(i)',vvel(i)

      write(nu_diag,*) ' '
      write(nu_diag,*) 'atm states and fluxes'
      write(nu_diag,*) '            uatm    = ',uatm (i)
      write(nu_diag,*) '            vatm    = ',vatm (i)
      write(nu_diag,*) '            potT    = ',potT (i)
      write(nu_diag,*) '            Tair    = ',Tair (i)
      write(nu_diag,*) '            Qa      = ',Qa   (i)
      write(nu_diag,*) '            rhoa    = ',rhoa (i)
      write(nu_diag,*) '            swvdr   = ',swvdr(i)
      write(nu_diag,*) '            swvdf   = ',swvdf(i)
      write(nu_diag,*) '            swidr   = ',swidr(i)
      write(nu_diag,*) '            swidf   = ',swidf(i)
      write(nu_diag,*) '            flw     = ',flw  (i)
      write(nu_diag,*) '            frain   = ',frain(i)
      write(nu_diag,*) '            fsnow   = ',fsnow(i)
      write(nu_diag,*) ' '
      write(nu_diag,*) 'ocn states and fluxes'
      write(nu_diag,*) '            frzmlt  = ',frzmlt (i)
      write(nu_diag,*) '            sst     = ',sst    (i)
      write(nu_diag,*) '            sss     = ',sss    (i)
      write(nu_diag,*) '            Tf      = ',Tf     (i)
      write(nu_diag,*) '            uocn    = ',uocn   (i)
      write(nu_diag,*) '            vocn    = ',vocn   (i)
      write(nu_diag,*) ' '
      write(nu_diag,*) 'srf states and fluxes'
      write(nu_diag,*) '            Tref    = ',Tref  (i)
      write(nu_diag,*) '            Qref    = ',Qref  (i)
      write(nu_diag,*) '            Uref    = ',Uref  (i)
      write(nu_diag,*) '            fsens   = ',fsens (i)
      write(nu_diag,*) '            flat    = ',flat  (i)
      write(nu_diag,*) '            evap    = ',evap  (i)
      write(nu_diag,*) '            flwout  = ',flwout(i)
      write(nu_diag,*) ' '
      write(nu_diag,*) 'shortwave'
      write(nu_diag,*) '            fsw          = ',fsw         (i)
      write(nu_diag,*) '            fswabs       = ',fswabs      (i)
      write(nu_diag,*) '            fswint_ai    = ',fswint_ai   (i)
      write(nu_diag,*) '            fswthru      = ',fswthru     (i)
      write(nu_diag,*) '            scale_factor = ',scale_factor(i)
      write(nu_diag,*) '            alvdr        = ',alvdr_ai    (i)
      write(nu_diag,*) '            alvdf        = ',alvdf_ai    (i)
      write(nu_diag,*) '            alidr        = ',alidr_ai    (i)
      write(nu_diag,*) '            alidf        = ',alidf_ai    (i)
      write(nu_diag,*) ' '

      call icepack_warnings_flush(nu_diag)

      end subroutine print_state

!=======================================================================
!
! Given the state variables (vicen, vsnon, zqin, etc.),
! compute variables needed for the vertical thermodynamic profiles
! (hin, hsn, zTin, etc.)
!
! The internal temperature is not a global variable, rather being computed 
! from other thermodynamic variables.
! This subroutine corresponds to subroutine init_vertical_profile, but without
! updating the tracers. 
! (See module icepack_therm_vertical, authors William H. Lipscomb, LANL ;
!      C. M. Bitz, UW)
! 
      subroutine calc_vertical_profile(nilyr,    nslyr,    &
                                       aicen,    vicen,    &
                                       vsnon,    zqin_out, &
                                       zqin,     zTin,     &
                                       zqsn,     zTsn,     &
                                       zSin )
                                       
      use icepack_therm_shared, only: l_brine
      use icepack_parameters, only: ktherm, depressT, puny, c1, heat_capacity
      use icepack_parameters, only: rhos, Lfresh, hs_min, cp_ice, min_salin
      use icepack_warnings, only: warnstr, icepack_warnings_add
      use icepack_warnings, only: icepack_warnings_setabort, icepack_warnings_aborted   
      use icepack_mushy_physics, only: temperature_mush
      use icepack_mushy_physics, only: liquidus_temperature_mush
      use icepack_mushy_physics, only: enthalpy_mush, enthalpy_of_melting
      use icepack_therm_shared, only: calculate_tin_from_qin, Tmin
      

                                                                                 
      integer (kind=int_kind), intent(in) :: &
         nilyr , & ! number of ice layers
         nslyr     ! number of snow layers

      real (kind=dbl_kind), intent(in) :: &
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)
 
      real (kind=dbl_kind) :: &
         hilyr       , & ! ice layer thickness
         hslyr        ! snow layer thickness
 
      real (kind=dbl_kind) :: &
         hin         , & ! ice thickness (m)
         hsn             ! snow thickness (m)

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         zqin            ! ice layer enthalpy (J m-3)        
      real (kind=dbl_kind), dimension (:), intent(out) :: &
         zTin        , & ! internal ice layer temperatures 
         zqin_out        ! internal ice layer enthalpy to print 
      real (kind=dbl_kind), dimension (:), intent(out) :: & 
         zTsn            ! snow temperature 
      real (kind=dbl_kind), dimension (:), intent(in) :: &
         zSin            ! internal ice layer salinities        
      real (kind=dbl_kind), dimension (:), intent(in) :: &
         zqsn            ! snow enthalpy

      ! local variables
      real (kind=dbl_kind), dimension(nilyr) :: &
         Tmlts          ! melting temperature

      integer (kind=int_kind) :: &
         k               ! ice layer index

      real (kind=dbl_kind) :: &
         rnslyr,       & ! real(nslyr)
         zqsn_dum,     & ! melting temperature
         zqin_dum,     & ! melting temperature
         Tmax            ! maximum allowed snow/ice temperature (deg C)

      logical (kind=log_kind) :: &   ! for vector-friendly error checks
         tsno_high   , & ! flag for zTsn > Tmax
         tice_high   , & ! flag for zTin > Tmlt
         tsno_low    , & ! flag for zTsn < Tmin
         tice_low        ! flag for zTin < Tmin

      character(len=*),parameter :: subname='(out_vertical_profile)'

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      rnslyr = real(nslyr,kind=dbl_kind)

      tsno_high = .false.
      tice_high = .false.
      tsno_low  = .false.
      tice_low  = .false.

      !-----------------------------------------------------------------
      ! Surface temperature, ice and snow thickness
      ! Initialize internal energy
      !-----------------------------------------------------------------
      if (aicen .eq. 0d0) then
	zTsn(:) = 0d0
	zTin(:) = 0d0            
	
	zqin_out(:) = 0d0
      else
      
	hin    = vicen / aicen
	hsn    = vsnon / aicen
	hilyr    = hin / real(nilyr,kind=dbl_kind)
	hslyr    = hsn / rnslyr
	    print *, 'hin, hsn, nilyr : ', hin, hsn, nilyr
      !-----------------------------------------------------------------
      ! Snow enthalpy and maximum allowed snow temperature
      ! If heat_capacity = F, zqsn and zTsn are used only for checking
      ! conservation.
      !-----------------------------------------------------------------

	do k = 1, nslyr
	    print *, 'zqsn, rho, Lfresh : ', zqsn(k), rhos, Lfresh
      !-----------------------------------------------------------------
      ! Tmax based on the idea that dT ~ dq / (rhos*cp_ice)
      !                             dq ~ q dv / v
      !                             dv ~ puny = eps11
      ! where 'd' denotes an error due to roundoff.
      !-----------------------------------------------------------------

	  if (hslyr > hs_min/rnslyr .and. heat_capacity) then
            ! zqsn < 0              
            Tmax = -zqsn(k)*puny*rnslyr / &
                 (rhos*cp_ice*vsnon)
            zqsn_dum = zqsn(k)                 
	  else
            zqsn_dum= -rhos * Lfresh
            Tmax = puny
	  endif

      !-----------------------------------------------------------------
      ! Compute snow temperatures from enthalpies.
      ! Note: zqsn <= -rhos*Lfresh, so zTsn <= 0.
      !-----------------------------------------------------------------
	  zTsn(k) = (Lfresh + zqsn_dum/rhos)/cp_ice

      !-----------------------------------------------------------------
      ! Check for zTsn > Tmax (allowing for roundoff error) and zTsn < Tmin.
      !-----------------------------------------------------------------
	  if (zTsn(k) > Tmax) then
            tsno_high = .true.
	  elseif (zTsn(k) < Tmin) then
            tsno_low  = .true.
	  endif

	enddo                     ! nslyr

      !-----------------------------------------------------------------
      ! If zTsn is out of bounds, print diagnostics and exit.
      !-----------------------------------------------------------------

	if (tsno_high .and. heat_capacity) then
	  do k = 1, nslyr

            if (hslyr > hs_min/rnslyr) then
               Tmax = -zqsn(k)*puny*rnslyr / &
                    (rhos*cp_ice*vsnon)
            else
               Tmax = puny
            endif
	    print *, 'zTsn, Tmax : ', zTsn(k), Tmax
            if (zTsn(k) > Tmax) then
               write(warnstr,*) ' '
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, 'output thermo, zTsn > Tmax'
               call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
               call icepack_warnings_add(subname//" out_vertical_profile: Starting thermo, zTsn > Tmax" ) 
               return
            endif

	  enddo                  ! nslyr
	endif                     ! tsno_high

	if (tsno_low .and. heat_capacity) then
	  do k = 1, nslyr

            if (zTsn(k) < Tmin) then ! allowing for roundoff error
               write(warnstr,*) ' '
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, 'output thermo, zTsn < Tmin'
               return
            endif

	  enddo                  ! nslyr
	endif                     ! tsno_low

	do k = 1, nslyr

	  if (zTsn(k) > c0) then   ! correct roundoff error
            zTsn(k) = c0
	  endif

	enddo                     ! nslyr

	do k = 1, nilyr

      !---------------------------------------------------------------------
      !  Use initial salinity profile for thin ice
      !---------------------------------------------------------------------

	  if (ktherm == 1 .and. zSin(k) < min_salin-puny) then
            write(warnstr,*) ' '
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, 'output zSin < min_salin, layer', k
            return
	  endif
         
	  if (ktherm == 2) then
            Tmlts(k) = liquidus_temperature_mush(zSin(k))
	  else
            Tmlts(k) = -zSin(k) * depressT
	  endif
	  

      !-----------------------------------------------------------------
      ! Compute ice enthalpy
      ! If heat_capacity = F, zqin and zTin are used only for checking
      ! conservation.
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! Compute ice temperatures from enthalpies using quadratic formula
      !-----------------------------------------------------------------
         
	  if (ktherm == 2) then
            zTin(k) = temperature_mush(zqin(k),zSin(k))
	        zqin_out(k) = zqin(k)
	  else
            zTin(k) = calculate_Tin_from_qin(zqin(k),Tmlts(k))
            zqin_out(k) = zqin(k)
            print *, 'Enthalpy out', zqin(k)
	  endif

	  if (l_brine) then
            Tmax = Tmlts(k)
	  else                ! fresh ice
            Tmax = -zqin(k)*puny/(rhos*cp_ice*vicen)
	  endif

      !-----------------------------------------------------------------
      ! Check for zTin > Tmax and zTin < Tmin
      !-----------------------------------------------------------------
	  if (zTin(k) > Tmax) then
            tice_high = .true.
	  elseif (zTin(k) < Tmin) then
            tice_low  = .true.
	  endif

      !-----------------------------------------------------------------
      ! If zTin is out of bounds, print diagnostics and exit.
      !-----------------------------------------------------------------

	  if (tice_high .and. heat_capacity) then

            if (l_brine) then
               Tmax = Tmlts(k)
            else             ! fresh ice
               Tmax = -zqin(k)*puny/(rhos*cp_ice*vicen)
            endif

            if (zTin(k) > Tmax) then               
               if (ktherm == 2) then
                  zqin_dum = enthalpy_of_melting(zSin(k)) - c1
                  zTin(k) = temperature_mush(zqin_dum,zSin(k))
                  zqin_out(k) = zqin_dum
                  write(warnstr,*) subname, 'Corrected quantities'
               else
                  call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
                  call icepack_warnings_add(subname//" init_vertical_profile: Starting thermo, T > Tmax, layer" ) 
                  return
               endif
            endif
	  endif                  ! tice_high

	  if (tice_low .and. heat_capacity) then
            if (zTin(k) < Tmin) then
               write(warnstr,*) ' '
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, 'Starting thermo T < Tmin, layer', k
               return
            endif
	  endif                  ! tice_low

      !-----------------------------------------------------------------
      ! correct roundoff error
      !-----------------------------------------------------------------

	  if (ktherm /= 2) then

            if (zTin(k) > c0) then
               zTin(k) = c0
            endif

	  endif

	enddo                     ! nilyr
      endif ! if aicen = 0
      
      
      end subroutine calc_vertical_profile
      
!=======================================================================     
      
      
      
!=======================================================================

      end module icedrv_diagnostics

!=======================================================================
