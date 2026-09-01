MODULE gemdrv_load_forcing

   !!======================================================================
   !!                       ***  MODULE  load GEM forcing  ***
   !!=====================================================================
   !! History :  1.0   !  10-18  created by M. Plante 
   !! 			            based on code in CONCEPTS : sbcblk_rpn.F90
   !!----------------------------------------------------------------------  
   !!----------------------------------------------------------------------
   !!   load_data    : Compute the atmospheric fluxes over ice
   !!   test_output_rpn : Make outputs of the computed field.
   !!			  This is only to verify that the data is unchanged
   !!----------------------------------------------------------------------


   USE gemdrv_read_rpn
   USE icedrv_kinds   
   USE icedrv_calendar   
 

   IMPLICIT NONE
   PUBLIC

   INTEGER , PARAMETER ::   jpfld   = 15          ! maximum number of files to read
   INTEGER , PARAMETER ::   jp_wndi = 1           ! index of 10m (or zzu) wind velocity (i-component) (m/s)    at T-point
   INTEGER , PARAMETER ::   jp_wndj = 2           ! index of 10m (or zzu) wind velocity (j-component) (m/s)    at T-point
   INTEGER , PARAMETER ::   jp_humi = 3           ! index of specific humidity               ( Kg/Kg )
   INTEGER , PARAMETER ::   jp_qsr  = 4           ! index of solar heat                      (W/m2)
   INTEGER , PARAMETER ::   jp_qlw  = 5           ! index of Long wave                       (W/m2)
   INTEGER , PARAMETER ::   jp_tair = 6           ! index of 10m (or zzu) air temperature             (Kelvin)
   INTEGER , PARAMETER ::   jp_prec = 7           ! index of total precipitation (rain+snow) (Kg/m2/s)
   INTEGER , PARAMETER ::   jp_snow = 8           ! index of snow (solid prcipitation)       (kg/m2/s)
   INTEGER , PARAMETER ::   jp_hght = 9           ! index of thermo   level (m)
   INTEGER , PARAMETER ::   jp_hghm = 10          ! index of momentum level (m)
   INTEGER , PARAMETER ::   jp_prst = 11          ! index of thermo   level pressure (Pa)
   INTEGER , PARAMETER ::   jp_prsm = 12          ! index of momentum level pressure (Pa)
   INTEGER , PARAMETER ::   jp_prsg = 13          ! index of sea (ground) level pressure (Pa)
   INTEGER , PARAMETER ::   jp_mcmx = 14          ! index of mixt coupling/forcing mask (forcing==>sf(jp_mcmx)=1)
   INTEGER , PARAMETER ::   jp_prsn = 15          ! index of PN (virtual ground pressure) (Pa)

   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf   ! structure of input fields (file informations, fields read)

!-------------------------------------------------------------------------------
!!BLOCK FROM sbc_rpn and sbc_oce. 
   !                                             !!! RPNE bulk parameters
   REAL(wp), PARAMETER  ::   albo =    0.066      ! ocean albedo assumed to be constant
   REAL(wp), PARAMETER  ::   zeps0 = 1.e-13       ! small number
   REAL(wp), PARAMETER  ::   rhor_ref = 1.22      ! reference air density
   REAL(wp), PARAMETER  ::   zu_diag = 10.        ! hard coded diagnostic level for winds
   REAL(wp), PARAMETER  ::   zt_diag = 1.5        ! hard coded diagnostic level for T and q
   REAL(wp), SAVE       ::   Cice_t               ! transfer coefficient over ice (thermo level)
   REAL(wp), SAVE       ::   Cice_m               ! transfer coefficient over ice (momentum level)
   
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::     tprecip           !: total precipitation                          [Kg/m2/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::     sprecip           !: solid precipitation                          [Kg/m2/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   cat_i_rpn !: potential temperature (K)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   tatm_ice       !: air temperature [K]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   qatm_ice           !: specific humidity
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   zlm_i_rpn !: atmos momentum level height (m)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   zlt_i_rpn !: atmos thermo   level height (m)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   wndi_ice           !: i wind at T point
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   wndj_ice           !: j wind at T point
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   qlw_ice            !: incoming long-wave
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   qsr_ice        !: solar heat flux over ice                      [W/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   fr0_i_rpn !: ice fraction if coupled with CICE (begining of time step)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   t4a_i_rpn !: radiative temperature [K**4]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   qsb_i_rpn !: sensible heat flux [W/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   qlw_i_rpn !: upward long wave [W/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   rhor_rpn             !: air density (km/m3)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   ustr_rpn             !: friction velocity U* (m/s)
      
!-------------------------------------------------------   

   public :: sf,jp_wndi,jp_wndj,jp_tair,jp_humi

   
CONTAINS

   SUBROUTINE load_data( kt,GEM_rpn_list)
        implicit none
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE load_data  ***
      !! Based on the routine atm_blk_rpn in the CONCEPTS code
      !! ** Purpose :   provide at each time step the surface atm fluxes
      !!
      !!----------------------------------------------------------------------
  
      INTEGER, INTENT(in) ::   kt   ! time step (---> in CONCEPTS, this is the ocean time step
      !
      INTEGER  ::   ierror      ! return error code
      INTEGER  ::   ifpr,ji,jj  ! dummy loop indice
      INTEGER  ::   jfld        ! dummy loop arguments
      INTEGER  ::   ios         ! Local integer output status for namelist read
      REAL(wp) ::   xmask	! i DO NOT KNOW WHAT IS THIS : LANDMASK??
      character*512 GEM_rpn_list
      !                                         ! ====================== !
      IF( kt == 1 ) THEN                   !  First call kt=nit000  !
         !                                      ! ====================== !
 
!---------------------------------------------------------------------------------------------
         !                                         ! store namelist information in an array
         jfld = jpfld
         !

         ALLOCATE( sf(jfld), STAT=ierror )         ! set sf structure
         ALLOCATE( tprecip(jpi,jpj), STAT=ierror )
         ALLOCATE( sprecip(jpi,jpj), STAT=ierror )
         ALLOCATE( cat_i_rpn(jpi,jpj), STAT=ierror )
         ALLOCATE( tatm_ice(jpi,jpj), STAT=ierror )
         ALLOCATE( qatm_ice(jpi,jpj), STAT=ierror )
         ALLOCATE( zlm_i_rpn(jpi,jpj), STAT=ierror )
         ALLOCATE( zlt_i_rpn(jpi,jpj), STAT=ierror )
         ALLOCATE( wndi_ice(jpi,jpj), STAT=ierror )
         ALLOCATE( wndj_ice(jpi,jpj), STAT=ierror )
         ALLOCATE( qlw_ice(jpi,jpj,1), STAT=ierror )
         ALLOCATE( qsr_ice(jpi,jpj,1), STAT=ierror )
       

            DO ifpr= 1, jfld
              ALLOCATE( sf(ifpr)%fnow(jpi,jpj,1) )
              ALLOCATE( sf(ifpr)%fdta(jpi,jpj,1,2) )
            END DO
            
         ! First estimation of air density and ustar
           rhor_rpn(:,:) = rhor_ref
           ustr_rpn(:,:) = 0.
         !
      ENDIF

          CALL fld_read_rpn( kt, kt, sf, GEM_rpn_list)

      tprecip(:,:)     = sf(jp_prec)%fnow(:,:,1)  !Total precipitation
      sprecip(:,:)     = sf(jp_snow)%fnow(:,:,1)  !Snow precipitation
      cat_i_rpn(:,:)   = sf(jp_tair)%fnow(:,:,1)  !Potential temperature
      tatm_ice(:,:)    = sf(jp_tair)%fnow(:,:,1)* &
                        (sf(jp_prsg)%fnow(:,:,1)/ &
                         sf(jp_prst)%fnow(:,:,1))**(-CAPPA)
                                                  !Temperature
      qatm_ice(:,:)    = sf(jp_humi)%fnow(:,:,1)  !Specific humidity
      zlm_i_rpn(:,:)   = MAX(0.001,sf(jp_hghm)%fnow(:,:,1))
                                                  !Atmospheric forcing height (momentum)
                                                  !MAX to avoid problems if has unset halo rows
      zlt_i_rpn(:,:)   = MAX(0.001,sf(jp_hght)%fnow(:,:,1))
                                                  !Atmospheric forcing height (temperature)
                                                  !MAX to avoid problems if has unset halo rows

      wndi_ice(:,:)    = sf(jp_wndi)%fnow(:,:,1)  !i-wind component
      wndj_ice(:,:)    = sf(jp_wndj)%fnow(:,:,1)  !j-wind component
      qlw_ice(:,:,1)   = sf(jp_qlw )%fnow(:,:,1)  !Long wave down
      qsr_ice(:,:,1)   = sf(jp_qsr )%fnow(:,:,1)  !Short wave down

        print *, 'the data is loaded for timestep : ', kt 
   
     return
   END SUBROUTINE load_data
   
END MODULE gemdrv_load_forcing

