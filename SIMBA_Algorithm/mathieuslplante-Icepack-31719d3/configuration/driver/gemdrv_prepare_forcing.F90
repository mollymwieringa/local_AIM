MODULE gemdrv_prepare_forcing

   !!======================================================================
   !!                       ***  MODULE prepare GEM forcing  ***
   !!=====================================================================
   !! History :  1.0   !  10-18  created by M. Plante 
   !! 			            based on code in CONCEPTS : sbcice_cice.F90
   !!----------------------------------------------------------------------  
   !!----------------------------------------------------------------------
   !!   prepare_cice_sfc     : write the forcing fields in a cice compatible array
   !!   test_input_cice      : Make outputs of the array. 
   !!			       This is only to test the code
   !!   alloc_col_forcing    : create the array for the data at column location
   !!   prep_col_forcing     : write the data at column location to be sent to icepack
   !!----------------------------------------------------------------------
   !!======================================================================
   !!                       ***  MODULE  sbcice_cice  ***
   !! To couple with Icepack, the column component of CICE
   !!=====================================================================
   !!----------------------------------------------------------------------
   !!   sbc_ice_cice  : sea-ice model time-stepping and update ocean sbc over ice-covered area
   !!   
   !!   
   !!----------------------------------------------------------------------

   USE gemdrv_load_forcing     ! Surface boundary condition: CORE bulk

   IMPLICIT NONE
   PUBLIC


      REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::  &
         zlvl_rpn    , & ! momentum atm level height (m)
         zlvt_rpn    , & ! thermodynamic atm level height (m)
         uatm_rpn   , & ! wind velocity components (m/s)
         vatm_rpn    , &
         wind_rpn   , & ! wind speed (m/s)
         potT_rpn    , & ! air potential temperature  (K)
         Tair_rpn    , & ! air temperature  (K)
         Qa_rpn     , & ! specific humidity (kg/kg)
         rhoa_rpn    , & ! air density (kg/m^3)
         swvdr_rpn   , & ! sw down, visible, direct  (W/m^2)
         swvdf_rpn   , & ! sw down, visible, diffuse (W/m^2)
         swidr_rpn  , & ! sw down, near IR, direct  (W/m^2)
         swidf_rpn   , & ! sw down, near IR, diffuse (W/m^2)
         flw_rpn     , & ! incoming longwave radiation (W/m^2)   
         fsw_rpn         ! incoming longwave radiation (W/m^2)   
 
       real (kind=dbl_kind), dimension (:,:), allocatable, public :: &
         frain_rpn   , & ! rainfall rate (kg/m^2 s)
         fsnow_rpn       ! snowfall rate (kg/m^2 s)

         !!---------------------------------------------------------------------
      ! THIS CAN BE REPLACED BY USE ICE-FORCING WHEN IN CICE
       real (kind=dbl_kind), parameter, public :: &
         frcvdr = 0.28_dbl_kind, & ! frac of incoming sw in vis direct band
         frcvdf = 0.24_dbl_kind, & ! frac of incoming sw in vis diffuse band
         frcidr = 0.31_dbl_kind, & ! frac of incoming sw in near IR direct band
         frcidf = 0.17_dbl_kind    ! frac of incoming sw in near IR diffuse band     
      
         
         !!---------------------------------------------------------------------
      
      REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::  &     
         zlvl_col   , & ! momentum atm level height (m)
         zlvt_col   , & ! thermodynamic atm level height (m)
         uatm_col   , & ! wind velocity components (m/s)
         vatm_col   , &
         wind_col   , & ! wind speed (m/s)
         potT_col   , & ! air potential temperature  (K)
         Tair_col   , & ! air temperature  (K)
         Qa_col     , & ! specific humidity (kg/kg)
         rhoa_col   , & ! air density (kg/m^3)
         swvdr_col  , & ! sw down, visible, direct  (W/m^2)
         swvdf_col  , & ! sw down, visible, diffuse (W/m^2)
         swidr_col  , & ! sw down, near IR, direct  (W/m^2)
         swidf_col  , & ! sw down, near IR, diffuse (W/m^2)
         flw_col    , & ! incoming longwave radiation (W/m^2)   
         fsw_col    , & ! incoming longwave radiation (W/m^2) 
         frain_col  , & ! rainfall rate (kg/m^2 s)
         fsnow_col       ! snowfall rate (kg/m^2 s)
      !!---------------------------------------------------------------------   

   
CONTAINS


   SUBROUTINE prepare_cice_sfc (kt)
      !!---------------------------------------------------------------------
      !!                    ***  from CONCEPTS ROUTINE cice_sbc_in  ***
      !! ** Purpose: Set coupling fields and pass to CICE
      !!---------------------------------------------------------------------
      
      INTEGER, INTENT(in   ) ::   kt   ! ocean time step
      INTEGER  ::   ierror      ! return error code
      
         ALLOCATE( zlvl_rpn(jpi,jpj), STAT=ierror )         ! set sf structure
         ALLOCATE( zlvt_rpn(jpi,jpj), STAT=ierror )
         ALLOCATE( uatm_rpn(jpi,jpj), STAT=ierror )
         ALLOCATE( vatm_rpn(jpi,jpj), STAT=ierror )
         ALLOCATE( wind_rpn(jpi,jpj), STAT=ierror )
         ALLOCATE( potT_rpn(jpi,jpj), STAT=ierror )
         ALLOCATE( Tair_rpn(jpi,jpj), STAT=ierror )
         ALLOCATE( Qa_rpn(jpi,jpj), STAT=ierror )
         ALLOCATE( rhoa_rpn(jpi,jpj), STAT=ierror )
         ALLOCATE( swvdr_rpn(jpi,jpj), STAT=ierror )
         ALLOCATE( swvdf_rpn(jpi,jpj), STAT=ierror )
         ALLOCATE( swidr_rpn(jpi,jpj), STAT=ierror )
         ALLOCATE( swidf_rpn(jpi,jpj), STAT=ierror )
         ALLOCATE( flw_rpn(jpi,jpj), STAT=ierror )
         ALLOCATE( fsw_rpn(jpi,jpj), STAT=ierror )
         ALLOCATE( fsnow_rpn(jpi,jpj), STAT=ierror )
         ALLOCATE( frain_rpn(jpi,jpj), STAT=ierror )

         
! Pass CORE forcing fields to CICE (which will calculate heat fluxes etc itself)
         uatm_rpn(:,:) = wndi_ice(:,:)
         vatm_rpn(:,:) = wndj_ice(:,:)       
         wind_rpn(:,:) = SQRT ( wndi_ice(:,:)**2 + wndj_ice(:,:)**2 )   
         fsw_rpn(:,:) = qsr_ice(:,:,1)
         flw_rpn(:,:) = qlw_ice(:,:,1)
         Tair_rpn(:,:) = tatm_ice(:,:)
         potT_rpn(:,:) = cat_i_rpn(:,:)    ! Potential temp (K)
        
! Following line uses MAX(....) to avoid problems if tatm_ice has unset halo rows
         rhoa_rpn(:,:) = 101000. / ( 287.04 * MAX(1.0,tatm_ice(:,:)) ) ! Constant (101000.) atm pressure assumed
                                                     
         Qa_rpn(:,:) = qatm_ice(:,:)        
         zlvl_rpn(:,:)=zlm_i_rpn(:,:)      
         zlvt_rpn(:,:)=zlt_i_rpn(:,:) 
        
! Divide shortwave into spectral bands (as in prepare_forcing)
         swvdr_rpn(:,:)=qsr_ice(:,:,1)*frcvdr       ! visible direct        
         swvdf_rpn(:,:)=qsr_ice(:,:,1)*frcvdf       ! visible diffuse        
         swidr_rpn(:,:)=qsr_ice(:,:,1)*frcidr       ! near IR direct
         swidf_rpn(:,:)=qsr_ice(:,:,1)*frcidf       ! near IR diffuse

! Snowfall
! Ensure fsnow is positive (as in CICE routine prepare_forcing)
      fsnow_rpn(:,:)=MAX(sprecip(:,:),0.0)

! Rainfall
      frain_rpn(:,:)=(tprecip(:,:)-sprecip(:,:))
     
! sample print of the data, to test the code.     
!      print *, 'sample uatm : ', uatm_rpn(50,50)
!      print *, 'sample vatm : ', vatm_rpn(50,50)
!      print *, 'sample wind : ', wind_rpn(50,50)
!      print *, 'sample fsw : ', fsw_rpn(50,50)
!      print *, 'sample flw : ', flw_rpn(50,50)
!      print *, 'sample Tair : ', Tair_rpn(50,50)
!      print *, 'sample potT : ', potT_rpn(50,50)
!      print *, 'sample rhoa : ', rhoa_rpn(50,50)
!      print *, 'sample Qa : ', Qa_rpn(50,50)
      print *, 'sample zlvl : ', zlvl_rpn(50,50)
      print *, 'sample zlvs : ', zlvt_rpn(50,50)
!      print *, 'sample swvdr : ', swvdr_rpn(50,50)
!      print *, 'sample swvdf : ', swvdf_rpn(50,50)
!      print *, 'sample swidr : ', swidr_rpn(50,50)
!      print *, 'sample swidf : ', swidf_rpn(50,50)
!      print *, 'sample fsnow : ', fsnow_rpn(50,50)      
!      print *, 'sample frain : ', frain_rpn(50,50)
                      
      return
   END SUBROUTINE prepare_cice_sfc
   
   SUBROUTINE alloc_col_forcing(ktmax)
	implicit none

      INTEGER, INTENT(in   ) ::   ktmax   ! ocean time step
!      INTEGER, INTENT(in   ) ::   ksbc ! surface forcing type

      INTEGER  ::   ierror      ! return error code
  
	print *, 'Allocating the column forcing vectors'      


         ALLOCATE( zlvl_col(ktmax),  STAT=ierror )         ! set sf structure
         ALLOCATE( zlvt_col(ktmax),  STAT=ierror )
         ALLOCATE( uatm_col(ktmax),  STAT=ierror )
         ALLOCATE( vatm_col(ktmax),  STAT=ierror )
         ALLOCATE( wind_col(ktmax),  STAT=ierror )
         ALLOCATE( potT_col(ktmax),  STAT=ierror )
         ALLOCATE( Tair_col(ktmax),  STAT=ierror )
         ALLOCATE( Qa_col(ktmax),    STAT=ierror )
         ALLOCATE( rhoa_col(ktmax),  STAT=ierror )
         ALLOCATE( swvdr_col(ktmax), STAT=ierror )
         ALLOCATE( swvdf_col(ktmax), STAT=ierror )
         ALLOCATE( swidr_col(ktmax), STAT=ierror )
         ALLOCATE( swidf_col(ktmax), STAT=ierror )
         ALLOCATE( flw_col(ktmax),   STAT=ierror )
         ALLOCATE( fsw_col(ktmax),   STAT=ierror )
         ALLOCATE( fsnow_col(ktmax), STAT=ierror )
         ALLOCATE( frain_col(ktmax), STAT=ierror )

      return
   END SUBROUTINE alloc_col_forcing  
   
   
   
   SUBROUTINE prep_col_forcing(i_col,j_col,kt)
	implicit none

      INTEGER, INTENT(in   ) ::   i_col, j_col, kt   ! ocean time step

      INTEGER  ::   ierror      ! return error code

         uatm_col(kt) = uatm_rpn(i_col,j_col) 
         vatm_col(kt) = vatm_rpn(i_col,j_col)
         wind_col(kt) = wind_rpn(i_col,j_col)
         fsw_col(kt) = fsw_rpn(i_col,j_col) 
         flw_col(kt) = flw_rpn(i_col,j_col) 
         Tair_col(kt) = Tair_rpn(i_col,j_col) 
         rhoa_col(kt) = rhoa_rpn(i_col,j_col)
         Qa_col(kt) = Qa_rpn(i_col,j_col) 
         zlvl_col(kt) = zlvl_rpn(i_col,j_col)
         zlvt_col(kt) = zlvt_rpn(i_col,j_col)
         potT_col(kt) = potT_rpn(i_col,j_col) ! Potential temp (K)
         swvdr_col(kt) = swvdr_rpn(i_col,j_col)      ! visible direct
         swvdf_col(kt) = swvdf_rpn(i_col,j_col)     ! visible diffuse
         swidr_col(kt) = swidr_rpn(i_col,j_col)    ! near IR direct
         swidf_col(kt) = swidf_rpn(i_col,j_col)   ! near IR diffuse
         fsnow_col(kt) = fsnow_rpn(i_col,j_col)
         frain_col(kt) = frain_rpn(i_col,j_col)
         
!! sample print of the data, to test the code.           
!        print *, 'sample uatm : ', uatm_col(kt)
!        print *, 'sample vatm : ', vatm_col(kt)
!        print *, 'sample wind : ', wind_col(kt)
!        print *, 'sample fsw : ', fsw_col(kt)
!        print *, 'sample flw : ', flw_col(kt)
!        print *, 'sample Tair : ', Tair_col(kt)
!        print *, 'sample rhoa : ', rhoa_col(kt)
!        print *, 'sample Qa : ', Qa_col(kt)
        print *, 'sample zlvl : ', zlvl_col(kt)
        print *, 'sample zlvs : ', zlvt_col(kt)
!        print *, 'sample potT : ', potT_col(kt)
!        print *, 'sample swvdr : ', swvdr_col(kt)
!        print *, 'sample swvdf : ', swvdf_col(kt)
!        print *, 'sample swidr : ', swidr_col(kt)
!        print *, 'sample swidf : ', swidf_col(kt)
!        print *, 'sample fsnow : ', fsnow_col(kt)
!        print *, 'sample frain : ', frain_col(kt)
        
      return
   END SUBROUTINE prep_col_forcing   
   
   
   
   

   !!======================================================================
END MODULE gemdrv_prepare_forcing


