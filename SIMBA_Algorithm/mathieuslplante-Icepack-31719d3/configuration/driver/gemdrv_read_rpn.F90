MODULE gemdrv_read_rpn
   !!======================================================================
   !!         ***  from CONCEPTS MODULE  fldread_rpn  ***
   !! Read input field for surface boundary condition 
   !!                 (rpn files)
   !!=====================================================================
   !! History CONCEPTS:  9.0  !  08-07  (F. Roy with original subroutines 
   !!                                           provided by  M. Desgagne) 
   !!  	      ICEPACK:   0.0  !  11-2018  (M. Plante, with subroutines from above) 
   !!----------------------------------------------------------------------
   !!   fld_read_rpn : 
   !!                   read input fields used for the computation of the
   !!                   surface boundary condition (written in rpn files)
   !!----------------------------------------------------------------------

    USE icedrv_kinds  
    USE icedrv_calendar
    USE gemdrv_sbc_rpntls

   IMPLICIT NONE
   PUBLIC   

   !! * Routine accessibility
   
   INTEGER :: iun_std(1000), &  ! Unit number vector
  &           nfile_std, &      ! Number of std files
  &           kt_sbc            ! Time step used for atm. data interpolation, or treatment of coupling 
   LOGICAL :: atm_file_L

   REAL(wp) ::   dt_gem_atm !: spacing between forcing fields (GEM) (hours at this stage)
   REAL(wp) ::   dt_gem_prc !: precip dt used for cumulation in GEM (hours at this stage)
                            !: equal to dt_gem_atm if no average in pre-process
                    
   LOGICAL  ::   ln_gem_intrp = .TRUE.   !: switch to activate linear time-interpolation
   INTEGER  ::   nh_gem_offs = 0          !: offset to match dates when searching forcing fields
                    !: example, if daily averages correspond to hour 3, nh_gem_offs = 3
                    !: must be positive
   LOGICAL  ::   ln_gem_avgrad = .FALSE.   !: switch to activate 24h averaging for radiation
   INTEGER, PUBLIC, DIMENSION(2) :: lev_gem_forc=(/ 28257976, 11950 /) 
                    !: Vertical level descriptor (equivalent to ip1), the first
                    !: one is the newstyle version, and the second one the
                    !: oldstyle. It allows to mix oldstyle and newstyle code
                    !: in atmospheric forcing std files.
                    !: The oldstyle code (second value) will only be used
                    !: if fields coded in newstyle are not found.
                    !: The default here is eta=0.995 
   TYPE, PUBLIC ::   FLD        !: Input field related variables
      REAL(wp) , ALLOCATABLE, DIMENSION(:,:,:) ::   fnow       ! input fields interpolated to now time step
      REAL(wp) , ALLOCATABLE, DIMENSION(:,:,:,:) ::   fdta       ! 2 consecutive record of input fields
   END TYPE FLD  
   
   REAL(wp), PUBLIC, SAVE ::   &
                 zlev_gem   = -9.       !: optional fixed forcing level from namelist
                                        !: allow to skip reading of PX,P0
                                        !: do nothing if left to -9.
  

   INTEGER, PUBLIC  ::   jpi = 528   ! = ( jpiglo-2*jpreci + (jpni-1) ) / jpni + 2*jpreci   !: first  dimension
   INTEGER, PUBLIC  ::   jpj = 735  ! = ( jpjglo-2*jprecj + (jpnj-1) ) / jpnj + 2*jprecj   !: second dimension  
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: lat_rpn
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: lon_rpn   
   
   REAL(wp), PUBLIC, PARAMETER ::   RGASV  =.46151e+3      ! J K-1 kg-1; gas constant 
                                                           ! for water vapour
   REAL(wp), PUBLIC, PARAMETER ::   CAPPA  =.28549121795   ! RGASD/CPD         ! ; Von Karman constant
   REAL(wp), PUBLIC, PARAMETER ::   RGASD  =.28705e+3    
   REAL(wp), PUBLIC, PARAMETER ::   DELTA  =.6077686814144 ! ; 1/EPS1 - 1
   REAL(wp), PUBLIC            ::   grav  = 9.80665_wp     !: gravity                            [m/s2]   
   
   INTEGER ::   numout          =    6      !: logical unit for output print; Set to stdout to ensure any early
                                            !  output can be collected; do not change
                                            
   INTEGER       ::   nn_date0_hour=6  !: initial hour of the calendar da 
   
   PUBLIC   fld_read_rpn   ! called by sbc... modules

CONTAINS

   SUBROUTINE fld_read_rpn( kt, kn_fsbc, sd, GEM_rpn_list )
      implicit none
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE fld_read_rpn  ***
      !!                   
      !! ** Purpose :   provide at each time step atmospheric variables
      !!                needed by sbc routine (core)
      !!                (UU,VV,TT, etc.) 
      !!
      !! ** Method  :   READ each input fields in STD files using RPN libs
      !!      and intepolate it to the model time-step.
      !!      Forcing files must be listed in "liste_inputfiles",located in
      !!      execution directory
      !!----------------------------------------------------------------------
      INTEGER  , INTENT(in   )               ::   kt        ! ocean time step
      INTEGER  , INTENT(in   )               ::   kn_fsbc   ! sbc computation period (in time step) 
      
   
      TYPE(FLD), INTENT(inout), DIMENSION(:) ::   sd        ! input field related variables
     
      !!
      INTEGER  ::   jf         ! dummy indices
      INTEGER  ::   ios        ! namelist error control
      
      character*512 GEM_rpn_list
      !!---------------------------------------------------------------------

      IF( kt == 1 ) THEN
      
        DO jf = 1, SIZE( sd )                     !    LOOP OVER FIELD    !
          sd(jf)%fnow(:,:,:)   = 0.0
          sd(jf)%fdta(:,:,:,:) = 0.0 ! tricky but working
        ENDDO
        
        kt_sbc=kt-1	
        call init_atm (sd, GEM_rpn_list) ! get the initial forcing for step 0
        
      ENDIF

      kt_sbc=kt
      call linear_tint (sd) ! get the forcing for current time

       return
      END SUBROUTINE fld_read_rpn

      
      
      SUBROUTINE init_atm (sd, GEM_rpn_list)
      implicit none
      
      TYPE(FLD), INTENT(inout), DIMENSION(:) ::   sd
      integer yy,mo,dd,hh,mm,ss,dum
      character*16 datev
      character*512 GEM_rpn_list
      real*8  dayfrac,one,sid,rsid
      parameter(one=1.0d0, sid=86400.0d0, rsid=one/sid)

     !getting the initial date (Mod_runstrt_S) :
        write(Mod_runstrt_S(1:8),'(i8.8)') idate0
        Mod_runstrt_S(9:9)='.'
        write(Mod_runstrt_S(10:11),'(i2.2)') nh_gem_offs + nn_date0_hour
        Mod_runstrt_S(12:16)='0000 '  !Already assumed in opa
        print *, 'we get the Mod_runstrt_S : ', Mod_runstrt_S
        print *, 'Going into atmopenf, with ksbc : ', kt_sbc

      call atm_openf(GEM_rpn_list)
      dayfrac = dble(kt_sbc)*max(dt,3600.0)*rsid ! current timestep into seconds from sim start      
      call incdatsd  (datev,Mod_runstrt_S,dayfrac) !getting the date of current time
      call prsdate   (yy,mo,dd,hh,mm,ss,dum,datev,1) !from date string to yy,mm,etc
      call pdfjdate2 (tforc_2,yy,mo,dd,hh,mm,ss) ! get current time (sec) into tforc_2
      print *, 'Going into (initial) atm_get_data, datev, dayfrac : ', datev, dayfrac, dt
      
      call atm_getdata (datev,.false.,.true.,sd) ! get the date for current time (datev)
       print *, 'end atm_getdata initial'
      return
      END SUBROUTINE init_atm
      
      
   
      SUBROUTINE atm_openf(GEM_rpn_list)
      implicit none
      
      integer maxnfile
      parameter ( maxnfile=1000 )
      character*512 filename(maxnfile),fn,pwd,GEM_rpn_list
      integer  fnom,fstouv,fstlnk
      external fnom,fstouv,fstlnk
      integer err,err1,err2,i,cnt,unf,wkoffit
      nfile_std = 0
      cnt       = 0
      unf       = 0
	call GETCWD(pwd)
	print *, trim(pwd)

      if (fnom(unf,GEM_rpn_list,'SEQ+OLD',0).lt.0) then
       call abort('atm_openf: GEM_atm_forcing.txt, error fnom')
      endif
 77   cnt=cnt+1
      if (cnt.gt.maxnfile) then 
        call abort('atm_openf: maxnfile not large enough ---ABORT---')
      endif
      read (unf, '(a)', end = 9120) filename(cnt)

      goto 77
 9120 nfile_std = cnt - 1
      close(unf)
      print *, nfile_std, ' files'
      do cnt = 1, nfile_std
         fn  = filename(cnt)
         print *, fn, cnt
         err = wkoffit(fn)

         if ((err.ne.1).and.(err.ne.33)) then
            filename(cnt) = '@#$%^&'
         endif
      end do
      i=0
      do cnt = 1, nfile_std
         if (filename(cnt).ne.'@#$%^&') then
            i = i+1
            filename(i) = filename(cnt)
         endif
      end do
      nfile_std = i
      if (nfile_std.lt.1) then 
        call abort('atm_openf: NO STD FILE available ---ABORT---')
      endif
      do cnt = 1, nfile_std
         iun_std(cnt) = 0
         err1 = FNOM  (iun_std(cnt),trim(filename(cnt)),'RND+OLD',0)
         err2 = FSTOUV(iun_std(cnt),'RND')
         if ((err1.lt.0).or.(err2.lt.0)) &
           call abort('atm_openf: err1... ---ABORT---')
      end do
      err = fstlnk (iun_std,nfile_std)
      atm_file_L = .true.
 1001 format ('Opening STD file: ',a,2x,'UNIT= ',i4)
 1002 format ('Opening BINARY file: ',a,2x,'UNIT= ',i4)
      return

      END SUBROUTINE atm_openf     
         
      
      

      SUBROUTINE linear_tint (sd)
      implicit none
      TYPE(FLD), INTENT(inout), DIMENSION(:) ::   sd     
      INTEGER  ::   jf         ! dummy indices

      character*16 datev,daten
      integer yy,mo,dd,hh,mm,ss,dum,datm
      real*8  tx,dayfrac,one,sid,rsid
      real(wp) ::  b
      parameter(one=1.0d0, sid=86400.0d0, rsid=one/sid)
          
      dayfrac = dble(kt_sbc)*max(dt,3600.0)*rsid 
      call incdatsd  (datev,Mod_runstrt_S,dayfrac)
      print *, 'datev, Mod_runstrt_S, dayfrac :', datev, Mod_runstrt_S, dayfrac
      if (datev.gt.current_atmf) then
      print *, 'Changing the forcing data'
      
      call prsdate   (yy,mo,dd,hh,mm,ss,dum,datev,1)      
      if (yy < 2011)  then 
 	  dt_gem_atm = 3600.
      else 
           dt_gem_atm = 3600.*3.0
      endif      
         dayfrac = dble(dt_gem_atm)*rsid 
         call incdatsd (daten,current_atmf,dayfrac)
           print *, 'Going into atm_get_data, datev : ', datev     
         call atm_getdata (daten,.true.,.true.,sd)
      endif
      call prsdate   (yy,mo,dd,hh,mm,ss,dum,datev,1)
      call pdfjdate2 (tx,yy,mo,dd,hh,mm,ss)
      
      if (ln_gem_intrp) then
        b=(tx-tforc_1)/(tforc_2-tforc_1)
! Using averaged fields with date stamp corresponding to the end of
! the average serie
      else
        b=1.
      endif

      print *, 'interpolating'
      DO jf = 1, SIZE( sd )
        call inter_field3 (b,sd(jf)%fdta(:,:,1,1), &
                             sd(jf)%fdta(:,:,1,2), &
                             sd(jf)%fnow(:,:,1),jpi,jpj)
      ENDDO

 1001 format (/' ----- ABORT ----- No ATM data valid at: ',a)
      return

      END SUBROUTINE linear_tint



 
      
      
      
      
      SUBROUTINE atm_getdata (datev,put,get,sd)
      implicit none
      character*(*) datev
      logical put,get
      TYPE(FLD), INTENT(inout), DIMENSION(:) ::   sd
      integer yy,mo,dd,hh,mm,ss,datm,dum,jf
      !print *, 'switching tforc : ', tforc_1, tforc_2
      tforc_1        = tforc_2
      !print *, 'switching tforc : ', tforc_1, tforc_2
      !print *, 'switching data T sample: ', sd(6)%fdta(50,50,1,1), sd(6)%fdta(50,50,1,2)   
      DO jf = 1, SIZE( sd )     
        sd(jf)%fdta(:,:,1,1) = sd(jf)%fdta(:,:,1,2)
      ENDDO
      !print *, 'Going into atm_read'
      call atm_read_rpn (datev,sd)
      current_atmf = datev
      call prsdate   (yy,mo,dd,hh,mm,ss,dum,datev,1)
      call pdfjdate2 (tforc_2,yy,mo,dd,hh,mm,ss)
      !print *, 'switching data T sample: ', sd(6)%fdta(50,50,1,1), sd(6)%fdta(50,50,1,2)
      return
      END SUBROUTINE atm_getdata

      SUBROUTINE atm_read_rpn (datev,sd)
      implicit none
      character*(*) datev
      character*16 datev2
      TYPE(FLD), INTENT(inout), DIMENSION(:) ::   sd
      integer datm, dum
      real KNAMS,kprec
      parameter (KNAMS=0.514791)
      integer i,j,iavg,navg,nivr,nivt,nivm
      integer, dimension(2) :: lev_nul=(/-1,-1/), lev_wrk
      integer yy,mo,dd,hh,mm,ss
    
      real*8  dayfrac,one,sid,rsid
      character*16 datew
      parameter(one=1.0d0, sid=86400.0d0, rsid=one/sid)
      real(wp) theta,pi,Tv


      call prsdate   (yy,mo,dd,hh,mm,ss,dum,datev,1)      
      if (yy < 2011)  then 
 	  dt_gem_prc = 3600.
      else 
           dt_gem_prc = 3600.*3.0
      endif    
      kprec=1000./dt_gem_prc
      print *, 'starting to write the rpn into arrays'
! CAREFUL Must be structured based on sbcblk_rpn.F90 indexing
! Assuming that winds and HU,TT are at same height, will need adaptation
! when using GEM staggered coordinates
      call readfstatm (  sd(1)%fdta (:,:,1,2),jpi,jpj,'UUOR',datev,lev_gem_forc,KNAMS,0.,nivr )
	print *, 'UUOR sample : ', sd(1)%fdta(50,50,1,2)
      call readfstatm (  sd(2)%fdta (:,:,1,2),jpi,jpj,'VVOR',datev,lev_gem_forc,KNAMS,0.,nivr )
	print *, 'VUOR sample : ', sd(2)%fdta(50,50,1,2)
	
      nivm=nivr
      
      call readfstatm (  sd(3)%fdta (:,:,1,2),jpi,jpj,'HU'  ,datev,lev_gem_forc,1.0  ,0.,nivr )
	print *, 'HU sample : ', sd(2)%fdta(50,50,1,2)
	
      nivt=nivr

        write(numout,*) 'atm_read_rpn: thermodynamic atm level: ', nivt
        write(numout,*) 'atm_read_rpn: momentum      atm level: ', nivm

      call readfstatm (sd(4)%fdta (:,:,1,2),jpi,jpj,'FB'  ,datev,lev_nul, 1.0  ,0.,nivr )
	print *, 'FB sample : ', sd(4)%fdta(50,50,1,2)        
      call readfstatm (sd(5)%fdta (:,:,1,2),jpi,jpj,'FI'  ,datev,lev_nul, 1.0  ,0.,nivr )
	print *, 'FI sample : ', sd(4)%fdta(50,50,1,2)        
      call readfstatm (  sd(6)%fdta (:,:,1,2),jpi,jpj,'TT'  ,datev,lev_gem_forc,1.0  ,273.16,nivr )
	print *, 'TT sample : ', sd(6)%fdta(50,50,1,2) 
	
      if (nivr /= nivt) call abort('atm_read_rpn: TT and HU should have same ip1 -- WILL ABORT')

      call readfstatm (  sd(7)%fdta (:,:,1,2),jpi,jpj,'PR'  ,datev,lev_nul, kprec,0.,nivr )
	print *, 'PR sample : ', sd(7)%fdta(50,50,1,2) 
        print *, 'zlev_gem :', zlev_gem
      if (zlev_gem.lt.0.) then
        lev_wrk(:) = nivt
        call readfstatm (sd(11)%fdta(:,:,1,2),jpi,jpj,'PX'  ,datev,lev_wrk, 100.0,0.,nivr )
	  print *, 'PX sample : ', sd(11)%fdta(50,50,1,2)
	  
        if ( nivm /= nivt ) then ! staggered case
          lev_wrk(:) = nivm
          call readfstatm (sd(12)%fdta(:,:,1,2),jpi,jpj,'PX'  ,datev,lev_wrk, 100.0,0.,nivr )
	    print *, 'PX-2 sample : ', sd(12)%fdta(50,50,1,2)
        else
          sd(12)%fdta(:,:,1,2)=sd(11)%fdta(:,:,1,2)
	    print *, 'PX sample-double : ', sd(12)%fdta(50,50,1,2)
        endif
        
        call readfstatm (sd(13)%fdta(:,:,1,2),jpi,jpj,'P0'  ,datev,lev_nul, 100.0,0.,nivr )
	  print *, 'P0 sample : ', sd(13)%fdta(50,50,1,2)
        !        if(ln_apr_dyn) &
!        call readfstatm (sd(15)%fdta(:,:,1,2),jpi,jpj,'PN'  ,datev,lev_nul, 100.0,0.,nivr )
!   ... Thickness of first atm. model layer
!   ... (hydrostatic relation approximating Tv_bar=Tv(layer=1))
        do j=1,jpj
        do i=1,jpi
          Tv= sd(6)%fdta(i,j,1,2) * (1.0D0 + DELTA*sd(3)%fdta(i,j,1,2))
          sd(9 )%fdta(i,j,1,2) = -RGASD*Tv/grav*log(sd(11)%fdta(i,j,1,2)/sd(13)%fdta(i,j,1,2))
          sd(10)%fdta(i,j,1,2) = -RGASD*Tv/grav*log(sd(12)%fdta(i,j,1,2)/sd(13)%fdta(i,j,1,2))
        enddo
        enddo
       else
        sd(9 )%fdta(:,:,1,2) = zlev_gem
        sd(10)%fdta(:,:,1,2) = zlev_gem
        sd(11)%fdta(:,:,1,2) = 100000.
        sd(12)%fdta(:,:,1,2) = 100000.
        sd(13)%fdta(:,:,1,2) = 100000.
        sd(15)%fdta(:,:,1,2) = 100000.
      endif


      do j=1,jpj 
      do i=1,jpi 
        sd(8)%fdta(i,j,1,2)=0.
        if (sd(6)%fdta(i,j,1,2).le.273.16) sd(8)%fdta(i,j,1,2)=sd(7)%fdta(i,j,1,2)
                ! Now using upper level, 8 corresponds to snow
      enddo
      enddo

! Converts sd(6)%fdta (:,:,1,2) to potential temperature
      sd(6)%fdta(:,:,1,2) = sd(6)%fdta(:,:,1,2)*    &
     &                (sd(13)%fdta(:,:,1,2)/sd(11)%fdta(:,:,1,2))**CAPPA
     
     
!-----extraction of the lat and lon arrays      
      if (datev .eq. Mod_runstrt_S) then
            call prsdate   (2001,01,10,00,00,00,dum,datev2,2)
      	  print *, 'DATEV2: ', datev2
          allocate(lat_rpn(jpi,jpj))
          allocate(lon_rpn(jpi,jpj))        
          call readfstatm ( lat_rpn(:,:),jpi,jpj,'^^',datev2,lev_gem_forc,1.,0.,nivr )
	  print *, 'LAT sample : ', lat_rpn(50,50)   
          call readfstatm ( lon_rpn(:,:),jpi,jpj,'>>',datev2,lev_gem_forc,1.,0.,nivr )
	  print *, 'LON sample : ', lon_rpn(50,50)     
      endif
 !-------------       
	print *, 'END get DATA!!!!'
	
      return
      
      END SUBROUTINE atm_read_rpn

      
      
      
      SUBROUTINE readfstatm (f,jpi,jpj,nomvar,dat,niv,factm,facta,nivr)


      implicit none

      character*(*) nomvar,dat
      integer jpi,jpj
      INTEGER niv(2)
      real factm,facta
      REAL(wp) , DIMENSION(jpi,jpj) :: f
      INTEGER, INTENT(out) :: nivr
      INTEGER fstinl,fstluk,ip1_all,errout

      integer ni1,nj1,nk1,err,nlis,lislon,i,j,datm,mi,mj
      parameter (nlis = 1024)
      integer liste (nlis)
      real*4 wrk(jpi,jpj)
      
	print *, 'In readfstatm for : ', nomvar
      nivr=0
      errout=0
        call datp2f (datm,dat)
        if (niv(1).eq.-1) then
          err = fstinl (iun_std,ni1,nj1,nk1,datm,' ',niv(1),-1,-1,' ',nomvar,  &
     &                                                   liste,lislon,nlis)
          nivr=niv(1)
         else
          err = fstinl (iun_std,ni1,nj1,nk1,datm,' ',niv(1),-1,-1,' ',nomvar,  &
     &                                                   liste,lislon,nlis)
          nivr=niv(1)
          if (lislon.ne.1) then
          err = fstinl (iun_std,ni1,nj1,nk1,datm,' ',niv(2),-1,-1,' ',nomvar,  &
     &                                                   liste,lislon,nlis)
          nivr=niv(2)
          endif
        endif
        
        if (nomvar.eq.'>>') then
	 print *, 'datm: ', datm
          err = fstinl (iun_std,ni1,nj1,nk1,datm,' ',1001,-1,-1,' ',nomvar,  &
     &                                                   liste,lislon,nlis)
          nivr=1001
          err = fstluk (wrk,liste(1),ni1,nj1,nk1)
		print *, 'lon_sample : ', wrk(50,50)
        elseif (nomvar.eq.'^^') then
          err = fstinl (iun_std,ni1,nj1,nk1,datm,' ',1001,-1,-1,' ',nomvar,  &
     &                                                   liste,lislon,nlis)
          nivr=1001
          err = fstluk (wrk,liste(1),ni1,nj1,nk1)
        endif       
        
        
        if (lislon.ne.1) then
            print*, 'Variable ',nomvar,' valid: ',dat,  &
     &              ' NOT found ---WILL ABORT---'
            errout=-1
        else
          err = fstluk (wrk,liste(1),ni1,nj1,nk1)
        endif

	print *, 'writing wrk into f'
	print *, 'jpi, jpj', jpi, jpj 
	do i = 1,jpi
          do j = 1,jpj
            f(i,j) = wrk(i,j)
          enddo
        enddo
        
        f   = f*factm + facta

      return
      END SUBROUTINE readfstatm
      
      
      SUBROUTINE get_col_position (lat_col,lon_col,i_col,j_col)
      implicit none

      real (kind=dbl_kind), intent(in):: lat_col, lon_col
      INTEGER, INTENT(out) :: i_col, j_col     
      INTEGER, DIMENSION(2) :: ind_col 
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: wrk_col   
      
      allocate(wrk_col(jpi,jpj))   
      
      wrk_col = ((lat_rpn - lat_col)**2d0 + (lon_rpn - lon_col)**2d0)**0.5d0
      ind_col = minloc(wrk_col)
      i_col = ind_col(1)
      j_col = ind_col(2)
      print *, 'POSITION: ', lon_col, lat_col, lat_rpn(i_col,j_col), lon_rpn(i_col,j_col)
      print *, 'POSITION2: ', lon_col, lat_col, lat_rpn(j_col,i_col), lon_rpn(j_col,i_col)  
      return
      END SUBROUTINE get_col_position      
      


   !!==============================================================================

END MODULE gemdrv_read_rpn
