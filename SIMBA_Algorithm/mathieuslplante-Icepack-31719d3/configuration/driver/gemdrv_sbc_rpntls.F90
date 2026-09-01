module gemdrv_sbc_rpntls

!copyright (C) 2001  MSC-RPN COMM  %%%MC2%%%

   !!======================================================================
   !!         ***  from CONCEPTS MODULE  sbc_rpntls  ***
   !!=====================================================================
   !! History ICEPACK:   0.0  !  11-2018  (M. Plante, with subroutines from above) 
   !!----------------------------------------------------------------------

   IMPLICIT NONE
   PUBLIC

   INTEGER, PUBLIC, PARAMETER ::   dp = SELECTED_REAL_KIND(12,307)   !: double precision (real 8)
   INTEGER, PUBLIC, PARAMETER ::   wp = dp                              !: working precision   
   INTEGER, PUBLIC, PARAMETER ::   sp = SELECTED_REAL_KIND( 6, 37)   !: single precision (real 4)
   
   REAL (kind=8), PUBLIC :: tforc_1, &  ! Julian day (decimal) atm field 1
  &                 tforc_2     ! Julian day, atm field 2
   CHARACTER (LEN=16), PUBLIC :: current_atmf, & ! RPN formatted date (atm field no 2)
  &                              Mod_runstrt_S   ! RPN formatted date (simulation starting time)

   PUBLIC incdatsd, datp2f, datf2p, prsdate, pdfjdate2, inter_field3, &
          d_rawfstw, r_rawfstw, r_rawfstwm
          

CONTAINS

      SUBROUTINE incdatsd(newdate, olddate, dt)
      IMPLICIT NONE
      character(LEN=16) :: newdate,olddate
      real*8 :: dt

      real*8 ::jolddate,jnewdate

      integer :: newyy,newmo,newdd,newhh,newmm,newss
      integer :: oldyy,oldmo,olddd,oldhh,oldmm,oldss,oldsign      
!----------------------------------------------------------------------------      
      call prsdate(oldyy,oldmo,olddd,oldhh,oldmm,oldss,oldsign,olddate,1)

      call pdfjdate(jolddate,oldyy,oldmo,olddd,oldhh,oldmm,oldss)
      jnewdate=jolddate+dt

      call pdfcdate(newyy,newmo,newdd,newhh,newmm,newss,jnewdate)

      write(newdate,12) newyy,newmo,newdd,newhh,newmm,newss
 12   format(i4.4,i2.2,i2.2,'.',i2.2,i2.2,i2.2)
      return

      END SUBROUTINE incdatsd


      SUBROUTINE prsdate(yy,mo,dd,hh,mm,ss,sign,date,mode)
      IMPLICIT NONE

      integer :: yy,mo,dd,hh,mm,ss,sign,mode
      character(LEN=*) :: date
      character(LEN=16) :: tmpdate
      character(LEN=4) :: cyy
      character(LEN=2) :: cmo,cdd,chh,cmm,css
!----------------------------------------------------------------

      if (mode == 1) then

      if (date(1:1) == '-') then
         sign = -1
         tmpdate=date(2:16)
      else
         if (date(1:1) == ' ') then
            sign = 1
            tmpdate=date(2:16)
         else
            sign = 1
            tmpdate=date(1:15)
         endif
      endif

      cyy=tmpdate(1:4)
      cmo=tmpdate(5:6)
      cdd=tmpdate(7:8)
      chh=tmpdate(10:11)
      cmm=tmpdate(12:13)
      css=tmpdate(14:15)

      read(cyy,'(I4)') yy
      read(cmo,'(I2)') mo
      read(cdd,'(I2)') dd
      read(chh,'(I2)') hh
      read(cmm,'(I2)') mm
      read(css,'(I2)') ss

      elseif (mode == 2) then

      write(date,10) yy,mo,dd,hh,mm,ss
 10   format(i4.4,i2.2,i2.2,'.',i2.2,i2.2,i2.2)

      else

      write(date,12) yy,mo,dd,hh
 12   format(i4.4,i2.2,i2.2,i2.2,'_000')

      endif
      
      RETURN

      END SUBROUTINE prsdate

      SUBROUTINE pdfjdate2 (jdate,yyyy,mo,dd,hh,mm,ss)
      implicit none
      real*8 jdate
      integer yyyy,mo,dd,hh,mm,ss
!
!  calculate julian calendar day
!  see cacm letter to editor by fliegel and flandern 1968
!  page 657
!
      integer jd,jyy,jmo,jdd
      real*8 one,sec_in_day
      parameter (one=1.0d0, sec_in_day=one/86400.0d0)

      jd(jyy,jmo,jdd)=jdd-32075+1461*(jyy+4800+(jmo-14)/12)/4     &
     &     +  367*(jmo-2-(jmo-14)/12*12)/12 - 3     &
     &     *((jyy+4900+(jmo-14)/12)/100)/4


      jdate = jd(yyyy,mo,dd) - 2433646 ! good from 1951 onwards
      if (jdate.lt.0) then
         print*, 'Negative Julian day in pdfjdate2 --- ABORT ---'
         stop
      endif

      jdate = jdate + dble(hh*3600+mm*60+ss)*sec_in_day

      return

      END SUBROUTINE pdfjdate2

      SUBROUTINE pdfjdate(jdate,yyyy,mo,dd,hh,mm,ss)
      IMPLICIT NONE

      real*8 :: jdate
      integer ::  yyyy,mo,dd,hh,mm,ss
!
!  calculate julian calendar day
!  see cacm letter to editor by fliegel and flandern 1968
!  page 657
!
      integer jd,jyy,jmo,jdd

! Old fortran: a  statement function
      jd(jyy,jmo,jdd)=jdd-32075+1461*(jyy+4800+(jmo-14)/12)/4    &
     &     +  367*(jmo-2-(jmo-14)/12*12)/12 - 3         &
     &     *((jyy+4900+(jmo-14)/12)/100)/4
!------------------------------------------------------------------------

      jdate = jd(yyyy,mo,dd)
      jdate = jdate + (hh*3600+mm*60+ss)/86400.0
      return

      END SUBROUTINE pdfjdate

      SUBROUTINE pdfcdate(yyyy,mo,dd,hh,mm,ss,jdate)
      IMPLICIT NONE
      real*8 jdate
      integer yyyy,mo,dd,hh,mm,ss,seconds

      real*8 :: f,rj
!--------------------------------------------------------------------

      rj = int(jdate)
      f = jdate - rj
      seconds = nint(f * 86400.0)
      
      ss = mod(seconds, 60)
      mm = mod(seconds - ss,3600)/60
      
      
      hh = (seconds-60*mm-ss) / 3600
      if (hh.eq.24) then
         hh = 0
         seconds = seconds - 86400
         rj = rj+1.0
      endif
      mm = (seconds - hh * 3600 - ss) / 60
      
      call datec(int(rj),yyyy,mo,dd)
      
      return

      END SUBROUTINE pdfcdate

      SUBROUTINE datp2f (fstdate,mc2date)
      IMPLICIT NONE

      integer :: fstdate
      character(LEN=*) :: mc2date
      integer :: yy,mo,dd,hh,mm,ss,dat2,dat3,newdate,err
      character(LEN=4) :: cyy
      character(LEN=2) :: cmo,cdd,chh,cmm,css
!-----------------------------------------------------------
      cyy=mc2date(1:4)
      cmo=mc2date(5:6)
      cdd=mc2date(7:8)
      chh=mc2date(10:11)
      cmm=mc2date(12:13)
      css=mc2date(14:15)

      read(cyy,'(I4)') yy
      read(cmo,'(I2)') mo
      read(cdd,'(I2)') dd
      read(chh,'(I2)') hh
      read(cmm,'(I2)') mm
      read(css,'(I2)') ss

      dat2= yy*10000 + mo*100 + dd
      dat3= hh*1000000 + mm*10000 + ss*100
      err = newdate(fstdate,dat2,dat3,3)

      RETURN
      END SUBROUTINE datp2f

      subroutine datf2p (mc2date,fstdate)
      implicit none
!
      character(LEN=*) :: mc2date
      integer :: fstdate
!
!ARGUMENTS 
!     NAMES     I/O  TYPE  A/S DESCRIPTION
!
!     mc2date    O     C    S  date encoded in mc2 format
!     fstdate    I     I    S  date encoded in RPN standard file format
!
!MODULES 
!
!
      integer :: yy,mo,dd,hh,mm,ss
      integer :: dat2,dat3,newdate,err
!     
      err= newdate(fstdate,dat2,dat3,-3)
!
      yy = dat2/10000
      mo = mod(dat2,10000)/100
      dd = mod(dat2,100)
      hh = dat3/1000000
      mm = mod(dat3,1000000)/10000
      ss = mod(dat3,10000)/100
!
      write(mc2date(1:16),10) yy,mo,dd,hh,mm,ss
 10   format(i4.2,i2.2,i2.2,'.',i2.2,i2.2,i2.2,' ')
!
      return
      END SUBROUTINE datf2p

      SUBROUTINE inter_field3 (b,fdta1,fdta2,fnow,ni,nj)
      implicit none
      real(wp), intent(in)                      :: b
      integer,  intent(in)                      :: ni,nj
      real(wp), INTENT(in),    DIMENSION(ni,nj) :: fdta1, fdta2
      real(wp), INTENT(inout), DIMENSION(ni,nj) :: fnow

      real(wp)                                  :: fac1

      fac1 = 1. - b
      fnow(:,:) = fac1*fdta1(:,:) + b*fdta2(:,:)
      return

      END SUBROUTINE inter_field3


      subroutine d_rawfstw (rs,nx,ny,varname,stepno,dt,date0,out_file)

      implicit none
!
      integer :: nx,ny,stepno,date0
      real(kind=wp) :: dt
      real(kind=wp) :: rs(nx,ny)
      character(LEN=*) :: varname,out_file
!
      integer fnom,fclos
      integer id_unit,err,ip2
!
      real(kind=sp) :: rl_work(nx,ny)
!----------------------------------------------------------------------
      id_unit=0
      err = fnom (id_unit, out_file, 'rnd', 0)
      call fstouv (id_unit,'rnd')
!
      ip2 =(stepno*dt)/3600
      rl_work=REAL(rs,sp)   ! conversion to single precision 
      call fstecr (rl_work,rl_work,-32,id_unit, date0, NINT(max(dt,3600.0)), stepno, nx,ny,1, &
                &  0,ip2,0,'P',varname,'rslt2','X',1,1,1,1,1,.false.)

      call fstfrm (id_unit)
      err = fclos (id_unit)
!
      return
      end subroutine d_rawfstw


      subroutine r_rawfstw (rs,nx,ny,varname,stepno,dt,date0,out_file)

      implicit none
!
      integer :: nx,ny,stepno,date0
      real(kind=wp) :: dt
      real(kind=sp) :: rs(nx,ny)
      character(LEN=*) :: varname,out_file
!
      integer fnom,fclos
      integer id_unit,err,ip2
!----------------------------------------------------------------------
      id_unit=0
      err = fnom (id_unit, out_file, 'rnd', 0)
      call fstouv (id_unit,'rnd')
!
      ip2 =(stepno*3600.0)/3600
      call fstecr (rs,rs,-32,id_unit, date0, NINT(dt), stepno, nx,ny,1, &
                &  0,ip2,0,'P',varname,'rslt2','X',1,1,1,1,1,.false.)

      call fstfrm (id_unit)
      err = fclos (id_unit)
!
      return
      end subroutine r_rawfstw

      subroutine r_rawfstwm (rs,ms,nx,ny,varname,stepno,dt,date0,out_file)

      implicit none
!
      integer :: nx,ny,stepno,date0
      real(kind=wp) :: dt
      real(kind=sp) :: rs(nx,ny),ms(nx,ny)
      character(LEN=*) :: varname,out_file
!
      integer fnom,fclos
      integer id_unit,err,ip2,i,j
      integer(kind=sp):: wms(nx,ny),nmean
      real(kind=sp) ::   wrs(nx,ny), rmean

!----------------------------------------------------------------------

      rmean=0.
      nmean=0
      do j=1,ny
      do i=1,nx
        wms(i,j)=nint(ms(i,j))
        if (wms(i,j).eq.1) then
           rmean=rmean+rs(i,j)
           nmean=nmean+1
        endif
      enddo
      enddo
      if (nmean.gt.0) rmean=rmean/nmean
      do j=1,ny
      do i=1,nx
        wrs(i,j)=rs(i,j)
        if (wms(i,j).ne.1) wrs(i,j)=rmean
      enddo
      enddo

      id_unit=0  
      err = fnom (id_unit, out_file, 'rnd', 0)
      call fstouv (id_unit,'rnd')
!
      ip2 =(stepno*3600.0)/3600
      call fstecr (wrs,wrs,-32,id_unit, date0, NINT(dt), stepno, nx,ny,1, &
                &  0,ip2,0,'P@',varname,'rslt2','X',1,1,1,1,1,.false.)
      call fstecr (wms,wms,-1,id_unit, date0, NINT(dt), stepno, nx,ny,1, &
                &  0,ip2,0,'@@',varname,'rslt2','X',1,1,1,1,2,.false.)

      call fstfrm (id_unit)
      err = fclos (id_unit)
!
      return
      end subroutine r_rawfstwm


end module gemdrv_sbc_rpntls
