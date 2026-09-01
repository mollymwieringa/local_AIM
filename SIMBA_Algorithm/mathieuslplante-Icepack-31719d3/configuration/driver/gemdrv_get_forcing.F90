   !!======================================================================
   !!                       ***  MODULE  get GEM forcing  ***
   !!=====================================================================
   !! History :  1.0   !  10-18  created by M. Plante 
   !! 			            based on code in CONCEPTS
   !!----------------------------------------------------------------------  
   !! Gets the forcing data from GEM, at a specific location, to force Icepack
   !!
   !!           If we initialise with buoy data, the location is taken 
   !!	        from the buoy data. Otherwise, the location is hardcoded
   !!           in icedrv_init.F90, to a location close to Nain.
   !!
   !!----------------------------------------------------------------------
   !!   prepare_gem_forcing     : Compute the atmospheric fluxes over ice
   !!----------------------------------------------------------------------



module gemdrv_get_forcing

     use gemdrv_load_forcing
     use gemdrv_prepare_forcing
     use icedrv_calendar

     
      implicit none

      integer :: kt 
      integer :: i_col, j_col


 contains      
      
      SUBROUTINE prepare_gem_forcing(idate0,npt,lat_col, lon_col,GEM_rpn_list)
           integer, intent(in) :: idate0, npt 
           real (kind=dbl_kind), intent(in) :: lat_col, lon_col	   
           character*512 GEM_rpn_list
          
      call alloc_col_forcing(npt)
	  DO kt = 1, npt

          call load_data( kt,GEM_rpn_list )
	      call prepare_cice_sfc(kt)

	      IF (kt .eq. 1) then 
	          call get_col_position(lat_col,lon_col,i_col,j_col)
	          print *, lat_col,lon_col,i_col,j_col
	      ENDIF

	      call prep_col_forcing(i_col,j_col,kt)
	  
	  ENDDO
	
      END SUBROUTINE prepare_gem_forcing	

    END MODULE gemdrv_get_forcing
