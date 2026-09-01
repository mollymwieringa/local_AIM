
import shutil
import os
import sys
import array
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime, date, time, timedelta
from contextlib import closing
import cmocean as cmo
from scipy import stats


class SfcRetrieval:
	
	
    def __init__(self,config = None, Data = None, Buoy = None):
		# Initialising the surface retrieval algorithm
		
		#Reading options from the namelist:
        self.STEP = int(config['STEP'])  #1 if there is a "plateau", 0 otherwise
        self.HEAT = int(config['HEAT']) #1 if there are heating cycles (SIMBAs), 0 otherwise 
        self.sonar = int(config['sonar']) #1 if IMB has sonar data (SIMB3 buoys)
        self.sno_ice= int(config['SnowIce']) #1 if we expect snow-ice, 0 to keep the ice surface immobile.
        
        self.dtemp = float(config['DeltaTmin_snow']) # Temperature threshold to detect snow
        self.dtemp_i = float(config['DeltaTmin_ice']) # Temperature threshold to detect ice
        self.temp = self.dtemp + self.dtemp_i # Total minimum temperature difference in profile
        
        self.crit = float(config['dTdZcrit']) 
        self.dTdzSnow = float(config['dTdZsnow'])
        
        self.Datamin = int(config['DataPositionMin'])
        self.AirLim = int(config['AirPositionLim'])
        self.SnowLim = int(config['SnowPositionLim'])
        self.IceLim = int(config['IcePositionLim'])
        self.OceLim = int(config['OceanPositionLim'])
        self.init_time = int(config['StartLag'])
        
        self.bottom_init = 0.0
        
        self.sechour = 60*60
        
        
        
        a,b = Data.shape
        if Buoy == None:
            self.nsensors = b
            self.tstep = 6
        else:
            self.nsensors = Buoy.nsensors
            self.tstep = Buoy.tstep
            
            
        self.z = np.zeros((self.nsensors+1,1))
        for ks in range(1, self.nsensors+1):
            self.z[ks] = ks*-2.0
        
        
        #----------------------------------------
        # Set analysis vectors and arrays : derivatives
        #----------------------------------------
            
        self.first_der = self.nanfill((a,self.nsensors+1))
        self.second_der = self.nanfill((a,self.nsensors+1))
        self.topsnow = self.nanfill((a,1))
        self.bottomsnow = self.nanfill((a,1))
        self.snowthick = self.nanfill((a,1))

        self.icetop = self.nanfill((a,1))
        self.icebottom = self.nanfill((a,1))
        self.icethick = self.nanfill((a,1))
        self.air_temperature = self.nanfill((a,1))
        self.hs_int = self.nanfill((a,1))
        self.bot_snow_int = self.nanfill((a,1))

        self.snowbot_mode = np.nan
        self.icetop_mode = np.nan
        self.LengthPlateau = 0
        

        #----------------------------------------
        # Set analysis vectors and arrays : minimisation
        #----------------------------------------
        
        self.topsnow_minim = self.nanfill((a,1))
        self.bottomsnow_minim = self.nanfill((a,1))
        self.step_minim = self.nanfill((a,1))
        self.snowthick_minim = self.nanfill((a,1))
        self.icetop_minim = self.nanfill((a,1))
        self.icebottom_minim = self.nanfill((a,1))
        self.icethick_minim = self.nanfill((a,1))
        self.snowice_minim = self.nanfill((a,1))
        
        self.hi_int = self.nanfill((a,1))
        self.hs_int = self.nanfill((a,1))
        self.snowice = self.nanfill((a,1))
        self.congelation = self.nanfill((a,1))
    
        self.Delt = int(config['Delt'])

        self.t = self.nanfill((int(a*self.tstep/24),1))
        self.indice = self.nanfill((int(a*self.tstep/24),1))
        
        # Temperature time series of snow layers, ice layers, interfaces
        self.Ttop = self.nanfill((a,1))
        self.Tsnow_layer = self.nanfill((a,1))
        self.Tinterface = self.nanfill((a,1))
        self.Tbottom = self.nanfill((a,1)) 
        self.Ndata = a






    def compute_inferfaces_minimisation(self,Buoy = None,Data=None):
    #-------------------------------------------------------------------
    #-------------------------------------------------------------------
    #        GET THE SNOW AND ICE INTERFACE: FROM MINIMATION
    #        
    #  This function finds the position of the material interfaces by the
    #  minimation of error functions on a theoretical temperature 
    #  profile.
    #  
    #  The vertical profiles are described by the following points to 
    #  be retrieved:
    #  
    #  Z1, T1 : Position and temperature of the air-snow interface
    #  Z2, T2 : Position and temperature of the snow-ice surface
    #  Zice0, Tice0 : Original ice surface (at deployment) 
    #  Zp, Tp : Other end of surface plateau (SAMS buoys)
    #  Zc, Tc : Sensor located inside the ice
    #  Zio, Tio : Position and temperature of the ice-ocean interface
    #-------------------------------------------------------------------

        #   Start by computing derivatives, and find initial guess on the 
        #   interface positions
        self.compute_derivatives(Buoy = Buoy, Data=Data)
        self.Initial_estimates_from_derivatives(Buoy = Buoy, Data=Data)

        a,b = Data.shape
        print('finding interfaces from minimisation')

        #-------------------------------------
        #  Define positions of the sensors that are laying flat on the ice
        #     based on the invariant kink in Z-derivatives
        #     They are defined as positions Za and Zb.
        #-------------------------------------
        
        # Take initial guess values from the derivatives
        Zice0 = int(np.squeeze(self.snowbot_mode))
        Zp = int(np.squeeze(self.icetop_mode))
        
        #Populate the time series of upper sensor of the flat section,
        #which remains constant.
        self.step_minim[:] = Zice0   

        # Get the ocean surface temperature from deeper sensors
        Tio, dummy = stats.mode(Data[0,self.OceLim-10:self.OceLim])
        
        #initialized previous itterate memory terms used later in the algo.
        Tice0_prev = -3.0
        Zice0_prev = int(Zice0)
        Zsi_init = int(Zice0)
        Tminim0 = 1
        Tmin_prev = -3.0
        kyesterday = 0
        
        
        # Starting loop over profiles 
        for k1 in range(0,int(a)):
			
			#initialise the interface location as nans
            self.topsnow_minim[k1] = np.nan
            self.bottomsnow_minim[k1] = np.nan
            self.icebottom_minim[k1] = np.nan  
            
            
            
            # k is an indice for daily data points,
            kday = int(np.min([int(k1*self.tstep/24),int(a*self.tstep/24)-1]))
            #set the data interval corresponding to the current day
            k_4days =np.max([0,(kday-4)*int(24/self.tstep)])
            # 
            #Find the time in the day with coldest temperature. Use as condition for ice to be cold in last 24h  
            #First guess is data close to 3-6am
            k_night = np.nanmin([int( kday*(24/self.tstep)+self.init_time),int(a)])             
            Tmin_init = -3.0
            for deltak in range (int(-24/(self.tstep*4)), int(24/(self.tstep*4)+1)):                    
                Tmin = np.nanmin(Data[int(k_night + deltak),
                                 self.Datamin:int(self.nsensors-10)])
                if (Tmin < Tmin_init):
                    k_minTemp = k_night  + deltak
                    Tmin_init = Tmin                        
            k_night = int(k_minTemp)           
            
 
            #from now on, k is the "day number" and k1 is the actual data row indice

            if (np.isnan(Data[k1,2])):
                #Skip rows with invalid data
                print('NaN at k1 = ', k1   )
            else:
                
                
                #-------------------------------------
                # INITIALIZING:
                # Defining initial reference points in the profiles
                # 
                # Get the temperature at the reference pts (Tmin,Tice0,Tp,Tocean)
                #-------------------------------------

                # Position and temperature of lowest temperature in profile. Used as initial guess                
                Tmin= np.nanmin(Data[k1,self.Datamin:self.nsensors-10])
                Zmin = np.argmin((Data[k1,self.Datamin:self.nsensors-10]-Tmin)**2.0)+self.Datamin
                
                # Air temperature
                Tair = np.mean(Data[k1,self.Datamin:self.Datamin+5])
            
                # Ocean temperature.
                Tocean_mode, dummy = stats.mode(Data[k1,self.nsensors-20:self.nsensors-5])
      
                # Wire plateau
                Tice0 = Data[k1,Zice0]
                Tp = Data[k1,Zp]
                
                #Tc: A random point in the ice: with a temperature between Tb and Tocean. 
                #Set to 1/3rd of the temperature difference so that we are closer to 
                #the ice bottom.
                Tc = np.max([(Tp-Tocean_mode)/3.0 + Tocean_mode, -5.0])
                Zc = int(np.argmin((Data[k1,Zp:self.nsensors-10]-Tc)**2.0))+Zp
                Tc = Data[k1,Zc]
    
                # Re-initialise profile ice temperature gradient 
                beta_ice = 0.0
                Tanalytic_rec = self.nanfill((self.nsensors+1))
    
    
                # Condition to continue: negative gradients in snow and ice
                if ((Zmin <= Zice0) & (Zc>Zp) ):
					      
					      					      
                    #-----------------------------------------------------------------
                    # 1. ICE BOTTOM: Found by minimizing about the bottom kink
                    #
                    #  Here, we make 20 test theoretical profiles, each with ice-ocean 
                    #  interface set at a given sensor position. The real 
                    #  interface is then set based on the profile that
                    #  minimes the error on a piece-wise linear profile.
                    #
                    #  We only scan from the sensor where T=-2.2C. 
                    #  This is to reduce computation cost while making sure 
                    #  that we start above the interface, but close to it. 
                    #-----------------------------------------------------------------
                    
                    
                    # Condition to continue: minimum deltaT in ice
                    if ((Tp < Tocean_mode - self.dtemp_i)):
						
                        #Initializng the theoretical profile:
                        Zio_init = int(np.argmin((Data[k1,Zp:self.nsensors-10]+2.2)**2.0))+Zp
                        Tio_init = Data[k1,Zio_init]
                        Tanalytic_rec = self.nanfill((self.nsensors+1))
                        Err_min = 9999.0
                        
                        #Scan the sensors to find theoretical profile that best match the obs. 
                        for dz in range(-10,10):
                            
                            #re-initialise the theoretical profiles
                            Tanalytic = self.nanfill((self.nsensors+1))   #reset the piece-wise linear profile
                            Zio_Test = np.nanmax([Zc,Zio_init + dz])     #Test sensor position, limited to below Zc
                            Tio_Test = Tocean_mode                       
                            beta_ice_test = (Tio_Test-Tc)/(Zio_Test-Zc)  # temperature slope based on Tc and the test position

                            # Compute the error between the test profile and the real obs:
                            Err = 0.0
                            for kz in range(int(Zio_Test-10),int(np.nanmin([Zio_Test+10, self.nsensors-10]))):
                                if (kz <= Zio_Test):
                                    Tanalytic[kz] = beta_ice_test*(kz-Zc) + Tc
                                    Err = Err + (Tanalytic[kz] - Data[k1,kz])**2.0
                                else :
                                    Tanalytic[kz] = Tio_Test
                                    Err = Err + (Tanalytic[kz] - Data[k1,kz])**2.0
                                    
                            #Record the position and temperature if the error is smaller than previous iterates        
                            if (Err < Err_min):
                                Err_min = Err
                                Zio = Zio_Test
                                Tio = Tio_Test
                                beta_ice = beta_ice_test
                                Tanalytic_rec = Tanalytic
    
             
                        self.icebottom_minim[k1] = Zio
                        self.Tbottom[k1] = Data[k1,Zio]
                        
                        #Record the first trustworthy bottom ice position (after 4 days, so that the hole is refrozen)
                        if self.bottom_init == 0 and k1>15:
                            self.bottom_init = Zio
                            
                
                
                #------------------------------------------------------------
                # 2. Air-snow interface.
                #      as the highest point where grad-T exceeds that in the ice
                #
                # We move the air-snow interface to test location if:
                #   - The inflection is largest
                #   - The slope in snow below is larger than a set threshold
                #   - The slope in air above is smaller then threshold*2
                # 
                #------------------------------------------------------------

                # Conditions to continue:  
                #  - minimum deltaT in the snow, for at least 6h,
			    #  - no temperature gradient reversal in snow (piecewise linear assumption)	 
			    # This is added to filter rare cases under rapid cooling weather, 
			    # when the air temperature conditions is met
				# but the snow layer remains warmer than the ice below
				
                    Tmax = np.max(Data[k1,Zmin:Zice0])
                    if (((Tmin < Tice0 - self.dtemp) & 
                         (Tmin_prev < Tice0_prev - self.dtemp)) & 
                         ((Tmax < Tice0+0.25) & (Tair <= Tice0))):

						#re-initialise parameters for new profile
                        Curv_init = 0.0
                        Zas = Zmin     #Initial guess: snow top is where Tmin is.
            
                        #Span from Tmin to the original ice surface to test if it is the interface
                        for z_test in range(Zmin,int(Zice0)):
							
							#Temperature at test location
                            Tas_test = Data[k1,z_test]
                            
                            #slope over 8cm
                            T_above = Data[k1,z_test-2]
                            T_below = Data[k1,z_test+2]
                            beta_snow = (T_below-T_above)/2
                            #slope over 4cm above and below the test location
                            beta_snow_above = (Tas_test-T_above)/(-4.0)
                            beta_snow_below = (T_below-Tas_test)/(-4.0)
                        
                            #Curvature over 8cm
                            Curv =(T_below-2.0*Tas_test+T_above)/(2.0)
                        
                            # Determining if the location is the interface
                            if (Curv > Curv_init) & (-beta_snow_below > self.dTdzSnow) & (-beta_snow_above < self.dTdzSnow*2.0):
                                Zas = int(np.nanmax([z_test,Zas]))
                                Tas = Data[k1,Zas]
                                beta_snow = (Tas-Tice0)/(Zas-Zice0)
                                Curv_init = Curv
                                for kz in range(Zas,Zice0):
                                    Tanalytic_rec[kz] =  Tas + beta_snow*(kz-Zas)
                           

                #------------------------------------------------------------
                # 3. Snow-ice (flooding)
                #       - based on where the slope exceeds that in the ice
                #------------------------------------------------------------
        
                        #Initializing
                        
                        beta_snow = (Tas-Tice0)/(Zas-Zice0)
                        beta_diff = beta_snow  - beta_ice
                        Zsi = Zas
                        Zsi_high = Zas
                        Zsi_low = Zas
                        
                        # Find interval with slopes between snow and ice
                        for Zsi_test in range(int(Zas+1),Zice0+1):
							
                            #Get the temperatures and gradients about the test position
                            Tsi_test = Data[k1,Zsi_test]
                            T_above = Data[k1,Zsi_test-1]
                            T_below = Data[k1,Zsi_test+1]
                            beta_above = (Tsi_test-T_above)/(2.0)
                            beta_below = (T_below-Tsi_test)/(2.0)
                            beta_test = (T_below-T_above)/(2.0)
                        
                            #Find the lowest position with a gradient close to that in snow
                            if ((beta_test > beta_snow - beta_diff*0.2)):
                                Zsi_high = np.nanmax([Zsi_high, Zsi_test])
                            #Find the lowest position with slope larger than in the ice
                            if ((beta_test > beta_snow - beta_diff*0.8)):
                                Zsi_low = np.nanmax([Zsi_low, Zsi_test])

                        
                        
                        #Define the interface by the point with largest curvature
                        Curv_si_init = 0.0
                        Zsi = Zsi_high
                        for z_test in range(Zsi_high,Zsi_low+2):
                            #Compute curvature at point
                            T_above = Data[k1,z_test-2]
                            T_test = Data[k1,z_test]
                            T_below = Data[k1,z_test+2]
                            Curv_si =(T_below-2.0*T_test+T_above)/(2.0)
          
                            if (Curv_si < Curv_si_init):
                                Zsi = np.nanmax([Zsi,z_test])
                                Curv_si_init = Curv_si 
                                
                        Zsi = np.min([Zsi,Zice0])
                        Zsi_prev = np.min([Zsi,Zice0])           
                        Tsi = Data[k1,Zsi]
                   
                        #We only accept snow-ice if it does not get back down in the following 4 days.
                        #I.e., reset the interface if it is not permanent
                        kyesterday = kday*int(24/self.tstep)
                        if kyesterday == k1:
                            kyesterday = int(np.max([0,(kday-1)*int(24/self.tstep)]))
                        if ((np.nanmin(self.bottomsnow_minim[kyesterday:k1])< Zsi_prev) & (Zsi <= Zice0)):
                            interval_4day = self.bottomsnow_minim[k_4days:k1] 
                            interval_4day[self.bottomsnow_minim[k_4days:k1]<Zsi_prev] = Zsi_prev
                            self.bottomsnow_minim[k_4days:k1] = interval_4day
                            self.icetop_minim[k_4days:k1] = interval_4day


                #-----------------------------------------------------------------
                # End of retrival.
                # Writing results in diagnostic time-series
                #-----------------------------------------------------------------
   
                        self.step_minim[k1] = Zice0
                        self.topsnow_minim[k1] = Zas
                        self.bottomsnow_minim[k1] = Zsi
                        self.snowthick_minim[k1] = (Zsi*self.sno_ice + Zice0*(1-self.sno_ice) - Zas)*2.0
                            
                        self.Ttop[k1] = Data[k1,Zas]
                        self.Tsnow_layer[k1] = Data[k1,int((Zsi+Zas)/2)]
                        self.Tinterface[k1] = Data[k1,Zsi]
                        self.Tbottom[k1] = Tocean_mode
                        
                        #Ice thickness, with the snow-ice and removing the plateau section.
                        if (self.STEP == 1):
                            self.icetop_minim[k1] = Zsi
                            if (Tp < Tocean_mode - self.dtemp_i):
                                self.icethick_minim[k1] = (Zio - Zp)*2.0 + (Zice0 - Zsi)*2.0*self.sno_ice
					
                        else:
                            self.icetop_minim[k1] = Zsi
                            if (Tp < Tocean_mode - self.dtemp_i):
                                self.icethick_minim[k1] = (Zio - Zice0)*2.0 + (Zice0 - Zsi)*2.0*self.sno_ice

                # Update the history parameters
                if k1 == k_night:
                    Tmin_prev = Tmin
                    Tice0_prev = Tice0 
                    kyesterday = k_night        


        #---------------------------------------------------------------
        #  Mask data in the arrays and perform daily running means
        #---------------------------------------------------------------

        #set snowice layer to nan if there is nan in the snow layer            
        self.bottomsnow_minim[np.isnan(self.topsnow_minim[:])] = np.nan
        
        #Prepare mask to hide points in the running means if no data in last or next 18 hours.
        topsnow_tmp =self.topsnow_minim.copy()
        icebottom_tmp=self.icebottom_minim.copy()
        bottomsnow_tmp=self.bottomsnow_minim.copy()
        
        topsnow_tmp2 =self.topsnow_minim.copy()
        icebottom_tmp2=self.icebottom_minim.copy()
        bottomsnow_tmp2=self.bottomsnow_minim.copy()
        
        topsnow_tmp3 =self.topsnow_minim.copy()
        icebottom_tmp3=self.icebottom_minim.copy()
        bottomsnow_tmp3=self.bottomsnow_minim.copy()
        
        topsnow_tmp[1:-1]= ~np.isnan(self.topsnow_minim[0:-2])+~np.isnan(self.topsnow_minim[1:-1])+~np.isnan(self.topsnow_minim[2:])
        icebottom_tmp[1:-1] = ~np.isnan(self.icebottom_minim[0:-2])+~np.isnan(self.icebottom_minim[1:-1])+~np.isnan(self.icebottom_minim[2:])
        bottomsnow_tmp[1:-1] = ~np.isnan(self.bottomsnow_minim[0:-2])+~np.isnan(self.bottomsnow_minim[1:-1])+~np.isnan(self.bottomsnow_minim[2:])

        topsnow_tmp2[0:-2]= ~np.isnan(self.topsnow_minim[0:-2])+~np.isnan(self.topsnow_minim[1:-1])+~np.isnan(self.topsnow_minim[2:])
        icebottom_tmp2[0:-2] = ~np.isnan(self.icebottom_minim[0:-2])+~np.isnan(self.icebottom_minim[1:-1])+~np.isnan(self.icebottom_minim[2:])
        bottomsnow_tmp2[0:-2] = ~np.isnan(self.bottomsnow_minim[0:-2])+~np.isnan(self.bottomsnow_minim[1:-1])+~np.isnan(self.bottomsnow_minim[2:])

        topsnow_tmp3[2:]= ~np.isnan(self.topsnow_minim[0:-2])+~np.isnan(self.topsnow_minim[1:-1])+~np.isnan(self.topsnow_minim[2:])
        icebottom_tmp3[2:] = ~np.isnan(self.icebottom_minim[0:-2])+~np.isnan(self.icebottom_minim[1:-1])+~np.isnan(self.icebottom_minim[2:])
        bottomsnow_tmp3[2:] = ~np.isnan(self.bottomsnow_minim[0:-2])+~np.isnan(self.bottomsnow_minim[1:-1])+~np.isnan(self.bottomsnow_minim[2:])


        #Transforming the data into daily means
        self.topsnow_minim[2:-1] = np.nanmean([self.topsnow_minim[0:-3],self.topsnow_minim[1:-2],self.topsnow_minim[2:-1],self.topsnow_minim[3:]],0)
        self.icebottom_minim[2:-1] = np.nanmean([self.icebottom_minim[0:-3],self.icebottom_minim[1:-2],self.icebottom_minim[2:-1],self.icebottom_minim[3:]],0)
        self.bottomsnow_minim[2:-1] = np.nanmean([self.bottomsnow_minim[0:-3],self.bottomsnow_minim[1:-2],self.bottomsnow_minim[2:-1],self.bottomsnow_minim[3:]],0)

        #Apply the mask
        self.topsnow_minim[topsnow_tmp[:]==False] = np.nan
        self.icebottom_minim[icebottom_tmp[:]==False] = np.nan     
        self.bottomsnow_minim[bottomsnow_tmp[:]==False] = np.nan   
        self.topsnow_minim[topsnow_tmp2[:]==False] = np.nan
        self.icebottom_minim[icebottom_tmp2[:]==False] = np.nan     
        self.bottomsnow_minim[bottomsnow_tmp2[:]==False] = np.nan   
        self.topsnow_minim[topsnow_tmp3[:]==False] = np.nan
        self.icebottom_minim[icebottom_tmp3[:]==False] = np.nan     
        self.bottomsnow_minim[bottomsnow_tmp3[:]==False] = np.nan   
        
        self.step_minim[:] = Zice0
        #---------------------------------------------------------------
        #  Record the results in ice mass balance time series 
        #---------------------------------------------------------------
                           
        if (self.STEP == 1): #if there is a plateau
            self.icetop_minim = self.bottomsnow_minim.copy()
            self.icethick_minim = (self.icebottom_minim - Zp)*2.0 + (Zice0 - self.bottomsnow_minim) *2.0*self.sno_ice
								
        else: # If there is no plateau, e.g. for SIMB3 buoys
            self.icetop_minim = self.bottomsnow_minim.copy()
            self.icethick_minim = (self.icebottom_minim - Zice0)*2.0 + (Zice0 - self.bottomsnow_minim) *2.0*self.sno_ice

        self.snowthick_minim = (self.bottomsnow_minim*self.sno_ice + Zice0*(1-self.sno_ice) - self.topsnow_minim)*2.0

        self.hs_int = self.snowthick_minim.copy()
        self.hi_int = self.icethick_minim.copy()        
        self.snowice = (self.step_minim - self.bottomsnow_minim)*2.0
        self.congelation = (self.icebottom_minim-self.bottom_init)*2.0
        
        # Making 1-day running means of snowice and congelation
        self.snowice[2:-2] = np.nanmean([self.snowice[0:-4],self.snowice[1:-3],self.snowice[2:-2],self.snowice[3:-1],self.snowice[4:]],0)
        self.congelation[2:-2] = np.nanmean([self.congelation[0:-4],self.congelation[1:-3],self.congelation[2:-2],self.congelation[3:-1],self.congelation[4:]],0)
        self.congelation[np.isnan(np.squeeze(self.hi_int))] = np.nan

 




    def compute_derivatives(self,Buoy = None,Data=None):
    #-------------------------------------------------------------------
    #-------------------------------------------------------------------     
    #  This function computes the first and second vertical derivatives
    #  of the temperature data.
    #  
    #  This method to detect the interfaces is not perfoming well, but
    #  it is useful to determine the location of the "thermistor string
    #  plateau" in the SAMS buoy, and the initial snow-ice interface in
    #  any buoy, which is needed in the more robust minimation retrieval
    #  algorithm .      
    #-------------------------------------------------------------------

        a,b = Data.shape
            
        #----------------------------------------
        # Set vertical profile (sensor positions)
        #----------------------------------------
        z = np.zeros((self.nsensors+1,1))
        for k in range(1, self.nsensors+1):
            z[k] = k*-2.0
        zf = len(z)
  
  
        #----------------------------------------
        # Calculate derivatives
        #----------------------------------------
  
        for k1 in range(0, a-1):
  
            # Compute the derivative (center difference)
            self.air_temperature[k1] = np.mean(Data[k1,5:15])
            for kz in range(1, zf-1):
                self.first_der[k1,kz] = (Data[k1,kz+1] - Data[k1,kz-1])/4.0
      
            # Compute the second derivative (center difference)
            for kz in range(1, zf-1):
                self.second_der[k1,kz] = (self.first_der[k1,kz+1] - self.first_der[k1,kz-1])/4.0



    def Initial_estimates_from_derivatives(self,Buoy = None,Data=None): 
	#-------------------------------------------------------------------
    #-------------------------------------------------------------------
    #        
    #  This function approximates the location of the 
    #  material interfaces from the derivative.
    #  
    #  This method to detect the interfaces is not perfoming well (it is noisy), but
    #  it is useful initial guess in the more precised minimation retrieval
    #  algorithm.    
    #
    #  Estimates for the snow and ice layers are are processed independently  
    #-------------------------------------------------------------------

        #-----------------------------------------------
        # Snow layer
        #-----------------------------------------------
		
        a,b = Data.shape
        for k1 in range(0, a-1):
                
            #------------------------------------------------------
            # Find the minimum temperature at the top of the chain
            #------------------------------------------------------
            # We find the indice of the minimum temperature, 
            # assuming that the snow layer is below this point
                
            # Note that Here, we skip the 10 top sensors (20 cm) 
            # and assume that the snow syrface is <140cm away... 
            # Needs improvement
                
            data_time_k = Data[k1,self.Datamin:self.AirLim]
            airTmin = np.amin(data_time_k)
                
            start_ind = np.argmin(data_time_k)+self.Datamin

            # If the air temperature is too high, there is no gradient
            # and we skip the computations
            if ((airTmin < -self.temp) and (np.isnan(airTmin)==0)):
                    
                #Not sure why we need this...
                if self.STEP == 1 and (start_ind > self.AirLim-self.Datamin):
                    print('Start is below 60, does not work for k1 = ', k1 )
                    nsnowmax = np.nan
                else:
                        
                    #-----------------------------------------------
                    # Find the position of air-snow interface
                    #-----------------------------------------------
                        
                    #The maximum derivative ~ middle of the snow layer
                    TopChainDeriv = self.first_der[k1,start_ind:self.SnowLim]
                    midsnow_loc = np.argmax(TopChainDeriv)+start_ind
          
                    #find the top of the snow as the nearest point above
                    # the largest slope with 
                    # a gradient <0.1 (self.dTdzSnow)
                    TopChainDeriv = self.first_der[k1,0:midsnow_loc]
                    N = len(TopChainDeriv)
                    indices = np.arange(N)
                    condition = indices[TopChainDeriv<self.dTdzSnow]
                    N = len(condition)   
                    if N == 0:
                        print("No snow detected!!!!")
                    else:    
                        highest_snow_sensor = np.nanmax(condition)
                        self.topsnow[k1] = highest_snow_sensor
        
                        #-----------------------------------------------
                        # Find the position of snow - ice interface
                        #-----------------------------------------------
         
                        #Method: assuming that the interface is at the largest curvature point.
                        if self.STEP == 1:
                            MaxLim = self.SnowLim
                        else: 
                            MaxLim = self.IceLim #we are spanning wider
                        SnowLayerInflection = self.second_der[k1,int(self.topsnow[k1]+1):MaxLim]
                        lowest_snow_sensor = np.argmin(SnowLayerInflection)+int(self.topsnow[k1]+1)
                            
                        self.bottomsnow[k1] = lowest_snow_sensor
                        self.snowthick[k1] = (self.bottomsnow[k1]-self.topsnow[k1])*2.0
                        
                        if self.STEP==0:
                            self.icetop[k1] = lowest_snow_sensor2
            
        # Interpolating to fill the gaps in data...
        self.bot_snow_int  = self.interp_data(Data = np.squeeze(self.bottomsnow))
        self.hs_int  = self.interp_data(Data=np.squeeze(self.snowthick))
    
    
        #-----------------------------------------------
        # Ice layer
        #-----------------------------------------------
        # We do an other loop so that if the snow 
        # retrieval fails, we can nonetheless find the ice based on 
        # interpolated values.
            
        for k1 in range(0, a):
            
            data_time_k = Data[k1,self.Datamin:self.AirLim]
            airTmin = np.amin(data_time_k)
            
            if ((airTmin < -self.temp) and (np.isnan(airTmin)==0)):
                
                if ((np.isnan(self.bot_snow_int[k1])==0)):
                 
                    #-----------------------------------------------
                    # Find the position of top ice interface
                    #-----------------------------------------------
                    #    This is done assuming that position snow bottom 
                    #    is below the top snow. As the bottom snow is
                    #    previously guessed by the inflection point, we 
                    #    assume that the maximum gradient now will be in
                    #    the ice, thus below the plateau.
                        
                    # finding the maximum gradient below the plateau
                        
                    
                    if self.STEP == 1:
                        start_ind = int(self.bot_snow_int[k1])+3  
                        if start_ind > self.SnowLim:
                            print('does not work for k1 = ', k1 )
                        else:
                            BelowSnowDeriv = self.first_der[k1,start_ind:self.IceLim]
                            topice_guess = np.argmax(BelowSnowDeriv)+start_ind
                            #print(start_ind,topice_guess)
                            Inflection = self.second_der[k1,start_ind:topice_guess]
                            N = len(Inflection)
                            if (N == 0):
                                print('Top ice does not work for k1 = ', k1 )
                            else:
                                self.icetop[k1] = np.argmax(Inflection)+start_ind
                    else:    
                        
                        start_ind = np.argmin(data_time_k)+self.Datamin  
                        BelowSnowDeriv = self.first_der[k1,start_ind:self.SnowLim]
                        topice_guess = np.argmax(BelowSnowDeriv)+start_ind
                        Inflection = self.second_der[k1,topice_guess:self.IceLim]
                        N = len(Inflection)
                        if (N == 0):
                            print('Top ice does not work for k1 = ', k1 )
                        else:
                            self.icetop[k1] = np.argmin(Inflection)+topice_guess
                            
                    #-----------------------------------------------
                    # Find the position of ice-ocean interface
                    #-----------------------------------------------
                
                if ((np.isnan(self.icetop[k1])==0)):    
                    IceDeriv = self.first_der[k1,int(self.icetop[k1]+1):self.OceLim]
                            
                    N = len(IceDeriv)
                    indices = np.arange(N)
                    condition = indices[IceDeriv<self.crit]
                    TopOcean_loc = np.nanmin(condition)
                            
                    #Ice bottom is directly above the upper ocean 
                    #sensor
                    
                    self.icebottom[k1] = TopOcean_loc+self.icetop[k1]
                    self.icethick[k1] = (self.icebottom[k1]-self.icetop[k1])*2.0
    
        self.snowbot_mode, dummy= stats.mode(self.bot_snow_int[0:int((24.0/self.tstep)*7)])
        self.icetop_mode, dummy = stats.mode(self.icetop[0:int((24.0/self.tstep)*7)])
        self.LengthPlateau = self.icetop_mode - self.snowbot_mode
        



    def nanfill(self,siz):
    #This creates an array of nans:
        A = np.zeros(siz)
        A.fill(np.nan)
        return A


    #--------------------------------------------------
    def interp_data(self,Data = None):
        a = len(Data)
        val_interp = Data.copy()
        val_dum = np.arange(a)
        val_nan = val_dum[np.isnan(Data)==1]
        val_1 = val_dum[np.isnan(Data)==0]
        for k_nan in range(0, len(val_nan)):
            for k_1 in range(0, len(val_1)-1):
                if ( val_1[k_1]<val_nan[k_nan]) and ( val_1[k_1+1]>val_nan[k_nan]):
                    indL = val_1[k_1]
                    indR = val_1[k_1+1]
                    val_interp[int(val_nan[k_nan])] = Data[indL] + (val_nan[k_nan]-indL)*(Data[indR]- Data[indL])/(indR-indL)
        return val_interp


   
    def prep_icepack_initfiles(self, Data = None, OutputFolder = None, Datelist = None):

        a = len(self.hi_int)
        for k1 in range(0, a):
            if (np.isnan(self.hi_int[k1])==0):

                if (np.isnan(self.hs_int[k1])==0) and (np.isnan(self.hi_int[k1])==0):
                    tops = int(self.topsnow_minim[k1])
                    bots = int(self.bottomsnow_minim[k1])
                    boti = int(self.icebottom_minim[k1])

                    Tprofile_snow_data = Data.data[k1,tops:bots]
                    ns = len(Tprofile_snow_data)
                    Tprofile_ice_data1 = Data.data[k1,bots:int(self.snowbot_mode)]
                    Tprofile_ice_data2 = Data.data[k1,int(self.icetop_mode):boti]
                    ni1 = len(Tprofile_ice_data1)
                    ni2 = len(Tprofile_ice_data2)
                    ni = len(Tprofile_ice_data1)  + len(Tprofile_ice_data2)
                    Tsfc = Data.data[k1,tops]
                    Data_out = np.squeeze(np.zeros((ni+ns+8,1)) )
                    Data_out[0] = Data.Lat[k1]
                    Data_out[1] = Data.Lon[k1]
                    Data_out[2] = Tsfc
                    Data_out[3] = self.hs_int[k1]
                    Data_out[4] = self.hi_int[k1]
                    Data_out[5] = ns
                    Data_out[6] = ni
                    Data_out[7:ns+7] = Tprofile_snow_data
                    Data_out[ns+7:ni1+ns+7] = Tprofile_ice_data1
                    Data_out[ni1+ns+7:ni2+ni1+ns+7] = Tprofile_ice_data2
                    
                    filename_out = "%s%s" % (OutputFolder, Datelist[k1])

                    np.savetxt(filename_out, Data_out, delimiter=",")
                    
                    


    def PrintMassBalance(self,Data = None, OutputFolder = None):
		
        DateMax = Data.datelist[np.squeeze(self.hi_int[:]==np.nanmax(self.hi_int[:]))]
        DateMax = DateMax[0] 
        Congel = self.congelation[np.squeeze(self.hi_int[:]==np.nanmax(self.hi_int))]
        Congel = Congel[0] 
        Snice = self.snowice[np.squeeze(self.hi_int[:]==np.nanmax(self.hi_int))]
        Snice = Snice[0]              
        with open('%sPrinted_MassBalance.txt' % OutputFolder, 'w') as f:			
            f.write('Initial ice thickness and snow depth with date: %s, %s, %s \n' % (self.hi_int[22],self.hs_int[22],Data.datelist[22])) 
            hfb = self.hi_int - (self.hs_int*330 + self.hi_int*917)/1026
            f.write('Initial freeboard is : %s \n' % hfb[22])
            f.write('Max ice thickness and date: %s, %s   \n' % (np.nanmax(self.hi_int),DateMax)) 
            f.write('Total ice congelation and snow-ice: %s, %s   \n' % (Congel, Snice)) 
            DateMax = Data.datelist[np.squeeze(self.hs_int[:]==np.nanmax(self.hs_int))]
            DateMax = DateMax[0] 
            f.write('Max snow thickness and date: %s, %s\n' % (np.nanmax(self.hs_int),DateMax)) 
            for k in range(0,len(hfb)):
                f.write('hi, hs, freeboard is : %s, %s, %s,,,,,%s, %s,              on %s \n' % (self.hi_int[k],self.hs_int[k],hfb[k], self.congelation[k],self.snowice[k],Data.datelist[k]))

        return
        
        
#End of programm        
        
