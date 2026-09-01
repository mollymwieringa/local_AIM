#!/usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from datetime import datetime, date, time, timedelta
import os
import sys
import cmocean as cmo

from SfcRetrieval import SfcRetrieval
import configparser


from TimeUtil import TimeUtil

from SAMSIMBdata import SAMSIMB_DAsetup
from SAMSIMBdata import SAMSIMBdata
from SAMSIMBdata import SAMSIMBheatdata

#------------------------------------------------------------
# Loading the IMB1 data
#------------------------------------------------------------

OutputFolder = 'Outputs/'
visuals = vis()

config_IMB1 = configparser.ConfigParser()
config_IMB1.read('Buoy_data/namelist_IMB1.ini')
BuoySetup = SAMSIMB_DAsetup(config=config_IMB1['IMB_setup'])
timeIMB1 = TimeUtil(config = config_IMB1['Time'])
IMB1Data = SAMSIMBdata(ExpSetup=BuoySetup, time = timeIMB1)
IMB1_heat = SAMSIMBheatdata(ExpSetup=BuoySetup, time = timeIMB1)

Rtrvl_IMB1 = SfcRetrieval(config=config_IMB1['Retrieval'],Data = IMB1Data.data, Buoy = BuoySetup)
Rtrvl_IMB1.compute_inferfaces_minimisation(Data= IMB1Data.data,Buoy= BuoySetup)
Rtrvl_IMB1.PrintMassBalance(Data = IMB1Data, OutputFolder = OutputFolder +'IMB1_MassBalance')

#------------------------------------------------------------
# Loading the IMB2 data
#------------------------------------------------------------

config_IMB2 = configparser.ConfigParser()
config_IMB2.read('Buoy_data/namelist_IMB2.ini')
BuoySetup = SAMSIMB_DAsetup(config=config_IMB2['IMB_setup'])
timeIMB2 = TimeUtil(config = config_IMB2['Time'])
IMB2Data = SAMSIMBdata(ExpSetup=BuoySetup, time = timeIMB2)
IMB2_heat = SAMSIMBheatdata(ExpSetup=BuoySetup, time = timeIMB2)

Rtrvl_IMB2 = SfcRetrieval(config=config_IMB2['Retrieval'],Data = IMB2Data.data, Buoy = BuoySetup)
Rtrvl_IMB2.compute_inferfaces_minimisation(Data= IMB2Data.data,Buoy= BuoySetup)
Rtrvl_IMB2.PrintMassBalance(Data = IMB2Data, OutputFolder = OutputFolder +'IMB2_MassBalance')


#------------------------------------------------------------
# Icepack analysis
#------------------------------------------------------------



#Standard comparison Bl99 vs. mushy simulations
#Figure 3, 4, 5,6
OutputFolder = 'Outputs/'
visuals = vis()
for krun in range(0, 2):

    if krun == 0: #For the IMB1 site

        BuoyData = IMB1Data
        Rtrvl = Rtrvl_IMB1
        time = timeIMB1
        Data_heat = IMB1_heat

        config = configparser.ConfigParser()
        config.read('Icepack_simulations/IMB1_Bl99_ctrl/namelist_Icepack.ini')

        config_mushy = configparser.ConfigParser()
        config_mushy.read('Icepack_simulations/IMB1_mushy_ctrl/namelist_Icepack.ini')


    if krun == 1:  #For the IMB2 site

        BuoyData = IMB2Data
        Rtrvl = Rtrvl_IMB2
        time = timeIMB2
        Data_heat = IMB2_heat

        config = configparser.ConfigParser()
        config.read('Icepack_simulations/IMB2_Bl99_ctrl/namelist_Icepack.ini')

        config_mushy = configparser.ConfigParser()
        config_mushy.read('Icepack_simulations/IMB2_mushy_ctrl/namelist_Icepack.ini')

    #--------------------------------------------------
    # Icepack BL99 data:
    Bl99_MetaData = IcepackDatasetup(config['Icepack_setup'],config_Labels = config['Diag_infos'])
    Data_Bl99 = IcepackData(meta = Bl99_MetaData)
    Data_Bl99.Interpolate_into_IMB_profile(nsensors=Rtrvl.nsensors-1,meta = Bl99_MetaData )

    #--------------------------------------------------
    # Icepack mushy data:
    mushy_MetaData = IcepackDatasetup(config_mushy['Icepack_setup'],config_Labels = config_mushy['Diag_infos'])
    Data_mushy = IcepackData(meta = mushy_MetaData)
    Data_mushy.Interpolate_into_IMB_profile(nsensors=Rtrvl.nsensors-1,meta = mushy_MetaData )


    visuals.Compute_Errors(Data_buoy = BuoyData,
                           Data_model = Data_Bl99,
                           variable = Data_Bl99.Hiavg*100.0,
                           reference = Rtrvl.hi_int,
                           krun = krun,
                           OutputFolder = Data_Bl99.OutputFolder + 'IMB%s_Bl99_Hi' % (krun+1))

    visuals.Compute_Errors(Data_buoy = BuoyData,
                           Data_model = Data_Bl99,
                           variable = Data_Bl99.Hsavg*100.0,
                           reference = Rtrvl.hs_int,
                           krun = krun,
                           OutputFolder = Data_Bl99.OutputFolder + 'IMB%s_Bl99_Hs' % (krun+1))


    visuals.Compute_Errors(Data_buoy = BuoyData,
                           Data_model = Data_Bl99,
                           variable = Data_Bl99.CongTotal*100.0,
                           reference = Rtrvl.congelation,
                           krun = krun,
                           OutputFolder = Data_Bl99.OutputFolder + 'IMB%s_Bl99_cong' % (krun+1))

    visuals.Compute_Errors(Data_buoy = BuoyData,
                           Data_model = Data_Bl99,
                           variable = Data_Bl99.CumulSnowice*100.0,
                           reference = Rtrvl.snowice,
                           krun = krun,
                           OutputFolder = Data_Bl99.OutputFolder + 'IMB%s_Bl99_snowice' % (krun+1))

    visuals.Compute_Errors(Data_buoy = BuoyData,
                           Data_model = Data_mushy,
                           variable = Data_mushy.Hiavg*100.0,
                           reference = Rtrvl.hi_int,
                           krun = krun,
                           OutputFolder = Data_mushy.OutputFolder + 'IMB%s_mushy_Hi' % (krun+1))

    visuals.Compute_Errors(Data_buoy = BuoyData,
                           Data_model = Data_mushy,
                           variable = Data_mushy.Hsavg*100.0,
                           reference = Rtrvl.hs_int,
                           krun = krun,
                           OutputFolder = Data_mushy.OutputFolder + 'IMB%s_mushy_Hs' % (krun+1))

    visuals.Compute_Errors(Data_buoy = BuoyData,
                           Data_model = Data_mushy,
                           variable = Data_mushy.CongTotal*100.0,
                           reference = Rtrvl.congelation,
                           krun = krun,
                           OutputFolder = Data_mushy.OutputFolder + 'IMB%s_mushy_cong' % (krun+1))

    visuals.Compute_Errors(Data_buoy = BuoyData,
                           Data_model = Data_mushy,
                           variable = Data_mushy.CumulSnowice*100.0,
                           reference = Rtrvl.snowice,
                           krun = krun,
                           OutputFolder = Data_mushy.OutputFolder + 'IMB%s_mushy_snowice' % (krun+1))

    # Figure 3
    visuals.figure_temperatures_weather(Data = BuoyData,
                                Data_model = Data_mushy,
                                Data_station = data_eccc,
                                Rtrvl = Rtrvl,
                                krun = krun,
                                OutputFolder = OutputFolder)

    # Figure 4
    visuals.figure_RtrvlValidation(Data=BuoyData,
                                time=time,
                                Rtrvl = Rtrvl,
                                Data_heat = Data_heat,
                                krun = krun,
                                OutputFolder =  OutputFolder)
    # Figure 5
    visuals.figure_thickness_obs(data_buoy = BuoyData,
                                       rtrvl = Rtrvl,
                                       krun = krun,
                                       OutputFolder = OutputFolder)

    # Figure 6
    visuals.figure_temperature_intercomp(Data_buoy = BuoyData,
                                       Data_model = Data_Bl99,
                                       Data_mushy = Data_mushy,
                                       Rtrvl_buoy = Rtrvl,
                                       krun = krun,
                                       OutputFolder = '%s%s' % (OutputFolder, 'Figure6'))


    visuals.Compute_Errors(Data_buoy = BuoyData,
                           Data_model = Data_Bl99,
                           variable = Data_Bl99.Tair,
                           reference = BuoyData.airT,
                           krun = krun,
                           OutputFolder = OutputFolder + 'IMB%s_GDPSbias' % (krun+1))

del visuals


#Figure 7: BL99 with and without snow ice
visuals = vis()
for krun in range(0, 2):

    if krun == 0:

        BuoyData = IMB1Data
        Rtrvl = Rtrvl_IMB1

        config = configparser.ConfigParser()
        config.read('Icepack_simulations/IMB1_Bl99_ctrl/namelist_Icepack.ini')

        config_nosni = configparser.ConfigParser()
        config_nosni.read('Icepack_simulations/IMB1_Bl99_noFlooding/namelist_Icepack.ini')

    if krun == 1:

        BuoyData = IMB2Data
        Rtrvl = Rtrvl_IMB2

        config = configparser.ConfigParser()
        config.read('Icepack_simulations/IMB2_Bl99_ctrl/namelist_Icepack.ini')

        config_nosni = configparser.ConfigParser()
        config_nosni.read('Icepack_simulations/IMB2_Bl99_noFlooding/namelist_Icepack.ini')

    #--------------------------------------------------
    # Load BL99 run with snowice:
    MetaData = IcepackDatasetup(config['Icepack_setup'],config_Labels = config['Diag_infos'])
    Data_sni = IcepackData(meta = MetaData)

    #--------------------------------------------------
    # Load BL99 run without snowice:
    MetaData = IcepackDatasetup(config_nosni['Icepack_setup'],config_Labels = config_nosni['Diag_infos'])
    Data_nosni = IcepackData(meta = MetaData)

    visuals.Compute_Errors(Data_buoy = BuoyData,
                           Data_model = Data_nosni,
                           variable = Data_nosni.Hiavg*100.0,
                           reference = Rtrvl.hi_int,
                           krun = krun,
                           OutputFolder = Data_nosni.OutputFolder + 'IMB%s_nosni_Hi' % (krun+1))

    visuals.Compute_Errors(Data_buoy = BuoyData,
                           Data_model = Data_nosni,
                           variable = Data_nosni.Hsavg*100.0,
                           reference = Rtrvl.hs_int,
                           krun = krun,
                           OutputFolder = Data_nosni.OutputFolder + 'IMB%s_nosni_Hs' % (krun+1))

    visuals.Compute_Errors(Data_buoy = BuoyData,
                           Data_model = Data_nosni,
                           variable = Data_nosni.CongTotal*100.0,
                           reference = Rtrvl.congelation,
                           krun = krun,
                           OutputFolder = Data_nosni.OutputFolder + 'IMB%s_nosni_cong' % (krun+1))

    visuals.Compute_Errors(Data_buoy = BuoyData,
                           Data_model = Data_nosni,
                           variable = Data_nosni.CumulSnowice*100.0,
                           reference = Rtrvl.snowice,
                           krun = krun,
                           OutputFolder = Data_nosni.OutputFolder + 'IMB%s_nosni_snowice' % (krun+1))


    visuals.F7_thickness_models(data_buoy = BuoyData,
                                    rtrvl = Rtrvl,
                                    krun = krun,
                                    data_1 = Data_sni,
                                    data_2 = Data_nosni,
                                    Figure = 7,
                                    OutputFolder=OutputFolder )
del visuals

#Figure 8: Mushy, with and without snow ice
visuals = vis()
for krun in range(0, 2):

    if krun == 0:

        BuoyData = IMB1Data
        Rtrvl = Rtrvl_IMB1

        config = configparser.ConfigParser()
        config.read('Icepack_simulations/IMB1_mushy_ctrl/namelist_Icepack.ini')

        config_nosni = configparser.ConfigParser()
        config_nosni.read('Icepack_simulations/IMB1_mushy_noFlooding/namelist_Icepack.ini')

    if krun == 1:

        BuoyData = IMB2Data
        Rtrvl = Rtrvl_IMB2

        config = configparser.ConfigParser()
        config.read('Icepack_simulations/IMB2_mushy_ctrl/namelist_Icepack.ini')

        config_nosni = configparser.ConfigParser()
        config_nosni.read('Icepack_simulations/IMB2_mushy_noFlooding/namelist_Icepack.ini')

    #--------------------------------------------------
    # Icepack data: with snow-ice
    MetaData = IcepackDatasetup(config['Icepack_setup'],config_Labels = config['Diag_infos'])
    Data_sni = IcepackData(meta = MetaData)

    #--------------------------------------------------
    # Icepack data: without snow-ice
    MetaData = IcepackDatasetup(config_nosni['Icepack_setup'],config_Labels = config_nosni['Diag_infos'])
    Data_nosni = IcepackData(meta = MetaData)


    visuals.Compute_Errors(Data_buoy = BuoyData,
                           Data_model = Data_nosni,
                           variable = Data_nosni.Hiavg*100.0,
                           reference = Rtrvl.hi_int,
                           krun = krun,
                           OutputFolder = Data_nosni.OutputFolder + 'IMB%s_nosni_Hi' % (krun+1))

    visuals.Compute_Errors(Data_buoy = BuoyData,
                           Data_model = Data_nosni,
                           variable = Data_nosni.Hsavg*100.0,
                           reference = Rtrvl.hs_int,
                           krun = krun,
                           OutputFolder = Data_nosni.OutputFolder + 'IMB%s_nosni_Hs' % (krun+1))

    visuals.Compute_Errors(Data_buoy = BuoyData,
                           Data_model = Data_nosni,
                           variable = Data_nosni.CongTotal*100.0,
                           reference = Rtrvl.congelation,
                           krun = krun,
                           OutputFolder = Data_nosni.OutputFolder + 'IMB%s_nosni_cong' % (krun+1))

    visuals.Compute_Errors(Data_buoy = BuoyData,
                           Data_model = Data_nosni,
                           variable = Data_nosni.CumulSnowice*100.0,
                           reference = Rtrvl.snowice,
                           krun = krun,
                           OutputFolder = Data_nosni.OutputFolder + 'IMB%s_nosni_snowice' % (krun+1))

    visuals.F7_thickness_models(data_buoy = BuoyData,
                                    rtrvl = Rtrvl,
                                    krun = krun,
                                    data_1 = Data_sni,
                                    data_2 = Data_nosni,
                                    Figure = 8,
                                    OutputFolder=OutputFolder )
del visuals

#Figure 9-10
#Adding a porosity criteria and rate to the flooding parameterization.
visuals = vis()
for krun in range(0, 2):

    if krun == 0:

        BuoyData = IMB1Data
        Rtrvl = Rtrvl_IMB1

        config_ctrl = configparser.ConfigParser()
        config_ctrl.read('Icepack_simulations/IMB1_mushy_noFlooding/namelist_Icepack.ini')

        config_phi = configparser.ConfigParser()
        config_phi.read('Icepack_simulations/IMB1_mushy_FloodCriteria_phi0.005/namelist_Icepack.ini')

        config_manual = configparser.ConfigParser()
        config_manual.read('Icepack_simulations/IMB1_mushy_Manual_Flooding/namelist_Icepack.ini')

    if krun == 1:

        BuoyData = IMB2Data
        Rtrvl = Rtrvl_IMB2

        config_ctrl = configparser.ConfigParser()
        config_ctrl.read('Icepack_simulations/IMB2_mushy_noFlooding/namelist_Icepack.ini')

        config_phi = configparser.ConfigParser()
        config_phi.read('Icepack_simulations/IMB2_mushy_FloodCriteria_phi0.005/namelist_Icepack.ini')

        config_manual = configparser.ConfigParser()
        config_manual.read('Icepack_simulations/IMB2_mushy_Manual_Flooding/namelist_Icepack.ini')

    #--------------------------------------------------
    # Icepack data: ctrl run
    MetaData = IcepackDatasetup(config_ctrl['Icepack_setup'],config_Labels = config_ctrl['Diag_infos'])
    Data_ctrl = IcepackData(meta = MetaData)
    Data_ctrl.Interpolate_into_IMB_profile(nsensors=Rtrvl.nsensors-1,meta = MetaData )


    #--------------------------------------------------
    # Icepack data: phi = 0.005 criteria
    MetaData = IcepackDatasetup(config_phi['Icepack_setup'],config_Labels = config_phi['Diag_infos'])
    Data_phi = IcepackData(meta = MetaData)
    Data_phi.Interpolate_into_IMB_profile(nsensors=Rtrvl.nsensors-1,meta = MetaData )

    #--------------------------------------------------
    # Icepack data: Manual flooding onset
    MetaData = IcepackDatasetup(config_manual['Icepack_setup'],config_Labels = config_manual['Diag_infos'])
    Data_manual = IcepackData(meta = MetaData)
    Data_manual.Interpolate_into_IMB_profile(nsensors=Rtrvl.nsensors-1,meta = MetaData)


    visuals.Compute_Errors(Data_buoy = BuoyData,
                           Data_model = Data_phi,
                           variable = Data_phi.Hiavg*100.0,
                           reference = Rtrvl.hi_int,
                           krun = krun,
                           OutputFolder = Data_phi.OutputFolder + 'IMB%s_phi005_Hi' % (krun+1))

    visuals.Compute_Errors(Data_buoy = BuoyData,
                           Data_model = Data_phi,
                           variable = Data_phi.Hsavg*100.0,
                           reference = Rtrvl.hs_int,
                           krun = krun,
                           OutputFolder = Data_phi.OutputFolder + 'IMB%s_phi005_Hs' % (krun+1))

    visuals.Compute_Errors(Data_buoy = BuoyData,
                           Data_model = Data_phi,
                           variable = Data_phi.CongTotal*100.0,
                           reference = Rtrvl.congelation,
                           krun = krun,
                           OutputFolder = Data_phi.OutputFolder + 'IMB%s_phi005_cong' % (krun+1))

    visuals.Compute_Errors(Data_buoy = BuoyData,
                           Data_model = Data_phi,
                           variable = Data_phi.CumulSnowice*100.0,
                           reference = Rtrvl.snowice,
                           krun = krun,
                           OutputFolder = Data_phi.OutputFolder + 'IMB%s_phi005_snowice' % (krun+1))

    visuals.Compute_Errors(Data_buoy = BuoyData,
                           Data_model = Data_manual,
                           variable = Data_manual.Hiavg*100.0,
                           reference = Rtrvl.hi_int,
                           krun = krun,
                           OutputFolder = Data_manual.OutputFolder + 'IMB%s_manual_Hi' % (krun+1))

    visuals.Compute_Errors(Data_buoy = BuoyData,
                           Data_model = Data_manual,
                           variable = Data_manual.Hsavg*100.0,
                           reference = Rtrvl.hs_int,
                           krun = krun,
                           OutputFolder = Data_manual.OutputFolder + 'IMB%s_manual_Hs' % (krun+1))

    visuals.Compute_Errors(Data_buoy = BuoyData,
                           Data_model = Data_manual,
                           variable = Data_manual.CongTotal*100.0,
                           reference = Rtrvl.congelation,
                           krun = krun,
                           OutputFolder = Data_manual.OutputFolder + 'IMB%s_manual_cong' % (krun+1))

    visuals.Compute_Errors(Data_buoy = BuoyData,
                           Data_model = Data_manual,
                           variable = Data_manual.CumulSnowice*100.0,
                           reference = Rtrvl.snowice,
                           krun = krun,
                           OutputFolder = Data_manual.OutputFolder + 'IMB%s_manual_snowice' % (krun+1))

    visuals.Figure9_temperature_intercomp(Data_1 = Data_ctrl,
                                Data_2 = Data_phi,
                                Data_3 = Data_manual,
                                Rtrvl= Rtrvl,
                                krun = krun,
                                OutputFolder = OutputFolder + 'Figure9_IMB%s_' % (krun+1))

    visuals.figure_layer_thermodynamics(Data_buoy = BuoyData,
                                Data_1 = Data_ctrl,
                                Data_2 = Data_phi,
                                Data_3 = Data_manual,
                                N = 0,
                                krun = krun,
                                OutputFolder = OutputFolder + 'Figure10_IMB%s_' % (krun+1))

#Figure 11
#Comparison F: Testing the brine drainage parameters, in standard code: w_brine and phi_i_mushy
visuals = vis() 
for krun in range(0, 2):

    if krun == 0:
		
        BuoyData = IMB1Data
        Rtrvl = Rtrvl_IMB1
        
        config_ctrl = configparser.ConfigParser()
        config_ctrl.read('Icepack_simulations/IMB1_mushy_noFlooding/namelist_Icepack.ini')
        
        config_2 = configparser.ConfigParser()
        config_2.read('Icepack_simulations/IMB1_mushy_w0.002/namelist_Icepack.ini')
        
        config_3 = configparser.ConfigParser()
        config_3.read('Icepack_simulations/IMB1_mushy_w0.001/namelist_Icepack.ini')
           

    elif krun == 1:
		
        BuoyData = IMB1Data
        Rtrvl = Rtrvl_IMB1
        
        config_ctrl = configparser.ConfigParser()
        config_ctrl.read('Icepack_simulations/IMB1_mushy_noFlooding/namelist_Icepack.ini')
        
        config_2 = configparser.ConfigParser()
        config_2.read('Icepack_simulations/IMB1_mushy_Modified_phi_init_0.65/namelist_Icepack.ini')
        
        config_3= configparser.ConfigParser()
        config_3.read('Icepack_simulations/IMB1_mushy_Modified_phi_init_0.45/namelist_Icepack.ini')
           
    #--------------------------------------------------
    # Icepack Ctrl run:
    MetaData = IcepackDatasetup(config_ctrl['Icepack_setup'],config_Labels = config_ctrl['Diag_infos'])
    Data_ctrl = IcepackData(meta = MetaData)


    #--------------------------------------------------
    # Icepack second run, lower salinity:
    MetaData = IcepackDatasetup(config_2['Icepack_setup'],config_Labels = config_2['Diag_infos'])
    Data_2 = IcepackData(meta = MetaData)


    #--------------------------------------------------
    # Icepack third run, lowest salinity:
    MetaData = IcepackDatasetup(config_3['Icepack_setup'],config_Labels = config_3['Diag_infos'])
    Data_3 = IcepackData(meta = MetaData)
  
    visuals.Compute_Errors(Data_buoy = BuoyData,  
                           Data_model = Data_2,  
                           variable = Data_2.Hiavg*100.0,
                           reference = Rtrvl.hi_int,
                           krun = krun,
                           OutputFolder = Data_2.OutputFolder + 'IMB1_Exp%s_sim2_Hi' % (krun+1))

    visuals.Compute_Errors(Data_buoy = BuoyData,  
                           Data_model = Data_2,  
                           variable = Data_2.Hsavg*100.0,
                           reference = Rtrvl.hs_int,
                           krun = krun,
                           OutputFolder = Data_2.OutputFolder + 'IMB1_Exp%s_sim2_Hs' % (krun+1))    
    
    visuals.Compute_Errors(Data_buoy = BuoyData,  
                           Data_model = Data_2,  
                           variable = Data_2.CongTotal*100.0,
                           reference = Rtrvl.congelation,
                           krun = krun,
                           OutputFolder = Data_2.OutputFolder + 'IMB1_Exp%s_sim2_cong' % (krun+1))    

    visuals.Compute_Errors(Data_buoy = BuoyData,  
                           Data_model = Data_2,  
                           variable = Data_2.CumulSnowice*100.0,
                           reference = Rtrvl.snowice,
                           krun = krun,
                           OutputFolder = Data_2.OutputFolder + 'IMB1_Exp%s_sim2_snowice' % (krun+1))    

    visuals.Compute_Errors(Data_buoy = BuoyData,  
                           Data_model = Data_3,  
                           variable = Data_3.Hiavg*100.0,
                           reference = Rtrvl.hi_int,
                           krun = krun,
                           OutputFolder = Data_3.OutputFolder + 'IMB1_Exp%s_sim3_Hi' % (krun+1))

    visuals.Compute_Errors(Data_buoy = BuoyData,  
                           Data_model = Data_3,  
                           variable = Data_3.Hsavg*100.0,
                           reference = Rtrvl.hs_int,
                           krun = krun,
                           OutputFolder = Data_3.OutputFolder + 'IMB1_Exp%s_sim3_Hs' % (krun+1))    
    
    visuals.Compute_Errors(Data_buoy = BuoyData,  
                           Data_model = Data_3,  
                           variable = Data_3.CongTotal*100.0,
                           reference = Rtrvl.congelation,
                           krun = krun,
                           OutputFolder = Data_3.OutputFolder + 'IMB1_Exp%s_sim3_cong' % (krun+1))    

    visuals.Compute_Errors(Data_buoy = BuoyData,  
                           Data_model = Data_3,  
                           variable = Data_3.CumulSnowice*100.0,
                           reference = Rtrvl.snowice,
                           krun = krun,
                           OutputFolder = Data_3.OutputFolder + 'IMB1_Exp%s_sim3_snowice' % (krun+1))    
                                                                                  
    visuals.Figure11_bottomlayer_thermodynamics(Data_buoy = BuoyData,
                                Data_1 = Data_ctrl,
                                Data_2 = Data_2, 
                                Data_3 = Data_3,
                                krun = krun,
                                OutputFolder = OutputFolder + 'Figure11_') 

del visuals

#Figure 12
#Results with mushy simulation with snow ice and congelation tuned to best represent the observations
visuals = vis() 
for krun in range(0, 2):

    if krun == 0:
		
        BuoyData = IMB1Data
        Rtrvl = Rtrvl_IMB1
        
        config_mushy = configparser.ConfigParser()
        config_mushy.read('Icepack_simulations/IMB1_mushy_bestFit/namelist_Icepack.ini')

        site = "IMB1"

    if krun == 1:
		
        BuoyData = IMB2Data
        Rtrvl = Rtrvl_IMB2
        
        config_mushy = configparser.ConfigParser()
        config_mushy.read('Icepack_simulations/IMB2_mushy_bestFit/namelist_Icepack.ini')

        site = "IMB2"
    
    #--------------------------------------------------
    # Icepack data: BestFit mushy simulation
    
    mushy_MetaData = IcepackDatasetup(config_mushy['Icepack_setup'],config_Labels = config_mushy['Diag_infos'])
    Data_mushy = IcepackData(meta = mushy_MetaData)

    visuals.Compute_Errors(Data_buoy = BuoyData,  
                           Data_model = Data_mushy,  
                           variable = Data_mushy.Hiavg*100.0,
                           reference = Rtrvl.hi_int,
                           krun = krun,
                           OutputFolder = Data_mushy.OutputFolder + 'IMB%s_bestFit_Hi' % (krun+1))

    visuals.Compute_Errors(Data_buoy = BuoyData,  
                           Data_model = Data_mushy,  
                           variable = Data_mushy.Hsavg*100.0,
                           reference = Rtrvl.hs_int,
                           krun = krun,
                           OutputFolder = Data_mushy.OutputFolder + 'IMB%s_bestFit_Hs' % (krun+1))
    
   
    visuals.Compute_Errors(Data_buoy = BuoyData,  
                           Data_model = Data_mushy,  
                           variable = Data_mushy.CongTotal*100.0,
                           reference = Rtrvl.congelation,
                           krun = krun,
                           OutputFolder = Data_mushy.OutputFolder + 'IMB%s_bestFit_cong' % (krun+1))
                            
    visuals.Compute_Errors(Data_buoy = BuoyData,  
                           Data_model = Data_mushy,  
                           variable = Data_mushy.CumulSnowice*100.0,
                           reference = Rtrvl.snowice,
                           krun = krun,
                           OutputFolder = Data_mushy.OutputFolder + 'IMB%s_bestFit_snowice' % (krun+1))
                                

    visuals.F7_thickness_models(data_buoy = BuoyData, 
                                    rtrvl = Rtrvl, 
                                    krun = krun, 
                                    data_1 = Data_mushy, 
                                    Figure = 12,
                                    OutputFolder=OutputFolder)

del visuals

#Figure Appendix 1
#Testing the influence of the Phi Mushy initial parameter, in the original code (as in A)
visuals = vis()
for krun in range(0, 1):

    if krun == 0:

        BuoyData = IMB1Data
        Rtrvl = Rtrvl_IMB1

        config_ctrl = configparser.ConfigParser()
        config_ctrl.read('Icepack_simulations/IMB1_mushy_noFlooding/namelist_Icepack.ini')

        config_phi = configparser.ConfigParser()
        config_phi.read('Icepack_simulations/IMB1_mushy_noFlooding_phi_init_0.65/namelist_Icepack.ini')

    #--------------------------------------------------
    # Icepack Ctrl simulation:
    MetaData = IcepackDatasetup(config_ctrl['Icepack_setup'],config_Labels = config_ctrl['Diag_infos'])
    Data_ctrl = IcepackData(meta = MetaData)

    #--------------------------------------------------
    # Icepack simulation with phi_i = 0.65:
    MetaData = IcepackDatasetup(config_phi['Icepack_setup'],config_Labels = config_phi['Diag_infos'])
    Data_phi65 = IcepackData(meta = MetaData)

    visuals.figure_congelation_intercomp(Data_buoy = BuoyData,
                                         Rtrvl = Rtrvl,
                                         Data_1 = Data_ctrl,
                                         Data_2 = Data_phi65,
                                         krun = krun,
                                         Figure = 1,
                                         OutputFolder=OutputFolder + 'Figure_Appendix1_')

del visuals

#Figure Appendix 2
#Testing the new method for congelation, integrating a mushy layer in the lowest
#              ice layer, instead of only water and then frazil.
visuals = vis()
for krun in range(0, 1):

    if krun == 0:

        BuoyData = IMB1Data
        Rtrvl = Rtrvl_IMB1

        config_ctrl = configparser.ConfigParser()
        config_ctrl.read('Icepack_simulations/IMB1_mushy_noFlooding/namelist_Icepack.ini')

        config_newCong = configparser.ConfigParser()
        config_newCong.read('Icepack_simulations/IMB1_mushy_Modified_phi_init_0.85/namelist_Icepack.ini')

    #--------------------------------------------------
    # Icepack ctrl simulation:
    MetaData = IcepackDatasetup(config_ctrl['Icepack_setup'],config_Labels = config_ctrl['Diag_infos'])
    Data_ctrl = IcepackData(meta = MetaData)

    #--------------------------------------------------
    # Icepack simulation, new congelation scheme:
    MetaData = IcepackDatasetup(config_newCong['Icepack_setup'],config_Labels = config_newCong['Diag_infos'])
    Data_newCong = IcepackData(meta = MetaData)

    visuals.figure_congelation_intercomp(Data_buoy = BuoyData,
                                         Rtrvl = Rtrvl,
                                         Data_1 = Data_ctrl,
                                         Data_2 = Data_newCong,
                                         krun = krun,
                                         Figure = 2,
                                         OutputFolder=OutputFolder + 'Figure_Appendix2_')

del visuals

#Figure Appendix 3
#Testing the influence of the Phi Mushy initial parameter, in the modified congelation scheme
visuals = vis()
for krun in range(0, 2):

    if krun == 0:

        BuoyData = IMB1Data
        Rtrvl = Rtrvl_IMB1

        config_phi85 = configparser.ConfigParser()
        config_phi85.read('Icepack_simulations/IMB1_mushy_Modified_phi_init_0.85/namelist_Icepack.ini')

        config_phi65 = configparser.ConfigParser()
        config_phi65.read('Icepack_simulations/IMB1_mushy_Modified_phi_init_0.65/namelist_Icepack.ini')

        config_phi45 = configparser.ConfigParser()
        config_phi45.read('Icepack_simulations/IMB1_mushy_Modified_phi_init_0.45/namelist_Icepack.ini')

    if krun == 1:

        BuoyData = IMB2Data
        Rtrvl = Rtrvl_IMB2

        config_phi85 = configparser.ConfigParser()
        config_phi85.read('Icepack_simulations/IMB2_mushy_Modified_phi_init_0.85/namelist_Icepack.ini')

        config_phi65 = configparser.ConfigParser()
        config_phi65.read('Icepack_simulations/IMB2_mushy_Modified_phi_init_0.65/namelist_Icepack.ini')

        config_phi45 = configparser.ConfigParser()
        config_phi45.read('Icepack_simulations/IMB2_mushy_Modified_phi_init_0.45/namelist_Icepack.ini')

    #--------------------------------------------------
    # Icepack BL99 data:
    MetaData = IcepackDatasetup(config_phi85['Icepack_setup'],config_Labels = config_phi85['Diag_infos'])
    Data_phi85 = IcepackData(meta = MetaData)

    #--------------------------------------------------
    # Icepack mushy data:
    MetaData = IcepackDatasetup(config_phi65['Icepack_setup'],config_Labels = config_phi65['Diag_infos'])
    Data_phi65 = IcepackData(meta = MetaData)

    #--------------------------------------------------
    # No drainage data:
    MetaData = IcepackDatasetup(config_phi45['Icepack_setup'],config_Labels = config_phi45['Diag_infos'])
    Data_phi45 = IcepackData(meta = MetaData)


    visuals.Compute_Errors(Data_buoy = BuoyData,
                           Data_model = Data_phi85,
                           variable = Data_phi85.Hiavg*100.0,
                           reference = Rtrvl.hi_int,
                           krun = krun,
                           OutputFolder = Data_phi85.OutputFolder + 'IMB%s_phi0.85_Hi' % (krun+1))

    visuals.Compute_Errors(Data_buoy = BuoyData,
                           Data_model = Data_phi85,
                           variable = Data_phi85.Hsavg*100.0,
                           reference = Rtrvl.hs_int,
                           krun = krun,
                           OutputFolder = Data_phi85.OutputFolder + 'IMB%s_phi0.85_Hs' % (krun+1))

    visuals.Compute_Errors(Data_buoy = BuoyData,
                           Data_model = Data_phi85,
                           variable = Data_phi85.CongTotal*100.0,
                           reference = Rtrvl.congelation,
                           krun = krun,
                           OutputFolder = Data_phi85.OutputFolder + 'IMB%s_phi0.85_cong' % (krun+1))

    visuals.Compute_Errors(Data_buoy = BuoyData,
                           Data_model = Data_phi85,
                           variable = Data_phi85.CumulSnowice*100.0,
                           reference = Rtrvl.snowice,
                           krun = krun,
                           OutputFolder = Data_phi85.OutputFolder + 'IMB%s_phi0.85_snowice' % (krun+1))

    visuals.Compute_Errors(Data_buoy = BuoyData,
                           Data_model = Data_phi65,
                           variable = Data_phi65.Hiavg*100.0,
                           reference = Rtrvl.hi_int,
                           krun = krun,
                           OutputFolder = Data_phi65.OutputFolder + 'IMB%s_phi0.65_Hi' % (krun+1))

    visuals.Compute_Errors(Data_buoy = BuoyData,
                           Data_model = Data_phi65,
                           variable = Data_phi65.Hsavg*100.0,
                           reference = Rtrvl.hs_int,
                           krun = krun,
                           OutputFolder = Data_phi65.OutputFolder + 'IMB%s_phi0.65_Hs' % (krun+1))

    visuals.Compute_Errors(Data_buoy = BuoyData,
                           Data_model = Data_phi65,
                           variable = Data_phi65.CongTotal*100.0,
                           reference = Rtrvl.congelation,
                           krun = krun,
                           OutputFolder = Data_phi65.OutputFolder + 'IMB%s_phi0.65_cong' % (krun+1))

    visuals.Compute_Errors(Data_buoy = BuoyData,
                           Data_model = Data_phi65,
                           variable = Data_phi65.CumulSnowice*100.0,
                           reference = Rtrvl.snowice,
                           krun = krun,
                           OutputFolder = Data_phi65.OutputFolder + 'IMB%s_phi0.65_snowice' % (krun+1))

    visuals.Compute_Errors(Data_buoy = BuoyData,
                           Data_model = Data_phi45,
                           variable = Data_phi45.Hiavg*100.0,
                           reference = Rtrvl.hi_int,
                           krun = krun,
                           OutputFolder = Data_phi45.OutputFolder + 'IMB%s_phi0.45_Hi' % (krun+1))

    visuals.Compute_Errors(Data_buoy = BuoyData,
                           Data_model = Data_phi45,
                           variable = Data_phi45.Hsavg*100.0,
                           reference = Rtrvl.hs_int,
                           krun = krun,
                           OutputFolder = Data_phi45.OutputFolder + 'IMB%s_phi0.45_Hs' % (krun+1))

    visuals.Compute_Errors(Data_buoy = BuoyData,
                           Data_model = Data_phi45,
                           variable = Data_phi45.CongTotal*100.0,
                           reference = Rtrvl.congelation,
                           krun = krun,
                           OutputFolder = Data_phi45.OutputFolder + 'IMB%s_phi0.45_cong' % (krun+1))

    visuals.Compute_Errors(Data_buoy = BuoyData,
                           Data_model = Data_phi45,
                           variable = Data_phi45.CumulSnowice*100.0,
                           reference = Rtrvl.snowice,
                           krun = krun,
                           OutputFolder = Data_phi45.OutputFolder + 'IMB%s_phi0.45_snowice' % (krun+1))
    if krun == 0:
        visuals.figure_congelation_intercomp(Data_buoy = BuoyData,
                                         Rtrvl = Rtrvl,
                                         Data_1 = Data_phi85,
                                         Data_2 = Data_phi65,
                                         Data_3 =  Data_phi45,
                                         Figure = 3,
                                         krun = krun,
                                         OutputFolder=OutputFolder+ 'Figure_Appendix3_')

del visuals

#-----------------------------------------------------
#The rest below are supplementary runs not included in the analysis
#-----------------------------------------------------

#Temperature profiles, snow-ice experiment
OutputFolder = 'ExtraMaterial/'
visuals = vis()
for krun in range(0, 2):

    if krun == 0: #For the IMB1 site

        BuoyData = IMB1Data
        Rtrvl = Rtrvl_IMB1
        time = timeIMB1
        Data_heat = IMB1_heat

        config = configparser.ConfigParser()
        config.read('Icepack_simulations/IMB1_Bl99_noFlooding/namelist_Icepack.ini')

        config_mushy = configparser.ConfigParser()
        config_mushy.read('Icepack_simulations/IMB1_mushy_noFlooding/namelist_Icepack.ini')


    if krun == 1:  #For the IMB2 site

        BuoyData = IMB2Data
        Rtrvl = Rtrvl_IMB2
        time = timeIMB2
        Data_heat = IMB2_heat

        config = configparser.ConfigParser()
        config.read('Icepack_simulations/IMB2_Bl99_noFlooding/namelist_Icepack.ini')

        config_mushy = configparser.ConfigParser()
        config_mushy.read('Icepack_simulations/IMB2_mushy_noFlooding/namelist_Icepack.ini')

    #--------------------------------------------------
    # Icepack BL99 data:
    Bl99_MetaData = IcepackDatasetup(config['Icepack_setup'],config_Labels = config['Diag_infos'])
    Data_Bl99 = IcepackData(meta = Bl99_MetaData)
    Data_Bl99.Interpolate_into_IMB_profile(nsensors=Rtrvl.nsensors-1,meta = Bl99_MetaData )

    #--------------------------------------------------
    # Icepack mushy data:
    mushy_MetaData = IcepackDatasetup(config_mushy['Icepack_setup'],config_Labels = config_mushy['Diag_infos'])
    Data_mushy = IcepackData(meta = mushy_MetaData)
    Data_mushy.Interpolate_into_IMB_profile(nsensors=Rtrvl.nsensors-1,meta = mushy_MetaData )
 
    # Figure 6
    visuals.figure_temperature_intercomp(Data_buoy = BuoyData,
                                       Data_model = Data_Bl99,
                                       Data_mushy = Data_mushy,
                                       Rtrvl_buoy = Rtrvl,
                                       krun = krun,
                                       OutputFolder = '%s%s' % (OutputFolder, 'Snowice_Experiment'))

del visuals

#Temperature profiles, manual flooding onset
visuals = vis()
for krun in range(0, 2):

    if krun == 0: #For the IMB1 site

        BuoyData = IMB1Data
        Rtrvl = Rtrvl_IMB1
        time = timeIMB1
        Data_heat = IMB1_heat

        config = configparser.ConfigParser()
        config.read('Icepack_simulations/IMB1_Bl99_Manual_Flooding/namelist_Icepack.ini')

        config_mushy = configparser.ConfigParser()
        config_mushy.read('Icepack_simulations/IMB1_mushy_Manual_Flooding/namelist_Icepack.ini')


    if krun == 1:  #For the IMB2 site

        BuoyData = IMB2Data
        Rtrvl = Rtrvl_IMB2
        time = timeIMB2
        Data_heat = IMB2_heat

        config = configparser.ConfigParser()
        config.read('Icepack_simulations/IMB2_Bl99_Manual_Flooding/namelist_Icepack.ini')

        config_mushy = configparser.ConfigParser()
        config_mushy.read('Icepack_simulations/IMB2_mushy_Manual_Flooding/namelist_Icepack.ini')

    #--------------------------------------------------
    # Icepack BL99 data:
    Bl99_MetaData = IcepackDatasetup(config['Icepack_setup'],config_Labels = config['Diag_infos'])
    Data_Bl99 = IcepackData(meta = Bl99_MetaData)
    Data_Bl99.Interpolate_into_IMB_profile(nsensors=Rtrvl.nsensors-1,meta = Bl99_MetaData )

    #--------------------------------------------------
    # Icepack mushy data:
    mushy_MetaData = IcepackDatasetup(config_mushy['Icepack_setup'],config_Labels = config_mushy['Diag_infos'])
    Data_mushy = IcepackData(meta = mushy_MetaData)
    Data_mushy.Interpolate_into_IMB_profile(nsensors=Rtrvl.nsensors-1,meta = mushy_MetaData )
 
    # Figure 6
    visuals.figure_temperature_intercomp(Data_buoy = BuoyData,
                                       Data_model = Data_Bl99,
                                       Data_mushy = Data_mushy,
                                       Rtrvl_buoy = Rtrvl,
                                       krun = krun,
                                       OutputFolder = '%s%s' % (OutputFolder, 'Manual_Onset_Experiment'))

del visuals

#Temperature profiles, phi_min experiment
visuals = vis()
for krun in range(0, 2):

    if krun == 0: #For the IMB1 site

        BuoyData = IMB1Data
        Rtrvl = Rtrvl_IMB1
        time = timeIMB1
        Data_heat = IMB1_heat

        config = configparser.ConfigParser()
        config.read('Icepack_simulations/IMB1_mushy_FloodCriteria_phi0.005/namelist_Icepack.ini')

        config_2 = configparser.ConfigParser()
        config_2.read('Icepack_simulations/IMB1_mushy_FloodCriteria_phi0.006/namelist_Icepack.ini')

        config_3= configparser.ConfigParser()
        config_3.read('Icepack_simulations/IMB1_mushy_FloodCriteria_phi0.007/namelist_Icepack.ini')

    if krun == 1:  #For the IMB2 site

        BuoyData = IMB2Data
        Rtrvl = Rtrvl_IMB2
        time = timeIMB2
        Data_heat = IMB2_heat

        config = configparser.ConfigParser()
        config.read('Icepack_simulations/IMB2_mushy_FloodCriteria_phi0.005/namelist_Icepack.ini')

        config_2 = configparser.ConfigParser()
        config_2.read('Icepack_simulations/IMB2_mushy_FloodCriteria_phi0.006/namelist_Icepack.ini')

        config_3= configparser.ConfigParser()
        config_3.read('Icepack_simulations/IMB2_mushy_FloodCriteria_phi0.007/namelist_Icepack.ini')
        
    #--------------------------------------------------
    # Icepack data1:
    MetaData = IcepackDatasetup(config['Icepack_setup'],config_Labels = config['Diag_infos'])
    Data_1 = IcepackData(meta = MetaData)
    Data_1.Interpolate_into_IMB_profile(nsensors=Rtrvl.nsensors-1,meta = MetaData )

    #--------------------------------------------------
    # Icepack data2:
    MetaData2 = IcepackDatasetup(config_2['Icepack_setup'],config_Labels = config_2['Diag_infos'])
    Data_2 = IcepackData(meta = MetaData2)
    Data_2.Interpolate_into_IMB_profile(nsensors=Rtrvl.nsensors-1,meta = MetaData2 )
 
    #--------------------------------------------------
    # Icepack data3:
    MetaData3 = IcepackDatasetup(config_3['Icepack_setup'],config_Labels = config_3['Diag_infos'])
    Data_3 = IcepackData(meta = MetaData3)
    Data_3.Interpolate_into_IMB_profile(nsensors=Rtrvl.nsensors-1,meta = MetaData3 )
 
    # Figure 6
    visuals.figure_temperature_intercomp_supp(Data_buoy = BuoyData,
                                       Data_model = Data_1,
                                       Data_mushy = Data_2,
                                       Data_3 = Data_3,
                                       Rtrvl_buoy = Rtrvl,
                                       krun = krun,
                                       OutputFolder = '%s%s' % (OutputFolder, 'Phi_min_Experiment'))

del visuals

#Temperature profiles, brine drainage experiment
visuals = vis()
for krun in range(0, 2):

    if krun == 0: #For the IMB1 site

        BuoyData = IMB1Data
        Rtrvl = Rtrvl_IMB1
        time = timeIMB1
        Data_heat = IMB1_heat

        config = configparser.ConfigParser()
        config.read('Icepack_simulations/IMB1_mushy_noFlooding/namelist_Icepack.ini')

        config_2 = configparser.ConfigParser()
        config_2.read('Icepack_simulations/IMB1_mushy_w0.002/namelist_Icepack.ini')

        config_3= configparser.ConfigParser()
        config_3.read('Icepack_simulations/IMB1_mushy_w0.001/namelist_Icepack.ini')

    if krun == 1:  #For the IMB2 site

        BuoyData = IMB2Data
        Rtrvl = Rtrvl_IMB2
        time = timeIMB2
        Data_heat = IMB2_heat

        config = configparser.ConfigParser()
        config.read('Icepack_simulations/IMB2_mushy_noFlooding/namelist_Icepack.ini')

        config_2 = configparser.ConfigParser()
        config_2.read('Icepack_simulations/IMB2_mushy_w0.002/namelist_Icepack.ini')

        config_3= configparser.ConfigParser()
        config_3.read('Icepack_simulations/IMB2_mushy_w0.001/namelist_Icepack.ini')
        
    #--------------------------------------------------
    # Icepack data1:
    MetaData = IcepackDatasetup(config['Icepack_setup'],config_Labels = config['Diag_infos'])
    Data_1 = IcepackData(meta = MetaData)
    Data_1.Interpolate_into_IMB_profile(nsensors=Rtrvl.nsensors-1,meta = MetaData )

    #--------------------------------------------------
    # Icepack data2:
    MetaData2 = IcepackDatasetup(config_2['Icepack_setup'],config_Labels = config_2['Diag_infos'])
    Data_2 = IcepackData(meta = MetaData2)
    Data_2.Interpolate_into_IMB_profile(nsensors=Rtrvl.nsensors-1,meta = MetaData2 )
 
    #--------------------------------------------------
    # Icepack data3:
    MetaData3 = IcepackDatasetup(config_3['Icepack_setup'],config_Labels = config_3['Diag_infos'])
    Data_3 = IcepackData(meta = MetaData3)
    Data_3.Interpolate_into_IMB_profile(nsensors=Rtrvl.nsensors-1,meta = MetaData3 )
 
    # Figure 6
    visuals.figure_temperature_intercomp_supp(Data_buoy = BuoyData,
                                       Data_model = Data_1,
                                       Data_mushy = Data_2,
                                       Data_3 = Data_3,
                                       Rtrvl_buoy = Rtrvl,
                                       krun = krun,
                                       OutputFolder = '%s%s' % (OutputFolder, 'Brine_Experiment'))

del visuals

#Temperature profiles, phi_init experiment, modified congelation.
visuals = vis()
for krun in range(0, 2):

    if krun == 0: #For the IMB1 site

        BuoyData = IMB1Data
        Rtrvl = Rtrvl_IMB1
        time = timeIMB1
        Data_heat = IMB1_heat

        config = configparser.ConfigParser()
        config.read('Icepack_simulations/IMB1_mushy_Modified_phi_init_0.85/namelist_Icepack.ini')

        config_2 = configparser.ConfigParser()
        config_2.read('Icepack_simulations/IMB1_mushy_Modified_phi_init_0.65/namelist_Icepack.ini')

        config_3= configparser.ConfigParser()
        config_3.read('Icepack_simulations/IMB1_mushy_Modified_phi_init_0.45/namelist_Icepack.ini')

    if krun == 1:  #For the IMB2 site

        BuoyData = IMB2Data
        Rtrvl = Rtrvl_IMB2
        time = timeIMB2
        Data_heat = IMB2_heat

        config = configparser.ConfigParser()
        config.read('Icepack_simulations/IMB2_mushy_Modified_phi_init_0.85/namelist_Icepack.ini')

        config_2 = configparser.ConfigParser()
        config_2.read('Icepack_simulations/IMB2_mushy_Modified_phi_init_0.65/namelist_Icepack.ini')

        config_3= configparser.ConfigParser()
        config_3.read('Icepack_simulations/IMB2_mushy_Modified_phi_init_0.45/namelist_Icepack.ini')
        
    #--------------------------------------------------
    # Icepack data1:
    MetaData = IcepackDatasetup(config['Icepack_setup'],config_Labels = config['Diag_infos'])
    Data_1 = IcepackData(meta = MetaData)
    Data_1.Interpolate_into_IMB_profile(nsensors=Rtrvl.nsensors-1,meta = MetaData )

    #--------------------------------------------------
    # Icepack data2:
    MetaData2 = IcepackDatasetup(config_2['Icepack_setup'],config_Labels = config_2['Diag_infos'])
    Data_2 = IcepackData(meta = MetaData2)
    Data_2.Interpolate_into_IMB_profile(nsensors=Rtrvl.nsensors-1,meta = MetaData2 )
 
    #--------------------------------------------------
    # Icepack data3:
    MetaData3 = IcepackDatasetup(config_3['Icepack_setup'],config_Labels = config_3['Diag_infos'])
    Data_3 = IcepackData(meta = MetaData3)
    Data_3.Interpolate_into_IMB_profile(nsensors=Rtrvl.nsensors-1,meta = MetaData3 )
 
    # Figure 6
    visuals.figure_temperature_intercomp_supp(Data_buoy = BuoyData,
                                       Data_model = Data_1,
                                       Data_mushy = Data_2,
                                       Data_3 = Data_3,
                                       Rtrvl_buoy = Rtrvl,
                                       krun = krun,
                                       OutputFolder = '%s%s' % (OutputFolder, 'Phi_init_Experiment'))
del visuals

OutputFolder = 'Outputs/'
#Other porosity criteria simulations.
visuals = vis()
for krun in range(0, 2):

    if krun == 0:

        BuoyData = IMB1Data
        Rtrvl = Rtrvl_IMB1

        config = configparser.ConfigParser()
        config.read('Icepack_simulations/IMB1_Bl99_Manual_Flooding/namelist_Icepack.ini')

        config_phi006 = configparser.ConfigParser()
        config_phi006.read('Icepack_simulations/IMB1_mushy_FloodCriteria_phi0.006/namelist_Icepack.ini')

        config_phi007 = configparser.ConfigParser()
        config_phi007.read('Icepack_simulations/IMB1_mushy_FloodCriteria_phi0.007/namelist_Icepack.ini')

    if krun == 1:
        BuoyData = IMB2Data
        Rtrvl = Rtrvl_IMB2

        config = configparser.ConfigParser()
        config.read('Icepack_simulations/IMB2_Bl99_Manual_Flooding/namelist_Icepack.ini')

        config_phi006 = configparser.ConfigParser()
        config_phi006.read('Icepack_simulations/IMB2_mushy_FloodCriteria_phi0.006/namelist_Icepack.ini')

        config_phi007 = configparser.ConfigParser()
        config_phi007.read('Icepack_simulations/IMB2_mushy_FloodCriteria_phi0.007/namelist_Icepack.ini')

    #--------------------------------------------------
    # Icepack BL99 data:
    MetaData = IcepackDatasetup(config['Icepack_setup'],config_Labels = config['Diag_infos'])
    Data_Bl99manual = IcepackData(meta = MetaData)


    #--------------------------------------------------
    # Icepack mushy data:
    MetaData = IcepackDatasetup(config_phi006['Icepack_setup'],config_Labels = config_phi006['Diag_infos'])
    Data_phi006 = IcepackData(meta = MetaData)

    #--------------------------------------------------
    # Icepack phi data:
    MetaData = IcepackDatasetup(config_phi007['Icepack_setup'],config_Labels = config_phi007['Diag_infos'])
    Data_phi007 = IcepackData(meta = MetaData)

    visuals.Compute_Errors(Data_buoy = BuoyData,
                           Data_model = Data_Bl99manual,
                           variable = Data_Bl99manual.Hiavg*100.0,
                           reference = Rtrvl.hi_int,
                           krun = krun,
                           OutputFolder = Data_Bl99manual.OutputFolder + 'IMB%s_manual_Hi' % (krun+1))

    visuals.Compute_Errors(Data_buoy = BuoyData,
                           Data_model = Data_Bl99manual,
                           variable = Data_Bl99manual.Hsavg*100.0,
                           reference = Rtrvl.hs_int,
                           krun = krun,
                           OutputFolder = Data_Bl99manual.OutputFolder + 'IMB%s_manual_Hs' % (krun+1))


    visuals.Compute_Errors(Data_buoy = BuoyData,
                           Data_model = Data_Bl99manual,
                           variable = Data_Bl99manual.CongTotal*100.0,
                           reference = Rtrvl.congelation,
                           krun = krun,
                           OutputFolder = Data_Bl99manual.OutputFolder + 'IMB%s_manual_cong' % (krun+1))

    visuals.Compute_Errors(Data_buoy = BuoyData,
                           Data_model = Data_Bl99manual,
                           variable = Data_Bl99manual.CumulSnowice*100.0,
                           reference = Rtrvl.snowice,
                           krun = krun,
                           OutputFolder = Data_Bl99manual.OutputFolder + 'IMB%s_manual_snowice' % (krun+1))

    visuals.Compute_Errors(Data_buoy = BuoyData,
                           Data_model = Data_phi006,
                           variable = Data_phi006.Hiavg*100.0,
                           reference = Rtrvl.hi_int,
                           krun = krun,
                           OutputFolder = Data_phi006.OutputFolder + 'IMB%s_phi006_Hi' % (krun+1))

    visuals.Compute_Errors(Data_buoy = BuoyData,
                           Data_model = Data_phi006,
                           variable = Data_phi006.Hsavg*100.0,
                           reference = Rtrvl.hs_int,
                           krun = krun,
                           OutputFolder = Data_phi006.OutputFolder + 'IMB%s_phi006_Hs' % (krun+1))

    visuals.Compute_Errors(Data_buoy = BuoyData,
                           Data_model = Data_phi006,
                           variable = Data_phi006.CongTotal*100.0,
                           reference = Rtrvl.congelation,
                           krun = krun,
                           OutputFolder = Data_phi006.OutputFolder + 'IMB%s_phi006_cong' % (krun+1))

    visuals.Compute_Errors(Data_buoy = BuoyData,
                           Data_model = Data_phi006,
                           variable = Data_phi006.CumulSnowice*100.0,
                           reference = Rtrvl.snowice,
                           krun = krun,
                           OutputFolder = Data_phi006.OutputFolder + 'IMB%s_phi006_snowice' % (krun+1))

    visuals.Compute_Errors(Data_buoy = BuoyData,
                           Data_model = Data_phi007,
                           variable = Data_phi007.Hiavg*100.0,
                           reference = Rtrvl.hi_int,
                           krun = krun,
                           OutputFolder = Data_phi007.OutputFolder + 'IMB%s_phi0.007_Hi' % (krun+1))

    visuals.Compute_Errors(Data_buoy = BuoyData,
                           Data_model = Data_phi007,
                           variable = Data_phi007.Hsavg*100.0,
                           reference = Rtrvl.hs_int,
                           krun = krun,
                           OutputFolder = Data_phi007.OutputFolder + 'IMB%s_phi0.007_Hs' % (krun+1))

    visuals.Compute_Errors(Data_buoy = BuoyData,
                           Data_model = Data_phi007,
                           variable = Data_phi007.CongTotal*100.0,
                           reference = Rtrvl.congelation,
                           krun = krun,
                           OutputFolder = Data_phi007.OutputFolder + 'IMB%s_phi0.007_cong' % (krun+1))

    visuals.Compute_Errors(Data_buoy = BuoyData,
                           Data_model = Data_phi007,
                           variable = Data_phi007.CumulSnowice*100.0,
                           reference = Rtrvl.snowice,
                           krun = krun,
                           OutputFolder = Data_phi007.OutputFolder + 'IMB%s_phi0.007_snowice' % (krun+1))

#Figure 11 but for IMB2
visuals = vis()

for krun in range(0, 2):

    if krun == 0:

        BuoyData = IMB2Data
        Rtrvl = Rtrvl_IMB2

        config_ctrl = configparser.ConfigParser()
        config_ctrl.read('Icepack_simulations/IMB2_mushy_noFlooding/namelist_Icepack.ini')

        config_2 = configparser.ConfigParser()
        config_2.read('Icepack_simulations/IMB2_mushy_w0.002/namelist_Icepack.ini')

        config_3 = configparser.ConfigParser()
        config_3.read('Icepack_simulations/IMB2_mushy_w0.001/namelist_Icepack.ini')

    elif krun == 1:

        BuoyData = IMB2Data
        Rtrvl = Rtrvl_IMB2

        config_ctrl = configparser.ConfigParser()
        config_ctrl.read('Icepack_simulations/IMB2_mushy_Modified_phi_init_0.85/namelist_Icepack.ini')

        config_2 = configparser.ConfigParser()
        config_2.read('Icepack_simulations/IMB2_mushy_Modified_phi_init_0.65/namelist_Icepack.ini')

        config_3= configparser.ConfigParser()
        config_3.read('Icepack_simulations/IMB2_mushy_Modified_phi_init_0.45/namelist_Icepack.ini')

    #--------------------------------------------------
    # Icepack Ctrl run:
    MetaData = IcepackDatasetup(config_ctrl['Icepack_setup'],config_Labels = config_ctrl['Diag_infos'])
    Data_ctrl = IcepackData(meta = MetaData)


    #--------------------------------------------------
    # Icepack second run, lower salinity:
    MetaData = IcepackDatasetup(config_2['Icepack_setup'],config_Labels = config_2['Diag_infos'])
    Data_2 = IcepackData(meta = MetaData)


    #--------------------------------------------------
    # Icepack third run, lowest salinity:
    MetaData = IcepackDatasetup(config_3['Icepack_setup'],config_Labels = config_3['Diag_infos'])
    Data_3 = IcepackData(meta = MetaData)
  
 
    visuals.Compute_Errors(Data_buoy = BuoyData,  
                           Data_model = Data_ctrl,  
                           variable = Data_ctrl.Hiavg*100.0,
                           reference = Rtrvl.hi_int,
                           krun = krun,
                           OutputFolder = Data_ctrl.OutputFolder + 'IMB2_Exp%s_ctrl_Hi' % (krun+1))

    visuals.Compute_Errors(Data_buoy = BuoyData,  
                           Data_model = Data_ctrl,  
                           variable = Data_ctrl.Hsavg*100.0,
                           reference = Rtrvl.hs_int,
                           krun = krun,
                           OutputFolder = Data_ctrl.OutputFolder + 'IMB2_Exp%s_ctrl_Hs' % (krun+1))


    visuals.Compute_Errors(Data_buoy = BuoyData,  
                           Data_model = Data_ctrl,  
                           variable = Data_ctrl.CongTotal*100.0,
                           reference = Rtrvl.congelation,
                           krun = krun,
                           OutputFolder = Data_ctrl.OutputFolder + 'IMB2_Exp%s_ctrl_cong' % (krun+1))

    visuals.Compute_Errors(Data_buoy = BuoyData,  
                           Data_model = Data_ctrl,  
                           variable = Data_ctrl.CumulSnowice*100.0,
                           reference = Rtrvl.snowice,
                           krun = krun,
                           OutputFolder = Data_ctrl.OutputFolder + 'IMB2_Exp%s_ctrl_snowice' % (krun+1))

    visuals.Compute_Errors(Data_buoy = BuoyData,  
                           Data_model = Data_2,  
                           variable = Data_2.Hiavg*100.0,
                           reference = Rtrvl.hi_int,
                           krun = krun,
                           OutputFolder = Data_2.OutputFolder + 'IMB2_Exp%s_sim2_Hi' % (krun+1))

    visuals.Compute_Errors(Data_buoy = BuoyData,  
                           Data_model = Data_2,  
                           variable = Data_2.Hsavg*100.0,
                           reference = Rtrvl.hs_int,
                           krun = krun,
                           OutputFolder = Data_2.OutputFolder + 'IMB2_Exp%s_sim2_Hs' % (krun+1))    
    
    visuals.Compute_Errors(Data_buoy = BuoyData,  
                           Data_model = Data_2,  
                           variable = Data_2.CongTotal*100.0,
                           reference = Rtrvl.congelation,
                           krun = krun,
                           OutputFolder = Data_2.OutputFolder + 'IMB2_Exp%s_sim2_cong' % (krun+1))    

    visuals.Compute_Errors(Data_buoy = BuoyData,  
                           Data_model = Data_2,  
                           variable = Data_2.CumulSnowice*100.0,
                           reference = Rtrvl.snowice,
                           krun = krun,
                           OutputFolder = Data_2.OutputFolder + 'IMB2_Exp%s_sim2_snowice' % (krun+1))    

    visuals.Compute_Errors(Data_buoy = BuoyData,  
                           Data_model = Data_3,  
                           variable = Data_3.Hiavg*100.0,
                           reference = Rtrvl.hi_int,
                           krun = krun,
                           OutputFolder = Data_3.OutputFolder + 'IMB2_Exp%s_sim3_Hi' % (krun+1))

    visuals.Compute_Errors(Data_buoy = BuoyData,  
                           Data_model = Data_3,  
                           variable = Data_3.Hsavg*100.0,
                           reference = Rtrvl.hs_int,
                           krun = krun,
                           OutputFolder = Data_3.OutputFolder + 'IMB2_Exp%s_sim3_Hs' % (krun+1))    
    
    visuals.Compute_Errors(Data_buoy = BuoyData,  
                           Data_model = Data_3,  
                           variable = Data_3.CongTotal*100.0,
                           reference = Rtrvl.congelation,
                           krun = krun,
                           OutputFolder = Data_3.OutputFolder + 'IMB2_Exp%s_sim3_cong' % (krun+1)) 

    visuals.Compute_Errors(Data_buoy = BuoyData,  
                           Data_model = Data_3,  
                           variable = Data_3.CumulSnowice*100.0,
                           reference = Rtrvl.snowice,
                           krun = krun,
                           OutputFolder = Data_3.OutputFolder + 'IMB2_Exp%s_sim3_snowice' % (krun+1))

    visuals.Figure11_bottomlayer_thermodynamics(Data_buoy = BuoyData,
                                Data_1 = Data_ctrl,
                                Data_2 = Data_2, 
                                Data_3 = Data_3,
                                krun = krun,
                                OutputFolder = OutputFolder) 

del visuals

print('END OF SCRIPT! Congrats, this ran perfectly! ')
print('Or did it?...')
# End of program
