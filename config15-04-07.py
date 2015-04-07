import numpy as np
import aolUtil


# Define the basic configuration
basicTofConfig = {
    "acqCh": 10, 
    "baselineSubtraction":"early", # 'early', 'roi', 'none' 
    "baselineEnd_us": 0.207, 
    "calib_file": ("/reg/neh/operator_amoopr/"
		    + "amoh5215/psana/tofCalibs/tof_calib_0.json"), 
    "filterMethod": "none", # wavelet, average, wienerDeconv 
    "filterWaveletLevels": 10, 
    "filterWinerSNR": 1,
    "filterWinerResponse":1,
    "filterAverageNumPoints":4,
    "detectorSource": "DetInfo(AmoETOF.0:Acqiris.0)", 
    "t_min_us": 0.20, 
    "t_max_us": 0.26, 
    "t_slice": True
}

lclsConfig = {
	'lcls_photonEnergyA':1,
	'lcls_photonEnergyB':1}
lclsConfig = aolUtil.struct(lclsConfig)

retardationPV = 'AMO:R14:IOC:10:VHS0:CH0:VoltageMeasure'

# Acqiris channel asignment
acqiris_setup = {
        4:['ACQ1', 6],
        3:['ACQ4', 0],
        2:['ACQ4', 1],
        1:['ACQ4', 2],
        0:['ACQ4', 3],
        15:['ACQ4', 4],
        14:['ACQ4', 5],
        13:['ACQ4', 6],
        12:['ACQ4', 7],
        11:['ACQ1', 0],
        10:['ACQ1', 1],
        9:['ACQ1', 2],
        8:['ACQ2', 0],
        7:['ACQ1', 3],
        6:['ACQ2', 1],
        5:['ACQ1', 5],
        }

# ROI for the photo line
# red
time_roi_0_us_common = [0.232, 0.245]
time_roi_0_us = [time_roi_0_us_common]*16
# Background for the photoline
time_roi_2_us_common = [0.23, 0.232]
time_roi_2_us = [time_roi_2_us_common]*16

# ROI for the auger line
# green
time_roi_1_us_common = [0.215, 0.225]
time_roi_1_us = [time_roi_1_us_common]*16

energy_roi_0_eV_common = [40, 60]
energy_roi_0_eV = [energy_roi_0_eV_common]*16



# Make copies of the basic configuration for each of the detectors
tof_config_list = [None] * 16

def makeTofConfigList(online=True):
    global tof_config_list
    for i in range(16):
        tof_config_list[i] = basicTofConfig.copy()
        tof_config_list[i]['calibFile'] = ('/reg/neh/operator/amoopr/'
                + 'amoi0114/psana/tofCalibs/tof{}Calib.json'.format(i+1)) 
        tof_config_list[i]['detectorSource'] = acqiris_setup[i][0]
        tof_config_list[i]['acqCh'] = acqiris_setup[i][1]

        tof_config_list[i]['time_roi_0_us'] = time_roi_0_us[i]
        tof_config_list[i]['time_roi_2_us'] = time_roi_2_us[i]
        tof_config_list[i]['time_roi_1_us'] = time_roi_1_us[i]
        tof_config_list[i]['energy_roi_0_eV'] = energy_roi_0_eV[i]
        #tof_config_list[i]['energy_roi_1_eV'] = energy_roi_1_eV[i]
        #tof_config_list[i]['energy_roi_2_eV'] = energy_roi_2_eV[i]
	
makeTofConfigList(online=True)

minE_eV = 10
maxE_eV = 200
nEnergyBins = 2**9

energyScaleBinLimits = np.linspace(minE_eV, maxE_eV, nEnergyBins + 1)


fitMask = np.array([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15])
#fitMask = np.array([0,1,2,3,4,5,6, 8,9,10,11, 13,14,15])
#fitMask = np.array([0,1,2,3,5,6,7,8,9,10,11,13,14,15])

boolFitMask = np.array([i in fitMask for i in range(16)])
nanFitMask = np.array([1 if b else np.nan for b in boolFitMask])

offlineSource = 'exp=amoh5215:run=216'

# For CookieBox class debugging
domainToDisplay = 'Time'


# Stuff below are used in the debugging of the tofData class

dataSource = 'exp=amoc8114:run=31'
nEvents = 10
useZmq = False
