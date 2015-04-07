import numpy as np
import aolUtil


# Define the basic configuration
basicTofConfig = {
    "acqCh": 10, 
    "baselineSubtraction":"early", # 'early', 'roi', 'none' 
    "baselineEnd_us": 0.05, 
    "calibFile": ("/reg/neh/operator_amoopr/"
		    + "amoh5215/psana/tofCalibs/tof_calib_0.json"), 
    "filterMethod": "none", # wavelet, average, wienerDeconv 
    "filterWaveletLevels": 10, 
    "filterWinerSNR": 1,
    "filterWinerResponse":1,
    "filterAverageNumPoints":4,
    "detectorSource": "DetInfo(AmoETOF.0:Acqiris.0)", 
    "tMin_us": 0.0, 
    "tMax_us": 0.4, 
    "tSlice": True 
}

lclsConfig = {
	'lcls_photonEnergyA':1,
	'lcls_photonEnergyB':1}
lclsConfig = aolUtil.struct(lclsConfig)

retardationPV = 'AMO:R14:IOC:10:VHS0:CH0:VoltageMeasure'

# Acqiris channel asignment
acqiris_setup = {
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
        8:['ACQ1', 3],
        7:['ACQ1', 4],
        6:['ACQ1', 5],
        5:['ACQ1', 6],
        4:['ACQ1', 7]
        }


timeRoi0_us_common = [0.2, 0.21]	#red
timeRoi0_us = [timeRoi0_us_common]*16

timeRoi0Bg_us_common = [0., 0.05]
timeRoi0Bg_us = [timeRoi0Bg_us_common]*16

timeRoi1_us_common = [0.3, 0.4]	#green
timeRoi1_us = [timeRoi1_us_common]*16

energyRoi0_eV_common = [40, 60]
energyRoi0_eV = [energyRoi0_eV_common]*16



# Make copies of the basic configuration for each of the detectors
tofConfigList = [None] * 16

def makeTofConfigList(online=True):
    global tofConfigList
    for i in range(16):
        tofConfigList[i] = basicTofConfig.copy()
	tofConfigList[i]['calibFile'] = ('/reg/neh/operator/amoopr/'
		    + 'amoi0114/psana/tofCalibs/tof{}Calib.json'.format(i+1)) 
	if online:
            tofConfigList[i]['detectorSource'] = acqiris_setup[i][0]
            tofConfigList[i]['acqCh'] = acqiris_setup[i][1]

makeTofConfigList(online=True)

minE_eV = 10
maxE_eV = 500
nEnergyBins = 256

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
