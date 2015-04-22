import sys
import os.path
sys.path.append(
                os.path.dirname(os.path.abspath(__file__)) +
                '/aolPyModules')
from mpi4py import MPI
import arguments
#import matplotlib.pyplot as plt
#plt.ion()
import numpy as np
import time
import platform
from collections import deque
import importlib
import lmfit

import aolUtil
import simplepsana
import cookie_box
import lcls

# Set up the mpi cpmmunication
world = MPI.COMM_WORLD
rank = world.Get_rank()
worldSize = world.Get_size()
workers = worldSize - 1


#############
# Data definitions
s = 0 # Data size tracker

dRank = s
s += 1
dFiducials = s
s += 1
dTime = s
s += 1
dIntRoi0 = slice(s, s+16)
s += 16
dIntRoi1 = slice(s, s+16)
s += 16
dPol = slice(s, s+8)
s += 8
dEnergy = slice(s, s+2)
s += 2
dEL3 = s
s += 1
dFEE = slice(s, s+4)
s += 4
dDeltaK = s
s += 1
dDeltaEnc = slice(s, s+4)
s += 4

dSize = s

dTraces = None
d_energy_trace = None

def connect_to_data_source(args, config, verbose=False):
    # If online
    if not args.offline:
        # make the shared memory string
        dataSource = 'shmem=AMO.0:stop=no'
    else:
        dataSource = ':'.join([config.offline_source, 'idx'])
        if verbose:
            print config

    if verbose:
        # check the host name
        host = platform.node()
        print 'rank {} (on {}) connecting to datasource: {}'.format(
                rank,
                host,
                dataSource)

    return simplepsana.get_data_source(dataSource)


def importConfiguration(args, verbose=False):
    # Import the correct configuration module    
    # Get the path and the name of the config file
    confPath, confFileName = os.path.split(args.configuration)
    # Get the current working directory to be able to get back
    working_dir = os.getcwd()
    if len(confPath) == 0:
        confPath = working_dir
    # Change dir to the config directory
    os.chdir(confPath)
    
    # Print someinformation
    if verbose:
        print 'Loading configuration from directory "{}"'.format(confPath)
        print 'File name is {}'.format(confFileName)
    confName, _ = os.path.splitext(confFileName)
    if verbose:
        print 'Module name is {}'.format(confName)
    
    # Import the configuration moudule
    conf = importlib.import_module(confName)

    # Change back to the working directory
    os.chdir(working_dir)

    # Update the config with some parameters from the command line
    if args.dataSource is not None:
        conf.offline_source = args.dataSource

    # Return the configuration
    return conf


def getDetectorCalibration(verbose=False, fileName=''):
    if fileName == '':
        detCalib = aolUtil.struct({'path':'detCalib',
            'name':'calib'})
        # Get the latest detector callibration values from file
        if not os.path.exists(detCalib.path):
            os.makedirs(detCalib.path)
            np.savetxt(detCalib.path + '/' + detCalib.name + '0.txt', [1]*16)
        
        detCalib.fileNumber = np.max([int(f[len(detCalib.name):-4]) for f in
            os.listdir(detCalib.path) if len(f) > len(detCalib.name) and
            f[:len(detCalib.name)]==detCalib.name])
    else:
        detCalib = aolUtil.struct()
        splitPath = fileName.split('/')
        if len(splitPath) > 1:
            detCalib.path = '/'.join(splitPath[:-1])
        else:
            detCalib.path = '.'
        detCalib.name = splitPath[-1]
        detCalib.fileNumber = np.nan


    if args.calibrate == -1:
        detCalib.factors = np.loadtxt(detCalib.path + '/' +
                detCalib.name.strip('.txt') +
                '{}.txt'.format( detCalib.fileNumber if
                    np.isfinite(detCalib.fileNumber) else '' ) ) 
    else:
        detCalib.factors = np.ones(16)
    if verbose:
        print 'Detector factors =', detCalib.factors


    return detCalib


def saveDetectorCalibration(master_loop, detCalib, config, verbose=False, beta=0):
    # Concatenate the list of all calibration data
    calibValues = np.concatenate(master_loop.calibValues, axis=0)
    # Check the data that is not NAN
    I = np.isfinite(calibValues.sum(axis=1))
    # average the finite values
    average = calibValues[I,:].mean(axis=0)
    #factors = average[config.boolFitMask].max()/average
    params = cookie_box.initial_params()
    params['A'].value = 1
    params['beta'].value = beta
    params['tilt'].value = 0
    params['linear'].value = 1
    factors = cookie_box.model_function(params, np.radians(np.arange(0, 360,
        22.5))) * float(average.max()) / average
    factors[~config.boolFitMask] = np.nan
    
    if verbose:
        print len(calibValues)
        print master_loop.calibValues[0].shape
        print calibValues.shape
        print average
        print 'Calibration factors:', factors

    calibFile = (detCalib.path + '/' + detCalib.name +
                    '{}.txt'.format( detCalib.fileNumber+1 if
                    np.isfinite(detCalib.fileNumber) else '' ) )

    np.savetxt(calibFile, factors)
          

def master_data_setup(args):
    # Container for the master data
    # This is mainly for the averaging buffers
    master_data = aolUtil.struct()
    #master_data.energyAmplitude = None
    #master_data.energyAmplitudeRoi0 = None
    #master_data.timeAmplitude = None
    #master_data.timeAmplitudeFiltered = None
    #master_data.timeAmplitudeRoi0 = None
    #master_data.timeAmplitudeRoi1 = None

    # Storage for the trace avareaging buffer
    master_data.time_trace_buffer = deque([], args.traceAverage)
    master_data.roi_0_buffer = deque([], args.roi0Average)
    master_data.energy_trace_buffer = deque([], args.traceAverage)

    # Storage for roi averaging
    master_data.roi_0_buffer = deque([], args.roi0Average)
    master_data.roi_1_buffer = deque([], args.roi1Average)

    # Buffers for the polarization averageing
    master_data.pol_degree_buffer = deque([], args.polAverage)
    master_data.pol_angle_buffer = deque([], args.polAverage)
    master_data.pol_roi0_buffer = deque([], args.polAverage)

    return master_data
    
def master_loop_setup(args):
    # Master loop data collector
    master_loop = aolUtil.struct()
    # Define the plot interval from the command line input
    master_loop.tPlot = args.plotInterval

    # Make template of the array that should be sent between the ranks
    master_loop.buf_template = np.empty((dSize,), dtype=float)
    
    # Initialize the stop time
    master_loop.tStop = time.time()
    
    # Calibration
    if args.calibrate > -1:
        master_loop.calibValues = []
    
    return master_loop


def get_scales(env, cb, verbose=False):
    global dSize
    global dTraces
    global d_energy_trace

    # A struct object to hold the scale information
    scales = aolUtil.struct()

    # Assume that all the tofs have the same energy scale and use only the first
    # one.
    scales.energy_eV = cb.get_energy_scales_eV()[0]
    scales.energyRoi0_eV = cb.get_energy_scales_eV(roi=0)[0]
    scales.energyRoi0_slice = slice(
            scales.energy_eV.searchsorted(np.min(scales.energyRoi0_eV)),
            scales.energy_eV.searchsorted(np.max(scales.energyRoi0_eV), side='right'))
    
    # Get all the time scales
    scales.time_us = cb.get_time_scales_us()
    scales.timeRoi0_us = cb.get_time_scales_us(roi=0)
    scales.timeRoi0_slice = [slice(full.searchsorted(np.min(part)),
                                   full.searchsorted(np.max(part),
                                       side='right')) for full, part in
                                   zip(scales.time_us, scales.timeRoi0_us)]
    scales.timeRoi2_us = cb.get_time_scales_us(roi=2)
    scales.timeRoi1_us = cb.get_time_scales_us(roi=1)
    scales.timeRoi1_slice = [slice(full.searchsorted(np.min(part)),
                                   full.searchsorted(np.max(part),
                                       side='right')) for full, part in
                                   zip(scales.time_us, scales.timeRoi1_us)]
    for i, scale_list in enumerate([scales.timeRoi0_us,
                                    scales.timeRoi1_us,
                                    scales.timeRoi2_us]):
        if np.any([len(s)==0 for s in scale_list]):
            print ('ERROR: Roi {} is empty for at leas one of' + \
                    'the detectors.').format(i)
            sys.exit(0)

    # Calculate the background factors
    scales.tRoi0BgFactors = np.array(
            [np.float(len(s))/len(bg) for
             s,bg in zip(scales.timeRoi0_us, scales.timeRoi2_us)])

    # Get some angle vectors
    scales.angles = cb.getAngles('rad')
    scales.anglesFit = np.linspace(0, 2*np.pi, 100)

    # Update the data size descriptions in the globals
    traces_size = 16 * np.max([len(t) for t in scales.time_us])
    dTraces = slice(dSize, dSize + traces_size)
    dSize += traces_size

    energy_trace_size = len(scales.energy_eV)
    d_energy_trace = slice(dSize, dSize+energy_trace_size)
    dSize += energy_trace_size


    return scales

def event_data_container(args, event=None):
    # Set up some data containers
    if event is None:
        event = aolUtil.struct()
    event.sender = []
    event.fiducials = []
    event.times = []
    event.intRoi0 = []
    event.intRoi0Bg = []
    event.intRoi1 = []
    event.pol = []

    event.ebEnergyL3 = []
    event.gasDet = []

    if args.photonEnergy != 'no':
        event.energy = []

    event.deltaK = []
    event.deltaEnc = []

    return event

def get_event_data(config, scales, detCalib, cb, args, epics, verbose=False):
    if verbose:
        print 'Rank {} grabbing data.'.format(rank)

    data = np.zeros(dSize, dtype=float)
    
    # Get the time amplitudes
    time_amplitudes = cb.get_time_amplitudes_filtered()
    # If there is a None in the data, some axqiris trace is missing
    if None in time_amplitudes:
        # Don't try to pull any more data from this event
        print 'Rank {} detected None in time amplitude.'.format(rank)
        return None

    # This construction can handle traces of different length
    data[dTraces] = np.nan
    length = (dTraces.stop - dTraces.start) / 16
    #print 'length =', length
    #print data.shape, dTraces.start, dTraces.stop
    start = dTraces.start
    for t_sig in time_amplitudes:
        if t_sig is None:
            return req
        #print len(t_sig)
        #print start, start+len(t_sig)
        #print data.shape, data[dTraces].shape
        data[start : start + len(t_sig)] = t_sig
        start += length


#def masterEvent_data(evt_data, cb, master_loop, verbose=False):
#    # Grab the y data
#    try:
#        evt_data.energyAmplitude = cb.get_energy_amplitudes()[7]
#        evt_data.energyAmplitudeRoi0 = cb.get_energy_amplitudes(roi=0)[7]
#        #evt_data.energyAmplitude = np.average(cb.get_energy_amplitudes(),
#        #    axis=0)
#        #evt_data.energyAmplitudeRoi0 = np.average(
#        #    cb.get_energy_amplitudes(roi=0), axis=0)
#        evt_data.timeAmplitude = cb.get_time_amplitudes()
#        evt_data.timeAmplitudeFiltered = cb.get_time_amplitudes_filtered()
#        evt_data.timeAmplitudeRoi0 = cb.get_time_amplitudes_filtered(roi=0)
#        evt_data.timeAmplitudeRoi1 = cb.get_time_amplitudes_filtered(roi=1)
#   except TypeError:
#        evt_data.timeAmplitude = None
#        return
#    for trace in evt_data.timeAmplitudeFiltered:
#        if trace is None:
#            return
#    evt_data.time_amplitude_filtered_buffer.append(evt_data.timeAmplitudeFiltered)
#
#    if verbose:
#        print 'Rank', rank, '(master) grabbed one event.'
#    # Update the event counters
#    master_loop.nProcessed += 1


    # Get the intensities
    # Roi 0 base information (photoline)
    data[dIntRoi0] = np.array([np.sum(trace) for trace
                               in cb.get_time_amplitudes_filtered(roi=0)])
    # Use roi 2 for backgorud subtraction, could be commented out to not do
    # subtraction...
    data[dIntRoi0] -= (np.array([np.sum(trace) for trace
                                 in cb.get_time_amplitudes_filtered(roi=2)]) *
                       scales.tRoi0BgFactors)
    # Rescale the data for roi 0
    data[dIntRoi0] *= config.nanFitMask * detCalib.factors

    # Get roi 1 imformation (auger line)
    data[dIntRoi1] = (np.sum(cb.get_time_amplitudes_filtered(roi=1), axis=1) *
                      detCalib.factors * config.nanFitMask)
        
    # Get the initial fit parameters
    #params = cookie_box.initial_params(evt_data.intRoi0[-1])
    params = cookie_box.initial_params()
    params['A'].value, params['linear'].value, params['tilt'].value = \
            cb.proj.solve(data[dIntRoi0], args.beta)
    # Lock the beta parameter
    params['beta'].value = args.beta
    params['beta'].vary = False

    #print params['A'].value, params['linear'].value, params['tilt'].value
    
    # Perform the fit
    #print scales.angles[config.boolFitMask]
    #print evt_data.intRoi0[-1][config.boolFitMask]
    res = lmfit.minimize(
            cookie_box.model_function,
            params,
            args=(scales.angles[config.boolFitMask],
                data[dIntRoi0][config.boolFitMask]),
            method='leastsq')

    #print params['A'].value, params['linear'].value, params['tilt'].value
    
    #lmfit.report_fit(params)

    # Store the values
    data[dPol] = np.array([params['A'].value, params['A'].stderr,
                           params['beta'].value, params['beta'].stderr,
                           params['tilt'].value, params['tilt'].stderr,
                           params['linear'].value, params['linear'].stderr])
    
    # Get the energy traces
    data[d_energy_trace] = np.mean([trace for trace, test in
        zip(cb.get_energy_amplitudes(), config.energy_spectrum_mask)
        if test], axis=0)

    # Get the photon energy center and width
    if args.photonEnergy != 'no':
        data[dEnergy] = cb.getPhotonEnergy(energyShift=args.energyShift)

    # Get lcls parameters
    data[dEL3] = lcls.getEBeamEnergyL3_MeV()
    data[dFEE] = lcls.getPulseEnergy_mJ()
                        
    # timing information
    data[dFiducials] = lcls.getEventFiducial()
    data[dTime] = lcls.getEventTime()

    data[dDeltaK] = epics.value('USEG:UND1:3350:KACT')
    data[dDeltaEnc] = np.array(
        [epics.value('USEG:UND1:3350:{}:ENC'.format(i)) for i in range(1,5)])

    return data
  
def send_data_to_master(data, req, buffer, verbose=False):
    # wait if there is an active send request
    if req != None:
        req.Wait()
    #copy the data to the send buffer
    buffer = data.copy()
    if verbose and 0:
        print 'rank', rank, 'sending data'
    req = world.Isend([buffer, MPI.FLOAT], dest=0, tag=0)

    return req

 
def merge_arrived_data(data, master_loop, args, scales, verbose=False):
    if verbose:
        print 'Merging master and worker data.'
        #print 'len(buf) =', len(master_loop.buf)
        #print 'buf[0] =', master_loop.buf[0]
    # Unpack the data
    data.sender = np.array([d[dRank] for d in master_loop.buf])
    data.fiducials = np.array([d[dFiducials] for d in master_loop.buf])
    data.times = np.array([d[dTime] for d in master_loop.buf])
    data.intRoi0 = np.array([d[dIntRoi0] for d in master_loop.buf])
    data.intRoi1 = np.array([d[dIntRoi1] for d in master_loop.buf])

    # traces
    data.timeSignals_V = []
    for event_data in master_loop.buf:
        data.timeSignals_V.append([d[:len(scale)] for d, scale in
            zip(event_data[dTraces].reshape(16, -1), scales.time_us)])

    data.energy_signal = [d[d_energy_trace] for d in master_loop.buf]

    if args.photonEnergy != 'no':
        data.energy = np.array([d[dEnergy] for d in master_loop.buf])
    if args.calibrate > -1:
        master_loop.calibValues.append(data.intRoi0 if args.calibrate==0 else
                                       data.intRoi1)

    data.ebEnergyL3 = np.array([d[dEL3] for d in master_loop.buf])
    data.gasDet = np.array([d[dFEE] for d in master_loop.buf])

    data.deltaK = np.array([d[dDeltaK] for d in master_loop.buf])
    data.deltaEnc = np.array([d[dDeltaEnc] for d in master_loop.buf])

    # FEE energy thresholding
    if args.feeTh is None:
        fee_accepted = range(len(master_loop.buf))
    else:
        fee_accepted = [i for i, fee in enumerate(data.gasDet)
                        if fee > args.feeTh]

    # Polarizatio information
    data.pol = np.array([d[dPol] for d in master_loop.buf])
    data.pol_roi0_int = []
    # Moving average over the polarization data
    for i in range(len(master_loop.buf)):
        if i not in fee_accepted:
            data.pol[i, 6] = np.nan
            data.pol[i, 4] = np.nan
            data.pol_roi0_int.append(np.nan)
            continue

        if np.isfinite(data.pol[i, 6]):
            data.pol_degree_buffer.append(data.pol[i, 6])
            data.pol[i, 6] = np.average(data.pol_degree_buffer)

        if np.isfinite(data.pol[i, 4]):
            data.pol_angle_buffer.append(data.pol[i, 4])
            data.pol[i, 4] = np.average(data.pol_angle_buffer)

        amp = data.intRoi0[i].sum()
        if np.isfinite(amp):
            data.pol_roi0_buffer.append(amp)
            data.pol_roi0_int.append(np.average(data.pol_roi0_buffer))
        else:
            data.pol_roi0_int.appens(np.nan)



    # Fill the buffers
    for i in fee_accepted:
        if data.gasDet[i].mean() > args.feeTh:
            data.time_trace_buffer.append(data.timeSignals_V[i])
            data.roi_0_buffer.append(data.intRoi0[i])
            data.roi_1_buffer.append(data.intRoi1[i])
            data.energy_trace_buffer.append(data.energy_signal[i])

    # and compute averages
    if len(data.time_trace_buffer) > 0:
        data.traceAverage = np.mean(data.time_trace_buffer, axis=0)
        data.roi_0_average = np.mean(data.roi_0_buffer, axis=0)
        data.energy_trace_average = np.mean(data.energy_trace_buffer, axis=0)
        data.roi_1_average = np.mean(data.roi_1_buffer, axis=0)
    else:
        data.traceAverage = None
        data.roi_0_average = None
        data.energy_trace_aveage = None
        data.roi_1_average = None
            
def sendPVs(data, scales,  pvHandler, args):
    pv_data = {}
    # Polarization data
    # Should contain degree of circular polarization
    pv_data['polarization'] = np.array( [
        data.fiducials,
        data.pol[:,6],
        data.pol[:,7],
        data.pol[:,4],
        data.pol[:,5]] ).T.reshape(-1)
    # Intensity information in the detectors
    pv_data['intensities'] =  np.concatenate(
        [data.fiducials.reshape(-1,1), data.intRoi0],
        axis=1).reshape(-1)
    # Photon energy information
    if args.photonEnergy != 'no':
        pv_data['energy'] = np.concatenate(
                [data.fiducials.reshape(-1,1), data.energy],
                axis=1).reshape(-1)
        pv_data['spectrum'] = np.concatenate(
                [np.array([scales.energyRoi0_eV[0] + args.energyShift,
                    scales.energyRoi0_eV[1] - scales.energyRoi0_eV[0]]),
                    data.energyAmplitudeRoi0])

    pv_data['ebeam'] = np.array( [data.fiducials, data.ebEnergyL3]
            ).T.flatten() 
    
    # Send the data
    pvHandler.assign_data(verbose=False, **pv_data)
    pvHandler.flush_data()

    return

def zmqPlotting(data, scales, zmq):
 
    plot_data = {}
    plot_data['polar'] = {
            'roi0':data.roi_0_average,
            'roi1': data.roi_1_average,
            'A':data.pol[-1][0],
            'beta':data.pol[-1][2],
            'tilt':data.pol[-1][4],
            'linear':data.pol[-1][6]}


    plot_data['strip'] = [data.fiducials,
                          data.pol,
                          data.pol_roi0_int]
    plot_data['traces'] = {}
    if data.traceAverage != None:
        plot_data['traces']['timeScale'] = scales.time_us
        plot_data['traces']['timeRaw'] = data.timeSignals_V[-1]
        plot_data['traces']['timeFiltered'] = data.traceAverage
        plot_data['traces']['timeScaleRoi0'] = scales.timeRoi0_us
        plot_data['traces']['timeRoi0'] = [trace[slice] for trace, slice in
                                           zip(data.traceAverage,
                                               scales.timeRoi0_slice)]
        plot_data['traces']['timeScaleRoi1'] = scales.timeRoi1_us
        plot_data['traces']['timeRoi1'] = [trace[slice] for trace, slice in
                                           zip(data.traceAverage,
                                               scales.timeRoi1_slice)]
    if args.photonEnergy != 'no':
        plot_data['energy'] = np.concatenate(
                [data.fiducials.reshape(-1,1), data.energy],
                axis=1).reshape(-1)
    try:
        plot_data['spectrum'] = {}
        plot_data['spectrum']['energyScale'] = scales.energy_eV
        plot_data['spectrum']['energyScaleRoi0'] = scales.energyRoi0_eV
        plot_data['spectrum']['energyAmplitude'] = data.energy_trace_average
        plot_data['spectrum']['energyAmplitudeRoi0'] = \
            data.energy_trace_average[scales.energyRoi0_slice]
    except:
        plot_data['spectrum'] = {}
        plot_data['spectrum']['energyScale'] = [0]
        plot_data['spectrum']['energyScaleRoi0'] = [0]
        plot_data['spectrum']['energyAmplitude'] = [0]
        plot_data['spectrum']['energyAmplitudeRoi0'] = [0]

    #print plot_data
            
    zmq.sendObject(plot_data)
                

def openSaveFile(format, online=False, config=None):
    fileName = '/reg/neh/home/alindahl/output/amoi0114/'
    if online is True:
        t = time.localtime()
        fileName += 'online{}-{}-{}_{}-{}-{}.{}'.format(t.tm_year, t.tm_mon,
                t.tm_mday, t.tm_hour, t.tm_min, t.tm_sec, format)
    else:
        fileCount = 0
        if config is None:
            fileName += 'outfile{}.' + format
        else:
            fileName += 'run' + config.offline_source.split('=')[-1] + '_{}.' + format
        while os.path.exists(fileName.format(fileCount)):
            fileCount += 1
        fileName = fileName.format(fileCount)
    if format == 'txt':
        file = open(fileName,'w')

        file.write('eventTime\tfiducial')
        file.write('\tebEnergyL3')
        file.write('\tfee11\tfee12\tfee21\tfee22')
        for i in range(16):
            file.write('\tauger_{}'.format(i))
        for i in range(16):
            file.write('\tphoto_{}'.format(i))
        file.write('\tI0\tI0_err\tbeta\tbeta_err\ttilt\ttilt_err\tlinDegree\tlinDegree_err')
        file.write('\tdeltaK')
        for i in range(4):
            file.write('\tdeltaEnc{}'.format(i+1))

        file.write('\n')
        file.flush()
        return file

def write_dataToFile(file, data, format):
    if format == 'txt':
        for i in range(len(data.sender)):
            line = ( repr( data.times[i] ) + '\t' +
                    repr( data.fiducials[i] ) )
            line += '\t' + repr(data.ebEnergyL3[i])
            for a in data.gasDet[i,:]:
                line += '\t' + repr(a)
            for a in data.intRoi1[i,:]:
                line += '\t' + repr(a)
            for a in data.intRoi0[i,:]:
                line += '\t' + repr(a)
            for a in data.pol[i,:]:
                line += '\t' + repr(a)
            line += '\t' + repr(data.deltaK[i])
            for a in data.deltaEnc[i,:]:
                line += '\t' + repr(a)

            line += '\n'

            file.write(line)

        file.flush()

def closeSaveFile(file):
    try:
        file.close()
    except:
        pass

def main(args, verbose=False):
    try:
        # Import the configuration file
        config = importConfiguration(args, verbose=verbose)
    
        # Make a cookie box object
        cb = cookie_box.CookieBox(config, verbose=False if rank==1 else False)
        cb.proj.setFitMask(config.boolFitMask)
        if args.bgAverage != 1:
            cb.set_baseline_subtraction_averaging(args.bgAverage)


        # Read the detector transmission calibrations
        detCalib = getDetectorCalibration(verbose=verbose,
                fileName=args.gainCalib)

        # Change the configuration fit masks according to the factors
        config.nanFitMask = config.nanFitMask.astype(float)
        config.nanFitMask[np.isnan(detCalib.factors)] = np.nan
        config.boolFitMask[np.isnan(detCalib.factors)] = False

        # Connect to the correct datasource
        ds = connect_to_data_source(args, config, verbose=False)
        if not args.offline:
            # For online use just get the events iterator
            events = ds.events()
        else:
            # For offline use the indexing capabilities are used to enable
            # event skipping and real multi core advantage
            run = ds.runs().next()
            times = run.times()
            event_counter = args.skip + rank

        # Get the epics store
        epics = ds.env().epicsStore()

        # Get the next event. The purpouse here is only to make sure the
        # datasource is initialized enough so that the env object is avaliable.
        if not args.offline:
            evt = events.next()
        else:
            evt = run.event(times[event_counter])
            event_counter += 1
           
        # Get the scales that we need
        cb.setup_scales(config.energy_scale_eV, ds.env())
        scales = get_scales(ds.env(), cb)

        # The master have some extra things to do
        if rank == 0:
            # Set up the plotting in AMO
            if args.sendPlots:
                from ZmqSender import zmqSender
                zmq = zmqSender()
    
            # Set up the PV handler
            if args.sendPV:
                import pv_handler
                pvHandler = pv_handler.PvHandler(timeout=1.0)
    
            master_loop = master_loop_setup(args) 

            master_data = master_data_setup(args)

            if args.save_data != 'no':
                saveFile = openSaveFile(args.save_data, not args.offline, config)
            
        else:
            # set an empty request for the mpi send to master
            req = None
            # and the corresponding buffer
            buffer = np.empty(dSize, dtype=float)
    
        event_data = None
    
        # The main loop that never ends...
        while 1:
            # An event data container
            event_data = event_data_container(args, event=event_data)
                
            # The master should receive data in a timed loop
            if rank == 0:
                # Empty the buffer list
                master_loop.buf = []
                # and enter the timed loop
                if verbose:
                    print 'Master enters the timed loop.',
                    print 't_stop - time.time() = {} s.'.format(master_loop.tStop -
                                                          time.time())
                while time.time() < master_loop.tStop:
                    # Append a buffer
                    #if verbose:
                    #    print 'Master appending a new buffer.'
                    master_loop.buf.append(master_loop.buf_template.copy())
                    # Make a blockign receive from anyone
                    #if verbose:
                    #    print 'Master waits for data.'
                    world.Recv([master_loop.buf[-1], MPI.FLOAT],
                               source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG)
                    #if verbose:
                    #    print 'Master got data.'

                # On loop exit increment the stop time
                master_loop.tStop += master_loop.tPlot

                # Check how many events arrived
                master_loop.nArrived = len(master_loop.buf)
                if verbose:
                    print 'Master received {} events'.format(
                            master_loop.nArrived)
                if master_loop.nArrived == 0:
                    continue

                merge_arrived_data(master_data, master_loop, args,
                        scales, verbose=False)
    
                # Send data as PVs
                if args.sendPV:
                    sendPVs(master_data, scales, pvHandler, args)
    
                # Send data for plotting
                if args.sendPlots:
                    zmqPlotting(master_data, scales, zmq)

                if args.save_data != 'no':
                    write_dataToFile(saveFile, master_data, args.save_data)

                if ((args.calibrate > -1) and
                    (len(master_loop.calibValues))):
                        saveDetectorCalibration(master_loop, detCalib, config,
                                                verbose=False,
                                                beta=args.calibBeta)
                   
            else:
                # The workers
                # Get the next event
                if not args.offline:
                    evt = events.next()
                else:
                    evt = run.event(times[event_counter])
                    event_counter += workers

                # Pass the event to the processing
                cb.set_raw_data(evt, newDataFactor=args.floatingAverage)
                lcls.setEvent(evt)
    
                # Randomize the amplitudes if this is requested 
                if args.randomize:
                    cb.randomize_amplitudes()
    
                # Get the data for the event
                data = get_event_data(config, scales, detCalib,
                        cb, args, epics, verbose=False)
                if data is None:
                    if verbose:
                        print 'Rank {} got empty event.'.format(rank)
                    continue
                
                # Send the data to master
                req = send_data_to_master(data, req, buffer,
                                          verbose=verbose)


    except KeyboardInterrupt:
        print "Terminating program."

        if rank == 0:
            if args.save_data != 'no':
                closeSaveFile(saveFile)



if __name__ == '__main__':
    # Start here
    # parset the command line
    args = arguments.parse()
    if args.verbose:
        print args
        if args.sendPV:
            print 'Will send PVs'
        else:
            print 'Will NOT send PVs'

    main(args, args.verbose)

