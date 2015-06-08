import argparse
import numpy as np

def parse():
    parser = argparse.ArgumentParser(
            description =
            'Argument parser for the cookie box DELTA comissioning support')

    run_controll = parser.add_argument_group('Execution controll.')
    averaging = parser.add_argument_group('Averaging',
        'Parameters for averaging of the data.' +
        '\nFloating average should preferably be used on its own.' +
        ' Otherwhise double averaging is performed.')
    calibration = parser.add_argument_group('Gain calibration')
    offline = parser.add_argument_group('Offline')
    debug = parser.add_argument_group('Debug information')

    offline.add_argument(
            '-o',
            '--offline',
            action='store_true',
            help=('Use this flag to run from dummy data.' +
                'Permissions to acess the data might cause problems.' +
                'Default is to run from shared memory.')
            )

    #group1 = parser.add_mutually_exclusive_group()

    #group1.add_argument(
    #        '--alignement',
    #        action='store_true',
    #        help='Display to use when aligning the cookie box.')

    #group1.add_argument(
    #parser.add_argument(
    #        '--calibrate',
    #        default = -1,
    #        type = int,
    #        metavar='ROI',
    #        help=('All data in the gven ROI will be used to make' +
    #              ' a gain calibration.'))

    #parser.add_argument(
    #        '--calibBeta',
    #        default = 0,
    #        type = float,
    #        metavar = 'BETA',
    #        help = 'Beta parameter to be used for the calibration')

    debug.add_argument(
            '-v',
            '--verbose',
            action='store_true',
            help = 'Print stuff describing what is going on.')

    debug.add_argument(
            '-r',
            '--randomize',
            action='store_true',
            help = ('Randimize the amplitudes. Might me usefull' +
                'when running dummy data.'))

    run_controll.add_argument(
            '-i',
            '--plotInterval',
            metavar = 'TIME',
            type = float,
            default = 0.5,
            help = 'Time in seconds to wait between plots. Default = 0.5 s')

    run_controll.add_argument(
            '-b',
            '--beta',
            type = float,
            default = 2,
            help = ('Set the beta parameter used in the fitting.' +
                'Default = 2'))

    run_controll.add_argument(
            '-c',
            '--configuration',
            type = str,
            metavar = 'FILENAME',
            default = '/reg/neh/home/alindahl/cookieBoxDefaultConfig.py',
            help = ('Name of the module containing the configuration' +
                'information for the analysis. Default = ' +
                '"cookieBoxDefaultConfig.py"'))

    
    calibration.add_argument(
            '--gainCalib',
            metavar = 'FILENAME',
            type = str,
            default = '',
            help = ('Path to the gain/transmission calibration file. ' \
                    + 'By default the file with the highest number ' \
                    + 'will be used.'))

    run_controll.add_argument(
            '-P',
            '--sendPV',
            action = 'store_true',
            help = ('Use this flag to send data as PVs over the EPICS' +
                ' system. Default = NO PV writing.'))
    

    averaging.add_argument(
            '-f',
            '--floatingAverage',
            type = float,
            metavar = 'WEIGHT',
            default = None,
            help = ('Moving average with weight on last point.' +
                    ' This averaging is performed each core when' +
                    ' the raw data is grabbed from the evet.'))

    averaging.add_argument(
            '--feeTh',
            type=float,
            metavar = 'mJ',
            default = None,
            help='FEE threshold for averaging. Default: no threshold.')

    averaging.add_argument(
            '-tA',
            '--traceAverage',
            default = 1,
            type = int,
            metavar = 'numShots',
            help=('Averaging for time and energy traces in the plots.' +
                  ' Does not effect any other calculations.' +
                  ' Number of shots to use. Default = 1.'))

    averaging.add_argument(
            '-0A',
            '--roi0Average',
            type = int,
            default = 1,
            metavar = 'numShots',
            help = ('Averaging for the intensities in roi 0.' +
                ' Number of shots to use.' +
                ' Default = 1.'))

    averaging.add_argument(
            '-1A',
            '--roi1Average',
            type = int,
            default = 1,
            metavar = 'numShots',
            help = ('Averaging for the intensities in roi 1.' +
                ' Number of shots to use.' +
                ' Default = 1.'))

    averaging.add_argument(
            '-pA',
            '--polAverage',
            type=int,
            metavar='numShots',
            default=1,
            help=('Averaging for the polarization parameters.' +
                  ' Default = 1.'))

    averaging.add_argument(
            '--bgAverage',
            type = float,
            default = 1,
            metavar = 'WEIGHT',
            help = ('Averaging the regions used for baseline and background'
                + ' subtraction. Default: shot to shot calculation.'))


#    energyShifts = {
#            'He1s':24.6,
#            'Ne1s':870.2,
#            'Ar2s':326.3
#            }

    #parser.add_argument(
    #        '-pE',
    #        '--photonEnergy',
    #        default = 'no',
    #        #metavar = 'STATE',
    #        type = str,
    #        choices = energyShifts.keys(),
    #        help = ('Calculate photon energy information.' +
    #            ' Passed value is the offset from the time of flight' +
    #            ' energy scale. Default = no.'))

    calibration.add_argument(
            '--calibrate',
            default = -1,
            type = int,
            metavar='ROI',
            help=('All data in the gven ROI will be used to make' +
                  ' a gain calibration.'))

    calibration.add_argument(
            '--calibBeta',
            default = 0,
            type = float,
            metavar = 'BETA',
            help = 'Beta parameter to be used for the calibration.' +
                    ' Default = 0') 

    run_controll.add_argument(
            '-s',
            '--save_data',
            type = str,
            default = 'no',
            choices = ['no', 'txt', 'h5'],
            help = ('Output the data to a file' + 
                ' Default = no'))

    offline.add_argument(
            '--skip',
            type = int,
            default = 0,
            help = ('Skip events at start of run, for offline use only.'))

    offline.add_argument(
            '--num-events',
            type = int,
            default = -1,
            help = ('Number of events to process in offline analysis ' +
                    'before exiting. NUM_EVENTS < 0 => all. Default = -1'))

    offline.add_argument(
            '-d',
            '--dataSource',
            type=str,
            default=None,
            help=('Data source string that overrides the one set in the' +
                  ' configuration file. Only for offline use.'))

    run_controll.add_argument(
            '--no-plot',
            dest='sendPlots',
            action='store_false',
            help=('Stop the script from sending plots.'))

    args = parser.parse_args()

    # Unused optinos that sort of live in the code but should not be used
    args.photonEnergy = 'no'
    #args.floatingAverage = None

    if args.photonEnergy != 'no':
        args.energyShift = energyShifts[args.photonEnergy]

    return args




if __name__ == '__main__':
    args =  parse()
    print args
