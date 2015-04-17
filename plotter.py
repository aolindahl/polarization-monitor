import sys
import os.path
sys.path.append( os.path.dirname(os.path.abspath(__file__)) +  '/aolPyModules')
from collections import deque
import ZmqSender
import matplotlib.animation as animation
import zmq
import matplotlib.pyplot as plt
import numpy as np
import time
import cookie_box
import argparse


def mainPlotter(args, verbose=False):
    """\
    Plotting function for zmq transmitted data\
    """

    fulladdress='tcp://'+args.host+':'+args.port
    if verbose:
        print 'Connecting to: ', fulladdress

    context = zmq.Context()
    socket = context.socket(zmq.SUB)
    socket.setsockopt(zmq.SUBSCRIBE, 'data')
    socket.set_hwm(2)
    socket.connect(fulladdress)

    raw = 0
    filt = 1
    roi0 = 2
    roi1 = 3

    angle = np.linspace(0, np.pi*2, 1000)
  
    nHistory = 1e4
    storage = {
            'fiducial':deque([], nHistory),
            'degree':deque([], nHistory),
            'err_degree':deque([], nHistory),
            'tilt':deque([], nHistory),
            'err_tilt':deque([], nHistory),
            'intFull':deque([], nHistory),
            'intRoi0':deque([], nHistory),
            'intRoi1':deque([], nHistory) }


    def updatePlot(figs):
        socket.recv()
        data = socket.recv_pyobj(zmq.NOBLOCK)
        # we've got data, see if it is the most recent
        while 1:
            try:
                socket.recv(zmq.NOBLOCK)
                data = socket.recv_pyobj()
                # if we get here with no exception there was
                # more data
                print '*** throw away data'
                continue
            except zmq.error.Again:
                print '*** got last data'
                break
        if verbose: 
            print 'Data contains:'
            for k, v in data.iteritems():
                print '\t{} of type {}'.format(k, type(v)),
                if isinstance(v, dict):
                    print 'with components:'
                    for k2, v2 in v.iteritems():
                        print '\t\t{} of type {}'.format(k2, type(v2)),

                        if isinstance(v2, list):
                            print 'of length {}.'.format(len(v2))
                            print '\t\t\tlengths = (',
                            for elem in v2:
                                try:
                                    print '{},'.format(len(elem)),
                                except:
                                    pass
                            print ')'

                        elif isinstance(v2, np.ndarray):
                            print 'of shape {}.'.format(v2.shape)
                else:
                    print ''


            #try:
            #    print 'polar:roi0 =', data['polar']['roi0']
            #    print 'polar:roi1 =', data['polar']['roi1']
            #except:
            #    pass

            #try:
            #    print 'err_degree = ', data['strip'][1][:,7]
            #except:
            #    pass

        plotKey = 'polar'
        if plotKey in data.keys():
            d = data[plotKey]
            ax = figs[1].axes[0]
            lines = ax.lines
            for l, roi in enumerate(['roi0', 'roi1']):
                if roi in d.keys():
                    try:
                        lines[l].set_ydata(d[roi]/d[roi][ np.isfinite(d[roi])
                            ].max())
                    except:
                        print 'Data error'

                    
            # The fit in figure 2
            params = cookie_box.initial_params()
            for k in params:
                params[k].value = d[k]
            params['A'].value /= d['roi0'][ np.isfinite(d['roi0']) ].max()
            
            lines[2].set_ydata(cookie_box.model_function(params, angle))
            
            if args.beta != None:
                params['beta'].value = args.beta
                params['linear'].value = 1.0
                params['tilt'].value = 0
                y = cookie_box.model_function(params, angle)
                lines[3].set_ydata(y/y.max())
            
            ax.relim()
            ax.autoscale_view()
            figs[1].canvas.draw()



        plotKey = 'strip'
        if plotKey in data.keys():
            ax = figs[2].axes[0]
            ax.cla()
            ax.set_title('Degree of linear polarization.')

            storage['fiducial'].extend( data[plotKey][0] )
            storage['degree'].extend( data[plotKey][1][:,6] )
            storage['err_degree'].extend( data[plotKey][1][:,7] )
            storage['tilt'].extend( data[plotKey][1][:,4] )
            storage['err_tilt'].extend( data[plotKey][1][:,5] )
            #storage['intFull'].extend(data[plotKey][2])
            storage['intRoi0'].extend(data[plotKey][2])
            #storage['intRoi1'].extend(data[plotKey][4])



            #ax.plot(storage['fiducial'], storage['degree'],'.')
            #plt.sca(ax)
            if args.errors:
                ax.errorbar(np.array(storage['fiducial']),
                            np.array(storage['degree']),
                            np.array(storage['err_degree']),
                            linestyle='', marker='.', capsize=0)
            else:
                ax.plot(np.array(storage['fiducial']), 
                        np.array(storage['degree']),
                        linestyle='', marker='.')

            if args.a:
                ax.relim()
                ax.autoscale_view()
            else:
                ax.set_ylim(-0.4, 1.4)
            #ax.plot(data[plotKey][0], data[plotKey][1][:,6], '.b')
            #ax.errorbar(data[plotKey][0], data[plotKey][1][:,6],
            #        data[plotKey][1][:,7], ls='None', color='b')
            #ax.relim()
            #ax.autoscale_view()
            #ax.ylim


            ax = figs[2].axes[1]
            ax.cla()
            ax.set_title('Polarization tilt')
            if args.errors:
                ax.errorbar(np.array(storage['fiducial']),
                            np.array(storage['tilt']),
                            np.array(storage['err_tilt']),
                            linestyle='', marker='.', capsize=0)
            else:
                ax.plot(np.array(storage['fiducial']),
                        np.array(storage['tilt'])*180/np.pi,
                        linestyle='', marker='.')

            if args.a:
                ax.relim()
                ax.autoscale_view()
            else:
                ax.set_ylim(-90, 90)
            #ax.errorbar(data[plotKey][0], data[plotKey][1][:,4],
            #        data[plotKey][1][:,5], ls='None', color='b')

            ax = figs[2].axes[2]
            ax.cla()
            ax.set_title('Photoline intensity')
            ax.plot(np.array(storage['fiducial']),
                    np.array(storage['intRoi0']),
                    '.')
            ax.relim()
            ax.autoscale_view()
            min_value, _ = ax.get_ylim()
            ax.set_ylim(bottom=np.minimum(min_value, 0))


            figs[2].canvas.draw()


        plotKey = 'traces'
        if plotKey in data.keys():
            plot = data[plotKey]
            for i, ax in enumerate(figs[0].axes):
                lines = ax.lines
                
                if 'timeScale' in plot.keys():
                    lines[raw].set_xdata(data[plotKey]['timeScale'][i])
                    lines[filt].set_xdata(data[plotKey]['timeScale'][i])
                if 'timeScaleRoi0' in plot.keys():
                    lines[roi0].set_xdata(data[plotKey]['timeScaleRoi0'][i])
                if 'timeScaleRoi1' in plot.keys():
                    lines[roi1].set_xdata(data[plotKey]['timeScaleRoi1'][i])

                if 'timeRaw' in plot.keys():
                    lines[raw].set_ydata(plot['timeRaw'][i])
                if 'timeFiltered' in plot.keys():
                    lines[filt].set_ydata(plot['timeFiltered'][i])
                if 'timeRoi0' in plot.keys():
                    lines[roi0].set_ydata(plot['timeRoi0'][i])
                if 'timeRoi1' in plot.keys():
                    lines[roi1].set_ydata(plot['timeRoi1'][i])

                ax.relim()
                ax.autoscale_view()

            figs[0].canvas.draw()


        plotKey = 'spectrum'
        if plotKey in data.keys():
            plotData = data[plotKey]
            ax = figs[3].axes[0]
            
            l = ax.lines[0]
            l.set_xdata(plotData['energyScale'])
            l.set_ydata(plotData['energyAmplitude'])

            l = ax.lines[1]
            l.set_xdata(plotData['energyScaleRoi0'])
            l.set_ydata(plotData['energyAmplitudeRoi0'])

            ax.relim()
            ax.autoscale_view()

            ax.figure.canvas.draw()
            


    def initializePlot():
        if verbose:
            print 'init'

        # tof traces
        if verbose:
            print 'Make the tof trace plot.'
        fig1, f1ax = plt.subplots(4,4, sharex=True, sharey=True, num='Traces',
                figsize = (24, 18))

        #f1ax = np.array( fig1.axes).reshape(4,4)
        for ax in f1ax[:,0]:
            ax.set_ylabel('signal [V]')
        for ax in f1ax[-1,:]:
            ax.set_xlabel('time [us]')
        for i, ax in enumerate(f1ax.flatten()):
           ax.plot([],[], '--m', label='{} single'.format(22.5*i))
           ax.plot([],[], 'b', label='average')
           ax.plot([],[], 'r', label='ROI 0')
           ax.plot([],[], 'g', label='ROI 1')
           ax.grid(True)
           ax.legend(loc='upper right', prop={'size': 8})


        # Polar plot
        if verbose:
            print 'Make the polar plot.'
        fig2 = plt.figure('ROI amplitudes')
        f2ax = fig2.add_subplot(111, polar=True)
        f2ax.plot(np.arange(0, 2*np.pi, np.pi/8), np.ones(16), 'or',
                label='Photoelectrons')
        f2ax.plot(np.arange(0, 2*np.pi, np.pi/8), np.ones(16), 'sg',
                label='Auger electrons')

        f2ax.plot(angle, np.ones_like(angle), '-b')
        if args.beta != None:
            f2ax.plot(angle, np.ones_like(angle), '--k')

        f2ax.set_ylim(0, 1.1)

        f2ax.legend()

        # Polarization strip shart
        if verbose:
            print 'Make the polarization strip chart.'
        fig3 = plt.figure('Polarization')
        fig3.add_subplot(311)
        fig3.axes[0].set_title('Degree of linear polarization')
        fig3.axes[0].set_ylim(-0.5, 1.5)
        fig3.axes[0].grid(True)
        fig3.add_subplot(312)
        fig3.axes[1].set_title('Angle of residual linear polarization')
        fig3.axes[1].set_ylim(-2, 2)
        fig3.axes[1].grid(True)
        fig3.add_subplot(313)
        #fig3.axes[1].set_title('Angle of residual linear polarization')
        #fig3.axes[1].set_ylim(-2, 2)
        fig3.axes[1].grid(True)

        #f3ax.plot([],[], [], '.', label='Degree of circular polarization')

        # Energy spectrum plot
        if verbose:
            print 'Make the energy spectrum plot.'
        fig4 = plt.figure('Energy spectrum')
        fig4.add_subplot(111)
        fig4.axes[0].plot([],[],'b', [],[],'r')

        # Show the plots
        if args.show_plots:
            if verbose:
                print 'Showing the plots.'
            fig1.show()
            fig2.show()
            fig3.show()
            fig4.show()

        if verbose:
            print 'init finished.'

        return [fig1, fig2, fig3, fig4]
        
    figs = initializePlot()
    
    try:
        while 1:
            updatePlot(figs)

    except:
        if verbose:
            print 'Disconnect from host'
        socket.close()
        raise

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
            description = 'Plotter for cookie box data sent over the network')

    parser.add_argument(
            '-H', '--host', type = str, default = 'localhost',
            help = 'Host to subscribe to. Default = localhost')

    parser.add_argument(
            '-p', '--port', type=str, default='19820809',
            help='Port number to connect to. Default = 19820809')

    parser.add_argument(
            '-b', '--beta', type=float, default=None,
            help= 'Beta parameter to be used in the referece plot.' + 
            'Default gives no reference plot.')

    parser.add_argument(
            '-v', '--verbose', action='store_true',
            help = 'Print stuff ot the prompt.')

    parser.add_argument(
            '-a', action='store_true',
            help='Autoscale strip chart y axis.')

    parser.add_argument(
            '-e', '--errors', action='store_true',
            help='Plot with error bars in the strip chart.')

    parser.add_argument(
            '--no-plots', action='store_false', dest='show_plots',
            help='For debuggin, does not show the plots.')

    args = parser.parse_args()
    
    mainPlotter(args, args.verbose)
