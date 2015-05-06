# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 15:32:02 2015

@author: Anton O Lindahl
"""
import h5py
import numpy as np
import matplotlib.pyplot as plt
import lmfit

import sys
import os.path
sys.path.append(os.path.dirname(os.path.abspath(__file__)) +
                '/aolPyModules')
import cookie_box
from BurningDetectors_V6 import projector

proj = projector()

plt.ion()

h5_names = ['data/amom0115_5_1.h5',
            'data/amom0115_16_0.h5',
            'data/amom0115_25_0.h5']
h5_names = ['data/amom0115_31_0.h5']


photo_roi = [[240, 250]]*16
photo_roi_1 = [[237, 241]]*16
photo_roi_2 = [[242.5, 250]]*16

auger_roi = [[215, 220]]*16

traces = {}
average_traces = {}
time_scales = {}
photo_roi_slices = {}
photo_roi_1_slices = {}
photo_roi_2_slices = {}
photo_bg_slices = {}
photo_bg_1_slices = {}
photo_bg_2_slices = {}
auger_roi_slices = {}
fee_mean = {}
fee_valid = {}
fee = {}
auger_sum = {}
auger_signals = {}
auger_signals_average = {}
auger_factors = {}
photo_signals_average_corrected = {}
photo_signals_1_average_corrected = {}
photo_signals_2_average_corrected = {}
bg_factors = {}
bg_1_factors = {}
bg_2_factors = {}
photo_signals = {}
photo_signals_1 = {}
photo_signals_2 = {}
photo_bg = {}
photo_bg_1 = {}
photo_bg_2 = {}
photo_signals_corrected = {}
photo_signals_1_corrected = {}
photo_signals_2_corrected = {}
bg_fit_coeffs = {}
runs = []

for h5_name in h5_names:
    run = int(h5_name.split('_')[1])
    runs.append(run)

    traces[run] = []
    average_traces[run] = []
    time_scales[run] = []
    photo_roi_slices[run] = []
    photo_roi_1_slices[run] = []
    photo_roi_2_slices[run] = []
    auger_roi_slices[run] = []
    photo_bg_slices[run] = []
    photo_bg_1_slices[run] = []
    photo_bg_2_slices[run] = []
    photo_bg[run] = []
    photo_bg_1[run] = []
    photo_bg_2[run] = []
    photo_signals[run] = []
    photo_signals_1[run] = []
    photo_signals_2[run] = []
    photo_signals_corrected[run] = []
    photo_signals_1_corrected[run] = []
    photo_signals_2_corrected[run] = []
    auger_signals[run] = []
    bg_factors[run] = []
    bg_1_factors[run] = []
    bg_2_factors[run] = []
    bg_fit_coeffs[run] = []

    # with h5py.File(h5_name, 'r+') as h5_file:
    h5_file = h5py.File(h5_name, 'r+')

    valid = np.zeros(h5_file['fee'].len(), dtype=bool)
    hits = h5_file.attrs.get('n_events_set')
    valid[:hits] = 1

    fee[run] = h5_file['fee'][:, 2:].mean(axis=1)
    valid *= np.isfinite(fee[run]) * (fee[run] > 0.005)
    fee_valid[run] = fee[run][valid]
    fee_mean[run] = fee_valid[run].mean()

    for i in range(16):
        traces[run].append(
            h5_file['time_amplitudes/det_{}'.format(i)].value[valid, :])

        average_traces[run].append(np.average(traces[run][i],
                                              axis=0) / np.mean(fee_mean[run]))
#        average_traces[run].append(np.average(
#            h5_file['time_amplitudes/det_{}'.format(i)].value[valid, :],
#            axis=0) * 1e3)

        time_scales[run].append(
            h5_file['time_scales/det_{}'.format(i)].value * 1e3)

        photo_roi_slices[run].append(
            slice(time_scales[run][i].searchsorted(photo_roi[i][0]),
                  time_scales[run][i].searchsorted(photo_roi[i][1],
                                                   side='right')))

        photo_roi_1_slices[run].append(
            slice(time_scales[run][i].searchsorted(photo_roi_1[i][0]),
                  time_scales[run][i].searchsorted(photo_roi_1[i][1],
                                                   side='right')))

        photo_roi_2_slices[run].append(
            slice(time_scales[run][i].searchsorted(photo_roi_2[i][0]),
                  time_scales[run][i].searchsorted(photo_roi_2[i][1],
                                                   side='right')))

        photo_bg_slices[run].append(slice(photo_roi_slices[run][i].start - 30,
                                          photo_roi_slices[run][i].start))

        photo_bg_1_slices[run].append(
            slice(photo_roi_1_slices[run][i].start - 10,
                  photo_roi_1_slices[run][i].start))

        photo_bg_2_slices[run].append(
            slice(photo_roi_2_slices[run][i].start - 5,
                  photo_roi_2_slices[run][i].start))

        bg_fit_coeffs[run].append(np.polyfit(
            time_scales[run][i][photo_bg_slices[run][i]],
            average_traces[run][i][photo_bg_slices[run][i]], 1))

        bg_factors[run].append(
            np.polyval(bg_fit_coeffs[run][i],
                       time_scales[run][i][photo_roi_slices[run][i]]).sum() /
            np.polyval(bg_fit_coeffs[run][i],
                       time_scales[run][i][photo_bg_slices[run][i]]).sum())

#        bg_factors[run].append((photo_roi_slices[run][i].stop -
#                                photo_roi_slices[run][i].start) /
#                               (photo_bg_slices[run][i].stop -
#                                photo_bg_slices[run][i].start))

        bg_1_factors[run].append((photo_roi_1_slices[run][i].stop -
                                  photo_roi_1_slices[run][i].start) /
                                 (photo_bg_1_slices[run][i].stop -
                                  photo_bg_1_slices[run][i].start))

        bg_2_factors[run].append((photo_roi_2_slices[run][i].stop -
                                  photo_roi_2_slices[run][i].start) /
                                 (photo_bg_2_slices[run][i].stop -
                                  photo_bg_2_slices[run][i].start))

        auger_roi_slices[run].append(
            slice(time_scales[run][i].searchsorted(auger_roi[i][0]),
                  time_scales[run][i].searchsorted(auger_roi[i][1],
                                                   side='right')))

        photo_signals[run].append(
            traces[run][i][:, photo_roi_slices[run][i]].sum(axis=1))
        photo_bg[run].append(
            traces[run][i][:, photo_bg_slices[run][i]].sum(axis=1) *
            bg_factors[run][i])
        photo_signals_corrected[run].append(
            photo_signals[run][i] - photo_bg[run][i])

        photo_signals_1[run].append(
            traces[run][i][:, photo_roi_1_slices[run][i]].sum(axis=1))
        photo_bg_1[run].append(
            traces[run][i][:, photo_bg_1_slices[run][i]].sum(axis=1) *
            bg_1_factors[run][i])
        photo_signals_1_corrected[run].append(
            photo_signals_1[run][i] - photo_bg_1[run][i])

        photo_signals_2[run].append(
            traces[run][i][:, photo_roi_2_slices[run][i]].sum(axis=1))
        photo_bg_2[run].append(
            traces[run][i][:, photo_bg_2_slices[run][i]].sum(axis=1) *
            bg_2_factors[run][i])
        photo_signals_2_corrected[run].append(
            photo_signals_2[run][i] - photo_bg_2[run][i])

        auger_signals[run].append(
            traces[run][i][:, auger_roi_slices[run][i]].sum(axis=1))

    photo_signals_average_corrected[run] = np.average(
        photo_signals_corrected[run], axis=1)
    photo_signals_1_average_corrected[run] = np.average(
        photo_signals_1_corrected[run], axis=1)
    photo_signals_2_average_corrected[run] = np.average(
        photo_signals_2_corrected[run], axis=1)

    auger_signals_average[run] = np.average(
        auger_signals[run], axis=1)
    auger_factors[run] = (auger_signals_average[run].max() /
                          auger_signals_average[run])
    auger_sum[run] = np.sum(auger_signals_average[run])

# %% Signal tests
#
#sig_min = -0.5
#sig_max = 4
#n_sig_bins = 2**7
#sig_ax = np.linspace(sig_min, sig_max, 2 * n_sig_bins + 1)[1::2]
#
#n_rows = int(np.floor(np.sqrt(float(len(runs)))))
#n_cols = int(np.ceil(float(len(runs))/n_rows))
#
#sig_level_fig = plt.figure('Signal levels')
#sig_level_fig.clf()
#sig_level_fig, sig_level_axis_array = plt.subplots(n_rows*2, n_cols,
#                                                   num='Signal levels')
#det = 4
#for i_run, run in enumerate(runs):
#    ax = sig_level_axis_array.flatten()[i_run]
#    ax.set_title('run {}'.format(run))
#    ax.hist(photo_signals[run][det], n_sig_bins, (sig_min, sig_max))
#    ax.hist(photo_bg[run][det], n_sig_bins, (sig_min, sig_max))
#    ax.hist(photo_signals_corrected[run][det], n_sig_bins, (sig_min, sig_max))
#
#    ax = sig_level_axis_array.flatten()[i_run + len(runs)]
#    ax.plot(fee_valid[run], photo_signals_corrected[run][det], '.')
#    ax.plot(fee_valid[run], auger_signals[run][det], '.')
#sig_level_fig.tight_layout()
# %% Trace plots
try:
    trace_plot = plt.figure('Trace plot')
    trace_plot.clf()
except:
    pass

trace_plot, trace_axis_array = plt.subplots(4, 4, sharex=True, sharey=True,
                                            num='Trace plot')

for i_run, run in enumerate(runs):
    for i, ax in enumerate(trace_axis_array.flatten()):
        ax.plot(time_scales[run][i], average_traces[run][i],
                '-{}'.format('byc'[i_run]),
                label='{} {} deg'.format(run, 22.5*i))

        ax.plot(time_scales[run][i][auger_roi_slices[run][i]],
                average_traces[run][i][auger_roi_slices[run][i]], '.g')

        if run != 31:
            ax.plot(time_scales[run][i][photo_roi_slices[run][i]],
                    average_traces[run][i][photo_roi_slices[run][i]], '.r')
            ax.plot(time_scales[run][i][photo_bg_slices[run][i]],
                    average_traces[run][i][photo_bg_slices[run][i]], '.m')
    #        ax.plot(time_scales[run][i],
    #                np.polyval(bg_fit_coeffs[run][i], time_scales[run][i]),
    #                'm')
        else:
            ax.plot(time_scales[run][i][photo_roi_1_slices[run][i]],
                    average_traces[run][i][photo_roi_1_slices[run][i]], '.r')
            ax.plot(time_scales[run][i][photo_bg_1_slices[run][i]],
                    average_traces[run][i][photo_bg_1_slices[run][i]], '.m')

            ax.plot(time_scales[run][i][photo_roi_2_slices[run][i]],
                    average_traces[run][i][photo_roi_2_slices[run][i]], '.y')
            ax.plot(time_scales[run][i][photo_bg_2_slices[run][i]],
                    average_traces[run][i][photo_bg_2_slices[run][i]], '.c')
        if i % 4:
            plt.setp(ax.get_yticklabels(), visible=False)
        else:
            ax.set_ylabel('signal (mV)')
        if i / 4 < 3:
            plt.setp(ax.get_xticklabels(), visible=False)
        else:
            ax.set_xlabel('time (ns)')

        ax.grid(True)
        ax.legend(loc='best', fontsize='x-small', ncol=1)

#        ax.set_title('Run {}'.format(run))

ax.set_xlim(200, 260)
plt.tight_layout()

# %%
try:
    angular_plot = plt.figure('Angular')
    angular_plot.clf()
except:
    pass

angular_plot, angular_axis_array = plt.subplots(1, 3,
                                                subplot_kw={'polar': True},
                                                num='Angular',
                                                figsize=(8, 10))
phi = cookie_box.phi_rad
phi_line = np.linspace(0, 2*np.pi, 2**8)
norm_params = cookie_box.initial_params()
norm_params['A'].value = 1
norm_params['beta'].value = 2
norm_params['beta'].vary = False
norm_params['tilt'].value = 0
norm_params['tilt'].vary = False
norm_params['linear'].value = 1
norm_params['linear'].vary = False

#lmfit.minimize(cookie_box.model_function, norm_params,
#               args=(phi,
#                     photo_signals_average_corrected[5] * auger_factors[5]))
#beta2_factors = (cookie_box.model_function(norm_params, phi) /
#                 (photo_signals_average_corrected[5] * auger_factors[5]))

I_fit = np.ones(16, dtype=bool)
#I_fit[[4, 5, 11, 12]] = 0
proj.setFitMask(I_fit)

full_falctors = {}

for ax, run in zip(angular_axis_array, runs):
#    signal_factors = beta2_factors
#    signal_factors = auger_factors[run] * beta2_factors
    signal_factors = auger_factors[run]
    ax.plot(phi, auger_signals_average[run], 'gx', label='auger raw')
    ax.plot(phi, auger_signals_average[run] * signal_factors, 'gs',
            label='auger scaled')
    ax.plot(phi, photo_signals_average_corrected[run], 'rx',
            label='photo raw')
    ax.plot(phi, photo_signals_average_corrected[run] * signal_factors,
            'ro', label='photo scaled')

    params = cookie_box.initial_params()
    params['beta'].vary = False
    params['A'].value, params['linear'].value, params['tilt'].value = \
        proj.solve(photo_signals_average_corrected[run] * signal_factors,
                   2)

    res = lmfit.minimize(cookie_box.model_function, params,
                         args=(phi[I_fit],
                               (photo_signals_average_corrected[run] *
                                signal_factors)[I_fit]))

    lmfit.report_fit(params)
    ax.plot(phi_line, cookie_box.model_function(params, phi_line), 'm',
            label='fit')

    ax.set_title('Run {}'.format(run))
    ax.legend(loc='upper right', fontsize='x-small')
angular_plot.tight_layout()
#    ax.plot(phi, photo_signals_average_corrected[run] * factors, 'ro')

# %% Run 31 polar
if 31 in traces.keys():
    run = 31
    signal_factors = auger_factors[run]
    fig31 = plt.figure(run)
    fig31.clf()
    ax = plt.subplot(131, polar=True)
    ax.plot(phi, auger_signals_average[run], 'gx', label='auger raw')
    ax.plot(phi, auger_signals_average[run] * signal_factors, 'gs',
            label='auger scaled')
    ax.legend(loc='best', fontsize='small')
    ax.set_title('run31 augers')

    params = cookie_box.initial_params()
    params['beta'].vary = False
    params['A'].value, params['linear'].value, params['tilt'].value = \
        proj.solve(photo_signals_1_average_corrected[run] * signal_factors,
                   2)
    lmfit.minimize(cookie_box.model_function, params,
                   args=(phi[I_fit],
                         (photo_signals_1_average_corrected[run] *
                         signal_factors)[I_fit]))
    ax = plt.subplot(132, polar=True)
    ax.plot(phi, photo_signals_1_average_corrected[run], 'rx',
            label='photo 1 raw')
    ax.plot(phi, photo_signals_1_average_corrected[run] * signal_factors,
            'ro', label='photo 1 scaled')
    ax.plot(phi_line, cookie_box.model_function(params, phi_line), 'm',
            label='fit {:.3} lin {:.3} circ'.format(
                params['linear'].value,
                np.sqrt(1 - params['linear'].value**2)))
    ax.legend(loc='best', fontsize='small')
    ax.set_title('run31 first photoline (high energy)')

    params['A'].value, params['linear'].value, params['tilt'].value = \
        proj.solve(photo_signals_2_average_corrected[run] * signal_factors,
                   2)
    lmfit.minimize(cookie_box.model_function, params,
                   args=(phi[I_fit],
                         (photo_signals_2_average_corrected[run] *
                         signal_factors)[I_fit]))
    ax = plt.subplot(133, polar=True)
    ax.plot(phi, photo_signals_2_average_corrected[run], 'yx',
            label='photo 2 raw')
    ax.plot(phi, photo_signals_2_average_corrected[run] * signal_factors,
            'yo', label='photo 2 scaled')
    ax.plot(phi_line, cookie_box.model_function(params, phi_line), 'm',
            label='fit {:.3} lin {:.3} circ'.format(
                params['linear'].value,
                np.sqrt(1 - params['linear'].value**2)))
    ax.legend(loc='best', fontsize='small')
    ax.set_title('run31 second photoline (low energy)')
    fig31.tight_layout()

# %%
#try:
#    pol_prog_fig = plt.figure('Polar progression')
#    pol_prog_fig.clf()
#except:
#    pass
#
#n_rows = int(np.floor(np.sqrt(float(len(runs)))))
#n_cols = int(np.ceil(float(len(runs))/n_rows))
#
#pol_prog_fig, pol_prog_axis_array = plt.subplots(n_rows, n_cols,
#                                                 subplot_kw={'polar': True},
#                                                 num='Polar progression')
#
#for ax, run in zip(pol_prog_axis_array.flatten(), runs):
#    ax.set_title('Run {}'.format(run))
#    ax.plot(phi, photo_signals_average_corrected[run], 'rx')
##    ax.plot(phi, photo_signals_average_corrected[run] * beta2_factors, 'ro')
#
#    params = cookie_box.initial_params()
#    params['beta'].vary = False
##    params['A'].value, params['linear'].value, params['tilt'].value = \
##        proj.solve(photo_signals_average_corrected[run] * beta2_factors, 2)
##    res = lmfit.minimize(cookie_box.model_function, params,
##                         args=(phi[I_fit],
##                               (photo_signals_average_corrected[run] *
##                                beta2_factors)[I_fit]))
#
#    ax.plot(phi_line, cookie_box.model_function(params, phi_line), 'm')
#
#    print 'Run {}'.format(run)
#    lmfit.report_fit(params)
#
#plt.tight_layout()
