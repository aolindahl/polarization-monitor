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

runs = [283, 285]

h5_names = ['data/amoi0114_283_2.h5',
            'data/amoi0114_285_2.h5',
            'data/amoi0114_275_0.h5',
            'data/amoi0114_276_0.h5',
            'data/amoi0114_277_0.h5',
            'data/amoi0114_278_0.h5',
            'data/amoi0114_279_0.h5',
            'data/amoi0114_280_0.h5',
            'data/amoi0114_281_0.h5',
            'data/amoi0114_282_0.h5']

h5_names = [name for name in h5_names if
            np.any([str(run) in name for run in runs])]

photo_roi = [[232.5, 239], [233.5, 241.5], [237.5, 245.5], [237.5, 245.5],
             [233.5, 241.5], [233.5, 241.5], [233, 241], [233, 241],
             [233, 241], [233, 241], [233, 241], [233, 241],
             [233, 241], [233, 242], [233.5, 241], [233, 241]]

auger_roi = [[224, 230.5], [225, 231.5], [226.5, 233], [226.5, 233],
             [225, 231.5], [225, 231.5], [225, 231.5], [224.5, 231],
             [225, 231.5], [225, 231.5], [224.5, 231], [224.5, 231],
             [225, 231.5], [225, 231.5], [225, 231.5], [224.5, 231]]

traces = {}
average_traces = {}
time_scales = {}
photo_roi_slices = {}
photo_bg_slices = {}
auger_roi_slices = {}
fee_mean = {}
fee_valid = {}
fee = {}
auger_sum = {}
auger_signals = {}
auger_signals_average = {}
auger_factors = {}
photo_signals_average_corrected = {}
bg_factors = {}
photo_signals = {}
photo_bg = {}
photo_signals_corrected = {}
bg_fit_coeffs = {}
runs = []

for h5_name in h5_names:
    run = int(h5_name.split('_')[1])
    runs.append(run)

    traces[run] = []
    average_traces[run] = []
    time_scales[run] = []
    photo_roi_slices[run] = []
    auger_roi_slices[run] = []
    photo_bg_slices[run] = []
    photo_bg[run] = []
    photo_signals[run] = []
    photo_signals_corrected[run] = []
    auger_signals[run] = []
    bg_factors[run] = []
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

        photo_bg_slices[run].append(slice(photo_roi_slices[run][i].start - 3,
                                          photo_roi_slices[run][i].start))

        bg_fit_coeffs[run].append(np.polyfit(
            time_scales[run][i][photo_bg_slices[run][i]],
            average_traces[run][i][photo_bg_slices[run][i]], 1))

#        bg_factors[run].append(
#            np.polyval(bg_fit_coeffs[run][i],
#                       time_scales[run][i][photo_roi_slices[run][i]]).sum() /
#            np.polyval(bg_fit_coeffs[run][i],
#                       time_scales[run][i][photo_bg_slices[run][i]]).sum())

        bg_factors[run].append((photo_roi_slices[run][i].stop -
                                photo_roi_slices[run][i].start) /
                               (photo_bg_slices[run][i].stop -
                                photo_bg_slices[run][i].start))

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

        auger_signals[run].append(
            traces[run][i][:, auger_roi_slices[run][i]].sum(axis=1))

    photo_signals_average_corrected[run] = np.average(
        photo_signals_corrected[run], axis=1)

    auger_signals_average[run] = np.average(
        auger_signals[run], axis=1)
    auger_factors[run] = (auger_signals_average[run].max() /
                          auger_signals_average[run])
    auger_sum[run] = np.sum(auger_signals_average[run])

# %% Signal tests

sig_min = -0.5
sig_max = 4
n_sig_bins = 2**7
sig_ax = np.linspace(sig_min, sig_max, 2 * n_sig_bins + 1)[1::2]

n_rows = int(np.floor(np.sqrt(float(len(runs)))))
n_cols = int(np.ceil(float(len(runs))/n_rows))

sig_level_fig = plt.figure('Signal levels')
sig_level_fig.clf()
sig_level_fig, sig_level_axis_array = plt.subplots(n_rows*2, n_cols,
                                                   num='Signal levels')
det = 4
for i_run, run in enumerate(runs):
    ax = sig_level_axis_array.flatten()[i_run]
    ax.set_title('run {}'.format(run))
    ax.hist(photo_signals[run][det], n_sig_bins, (sig_min, sig_max))
    ax.hist(photo_bg[run][det], n_sig_bins, (sig_min, sig_max))
    ax.hist(photo_signals_corrected[run][det], n_sig_bins, (sig_min, sig_max))

    ax = sig_level_axis_array.flatten()[i_run + len(runs)]
    ax.plot(fee_valid[run], photo_signals_corrected[run][det], '.')
    ax.plot(fee_valid[run], auger_signals[run][det], '.')
plt.tight_layout()
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
        ax.plot(time_scales[run][i][photo_roi_slices[run][i]],
                average_traces[run][i][photo_roi_slices[run][i]], '.r')
        ax.plot(time_scales[run][i][auger_roi_slices[run][i]],
                average_traces[run][i][auger_roi_slices[run][i]], '.g')
        ax.plot(time_scales[run][i][photo_bg_slices[run][i]],
                average_traces[run][i][photo_bg_slices[run][i]], '.m')
        ax.plot(time_scales[run][i],
                np.polyval(bg_fit_coeffs[run][i], time_scales[run][i]),
                'm')

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

ax.set_xlim(220, 245)
plt.tight_layout()

# %%
try:
    angular_plot = plt.figure('Angular')
    angular_plot.clf()
except:
    pass

angular_plot, angular_axis_array = plt.subplots(1, 2,
                                                subplot_kw={'polar': True},
                                                num='Angular')
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

lmfit.minimize(cookie_box.model_function, norm_params,
               args=(phi,
                     photo_signals_average_corrected[285] *
                     auger_factors[285]))
beta2_factors = (cookie_box.model_function(norm_params, phi) /
                 (photo_signals_average_corrected[285] * auger_factors[285]))

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
    ax.legend(loc='best', fontsize='small')
plt.tight_layout()
#    ax.plot(phi, photo_signals_average_corrected[run] * factors, 'ro')

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
