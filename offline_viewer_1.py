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
fee = {}
auger_sum = {}
auger_signals_average = {}
photo_signals = {}
photo_signals_average = {}
photo_signals_average_corrected = {}
photo_bg = {}
bg_factors = {}
photo_signals_corrected = {}
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
    photo_signals[run] = []
    bg_factors[run] = []

    # with h5py.File(h5_name, 'r+') as h5_file:
    h5_file = h5py.File(h5_name, 'r+')

    valid = np.zeros(h5_file['fee'].len(), dtype=bool)
    hits = h5_file.attrs.get('n_events_set')
    valid[:hits] = 1

    fee[run] = h5_file['fee'][:, 2:].mean(axis=1)
    valid *= np.isfinite(fee[run]) * (fee[run] > 0.01)
    fee_mean[run] = np.mean(fee[run][valid], axis=0)

    for i in range(16):
        traces[run].append(h5_file['time_amplitudes/det_{}'.format(i)].value)

        average_traces[run].append(np.average(traces[run][i][valid, :],
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

        bg_factors[run].append((photo_roi_slices[run][i].stop -
                                photo_roi_slices[run][i].start) /
                               (photo_bg_slices[run][i].stop -
                                photo_bg_slices[run][i].start))

        auger_roi_slices[run].append(
            slice(time_scales[run][i].searchsorted(auger_roi[i][0]),
                  time_scales[run][i].searchsorted(auger_roi[i][1],
                                                   side='right')))

        photo_signals[run].append(
            traces[run][i][photo_roi_slices[run][i]].mean() -
            traces[run][i][photo_bg_slices[run][i]].mean() *
            bg_factors[run][i])

    auger_signals_average[run] = np.array([np.sum(tr[sl]) for tr, sl in
                                           zip(average_traces[run],
                                           auger_roi_slices[run])])
    auger_sum[run] = np.sum(auger_signals_average[run])

    photo_signals_average[run] = np.array([np.sum(tr[sl]) for tr, sl in
                                          zip(average_traces[run],
                                              photo_roi_slices[run])])

    photo_bg[run] = np.array([np.sum(tr[sl]) for tr, sl in
                             zip(average_traces[run],
                                 photo_bg_slices[run])])

    photo_signals_average_corrected[run] = photo_signals_average[run] - \
        photo_bg[run]

# Renormalization
norm_factor_285 = auger_sum[283] / auger_sum[285]
for i in range(16):
    average_traces[285][i] *= norm_factor_285

# %%
try:
    trace_plot = plt.figure('Trace plot')
    trace_plot.clf()
except:
    pass

trace_plot, trace_axis_array = plt.subplots(4, 4, sharex=True, sharey=True,
                                            num='Trace plot')

for i_run, run in enumerate([285, 283]):
    for i, ax in enumerate(trace_axis_array.flatten()):
        ax.plot(time_scales[run][i], average_traces[run][i],
                '-{}'.format('by'[i_run]),
                label='{} {} deg'.format(run, 22.5*i))
        ax.plot(time_scales[run][i][photo_roi_slices[run][i]],
                average_traces[run][i][photo_roi_slices[run][i]], '.r')
        ax.plot(time_scales[run][i][auger_roi_slices[run][i]],
                average_traces[run][i][auger_roi_slices[run][i]], '.g')

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

        ax.set_title('Run {}'.format(run))

ax.set_xlim(220, 250)
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
               args=(phi, photo_signals_average_corrected[285]))
beta2_factors = (cookie_box.model_function(norm_params, phi) /
                 photo_signals_average_corrected[285])

I_fit = np.ones(16, dtype=bool)
I_fit[4] = 0
I_fit[8] = 0
proj.setFitMask(I_fit)

full_falctors = {}
auger_factors = {}

for ax, run in zip(angular_axis_array, [283, 285]):
    auger_factors[run] = (auger_signals_average[run] /
                          auger_signals_average[285])
    auger_factors[run] /= auger_factors[run].max()
    full_falctors[run] = beta2_factors * auger_factors[run]
    full_falctors[run] /= full_falctors[run].max()
#    factors = np.max(auger_signals[run]) / auger_signals[run]
    ax.plot(phi, auger_signals_average[run], 'gx')
    ax.plot(phi, auger_signals_average[run] * auger_factors[run], 'gs')
    ax.plot(phi, photo_signals_average_corrected[run], 'rx',
            label='raw')
    ax.plot(phi, photo_signals_average_corrected[run] * full_falctors[run],
            'ro', label='normalized')

    params = cookie_box.initial_params()
    params['beta'].vary = False
    params['A'].value, params['linear'].value, params['tilt'].value = \
        proj.solve(photo_signals_average_corrected[run] * full_falctors[run],
                   2)

    res = lmfit.minimize(cookie_box.model_function, params,
                         args=(phi[I_fit],
                               (photo_signals_average_corrected[run] *
                                full_falctors[run])[I_fit]))

    lmfit.report_fit(params)
    ax.plot(phi_line, cookie_box.model_function(params, phi_line), 'm',
            label='fit')

    ax.set_title('Run {}'.format(run))
    ax.legend(loc='best', fontsize='small')
plt.tight_layout()
#    ax.plot(phi, photo_signals_average_corrected[run] * factors, 'ro')

# %%
try:
    pol_prog_fig = plt.figure('Polar progression')
    pol_prog_fig.clf()
except:
    pass

n_rows = int(np.floor(np.sqrt(float(len(runs)))))
n_cols = int(np.ceil(float(len(runs))/n_rows))

pol_prog_fig, pol_prog_axis_array = plt.subplots(n_rows, n_cols,
                                                 subplot_kw={'polar': True},
                                                 num='Polar progression')

for ax, run in zip(pol_prog_axis_array.flatten(), runs):
    ax.set_title('Run {}'.format(run))
    ax.plot(phi, photo_signals_average_corrected[run], 'rx')
    ax.plot(phi, photo_signals_average_corrected[run] * beta2_factors, 'ro')

    params = cookie_box.initial_params()
    params['beta'].vary = False
    params['A'].value, params['linear'].value, params['tilt'].value = \
        proj.solve(photo_signals_average_corrected[run] * beta2_factors, 2)
    res = lmfit.minimize(cookie_box.model_function, params,
                         args=(phi[I_fit],
                               (photo_signals_average_corrected[run] *
                                beta2_factors)[I_fit]))

    ax.plot(phi_line, cookie_box.model_function(params, phi_line), 'm')

    print 'Run {}'.format(run)
    lmfit.report_fit(params)

plt.tight_layout()
