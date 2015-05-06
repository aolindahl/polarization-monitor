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

h5_names = ['data/amom0115_5_0.h5',
            'data/amom0115_16_0.h5',
            'data/amom0115_17_0.h5',
            'data/amom0115_22_0.h5',
            'data/amom0115_23_0.h5',
            'data/amom0115_24_0.h5',
            'data/amom0115_25_0.h5']
h5_names = ['data/amom0115_31_0.h5']


photo_roi = [[240, 250]]*16

auger_roi = [[215, 220]]*16

phi = cookie_box.phi_rad

for h5_name in h5_names:
    run = int(h5_name.split('_')[1])

    traces = []
    average_traces = []
    time_scales = []
    photo_roi_slices = []
    auger_roi_slices = []
    photo_bg_slices = []
    photo_bg = []
    photo_signals = []
    photo_signals_corrected = []
    auger_signals = []
    bg_factors = []
    bg_fit_coeffs = []

    # with h5py.File(h5_name, 'r+') as h5_file:
    h5_file = h5py.File(h5_name, 'r+')

    valid = np.zeros(h5_file['fee'].len(), dtype=bool)
    hits = h5_file.attrs.get('n_events_set')
    valid[:hits] = 1

    fee = h5_file['fee'][:, 2:].mean(axis=1)
    valid *= np.isfinite(fee) * (fee > 0.001)
    time = h5_file['event_time'][valid]
    valid = np.where(valid)[0][time.argsort()]
    time.sort()
    time -= time[0]
    fee_valid = fee[valid]
    fee_mean = fee_valid[run].mean()

    n_valid = len(valid)

    fiducial = h5_file['fiducial'].value[valid]
    enc = h5_file['delta_encoders'].value[valid, :]

    for det in range(16):
        traces.append(
            h5_file['time_amplitudes/det_{}'.format(det)].value[valid, :])

        average_traces.append(
            np.average(traces[det], axis=0) / np.mean(fee_mean))

        time_scales.append(
            h5_file['time_scales/det_{}'.format(det)].value * 1e3)

        photo_roi_slices.append(
            slice(time_scales[det].searchsorted(photo_roi[det][0]),
                  time_scales[det].searchsorted(photo_roi[det][1],
                                                side='right')))

        photo_bg_slices.append(slice(photo_roi_slices[det].start - 30,
                                     photo_roi_slices[det].start))

        bg_fit_coeffs.append(np.polyfit(
            time_scales[det][photo_bg_slices[det]],
            average_traces[det][photo_bg_slices[det]], 1))

        bg_factors.append(
            np.polyval(bg_fit_coeffs[det],
                       time_scales[det][photo_roi_slices[det]]).sum() /
            np.polyval(bg_fit_coeffs[det],
                       time_scales[det][photo_bg_slices[det]]).sum())

#        bg_factors.append((photo_roi_slices[det].stop -
#                                photo_roi_slices[det].start) /
#                               (photo_bg_slices[det].stop -
#                                photo_bg_slices[det].start))

        auger_roi_slices.append(
            slice(time_scales[det].searchsorted(auger_roi[det][0]),
                  time_scales[det].searchsorted(auger_roi[det][1],
                                                side='right')))

        photo_signals.append(
            traces[det][:, photo_roi_slices[det]].sum(axis=1))
        photo_bg.append(
            traces[det][:, photo_bg_slices[det]].sum(axis=1) * bg_factors[det])

        photo_signals_corrected.append(photo_signals[det] - photo_bg[det])

        auger_signals.append(traces[det][:, auger_roi_slices[det]].sum(axis=1))

    photo_signals_corrected = np.array(photo_signals_corrected)
    photo_signals_average_corrected = np.average(
        photo_signals_corrected, axis=1)

    auger_signals_average = np.average(auger_signals, axis=1)
    auger_factors = (auger_signals_average.max() / auger_signals_average)
    auger_sum = np.sum(auger_signals_average)

    params = cookie_box.initial_params()
    params['beta'].value = 2
    params['beta'].vary = False

    n_points = 1000
    n_in_average = 2**3
    stepsize = np.max((1, n_valid / n_points))
    n_points = np.min((n_points, n_valid))
    idx_list = np.arange(n_points) * stepsize

    tilt = np.empty_like(idx_list, dtype=float)
    tilt_error = np.empty_like(idx_list, dtype=float)
    lin = np.empty_like(idx_list, dtype=float)
    lin_error = np.empty_like(idx_list, dtype=float)

    fee_select = np.empty_like(idx_list, dtype=float)

    for point, idx in enumerate(idx_list):
        s = slice(idx - int((n_in_average-1.0)/2),
                  idx + 1 + int((n_in_average-1.0)/2))

        data = photo_signals_corrected[:, s].mean(axis=1) * auger_factors

        fee_select[point] = fee_valid[s].mean()

        params['A'].value, params['linear'].value, params['tilt'].value = \
            proj.solve(data, 2)

        res = lmfit.minimize(cookie_box.model_function, params,
                             args=(phi, data))
        tilt[point] = params['tilt'].value
        tilt_error[point] = params['tilt'].stderr
        lin[point] = params['linear'].value
        lin_error[point] = params['linear'].stderr

    tilt = np.rad2deg(tilt)
    tilt_error = np.rad2deg(tilt_error)
    time_select = time[idx_list]

    pol_plot = plt.figure('Polarization run {}'.format(run), figsize=(8, 10))
    pol_plot.clf()
    pol_plot.suptitle('Run {}, {} shot average'.format(run, n_in_average))

    lin_degree_ax = pol_plot.add_subplot(3, 2, 1)
    lin_degree_ax.errorbar(time_select, lin*100, lin_error*100, fmt='.',
                           capsize=0, errorevery=n_points/100)
    lin_degree_ax.set_ylabel('degree linear (%)')
    lin_degree_ax.set_title('Linear polarization')
    ybounds = lin_degree_ax.get_ybound()
    lin_degree_ax.set_ylim(np.max((ybounds[0], -10)),
                           np.min((ybounds[1], 110)))

    circ_degree_ax = pol_plot.add_subplot(3, 2, 2, sharex=lin_degree_ax)
    circ_degree_ax.errorbar(time_select, np.sqrt(1-lin**2)*100,
                            np.sqrt(lin**2/(1.0-lin**2)) * lin_error*100,
                            fmt='.', capsize=0, errorevery=n_points/100)
    circ_degree_ax.set_ylabel('degree circular (%)')
    circ_degree_ax.set_title('Circular polarization')
    ybounds = circ_degree_ax.get_ybound()
    circ_degree_ax.set_ylim(np.max((ybounds[0], -10)),
                            np.min((ybounds[1], 110)))

    tilt_ax = pol_plot.add_subplot(3, 2, 3, sharex=lin_degree_ax)
    tilt_ax.errorbar(time_select, tilt, tilt_error, fmt='.',
                     capsize=0, errorevery=n_points/100)
    tilt_ax.set_ylabel('tilt angle (degrees)')
    tilt_ax.set_title('Polarization tilt')
    ybounds = tilt_ax.get_ybound()
    tilt_ax.set_ylim(np.max((ybounds[0], -95)),
                     np.min((ybounds[1], 95)))

    enc_ax = pol_plot.add_subplot(3, 2, 5, sharex=lin_degree_ax)
    for i_enc in range(4):
        enc_ax.plot(time[idx_list], enc[idx_list, i_enc],
                    label='encoder {}'.format(i_enc+1))
    enc_ax.set_xlabel('time (s)')
    enc_ax.set_ylabel('encoder values (?)')
    enc_ax.set_title('Delta position encoders')
    enc_ax.legend(loc='best', fontsize='small', ncol=2)

    fee_ax = pol_plot.add_subplot(3, 2, 4, sharex=lin_degree_ax)
    fee_ax.plot(time_select, fee_select, '.')
    fee_ax.set_ylabel('energy (mJ)')
    fee_ax.set_title('Downstream gas detectors')

    pol_plot.tight_layout(rect=(0, 0, 1, 0.95))
    pol_plot.savefig('figs/' + pol_plot.get_label().replace(' ', '_') + '.pdf')
