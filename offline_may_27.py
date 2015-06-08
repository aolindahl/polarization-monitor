# -*- coding: utf-8 -*-
"""
Created on Wed May 27 18:49:48 2015

@author: antlin
"""

import offline_handler
import numpy as np
import matplotlib.pyplot as plt
from aolPyModules import cookie_box
import lmfit

plt.ion()

phi = np.linspace(0, 2 * np.pi, 16, endpoint=False)
phi_line = np.linspace(0, 2 * np.pi, 1000)

auger_roi_720 = [226, 234]
photo_roi_720 = [[236.5, 250],
                 [236.5, 250],
                 [242.0, 260],
                 [242.0, 260],
                 [236.5, 250],
                 [236.5, 250],
                 [236.5, 250],
                 [236.5, 250],
                 [236.5, 250],
                 [236.5, 250],
                 [236.5, 250],
                 [236.5, 250],
                 [236.5, 250],
                 [236.5, 250],
                 [236.5, 250],
                 [236.5, 250]]

auger_roi_1200 = [217, 220]
photo_roi_1200 = [
                    [231, 240],
                    [231, 240],
                    [234, 244],
                    [234, 244],
                    [231, 240],
                    [231, 240],
                    [231, 240],
                    [231, 240],
                    [231, 240],
                    [231, 240],
                    [231, 240],
                    [231, 240],
                    [231, 240],
                    [231, 240],
                    [231, 240],
                    [231, 240],
                    [231, 240],
                    ]

fit_mask = np.ones(16, dtype=bool)
fit_mask[2:4] = 0

ref_720 = offline_handler.DataSet('data/amom0115_35_1.h5', 'ref 720')
ref_1200 = offline_handler.DataSet('data/amom0115_41_0.h5', 'ref 1200')

data_list = []
data_list.append(offline_handler.DataSet('data/amom0115_35_1.h5',
                                         'run 35 720 eV lin (ref)'))
data_list.append(offline_handler.DataSet('data/amom0115_36_1.h5',
                                         'run 36 720 eV c1'))
data_list.append(offline_handler.DataSet('data/amom0115_38_0.h5',
                                         'run 38 720 eV c2'))
data_list.append(offline_handler.DataSet('data/amom0115_41_0.h5',
                                         'run 41 1200 eV lin (ref)'))
data_list.append(offline_handler.DataSet('data/amom0115_46_0.h5',
                                         'run 46 1200 eV c1'))
data_list.append(offline_handler.DataSet('data/amom0115_52_0.h5',
                                         'run 52 1200 eV c2'))


ref_720.auger_slices = auger_roi_720
ref_720.photo_slices = photo_roi_720
ref_1200.auger_slices = auger_roi_1200
ref_1200.photo_slices = photo_roi_1200
# offline_handler.plot_traces('ref_traces', ref)

for ref in [ref_720, ref_1200]:
    photo_ref_raw = ref.photo_amplitudes.mean(axis=0)
    auger_ref_raw = ref.auger_amplitudes.mean(0)
    
    auger_factors = auger_ref_raw.max() / auger_ref_raw
    
    auger_ref = auger_ref_raw * auger_factors
    photo_ref = photo_ref_raw * auger_factors
    
    params = cookie_box.initial_params(photo_ref)
    params['linear'].value = 1
    params['linear'].vary = False
    params['beta'].value = 2
    params['beta'].vary = False
    params['tilt'].value = 0
    params['tilt'].vary = False
    
    lmfit.minimize(cookie_box.model_function, params, args=(phi, photo_ref))
    
    ref_factors = cookie_box.model_function(params, phi) / photo_ref_raw
    ref_factors[[3, 4, 5, 11, 12, 13]] = auger_factors[[3, 4, 5, 11, 12, 13]]
    
    photo_ref_corr = photo_ref_raw * ref_factors
    
    plt.figure('{} polar'.format(ref.name))
    plt.clf()
    plt.subplot(111, polar=False)
    plt.plot(phi, photo_ref_raw, 'rx')
    plt.plot(phi, photo_ref, 'ro')
    plt.plot(phi, photo_ref_corr, 'r.')
    plt.plot(phi_line, cookie_box.model_function(params, phi_line), '-')
    
    if '720' in ref.name:
        ref_factors_720 = ref_factors.copy()
    else:
        ref_factors_1200 = ref_factors.copy()


pol_fig = plt.figure('720 eV polar')
pol_fig.clf()
for i, data in enumerate(data_list):
    data.auger_slices = auger_roi_720 if i < 3 else auger_roi_1200
    data.photo_slices = photo_roi_720 if i < 3 else photo_roi_1200
    data.det_factors = ref_factors_720 if i < 3 else ref_factors_1200
    tr_fig = offline_handler.plot_traces('traces {}'.format(data.name), data)
    ax = pol_fig.add_subplot(2, 3, i+1, polar=True)
    ax.set_title(format(data.name))
    pol_fig = offline_handler.polar_plot(data, ax=ax,
                                         reset_scaling=False,
                                         fit_mask=fit_mask)

