# -*- coding: utf-8 -*-
"""
Created on Wed May 27 12:59:34 2015

@author: antlin
"""

import h5py
import numpy as np
import sys
import matplotlib.pyplot as plt
from aolPyModules import cookie_box
import lmfit

photo_roi = [[236.5, 250],
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

auger_roi = [226, 234]


def list_hdf5_content(group, indent='  '):
    for k, v in group.iteritems():
        print '{}"{}"'.format(indent, k),
        if isinstance(v, h5py.Group):
            print 'group with members:'
            list_hdf5_content(v, indent=indent + '  ')
        elif isinstance(v, h5py.Dataset):
            print '\t{} {}'.format(v.shape, v.dtype)


def get_average(obj, name=None, selection=slice(None)):
    if name is None:
        dset = obj
    else:
        try:
            dset = obj[name]
        except KeyError as e:
            print e.message
            sys.exit(1)
        return dset[selection].mean(axis=0)


def plot_traces(fig_num, dataset):
    fig = plt.figure(fig_num)
    fig.clf()
    ax = None
    for i in range(16):
        ax = fig.add_subplot(4, 4, i+1, sharex=ax, sharey=ax)
        t_axis = dataset.time_scales[i]
        ax.plot(t_axis, dataset.time_amplitudes_averaged[i, :],
                label='{} deg'.format(i*22.5))

        photo_sl = dataset.photo_slices[i]

        auger_sl = dataset.auger_slices[i]

        ax.plot(t_axis[photo_sl],
                dataset.time_amplitudes_averaged[i, photo_sl],
                'r',
                label='photo')
        ax.plot(t_axis[auger_sl],
                dataset.time_amplitudes_averaged[i, auger_sl],
                'g',
                label='auger')
        ax.legend(loc='best', fontsize='small')
        ax.grid(True)

    ax.set_xbound(upper=280)
    plt.tight_layout()

    return fig


def polar_plot(dataset, fig_name=None, polar=True, ax=None,
               reset_scaling=False,
               fit_mask=np.ones(16, dtype=bool)):
    if ax is None:
        fig = plt.figure(fig_name)
        fig.clf()
        ax = fig.add_subplot(111, polar=polar)
    else:
        fig = ax.figure
    phi = np.linspace(0, 2 * np.pi, 16, endpoint=False)
    phi_line = np.linspace(0, 2 * np.pi, 2**10)

    if reset_scaling:
        det_factors = dataset.det_factors.copy()
        dataset.det_factors = np.ones_like(det_factors)

    auger = dataset.auger_amplitudes.mean(axis=0)
    photo = dataset.photo_amplitudes.mean(axis=0)

    if reset_scaling:
        det_calib = auger.max() / auger
        ax.plot(phi, auger, 'gx')
        auger *= det_calib

        ax.plot(phi, photo, 'rx')
        photo *= det_calib

    ax.plot(phi, auger, 'gs', label='auger')
    ax.plot(phi, photo, 'ro', label='photo')

    params = cookie_box.initial_params(photo)
    params['beta'].vary = False
    params['beta'].value = 2
    lmfit.minimize(cookie_box.model_function, params,
                   args=(phi[fit_mask], photo[fit_mask]))
    lmfit.report_fit(params)

    ax.plot(phi_line, cookie_box.model_function(params, phi_line),
            '-m', label='{:.1f} % lin {:.1f} deg'.format(
            params['linear'].value*100,
            np.rad2deg(params['tilt'].value)))

    ax.grid(True)
    ax.legend(loc='center', bbox_to_anchor=(0, 0), fontsize='medium')
    plt.tight_layout()

    if reset_scaling:
        dataset.det_factors = det_factors

    return fig


def get_bg_index_list(sl, nr_points):
    return (range(sl.start - nr_points, sl.start) +
            range(sl.stop, sl.stop + nr_points))


class DataSet(object):
    """Class to handle the a dataset contnained in an hdf5 file."""

    @property
    def h5_file(self):
        return self._h5_file

    @h5_file.setter
    def h5_file(self, x):
        print 'WARNING: The "h5_file" property can only be set at creation.'

    def __init__(self, file_name, name=None):
        """Initialize the instance based on a file name."""

        self._h5_name = file_name
        self._h5_file = h5py.File(file_name, 'r')

        self._name = file_name if name is None else name

        self._det_factors = np.ones(16, dtype=float)

        self._time_scales = np.array(())
        self._time_amplitudes_averaged = np.array(())
        self._time_amplitudes_averaged_selection = slice(None)

        self._photo_amplitudes = np.array(())
        self._auger_amplitudes = np.array(())

        self.event_selection = slice(None)

    def __del__(self):
        self._h5_file.close()

    @property
    def name(self):
        return self._name

    def print_content(self, indent=''):
        list_hdf5_content(self.h5_file, indent=indent)

    def get_average(self, path, selection=slice(None)):
        return get_average(self.h5_file, path, selection=selection)

    @property
    def time_amplitudes_averaged(self):
        if (len(self._time_amplitudes_averaged) == 0 or
            self._time_amplitudes_averaged_selection != self.event_selection):

            group = self.h5_file['time_amplitudes']
            n_detectors = len(group.keys())
            self._time_amplitudes_averaged = np.array([
                get_average(group, 'det_{}'.format(i),
                            selection=self.event_selection) for
                i in range(n_detectors)])
        self._time_amplitudes_averaged_selection = self.event_selection
        return self._time_amplitudes_averaged

    @property
    def time_scales(self):
        if len(self._time_scales) == 0:
            group = self.h5_file['time_scales']
            n_detectors = len(group.keys())
            self._time_scales = np.array(
                [group['det_{}'.format(i)].value for i in range(n_detectors)]
                ) * 1e3

        return self._time_scales

    @property
    def photo_slices(self):
        try:
            return self._photo_slices
        except:
            print 'ERROR: ROI for photoline not set yet.'
            sys.exit(1)

    @photo_slices.setter
    def photo_slices(self, x):
        try:
            iter(x[0])
        except:
            x = [x]*16

        self._photo_slices = []
        self._photo_bg_index_list = []
        self._photo_bg_factors = []
        for limits, t_axis in zip(x, self.time_scales):
            sl = slice(t_axis.searchsorted(limits[0]),
                       t_axis.searchsorted(limits[1], side='right'))
            bg_I = get_bg_index_list(sl, 5)
            self._photo_slices.append(sl)
            self._photo_bg_index_list.append(bg_I)
            self._photo_bg_factors.append(float(sl.stop-sl.start) / len(bg_I))
        self._photo_bg_factors = np.array(self._photo_bg_factors)

    @property
    def auger_slices(self):
        try:
            return self._auger_slices
        except:
            print 'ERROR: ROI for auger line not set yet.'
            sys.exit(1)

    @auger_slices.setter
    def auger_slices(self, x):
        try:
            iter(x[0])
        except:
            x = [x]*16

        self._auger_slices = []
        for limits, t_axis in zip(x, self.time_scales):
            self._auger_slices.append(
                slice(t_axis.searchsorted(limits[0]),
                      t_axis.searchsorted(limits[1], side='right')))

    def get_time_amplitudes(self, selections=[slice(None)]*16):
        names = map(lambda x: 'time_amplitudes/det_{}'.format(x), range(16))
        return [self.h5_file[name][:, selections[i]] for i, name in
                enumerate(names)]

    @property
    def photo_amplitudes(self):
        if len(self._photo_amplitudes) == 0:
            raw = np.array([det.sum(1) for det in
                            self.get_time_amplitudes(self._photo_slices)]).T
            bg = np.array(
                [det.sum(1) for det in
                 self.get_time_amplitudes(self._photo_bg_index_list)]).T
#            print 'raw:', raw.shape
#            print 'gb:', bg.shape
#            print 'factors:', self._photo_bg_factors.shape
            self._photo_amplitudes = ((raw - bg * self._photo_bg_factors) *
                                      self._det_factors)

        return self._photo_amplitudes

    @property
    def auger_amplitudes(self):
        if len(self._auger_amplitudes) == 0:
            self._auger_amplitudes = (np.array(
                [det.sum(1) for det in
                 self.get_time_amplitudes(self._auger_slices)]).T *
                self._det_factors)

        return self._auger_amplitudes

    @property
    def fee(self):
        return self.h5_file['fee'].value

    @property
    def det_factors(self):
        return self._det_factors

    @det_factors.setter
    def det_factors(self, new_val):
        self._det_factors = new_val

if __name__ == '__main__':
    dataset = DataSet('data/amom0115_35_1.h5')

    dataset.photo_slices = photo_roi
    dataset.auger_slices = auger_roi
#    dataset.print_content()

    plot_traces('Average traces', dataset)

    polar_plot(dataset, 'Polar', reset_scaling=True)
