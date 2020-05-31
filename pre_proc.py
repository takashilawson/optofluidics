import h5py
import pandas as pd
import os
import numpy as np
import math as math
from datetime import datetime, timedelta

class Datafile:

    def __init__(self, file_path):
        self.file_path = file_path
        measurement_list = []
        key_list = []
        if self.file_path.endswith(".h5"):
            data = h5py.File(self.file_path, 'r')
            for counter, measurement in enumerate(data.keys()):
                measurement_list.append(measurement)
                key_list.append(counter)
            self.exp_labels_list = measurement_list
            self.exp_key_list = key_list
            print('datafile intialised successfully \n')
        else:
            self.exp_labels_list = measurement_list
            self.exp_key_list = key_list
            print('error: the file is not a .h5 file \n')

    def list_groups(self):
        for counter, label in enumerate(self.exp_labels_list):
            print('Key {}: {} \n'.format(str(counter), label))

    def __repr__(self):
        return self.exp_key_list
        return self.exp_labels_list

class Dataset:

    def __init__(self, datafile, exp_key = 0, timelapse_key = 0):
        self.datafile = datafile
        self.exp_key = exp_key
        self.timelapse_key = timelapse_key
        self.exp_label = self.datafile.exp_labels_list[self.exp_key]
        self.h5_loc = r'/{}/timelapse_{}'.format(self.exp_label, str(self.timelapse_key))

        if self.exp_key >= max(self.datafile.exp_key_list):
            pass
        else:
            print('error: key out of range \n')

        hf = h5py.File(self.datafile.file_path, 'r')
        self.raw_data = hf.get(self.h5_loc)
        print ('{} dataset loaded successfully \n'.format(self.h5_loc))

    def pre_process(self):
        wav_arr_raw = np.array(self.raw_data['spectrum_0'].attrs['wavelengths'])
        self.wavelengths = wav_arr_raw
        self.back_spectra_arr = np.array(self.raw_data['spectrum_0'].attrs['background'])
        ref_spectra_raw = np.array(self.raw_data['spectrum_0'].attrs['reference'])
        self.ref_spectra_arr = np.subtract(ref_spectra_raw,self.back_spectra_arr)

        corr_data = []
        times_proc = []
        time_ref = str(self.raw_data['spectrum_0'].attrs['creation_timestamp'])

        try:
            time_ref = datetime.strptime((time_ref.replace('b','')).replace('\'',''),"%Y-%m-%dT%H:%M:%S.%f")
        except ValueError:
            time_ref = datetime.strptime((time_ref.replace('b','')).replace('\'',''),"%Y-%m-%dT%H:%M:%S")

        print('measurement was started at {}, \n normalising times and applying a background correction \n'.format(time_ref))

        for counter, spectra in enumerate(self.raw_data.keys()):
            corr_data.append(self.raw_data[spectra]-self.back_spectra_arr)
            time = str(self.raw_data[spectra].attrs['creation_timestamp'])
            try:
                time = datetime.strptime((time.replace('b','')).replace('\'',''),"%Y-%m-%dT%H:%M:%S.%f")
            except ValueError:
                time = datetime.strptime((time.replace('b','')).replace('\'',''),"%Y-%m-%dT%H:%M:%S")
            deltatime = time - time_ref
            times_proc.append(deltatime.total_seconds())

        self.times = np.array(times_proc)
        print('measurement contains {} spectra with {} wavelengths \n'.format(len(self.times),len(self.wavelengths)))
        pre_proc_data = pd.DataFrame(corr_data, index = self.times, columns = self.wavelengths)
        self.pre_proc_data = pre_proc_data.sort_index(axis=0)
        self.times = np.sort(self.times)
        return self.pre_proc_data

    def max_counts(self):
        return np.nanmax(self.pre_proc_data)

    def export_pre_proc(self, file_path, export_index = True, export_header = True):
        self.pre_proc_data.to_csv(file_path+'.csv', index = export_index, header = export_header)
        print('Pre-processed data saved to .csv succesfully. \n')

    def calculate_abs(self):
        abs=-np.log10(self.pre_proc_data.div(self.ref_spectra_arr))
        self.abs_data=abs
        return self.abs_data

    def max_abs(self):
        return np.nanmax(self.abs_data)

    def export_abs(self, file_path, export_index = True, export_header = True):
        self.abs_data.to_csv(os.path.join(file_path,'.csv'), index = export_index, header = export_header)
        print('Absorbance data saved to .csv succesfully. \n')

    def find_nearest_wav(self, wavelength):
        idx = np.searchsorted(self.wavelengths, wavelength, side="left")
        if idx > 0 and (idx == len(self.wavelengths) or math.fabs(wavelength - self.wavelengths[idx-1]) < math.fabs(wavelength - self.wavelengths[idx])):
            return self.wavelengths[idx-1]
        else:
            return self.wavelengths[idx]

    def find_nearest_time(self, time):
        idx = np.searchsorted(self.times, time, side="left")
        if idx > 0 and (idx == len(self.times) or math.fabs(time - self.times[idx-1]) < math.fabs(time - self.times[idx])):
            return self.times[idx-1]
        else:
            return self.times[idx]
