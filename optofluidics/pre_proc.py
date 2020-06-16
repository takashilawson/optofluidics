from datetime import datetime, timedelta
import h5py
import numpy as np
import math as math
import os
import pandas as pd
from time import perf_counter

class Datafile:

    """ Datafile class for viewing and selecting Datasets in a HDF5 file for
            processing.

	Attributes:
        exp_labels_list representing the list of experiment labels (list)
        file_path representing the file path to the HDF5 file of interest
            (string)

    """

    def __init__(self, file_path):

        """ Intialisation checks whether the file is a HDF5 file.

        """

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
            print('Datafile intialised successfully \n')
        else:
            self.exp_labels_list = measurement_list
            self.exp_key_list = key_list
            print('Error: the file is not a .h5 file \n')

    def list_groups(self):

        """ Function to returns the list of measurements in the datafile

        Args:
    		None

    	Returns:
    		string: list of measurements

        """

        for counter, label in enumerate(self.exp_labels_list):
            print('Key {}: {} \n'.format(str(counter), label))

    def __repr__(self):

        """Function to output the characteristics of the Datafile instance

		Args:
			None

		Returns:
			string: characteristics of the Datafile

		"""

        return self.exp_key_list
        return self.exp_labels_list

class Dataset:

    """ Dataset class for processing spectral timelapses.

	Attributes:
        abs_data representing the absorbance values (pd DataFrame)
        back_spectra_arr representing the background counts (array)
        datafile representing a Datafile object specified above
        exp_key representing the measurement you wish to access in the datafile
            (float)
        h5_loc representing the h5 file internal location of the timelapse
            (string)
        pre_proc_data representing the background-corrected counts
            (pd DataFrame)
        raw_data representing the counts of the spectra (HDF5 object)
        ref_spectra_arr representing the reference counts (array)
        timelapse_key representing the timelapse you wish to access (usually 0)
            (float)
        times representing the elapsed times of the spectra (array)
        wavelengths representing the wavelengths of the spectra (array)

    """

    def __init__(self, datafile, exp_key = 0, timelapse_key = 0):

        """Initialisation checks the validity of the key specified.

		Args:
			None

		Returns:
			None

		"""

        self.datafile = datafile
        self.exp_key = exp_key
        self.timelapse_key = timelapse_key
        self.exp_label = self.datafile.exp_labels_list[self.exp_key]

        # h5_loc is the internal HDF5 location reference to the timelapse of interest
        self.h5_loc = r'/{}/timelapse_{}'.format(self.exp_label, str(self.timelapse_key))

        # validation on whether the specified key exists in the datafile
        if self.exp_key <= max(self.datafile.exp_key_list):
            pass
        else:
            print('Error: key out of range \n')

        hf = h5py.File(self.datafile.file_path, 'r')
        self.raw_data = hf.get(self.h5_loc)
        print ('{} dataset loaded successfully \n'.format(self.h5_loc))

    def pre_process(self):

        """Function to pre-process the data to background correct the counts and
                calculate elapsed times from spectral timestamps.

		Args:
			None

		Returns:
			pd Dataframe: background-corrected data with elapsed times as
                indices and wavelength as columns

		"""
        t1_start = perf_counter()
        wav_arr_raw = np.array(self.raw_data['spectrum_0'].attrs['wavelengths'])
        self.wavelengths = wav_arr_raw
        self.back_spectra_arr = np.array(self.raw_data['spectrum_0'].attrs['background'])

        corr_data = []
        times_proc = []

        # extract reference point for 0 seconds
        time_ref = str(self.raw_data['spectrum_0'].attrs['creation_timestamp'])

        # spectrometer adds 'b' and quotation marks to timestamps that must be removed
        # some spectra are taken on X.000000s which does not have a .%f component - use try and except
        try:
            time_ref = datetime.strptime((time_ref.replace('b','')).replace('\'',''),"%Y-%m-%dT%H:%M:%S.%f")
        except ValueError:
            time_ref = datetime.strptime((time_ref.replace('b','')).replace('\'',''),"%Y-%m-%dT%H:%M:%S")

        print('Measurement was started at {}, \n normalising times and applying a background correction \n'.format(time_ref))

        # applies background correction
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
        print('Measurement contains {} spectra with {} wavelengths \n'.format(len(self.times),len(self.wavelengths)))

        # data is stored as a pd Dataframe with elapsed times as indices and wavelengths as columns
        pre_proc_data = pd.DataFrame(corr_data, index = self.times, columns = self.wavelengths)

        # data may be disordered in time when iterated through
        # sort the data by elapsed time
        self.pre_proc_data = pre_proc_data.sort_index(axis=0)
        self.times = np.sort(self.times)

        t1_stop = perf_counter()
        print("Elapsed time for pre-processing:", t1_stop-t1_start)

        return self.pre_proc_data

    def max_counts(self):

        """Function to return the maximum counts in the timelapse

		Args:
			None

		Returns:
			float: maximum counts

		"""

        return np.nanmax(self.pre_proc_data)

    def preproc2csv(self, file_path, export_index = True, export_header = True):

        """Function to export pre-processed data to a csv file

		Args:
			file_path representing the csv file location to save to
            export_index representing a Boolean on whether to export the times
                as the first column
            export_header representing a Boolean on whether to export the
                wavelengths as the first row

		Returns:
			None

		"""

        self.pre_proc_data.to_csv(file_path+'.csv', index = export_index, header = export_header)
        print('Pre-processed data saved to .csv succesfully. \n')

    def calculate_abs(self):

        """Function to extract reference spectrum and calculate absorbances

		Args:
			None

		Returns:
			pd Dataframe: absorbance data with elapsed times as indices and
                wavelength as columns

		"""
        ref_spectra_raw = np.array(self.raw_data['spectrum_0'].attrs['reference'])
        self.ref_spectra_arr = np.subtract(ref_spectra_raw,self.back_spectra_arr)
        abs=-np.log10(self.pre_proc_data.div(self.ref_spectra_arr))
        self.abs_data=abs
        return self.abs_data

    def max_abs(self):

        """Function to return the maximum absorbance in the timelapse

		Args:
			None

		Returns:
			float: maximum absorbance value

		"""

        return np.nanmax(self.abs_data)

    def abs2csv(self, file_name, file_path, export_index = True, export_header = True):

        """Function to export absorbance data to a csv file

		Args:
			file_name representing the file name with file extension (string)
            file_path representing the folder to save tgo (string)
            export_index representing a Boolean on whether to export the times
                as the first column
            export_header representing a Boolean on whether to export the
                wavelengths as the first row

		Returns:
			None

		"""

        self.abs_data.to_csv(os.path.join(file_path,file_name), index = export_index, header = export_header)
        print('Absorbance data saved to .csv succesfully. \n')

    def find_nearest_wav(self, wavelength):

        """Function to find the nearest wavelength in the timelapse to the
                wavelength specified (float)

		Args:
		      wavelength representing the wavelength you wish to search for

		Returns:
			Float: closest wavelength in timelapse dataset

		"""

        idx = np.searchsorted(self.wavelengths, wavelength, side="left")
        if idx > 0 and (idx == len(self.wavelengths) or math.fabs(wavelength - self.wavelengths[idx-1]) < math.fabs(wavelength - self.wavelengths[idx])):
            return self.wavelengths[idx-1]
        else:
            return self.wavelengths[idx]

    def find_nearest_time(self, time):

        """Function to find the nearest time in the timelapse to the
                time specified

		Args:
		      time representing the time you wish to search for (float)

		Returns:
			Float: closest time in timelapse dataset

		"""

        idx = np.searchsorted(self.times, time, side="left")
        if idx > 0 and (idx == len(self.times) or math.fabs(time - self.times[idx-1]) < math.fabs(time - self.times[idx])):
            return self.times[idx-1]
        else:
            return self.times[idx]

    def __repr__(self):

        """Function to output the characteristics of the Dataset instance

		Args:
			None

		Returns:
			string: characteristics of the Dataset

		"""

        return self.datafile.file_path
