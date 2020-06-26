from .pre_proc import Datafile, Dataset

from datetime import datetime, timedelta
import numpy as np
import math as math
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.pyplot import rc_context
import os
import pandas as pd
from time import perf_counter

# specify the three colours to use in plots
OPT=['#347f3a','#36358b','#e47327']

# specify a rc_file which contains plotting parameters
plot_config_file = 'plotting_params.txt'
dir = os.path.abspath(os.path.dirname(__file__))
rc_fname=os.path.join(dir, plot_config_file)

class SV:

    """ Stern Volmer class for converting spectral information for Stern Volmer
            analysis.

	Attributes:
        abs representing the absorbance of the solution at the specified
            wavelength
        counts representing the non-normalised fluorescence counts (float)
        dataset_f representing a Dataset object specified in pre_proc.py for
            the sample fluorescence.
        dataset_a representing a Dataset object specified in pre_proc.py for
            the sample absorbance.
        ncounts representing the normalised fluorescence counts (float)
        normalise representing a Boolean on whether to normalise counts or not
        quencher representing the quencher concentration (float)
        state representing the status of the data (string)
        wavelength representing the peak fluorescence wavelength (float)

    """

    def __init__(self,dataset_f,dataset_a,wavelength,quencher,normalise=True):

        self.dataset_f = dataset_f
        self.dataset_a = dataset_a
        self.normalise = normalise
        self.quencher = quencher
        self.state = 'pre_proc'
        self.wavelength = wavelength

    def nearest_value(self,value,value_arr):

        """Function to find the nearest value in a value array

		Args:
			value representing the value you wish to find (float)
            value_array representing the array you wish to search through (array)

		Returns:
			float: nearest value to that specfied

		"""

        idx = np.searchsorted(value_arr, value, side="left")
        if idx > 0 and (idx == len(value_arr) or math.fabs(value- value_arr[idx-1]) < math.fabs(value - value_arr[idx])):
            return value_arr[idx-1]
        else:
            return value_arr[idx]

    def widen_range(self,value,value_arr,range):

        """Function to return the nearest values for the specified
                range defined by [Value-Range/2, Value + Range/2]

		Args:
			value representing the value you wish to find (float)
            value_array representing the array you wish to search through (array)
            range representing the range you wish to expand to (float)

		Returns:
			Two floats: lowest value and highest value for range

		"""

        middle = self.nearest_value(value, value_arr)
        return self.nearest_value(value - range/2 , value_arr), self.nearest_value(value + range/2 , value_arr)

    def calculate_counts(self,range):

        """Function to return the absorbance-normalised counts

		Args:
			range representing the range either side of the central wavelength to
                average counts over (float)

		Returns:
			float: averaged and absorbance-normalised counts

		"""

        counts_f_arr=self.dataset_f.pre_proc_data
        low_limit, high_limit = self.widen_range(self.wavelength,self.dataset_f.wavelengths,range)

        print('Calculating average counts using range {} nm to {} nm.'.format(low_limit, high_limit))

        # temp stores the count values in the wavelength range
        temp = counts_f_arr.loc[:,low_limit:high_limit]
        print('Counts will be averaged over {} data points'.format((np.array(temp)).shape[1]))

        # average counts over the wavelength range and over time
        counts_f = temp.mean()
        self.counts = counts_f

        # check whether counts should be normalised
        if self.normalise is True:

            # apply similar process to absorbance data
            abs_a_arr=self.dataset_a.abs_data
            low_limit, high_limit = self.widen_range(self.wavelength,self.dataset_a.wavelengths,range)
            temp = abs_a_arr.loc[:,low_limit:high_limit]
            abs_a = temp.mean()
            self.abs = abs_a

        else:
            self.abs = 1

        self.ncounts = self.counts/self.abs
        self.state = 'calculated'

        return self.ncounts

    def __repr__(self):

        """Function to output the characteristics of the Stern Volmer instance

		Args:
			None

		Returns:
			string: characteristics of the Stern Volmer instance

		"""

        return str(self.dataset.exp_label)
