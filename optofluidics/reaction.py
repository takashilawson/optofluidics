from .pre_proc import Datafile, Dataset

from datetime import datetime, timedelta
from lmfit import minimize, Parameters, fit_report
import numpy as np
import math as math
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.pyplot import rc_context
import os
import pandas as pd
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from scipy.interpolate import interp1d
import scipy.special as scs
from time import perf_counter

# specify the three colours to use in plots
OPT=['#347f3a','#36358b','#e47327']

# specify a rc_file which contains plotting parameters
plot_config_file = 'plotting_params.txt'
dir = os.path.abspath(os.path.dirname(__file__))
rc_fname=os.path.join(dir, plot_config_file)

class Reaction:

    """ Reaction class for converting spectral information into chemical
            information and for post-processing procedures.

	Attributes:
        conc_profile representing the concentration profile (pd DataFrame)
        conc_profile_d representing the drift-corrected concentration profile
            (pd DataFrame)
        cycles representing the number of UV on-off periods (float)
        dataset representing a Dataset object specified in pre_proc.py
        epsilon representing the absorption coefficient for Beer-Lambert Law
            calculations (float)
        model representing an initial model (pd DataFrame)
        path_length representing the cm path length for Beer-Lambert Law
            calculations (float)
        rates representing the time-dependent rates (pd DataFrame)
        state representing the status of the data (string)
        turning_points representing the times of turning points (list)
        wavelength representing the wavelength for Beer-Lambert Law calculations
            (float)

    """

    def __init__(self,dataset,wavelength,epsilon,path_length):

        self.dataset = dataset
        self.wavelength = wavelength
        self.epsilon = epsilon
        self.path_length = path_length
        self.state='pre_proc'

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

    def nearest_value_idx(self,value,value_arr):

        """Function to find the nearest value in a value array

		Args:
			value representing the value you wish to find (float)
            value_array representing the array you wish to search through (array)

		Returns:
			float: nearest value to that specfied

		"""

        idx = np.searchsorted(value_arr, value, side="left")
        if idx > 0 and (idx == len(value_arr) or math.fabs(value- value_arr[idx-1]) < math.fabs(value - value_arr[idx])):
            return idx-1
        else:
            return idx

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

    def calculate_conc(self,range):

        """Function to return the concentration profile

		Args:
			range representing the range either side of the central wavelength to
                apply Beer-Lambert Law calculations to (float)

		Returns:
			pd Dataseries: concentration profile with elapsed times as indices

		"""

        abs_data=self.dataset.abs_data
        low_limit, high_limit = self.widen_range(self.wavelength,self.dataset.wavelengths,range)

        print('Calculating concentration using range {} nm to {} nm.'.format(low_limit, high_limit))

        # temp stores the absorbance values in the wavelength range
        temp = abs_data.loc[:,low_limit:high_limit]
        print('Absorbance will be averaged over {} data points'.format((np.array(temp)).shape[1]))

        # average Beer-Lambert calculations at each wavelength over the
        # wavelength range
        temp_avg = temp.mean(axis=1)
        self.conc_profile = temp_avg.div(self.epsilon*self.path_length*1e-6)
        return self.conc_profile

    def linear_drift(self,fit_times):

        """Function to return a drift corrected concentration profile (linear)

		Args:
			fit_times representing the times at which the absorbance should be
                nil (array)

		Returns:
			pd Dataseries: concentration profile with elapsed times as indices

		"""

        temp = self.conc_profile
        fit_c_arr=[]
        fit_t_arr=[]

        for time in fit_times:
            fit_t_arr.append(self.nearest_value(time,self.dataset.times))

        fit_c_arr=self.conc_profile.loc[fit_t_arr]

        # fit to linear function
        def drift_func(x,a,c):
            return a*x+c

        popt, pcov = curve_fit(drift_func, fit_t_arr, fit_c_arr)
        drift_delta = pd.Series(drift_func(self.dataset.times,popt[0],popt[1]),index=self.dataset.times)
        self.conc_profile_d = self.conc_profile.subtract(drift_delta)
        self.state = 'drift_corrected'
        return self.conc_profile_d

    def plot_conc(self):

        """Function to return a plot of the concentration profile. Drift corrected
                data will be shown if calculated.

		Args:
			None

		Returns:
			plot

		"""

        with rc_context(fname=rc_fname):
            fig = plt.figure(figsize=(8, 4))
            plt.plot(self.conc_profile,color=OPT[0],label='data')

            # check to see if drift-corrected data exists
            if self.state=='drift_corrected':
                plt.plot(self.conc_profile_d,color=OPT[1],label='data corrected')
                plt.legend(frameon='True')
            else:
                pass

            plt.xlim(0, max(self.dataset.times))
            plt.ylim(0, max(self.conc_profile.values)+2)
            plt.axhline(0,color='gray',ls='--')
            plt.xlabel('time / s')
            plt.ylabel('concentration / $\mu$M')
            plt.title('Concentration Plot \n {}'.format(self.dataset.exp_label))
            plt.show()

    def export_conc(self):

        """Function to export concentration profile to xlsx. Drift corrected data
                will be exported if it exists.

		Args:
			None

		Returns:
			None

		"""

        today = datetime.today()
        d1 = today.strftime("%Y-%b-%d-%H-%M")

        # save xlsx with today's date and time in file name
        writer = pd.ExcelWriter('conc_export_{}.xlsx'.format(d1), engine='xlsxwriter')
        self.conc_profile.to_excel(writer, sheet_name='Sheet1')

        # export drift-corrected data if exists
        if self.state=='drift_corrected':
            self.conc_profile_d.to_excel(writer, sheet_name='Sheet2')
        else:
            pass
        writer.save()

    def find_turn_points(self, first, prom, translate=50, plot=False):

        """Function to find turning points in concentration profile

		Args:
			first representing the first UV on time (first turning point doesn't
                account for delay time, so this must be manually specified)
            prom representing the prominence of peaks (float)
            translate representing a shift correction in time for identified
                peaks (float)
            plot representing a Boolean on whether a plot should be returned

		Returns:
			array representing the iloc indices for the turning points

		"""
        # perform turning point identification on drift-corrected data if it
        # exists
        if self.state=='drift_corrected':
            dx = self.conc_profile_d.diff()
        else:
            dx = self.conc_profile.diff()

        # apply exponentially weighted mean (ewm) with 50 data point span to
        # smoothen data before peak identification
        self.rate = dx.ewm(span = 50).mean()
        indices_ON = find_peaks(self.rate,prominence=prom)[0]
        indices_OFF = find_peaks(-self.rate,prominence=prom)[0]
        indices_all = sorted(np.concatenate((indices_ON, indices_OFF), axis=None))
        indices_all = np.array(indices_all)-translate
        cycles = len(indices_all)/2
        self.cycles = cycles
        print('{} cycles detected'.format(self.cycles))

        if len(indices_all) % 2 == 0:
            pass
        else:
            print('Warning: Number of cycles invalid. Adjust prominence.')
        if cycles > 3:
            print('Warning: More than three cycles.')
        else:
            pass

        indices_all = np.delete(indices_all,0)

        if plot==True:
            with rc_context(fname=rc_fname):
                fig = plt.figure(figsize=(8, 4))
                plt.plot(self.conc_profile,label='data',color=OPT[0],linewidth=3)
                plt.plot(self.conc_profile.iloc[indices_all],marker='o',ls='none',ms=8,mfc='none',mec=OPT[2],mew=2,label='turning points')
                plt.legend(frameon='True')
                plt.xlim(0,max(self.dataset.times))
                plt.ylim(0,max(self.conc_profile)*1.05)
                plt.ylabel('Concentration / $\mu$M')
                plt.xlabel('Time / s')
                plt.show()
        else:
            pass

        times_all = []
        for iloc in indices_all:
            times_all.append(int(self.dataset.times[iloc]))

        times_all = np.insert(times_all, 0, first)

        self.turning_points = times_all
        return times_all

    def update_turn_points(self,times_arr):

        """Function to update turning points array

		Args:
            times_arr representing the times at which turning points exist

        Returns:
			None

		"""
        self.turning_points = times_arr

    def create_model(self,pars):

        """Function to create model based on rate constants

		Args:
            params representing the model parameters consisting of
                {
                t_erf representing the delay time in Erf model
                kc representing the complexation rate constant
                kbr representing the product bleaching rate constant
                k representing the saturation rate constant in the Erf model
                kr representing the back-reaction rate constant
                o representing the standard deviation in the Erf model
                c0 representing the initial concentration of reactant
                K representing the equilibrium constant for complexation
                end representing the total time for the model
                }

        Returns:
			pd DataFrame

		"""

        if hasattr(self, 'turning_points'):
            pass
        else:
            print("Run find_turn_points or update_turn_points first.")

        params = pars.valuesdict()
        t_erf = params['t_erf']
        kc = params['kc']
        kbr = params['kbr']
        k = params['k']
        kr = params['kr']
        o = params['o']
        c0 = params['c0']
        K = params['K']
        end = params['end']

        times_arr = self.turning_points

        if self.cycles == len(times_arr)/2:
            pass
        else:
            print('Error: Number of cycles does not match up with time array.')

        # creation of dictionary with all concentration values
        Concentrations = {}

        # setting concentration start values at t=0
        Concentrations["MV"] = [c0]
        Concentrations["R"] = [0]
        Concentrations["C"] = [0]
        Concentrations["B"] = [0]

        # create time array for model
        t_model = np.linspace(0,end,end+1,dtype='int64')

        # create arrays to track UV on and off periods
        t_uv = [0]
        uv_bool = [0]

        # populate arrays to track UV on and off periods
        uv_bool_test = False # UV is off initially
        for time in range(1,len(t_model)):
            # if time matches times_arr we are at a switching point
            # update uv_bool_test to opposite
            if time in times_arr:
                uv_bool_test = not uv_bool_test
            else:
                pass
            # when UV is off, t_uv is constant and uv_bool is 0 (off)
            if uv_bool_test is False:
                t_uv.append(t_uv[time-1])
                uv_bool.append(0)

            # when UV is on, t_uv increases by 1 and uv_bool is 1 (on)
            else:
                t_uv.append(t_uv[time-1]+1)
                uv_bool.append(1)

        # define functions for rate constants
        # plus t_erf here as t_uv is referenced to turning points NOT UV on/off
        # periods
        def rate_k(time):
            return uv_bool[time]*(0.5*k*(1+scs.erf((t_uv[time]-t_erf)/(np.sqrt(2)*o))))

        def rate_kr(time):
            return kr

        def rate_kc(time):
            # need to equilibrate complexation at start, 0.1 returns no errors.
            if time < 60 and K !=0:
                return 0.1
            elif time >= 60 and K !=0:
                return kc
            else:
                return 0

        def rate_kcr(time):
            # need to equilibrate complexation at start, 0.1 returns no errors.
            if time < 60 and K != 0:
                return 0.1*(1/K)
            elif time >= 60 and K != 0:
                return rate_kc(time)*(1/K)
            else:
                return 0

        def rate_kbr(time):
            return uv_bool[time]*kbr

        # creation of dictionary with all rate values
        Rates = {}

        # intialise rates dictionary
        Rates['k(t)'] = []
        Rates['kr'] = []
        Rates['kc'] = []
        Rates['kcr'] = []
        Rates['kbr'] = []

        # populate rates dictionary
        for time in t_model:
            Rates['k(t)'].append(rate_k(time))
            Rates['kr'].append(rate_kr(time))
            Rates['kbr'].append(rate_kbr(time))
            Rates['kcr'].append(rate_kcr(time))
            Rates['kc'].append(rate_kc(time))

        self.rates = pd.DataFrame(Rates, index=t_model)

        t_model=[0]
        for i in range(1,end):
            # the concentration values of the prevous time step are read-out and
            # temporarliy saved
            MV = Concentrations["MV"][i-1]
            R = Concentrations["R"][i-1]
            B = Concentrations["B"][i-1]
            C = Concentrations["C"][i-1]

            # the rate constants of the prevous time step are read-out and
            # temporarliy saved
            kt1 = Rates['k(t)'][i-1]
            kr1 = Rates['kr'][i-1]
            kc1 = Rates['kc'][i-1]
            kcr1 = Rates['kcr'][i-1]
            kbr1 = Rates['kbr'][i-1]

            # kinetic equations are used to calculate the next concentration
            # value from past concentration values and past rate constants
            Concentrations["MV"].append(MV-kt1*MV-kc1*MV+kr1*R+kcr1*C)
            Concentrations["C"].append(C+kc1*MV-kcr1*C)
            Concentrations["R"].append(R-kr1*R-kbr1*R+kt1*MV)
            Concentrations["B"].append(B+kbr1*R)
            t_model.append(i)

        self.model = pd.DataFrame(Concentrations, index=t_model)
        return pd.DataFrame(Concentrations, index=t_model)

    def get_rates(self):

        """Function to return all time-dependent rates

		Args:
            None

		Returns:
			pd DataFrame containing the time-dependent rates

		"""

        return self.rates

    def model2csv(self, file_name, file_path, export_index = True, export_header = True):

        """Function to export model to a csv file

		Args:
			file_name representing the file name with csv file extension (string)
            file_path representing the folder to save to (string)
            export_index representing a Boolean on whether to export the times
                as the first column
            export_header representing a Boolean on whether to export the
                components as the first row

		Returns:
			None

		"""

        self.model.to_csv(os.path.join(file_path,file_name), index = export_index, header = export_header)
        print('Model data saved to .csv succesfully. \n')

    def residual(self,params):

        """Function to return residual to model

		Args:
            params representing the model parameters consisting of
                {
                t_erf representing the delay time in Erf model
                kc representing the complexation rate constant
                kbr representing the product bleaching rate constant
                k representing the saturation rate constant in the Erf model
                kr representing the back-reaction rate constant
                o representing the standard deviation in the Erf model
                c0 representing the initial concentration of reactant
                K representing the equilibrium constant for complexation
                end representing the total time for the model
                }
            complex_bool representing a Boolean on whether to account for
                complexation

        Returns:
			pd Series representing the residuals

		"""
        model_all = self.create_model(params)
        model_R = model_all["R"]

        f = interp1d(model_R.index.values, model_R.values)
        model_inter1pd = f(self.dataset.times)

        # create interpolation function for model
        if self.state == 'drift_corrected':
            exp_data = self.conc_profile_d
        else:
            exp_data = self.conc_profile

        residual = model_inter1pd-exp_data
        return np.array(residual.values)

    def fit_model(self,params,meth,save_path):

        """Function to fit experimental data to the model

		Args:
            params representing the initial fitting parameters consisting of
                {
                t_erf representing the delay time in Erf model
                kc representing the complexation rate constant
                kbr representing the product bleaching rate constant
                k representing the saturation rate constant in the Erf model
                kr representing the back-reaction rate constant
                o representing the standard deviation in the Erf model
                c0 representing the initial concentration of reactant
                K representing the equilibrium constant for complexation
                end representing the total time for the model
                }
            method representing the lmfit algorithm to use (string)
            save_path representing the file path to save the fit report to.

        Returns:
			lmfit minimise object

		"""
        t1_start = perf_counter()
        out = minimize(self.residual, params, method=meth)
        t1_stop = perf_counter()

        with open(os.path.join(save_path,'{}-fit-{}-report.txt'.format(self.dataset.exp_label,meth)), 'w') as fh:
            fh.write(fit_report(out))
        fh.close()

        print("Elapsed time for fitting:", t1_stop-t1_start,"\n")
        print(fit_report(out))
        print("\n The model attribute has been updated with these parameters. ")
        return out

    def update_model_range(self,range):
        i1=self.nearest_value_idx(range[0],self.dataset.times)
        i2=self.nearest_value_idx(range[1],self.dataset.times)
        self.i1=i1
        self.i2=i2
        return(i1,i2)

    def residual_range(self,params):

        """Function to return residual to model

		Args:
            range representing the fitting range (array)
            params representing the model parameters consisting of
                {
                t_erf representing the delay time in Erf model
                kc representing the complexation rate constant
                kbr representing the product bleaching rate constant
                k representing the saturation rate constant in the Erf model
                kr representing the back-reaction rate constant
                o representing the standard deviation in the Erf model
                c0 representing the initial concentration of reactant
                K representing the equilibrium constant for complexation
                end representing the total time for the model
                }
            complex_bool representing a Boolean on whether to account for
                complexation

        Returns:
			pd Series representing the residuals

		"""

        model_all = self.create_model(params)
        model_R = model_all["R"]

        f = interp1d(model_R.index.values, model_R.values)
        model_inter1pd = f(self.dataset.times[self.i1:self.i2])

        # create interpolation function for model
        if self.state == 'drift_corrected':
            exp_data = self.conc_profile_d.iloc[self.i1:self.i2]
        else:
            exp_data = self.conc_profile.iloc[self.i1:self.i2]

        residual = model_inter1pd-exp_data
        return np.array(residual.values)

    def fit_model_range(self,params,meth,save_path):

        """Function to fit experimental data to the model

		Args:
            range representing the fitting range (array)
            params representing the initial fitting parameters consisting of
                {
                t_erf representing the delay time in Erf model
                kc representing the complexation rate constant
                kbr representing the product bleaching rate constant
                k representing the saturation rate constant in the Erf model
                kr representing the back-reaction rate constant
                o representing the standard deviation in the Erf model
                c0 representing the initial concentration of reactant
                K representing the equilibrium constant for complexation
                end representing the total time for the model
                }
            method representing the lmfit algorithm to use (string)
            save_path representing the file path to save the fit report to.

        Returns:
			lmfit minimise object

		"""

        t1_start = perf_counter()
        out = minimize(self.residual_range, params, method=meth)
        t1_stop = perf_counter()

        with open(os.path.join(save_path,'{}-fit-{}-report.txt'.format(self.dataset.exp_label,meth)), 'w') as fh:
            fh.write(fit_report(out))
        fh.close()

        print("Elapsed time for fitting:", t1_stop-t1_start,"\n")
        print(fit_report(out))
        print("\n The model attribute has been updated with these parameters. ")
        return out

    def __repr__(self):

        """Function to output the characteristics of the Reaction instance

		Args:
			None

		Returns:
			string: characteristics of the Reaction

		"""

        return str(self.dataset.exp_label)
