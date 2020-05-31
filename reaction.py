from pre_proc import Datafile, Dataset
import pandas as pd
import numpy as np
import math as math
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.pyplot import rc_context
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

OPT=['#347f3a','#36358b','#e47327']
rc_fname='plotting_params'

class Reaction:

    def __init__(self,dataset,wavelength,epsilon,path_length):

        self.dataset = dataset
        self.wavelength = wavelength
        self.epsilon = epsilon
        self.path_length = path_length
        self.state='pre_proc'

    def nearest_value(self,value,value_arr):
        idx = np.searchsorted(value_arr, value, side="left")
        if idx > 0 and (idx == len(value_arr) or math.fabs(value- value_arr[idx-1]) < math.fabs(value - value_arr[idx])):
            return value_arr[idx-1]
        else:
            return value_arr[idx]

    def widen_range(self,value,value_arr,range):
        middle = self.nearest_value(value, value_arr)
        return self.nearest_value(value - range/2 , value_arr), self.nearest_value(value + range/2 , value_arr)

    def calculate_conc(self,range):

        abs_data=self.dataset.abs_data
        low_limit, high_limit = self.widen_range(self.wavelength,self.dataset.wavelengths,range)

        print('Calculating concentration using range {} nm to {} nm.'.format(low_limit, high_limit))

        temp = abs_data.loc[:,low_limit:high_limit]
        temp_avg = temp.mean(axis=1)
        self.conc_profile = temp_avg.div(self.epsilon*self.path_length*1e-6)
        return self.conc_profile

    def linear_drift(self,fit_times):
        temp = self.conc_profile
        fit_c_arr=[]
        fit_t_arr=[]

        for time in fit_times:
            fit_t_arr.append(self.nearest_value(time,self.dataset.times))

        fit_c_arr=self.conc_profile.loc[fit_t_arr]

        def drift_func(x,a,c):
            return a*x+c

        popt, pcov = curve_fit(drift_func, fit_t_arr, fit_c_arr)
        drift_delta = pd.Series(drift_func(self.dataset.times,popt[0],popt[1]),index=self.dataset.times)
        self.conc_profile_d = self.conc_profile.subtract(drift_delta)
        self.state = 'drift_corrected'
        return self.conc_profile_d

    def plot_conc(self):
        with rc_context(fname=rc_fname):
            fig = plt.figure(figsize=(8, 4))
            plt.plot(self.conc_profile,color=OPT[0],label='data')
            if self.state=='drift_corrected':
                plt.plot(self.conc_profile_d,color=OPT[1],label='data corrected')
                plt.legend()
            else:
                pass
            plt.xlim(0, max(self.dataset.times))
            plt.ylim(bottom=0)
            plt.axhline(0,color='gray',ls='--')
            plt.xlabel('time / s')
            plt.ylabel('concentration / $\mu$M')
            plt.title('Concentration Plot \n {}'.format(self.dataset.exp_label))
            plt.show()

    def export_conc(self):
        today = datetime.today()
        d1 = today.strftime("%Y-%b-%d-%H-%M")
        writer = pd.ExcelWriter('conc_export_{}.xlsx'.format(d1), engine='xlsxwriter')
        self.conc_profile.to_excel(writer, sheet_name='Sheet1')
        if self.state=='drift_corrected':
            self.conc_profile_d.to_excel(writer, sheet_name='Sheet2')
        else:
            pass
        writer.save()

    def find_turn_points(self, prom, translate=50, plot=False):
        if self.state=='drift_corrected':
            dx = self.conc_profile_d.diff()
        else:
            dx = self.conc_profile.diff()
        self.rate = dx.ewm(span = 50).mean()
        indices_ON = find_peaks(self.rate,prominence=prom)[0]
        indices_OFF = find_peaks(-self.rate,prominence=prom)[0]
        indices_all = sorted(np.concatenate((indices_ON, indices_OFF), axis=None))
        indices_all = np.array(indices_all)-translate

        if plot==True:
            with rc_context(fname=rc_fname):
                fig = plt.figure(figsize=(8, 4))
                plt.plot(self.conc_profile,label='data',color=OPT[0],linewidth=3)
                plt.plot(self.conc_profile.iloc[indices_all],marker='o',ls='none',ms=8,mfc='none',mec=OPT[2],mew=2,label='turning points')
                plt.legend()
                plt.xlim(0,max(self.dataset.times))
                plt.ylim(0,max(self.conc_profile)*1.05)
                plt.ylabel('Concentration / $\mu$M')
                plt.xlabel('Time / s')
                plt.show()
        else:
            pass
        print('returning the iloc indices for conc_profile')
        return indices_all
