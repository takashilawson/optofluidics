from .pre_proc import Datafile, Dataset
from .reaction import Reaction

from datetime import datetime, timedelta
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.pyplot import rc_context
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import os
import pandas as pd

# specify three colours to use in plots
OPT=['#347f3a','#36358b','#e47327']

# specify matplotlib configuration file name
plot_config_file = 'plotting_params.txt'
dir = os.path.abspath(os.path.dirname(__file__))
rc_fname=os.path.join(dir, plot_config_file)

def produce_colour_arr(arr):

    """Function to produce a corresponding colour array for a dataset.

	Args:
		arr representing the data array you want to produce a colour array for

	Returns:
		array representing the colour array

	"""

    evenly_spaced_interval = np.linspace(0, 1, len(arr))

    #cividis colour map used for plotting
    colors = [cm.cividis(x) for x in evenly_spaced_interval]

    return colors

def plot_counts(dataset):

    """Function to plot background-corrected counts at all times

	Args:
		dataset representing the Dataset object

	Returns:
		none

	"""

    temp = dataset.pre_proc_data
    temp = temp.transpose()
    colors = produce_colour_arr(dataset.times)

    with rc_context(fname=rc_fname):
        for i, time in enumerate(temp):
            plt.plot(temp[time], color = colors[i], linewidth = 0.5)

        plt.xlim(400,900)
        plt.ylim(0,dataset.max_counts()*1.05)
        plt.xlabel('wavelength / nm')
        plt.ylabel('counts')
        plt.title('Spectral Plot \n {}'.format(dataset.exp_label))
        plt.show()

def plot_trace(dataset):

    """Function to plot referenced counts at all times

	Args:
		dataset representing the Dataset object

	Returns:
		none

	"""

    temp = dataset.pre_proc_data
    trace = temp.div(dataset.ref_spectra_arr)
    trace = trace.transpose()

    colors = produce_colour_arr(dataset.times)

    with rc_context(fname=rc_fname):
        for i, time in enumerate(trace):
            plt.plot(trace[time], color = colors[i], linewidth = 0.5)

        plt.xlim(400,900)
        plt.ylim(0,2)
        plt.xlabel('wavelength / nm')
        plt.ylabel('counts')
        plt.title('Referenced Spectral Plot \n {}'.format(dataset.exp_label))
        plt.show()

def plot_abs(dataset,absMIN,absMAX):

    """Function to plot absorbance at all times

	Args:
		dataset representing the Dataset object
        absMIN representing the minimum absorbance to plot
        absMAX representing the maximum absorbance to plot

	Returns:
		none

	"""

    temp = dataset.abs_data
    temp = temp.transpose()
    colors = produce_colour_arr(dataset.times)

    with rc_context(fname=rc_fname):
        for i, time in enumerate(temp):
            plt.plot(temp[time], color = colors[i], linewidth = 0.5)

        plt.xlim(400,900)
        plt.ylim(absMIN,absMAX)
        plt.xlabel('wavelength / nm')
        plt.ylabel('absorbance')
        plt.title('Absorbance Plot \n {}'.format(dataset.exp_label))
        plt.show()

def wav_counts(dataset,wavelength_list,plot=False):

    """Function to return background-corrected counts for specified wavelengths
        over time

	Args:
		dataset representing the Dataset object
        wavelength_list representing a list of wavelengths (list)
        plot representing a Boolean on whether to create a plot (Boolean)

	Returns:
		pd DataFrame with background-corrected counts at each wavelength

	"""

    temp = dataset.pre_proc_data
    wav_plt = []

    for wavelength in wavelength_list:
        nearest_wav = dataset.find_nearest_wav(wavelength)
        wav_plt.append(nearest_wav)
        print('Will return {} nm in dataset'.format(np.round(nearest_wav,2)))

    if plot==True:
        colors = produce_colour_arr(wav_plt)
        with rc_context(fname=rc_fname):
            for i, wav in enumerate(wav_plt):
                plt.plot(temp.loc[:,wav],color=colors[i],label='{} nm'.format(round(wav)))

            plt.legend(frameon='True')
            plt.xlim(0,max(dataset.times))
            plt.ylim(bottom=0)
            plt.xlabel('time / s')
            plt.ylabel('counts')
            plt.title('Spectral Plot \n {}'.format(dataset.exp_label))
            plt.show()

    else:
        pass

    return temp.loc[:,wav_plt]

def wav_trace(dataset,wavelength_list,plot=False):

    """Function to return referenced counts for specified wavelengths
        over time

	Args:
		dataset representing the Dataset object
        wavelength_list representing a list of wavelengths (list)
        plot representing a Boolean on whether to create a plot (Boolean)

	Returns:
		pd DataFrame with referenced counts at each wavelength

	"""

    temp = dataset.pre_proc_data
    trace = temp.div(dataset.ref_spectra_arr)

    wav_plt = []

    for wavelength in wavelength_list:
        nearest_wav = dataset.find_nearest_wav(wavelength)
        wav_plt.append(nearest_wav)
        print('Will return {} nm in dataset'.format(np.round(nearest_wav,2)))

    if plot==True:
        colors = produce_colour_arr(wav_plt)
        with rc_context(fname=rc_fname):
            for i, wav in enumerate(wav_plt):
                plt.plot(trace.loc[:,wav],color=colors[i],label='{} nm'.format(round(wav)))

            plt.legend(frameon='True')
            plt.xlim(0,max(dataset.times))
            plt.ylim(bottom=0)
            plt.xlabel('wavelength / nm')
            plt.ylabel('counts')
            plt.title('Referenced Spectral Plot \n {}'.format(dataset.exp_label))
            plt.show()

    else:
        pass

    return trace.loc[:,wav_plt]

def wav_abs(dataset,wavelength_list,absMIN,absMAX,plot=False):

    """Function to return absorbance for specified wavelengths over time

	Args:
		dataset representing the Dataset object
        wavelength_list representing a list of wavelengths (list)
        absMIN representing the minimum absorbance to plot
        absMAX representing the maximum absorbance to plot
        plot representing a Boolean on whether to create a plot (Boolean)

	Returns:
		pd DataFrame with absorbance at each wavelength

	"""

    temp=dataset.abs_data
    wav_plt = []

    for wavelength in wavelength_list:
        nearest_wav = dataset.find_nearest_wav(wavelength)
        wav_plt.append(nearest_wav)
        print('Will return {} nm in dataset'.format(np.round(nearest_wav,2)))

    if plot==True:
        colors = produce_colour_arr(wav_plt)

        with rc_context(fname=rc_fname):
            for i, wav in enumerate(wav_plt):
                plt.plot(temp.loc[:,wav],color=colors[i],label='{} nm'.format(round(wav)))

            plt.legend(frameon='True')
            plt.xlim(0,max(dataset.times))
            plt.ylim(absMIN,absMAX)
            plt.xlabel('time / s')
            plt.ylabel('absorbance')
            plt.title('Absorbance Plot \n {}'.format(dataset.exp_label))
            plt.show()

    else:
        pass

    return temp.loc[:,wav_plt]

def time_counts(dataset,time_list,plot=False):

    """Function to return background-corrected spectra at specified times

	Args:
		dataset representing the Dataset object
        time_list representing a list of times (list)
        plot representing a Boolean on whether to create a plot (Boolean)

	Returns:
		pd DataFrame with background-corrected spectra at each time

	"""

    temp=dataset.pre_proc_data
    time_plt = []

    for time in time_list:
        nearest_time = dataset.find_nearest_time(time)
        time_plt.append(nearest_time)
        print('Will plot {} s in dataset'.format(np.round(nearest_time)))

    colors = produce_colour_arr(time_plt)

    if plot==True:
        with rc_context(fname=rc_fname):
            for i, time in enumerate(time_plt):
                plt.plot(temp.loc[time,:],color=colors[i],label='{} s'.format(round(time)))

            plt.legend(frameon='True')
            plt.xlim(400,900)
            plt.ylim(bottom=0)
            plt.xlabel('wavelength / nm')
            plt.ylabel('counts')
            plt.title('Spectral Plot \n {}'.format(dataset.exp_label))
            plt.show()

    else:
        pass

    return temp.loc[time_plt,:]

def time_trace(dataset,time_list,plot=False):

    """Function to return referenced spectra at specified times

	Args:
		dataset representing the Dataset object
        time_list representing a list of times (list)
        plot representing a Boolean on whether to create a plot (Boolean)

	Returns:
		pd DataFrame with referenced spectra at each time

	"""

    temp=dataset.pre_proc_data
    trace = temp.div(dataset.ref_spectra_arr)

    time_plt = []

    for time in time_list:
        nearest_time = dataset.find_nearest_time(time)
        time_plt.append(nearest_time)
        print('Will plot {} s in dataset'.format(np.round(nearest_time)))

    if plot==True:
        colors = produce_colour_arr(time_plt)

        with rc_context(fname=rc_fname):
            for i, time in enumerate(time_plt):
                plt.plot(trace.loc[time,:],color=colors[i],label='{} s'.format(round(time)))

            plt.legend(frameon='True')
            plt.xlim(400,900)
            plt.ylim(bottom=0)
            plt.xlabel('wavelength / nm')
            plt.ylabel('counts')
            plt.title('Referenced Spectral Plot \n {}'.format(dataset.exp_label))
            plt.show()
    else:
        pass

    return trace.loc[time_plt,:]

def time_abs(dataset,time_list,absMIN,absMAX,plot=False):

    """Function to return absorbance spectra at specified times

	Args:
		dataset representing the Dataset object
        time_list representing a list of times (list)
        absMIN representing the minimum absorbance to plot
        absMAX representing the maximum absorbance to plot
        plot representing a Boolean on whether to create a plot (Boolean)

	Returns:
		pd DataFrame with absorbance spectra at each time

	"""

    temp=dataset.abs_data
    time_plt = []

    for time in time_list:
        nearest_time = dataset.find_nearest_time(time)
        time_plt.append(nearest_time)
        print('Will plot {} s in dataset'.format(np.round(nearest_time)))

    if plot==True:
        colors = produce_colour_arr(time_plt)

        with rc_context(fname=rc_fname):
            for i, time in enumerate(time_plt):
                plt.plot(temp.loc[time,:],color=colors[i],label='{} s'.format(round(time)))

            plt.legend(frameon='True')
            plt.xlim(400,900)
            plt.ylim(absMIN,absMAX)
            plt.xlabel('wavelength / nm')
            plt.ylabel('absorbance')
            plt.title('Absorbance Plot \n {}'.format(dataset.exp_label))
            plt.show()
    else:
        pass

    return temp.loc[time_plt,:]

def colourplot(dataset,times,wavs,absorbanceMIN,absorbanceMAX):

    """Function to plot an absorbance colourmap

	Args:
		dataset representing the Dataset object
        times representing a list of times (list)
        wavs representing a list of wavelengths (list)
        absorbanceMIN representing the minimum on the abs scale (float)
        absorbanceMAX representing the maximum on the abs scale (float)

	Returns:
		None

	"""

    # validation to ensure the colour map creation returns no errors
    if absorbanceMIN > 0:
        print('Miniumum absorbance value must be below 0')
    else:
        pass

    # create a colour map that is blue for negative and red for positive
    levels = MaxNLocator(nbins=25).tick_values(absorbanceMIN,absorbanceMAX)
    OD_n=int(np.floor((((absorbanceMAX-0)/(absorbanceMAX-absorbanceMIN))*256)))
    OD_p=256-OD_n
    cm_t=cm.get_cmap('RdBu', OD_p)
    cm_b=cm.get_cmap('RdBu', OD_n)
    newcolors=np.vstack((cm_t(np.linspace(1, 0.5, OD_p)),cm_b(np.linspace(0.5, 0, OD_n))))
    newcmp=ListedColormap(newcolors, name='shifted_cmap')

    with rc_context(fname=rc_fname):
        fig = plt.figure(figsize=(8, 4))
        grid = plt.GridSpec(4, 8, hspace=0.3, wspace=0.3, right=0.75, bottom=0.2)
        main_ax = fig.add_subplot(grid[:-1, 1:])
        y_plot = fig.add_subplot(grid[:-1, 0], sharey=main_ax)
        x_plot = fig.add_subplot(grid[-1, 1:], sharex=main_ax)

    # main axes
    with rc_context(fname=rc_fname):
        cb=main_ax.contourf(dataset.times, dataset.wavelengths, dataset.abs_data.transpose(), levels=levels, cmap=newcmp)
        main_ax.set_ylim(450,900)
        main_ax.set_xlim(0,max(dataset.times))

        for time in times:
            main_ax.axvline(time,color='black',linewidth=1,linestyle=':')

        for wav in wavs:
            main_ax.axhline(wav,color='black',linewidth=1,linestyle=':')

        main_ax.tick_params(labelcolor='none', top=True, bottom=True, left=True, right=True)

    # time trace at specific wavelength
    with rc_context(fname=rc_fname):
        x_plot.plot(wav_abs(dataset,wavs,absorbanceMIN,absorbanceMAX))
        x_plot.set_ylim(absorbanceMIN,absorbanceMAX)
        x_plot.set_xlim(0,max(dataset.times))
        x_plot.set_xlabel('Time / s')

    # wavelength trace at specific time
    with rc_context(fname=rc_fname):
        y_plot.plot((time_abs(dataset,times,absorbanceMIN,absorbanceMAX).transpose()),dataset.wavelengths)
        y_plot.axvline(0,color='lightgrey',linewidth=1,linestyle='--')
        y_plot.set_xlim(absorbanceMIN,absorbanceMAX)
        y_plot.set_ylim(450,900)
        y_plot.set_ylabel('Wavelength / nm')

    # colour bar
    with rc_context(fname=rc_fname):
        fig.add_subplot(111, frameon=False)
        plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
        plt.grid(False)
        cbax = fig.add_axes([0.8, 0.41, 0.03, 0.45])
        cbar=fig.colorbar(cb,cax=cbax, orientation='vertical')
        cbar.set_label('Absorbance', rotation=270,labelpad=15)

    plt.show()

def plot_model(reaction):

    """Function to plot the modelled concentration profile

	Args:
		reaction representing a Reaction object

	Returns:
		None

	"""

    if hasattr(reaction, 'model'):
        with rc_context(fname=rc_fname):
            plt.plot(reaction.model["R"],color=OPT[0],label='Radical Cation')
            plt.legend(frameon='True')

            for time in reaction.turning_points:
                plt.axvline(time,color='grey',linestyle='--',linewidth='1')

            plt.xlim(0,max(reaction.dataset.times)+100)
            plt.ylim(0,max(reaction.model["R"].values)+2)
            plt.xlabel('Time /s')
            plt.ylabel('Concentration / $\mu$M')
            plt.show()

    else:
        print('Calculate model first.')

def plot_rates(reaction):

    """Function to plot k(t) as a function of time

	Args:
		reaction representing a Reaction object

	Returns:
		None

	"""

    if hasattr(reaction, 'model'):

        rates = reaction.rates

        with rc_context(fname=rc_fname):

            for i, column in enumerate(rates):
                plt.plot(rates[column],label='{}'.format(column))

            for time in reaction.turning_points:
                plt.axvline(time,color='grey',linestyle='--',linewidth='1')

            plt.xlim(0,max(reaction.dataset.times)+100)
            plt.ylim(bottom=0)
            plt.xlabel('Time /s')
            plt.ylabel('Rate / $\mu$M$s^{-1}$')
            plt.legend(frameon='True')
            plt.show()

    else:
        print('Calculate model first.')

def compare_model(reaction):

    """Function to plot both the experimental data and model side-by-side

	Args:
		reaction representing a Reaction object

	Returns:
		None

	"""
    # create 1um error bars
    error = np.ones(len(reaction.model.index))

    if hasattr(reaction, 'model'):
        with rc_context(fname=rc_fname):
            mod = reaction.model["R"]
            plt.plot(mod, color='tomato',label='Model')

            if reaction.state=='drift_corrected':
                y = reaction.conc_profile_d
            else:
                y = reaction.conc_profile

            plt.plot(y,color=OPT[0],label='Data')
            plt.fill_between(reaction.model.index, mod.values - error, mod.values + error, color='tomato', alpha=0.3, label='$\pm 1\mu$M interval')

            for time in reaction.turning_points:
                plt.axvline(time,color='grey',linestyle='--',linewidth='1')

            plt.xlim(0,max(reaction.dataset.times)+100)
            plt.ylim(0,max(y.values)+2)
            plt.xlabel('Time /s')
            plt.ylabel('Concentration / $\mu$M')
            plt.legend(frameon='True')
            plt.show()

    else:
        print('Calculate model first.')
