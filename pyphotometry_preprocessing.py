# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 17:42:57 2024

@author: Manu
"""

import os
import numpy as  np
import pylab as plt
from scipy.signal import medfilt, butter, filtfilt
from scipy.stats import linregress
from scipy.optimize import curve_fit, minimize
from data_import import import_ppd
from barcodes import extract_barcodes_from_times
import scipy

#set default plot properties
plt.rcParams['figure.figsize'] = [14, 12] # Make default figure size larger.
plt.rcParams['axes.xmargin'] = 0          # Make default margin on x axis zero.
plt.rcParams['axes.labelsize'] = 12     #Set default axes label size 
plt.rcParams['axes.titlesize']=15
plt.rcParams['axes.titleweight']='heavy'
plt.rcParams['ytick.labelsize']= 10
plt.rcParams['xtick.labelsize']= 10
plt.rcParams['legend.fontsize']=12
plt.rcParams['legend.markerscale']=2

preprocessed_data = {}


data_folder = r'\\research-cifs.nyumc.org\research\buzsakilab\Buzsakilabspace\LabShare\ZutshiI\patchTask\test_random_pulses'
file = 'N7_striatum-2024-10-22-172942' 
data_filename = file + '.ppd'
data = import_ppd(os.path.join(data_folder, data_filename))
dLight_raw = data['analog_1'] #green
TdTom_raw = data['analog_2'] #red
barcode = data['digital_1']
time_seconds = data['time']/1000
sampling_rate = data['sampling_rate']

preprocessed_data['sampling_rate'] = sampling_rate
#%%
""" GET BARCODE TIMESTAMPS """


timestampsOn = []
timestampsOff = []
for i,signal in enumerate(barcode):
    if signal == 1 and barcode[i-1]==0:
        timestamp_On = i
        timestampsOn.append(i)
    if signal == 0 and barcode[i-1]==1:
        timestamp_Off = i
        timestampsOff.append(i)

timestampsOnOff = []
for On,Off in zip(timestampsOn,timestampsOff):
    timestampsOnOff.append([On, Off])

timestampsOn = np.array(timestampsOn)
timestampsOff = np.array(timestampsOff)

preprocessed_data['timestampsOn'] = timestampsOn
preprocessed_data['timestampsOnOff'] = timestampsOnOff
preprocessed_data['highLow'] = barcode

#%%
"""RAW SIGNALS"""
# Plot signals

fig,ax1=plt.subplots()  # create a plot to allow for dual y-axes plotting
plot1=ax1.plot(time_seconds, dLight_raw, 'g', label='dLight') #plot dLight on left y-axis
ax2=plt.twinx()# create a right y-axis, sharing x-axis on the same plot
plot2=ax2.plot(time_seconds, TdTom_raw, 'r', label='TdTomato') # plot TdTomato on right y-axis

# Plot rewards times as ticks.
#reward_ticks = ax1.plot(reward_cue_times, np.full(np.size(reward_cue_times), 1.625), label='Reward Cue', color='w', marker="|", mec='k')


#ax1.set_ylim(1.25, 1.65)
#ax2.set_ylim(1.35, 1.75)
ax1.set_xlabel('Time (seconds)')
ax1.set_ylabel('dLight Signal (V)', color='g')
ax2.set_ylabel('TdTomato Signal (V)', color='r')
ax1.set_title('Raw signals')

lines = plot1 + plot2  #line handle for legend
labels = [l.get_label() for l in lines]  #get legend labels
legend = ax1.legend(lines, labels, loc='upper right', bbox_to_anchor=(0.98, 0.93)) #add legend

#%%
"""DENOISING"""
# Lowpass filter - zero phase filtering (with filtfilt) is used to avoid distorting the signal.
b,a = butter(2, 10, btype='low', fs=sampling_rate)
dLight_denoised = filtfilt(b,a, dLight_raw)
TdTom_denoised = filtfilt(b,a, TdTom_raw)

fig,ax1=plt.subplots()
plot1=ax1.plot(time_seconds, dLight_denoised, 'g', label='dLight denoised')
ax2=plt.twinx()
plot2=ax2.plot(time_seconds, TdTom_denoised, 'r', label='TdTomato denoised')
#reward_ticks = ax1.plot(reward_cue_times, np.full(np.size(reward_cue_times), 1.625), label='Reward Cue', color='w', marker="|", mec='k', ms=10)

#ax1.set_ylim(1.25, 1.65)
#ax2.set_ylim(1.35, 1.75)
ax1.set_xlabel('Time (seconds)')
ax1.set_ylabel('dLight Signal (V)', color='g')
ax2.set_ylabel('TdTomato Signal (V)', color='r')
ax1.set_title('Denoised signals')

lines = plot1+plot2 #line handle for legend
labels = [l.get_label() for l in lines]  #get legend labels
legend = ax1.legend(lines, labels, loc='upper right', bbox_to_anchor=(0.98, 0.92)) #add legend

fig,ax1=plt.subplots()  
plot1=ax1.plot(time_seconds, dLight_raw, color='g', alpha=0.3, label='dLight raw')
ax2=plt.twinx()
plot2=ax2.plot(time_seconds, TdTom_raw, color='r', alpha=0.3, label='TdTomato raw') 
plot3=ax1.plot(time_seconds, dLight_denoised, color='g', label='dLight denoised') 
plot4=ax2.plot(time_seconds, TdTom_denoised, color='r', label='TdTomato denoised') 
#reward_ticks = ax1.plot(reward_cue_times, np.full(np.size(reward_cue_times), 1.59), label='Reward Cue',color='w', marker="v", mfc='k', mec='k', ms=8)

ax1.set_xlabel('Time (seconds)')
ax1.set_ylabel('dLight Signal (V)', color='g')
ax2.set_ylabel('TdTomato Signal (V)', color='r')
ax1.set_title('Denoised signals')

lines = plot1+plot2 + plot3 + plot4 
labels = [l.get_label() for l in lines]
legend = ax1.legend(lines, labels, loc='upper right', bbox_to_anchor=(0.93, 0.99))
ax1.set_xlim(1000, 1060) # 60 sec window
#ax1.set_ylim(1.4, 1.625)
#ax2.set_ylim(1.4, 1.625);

#%%
"""PHOTOBLEACHING CORRECTION"""

#METHOD 1
# The double exponential curve we are going to fit.
def double_exponential(t, const, amp_fast, amp_slow, tau_slow, tau_multiplier):
    '''Compute a double exponential function with constant offset.
    Parameters:
    t       : Time vector in seconds.
    const   : Amplitude of the constant offset. 
    amp_fast: Amplitude of the fast component.  
    amp_slow: Amplitude of the slow component.  
    tau_slow: Time constant of slow component in seconds.
    tau_multiplier: Time constant of fast component relative to slow. 
    '''
    tau_fast = tau_slow*tau_multiplier
    return const+amp_slow*np.exp(-t/tau_slow)+amp_fast*np.exp(-t/tau_fast)

# Fit curve to dLight signal.
max_sig = np.max(dLight_denoised)
inital_params = [max_sig/2, max_sig/4, max_sig/4, 3600, 0.1]
bounds = ([0      , 0      , 0      , 600  , 0],
          [max_sig, max_sig, max_sig, 36000, 1])
dLight_parms, parm_cov = curve_fit(double_exponential, time_seconds, dLight_denoised, 
                                  p0=inital_params, bounds=bounds, maxfev=1000)
dLight_expfit = double_exponential(time_seconds, *dLight_parms)

# Fit curve to TdTomato signal.
max_sig = np.max(TdTom_denoised)
inital_params = [max_sig/2, max_sig/4, max_sig/4, 3600, 0.1]
bounds = ([0      , 0      , 0      , 600  , 0],
          [max_sig, max_sig, max_sig, 36000, 1])
TdTom_parms, parm_cov = curve_fit(double_exponential, time_seconds, TdTom_denoised, 
                                  p0=inital_params, bounds=bounds, maxfev=1000)
TdTom_expfit = double_exponential(time_seconds, *TdTom_parms)

#plot fits over denoised data
fig,ax1=plt.subplots()  
plot1=ax1.plot(time_seconds, dLight_denoised, 'g', label='dLight')
plot3=ax1.plot(time_seconds, dLight_expfit, color='k', linewidth=1.5, label='Exponential fit') 
ax2=plt.twinx()
plot2=ax2.plot(time_seconds, TdTom_denoised, color='r', label='TdTomato') 
plot4=ax2.plot(time_seconds, TdTom_expfit,color='k', linewidth=1.5) 


dLight_detrended = dLight_denoised - dLight_expfit
TdTom_detrended = TdTom_denoised - TdTom_expfit

fig,ax1=plt.subplots()  
plot1=ax1.plot(time_seconds, dLight_detrended, 'g', label='dLight')
ax2=plt.twinx()
plot2=ax2.plot(time_seconds, TdTom_detrended, color='r', label='TdTomato') 

ax1.set_xlabel('Time (seconds)')
ax1.set_ylabel('dLight Signal (V)', color='g')
ax2.set_ylabel('TdTomato Signal (V)', color='r')
ax1.set_title('Bleaching Correction by Double Exponential Fit')

lines = plot1+plot2 
labels = [l.get_label() for l in lines]  
legend = ax1.legend(lines, labels, loc='upper right'); 
ax1.set_ylim(-0.18, 0.12)
ax2.set_ylim(-0.1, 0.2);


ax1.set_xlabel('Time (seconds)')
ax1.set_ylabel('dLight Signal (V)', color='g')
ax2.set_ylabel('TdTomato Signal (V)', color='r')
ax1.set_title('Denoised signals with double exponential fits')

lines = plot1 + plot2 + plot3
labels = [l.get_label() for l in lines]  
legend = ax1.legend(lines, labels, loc='upper right'); 
#ax1.set_ylim(1.27, 1.62)
#ax2.set_ylim(1.35, 1.7);

# #METHOD 2
# b,a = butter(2, 0.001, btype='high', fs=sampling_rate)
# dLight_highpass = filtfilt(b,a, dLight_denoised, padtype='even')
# TdTom_highpass = filtfilt(b,a, TdTom_denoised, padtype='even')


# fig,ax1=plt.subplots()  
# plot1=ax1.plot(time_seconds, dLight_highpass, 'g', label='dLight')
# ax2=plt.twinx()
# plot2=ax2.plot(time_seconds, TdTom_highpass, color='r', label='TdTomato') 

# ax1.set_xlabel('Time (seconds)')
# ax1.set_ylabel('dLight Signal (V)', color='g')
# ax2.set_ylabel('TdTomato Signal (V)', color='r')
# ax1.set_title('Bleaching Correction by Highpass Filtering')

# lines = plot1+plot2 
# labels = [l.get_label() for l in lines]  
# legend = ax1.legend(lines, labels, loc='upper right'); 
# ax1.set_ylim(-0.18, 0.12)
# ax2.set_ylim(-0.1, 0.2);

#%%
"""MOTION CORRECTION"""
slope, intercept, r_value, p_value, std_err = linregress(x=TdTom_detrended, y=dLight_detrended)

plt.scatter(TdTom_detrended[::5], dLight_detrended[::5],alpha=0.1, marker='.')
x = np.array(plt.xlim())
plt.plot(x, intercept+slope*x)
plt.xlabel('TdTomato')
plt.ylabel('dLight')
plt.title('TdTomato - dLight correlation.')

print('Slope    : {:.3f}'.format(slope))
print('R-squared: {:.3f}'.format(r_value**2))


dLight_est_motion = intercept + slope * TdTom_detrended
dLight_corrected = dLight_detrended - dLight_est_motion

fig,ax1=plt.subplots()  
plot1=ax1.plot(time_seconds, dLight_detrended, 'b' , label='dLight - pre motion correction', alpha=0.5)
plot3=ax1.plot(time_seconds, dLight_corrected, 'g', label='dLight - motion corrected', alpha=0.5)
plot4=ax1.plot(time_seconds, dLight_est_motion - 0.05, 'y', label='estimated motion')
#reward_ticks = ax1.plot(reward_cue_times, np.full(np.size(reward_cue_times), 0.08), label='Reward Cue',color='w', marker="v", mfc='k', mec='k', ms=8)

ax1.set_xlabel('Time (seconds)')
ax1.set_ylabel('dLight Signal (V)', color='g')
ax1.set_title('Motion Correction')

lines = plot1+plot3+plot4
labels = [l.get_label() for l in lines]  
legend = ax1.legend(lines, labels, loc='upper right', bbox_to_anchor=(0.95, 0.98))

#ax1.set_xlim(1000, 1060)  # 60 sec window
#ax1.set_ylim(-0.075, 0.1);

#%%
"""NORMALIZATION"""
# #METHOD 1: dF/F
# dLight_dF_F = 100*dLight_corrected/dLight_expfit

# fig,ax1=plt.subplots()  
# plot1=ax1.plot(time_seconds, dLight_dF_F, 'g', label='dLight dF/F')
# #reward_ticks = ax1.plot(reward_cue_times, np.full(np.size(reward_cue_times), 6), label='Reward Cue',color='w', marker="v", mfc='k', mec='k', ms=8)

# ax1.set_xlabel('Time (seconds)')
# ax1.set_ylabel('dLight dF/F (%)')
# ax1.set_title('dLight dF/F')

# lines = plot1
# labels = [l.get_label() for l in lines]  
# legend = ax1.legend(lines, labels, loc='upper right', bbox_to_anchor=(0.95, 0.98))

# ax1.set_xlim(1000, 1060)
# ax1.set_ylim(-3, 7);

#METHOD 2: Z-SCORE
dLight_zscored = (dLight_corrected-np.mean(dLight_corrected))/np.std(dLight_corrected)


fig,ax1=plt.subplots()  
plot1=ax1.plot(time_seconds, dLight_zscored, 'g', label='dLight z-score')
#reward_ticks = ax1.plot(reward_cue_times, np.full(np.size(reward_cue_times), 6), label='Reward Cue',color='w', marker="v", mfc='k', mec='k', ms=8)

ax1.set_xlabel('Time (seconds)')
ax1.set_ylabel('dLight z-score')
ax1.set_title('dLight z-scored')

lines = plot1
labels = [l.get_label() for l in lines]  
legend = ax1.legend(lines, labels, loc='upper right', bbox_to_anchor=(0.95, 0.98))

# ax1.set_xlim(1000, 1060)
#ax1.set_ylim(-3, 7);
preprocessed_data['dlight'] = dLight_zscored

#%% SAVE DATA IN MATLAB FORMAT
scipy.io.savemat(os.path.join(data_folder, f'{file}_photometry.mat'), {'photometryData': preprocessed_data})
