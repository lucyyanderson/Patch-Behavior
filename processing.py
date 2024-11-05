
# For analyzing pyphotometry data. Similar to pyphotometry_preprocessing.py, but plots fewer things 
# and uses the pyphotometry function preprocess_data. 

import os
import numpy as  np
import matplotlib.pyplot as plt
from scipy.signal import medfilt, butter, filtfilt
from scipy.stats import linregress
from scipy.optimize import curve_fit, minimize
import scipy


from data_import import import_ppd, preprocess_data

#data_folder = r'\\research-cifs.nyumc.org\research\buzsakilab\Buzsakilabspace\LabShare\ZutshiI\patchTask\test_random_pulses'
data_folder = '/Users/lucyanderson/Library/Mobile Documents/com~apple~CloudDocs/Buzsaki/pyphotometry/N7/10.29/'
file = 'N7_striatum-2024-10-29-130040' 
data_filename = file + '.ppd'
data = import_ppd(os.path.join(data_folder, data_filename))   #, low_pass=20, high_pass=0.001)
dLight_raw = data['analog_1'] # green
cherry_raw = data['analog_2'] # red
barcode = data['digital_1'] 
time_seconds = data['time']/1000
sampling_rate = data['sampling_rate']

preprocessed_data = {}


processed_signal = preprocess_data(data_dict=data, 
                                   signal="analog_1", 
                                   control="analog_2", 
                                   low_pass=10,
                                   normalisation="dF/F",
                                   plot=True,
                                   fig_path = os.path.join(os.path.expanduser('~'), 'Downloads', 'N7_striatum.png'))  # change name you want to save as

# Save plot to the Downloads folder
#downloads_path = os.path.join(os.path.expanduser('~'), 'Downloads', 'N7_striatum_plot.png')
#plt.gcf().savefig(downloads_path)

# Plot raw signals
fig,ax1=plt.subplots()  # create a plot to allow for dual y-axes plotting
plot1=ax1.plot(time_seconds, dLight_raw, 'g', label='dLight') #plot dLight on left y-axis
ax2=plt.twinx()# create a right y-axis, sharing x-axis on the same plot
plot2=ax2.plot(time_seconds, cherry_raw, 'r', label='mCherry') # plot mCherry on right y-axis

# Add labels for both axes
ax1.set_xlabel('Time (seconds)')
ax1.set_ylabel('dLight', color='g')
ax2.set_ylabel('mCherry', color='r')

plt.show()


# Plot rewards times as ticks.
#reward_ticks = ax1.plot(reward_cue_times, np.full(np.size(reward_cue_times), 1.625), label='Reward Cue', color='w', marker="|", mec='k')


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

preprocessed_data['sampling_rate'] = sampling_rate
preprocessed_data['barcodesOn'] = timestampsOn # formerly timestampsOn
preprocessed_data['barcodesOnOff'] = timestampsOnOff # formerly timestampsOnOff
preprocessed_data['highLow'] = barcode
preprocessed_data['timestamps'] = time_seconds


#%% normalisation to save to matlab

"""DENOISING"""
# Lowpass filter - zero phase filtering (with filtfilt) is used to avoid distorting the signal.
b,a = butter(2, 10, btype='low', fs=sampling_rate)
dLight_denoised = filtfilt(b,a, dLight_raw)
cherry_denoised = filtfilt(b,a, cherry_raw)


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

# Fit curve to mCherry signal.
max_sig = np.max(cherry_denoised)
inital_params = [max_sig/2, max_sig/4, max_sig/4, 3600, 0.1]
bounds = ([0      , 0      , 0      , 600  , 0],
          [max_sig, max_sig, max_sig, 36000, 1])
cherry_parms, parm_cov = curve_fit(double_exponential, time_seconds, cherry_denoised, 
                                  p0=inital_params, bounds=bounds, maxfev=1000)
cherry_expfit = double_exponential(time_seconds, *cherry_parms)

dLight_detrended = dLight_denoised - dLight_expfit
cherry_detrended = cherry_denoised - cherry_expfit


"""MOTION CORRECTION"""
slope, intercept, r_value, p_value, std_err = linregress(x=cherry_detrended, y=dLight_detrended)

#print('Slope    : {:.3f}'.format(slope))
#print('R-squared: {:.3f}'.format(r_value**2))


dLight_est_motion = intercept + slope * cherry_detrended
dLight_corrected = dLight_detrended - dLight_est_motion


"""NORMALIZATION"""
 #METHOD 1: dF/F
dLight_dF_F = 100*dLight_corrected/dLight_expfit

#METHOD 2: Z-SCORE
dLight_zscored = (dLight_corrected-np.mean(dLight_corrected))/np.std(dLight_corrected)

#%% SAVE DATA IN MATLAB FORMAT
preprocessed_data['dlight_z'] = dLight_zscored
preprocessed_data['dlight_df'] = dLight_dF_F
scipy.io.savemat(os.path.join(data_folder, f'{file}_photometry.mat'), {'photometryData': preprocessed_data})

