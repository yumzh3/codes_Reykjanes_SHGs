## conduct spectral analysis of SHGs MgO time series

## written by Mingzhen Yu
## Jun 17, 2024
## last modified: Sep 30, 2025
    
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import correlate
from scipy.signal import periodogram
from scipy.signal.windows import dpss
from scipy.signal import detrend
from statsmodels.tsa.ar_model import AutoReg
from scipy.signal import find_peaks
import matplotlib.gridspec as gridspec
from matplotlib.font_manager import FontProperties
import matplotlib.patheffects as path_effects


# %% spectral analysis
# load MgO time series data for core 22, e.g., data_22
# get the age data as time, and MgO data as signal

signal_nan_index = signal[signal[:].isna()].index
signal = pd.to_numeric(signal,errors='coerce')
signal = signal.interpolate(limit_direction='both')

detrended = detrend(signal)  # detrend improves clarity of cyclic signals
ar_model = AutoReg(detrended,lags=1,old_names=False).fit()  # prewhitening with AR(1) to flatten red-noise bias, 1 means lags=1
ar_coef = ar_model.params[1]
prewhitened = detrended[1:]-ar_coef*detrended[:-1]

signal_pretreated = prewhitened

# === use single-taper spectral analysis
fs = 1/2
frequencies, power = periodogram(signal_pretreated , fs=0.5) 
frequencies = frequencies[1:] # remove zero frequency at the beginning to avoid division by zero
periods = 1/frequencies  
power = power[1:]

np.random.seed(0)
n_simulations = 1000
simulated_max = []
for _ in range(n_simulations):
    shuffled = np.random.permutation(signal_pretreated)
    _, sim_power = periodogram(shuffled, fs=fs)
    simulated_max.append(sim_power.max())
p_value = np.mean(np.array(simulated_max) >= power.max())

sorted_indices_power = np.argsort(power)[::-1]
np.random.seed(0)
target_freq = frequencies[sorted_indices_power[0]]  # target frequency for p-value test which is the frequency with highest power
n_simulations = 1000
simulated_power = []
for _ in range(n_simulations):
    shuffled = np.random.permutation(signal_pretreated)
    f, sim_power = periodogram(shuffled, fs=fs)
    target_idx = np.argmin(np.abs(f-target_freq))
    simulated_power.append(sim_power[target_idx])
p_value_target = np.mean(np.array(simulated_power) >= power.max())

np.random.seed(0)
target_freq_2 = frequencies[sorted_indices_power[1]]  # target frequency for p-value test which is the frequency with highest power
n_simulations = 1000
simulated_power = []
for _ in range(n_simulations):
    shuffled = np.random.permutation(signal_pretreated)
    f, sim_power = periodogram(shuffled, fs=fs)
    target_idx = np.argmin(np.abs(f-target_freq_2))
    simulated_power.append(sim_power[target_idx])
p_value_target_2 = np.mean(np.array(simulated_power) >= power[sorted_indices_power[1]])

np.random.seed(0)
target_freq_3 = frequencies[sorted_indices_power[2]]  # target frequency for p-value test which is the frequency with highest power
n_simulations = 1000
simulated_power = []
for _ in range(n_simulations):
    shuffled = np.random.permutation(signal_pretreated)
    f, sim_power = periodogram(shuffled, fs=fs)
    target_idx = np.argmin(np.abs(f-target_freq_3))
    simulated_power.append(sim_power[target_idx])
p_value_target_3 = np.mean(np.array(simulated_power) >= power[sorted_indices_power[2]])

# === plot
plot_params = {'figsize': (16,13),
               'linecolor_detrend': 'gray',
               'linecolor_prewhitened': 'blue',
               'linecolor_spectrum': 'k',
               'linecolor_power': 'red',
               'linewidth_detrend': 2,
               'linewidth_prewhitened': 2,
               'linewidth_spectrum': 2,
               'linewidth_power': 1.5,
               'linestyle_detrend':'--',
               'linestyle_prewhitened':'-',
               'linestyle_spectrum':'-',
               'linestyle_power': '--',
               'bin_face_color': 'lightgray',
               'bin_edge_color': 'k',
               'fontsize_tickmark': 20,
               'fontsize_label': 22,
               'fontsize_letter': 22,
               'fontname_all': 'Arial',
               'textcolor': 'k',
               'textstyle': 'italic'}
fig = plt.figure(figsize=plot_params['figsize'])
gs = gridspec.GridSpec(2, 2)
# == plot detrend data and prewhitened data
ax0 = fig.add_subplot(gs[0])
ax0.plot(time, detrended,color=plot_params['linecolor_detrend'],
         linewidth=plot_params['linewidth_detrend'],linestyle=plot_params['linestyle_detrend'],
         label='detrended data')
ax0.plot(time[1:], prewhitened,color=plot_params['linecolor_prewhitened'],
         linewidth=plot_params['linewidth_prewhitened'],linestyle=plot_params['linestyle_prewhitened'],
         label='prewhitened data')
ax0.tick_params(labelsize=plot_params['fontsize_tickmark'])  
for label in ax0.get_xticklabels() + ax0.get_yticklabels():
    label.set_fontname(plot_params['fontname_all'])
ax0.set_xlim(0,140)
ax0.set_ylim(-0.45,0.35)
ax0.set_xlabel("Age (ka)", fontsize=plot_params['fontsize_label'], font=plot_params['fontname_all'])
ax0.set_ylabel("Pretreated MgO", fontsize=plot_params['fontsize_label'], font=plot_params['fontname_all'])
leg = ax0.legend(loc='upper left',prop=FontProperties(family=plot_params['fontname_all'],size=plot_params['fontsize_tickmark']),
           facecolor='white',frameon=True,markerscale=1,
           handletextpad=0.05, borderaxespad=0.1, labelspacing=0.2)
leg.get_frame().set_edgecolor('black')
fig.text(0.016, 0.97, 'A', fontsize=plot_params['fontsize_letter'], fontweight='bold', fontname=plot_params['fontname_all'])
# == plot power spectrum with top peaks
ax1 = fig.add_subplot(gs[1])
ax1.plot(frequencies, power, color=plot_params['linecolor_spectrum'],
         linewidth=plot_params['linewidth_spectrum'],linestyle=plot_params['linestyle_spectrum'],
         label='power spectrum')
for i in range(3):
    ax1.axvline(frequencies[sorted_indices_power[i]], color='red', linestyle='--')
ax1.text(frequencies[sorted_indices_power[0]]-0.013, -0.004,
         str(round(periods[sorted_indices_power[0]])),
         color=plot_params['textcolor'], ha='left', va='bottom', 
         fontsize=plot_params['fontsize_label'],font=plot_params['fontname_all'],fontstyle=plot_params['textstyle'])
ax1.text(frequencies[sorted_indices_power[1]], -0.004,
         str(round(periods[sorted_indices_power[1]])),
         color=plot_params['textcolor'], ha='left', va='bottom', 
         fontsize=plot_params['fontsize_label'],font=plot_params['fontname_all'],fontstyle=plot_params['textstyle'])
ax1.text(frequencies[sorted_indices_power[2]]-0.013, -0.004,
         str(round(periods[sorted_indices_power[2]])),
         color=plot_params['textcolor'], ha='left', va='bottom', 
         fontsize=plot_params['fontsize_label'],font=plot_params['fontname_all'],fontstyle=plot_params['textstyle'])
ax1.tick_params(labelsize=plot_params['fontsize_tickmark'])  
for label in ax1.get_xticklabels() + ax1.get_yticklabels():
    label.set_fontname(plot_params['fontname_all'])
ax1.set_xlabel('Frequency (1/kyr)', fontsize=plot_params['fontsize_label'], font=plot_params['fontname_all'])
ax1.set_ylabel('Spectral Power', fontsize=plot_params['fontsize_label'], font=plot_params['fontname_all'])
leg = ax1.legend(loc='upper right',prop=FontProperties(family=plot_params['fontname_all'],size=plot_params['fontsize_tickmark']),
           facecolor='white',frameon=True,markerscale=1,
           handletextpad=0.05, borderaxespad=0.1, labelspacing=0.2)
leg.get_frame().set_edgecolor('black')
fig.text(0.51, 0.97, 'B', fontsize=plot_params['fontsize_letter'], fontweight='bold', fontname=plot_params['fontname_all'])
# === Monte Carlo shuffle test histogram for global
ax2 = fig.add_subplot(gs[2])
ax2.hist(simulated_max, bins=30, color=plot_params['bin_face_color'], edgecolor=plot_params['bin_edge_color'])
ax2.axvline(power.max(), color=plot_params['linecolor_power'], linestyle=plot_params['linestyle_power'],
            linewidth=plot_params['linewidth_power'],
            label=f'observed max power\n$p$-value = {np.mean(np.array(simulated_max) >= power.max()):.3f}')
ax2.tick_params(labelsize=plot_params['fontsize_tickmark'])  
for label in ax2.get_xticklabels() + ax2.get_yticklabels():
    label.set_fontname(plot_params['fontname_all'])
ax2.set_xlabel('Maximum Power in Shuffled Spectra', fontsize=plot_params['fontsize_label'], font=plot_params['fontname_all'])
ax2.set_ylabel('Frequency', fontsize=plot_params['fontsize_label'], font=plot_params['fontname_all'])
leg = ax2.legend(loc='upper right',prop=FontProperties(family=plot_params['fontname_all'],size=plot_params['fontsize_tickmark']),
           facecolor='white',frameon=True,markerscale=1,
           handletextpad=0.05, borderaxespad=0.1, labelspacing=0.2)
leg.get_frame().set_edgecolor('black')
fig.text(0.016, 0.47, 'C', fontsize=plot_params['fontsize_letter'], fontweight='bold', fontname=plot_params['fontname_all'])
# === Bar Plot of Top Peak p-values
ax3 = fig.add_subplot(gs[3])
bars = ax3.bar([f'{round(periods[sorted_indices_power[i]])} kyr' for i in range(3)],
               [-np.log10(p) for p in [p_value_target,p_value_target_2,p_value_target_3]],
               color=['darkgreen' if p < 0.05 else 'gray' for p in [p_value_target,p_value_target_2,p_value_target_3]])
ax3.axhline(-np.log10(0.05), color='red', linestyle='--', linewidth=1.5,label='$p$-value = 0.05')
ax3.tick_params(labelsize=plot_params['fontsize_tickmark'])  
for label in ax3.get_xticklabels() + ax3.get_yticklabels():
    label.set_fontname(plot_params['fontname_all'])
ax3.set_ylabel('-log10($p$-value)', fontsize=plot_params['fontsize_label'], font=plot_params['fontname_all'])
leg = ax3.legend(loc='upper right',prop=FontProperties(family=plot_params['fontname_all'],size=plot_params['fontsize_tickmark']),
           facecolor='none',frameon=True,markerscale=1,
           handletextpad=0.05, borderaxespad=0.1, labelspacing=0.2)
leg.get_frame().set_edgecolor('black')
fig.text(0.51, 0.47, 'D', fontsize=plot_params['fontsize_letter'], fontweight='bold', fontname=plot_params['fontname_all'])
plt.tight_layout()



