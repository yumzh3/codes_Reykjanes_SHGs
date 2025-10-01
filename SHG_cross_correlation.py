## conduct cross-correlation between sediment-hosted glasses (SHGs) geochemical time seris and d18O stack

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


# %% cross-correlation
## load time series data for both SHGs and d18O.
""" 
load moving-average time series data of SHGs from core 22 provided in Supplementary DataS3.
load global benthic d18O stack from Lisiecki and Raymo, 2005 (LR04).
make sure ages in two data sets are the same and d18O value is reversed.
"""
data = pd.read_excel('the SHG data file')  # load SHG data for elements
data_Mg8 = pd.read_excel('the SHG data file')  # load SHG data for elements that are normalized to 8 wt% MgO for differentiation correction
data = data[:-6]  # truncated the core 22 SHGs data to 115 ka
data_Mg8 = data_Mg8[:-6]  # truncated the core 22 SHGs data to 115 ka
data_LR04 = pd.read_excel('LR04 d18O stack data file')  # load LR04 d18O data. The data can be downloaded from their website.
data_LR04_sub = data_LR04[data_LR04['Time (ka)'].isin(data['movingAve_Age'])]  # get the data from LR04 that has the same age of the core data
data_LR04_sub['Benthic d18O (per mil)'] = -data_LR04_sub['Benthic d18O (per mil)']  # convert d18O of LR04 to reverse order

## components for cross-correlation
isotope_for_lag_trace = ['La','K','Rb','Ba','Sm','Th','La_Sm']  
isotope_for_lag_trace8 = ['La8','K8','Rb8','Ba8','Th8']

## parameter for bootstrap and montecarlo
n_bootstrap = 1000
n_montecarlo = 1000
random_seed = 42  # with this seed, results can be replicated if using the same data sets.

# === FUNCTION: Bootstrap Confidence Interval for Lag
"""
Bootstrap d18O and SHG time series scaled by their own standard errors and conduct the cross-correlation to get the best lag value.
Repeat 1000 times, and calculate the mean and 95% confidence intervals, and that is the bootstrap mean lag shown in supplement figures.
Standard error of LR04 d18O is provided in their data set downloaded from their website.
Standard error of core 22 SHGs is provided with the moving-average time seris in Supplementary DataS3.
"""
def bootstrap_lag_confidence(d18O, element_data, d18O_sem, element_sem, n_iter, seed,element_nan_index):
    np.random.seed(seed)
    best_lags = []
    n = len(d18O)   
    for _ in range(n_iter):
        d18O_synthetic = np.random.normal(loc=d18O, scale=d18O_sem)
        element_synthetic = np.random.normal(loc=element_data, scale=element_sem)
        element_synthetic[element_nan_index] = np.nan
        element_synthetic = pd.to_numeric(element_synthetic,errors='coerce')
        element_synthetic = pd.Series(element_synthetic)
        element_synthetic = element_synthetic.interpolate(limit_direction='both')
        d18O_norm = (d18O_synthetic - np.mean(d18O_synthetic)) / np.std(d18O_synthetic,ddof=1)
        elem_norm = (element_synthetic - np.mean(element_synthetic)) / np.std(element_synthetic,ddof=1)
        corr = correlate(d18O_norm, elem_norm, mode='full')
        lags = np.arange(-n + 1, n)
        best_lags.append(lags[np.argmax(corr)])
    ci_lower = np.percentile(best_lags, 2.5)
    ci_upper = np.percentile(best_lags, 97.5)
    return np.mean(best_lags), (ci_lower, ci_upper)
# === FUNCTION: Monte Carlo Shuffle Test for Correlation
"""
Suffule the SHG time series data and conduct cross-correlation with the fixed d18O record for 1000 times.
Check the statistical significance of the observed correlation from the original unshuffled data set.
"""
def monte_carlo_significance(d18O, element_data, n_iter, seed):
    np.random.seed(seed)
    d18O_norm = (d18O - np.mean(d18O)) / np.std(d18O,ddof=1)
    elem_norm = (element_data - np.mean(element_data)) / np.std(element_data,ddof=1)
    observed_corr = correlate(d18O_norm, elem_norm, mode='full')
    observed_max = np.max(np.abs(observed_corr))
    null_max_corrs = []
    for _ in range(n_iter):
        shuffled_element_data= np.random.permutation(element_data)
        shuffled_elem_norm = (shuffled_element_data - np.mean(shuffled_element_data)) / np.std(shuffled_element_data,ddof=1)
        null_corr = correlate(d18O_norm, shuffled_elem_norm, mode='full')
        null_max_corrs.append(np.max(np.abs(null_corr)))
    p_value = np.mean(np.array(null_max_corrs) >= observed_max)
    return observed_max, p_value

corr_observed_trace = {}
results_trace = []
d18O = data_LR04_sub['Benthic d18O (per mil)']
d18O = d18O.values
d18O_sem = data_LR04_sub['Standard error (per mil)']
d18O_sem = d18O_sem.values
for element in isotope_for_lag_trace:
    # == use moving-average time series
    element_data = data['movingAve_'+element]
    element_data = element_data.values
    element_sem = data['SEM_'+element]
    element_nan_index = element_sem[element_sem[:].isna()].index
    element_sem = np.nan_to_num(element_sem,nan=0)
    
    d18O_norm = (d18O - np.mean(d18O)) / np.std(d18O,ddof=1)
    elem_norm = (element_data - np.mean(element_data)) / np.std(element_data,ddof=1)
    corr = correlate(d18O_norm, elem_norm, mode='full')
    lags = np.arange(-len(d18O_norm) + 1, len(d18O_norm))
    best_lags_observed = lags[np.argmax(corr)]
    corr_observed_trace[element] = corr
    mean_lag,(ci_low,ci_high) = bootstrap_lag_confidence(d18O, element_data, d18O_sem, element_sem, n_bootstrap, random_seed,element_nan_index)
    max_corr,p_val = monte_carlo_significance(d18O, element_data, n_montecarlo, random_seed)
    results_trace.append({
        'Element': element,
        'Best_lag_observed': best_lags_observed,
        'Best_lag': mean_lag,
        'Lag_95%CI_low': ci_low,
        'Lag_95%CI_high': ci_high,
        'Max_correlation': max_corr,
        'Correlation_p-value': p_val})
df_results_trace = pd.DataFrame(results_trace)  # an output data frame for results.

corr_observed_trace8 = {}
results_trace8 = []
d18O = data_LR04_sub['Benthic d18O (per mil)']
d18O = d18O.values
d18O_sem = data_LR04_sub['Standard error (per mil)']
d18O_sem = d18O_sem.values
for element in isotope_for_lag_trace8:
    # == use moving average
    element_data = data_Mg8['movingAve_'+element]
    element_data = element_data.values
    element_sem = data_Mg8['SEM_'+element]
    element_sem = np.nan_to_num(element_sem,nan=0)
      
    d18O_norm = (d18O - np.mean(d18O)) / np.std(d18O,ddof=1)
    elem_norm = (element_data - np.mean(element_data)) / np.std(element_data,ddof=1)
    corr = correlate(d18O_norm, elem_norm, mode='full')
    lags = np.arange(-len(d18O_norm) + 1, len(d18O_norm))
    best_lags_observed = lags[np.argmax(corr)]
    corr_observed_trace8[element] = corr
    mean_lag,(ci_low,ci_high) = bootstrap_lag_confidence(d18O, element_data, d18O_sem, element_sem, n_bootstrap, random_seed, element_nan_index)
    max_corr,p_val = monte_carlo_significance(d18O, element_data, n_montecarlo, random_seed)
    results_trace8.append({
        'Element': element,
        'Best_lag_observed': best_lags_observed,
        'Best_lag': mean_lag,
        'Lag_95%CI_low': ci_low,
        'Lag_95%CI_high': ci_high,
        'Max_correlation': max_corr,
        'Correlation_p-value': p_val})
df_results_trace8 = pd.DataFrame(results_trace8)  # an output data frame for results.

# === plot
plot_params = {'figsize': (18,26),
               'marker_corr': 'o',
               'marker_fc_corr':'k',
               'linecolor_corr': 'k',
               'linecolor_observed_lag': 'red',
               'linecolor_bstp_lag': 'dimgray',
               'linewidth_corr': 1.5,
               'linewidth_observed_lag': 2,
               'linewidth_bstp_lag': 1.5,
               'linestyle_corr':'-',
               'linestyle_observed_lag':'--',
               'linestyle_bstp_lag':'--',
               'shade_color_CI': 'gray',
               'fontsize_tickmark': 20,
               'fontsize_label': 20,
               'fontsize_letter': 22,
               'fontname_all': 'Arial',
               'textcolor': 'k',
               'textstyle': 'italic'}
fig,axes = plt.subplots(nrows=4,ncols=2,sharex='col',   
                       figsize=plot_params['figsize'],
                       gridspec_kw={'hspace':0}) 

ax1 = axes[0][0]
ax1.plot(lags,corr_observed_trace['La'],color=plot_params['linecolor_corr'],
          marker=plot_params['marker_corr'],markerfacecolor=plot_params['marker_fc_corr'],
          linewidth=plot_params['linewidth_corr'],linestyle=plot_params['linestyle_corr'],zorder=0)
ax1.axvline(df_results_trace.loc[df_results_trace['Element']=='La','Best_lag'].values[0],
                    color=plot_params['linecolor_bstp_lag'],linewidth=plot_params['linewidth_bstp_lag'],
                    linestyle=plot_params['linestyle_bstp_lag'],label=f"bootstrap mean lag = {df_results_trace.loc[df_results_trace['Element']=='La','Best_lag'].values[0]:.1f}",zorder=1)
ax1.axvspan(df_results_trace.loc[df_results_trace['Element']=='La','Lag_95%CI_low'].values[0],
                    df_results_trace.loc[df_results_trace['Element']=='La','Lag_95%CI_high'].values[0],
                    color=plot_params['shade_color_CI'],alpha=0.3,zorder=1.5)
ax1.axvline(df_results_trace.loc[df_results_trace['Element']=='La','Best_lag_observed'].values[0],
                    color=plot_params['linecolor_observed_lag'],linewidth=plot_params['linewidth_observed_lag'],
                    linestyle=plot_params['linestyle_observed_lag'],label=f"best lag observed = {df_results_trace.loc[df_results_trace['Element']=='La','Best_lag_observed'].values[0]}",zorder=2)
ax1.set_ylim(-30,65)
ax1.tick_params(labelsize=plot_params['fontsize_tickmark'])  
for label in ax1.get_xticklabels() + ax1.get_yticklabels():
    label.set_fontname(plot_params['fontname_all'])
ax1.set_ylabel('Correlation',fontsize=plot_params['fontsize_label'],font=plot_params['fontname_all'])
leg = ax1.legend(loc='upper left',
            prop={'family':plot_params['fontname_all'],'size':plot_params['fontsize_label']},
            facecolor='white',frameon=True,markerscale=0.3,
            handletextpad=0.1, borderaxespad=0.05, labelspacing=0.15, columnspacing=0.5)
leg.get_frame().set_edgecolor('black')
ax1.text(52,57,'La',fontsize=plot_params['fontsize_label'],fontname=plot_params['fontname_all'])
fig.text(0.08, 0.875, 'A', fontsize=plot_params['fontsize_letter'], fontweight='bold', fontname=plot_params['fontname_all'])

ax2 = axes[0][1]
ax2.plot(lags,corr_observed_trace['La_Sm'],color=plot_params['linecolor_corr'],
          marker=plot_params['marker_corr'],markerfacecolor=plot_params['marker_fc_corr'],
          linewidth=plot_params['linewidth_corr'],linestyle=plot_params['linestyle_corr'],zorder=0)
ax2.axvline(df_results_trace.loc[df_results_trace['Element']=='La_Sm','Best_lag'].values[0],
                    color=plot_params['linecolor_bstp_lag'],linewidth=plot_params['linewidth_bstp_lag'],
                    linestyle=plot_params['linestyle_bstp_lag'],zorder=1)
ax2.axvspan(df_results_trace.loc[df_results_trace['Element']=='La_Sm','Lag_95%CI_low'].values[0],
                    df_results_trace.loc[df_results_trace['Element']=='La_Sm','Lag_95%CI_high'].values[0],
                    color=plot_params['shade_color_CI'],alpha=0.3,zorder=1.5)
ax2.axvline(df_results_trace.loc[df_results_trace['Element']=='La_Sm','Best_lag_observed'].values[0],
                    color=plot_params['linecolor_observed_lag'],linewidth=plot_params['linewidth_observed_lag'],
                    linestyle=plot_params['linestyle_observed_lag'],zorder=2)
ax2.set_ylim(-38,58)
ax2.tick_params(labelsize=plot_params['fontsize_tickmark'])  
for label in ax2.get_xticklabels() + ax2.get_yticklabels():
    label.set_fontname(plot_params['fontname_all'])
ax2.set_ylabel('Correlation',fontsize=plot_params['fontsize_label'],font=plot_params['fontname_all'])
ax2.text(41,50,"La/Sm",
                fontsize=plot_params['fontsize_label'],fontname=plot_params['fontname_all'])
txt = ax2.text(-59,46,f"bootstrap mean lag = {df_results_trace.loc[df_results_trace['Element']=='La_Sm','Best_lag'].values[0]:.1f}\nbest lag observed = {df_results_trace.loc[df_results_trace['Element']=='La_Sm','Best_lag_observed'].values[0]}",
                fontsize=plot_params['fontsize_label'],fontname=plot_params['fontname_all'],
                bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.2'))
fig.text(0.51, 0.875, 'B', fontsize=plot_params['fontsize_letter'], fontweight='bold', fontname=plot_params['fontname_all'])

ax3 = axes[1][0]
ax3.plot(lags,corr_observed_trace['K'],color=plot_params['linecolor_corr'],
          marker=plot_params['marker_corr'],markerfacecolor=plot_params['marker_fc_corr'],
          linewidth=plot_params['linewidth_corr'],linestyle=plot_params['linestyle_corr'],zorder=0)
ax3.axvline(df_results_trace.loc[df_results_trace['Element']=='K','Best_lag'].values[0],
                    color=plot_params['linecolor_bstp_lag'],linewidth=plot_params['linewidth_bstp_lag'],
                    linestyle=plot_params['linestyle_bstp_lag'],zorder=1)
ax3.axvspan(df_results_trace.loc[df_results_trace['Element']=='K','Lag_95%CI_low'].values[0],
                    df_results_trace.loc[df_results_trace['Element']=='K','Lag_95%CI_high'].values[0],
                    color=plot_params['shade_color_CI'],alpha=0.3,zorder=1.5)
ax3.axvline(df_results_trace.loc[df_results_trace['Element']=='K','Best_lag_observed'].values[0],
                    color=plot_params['linecolor_observed_lag'],linewidth=plot_params['linewidth_observed_lag'],
                    linestyle=plot_params['linestyle_observed_lag'],zorder=2)
ax3.set_ylim(-30,55)
ax3.tick_params(labelsize=plot_params['fontsize_tickmark'])  
for label in ax3.get_xticklabels() + ax3.get_yticklabels():
    label.set_fontname(plot_params['fontname_all'])
ax3.set_ylabel('Correlation',fontsize=plot_params['fontsize_label'],font=plot_params['fontname_all'])
ax3.text(48,45,"K$_2$O",
                fontsize=plot_params['fontsize_label'],fontname=plot_params['fontname_all'])
ax3.text(-59,45,f"bootstrap mean lag = {df_results_trace.loc[df_results_trace['Element']=='K','Best_lag'].values[0]:.1f}\nbest lag observed = {df_results_trace.loc[df_results_trace['Element']=='K','Best_lag_observed'].values[0]}",
                fontsize=plot_params['fontsize_label'],fontname=plot_params['fontname_all'],
                bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.2'))
fig.text(0.08, 0.683, 'C', fontsize=plot_params['fontsize_letter'], fontweight='bold', fontname=plot_params['fontname_all'])

ax4 = axes[1][1]
ax4.plot(lags,corr_observed_trace['Th'],color=plot_params['linecolor_corr'],
          marker=plot_params['marker_corr'],markerfacecolor=plot_params['marker_fc_corr'],
          linewidth=plot_params['linewidth_corr'],linestyle=plot_params['linestyle_corr'],zorder=0)
ax4.axvline(df_results_trace.loc[df_results_trace['Element']=='Th','Best_lag'].values[0],
                    color=plot_params['linecolor_bstp_lag'],linewidth=plot_params['linewidth_bstp_lag'],
                    linestyle=plot_params['linestyle_bstp_lag'],zorder=1)
ax4.axvspan(df_results_trace.loc[df_results_trace['Element']=='Th','Lag_95%CI_low'].values[0],
                    df_results_trace.loc[df_results_trace['Element']=='Th','Lag_95%CI_high'].values[0],
                    color=plot_params['shade_color_CI'],alpha=0.3,zorder=1.5)
ax4.axvline(df_results_trace.loc[df_results_trace['Element']=='Th','Best_lag_observed'].values[0],
                    color=plot_params['linecolor_observed_lag'],linewidth=plot_params['linewidth_observed_lag'],
                    linestyle=plot_params['linestyle_observed_lag'],zorder=2)
ax4.set_ylim(-30,55)
ax4.tick_params(labelsize=plot_params['fontsize_tickmark'])  
for label in ax4.get_xticklabels() + ax4.get_yticklabels():
    label.set_fontname(plot_params['fontname_all'])
ax4.set_ylabel('Correlation',fontsize=plot_params['fontsize_label'],font=plot_params['fontname_all'])
ax4.text(52,45,"Th",
                fontsize=plot_params['fontsize_label'],fontname=plot_params['fontname_all'])
ax4.text(-59,45,f"bootstrap mean lag = {df_results_trace.loc[df_results_trace['Element']=='Th','Best_lag'].values[0]:.1f}\nbest lag observed = {df_results_trace.loc[df_results_trace['Element']=='Th','Best_lag_observed'].values[0]}",
                fontsize=plot_params['fontsize_label'],fontname=plot_params['fontname_all'],
                bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.2'))
fig.text(0.51, 0.683, 'D', fontsize=plot_params['fontsize_letter'], fontweight='bold', fontname=plot_params['fontname_all'])

ax5 = axes[2][0]
ax5.plot(lags,corr_observed_trace8['La8'],color=plot_params['linecolor_corr'],
          marker=plot_params['marker_corr'],markerfacecolor=plot_params['marker_fc_corr'],
          linewidth=plot_params['linewidth_corr'],linestyle=plot_params['linestyle_corr'],zorder=0)
ax5.axvline(df_results_trace8.loc[df_results_trace8['Element']=='La8','Best_lag'].values[0],
                    color=plot_params['linecolor_bstp_lag'],linewidth=plot_params['linewidth_bstp_lag'],
                    linestyle=plot_params['linestyle_bstp_lag'],zorder=1)
ax5.axvspan(df_results_trace8.loc[df_results_trace8['Element']=='La8','Lag_95%CI_low'].values[0],
                    df_results_trace8.loc[df_results_trace8['Element']=='La8','Lag_95%CI_high'].values[0],
                    color=plot_params['shade_color_CI'],alpha=0.3,zorder=1.5)
ax5.axvline(df_results_trace8.loc[df_results_trace8['Element']=='La8','Best_lag_observed'].values[0],
                    color=plot_params['linecolor_observed_lag'],linewidth=plot_params['linewidth_observed_lag'],
                    linestyle=plot_params['linestyle_observed_lag'],zorder=2)
ax5.set_ylim(-35,52)
ax5.tick_params(labelsize=plot_params['fontsize_tickmark'])  
for label in ax5.get_xticklabels() + ax5.get_yticklabels():
    label.set_fontname(plot_params['fontname_all'])
ax5.set_ylabel('Correlation',fontsize=plot_params['fontsize_label'],font=plot_params['fontname_all'])
ax5.text(46,40,"La$_{{\mathrm{{8.0}}}}$",
                fontsize=plot_params['fontsize_label'],fontname=plot_params['fontname_all'])
ax5.text(-59,40,f"bootstrap mean lag = {df_results_trace8.loc[df_results_trace8['Element']=='La8','Best_lag'].values[0]:.1f}\nbest lag observed = {df_results_trace8.loc[df_results_trace8['Element']=='La8','Best_lag_observed'].values[0]}",
                fontsize=plot_params['fontsize_label'],fontname=plot_params['fontname_all'],
                bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.2'))
fig.text(0.08, 0.495, 'E', fontsize=plot_params['fontsize_letter'], fontweight='bold', fontname=plot_params['fontname_all'])

ax6 = axes[2][1]
ax6.plot(lags,corr_observed_trace8['K8'],color=plot_params['linecolor_corr'],
          marker=plot_params['marker_corr'],markerfacecolor=plot_params['marker_fc_corr'],
          linewidth=plot_params['linewidth_corr'],linestyle=plot_params['linestyle_corr'],zorder=0)
ax6.axvline(df_results_trace8.loc[df_results_trace8['Element']=='K8','Best_lag'].values[0],
                    color=plot_params['linecolor_bstp_lag'],linewidth=plot_params['linewidth_bstp_lag'],
                    linestyle=plot_params['linestyle_bstp_lag'],zorder=1)
ax6.axvspan(df_results_trace8.loc[df_results_trace8['Element']=='K8','Lag_95%CI_low'].values[0],
                    df_results_trace8.loc[df_results_trace8['Element']=='K8','Lag_95%CI_high'].values[0],
                    color=plot_params['shade_color_CI'],alpha=0.3,zorder=1.5)
ax6.axvline(df_results_trace8.loc[df_results_trace8['Element']=='K8','Best_lag_observed'].values[0],
                    color=plot_params['linecolor_observed_lag'],linewidth=plot_params['linewidth_observed_lag'],
                    linestyle=plot_params['linestyle_observed_lag'],zorder=2)
ax6.set_ylim(-35,52)
ax6.tick_params(labelsize=plot_params['fontsize_tickmark'])  
for label in ax6.get_xticklabels() + ax6.get_yticklabels():
    label.set_fontname(plot_params['fontname_all'])
ax6.set_ylabel('Correlation',fontsize=plot_params['fontsize_label'],font=plot_params['fontname_all'])
ax6.text(41,40,"K$_2$O$_{{\mathrm{{8.0}}}}$",
                fontsize=plot_params['fontsize_label'],fontname=plot_params['fontname_all'])
ax6.text(-59,40,f"bootstrap mean lag = {df_results_trace8.loc[df_results_trace8['Element']=='K8','Best_lag'].values[0]:.1f}\nbest lag observed = {df_results_trace8.loc[df_results_trace8['Element']=='K8','Best_lag_observed'].values[0]}",
                fontsize=plot_params['fontsize_label'],fontname=plot_params['fontname_all'],
                bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.2'))
fig.text(0.51, 0.495, 'F', fontsize=plot_params['fontsize_letter'], fontweight='bold', fontname=plot_params['fontname_all'])

ax7 = axes[3][0]
ax7.plot(lags,corr_observed_trace8['Rb8'],color=plot_params['linecolor_corr'],
          marker=plot_params['marker_corr'],markerfacecolor=plot_params['marker_fc_corr'],
          linewidth=plot_params['linewidth_corr'],linestyle=plot_params['linestyle_corr'],zorder=0)
ax7.axvline(df_results_trace8.loc[df_results_trace8['Element']=='Rb8','Best_lag'].values[0],
                    color=plot_params['linecolor_bstp_lag'],linewidth=plot_params['linewidth_bstp_lag'],
                    linestyle=plot_params['linestyle_bstp_lag'],zorder=1)
ax7.axvspan(df_results_trace8.loc[df_results_trace8['Element']=='Rb8','Lag_95%CI_low'].values[0],
                    df_results_trace8.loc[df_results_trace8['Element']=='Rb8','Lag_95%CI_high'].values[0],
                    color=plot_params['shade_color_CI'],alpha=0.3,zorder=1.5)
ax7.axvline(df_results_trace8.loc[df_results_trace8['Element']=='Rb8','Best_lag_observed'].values[0],
                    color=plot_params['linecolor_observed_lag'],linewidth=plot_params['linewidth_observed_lag'],
                    linestyle=plot_params['linestyle_observed_lag'],zorder=2)
ax7.set_ylim(-40,55)
ax7.tick_params(labelsize=plot_params['fontsize_tickmark'])  
for label in ax7.get_xticklabels() + ax7.get_yticklabels():
    label.set_fontname(plot_params['fontname_all'])
ax7.set_ylabel('Correlation',fontsize=plot_params['fontsize_label'],font=plot_params['fontname_all'])
ax7.set_xlabel('Lag (1 lag = 2 kyr)',fontsize=plot_params['fontsize_label'],font=plot_params['fontname_all'])
ax7.text(45,42,"Rb$_{{\mathrm{{8.0}}}}$",
                fontsize=plot_params['fontsize_label'],fontname=plot_params['fontname_all'])
ax7.text(-59,42,f"bootstrap mean lag = {df_results_trace8.loc[df_results_trace8['Element']=='Rb8','Best_lag'].values[0]:.1f}\nbest lag observed = {df_results_trace8.loc[df_results_trace8['Element']=='Rb8','Best_lag_observed'].values[0]}",
                fontsize=plot_params['fontsize_label'],fontname=plot_params['fontname_all'],
                bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.2'))
fig.text(0.08, 0.307, 'G', fontsize=plot_params['fontsize_letter'], fontweight='bold', fontname=plot_params['fontname_all'])

ax8 = axes[3][1]
ax8.plot(lags,corr_observed_trace8['Th8'],color=plot_params['linecolor_corr'],
          marker=plot_params['marker_corr'],markerfacecolor=plot_params['marker_fc_corr'],
          linewidth=plot_params['linewidth_corr'],linestyle=plot_params['linestyle_corr'],zorder=0)
ax8.axvline(df_results_trace8.loc[df_results_trace8['Element']=='Th8','Best_lag'].values[0],
                    color=plot_params['linecolor_bstp_lag'],linewidth=plot_params['linewidth_bstp_lag'],
                    linestyle=plot_params['linestyle_bstp_lag'],zorder=1)
ax8.axvspan(df_results_trace8.loc[df_results_trace8['Element']=='Th8','Lag_95%CI_low'].values[0],
                    df_results_trace8.loc[df_results_trace8['Element']=='Th8','Lag_95%CI_high'].values[0],
                    color=plot_params['shade_color_CI'],alpha=0.3,zorder=1.5)
ax8.axvline(df_results_trace8.loc[df_results_trace8['Element']=='Th8','Best_lag_observed'].values[0],
                    color=plot_params['linecolor_observed_lag'],linewidth=plot_params['linewidth_observed_lag'],
                    linestyle=plot_params['linestyle_observed_lag'],zorder=2)
ax8.set_ylim(-35,52)
ax8.tick_params(labelsize=plot_params['fontsize_tickmark'])  
for label in ax8.get_xticklabels() + ax8.get_yticklabels():
    label.set_fontname(plot_params['fontname_all'])
ax8.set_ylabel('Correlation',fontsize=plot_params['fontsize_label'],font=plot_params['fontname_all'])
ax8.set_xlabel('Lag (1 lag = 2 kyr)',fontsize=plot_params['fontsize_label'],font=plot_params['fontname_all'])
ax8.text(45,40,"Th$_{{\mathrm{{8.0}}}}$",
                fontsize=plot_params['fontsize_label'],fontname=plot_params['fontname_all'])
ax8.text(-59,40,f"bootstrap mean lag = {df_results_trace8.loc[df_results_trace8['Element']=='Th8','Best_lag'].values[0]:.1f}\nbest lag observed = {df_results_trace8.loc[df_results_trace8['Element']=='Th8','Best_lag_observed'].values[0]}",
                fontsize=plot_params['fontsize_label'],fontname=plot_params['fontname_all'],
                bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.2'))
fig.text(0.51, 0.307, 'H', fontsize=plot_params['fontsize_letter'], fontweight='bold', fontname=plot_params['fontname_all'])

