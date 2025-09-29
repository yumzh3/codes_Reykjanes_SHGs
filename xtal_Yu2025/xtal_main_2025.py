# modeling for fractional and equilibrium crystallization
# calculate the liquid fraction, phase fractions, liquid and phase compositions in wt.% per each step during crystallization
# calculte either fractional or equilibrium crystallization
# results are saved in data frames "LLD_df_frac" and "LLD_df_equ"
# written by: Mingzhen Yu
# last modified: Sep 26, 2025

from wl1989stoich_2023 import *
from wl1989kdcalc_2023 import *
from wl1989models_2023 import *
from wlState_2023 import *
import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt
import random


# %% load starting compositions
LLDinput = pd.read_csv('your location/parental_magma_data.csv')
pm_location = 'Reykjanes_Murton2002'  # can add more lines of different starting compositions in the data file
index = np.where(LLDinput['location'] == pm_location)[0][0]
P = 1 # pressure for crystallization, unit is bar  

# %% pure fractional crystallization 
magma = LLDinput.loc[index,['SiO2','TiO2','Al2O3','FeO','MgO','K2O','MnO','Na2O','P2O5','CaO','NiO']]
magma['FeO'] = magma['FeO']*0.9  # assume FeO/FeO(t) = 0.9 for MORB
magma_trace = LLDinput.loc[index,['La','Ce','Nd','Sm','Eu','Gd','Dy','Er','Yb','Lu','Sr','Th','Ba','Rb','Ta','Nb','D=0','Pr','Pb','U']]
system_components = magma.to_dict()
T_system_components = oxideToComponent(system_components)
t_start = get_first_T(T_system_components, P, kdCalc = kdCalc_langmuir1992)
t_stop = t_start - 250
tstep = 0.5
trace_components = magma_trace.to_dict()
fl,fa_dict,major_oxide_dict,major_phase_oxide_dict,trace_dict = frac_model_trange(t_start, t_stop, tstep, system_components,trace_components,kd_trace,P,kdCalc = kdCalc_langmuir1992) 

T_df = pd.DataFrame(np.arange(t_start,t_stop,-tstep))
T_df = T_df-273.15
T_df.columns = ['T_C']
fl_dict = {'fl':fl}
fl_df = pd.DataFrame(fl_dict)
fa_df = pd.DataFrame(fa_dict)
major_oxide_df = pd.DataFrame(major_oxide_dict)
major_ol_oxide_df = pd.DataFrame(major_phase_oxide_dict['ol'])
major_cpx_oxide_df = pd.DataFrame(major_phase_oxide_dict['cpx'])
major_plg_oxide_df = pd.DataFrame(major_phase_oxide_dict['plg'])
trace_dict_df = pd.DataFrame(trace_dict)

LLD_df_frac = pd.concat([T_df,fl_df,fa_df,major_oxide_df,major_ol_oxide_df,major_cpx_oxide_df,major_plg_oxide_df,trace_dict_df],axis=1)
LLD_df_frac.columns = ['T_C','f_liq','f_plg','f_cpx','f_ol','liq_SiO2','liq_TiO2','liq_Al2O3','liq_FeO',\
               'liq_MgO','liq_K2O','liq_MnO','liq_Na2O','liq_P2O5','liq_CaO','liq_NiO','olSiO2',\
                   'olTiO2','olAl2O3','olFeO','olMgO','olK2O','olMnO','olNa2O','olP2O5','olCaO','olNiO',\
                       'cpxSiO2','cpxTiO2','cpxAl2O3','cpxFeO','cpxMgO','cpxK2O','cpxMnO','cpxNa2O',\
                           'cpxP2O5','cpxCaO','cpxNiO','plgSiO2','plgTiO2','plgAl2O3','plgFeO','plgMgO',\
                               'plgK2O','plgMnO','plgNa2O','plgP2O5','plgCaO','plgNiO','liq_La','liq_Ce',\
                                   'liq_Nd','liq_Sm','liq_Eu','liq_Gd','liq_Dy','liq_Er','liq_Yb','liq_Lu','liq_Sr',\
                                       'liq_Th','liq_Ba','liq_Rb','liq_Ta','liq_Nb','liq_D=0','liq_Pr','liq_Pb','liq_U']
LLD_df_frac['liq_FeOt'] = LLD_df_frac['liq_FeO']/0.9  # assume FeO/FeO(t) = 0.9 for MORB
LLD_df_frac['Fo'] = 100/(1+LLD_df_frac['olFeO']/LLD_df_frac['olMgO']*40.3/71.84)
LLD_df_frac['olNippm'] = LLD_df_frac['olNiO']*58.6934/74.69*10**4
LLD_df_frac['olMnppm'] = LLD_df_frac['olMnO']*54.938/70.94*10**4
LLD_df_frac['liq_Ni'] = LLD_df_frac['liq_NiO']*58.6934/74.69*10**4
LLD_df_frac['liq_FeOtMnO'] = LLD_df_frac['liq_FeOt']/LLD_df_frac['liq_MnO']


# %% pure equilibrium crystallization 
magma = LLDinput.loc[index,['SiO2','TiO2','Al2O3','FeO','MgO','K2O','MnO','Na2O','P2O5','CaO','NiO']]
magma['FeO'] = magma['FeO']*0.9  # assume FeO/FeO(t) = 0.9 for MORB
magma_trace = LLDinput.loc[index,['La','Ce','Nd','Sm','Eu','Gd','Dy','Er','Yb','Lu','Sr','Th','Ba','Rb','Ta','Nb','D=0','Pr','Pb','U']]
system_components = magma.to_dict()
T_system_components = oxideToComponent(system_components)
t_start = get_first_T(T_system_components, P, kdCalc = kdCalc_langmuir1992)
t_stop = t_start -250
tstep = 0.5
trace_components = magma_trace.to_dict()
fl,fa_dict,major_oxide_dict,major_phase_oxide_dict,trace_dict = eq_model_trange(t_start, t_stop,tstep,system_components,trace_components,kd_trace,P,kdCalc = kdCalc_langmuir1992) 

T_df = pd.DataFrame(np.arange(t_start,t_stop,-tstep))
T_df = T_df-273.15
T_df.columns = ['T_C']
fl_dict = {'fl':fl}
fl_df = pd.DataFrame(fl_dict)
fa_df = pd.DataFrame(fa_dict)
major_oxide_df = pd.DataFrame(major_oxide_dict)
major_ol_oxide_df = pd.DataFrame(major_phase_oxide_dict['ol'])
major_cpx_oxide_df = pd.DataFrame(major_phase_oxide_dict['cpx'])
major_plg_oxide_df = pd.DataFrame(major_phase_oxide_dict['plg'])
trace_dict_df = pd.DataFrame(trace_dict)

LLD_df_equ = pd.concat([T_df,fl_df,fa_df,major_oxide_df,major_ol_oxide_df,major_cpx_oxide_df,major_plg_oxide_df,trace_dict_df],axis=1)
LLD_df_equ.columns = ['T_C','f_liq','f_plg','f_cpx','f_ol','liq_SiO2','liq_TiO2','liq_Al2O3','liq_FeO',\
               'liq_MgO','liq_K2O','liq_MnO','liq_Na2O','liq_P2O5','liq_CaO','liq_NiO','olSiO2',\
                   'olTiO2','olAl2O3','olFeO','olMgO','olK2O','olMnO','olNa2O','olP2O5','olCaO','olNiO',\
                       'cpxSiO2','cpxTiO2','cpxAl2O3','cpxFeO','cpxMgO','cpxK2O','cpxMnO','cpxNa2O',\
                           'cpxP2O5','cpxCaO','cpxNiO','plgSiO2','plgTiO2','plgAl2O3','plgFeO','plgMgO',\
                               'plgK2O','plgMnO','plgNa2O','plgP2O5','plgCaO','plgNiO','liq_La','liq_Ce',\
                                   'liq_Nd','liq_Sm','liq_Eu','liq_Gd','liq_Dy','liq_Er','liq_Yb','liq_Lu','liq_Sr',\
                                       'liq_Th','liq_Ba','liq_Rb','liq_Ta','liq_Nb','liq_D=0','liq_Pr','liq_Pb','liq_U']
LLD_df_equ['liq_FeOt'] = LLD_df_equ['liq_FeO']/0.9  # assume FeO/FeO(t) = 0.9 for MORB
LLD_df_equ['Fo'] = 100/(1+LLD_df_equ['olFeO']/LLD_df_equ['olMgO']*40.3/71.84)
LLD_df_equ['olNippm'] = LLD_df_equ['olNiO']*58.6934/74.69*10**4
LLD_df_equ['olMnppm'] = LLD_df_equ['olMnO']*54.938/70.94*10**4
LLD_df_equ['liq_Ni'] = LLD_df_equ['liq_NiO']*58.6934/74.69*10**4
LLD_df_equ['liq_FeOtMnO'] = LLD_df_equ['liq_FeOt']/LLD_df_equ['liq_MnO']

    


