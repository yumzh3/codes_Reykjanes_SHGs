# calculate the liquidus temperature to start crystallization
# temperature decrease by 1 Celsius per step
# calculate the liquid fraction, phase fractions, liquid and phase compositions in wt.% per each step during crystallization
# calculte either fractional or equilibrium crystallization
# detailed algorithm is introduced in Weaver and Langmuir 1990
# originally written by Jocelyn Fuentes 2016
# modified by Mingzhen Yu 2021: add Ni and Mn in the system
# last modified: Sep 26, 2025

from wl1989stoich_2023 import *
from wl1989kdcalc_2023 import *
from wlState_2023 import *
import numpy as np
import math


# default values used to calculate the phase proportions called by function 'wlState_2023'
ta = {'cpx':1., 'plg':1., 'ol':2./3.}
uaj_plg = {'CaAl2O4':5./3., 'NaAlO2':2.5, 'MgO':0., 'FeO':0., 'CaSiO3':0., 'TiO2':0., 'KAlO2':0., 'PO52':0., 'MnO':0.,'NiO':0.} # CaAl2Si2O8, NaAlSi3O8
uaj_ol = {'CaAl2O4':0., 'NaAlO2':0., 'MgO':1, 'FeO':1, 'CaSiO3':0., 'TiO2':0., 'KAlO2':0., 'PO52':0., 'MnO':1., 'NiO':1.} # (Mg,Fe,Ni,Mn)2SiO4
uaj_cpx = {'CaAl2O4':4./3., 'NaAlO2':2., 'MgO':2., 'FeO':2., 'CaSiO3':1., 'TiO2':1., 'KAlO2':0., 'PO52':0., 'MnO':2., 'NiO':2.} # CaAl2SiO6, NaAlSi2O6, MgSiO3, FeSiO3, CaSiO3, MnSiO3, NiSiO3, CaTiO3
uaj = {'ol':uaj_ol, 'plg':uaj_plg, 'cpx':uaj_cpx}
kd_ol_trace = {'La':1*10**(-5),'Ce':2*10**(-5),'Nd':1*10**(-4),'Sm':4*10**(-4),'Eu':7*10**(-4),'Gd':0.001,'Dy':0.004,'Er':0.01,'Yb':0.021,'Lu':0.031,'Sr':0,'Th':0,'Ba':0,'Rb':0,'Nb':0,'Ta':0,'D=0':0,'Pr':5*10**(-5),'Pb':0.003,'U':0.00038}  ## Kd for ol from O'Neill and Jenner 2012 Nature supplement and references in. Pb and U from Salters and Stracke (2004)
kd_plg_trace = {'La':0.0348,'Ce':0.0278,'Nd':0.0179,'Sm':0.0132,'Eu':0.20,'Gd':0.0125,'Dy':0.0112,'Er':0.0116,'Yb':0.0155,'Lu':0.012,'Sr':2,'Th':0.01,'Ba':0.16,'Rb':0.016,'Nb':0.003,'Ta':0.003,'D=0':0,'Pr':0.05,'Pb':0.144,'U':0.012}  # Kd for pl from O'Neill and Jenner 2012 Nature supplement and references in. Eu and Th to U are from Bedard (2023) with An=0.84
kd_cpx_trace = {'La':0.081,'Ce':0.129,'Nd':0.26,'Sm':0.417,'Eu':0.479,'Gd':0.531,'Dy':0.596,'Er':0.597,'Yb':0.558,'Lu':0.532,'Sr':0.06,'Th':0.01,'Ba':0.02,'Rb':0,'Nb':0.019,'Ta':0.06,'D=0':0,'Pr':0.26,'Pb':0.01,'U':0.005}  # Kd for cpx from O'Neill and Jenner 2012 Nature supplement and references in. Ba, Rb, Ta, Nb and Pr from Hill et al. (2000). Th and U are from LaTourretter and Burnett (1992). Pb from Salters and Stracke (2004).
kd_trace = {'ol':kd_ol_trace, 'plg':kd_plg_trace, 'cpx':kd_cpx_trace}

# calculate the liquidus T (in Kelvin)
def get_first_T(system_components, P, kdCalc = kdCalc_langmuir1992):
    firstT = 2000.  # a guess for liquidus T
    deltaT = 100.
    qa, fa, major_liquid_components, solid_phase_components, num_iter = state(system_components,firstT,uaj, ta, P=P, kdCalc= kdCalc)
    fl = 1-sum(fa.values()) # liquid fraction in the system
    if num_iter == 3000:
        print('MAX ITERATION!')
    while (fl == 1.) or (deltaT > 1.):
        if fl == 1.:
            firstT = firstT-deltaT
        elif (fl < 1.) and (deltaT > 1.): # fl<1 means minerals start to xtalized already, means that deltaT is too large now
            firstT = firstT+deltaT
            deltaT = deltaT/10.
            firstT=firstT-deltaT
        qa, fa, major_liquid_components, solid_phase_components, num_iter = state(system_components,firstT,uaj, ta, P=P, kdCalc= kdCalc)
        fl = 1-sum(fa.values())
        if num_iter == 3000:
            print('MAX ITERATION!')
            firstT = 2000.
    return firstT

# calculate pure fractional xtalization including liquid fraction, phase fractions, liquid and phase compositions in wt.% and ppm
def frac_model_trange(t_start, t_stop, tstep, major_start_comp, trace_start_comp, kd_trace, P, kdCalc = kdCalc_langmuir1992):
    # tstep = 0.2
    bulk_d = {key:0. for key in trace_start_comp}
    trange = np.arange(t_start,t_stop, -tstep)
    system_components = oxideToComponent(major_start_comp)  # input and output are dictionary
    major_liquid_components = system_components.copy()
    trace_liquid_comp = trace_start_comp.copy()
    major_oxide_dict = {key:[] for key in major_start_comp}
    major_phase_oxide_dict = {phase:{key:[] for key in major_start_comp} for phase in ['ol','cpx','plg']}
    major_phase_oxides = {phase:[] for phase in ['ol','cpx','plg']}
    trace_dict = {key:[] for key in trace_start_comp}
    fl = []
    fa_dict = {phase:[] for phase in ['plg', 'cpx', 'ol']}
    for i in range(len(trange)):
        ## Major Elements
        if i == 0:
            qa, fa, major_liquid_components, major_phase_components, num_iter = state(major_liquid_components,trange[i],uaj, ta, P=P, kdCalc = kdCalc)
            for phase in fa:
                fa_dict[phase].append(fa[phase])
        else:
            major_liquid_components = oxideToComponent(major_oxides)
            qa, fa, major_liquid_components, major_phase_components, num_iter = state(major_liquid_components,trange[i],uaj, ta, P=P, kdCalc = kdCalc)
            for phase in fa:
                solid_phase = fa[phase]*fl[-1]+fa_dict[phase][-1]
                fa_dict[phase].append(solid_phase)
        major_oxides = cationFracToWeight(major_liquid_components)
        for phase in major_phase_oxides:
            major_phase_oxides[phase] = cationFracToWeight(major_phase_components[phase])
            for key in major_phase_oxides[phase]:
                major_phase_oxide_dict[phase][key].append(major_phase_oxides[phase][key])
        liq = (1. - sum(fa.values()))
        if i == 0:
            fl.append(liq)
        else:
            fl.append(liq*fl[-1])
        fa_tot = sum(fa.values())
        for key in major_oxides:
            major_oxide_dict[key].append(major_oxides[key])
        ## Trace Elements
        for elem in trace_start_comp:
            bulk_d[elem] = 0.  # Calculate Bulk D
            if fa_tot != 0.:
                for phase in fa:
                    bulk_d[elem] += (fa[phase]/fa_tot)*kd_trace[phase][elem]
            trace_liquid_comp[elem] = trace_liquid_comp[elem]/(liq +(1.-liq)*bulk_d[elem])  # Add erupted composition to eruption dictionary
            trace_dict[elem].append(trace_liquid_comp[elem]) 
    return fl, fa_dict, major_oxide_dict, major_phase_oxide_dict, trace_dict   

# calculate pure equilibrium xtalization including liquid fraction, phase fractions, liquid and phase compositions in wt.% and ppm
def eq_model_trange(t_start, t_stop, tstep, major_start_comp, trace_start_comp, kd_trace, P, kdCalc = kdCalc_langmuir1992):
    # tstep = 0.5
    bulk_d = {key:0. for key in trace_start_comp}
    trange = np.arange(t_start,t_stop, -tstep)
    system_components = oxideToComponent(major_start_comp)
    major_oxide_dict = {key:[] for key in major_start_comp}
    major_phase_oxide_dict = {phase:{key:[] for key in major_start_comp} for phase in ['ol','cpx','plg']}
    major_phase_oxides = {phase:[] for phase in ['ol','cpx','plg']}
    trace_dict = {key:[] for key in trace_start_comp}
    fl = []
    fa_dict = {phase:[] for phase in ['plg', 'cpx', 'ol']}
    for i in range(len(trange)):
        ## Major Elements
        qa, fa, major_liquid_components, major_phase_components, num_iter = state(system_components,trange[i],uaj, ta, P = P, kdCalc = kdCalc)
        for phase in fa:
            fa_dict[phase].append(fa[phase])
        major_oxides = cationFracToWeight(major_liquid_components)
        liq = 1.-sum(fa.values())
        fl.append(liq)
        fa_tot = sum(fa.values())
        for key in major_oxides:
            major_oxide_dict[key].append(major_oxides[key])
        for phase in major_phase_oxides:
            major_phase_oxides[phase] = cationFracToWeight(major_phase_components[phase])
            for key in major_phase_oxides[phase]:
                major_phase_oxide_dict[phase][key].append(major_phase_oxides[phase][key])
        ## Trace Elements
        for elem in trace_start_comp:
            bulk_d[elem] = 0.   # Calculate Bulk D
            if fa_tot != 0.:
                for phase in fa:
                    bulk_d[elem] += (fa[phase]/fa_tot)*kd_trace[phase][elem]
            trace_dict[elem].append(trace_start_comp[elem]/(liq +(1.-liq)*bulk_d[elem]))  # Add erupted composition to eruption dictionary
    return fl, fa_dict, major_oxide_dict, major_phase_oxide_dict, trace_dict

# Fractional Crystallization for PRMC- System Components never Change
# calculate liquid fraction, phase fractions, liquid and phase compositions in wt.% and ppm
def frac_model_trange_PRMC(f_target,t_start, tstep, major_start_comp, trace_start_comp, kd_trace, P, kdCalc = kdCalc_langmuir1992):
    # tstep = 0.2
    bulk_d = {key:0. for key in trace_start_comp}
    # trange = np.arange(t_stop,t_start, tstep)
    system_components = oxideToComponent(major_start_comp)
    major_liquid_components = system_components.copy()
    trace_liquid_comp = trace_start_comp.copy()
    major_oxide_dict = {key:[] for key in major_start_comp}
    major_phase_oxide_dict = {phase:{key:[] for key in major_start_comp} for phase in ['ol','cpx','plg']}
    major_phase_oxides = {phase:[] for phase in ['ol','cpx','plg']}
    trace_dict = {key:[] for key in trace_start_comp}
    fl = [1]
    fa_dict = {phase:[] for phase in ['plg', 'cpx', 'ol']}
    liq = 1
    T = t_start
    T_LLD = []
    while fl[-1] >= f_target:
        ## Major Elements
        if liq == 1:
            qa, fa, major_liquid_components, major_phase_components, num_iter = state(major_liquid_components,T,uaj, ta, P=P, kdCalc = kdCalc)
            for phase in fa:
                fa_dict[phase].append(fa[phase])
        else:
            major_liquid_components = oxideToComponent(major_oxides)
            qa, fa, major_liquid_components, major_phase_components, num_iter = state(major_liquid_components,T,uaj, ta, P=P, kdCalc = kdCalc)
            for phase in fa:
                solid_phase = fa[phase]*fl[-1]+fa_dict[phase][-1]
                fa_dict[phase].append(solid_phase)
        major_oxides = cationFracToWeight(major_liquid_components)
        for phase in major_phase_oxides:
            major_phase_oxides[phase] = cationFracToWeight(major_phase_components[phase])
            for key in major_phase_oxides[phase]:
                major_phase_oxide_dict[phase][key].append(major_phase_oxides[phase][key])
        if liq == 1:
            liq = (1. - sum(fa.values()))
            fl[0] = liq
        else:
            liq = (1. - sum(fa.values())) 
            fl.append(liq*fl[-1])
        fa_tot = sum(fa.values())
        for key in major_oxides:
            major_oxide_dict[key].append(major_oxides[key])
        ## Trace Elements
        for elem in trace_start_comp:
            bulk_d[elem] = 0.  # Calculate Bulk D
            if fa_tot != 0.:
                for phase in fa:
                    bulk_d[elem] += (fa[phase]/fa_tot)*kd_trace[phase][elem]
            trace_liquid_comp[elem] = trace_liquid_comp[elem]/(liq +(1.-liq)*bulk_d[elem])  # Add erupted composition to eruption dictionary
            trace_dict[elem].append(trace_liquid_comp[elem])
        T_LLD.append(T)
        if fl[-1]-f_target <= 0.003:
            T = T-0.01
        else:
            T = T-tstep
    return fl, fa_dict, major_oxide_dict, major_phase_oxide_dict, trace_dict, T_LLD

# Equilibrium Crystallization for PRMC- System Components never Change 
# calculate liquid fraction, phase fractions, liquid and phase compositions in wt.% and ppm
def eq_model_trange_PRMC(f_target, t_start, tstep, major_start_comp, trace_start_comp, kd_trace, P, kdCalc = kdCalc_langmuir1992):
    # tstep = 0.5
    bulk_d = {key:0. for key in trace_start_comp}
    # trange = np.arange(t_stop,t_start, tstep)
    system_components = oxideToComponent(major_start_comp)
    major_oxide_dict = {key:[] for key in major_start_comp}
    major_phase_oxide_dict = {phase:{key:[] for key in major_start_comp} for phase in ['ol','cpx','plg']}
    major_phase_oxides = {phase:[] for phase in ['ol','cpx','plg']}
    trace_dict = {key:[] for key in trace_start_comp}
    fl = []
    fa_dict = {phase:[] for phase in ['plg', 'cpx', 'ol']}
    liq = 1
    T = t_start
    T_LLD = []
    while liq >= f_target:
        ## Major Elements
        qa, fa, major_liquid_components, major_phase_components, num_iter = state(system_components,T,uaj, ta, P = P, kdCalc = kdCalc)
        for phase in fa:
            fa_dict[phase].append(fa[phase])
        major_oxides = cationFracToWeight(major_liquid_components)
        liq = 1.-sum(fa.values())
        fl.append(liq)
        fa_tot = sum(fa.values())
        for key in major_oxides:
            major_oxide_dict[key].append(major_oxides[key])
        for phase in major_phase_oxides:
            major_phase_oxides[phase] = cationFracToWeight(major_phase_components[phase])
            for key in major_phase_oxides[phase]:
                major_phase_oxide_dict[phase][key].append(major_phase_oxides[phase][key])
        ## Trace Elements
        for elem in trace_start_comp:
            bulk_d[elem] = 0.   # Calculate Bulk D
            if fa_tot != 0.:
                for phase in fa:
                    bulk_d[elem] += (fa[phase]/fa_tot)*kd_trace[phase][elem]
            trace_dict[elem].append(trace_start_comp[elem]/(liq +(1.-liq)*bulk_d[elem]))  # Add erupted composition to eruption dictionary
        T_LLD.append(T)
        if liq-f_target <= 0.003:
            T = T-0.01
        else:
            T = T-tstep
    return fl, fa_dict, major_oxide_dict, major_phase_oxide_dict, trace_dict, T_LLD


