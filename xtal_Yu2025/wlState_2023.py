# determine the phase saturation status and proportions in the system 
# detailed algorithm is introduced in Weaver and Langmuir 1990
# originally written by Jocelyn Fuentes 2016
# modified by Mingzhen Yu 2021: add Ni and Mn in the system
# last modified: Sep 30, 2025

from wl1989stoich_2023 import *
from wl1989kdcalc_2023 import *
import numpy as np
import math
 

fa_guess = {'plg':0., 'ol':0., 'cpx':0.}    
def state(system_components,T, uaj, ta, P, kdCalc = kdCalc_langmuir1992):
    """State determines the liquid composition and phases present in the system
    at a given temperature and possible pressure (depending on the Kd formula).
    It is possible to also pass a guess or liquid components. If none are given,
    then it is assumed there are no phases present and that the liquid components
    are the system components. T is in Kelvin and P is in bars.
    
    This is used for all of the major elements.
    System components must include SiO2, TiO2, Na2O, MgO, FeO, CaO, Al2O3, K2O,
    MnO, and P2O5, NiO.
    """
    liquid_components = system_components.copy()
    max_iter = 3000
    qa = {'plg':0., 'ol':0., 'cpx':0.}
    fa = {'plg':0., 'ol':0., 'cpx':0.}
    dfa = {}
    rj = {}
    tolerance = np.power(10.,-5.)  # 10**(-5)
    a = True  # Variable 'a' tells it whether or not to break the loop
    kdaj = kdCalc(liquid_components, T, P)  # Calculate Kd using liquid components
    phase_list = []
    liquid_components = {component:0. for component in kdaj['cpx']}
    for component in kdaj['cpx']:  # Calculate new solid fractions in equilibrium with the current state
        rj[component] = calculate_Rj(fa, kdaj, component)
        liquid_components[component] = rj[component]*system_components[component]
    for phase in fa:  # Calculate the initial saturation for each phase
        qa[phase] = calculate_Qa(liquid_components, kdaj, phase, ta, uaj)
        if (qa[phase]>0) or (fa[phase]>0):
            phase_list.append(phase)
    if len(phase_list) == 0:  # If there is no phase that is saturated, no need to enter the loop
        a = False
        solid_phase_components = {phase:{key:0 for key in kdaj['cpx']} for phase in fa}
    i = 0
    while (a == True) and (i<max_iter):
        i += 1
        if len(phase_list) !=0:  # Use Newton Method to find new Fa if there are phases present
            # THE FOLLOWING CODE USES THE MATRIX METHOD
            pab_dict = create_Pab_dict(rj, kdaj, liquid_components, uaj, phase_list)
            dfa = solve_matrix(pab_dict, qa, phase_list)
            if dfa == 'Singular':
                print('Singular')
            fa_new = {}
            tst = 0.
            for phase in phase_list:
                fa_new[phase] = fa[phase] + dfa[phase]
                if fa_new[phase]<0:  # Check to make sure the new Fa is greater than 0 and less than 1
                    fa_new[phase] = 0.1*fa[phase]
                elif fa_new[phase]>1:
                    fa_new[phase]= 0.9 + .1*fa[phase]
                #tst += abs(fa[phase] - fa_new[phase])
                fa[phase] = fa_new[phase]
            # Recalculate Liquid Percent
#            x = tst/len(phase_list)
#            if x <= tolerance:
#                a = False
            for component in kdaj['cpx']:
                rj[component] = calculate_Rj(fa, kdaj, component)
                liquid_components[component] = rj[component]*system_components[component]
            phase_list = []  # Recalculate Saturation
            solid_phase_components = {phase:{key:0 for key in kdaj['cpx']} for phase in fa}
            kdaj = kdCalc(liquid_components, T, P)
            for component in kdaj['cpx']:
                rj[component] = calculate_Rj(fa, kdaj, component)
                liquid_components[component] = rj[component]*system_components[component]
            for phase in fa:
                qa[phase] = calculate_Qa(liquid_components, kdaj, phase, ta, uaj)
                if (qa[phase]>0) or (fa[phase]>0):
                    phase_list.append(phase)
                    for component in kdaj['cpx']:
                        solid_phase_components[phase][component] = liquid_components[component]*kdaj[phase][component]
            qa_new = [qa[x] for x in phase_list]  # Qa should be very closs to zero
            if all(np.abs(value) <= tolerance for value in qa_new):
                a = False
            fl = 1.-sum(fa.values())
            if (1.-fl)>=1:
                a = True
            elif (1.-fl)<0:
                a = True
        else:
            a = False
    return qa, fa,liquid_components, solid_phase_components, i
            
            
# The following functions are called in the codes above             
def calculate_Rj(fa, kd, component):  # Called by calculate_Pab and calculate_Qa
    temp = 0.
    for p in fa:
        temp += fa[p]*(kd[p][component]-1)
    rj = 1./(1.+temp)
    return rj
    
def calculate_Qa(clj, kd, phase, ta, uaj):  # Called by State
    # Initialize rj, clj, and caj
    caj = {'plg':{key:0 for key in kd['cpx']}, 
    'cpx':{key:0 for key in kd['cpx']}, 
    'ol':{key:0 for key in kd['cpx']}}
    # Given composition in component form and the Temperature,
    # calculate the saturation of a given phase.
    qa = -ta[phase]
    for component in kd['cpx']:  # Calculate the liquid composition and the composition of the phases
        caj[phase][component] = clj[component]*kd[phase][component]
        qa += uaj[phase][component]*caj[phase][component]
    return qa

def calculate_Pab(rj, kd, phase1, phase2, liquid_components,uaj):  # Called by create_pab_dict
    pab = 0.
    for component in kd['cpx']:
        pab += uaj[phase1][component]*liquid_components[component]*kd[phase1][component]*(kd[phase2][component]-1)*rj[component]
    return pab

def create_Pab_dict(rj, kd, liquid_components, uaj, phase_list):  # Called by State
    pab = {'plg':{}, 'cpx':{}, 'ol':{}}
    for phase1 in phase_list:
        for phase2 in phase_list:
            pab[phase1][phase2] = calculate_Pab(rj, kd, phase1, phase2, liquid_components,uaj)
    return pab
    
def solve_matrix(pab, qa, phase_list):  # Called by State
    # NEED TO ADD PART IN CASE OF SINGULAR MATRIX
    pab_array = np.zeros([len(phase_list),len(phase_list)])
    qa_array = np.zeros([len(phase_list),1])
    dfa = {}
    for i in range(0,len(phase_list)):
        qa_array[i] = qa[phase_list[i]]
        for j in range(0,len(phase_list)):
            pab_array[i,j] = pab[phase_list[i]][phase_list[j]]
    det = np.linalg.det(pab_array)
    if det == 0:
         dfa = 'Singular'
    # Solve Matrix
    else:
        dfa_array = np.dot(np.linalg.inv(pab_array),(qa_array))
        for k in range(0,len(phase_list)):  # Convert back to Dicitonaries
            dfa[phase_list[k]] = dfa_array[k][0]
    return dfa       
    
def newton(qa, fa, kd, system_components, uaj):
    epsilon = np.power(10, -14)  # Don't want to divide by a number smaller than epsilon
    fa_new = {}
    for phase1 in fa:
        qa_prime = 0
        for component in system_components:
            if component == 'SiO2':
                pass
            else:
                kb = 0
                for phase2 in fa:
                    kb+=kd[phase2][component]-1
                    rj = calculate_Rj(fa, kd, component)
                    qa_prime -= kd[phase1][component]*system_components[component]*rj*rj*kb
        #if np.abs(qa_prime) < epsilon:
            #return "Qa_prime Too Small"
        #else:
        fa_new[phase1] = fa[phase1]-(qa[phase1]/qa_prime)
    return fa_new       




