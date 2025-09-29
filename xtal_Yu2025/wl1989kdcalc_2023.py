# calculate partition coefficients for different components between minerals (ol. pl, cpx) and liquids used in crystallization calculations
# calculate based on the unit of cation mole
# crystallization involves olivine, plagioclase, and clinopyroxene
# equations are based on Weaver and Langmuir 1990
# originally written by Jocelyn Fuentes 2016
# modified by Mingzhen Yu 2021: add Ni and Mn in the system, and update KDFeMg(ol/l) with equation from Toplis 2005
# last modified: Sep 26, 2025

from wl1989stoich_2023 import *
import numpy as np
import math


# calculate partition coefficients for different components between minerals (ol. pl, cpx) and liquids
def kdCalc_langmuir1992(components, T, P, H = None):
    """
    This uses Langmuir et al. 1992 equations.
    P in bars, T in Kelvin.
    """
    keys = ['MgO', 'FeO', 'TiO2', 'PO52', 'MnO', 'CaAl2O4', 'NaAlO2', 'KAlO2', 'CaSiO3','NiO']
    cpx = {key:0. for key in keys}
    plg = {key:0. for key in keys}
    ol = {key:0. for key in keys}
    oxide_wt = cationFracToWeight(components)
    SiO2_molar = oxide_wt['SiO2']/60.083/1.625
    Na2O_molar = components['NaAlO2']/2*1.833*30.99/(30.99*2)/1.625
    K2O_molar = components['KAlO2']/2*1.833*47.1/(47.1*2)/1.625
    if SiO2_molar > 0.6:
        SiO2_adj = 100*SiO2_molar + 100*(Na2O_molar+K2O_molar)*(11-550/(100-100*SiO2_molar))*np.exp(-0.13*100*(Na2O_molar+K2O_molar))
    else:
        SiO2_adj = 100*SiO2_molar + 100*(Na2O_molar+K2O_molar)*((46/(100-100*SiO2_molar)-0.93)*100*(Na2O_molar+K2O_molar)-533/(100-100*SiO2_molar)+9.69)
    anorthite=components['CaAl2O4']/(components['CaAl2O4']+1.5*components['NaAlO2'])
    plg['CaAl2O4'] = np.power(10.,(2446./T) - (1.122  + 0.2562*anorthite))
    plg['NaAlO2'] = np.power(10.,((3195. + (3283.*anorthite) + (0.0506*P))/T) - (1.885*anorthite) -2.3715)
    plg['KAlO2'] = 0.15
    plg['MnO'] = 0.031
    plg['PO52'] = 0.1
    cpx['MgO'] = np.power(10.,(((3798. + (0.021*P))/T) - 2.28))
    cpx['FeO'] = cpx['MgO']*np.power(10.,-0.6198)
    cpx['CaSiO3'] = np.power(10.,((1783. + (0.0038*P))/T) -0.753)
    cpx['CaAl2O4'] = np.power(10.,((2418. + (0.068*P))/T) -2.3)
    cpx['NaAlO2'] = np.power(10.,((5087. + (0.073*P))/T) - 4.48)
    cpx['TiO2'] = np.power(10.,((1034. + (0.053*P))/T) - 1.27)
    cpx['KAlO2'] = 0.007
    cpx['PO52'] = 0.05
    ol['KAlO2'] = 0.001
    ol['PO52'] = 0.2
    ol['MgO'] = np.exp((6921./T) + (3.4*components['NaAlO2']/2.) + 
                (6.3*components['KAlO2']/2.) + (0.00001154*P) - 3.27)
    KDFeMg = np.exp(-6766./(8.3144*T)-7.34/8.3144+math.log(0.036*SiO2_adj-0.22)+3000*(1-2*components['MgO']*100*ol['MgO']/66.67)/(8.3144*T))  # from Toplis 2005
    ol['FeO'] = ol['MgO']*KDFeMg  # KdFe2(ol/l) with the unit of cation mole
    ol['MnO'] = 0.79*ol['FeO']  # KDMnFe(ol/l) from Davis et al. (2013)
    ol['NiO'] = np.exp(4272/T+0.01582*oxide_wt['SiO2']-2.7622)*ol['MgO'] # fitted by MPN+Hzb dataset (Eqn.3 in the paper)
    cpx['NiO'] = 0.24*1.08 * ol['NiO'] # average KdNi(cpx/ol) from Sobolev et al. (2005) Table S1, 1.08 is a factor to convert wtKdNi(cpx/ol) to cmKdNi(cpx/ol), observed from Walter 1998 data
    cpx['MnO'] = 0.85*0.98  # KdMn(cpx/l) modified after Le Roux et al. (2011) Table 3, 0.98 is a factor to convert wtKdMn(cpx/l) to cmKdMn(cpx/l), observed from Walter 1998 data
    kd = {'cpx':cpx, 'ol':ol, 'plg':plg}
    return kd



