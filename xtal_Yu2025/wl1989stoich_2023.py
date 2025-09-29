# switch the unit of oxides between wt% and cation mole fraction
# calculate elements either in oxides or components, e.g., K2O, KO0.5, KAlO2
# use stoichiometry, based on Weaver and Langmuir 1990
# originally written by Jocelyn Fuentes 2016
# modified by Mingzhen Yu 2021: add Ni and Mn in the system
# last modified: Sep 26, 2025

import numpy as np
import math     


## from wl1989stoich2021.py
mass = {'Si':28.0855, 'Ti':47.867, 'Al':26.9815, 'Fe':55.845, 'Mg':24.305, 
        'Ca':40.078, 'Na':22.98977,'O':15.999,'K':39.0987, 'P':30.973, 'Mn':54.938, 'Ni': 58.6934}

# convert oxide in wt% to cation mole fraction
def oxideToMolFracElement(oxides):
    mol = {}
    oxide_norm = {} # oxides in wt.% after being normalized to 100
    cationMolFrac = {}
    #Normalize oxides to 100wt%
    # tot = sum(oxides.values())
    for element in oxides:
        # oxide_norm[element] = (oxides[element])/tot
        oxide_norm[element] = oxides[element]
    #Get moles of each oxides
    mol['SiO2'] = oxide_norm['SiO2']/(mass['Si']+2*mass['O']) 
    mol['TiO2'] = oxide_norm['TiO2']/(mass['Ti']+2*mass['O']) 
    mol['AlO15'] = (oxide_norm['Al2O3']*2.)/(2.*mass['Al']+3.*mass['O']) 
    mol['CaO'] = oxide_norm['CaO']/(mass['Ca']+mass['O']) 
    mol['NaO5'] = (oxide_norm['Na2O']*2.)/(2.*mass['Na']+mass['O']) 
    mol['MgO'] = oxide_norm['MgO']/(mass['Mg']+mass['O']) 
    mol['FeO'] = oxide_norm['FeO']/(mass['Fe']+mass['O']) 
    mol['KO5'] = (oxide_norm['K2O']*2.)/(2.*mass['K']+mass['O'])
    mol['PO52']=(oxide_norm['P2O5']*2.)/(2.*mass['P']+5.*mass['O'])
    mol['MnO'] = oxide_norm['MnO']/(mass['Mn']+mass['O']) 
    mol['NiO'] = oxide_norm['NiO']/(mass['Ni']+mass['O']) 
    #mol['HO5'] = oxide_norm['H2O']/(mass['H']+0.5*mass['O'])
    #Calculate mole fraction of each element
    tot1 = sum(mol.values())
    cationMolFrac = {element: mol[element]/tot1 for element in mol}
    return cationMolFrac

# convert oxide in cation mole fraction to component      
def molFractoComponent(cationMolFrac):
    compCationMolFrac = {}
    #Components: CaAl2O4, NaAlO2, MgO, FeO, CaSiO3, TiO2, KAlO2, PO52, MnO, NiO
    compCationMolFrac['MgO'] = cationMolFrac['MgO']
    compCationMolFrac['FeO'] = cationMolFrac['FeO']
    compCationMolFrac['KAlO2'] = 2*cationMolFrac['KO5']
    compCationMolFrac['PO52'] = cationMolFrac['PO52']
    compCationMolFrac['MnO'] = cationMolFrac['MnO']
    compCationMolFrac['NaAlO2'] = 2.*cationMolFrac['NaO5']
    al_CaAl2O4 = (cationMolFrac['AlO15'] - cationMolFrac['NaO5'] - cationMolFrac['KO5'])
    compCationMolFrac['CaAl2O4'] = 1.5*al_CaAl2O4
    compCationMolFrac['CaSiO3'] = 2.*(cationMolFrac['CaO'] - (compCationMolFrac['CaAl2O4']/3.))
    compCationMolFrac['TiO2'] = cationMolFrac['TiO2']
    compCationMolFrac['NiO'] = cationMolFrac['NiO']
    #compCationMolFrac['SiO2'] = cationMolFrac['SiO2'] - (compCationMolFrac['CaSiO3']/2.)
    #compCationMolFrac['HO5'] = cationMolFrac['HO5']
    return compCationMolFrac

# convert oxide in wt% to component, combine the two functions above
def oxideToComponent(oxides):  # oxides is dictonary, e.g., SiO2: 50, MgO: 9, ......
    cationMolFrac = oxideToMolFracElement(oxides)
    components = molFractoComponent(cationMolFrac)
    return components

# convert componenet back to oxide in wt%
def cationFracToWeight(components):
    oxide = {}
    oxide_dict = {}
    oxide['Na2O'] = (components['NaAlO2']/2.)*(((2.*mass['Na'])+mass['O'])/2.)
    oxide['TiO2'] = components['TiO2']*(mass['Ti'] + (2*mass['O']))
    oxide['Al2O3'] = ((components['CaAl2O4']*(2./3.)) +(components['NaAlO2']/2.)+(components['KAlO2']/2.))*(((2.*mass['Al'])+(3.*mass['O']))/2.)
    oxide['FeO'] = components['FeO']*(mass['Fe']+mass['O'])
    oxide['MgO'] = components['MgO']*(mass['Mg']+mass['O'])
    oxide['CaO'] = ((components['CaAl2O4']/3.) + (components['CaSiO3']/2.))*(mass['Ca']+mass['O'])
    oxide['K2O'] = (components['KAlO2']/2.)*(2.*mass['K']+mass['O'])/2.
    oxide['P2O5'] = components['PO52']*(2.*mass['P']+5.*mass['O'])/2.
    oxide['MnO'] = components['MnO']*(mass['Mn']+mass['O'])
    oxide['NiO'] = components['NiO']*(mass['Ni']+mass['O'])
    oxide['SiO2'] = (1.-sum(components.values())+(components['CaSiO3']/2.))*(mass['Si']+(2*mass['O']))
    tot = (sum(oxide.values())-oxide['FeO']+oxide['FeO']/0.9)/100.
    oxide_dict = {element: oxide[element]/tot for element in oxide}
    return oxide_dict

