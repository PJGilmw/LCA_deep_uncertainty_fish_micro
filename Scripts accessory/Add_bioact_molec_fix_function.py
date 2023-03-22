# -*- coding: utf-8 -*-
"""
Created on Mon Jan  9 15:58:28 2023

@author: Pierre Jouannais, Department of Planning, DCEA, Aalborg University
pijo@plan.aau.dk

"""
'''
Script containing the function which calculate the compound content in the microalgal biomass for 1 data point in the raw output dataframe.



'''

import pandas as pd
import decimal
from random import *
import pstats
from itertools import *
from math import*
import csv
import copy
import numpy as np
import random

import datetime
from time import *




def biomasscompo(lipid_af_dw,
                 ash_dw,
                 water_after_drying,
                phospholipid_fraction,
                 elemental_contents):
    """Returns the macronutrient and elemental biomass composition for :
        the ash-free dry weight (ashfree dw) , total dry weight (dw),
        and once the biomass is harvested and dried, with remaining water (dbio).

    #Inputs:

        #lipid_af_dw: lipid content of the ash-free, dry biomass ; g.g ashfree dw -1
        #ash_dw: ash content in dry biomass ; g.g-1 dw
        #water_after_drying: remaining water in dried biomass ; g.g dbio -1
        #phospholipid_fraction: share of phospholipds among lipids ; g.g lipids-1
        #elemental_contents: table containing the elemental composition of 
                             each macronutrient ; g.g macronutrient-1

    #Outputs:

        #prot_af_dw:  protein content ash-free dw ; g.g ashfree dw -1
        #carb_af_dw: carbohydrate content ash-free dw ; g.g ashfree dw-1
        #C_af_dw: C content ash-free dw ; g.g ashfree dw-1
        #N_af_dw: N content ash-free dw ; g.g ashfree dw-1
        #P_af_dw: P content ash-free dw ; g.g ashfree dw-1
        #K_af_dw: K content ash-free dw ; g.g ashfree dw-1
        #Mg_af_dw: Mg content ash-free dw ; g.g ashfree dw-1
        #S_af_dw: S content ash-free dw ; g.g ashfree dw-1

        #lip_dw: lipid content  dw ; g.g dw -1
        #prot_dw: protein content dw ; g.g dw -1
        #carb_dw: carb content dw ; g.g dw -1

        #C_dw: C content dw ; g.g dw -1
        #N_dw: N content dw ; g.g dw -1
        #P_dw: P content dw ; g.g dw -1
        #K_dw: K content dw ; g.g dw -1
        #Mg_dw: Mg content dw ; g.g dw -1
        #S_dw: S content dw ; g.g dw -1

        #lip_dried_biomass: lipid content dried biomass ; g.g dbio -1
        #prot_dried_biomass: prot content dried biomass ;  g.g dbio -1
        #carb_dried_biomass: carb content dried biomass ; g.g dbio -1
        #ash_dried_biomass: ash content dried biomass ; g.g dbio -1
        """

    # ash-free dry weight

    prot_af_dw = (1 - lipid_af_dw)/(5/3)
    carb_af_dw = (2/3)*prot_af_dw

    # dry weight

    lip_dw = lipid_af_dw*(1 - ash_dw)
    prot_dw = prot_af_dw*(1 - ash_dw)
    carb_dw = carb_af_dw*(1 - ash_dw)

    # dried biomass (with remaining water)

    lip_dried_biomass = lip_dw*(1 - water_after_drying)
    prot_dried_biomass = prot_dw*(1 - water_after_drying)
    carb_dried_biomass = carb_dw*(1 - water_after_drying)
    ash_dried_biomass = ash_dw*(1 - water_after_drying)


    # Elemental composition

    C_lip = lipid_af_dw*(elemental_contents.iloc[1, 0]
                         * (1 -phospholipid_fraction)
                         + elemental_contents.iloc[2, 0]*phospholipid_fraction)

    C_prot = prot_af_dw * elemental_contents.iloc[0, 0]
    C_carb = carb_af_dw * elemental_contents.iloc[3, 0]

    C_af_dw = C_lip + C_prot + C_carb
    C_dw = C_af_dw*(1 - ash_dw)

    N_lip = lipid_af_dw*(elemental_contents.iloc[1, 1]*(
        1-phospholipid_fraction)+elemental_contents.iloc[2, 1]*phospholipid_fraction)

    N_prot = prot_af_dw * elemental_contents.iloc[0, 1]
    N_carb = carb_af_dw * elemental_contents.iloc[3, 1]

    N_af_dw = N_lip + N_prot + N_carb
    N_dw = N_af_dw*(1 - ash_dw)

    P_lip = lipid_af_dw*(elemental_contents.iloc[1, 2]*(
        1-phospholipid_fraction)+elemental_contents.iloc[2, 2]*phospholipid_fraction)

    P_prot = prot_af_dw * elemental_contents.iloc[0, 2]
    P_carb = carb_af_dw * elemental_contents.iloc[3, 2]

    P_af_dw = P_lip + P_prot + P_carb
    P_dw = P_af_dw*(1 - ash_dw)

    # ratioC_N = C_af_dw / N_af_dw
    # ratioN_P = N_af_dw / P_af_dw


    # Other components according to Chlorella's observed ratios

    K_af_dw = N_af_dw*0.18
    Mg_af_dw = N_af_dw*0.08
    S_af_dw = N_af_dw*0.048

    K_dw = N_dw*0.18
    Mg_dw = N_dw*0.08
    S_dw = N_dw*0.048

    return [prot_af_dw,
            carb_af_dw,
            lip_dw,
            prot_dw,
            carb_dw,
            lip_dried_biomass,
            prot_dried_biomass,
            carb_dried_biomass,
            ash_dried_biomass,
            C_af_dw,
            N_af_dw,
            P_af_dw,
            K_af_dw,
            Mg_af_dw,
            S_af_dw,
            C_dw,
            N_dw,
            P_dw,
            K_dw,
            Mg_dw,
            S_dw]




    
def add_bioact_molec_dbio(datafr,elemental_contents):
    
        
    
    water_after_drying = 0.05 # unique
    new_column=[]
    for row_index in range(datafr.shape[0]):
        
        lipid_af_dw = datafr.loc[row_index,"lipid_af_dw"]
        ash_dw = datafr.loc[row_index,"ash_dw"]
        phospholipid_fraction =  datafr.loc[row_index,"phospholipid_fraction"]
        
        bioact_fraction_molec = datafr.loc[row_index,"bioact_fraction_molec"]
        
        random_bioclass =  datafr.loc[row_index,"random_bioclass"]

        if random_bioclass < (1/3):
            #print('okC1')
            biochemicalclass='lip'
            
        elif (1/3)<random_bioclass<(2/3)  :
            #print('okC2')
            biochemicalclass='prot'
            
        else:
            #print('okC2')
            biochemicalclass='carb'
    
        
        
        
        
        
        biomass_composition = biomasscompo(
            lipid_af_dw,
            ash_dw,
            water_after_drying,
            phospholipid_fraction,
            elemental_contents)
        
        
        
        # After harvesting and drying  (dbio)
        lip_dbio = biomass_composition[5]
        prot_dbio = biomass_composition[6]
        carb_dbio = biomass_composition[7]
        ash_dbio = biomass_composition[8]
        
        # Elementary composition
        
        
        # Calculating the absolute bioactive molecule content in the dried biomass
        if biochemicalclass == 'lip':
            bioact_molec_dbio = bioact_fraction_molec * lip_dbio
        
        elif biochemicalclass == 'carb':
            bioact_molec_dbio = bioact_fraction_molec * carb_dbio
        
        elif biochemicalclass == 'prot':
            bioact_molec_dbio = bioact_fraction_molec * prot_dbio
        
        
        
        new_column.append(bioact_molec_dbio)
    
    datafr_fixed = datafr.copy()
    datafr_fixed["bioact_molec_dbio"] = new_column
        
    
    return datafr_fixed

