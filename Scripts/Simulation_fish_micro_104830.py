# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 09:38:46 2022

@author: Pierre Jouannais, Department of Planning, DCEA, Aalborg University
pijo@plan.aau.dk

"""
'''
Script which produces the raw output table  with 1 data point being 1 LCA for Trout with and without compound.
It calls functions to perform the Montecarlo iterations in the foreground and combines with the background iterations imported as a pickle file.

This current script was parameterized to simulate 104830 data points (size of the background MC sample) but will adapt any size.

Requires instances with multiple CPUS
It was run on a remote server with 64 CPUs and took some days.


-1) Define parameters/factors (l.100))

-2) Initialize (l.645)

-3) Function for fast parallel montecarlo iterations in the background (l.487)

-3) RUN  (l.1369)



'''


import datetime
from time import *
import requests
import pickle
import cProfile
from scipy.integrate import odeint

import os
import sys


import pandas as pd
import decimal
from random import *
import pstats
from itertools import *
from math import*
import csv
import copy
import numpy as np
import seaborn as sns
import random


import bw2data
import bw2io
from bw2data.parameters import *
import brightway2 as bw


from SALib.test_functions import Ishigami
import math
from SALib.sample import saltelli
from SALib.sample import fast_sampler
from SALib.analyze import sobol
from SALib.analyze import fast
from SALib.sample import latin
from SALib.analyze import delta

import SALib

import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits import mplot3d

import shapefile as shp
import geopandas as gpd
import pysal as ps
from shapely.geometry import Polygon, mapping, Point


#import Cultivation_simul_Night_Harvest_1_2nd_faster_noreinjection as cultsimul
import Functions_for_physical_and_biological_calculations_3rd as functions
import Main_functions_Fish_micro as mainfunc

import Retrieving_solar_and_climatic_data_2nd as solardata

import Map_3rd as geo

import ray


import pert


import scipy
import numpy as np
from scipy import linalg

from ast import literal_eval



from itertools import compress



"""1) Define parameters/factors """

#  Parameters/ Factors dictionnaries

# Description of the parameters are given in the appendix.
# Here values can be changed for parameters with unique values (no distributions)
# Values with distributions will be overwritten by the distribution dictionnaries

#Biological parameters


Biodict = {'rhoalgae': 1070, # kg.m-3
           'lipid_af_dw': 0.3,   # .
           'ash_dw': 0.05,  # .
           'MJ_kglip': 36.3,  # MJ.kg-1
           'MJ_kgcarb': 17.3,  # MJ.kg-1
           'MJ_kgprot': 23.9,  # MJ.kg-1
           'PAR': 0.45,  # .
           'losspigmentantenna': 0.21,  # .
           'quantumyield': 8,  # molphotons.moloxygen-1
           'lossvoltagejump': 1-1267/1384,  # .
           'losstoATPNADPH': 1-590/1267,  # .
           'losstohexose': 1-469/590,  # .
           'lossrespiration': 0.20,  # .
           'bioact_fraction_molec': 0.1,  # .
           'random_no3': 1,   # .
           'Topt': 25,  # °C
           'T_plateau': 10,  # °C
           'dcell': 5*10**-6,  # m
           'incorporation_rate': 0.10,  # .
           'nutrient_utilisation': 0.9, # .
           'co2_utilisation': 0.85,  # .
           'phospholipid_fraction': 0.20,  # .

            'random_market_subst': 0.1,  # .
            
            # A value of 1 for Thermoregulation at night
            
            'random_night_monitoring': 0,
            
                         
            'random_bioclass': 0, 
            
            'fertiliser_substitutability_AD': 0.8,
                                 
            'carbon_degradibilidy_AD':0.7
           }

# Physical parameters

Physicdict = {'hconv': 6.6189,  # W.m−2 .K−1
              'Cp': 4.186,  # kJ.(kg.K)-1
              'rhomedium': 1000,  # kg.m-3
              'rhowater': 1000,  # kg.m-3
              'Cw': 2256,     # kJ.kg-1
              'CH4_LHV':50} # MJ.m-3



Locationdict = {'lat': 37.189, #°
                'long': -3.572, #°
                'azimuthfrontal': 90} # °


#Techno-operational parameters

Tech_opdict = {'height': 1.5,  # m
                   'tubediameter': 0.03,  # m
                   'gapbetweentubes': 0.01,  # m
                   'horizontaldistance': 0.2,  # m
                   'length_of_PBRunit': 30,    # m
                   'width_of_PBR_unit': 30,   # m
                   'biomassconcentration': 1.4,  # kg.m-3
                   'flowrate': 0.38,     # m.s-1
                   'centrifugation_efficiency': 0.98,  # .
                   'pumpefficiency': 0.90,  # .
                   'slurry_concentration': 0.15,  # .
                   'water_after_drying': 0.05,  # gwater.g dbio-1
                   'recyclingrateaftercentrifuge': 0.3, # gwater.g dbio-1
                   'rhosuspension': 1015,  # kg.m-3 
                   'roughness': 0.0000015,  # m
                   'cleaningvolumeVSfacilityvolume': 4,  # .
                   'concentration_hypo': 2*10**-3,  # kg.m-3
                   'concentration_hydro': 30,  # kg.m-3
                   'boilerefficiency': 0.75,  # .
                   'glass_life_expectancy': 50, # years
                   

                   'extraction': 'yes',
                   
                   'heat_pump': 'yes',
                   
                   'COP': 3,
                   
                   'heat_input_per_kg_biomass_AD':0.68, # kWh.kg dw -1 
                                 
                   "elec_input_per_kg_biomass_AD":0.108, # kWh.kg dw -1
                   
                   "fraction_maxyield":0.3,
                   
                   "depth_well": 10,
                   
                   "hconv": 9


                   }

# Fish farm 



Fish_farm_and_compound_effect_dict= {'HAIT_loss_lev': 0,   
                              
                                 'HAIT_loss_red': 0,
                                 
                                 'FRIT_loss_lev': 0,
                                 
                                 'FRIT_loss_red': 0,
                                 
                                 'GOIT1_loss_lev': 0, 
                                 
                                 'GOIT1_loss_red': 0, 
                                 
                                 'GOIT2_loss_red': 0, 
                                 
                                 'GOIT2_loss_lev': 0,
                                 
                                 'GOIT1bis_loss_lev': 0, 
                                 
                                 'GOIT1bis_loss_red': 0, 
                                                                  
                                 'GODK_loss_lev': 0, 
                                 
                                 'GODK_loss_red': 0, 
                                                                  
                                 'SFDK1_loss_lev': 0, 
                                 
                                 'SFDK2_loss_lev': 0, 
                                 
                                 'SFDK1_loss_red': 0, 
                                 
                                 'SFDK2_loss_red': 0, 
                                 
                                 'FRIT_micro_dose': 0, 

                                 'GODK_micro_dose': 0,  
                                                                  
                                 'GOIT2_micro_dose': 0,
                                 
                                 'GOIT1_micro_dose': 0,
                   
                                 'GOIT1bis_micro_dose': 0,
                                 
                                 
                                 'HAIT_para_loss_lev': 0,   
                              
                                 
                                 'FRIT_para_loss_lev': 0,
                                 
                                 
                                 'GOIT1_para_loss_lev': 0, 
                                 
                                                                  
                                 'GOIT2_para_loss_lev': 0,
                                 

                                 
                                 'HAIT_micro_dose': 0, # kWh.kg dw -1
                                 
                                 'SFDK2_micro_dose': 0, # kWh.kg dw -1
                                 
                                 'SFDK1_micro_dose': 0,
                                 
                                 
                                 "HAIT_FCR_red_ratio_frac":0,
                                 "FRIT_FCR_red_ratio_frac":0,
                                 "GOIT1_FCR_red_ratio_frac":0,
                                 "GOIT2_FCR_red_ratio_frac":0,
                                 "GOIT1bis_FCR_red_ratio_frac":0,
                                 "GODK_FCR_red_ratio_frac":0,
                                 "SFDK1_FCR_red_ratio_frac":0,
                                 "SFDK2_FCR_red_ratio_frac":0,
                                 "ratio_CO2_CH4_biogas_fish":0.465, 
                                 "CH4_volume_biogas_fish":200, # ml
                                 "P_in_fish": 0.00415,  # kg
                                 "N_in_fish": 0.022, # kg 
                                  "K_in_fish":0.0046,
                                  "Mg_in_fish":0.00038,
                                  "water_in_fish": 0.7,
                                  "gvs_gts_in_fish": 0.89,
                                  "fraction_non_ingested" : 0.05,
                                  "excretion_N_removal_efficiency" :0.742095258,
                                  "excretion_P_removal_efficiency" : 0.656726614,
                                  "solid_filter_efficiency": 0.3485,
                                 
    
                                  "CH4_volume_sludge_and_manure": 350,
                                  "share_fish_sludge_in_substrate": 0.33,
                                  "ratio_CO2_CH4_biogas_sludge": 0.465,
                                  "N_manure": 0.05,
                                  "P_manure": 0.01,
                                  "fertilizer_substi_digest_N":0.9,
                                  "fertilizer_substi_digest_P":0.5,
                                  "fertilizer_substi_digest_K": 0.5,
                                  "fertilizer_substi_digest_Mg": 0.5,
                                  

                                  "fertilizer_substi_manure_field_N": 0.54,
                                  "fertilizer_substi_manure_field_P": 0.3}



Production_mix_dict= {'mix_size':5,   
                              
                      'mode_lat': 49.55}





# Parameters distribution directories 

# Exact same parameters dictionnaries the normal ones but instead of one value,
# each parameter is assigned a list containing  :
   # [Distribution,min,max,mode or mean,sd]

# Distribution :
#   - 'unique' if no distribution. The value indicated indicated in
#     the normal dictionnary will be considered.
#   - 'unif' for uniform, uses min and max
#   - 'triang, uses mim max and mode with mode as a fracion of max-min

#

# Biological parameters



Biodict_distributions = {'rhoalgae': ['unif', [0, 1020, 1250, 0, 0]], 
                           
                         'lipid_af_dw': ['triang', [0, 0.1, 0.7, 0.2, 0]],    # see dist_lip_williams
                         
                         'ash_dw': ['unif', [0, 0.01, 0.10, 0.10, 0]],
                         
                         'MJ_kglip': ['unique', [36.3, 0, 0, 0, 0]],
                       
                         'MJ_kgcarb': ['unique', [17.3, 0, 0, 0, 0]],
                        
                         'MJ_kgprot': ['unique', [23.9, 0, 0, 0, 0]],
                         
                         'PAR': ['unique', [0.45, 0, 0, 0, 0]],
                        
                         'losspigmentantenna': ['unif', [0, 0.16, 0.24, 0, 0]],
                        
                         'quantumyield': ['unif', [0, 8, 11, 0, 0]],
                         
                         'lossvoltagejump': ['unique', [1-1267/1384, 0, 0, 0, 0]],
                       
                         'losstoATPNADPH': ['unique', [1-590/1267, 0, 0, 0, 0]],
                 
                         'losstohexose': ['unique', [1-469/590, 0, 0, 0, 0]],
                      
                         'lossrespiration': ['unif', [0, 0.16, 0.24, 0, 0]],
                    
                         'bioact_fraction_molec': ['unif', [0, 0.01, 1, 0, 0]],
                         
                         'random_no3': ['unif', [0, 0, 1, 0, 0]],
                         
                         'Topt': ['unif', [0, 15, 35, 25, 0]],
                         
                         'T_plateau': ['unif', [0, 5, 10, 0, 0]],
                      
                         'dcell': ['unif', [0, 0.000001, 0.00005, 0.000005, 0]],
                
                         'incorporation_rate': ['unique', [0, 0.01, 0.15, 0, 0]],

                         'nutrient_utilisation': ['unif', [0, 0.75, 0.90, 0, 0]],
             
                         'co2_utilisation': ['triang', [0, 0.5, 0.95, 0.9, 0]],
            
                         'phospholipid_fraction': ['unif', [0, 0.1, 0.6, 0, 0]],

                         'random_market_subst': ['unif', [0, 0, 1, 0, 0]],
                         
                         'random_night_monitoring': ['unif', [0, 0, 1, 0, 0]],
                         
                         'random_bioclass': ['unif', [0, 0, 1, 0, 0]], 

                         'fertiliser_substitutability_AD': ['unif', [0.8, 0.6, 0.9, 0, 0]],
                                 
                         'carbon_degradibilidy_AD': ['unif', [0.7, 0.5, 0.6, 0, 0]]
                         }


# Physical parameters

Physicdict_distributions = {
    
    'Cp': ['unique', [4.186, 0, 0, 0, 0]],
     
    'rhomedium': ['unique', [1000, 0, 0, 0, 0]],
    
    'rhowater': ['unique', [1000, 0, 0, 0, 0]],
    
    'Cw': ['unique', [2256, 0, 0, 0, 0]], #kJ.kg-1
    
    'CH4_LHV':['unique', [50, 0, 0, 0, 0]]} # MJ.kg-1
    






Locationdict_distributions = {'lat': ['unique', [43.695, 0, 0, 0, 0]],
                              
                              'long': ['unique', [1.922, 0, 0, 0, 0]],
                                                            
                              
                              'azimuthfrontal': ['unique', [90, 0, 0, 0, 0]]}



# Mode value of a PERT or triangular distribution

Tech_opdict_distributions = {'height': ['unique', [1.5, 0, 0, 0, 0]],   
                              
                                 'tubediameter': ['pert_nested', [0, 0.03, 0.1, 0., 0]],
                                 
                                 'gapbetweentubes': ['pert_nested', [0, 0.01, 0.1, 0, 0]],
                                 
                                 'horizontaldistance': ['pert_nested', [0, 0.2, 0.8, 0, 0]],
                                 
                                 'length_of_PBRunit': ['unique', [20, 10, 35, 0, 0]], 
                                 
                                 'width_of_PBR_unit': ['unique', [20, 3, 10, 0, 0]], 
                                 
                                 'biomassconcentration': ['pert_nested', [0, 1, 7, 0, 0]],  #lower,upper, mode
                                 
                                 'flowrate': ['pert_nested', [0, 0.2, 1, 0, 0]],
                                 
                                 'centrifugation_efficiency': ['unique', [0.98, 0, 0, 0, 0]], 
                                 
                                 'pumpefficiency': ['unique', [0.9, 0, 0, 0, 0]], 
                                                                  
                                 'slurry_concentration': ['unique', [0, 0.10, 0.20, 0, 0]],  
                                 
                                 'water_after_drying': ['unique', [0, 0.02, 0.07, 0, 0]], 
                                 
                                 'recyclingrateaftercentrifuge': ['unique', [0, 0.2, 0.4, 0, 0]],
                                 
                                 'rhosuspension': ['pert_nested', [0, 1000, 1200, 0, 0]], 
                                 
                                 'roughness': ['unique', [0.0000015, 0, 0, 0, 0]], 
                                 
                                 'cleaningvolumeVSfacilityvolume': ['unique', [0, 2, 6, 0, 0]], 
                                 
                                 'concentration_hypo': ['unique', [2*10**-3, 0, 0, 0, 0]], 
                                 
                                 'concentration_hydro': ['unique', [30, 0, 0, 0, 0]],
                                 
                                 'glass_life_expectancy': ['unique', [50, 0, 0, 0, 0]], 

                                 'boilerefficiency': ['unique', [0.75, 0, 0, 0, 0]],  
                                                                  
                                 'extraction': ['binary', ['yes', 0, 0, 0, 0]],
                                 
                                 'heat_pump': ['binary', ['yes', 0, 0, 0, 0]],
                   
                                 'COP': ['unique', [3, 2, 6, 0, 0]],
                                 
                                 'heat_input_per_kg_biomass_AD': ['unique', [0.68, 0, 0, 0, 0]], # kWh.kg dw -1
                                 
                                 'elec_input_per_kg_biomass_AD': ['unique', [0.108, 0, 0, 0, 0]],
                                 
                                 "fraction_maxyield": ['unique', [0.3, 0, 0, 0, 0]],
                                 
                                 'depth_well': ['pert_nested', [0, 5, 25, 0, 0]],
                                 'hconv': ['pert_nested', [0, 5, 10, 0, 0]]


                                 
                                 }
                                 

nested_list_techno_op=[]

to_be_deleted_op=[]
for param in Tech_opdict_distributions:
    

    if Tech_opdict_distributions[param][0] == 'pert_nested':
            nested_list_techno_op.append([param,Tech_opdict_distributions[param][1]])
                


# Fish farm 




Fish_farm_and_compound_effect_dict_distributions= {'HAIT_loss_lev': ['unif', [0, 0, 0.15, 0, 0]],   
                              
                                 'HAIT_loss_red': ['unif', [0, 0, 1, 0, 0]],
                                 
                                 'FRIT_loss_lev': ['unif', [0, 0, 0.15, 0, 0]],
                                 
                                 'FRIT_loss_red': ['unif', [0, 0, 1, 0, 0]],
                                 
                                 'GOIT1_loss_lev': ['unique', [0, 0, 0, 0, 0]],  # Never modified
                                 
                                 'GOIT1_loss_red': ['unique', [0, 0, 0, 0, 0]],  # Never modified
                                 
                                 'GOIT2_loss_red': ['unique', [0, 0, 0, 0, 0]],  # Never modified
                                 
                                 'GOIT2_loss_lev': ['unique', [0, 0, 0, 0, 0]],  # Never modified
                                 
                                 'GOIT1bis_loss_lev': ['unif', [0, 0, 0.15, 0, 0]],   # Never modified
                                 
                                 'GOIT1bis_loss_red': ['unif', [0, 0, 1, 0, 0]],  # Never modified
                                                                  
                                 'GODK_loss_lev': ['unif', [0, 0, 0.15, 0, 0]], 
                                 
                                 'GODK_loss_red': ['unif', [0, 0, 1, 0, 0]], 
                                 
                                 'SFDK1_loss_lev': ['unif', [0, 0, 0.15, 0, 0]], 
                                 
                                 'SFDK2_loss_lev': ['unif', [0, 0, 0.15, 0, 0]], 
                                 
                                 'SFDK1_loss_red': ['unif', [0, 0, 1, 0, 0]], 
                                 
                                 'SFDK2_loss_red': ['unique', [0, 0, 1, 0, 0]], 
                                 
                                 'FRIT_micro_dose': ['unif', [0, 0.00000017, 0.0047, 0, 0]], 

                                 'GODK_micro_dose': ['unif', [0, 0.00000017, 0.0047, 0, 0]],    
                                                                  
                                 'GOIT2_micro_dose': ['unique', [0, 0, 0, 0, 0]], # Always 0
                                 
                                 'GOIT1_micro_dose': ['unique', [0, 0, 0, 0, 0]],  # Always 0
                   
                                 'GOIT1bis_micro_dose': ['unif', [0, 0.00000017, 0.0047, 0, 0]], # Always 0
                                 
                                 
                                 'HAIT_para_loss_lev': ['unique', [0, 0, 0, 0, 0]],   
                              
                                 
                                 'FRIT_para_loss_lev': ['unique', [0, 0, 0, 0, 0]],
                                 
                                 
                                 'GOIT1_para_loss_lev': ['unique', [0, 0, 0, 0, 0]], 
                                 
                                                                  
                                 'GOIT2_para_loss_lev': ['unique', [0, 0, 0, 0, 0]],
                                 
                                 
                                 'HAIT_micro_dose': ['unif', [0, 0.00000017, 0.0047, 0, 0]], 
                                 
                                 'SFDK2_micro_dose': ['unif', [0, 0.00000017, 0.0047, 0, 0]], 
                                 
                                 'SFDK1_micro_dose': ['unif', [0, 0.00000017, 0.0047, 0, 0]], 
                                 
               
                                 "HAIT_FCR_red_ratio_frac":['unif', [0, 0, 1, 0, 0]],
                                 "FRIT_FCR_red_ratio_frac":['unif', [0, 0, 1, 0, 0]],
                                 "GOIT1_FCR_red_ratio_frac":['unique', [0,  0, 0, 0, 0]], # Never modified
                                 "GOIT2_FCR_red_ratio_frac":['unique', [0,  0, 0, 0, 0]],  # Never modified
                                 "GOIT1bis_FCR_red_ratio_frac":['unif', [0,  0, 1, 0, 0]],  # Never modified
                                 "GODK_FCR_red_ratio_frac":['unif', [0, 0, 1, 0, 0]],
                                 "SFDK1_FCR_red_ratio_frac":['unif', [0, 0, 1, 0, 0]],
                                 "SFDK2_FCR_red_ratio_frac":['unif', [0, 0, 1, 0, 0]],
                                 
                                 "ratio_CO2_CH4_biogas_fish":['unique', [0.465, 0.6, 0.9, 0, 0]], 
                                 "CH4_volume_biogas_fish": ['triang', [0, 390, 920, 550, 0]]  , # ml
                                 "P_in_fish":['unique', [0.00415, 0.0026, 0.005513, 0, 0]] ,
                                 "N_in_fish":['unique', [0.021627, 0, 0, 0, 0]]  ,
                                  "K_in_fish": ['unique', [0.004618313, 0.002954, 0.005513, 0, 0]]  ,
                                  "Mg_in_fish": ['unique', [0.00038425, 0.000301, 0.000431, 0, 0]],
                                  "water_in_fish": ['unique', [0.7, 0, 0, 0, 0]],
                                  "gvs_gts_in_fish": ['unique', [0.89, 0, 0, 0, 0]],
                                  "fraction_non_ingested" : ['unif', [0.05, 0.00, 0.05, 0, 0]],
                                  "excretion_N_removal_efficiency" :['unique', [0, 0, 0, 0, 0]],
                                  "excretion_P_removal_efficiency" : ['unique', [0, 0, 0, 0, 0]],
                                  "solid_filter_efficiency": ['unique', [0.3485, 0, 0, 0, 0]],
                                  "CH4_volume_sludge_and_manure": ['unif', [350, 300, 400, 0, 0]],
                                  "share_fish_sludge_in_substrate": ['unif', [0.33, 0.30, 0.40, 0, 0]],
                                  "ratio_CO2_CH4_biogas_sludge": ['unique', [0.465, 0, 0, 0, 0]],
                                  "N_manure": ['unique', [0.05, 0, 0, 0, 0]],
                                  "P_manure": ['unique', [0.01, 0, 0, 0, 0]],
                                  "fertilizer_substi_digest_N": ['unif', [0.9, 0.8, 1, 0, 0]],
                                  "fertilizer_substi_digest_P": ['unif', [0.6, 0.2, 0.8, 0, 0]],
                                  "fertilizer_substi_digest_K": ['unif', [0.6, 0.2, 0.8, 0, 0]],
                                  "fertilizer_substi_digest_Mg": ['unif', [0.6, 0.2, 0.8, 0, 0]],
                                  
                                  "fertilizer_substi_manure_field_N": ['unif', [0.5, 0.53, 0.55, 0, 0]],
                                  "fertilizer_substi_manure_field_P": ['unif', [0.28, 0.20, 0.385, 0, 0]]}

   
    

Production_mix_dict_distributions= {'mix_size': ['unif', [0, 1, 25, 0, 0]],   
                              
                                 'mode_lat': ['pert_nested', [0, 37, 60, 0, 0]]}








# Accessory functions

def createFolder(directory):
    '''Creates a folder/directory with the path given as input'''
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)


def export_pickle(var, name_var):
    '''Saves a pickle in the working driectory and
    saves the object in input across python sessions'''

    with open(name_var+'.pkl', 'wb') as pickle_file:
        pickle.dump(var, pickle_file, protocol=pickle.HIGHEST_PROTOCOL)

def export_pickle_2(var, name_var, namefolder_in_root):
    '''Saves a pickle in the working directory and
    saves the object in input across python sessions'''

    path_object = "../"+namefolder_in_root+"/"+name_var+".pkl"
    with open(path_object, 'wb') as pickle_file:
        pickle.dump(var, pickle_file, protocol=pickle.HIGHEST_PROTOCOL)



def importpickle(path):
    with open(path, 'rb') as pickle_load:
        obj = pickle.load(pickle_load)
    return obj    



""" 1 ) Initialize """


# Fish feed composition

fishfeed_table_withNP = pd.read_csv("../Data/Ingredient_composition_withNP.csv", sep=";",
                             header=0, encoding='unicode_escape', engine='python')

# Cleaning

fishfeed_table_withNP = fishfeed_table_withNP[0:8][[
    'Ingredient', 'kg.kg feed-1', 'Lipid', 'Protein', 'Carb', 'Ash', 'Water',"N","P"]]



# Calculaing incumbent fish feed biochmical profile

biochem_profile_feed_incumbent = functions.biochemprofile(fishfeed_table_withNP)


N_P_profile_feed_incumbent = [sum([fishfeed_table_withNP["N"][a] * fishfeed_table_withNP['kg.kg feed-1'][a] for a in range(len(fishfeed_table_withNP["N"]))]),
                              sum([fishfeed_table_withNP["P"][a] * fishfeed_table_withNP['kg.kg feed-1'][a] for a in range(len(fishfeed_table_withNP["P"]))])]
    
ingredient_profile_incumbent = list(fishfeed_table_withNP['kg.kg feed-1']) 



# Elemental composition of macronutrients

elemental_contents = pd.read_csv("../Data/elemental_contents.csv",
                                 sep=";",
                                 header=0,
                                 encoding='unicode_escape',
                                 engine='python')

# Cleaning

elemental_contents = elemental_contents.iloc[:, 1:]




digestibility_list = [0.95,0.92,0.71,0.5,0.5] # lip, prot, carb, ash, phosphorus



# Import foreground Technosphere Fish Production


Techno_Matrix_Fish_loaded = pd.read_excel('../Data/Technosphere_Fish_micro_17_03.xlsx',header=None)


Techno_Matrix_Fish_with_names =  np.array(Techno_Matrix_Fish_loaded)

Techno_Matrix_Fish_activities = Techno_Matrix_Fish_with_names[0,1:]

Techno_Matrix_Fish_products = Techno_Matrix_Fish_with_names[1:,0]

# Only amounts, no names
Techno_Matrix_Fish = Techno_Matrix_Fish_with_names[1:,1:]

Techno_Matrix_Fish.shape
if Techno_Matrix_Fish.shape[0]!=Techno_Matrix_Fish.shape[1]:
    
     sys.exit("The Technosphere matrix is not square") 



# Keep the names and the positions of the activities which are ONLY connected to background activites
# This means excluding activites starting with "FG_"

filter_names_ok_activities_fish = [name[0:3] !="FG_" for name in Techno_Matrix_Fish_activities]
filter_names_exclude = [name[0:3] =="FG_" for name in Techno_Matrix_Fish_activities]


activities_fish_background = list(compress(Techno_Matrix_Fish_activities, filter_names_ok_activities_fish))

activities_fish_foreground = list(compress(Techno_Matrix_Fish_activities, filter_names_exclude))

positions_background = [np.where(Techno_Matrix_Fish_activities == act)[0][0] for act in activities_fish_background]

positions_foreground = [np.where(Techno_Matrix_Fish_activities == act)[0][0] for act in activities_fish_foreground]


index_growth_stages= positions_foreground[:-2]

names_growth_stages_no_solid_filter_laguna = ['FG_Seafarm_DK_1','FG_Seafarm_DK_2']


index_growth_stages_no_filter_laguna = [np.where(Techno_Matrix_Fish_activities == act)[0][0] for act in names_growth_stages_no_solid_filter_laguna]

index_roe = np.where(Techno_Matrix_Fish_activities == "Roe_prod")[0][0]

index_feed = np.where(Techno_Matrix_Fish_activities == "Fish_feed")[0][0]  


index_dead_fish = np.where(Techno_Matrix_Fish_activities == "FG_Dead_fish_management_supply")[0][0]


index_Nemissions = np.where(Techno_Matrix_Fish_activities == "N_emissions_background")[0][0]

index_Pemissions = np.where(Techno_Matrix_Fish_activities == "P_emissions_background")[0][0]

index_sludge = np.where(Techno_Matrix_Fish_activities == "FG_Sludge management")[0][0]

index_growing_DK = np.where(Techno_Matrix_Fish_activities == "FG_Growing_out_DK")[0][0]

index_300g = np.where(Techno_Matrix_Fish_activities == "FG_Growing_out_IT_2_para")[0][0]

index_micro_compound = np.where(Techno_Matrix_Fish_activities == "Micro_compound_prod")[0][0]


index_biogas_updgrade = np.where(Techno_Matrix_Fish_activities == "Heat market substitution Fish")[0][0]
index_N_substitution = np.where(Techno_Matrix_Fish_activities == "N source production fish")[0][0]
index_P_substitution = np.where(Techno_Matrix_Fish_activities == "P source production fish")[0][0]
index_heat_substitution = np.where(Techno_Matrix_Fish_activities == "Heat market substitution Fish")[0][0]

index_drug1=np.where(Techno_Matrix_Fish_activities == "Drug_1_prod")[0][0]
index_drug2=np.where(Techno_Matrix_Fish_activities == "Drug_2_prod")[0][0]
index_drug3=np.where(Techno_Matrix_Fish_activities == "Drug_3_prod")[0][0]
index_drug4=np.where(Techno_Matrix_Fish_activities == "Drug_4_prod")[0][0]


# Prepare demand vector
demand_vector= np.zeros(Techno_Matrix_Fish.shape[0])

index_FU = np.where(Techno_Matrix_Fish_activities == "FG_Seafarm_DK_2")[0][0]

demand_vector[index_FU] = 1





                             
                             



def initial(Techno_Matrix_Fish,
            demand_vector,
            index_300g,
            index_growing_DK,
            index_roe,
            index_FU,
            index_feed,
            index_dead_fish):
    """Function to calculate the current economic indicators,
    before modification"""    
    
    # Before modification

     
    HAIT_loss_lev=0  
                                     
    FRIT_loss_lev=0
                                     
                                     
    GOIT1_loss_lev=0
                                     
                                                             
    GOIT2_loss_lev=0
                                     
    GOIT1bis_loss_lev=0
                                     
    
    GODK_loss_lev=0
                                     
                                                                      
    SFDK1_loss_lev= 0
                                     
    SFDK2_loss_lev= 0 
                                     
    SFDK1_loss_red= 0 
    SFDK2_loss_red=0
    FRIT_loss_red=0
    GODK_loss_red=0
    GOIT2_loss_red =0   
    HAIT_loss_red=0
    GOIT1bis_loss_red=0
            
    HAIT_para_loss_lev=0
    FRIT_para_loss_lev=0
    GOIT1_para_loss_lev=0
    GOIT2_para_loss_lev=0
         
    
                 
    FRIT_micro_dose= 0 
    
    GODK_micro_dose= 0  
                                                                      
    GOIT2_micro_dose= 0
                                     
    GOIT1_micro_dose= 0
                       
    GOIT1bis_micro_dose= 0
                              
    HAIT_micro_dose= 0 
                                     
    SFDK2_micro_dose= 0 
                                     
    SFDK1_micro_dose= 0
            
    
    HAIT_FCR_red_ratio=1
    FRIT_FCR_red_ratio=1
    GOIT1_FCR_red_ratio=1
    GOIT2_FCR_red_ratio=1
    GOIT1bis_FCR_red_ratio=1
    GODK_FCR_red_ratio=1
    SFDK1_FCR_red_ratio=1
    SFDK2_FCR_red_ratio=1
    
    
    ratio_CO2_CH4_biogas_fish = 1
    CH4_volume_biogas_fish = 1
    
    
    # Just initialize with dummy values for these parameters here


    P_in_fish= 0.00415  # kg
    N_in_fish= 0.021627 # kg
    K_in_fish=0.004618313
    Mg_in_fish=0.00038425
    water_in_fish= 0.7
    gvs_gts_in_fish= 0.89
    
    gvs_gts_in_fish= 0.89
    
    
    
    
    biogenic_CO2_emitted = 10
    amount_of_combustion_natural_gas = 23
    m3_natural_gas_substituted = 9
    amount_of_upgrading_act = 2

    excretion_N_removal_efficiency = 1
    excretion_P_removal_efficiency = 1
    solid_filter_efficiency =1

    fertilizer_substi_digest_P=1
    fertilizer_substi_digest_N=1
    fertilizer_substi_digest_K=1
    fertilizer_substi_digest_Mg=1
    
    MJ_substituted_per_kilo_deadfish = 1
    
    # Initialize the numerical Tehcnosphere
    Techno_fish_num=np.zeros((Techno_Matrix_Fish.shape[0],Techno_Matrix_Fish.shape[0]))
    
    for i in range(Techno_fish_num.shape[0]):
    
        for j in range(Techno_fish_num.shape[1]):
            
            if type(Techno_Matrix_Fish[i,j])==str: # Then we evaluate the string
            
                Techno_fish_num[i,j] = eval(Techno_Matrix_Fish[i,j])
                
            else:        # Then we just take the float
                Techno_fish_num[i,j] = Techno_Matrix_Fish[i,j]
    

    [Tech_Matrix,
     eco_FCR_0, 
     bio_FCR_0,
     Loss_ratio_INC_0, 
     ratio_loss_biological_INC_0,_] = mainfunc.calculate_economic_indicators(Techno_fish_num,
                            index_300g,
                            index_growing_DK,
                            index_roe,
                            index_FU,
                            index_feed,
                            index_dead_fish,
                            index_micro_compound) 
    




    return Tech_Matrix,eco_FCR_0,bio_FCR_0, Loss_ratio_INC_0, ratio_loss_biological_INC_0





[Tech_Matrix,ECO_FCR_0,bio_FCR_0, Loss_ratio_INC_0, ratio_loss_biological_INC_0]=initial(Techno_Matrix_Fish,
            demand_vector,
            index_300g,
            index_growing_DK,
            index_roe,
            index_FU,
            index_feed,
            index_dead_fish)

#Incumbent dose for drugs


supply_ini=linalg.solve(Tech_Matrix, demand_vector)

ini_dose_drug1 = supply_ini[index_drug1]
ini_dose_drug2 = supply_ini[index_drug2]
ini_dose_drug3 = supply_ini[index_drug3]
ini_dose_drug4 = supply_ini[index_drug4]



""" 2 )Initialize brightway and calculate LCIA for 1 unit of each technosphere input """


bw.projects.set_current('Fish_and_micro_2103') 



# Loading Ecoinvent
Ecoinvent = bw.Database('ecoinvent 3.8 conseq')


# Loading foreground database

FISHMIC = bw.Database('AH_combi_1')



biosphere=bw.Database('biosphere3')

# Loading foreground database


MICAH = bw.Database('Micro_for_combi')

# Collecting the foreground activities which are only connected to a background activity for the microalgal compound

list_micro_algae_FU = []

list_micro_algae_names = []

        


for act in MICAH:
    # The original Molecule production activity
    if act['name'] == 'Molecule production PBR':
        # print('oik')
        Molprod = MICAH.get(act['code'])   


for act in MICAH:
    # The original wastewater treatment activity
    if act['name'] == 'Wastewater treatment (without electricity) PBR':
        wastewater = MICAH.get(act['code'])


LCIdict_micro={}

for exc in list(Molprod.exchanges()):

        if exc['type']!='production':
            
            exchange1 = MICAH.get(exc['input'][1])  
            
            name_exchange = exchange1['name']  
            
            list_micro_algae_names.append(name_exchange)
        
            list_micro_algae_FU.append({exchange1 : 1})
            LCIdict_micro[name_exchange]=0
    

# Electricity for europe


list_electricity_CNTR_FU = []

list_electricity_CNTR_names = []



for act in Ecoinvent:
            if 'market for electricity' in act['name']:
                if 'medium voltage' in act['name']:
                    if 'municipal' not in act['name']:
                        if 'label-certified' not in act['name']:
                            if act['location'] in ['ES','SE','IT','PT','GB','FR','GR','CY','IE','DE']:  
                                print(act)  
                                list_electricity_CNTR_FU.append({act:1})
                                list_electricity_CNTR_names.append(act["location"])




### Additional variables for the regionalization of the wastewater treatment of microalgae

electricity_low_voltage_input_per_m3_wastewater = 0.20571
electricity_high_voltage_input_per_m3_wastewater = 0.0086946





list_meth =[('ReCiPe Midpoint (H)', 'climate change', 'GWP100'), 
            ('ReCiPe Midpoint (H)', 'human toxicity', 'HTPinf'),
            ('ReCiPe Midpoint (H)', 'freshwater ecotoxicity', 'FETPinf'),
             ('TRACI', 'environmental impact', 'eutrophication'),
             ('ReCiPe Midpoint (H)', 'terrestrial ecotoxicity', 'TETPinf'),
             ('ReCiPe Midpoint (H)', 'fossil depletion', 'FDP'),
             ('ReCiPe Midpoint (H)', 'terrestrial acidification', 'TAP100'),
             ('ReCiPe Midpoint (H)', 'freshwater eutrophication', 'FEP'),
             ('ReCiPe Midpoint (H)', 'ozone depletion', 'ODPinf')]


list_cfs= [bw.Method((meth)).load() for meth in list_meth]


# Add uncertain parameters corresponding to the impacts per kg of drugs
for meth in list_meth:
    
    name_param = "impact_Drug_prod" + meth[-1]

    
    Fish_farm_and_compound_effect_dict[name_param]=0

    
    Fish_farm_and_compound_effect_dict_distributions[name_param]=["unif",[0,10,1000,0,0]]  

    
    


                      
#activities_fish_background

list_fish_FU = []
for name_act in activities_fish_background:
    
    for act in FISHMIC:
        
        if act["name"] == name_act:
            
            list_fish_FU.append({act:1})


"""Calculating impacts for necessary background activites"""



# Collect Characterization matrices in a dictionnary



C_matrixes ={}
# Initialize a LCA object, whatever object.

Lca=bw.LCA(list_fish_FU[2],('ReCiPe Midpoint (H)', 'climate change', 'GWP100'))
Lca.lci()
Lca.lcia()


for meth in list_meth:
    Lca.switch_method(meth)
    C_matrixes[meth]=Lca.characterization_matrix

    



# Function to modify the electricty impacts according to the national mix

# Necessary
list_processes_electricity_micro = []
for exc in list(Molprod.exchanges()):
        if exc['type']!='production':
            
            exchange1 = MICAH.get(exc['input'][1])  
            
            exchange1_name=exchange1['name']
            print(exchange1)
            for exc in list(exchange1.exchanges()):
                
               # "MJ heat biogas AD PBR" has an input of the modified anaerobic digestion actvitiy 
               #from the foreground database, it is not electricity and it can't be found with Ecoinvent.get() 
                
                if exc['type']=='technosphere' and exchange1_name!="MJ heat biogas" and exchange1_name!='anaerobic digestion of manure Europe':
                    
                    act_background = Ecoinvent.get(exc['input'][1])  
                    
                    name_background = act_background['name']
                    
                    #print(name_background)
                    if 'electricity' in name_background:
                        print('ok')
                        list_processes_electricity_micro.append(exchange1_name)
          
 

# For natural_gas


   
Dict_FU_NG = {}

list_skarka=['ES','SE','IT','PT','GB','FR','GR','CY','IE','DE'] 


for loc in list_skarka:
    
    if loc != "CY" and loc!="PT":
        for act in Ecoinvent:
            if "market for natural gas, high pressure" in act["name"] and act["location"] == loc:
                #print(act)
                Dict_FU_NG[loc]={Ecoinvent.get(act["code"]):1}
                
    else:  
        for act in Ecoinvent:
            
            if "market group for natural gas, high pressure" in act["name"]:
                #print(act["location"])
            
                if act["location"] == "Europe without Switzerland":
                   # print(act)
                    Dict_FU_NG[loc]={Ecoinvent.get(act["code"]):1}        
                    
            
            
list_natural_gas_CNTR_FU = [Dict_FU_NG[a] for a in Dict_FU_NG]
        



list_natural_gas_CNTR_names=[a+"gas" for a in list_skarka]





# Create dictionnaries with incumbent production outputs for the fish farm growth stages
# Needed to ease the calculations later


# the key is the name of the associated loss level parameter
Dict_incumbent_outputs_growth_stages_loss_level = {'HAIT_loss_lev': 2400,   
                                                   'FRIT_loss_lev': 12000,
                                                   'GOIT1_loss_lev': 64000, 
                                                   'GOIT2_loss_lev': 96000,
                                                   'GOIT1bis_loss_lev': 32000, 
                                                   'GODK_loss_lev': 319822+750100+194750, 
                                                   'SFDK1_loss_lev': 3512877, 
                                                   'SFDK2_loss_lev': 4934867}
                                 
                               
                                
                                 



# the key is the name of the associated loss reduction parameter

Dict_incumbent_outputs_growth_stages_loss_red = {'HAIT_loss_red': 2400,   
                                                   'FRIT_loss_red': 12000,
                                                   'GOIT1_loss_red': 64000, 
                                                   'GOIT2_loss_red': 96000,
                                                   'GOIT1bis_loss_red': 32000, 
                                                   'GODK_loss_red': 319822+750100+194750, 
                                                   'SFDK1_loss_red': 3512877, 
                                                   'SFDK2_loss_red': 4934867}




# Create dictionnaries with incumbent losses outputs for the fish farm growth stages


# the key is the name of the associated loss reduction parameter

Dict_incumbent_losses_growth_stages_loss_red = {'HAIT_loss_red': 400,   
                                                   'FRIT_loss_red': 0,
                                                   'GOIT1_loss_red': 16000, 
                                                   'GOIT2_loss_red': 0,
                                                   'GOIT1bis_loss_red': 0, 
                                                   'GODK_loss_red': 31053, 
                                                   'SFDK1_loss_red': 266783, 
                                                   'SFDK2_loss_red': 0}




# the key is the name of the associated loss level parameter
Dict_incumbent_losses_growth_stages_loss_level = {'HAIT_loss_lev': 400,   
                                                   'FRIT_loss_lev': 0,
                                                   'GOIT1_loss_lev': 16000, 
                                                   'GOIT2_loss_lev': 0,
                                                   'GOIT1bis_loss_lev': 0, 
                                                   'GODK_loss_lev': 31053, 
                                                   'SFDK1_loss_lev': 266783, 
                                                   'SFDK2_loss_lev': 0}
                                 

Dict_FCR_bio = {'HAIT_FCR': 2700/(2400+400),   
            'FRIT_FCR': 12400/(12000-2400),
            'GOIT1bis_FCR': 33333/(32000-12000), 
            'GODK_FCR': 1135494/(319822+750100+194750+31053-162266), 
            'SFDK1_FCR': 2345650/(3512877+266783-1926316), 
            'SFDK2_FCR': 3152350/(4934867+710210-3512877)}
                        


       

""" Generate the grid of locations for the microalgal production"""


def pre_download_climatic_data(latsec,methods_selected,biosph,MICAH,Ecoinvent,cultivationperiod ):
    
    """Functio nwhich creates a grid of random location over the european countries 
    and download the climatic data associated to each location   """
    
    
    geodataframe_metrop, gdf_points, list_points_grid = geo.geodataframe_initialize_2(latsec,methods_selected,biosph,MICAH,Ecoinvent )

    
    for country_index in range(geodataframe_metrop.shape[0]):
        
        list_points = geodataframe_metrop['random_points'].iloc[country_index]
    
        for point_index in range(len(list_points)):
        
        # corresponding rank (row) of this point in the geodataframe containing points
        
            
        
            point = list_points[point_index]
            print("point",point)
            

            for month in cultivationperiod:
        
                solardata.Qreceived_bym2PBR_month(round(point[1],3),     # Longitude (x) is given before latitude (y) in geodataframe)
                                                        round(point[0],3),
                                                        month,
                                                        90,
                                                        1.5,
                                                        0.3,
                                                        0.05,
                                                        0.3,
                                                        30,
                                                        30)
    
    x = datetime.datetime.now()
    
    month=str(x.month)
    day=str(x.day)
    microsec=str(x.strftime("%f"))
                 
    name_geo_countries ='geodataframe_input'+"_"+month+"_"+day+"_"+microsec
    
    name_geo_points ='gdfpoints_input'+"_"+month+"_"+day+"_"+microsec

    name_file_list_points_grid = 'list_points_grid_input'+"_"+month+"_"+day+"_"+microsec

    export_pickle_2(geodataframe_metrop, name_geo_countries, "geo objects")
    
    export_pickle_2(gdf_points, name_geo_points, "geo objects")

    export_pickle_2(list_points_grid, name_file_list_points_grid, "geo objects")

    
    return geodataframe_metrop, gdf_points  , list_points_grid             
    



# Latitude distance between two random locations on the grid
latsec= 1

# From April to September
cultivationperiod=[4,5,6,7,8,9] 


# Generate grid of random locations and downolad climatic data
geodataframe_metrop,gdf_points, list_points_grid = pre_download_climatic_data(latsec,list_meth,biosphere,MICAH,Ecoinvent,cultivationperiod )



min_lat= min([coord[0][1] for coord in list_points_grid])  # get minimum latitude
max_lat= max([coord[0][1] for coord in list_points_grid])  # get maximum latitude


# The mode latitude will be the mode of a triangular distribution over Europe which detemrines a production mix
# A produciton mix with a certain mode is therefore more concentrated around this latitude.
Production_mix_dict_distributions["mode_lat"]= ["pert_nested",[0, min_lat+0.5, max_lat-0.5, 0, 0]]


# 35 Techno-operational setup per location and strain in a mix
size_nested_sample = 35


      

# Import the background mc file. Replace with the desired file
list_array_total_mc_sorted = importpickle("../Background_mc/chunk_Background_mc_12_15_919097_size=104832.pkl")

# We will simulate as many iterations as there is in the background
size=len(list_array_total_mc_sorted)





list_FU_combined_mc = list_micro_algae_FU + list_fish_FU + list_electricity_CNTR_FU + list_natural_gas_CNTR_FU
list_FU_combined_names_mc = list_micro_algae_names + activities_fish_background + list_electricity_CNTR_names + list_natural_gas_CNTR_names


drug_inputs_names=['Drug_1_prod','Drug_2_prod','Drug_3_prod','Drug_4_prod','Drug_5_prod']






"""3) RUN"""

random.seed()
np.random.seed()


dict_correspondance_techno_growth_stagenames={0:"HAIT",
                                              1:"FRIT",
                                              5:"GOIT2_para",
                                              6:"GOIT1bis",
                                              7:"GODK",
                                              8:"SFDK1",
                                              9:"SFDK2"}

index_growth_stages_to_modif = [0, 1, 6, 7, 8, 9]



time3=time()
res,r= mainfunc.simulations_fish_micro(Tech_opdict_distributions,
                           Biodict_distributions,
                           Locationdict_distributions,
                           Physicdict_distributions,
                           Fish_farm_and_compound_effect_dict_distributions,
                           Production_mix_dict_distributions,
                           size, 
                           Tech_opdict,
                           Biodict,
                           Locationdict,
                           Physicdict,
                           Fish_farm_and_compound_effect_dict,
                           Production_mix_dict,
                           LCIdict_micro,
                           cultivationperiod,
                           fishfeed_table_withNP,
                           elemental_contents,
                           size_nested_sample,
                           list_points_grid,
                           demand_vector,
                           Techno_Matrix_Fish,
                           list_FU_combined_names_mc,
                           list_array_total_mc_sorted,
                           activities_fish_background,
                           filter_names_ok_activities_fish,
                           list_processes_electricity_micro,
                           list_meth,
                           Dict_incumbent_losses_growth_stages_loss_level,
                           Dict_incumbent_losses_growth_stages_loss_red,
                           Dict_incumbent_outputs_growth_stages_loss_red,
                           Dict_incumbent_outputs_growth_stages_loss_level,
                           biochem_profile_feed_incumbent,
                           N_P_profile_feed_incumbent,
                           ingredient_profile_incumbent,
                           index_growth_stages,
                           index_feed,
                           index_dead_fish,
                           digestibility_list,
                           index_Nemissions,
                           index_Pemissions,
                           index_growth_stages_no_filter_laguna,
                           index_sludge,
                           electricity_low_voltage_input_per_m3_wastewater,
                           electricity_high_voltage_input_per_m3_wastewater,
                           index_roe,
                           Dict_FCR_bio,
                           list_cfs,
                           ECO_FCR_0,
                           ratio_loss_biological_INC_0,
                           index_micro_compound,
                           index_growing_DK,
                           index_300g,
                           drug_inputs_names,
                           index_FU,
                           dict_correspondance_techno_growth_stagenames,
                           index_growth_stages_to_modif,
                           bio_FCR_0,
                           index_biogas_updgrade,
                           index_N_substitution,
                           index_P_substitution,
                           index_heat_substitution)




timetot=time()-time3

print("tiiiime",timetot)




x = datetime.datetime.now()

month=str(x.month)
day=str(x.day)
microsec=str(x.strftime("%f"))
             


name_file='result_AHmicro'+"_"+month+"_"+day+"_"+microsec
name_file_contribution='result_contribution_AHmicro'+"_"+month+"_"+day+"_"+microsec


res.to_csv("../Outputs/"+str(name_file)+'.csv', sep=';', encoding='utf-8')

export_pickle_2(r, name_file_contribution, "Outputs")




