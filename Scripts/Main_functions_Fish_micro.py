# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 09:33:34 2022

@author: Pierre Jouannais, Department of Planning, DCEA, Aalborg University
pijo@plan.aau.dk

"""
'''
Script which produces the raw output with 1 data point being 1 LCA for Trout with and without compound.
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
# Set working directory to file location 
# (works only when executing the whole file and not only sections (Run Current cell))

# currentfolder=os.path.dirname(os.path.realpath(__file__))
# os.chdir(currentfolder)

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
from SALib.analyze import delta

import SALib

import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits import mplot3d

import shapefile as shp
import geopandas as gpd
import pysal as ps
from shapely.geometry import Polygon, mapping, Point


import Cultivation_simul_2nd as cultsimul
import Functions_for_physical_and_biological_calculations_3rd as functions
import Map_3rd as geo
import Retrieving_solar_and_climatic_data_2nd as solardata


import Technosphere_matrix_modifications_3rd as techno_modif

import ray

import pert


import scipy
import numpy as np
from scipy import linalg

from ast import literal_eval











"""Sampling Functions """






def sampling_func_nested(Tech_opdict,
                         nested_list_techno_op,
                         size_nested_sample):

    """Function which creates a Monte Carlo Sample of the "nested" uncertain parameters.
      In the current model, It is used to generate samples of uncertain Techno-operational parameters"""   
    
    
    sample_array_nested = np.array([[0]*len(nested_list_techno_op)]*size_nested_sample,dtype="float32")
    list_bounds_nested = [a[1][1:3]+[Tech_opdict[a[0]]] for a in nested_list_techno_op ]  
                       
    
    for row_index in range(sample_array_nested.shape[0]):

        sample_row =[pert(bounds_param[0],bounds_param[2],bounds_param[2],size=1)[0] for bounds_param in list_bounds_nested]
        
        sample_array_nested[row_index,:]=sample_row
   
    
   
    names_param_nested_techno_op = [a[0] for a in nested_list_techno_op]
   
    return sample_array_nested, names_param_nested_techno_op
   




def pert(a, b, c, *, size=1, lamb=4):
    """Random draw in a pert distribution """
    
    r = c - a
    alpha = 1 + lamb * (b - a) / r
    beta = 1 + lamb * (c - b) / r
    
    return a + np.random.beta(alpha, beta, size=size) * r
    









    
def sampling_func_total_montecarlo(Tech_opdict_distributions,
                  Biodict_distributions,
                  Locationdict_distributions, 
                  Physicdict_distributions,
                  Fish_farm_and_compound_effect_dict_distributions,
                  Production_mix_dict_distributions,
                  size):
    '''Function which returns a random sample for the input space of 
    non-constant parameters via Montecarlo sampling
    Inputs:
        
        # distributions dictionnaries : 
                  Biodict_distributions,
                  Locationdict_distributions, 
                  Physicdict_distributions
                  
        #size : Size of the sample 

    Outputs :
        
        -sample :  Generated sample. Array with 1 row = 1 combination of uncertain parameters
        -names_param :  List of names of the uncertain parameters
        -names_param_op :  List of names of the uncertain parameters from Tech_opdict_distributions
        -names_param_bio :  List of names of the uncertain parameters from Biodict_distributions
        -names_param_geo :  List of names of the uncertain parameters from Locationdict_distributions
        -names_param_phy :  List of names of the uncertain parameters from Physicdict_distributions
        - names_param_fish_farm_compound_effect :  List of names of the uncertain parameters from fish_farm_compound_effect_distributions
        -names_param_prod_mix : List of names of the uncertain parameters from prod_mix_distributions
        -nested_list_techno_op :  List of names of the uncertain parameters from Techno_op_distributions
         
        
        
        
    '''


    # Creating a copy of the distribtutions dictionnaires containing only the
    # parameters with distributions(not the constant ones)

    # Tech_opdict
    Tech_opdict_distributions_input = Tech_opdict_distributions.copy()
    
    to_be_deleted_op = []
    
    nested_list_techno_op = []
    nested_list_prod_mix = []
    
    # Collecting the parameters that are constant
    
    for param in Tech_opdict_distributions_input:
        
                
        if Tech_opdict_distributions_input[param][0] == 'unique' or Tech_opdict_distributions_input[param][0] == 'binary':
            to_be_deleted_op.append(param)
    
        elif Tech_opdict_distributions_input[param][0] == 'pert_nested':
                nested_list_techno_op.append([param,Tech_opdict_distributions_input[param][1]])
                
                
    # Deleting the parameters that are constant
    
    for a in to_be_deleted_op:
        Tech_opdict_distributions_input.pop(a)

    # Biodict

    Biodict_distributions_input = Biodict_distributions.copy()

    to_be_deleted_bio = []

    for param in Biodict_distributions_input:
        if Biodict_distributions_input[param][0] == 'unique' or Biodict_distributions_input[param][0] == 'binary':
            to_be_deleted_bio.append(param)

    for a in to_be_deleted_bio:
        Biodict_distributions_input.pop(a)

    # Geography

    Locationdict_distributions_input = Locationdict_distributions.copy()

    to_be_deleted_geo = []

    for param in Locationdict_distributions_input:
        if Locationdict_distributions_input[param][0] == 'unique' or Locationdict_distributions_input[param][0] == 'binary':
            to_be_deleted_geo.append(param)

    for a in to_be_deleted_geo:
        Locationdict_distributions_input.pop(a)

    # Physics

    Physicdict_distributions_input = Physicdict_distributions.copy()

    to_be_deleted_phy = []

    for param in Physicdict_distributions_input:
        if Physicdict_distributions_input[param][0] == 'unique' or Physicdict_distributions_input[param][0] == 'binary':
            to_be_deleted_phy.append(param)

    for a in to_be_deleted_phy:
        Physicdict_distributions_input.pop(a)

    # Fish farm and compound effect

    Fish_farm_and_compound_effect_dict_distributions_input = Fish_farm_and_compound_effect_dict_distributions.copy()

    to_be_deleted_fish = []

    for param in Fish_farm_and_compound_effect_dict_distributions_input:
        
        if Fish_farm_and_compound_effect_dict_distributions_input[param][0] == 'unique' or Fish_farm_and_compound_effect_dict_distributions_input[param][0] == 'binary':
            
            to_be_deleted_fish.append(param)

    for a in to_be_deleted_fish:
        Fish_farm_and_compound_effect_dict_distributions_input.pop(a)
        
        
    # Production mix

    Production_mix_dict_distributions_input = Production_mix_dict_distributions.copy()

    to_be_deleted_prod_mix = []

    for param in Production_mix_dict_distributions_input:
        
        if Production_mix_dict_distributions_input[param][0] == 'unique' or Production_mix_dict_distributions_input[param][0] == 'binary':
            
            to_be_deleted_prod_mix.append(param)
            
            
        elif Production_mix_dict_distributions_input[param][0] == 'pert_nested':
            
            nested_list_prod_mix.append([param,Production_mix_dict_distributions_input[param][1]])
                
              

    for a in to_be_deleted_prod_mix:
        Production_mix_dict_distributions_input.pop(a)        




    # Collecting names, bounds , dists 
    names_param = []
    bounds = []
    dists = []
    

    # 1) Operational
    
    names_param_op = []

    count_col = -1
    col_with_triang = [] # Triangular distributions will need post processing
    actual_lower_bounds = []

    for param in Tech_opdict_distributions_input:
        
        count_col += 1
        
        names_param_op.append(param)
        
        distrib = Tech_opdict_distributions_input[param][0]
        
        if distrib == 'unif' or distrib =='pert_nested':  # then bounds = upper bound, lower bound

            bounds.append([Tech_opdict_distributions_input[param]
                           [1][1], Tech_opdict_distributions_input[param][1][2]])
            
            dists.append('unif')


        elif distrib == 'norm':  # then bounds = mean, sd

            bounds.append([Tech_opdict_distributions_input[param]
                           [1][3], Tech_opdict_distributions_input[param][1][4]])
            
            dists.append(distrib)



        elif distrib == 'triang':  # then bounds = width, location of the peak in % of width, Assume lower bound = 0


            dists.append(distrib)


            
            bounds.append([Tech_opdict_distributions_input[param][1][1],
               Tech_opdict_distributions_input[param][1][3],
               Tech_opdict_distributions_input[param][1][2]])            

    # 2) Biological

    names_param_bio = []

    for param in Biodict_distributions_input:

        count_col += 1
        
        names_param_bio.append(param)
        
        distrib = Biodict_distributions_input[param][0]
        
        dists.append(distrib)
        
        if distrib == 'unif':  # then bounds = upper bound, lower bound

            bounds.append([Biodict_distributions_input[param][1]
                           [1], Biodict_distributions_input[param][1][2]])

        elif distrib == 'norm':  # then bounds = mean, sd

            bounds.append([Biodict_distributions_input[param][1]
                           [3], Biodict_distributions_input[param][1][4]])

        elif distrib == 'triang':  # then bounds = width, location of the peak in % of width, Assume lower bound = 0


            
            bounds.append([Biodict_distributions_input[param][1][1],
               Biodict_distributions_input[param][1][3],
               Biodict_distributions_input[param][1][2]])

    # 3) Geography
    
    names_param_geo = []

    for param in Locationdict_distributions_input:

        count_col += 1

        names_param_geo.append(param)

        distrib = Locationdict_distributions_input[param][0]

        dists.append(distrib)

        if distrib == 'unif':  # then bounds = upper bound, lower bound

            bounds.append([Locationdict_distributions_input[param]
                           [1][1], Locationdict_distributions_input[param][1][2]])

        elif distrib == 'norm':  # then bounds = mean, sd

            bounds.append([Locationdict_distributions_input[param]
                           [1][3], Locationdict_distributions_input[param][1][4]])

        elif distrib == 'triang':  # then bounds = width, location of the peak in % of width, Assume lower bound = 0


            bounds.append([Locationdict_distributions_input[param][1][1],
                           Locationdict_distributions_input[param][1][3],
                           Locationdict_distributions_input[param][1][2]])



    # 3) Physics
    
    names_param_phy = []

    for param in Physicdict_distributions_input:
        
        distrib = Physicdict_distributions_input[param][0]
        
        dists.append(distrib)

        count_col += 1
        
        names_param_phy.append(param)

        if distrib == 'unif':  # then bounds = upper bound, lower bound

            bounds.append([Physicdict_distributions_input[param]
                           [1][1], Physicdict_distributions_input[param][1][2]])

        elif distrib == 'norm':  # then bounds = mean, sd

            bounds.append([Physicdict_distributions_input[param]
                           [1][3], Physicdict_distributions_input[param][1][4]])

        elif distrib == 'triang':  # then bounds = width, location of the peak in % of width, Assume lower bound = 0


            
            bounds.append([Physicdict_distributions_input[param][1][1],
                       Physicdict_distributions_input[param][1][3],
                       Physicdict_distributions_input[param][1][2]])
            




    # 4) Fish farm and compound effect
    
    names_param_fish_farm_compound_effect = []

    for param in Fish_farm_and_compound_effect_dict_distributions_input:
        
        distrib = Fish_farm_and_compound_effect_dict_distributions_input[param][0]
        
        dists.append(distrib)

        count_col += 1
        
        names_param_fish_farm_compound_effect.append(param)

        if distrib == 'unif':  # then bounds = upper bound, lower bound

            bounds.append([Fish_farm_and_compound_effect_dict_distributions_input[param]
                           [1][1], Fish_farm_and_compound_effect_dict_distributions_input[param][1][2]])

        elif distrib == 'norm':  # then bounds = mean, sd

            bounds.append([Fish_farm_and_compound_effect_dict_distributions_input[param]
                           [1][3], Fish_farm_and_compound_effect_dict_distributions_input[param][1][4]])
            

        elif distrib == 'triang':  # then bounds = width, location of the peak in % of width, Assume lower bound = 0


            bounds.append([Fish_farm_and_compound_effect_dict_distributions_input[param][1][1],
                       Fish_farm_and_compound_effect_dict_distributions_input[param][1][3],
                       Fish_farm_and_compound_effect_dict_distributions_input[param][1][2]])            
                        
            



    # 5) Prod mix
    
    names_param_prod_mix = []

    for param in Production_mix_dict_distributions_input:
        
        distrib = Production_mix_dict_distributions_input[param][0]
        

        count_col += 1
        
        names_param_prod_mix.append(param)

        if distrib == 'unif' or distrib =='pert_nested':  # then bounds = upper bound, lower bound

            bounds.append([Production_mix_dict_distributions_input[param]
                           [1][1], Production_mix_dict_distributions_input[param][1][2]])
            
            dists.append("unif")


        elif distrib == 'norm':  # then bounds = mean, sd

            bounds.append([Production_mix_dict_distributions_input[param]
                           [1][3], Production_mix_dict_distributions_input[param][1][4]])
            
            dists.append(distrib)

        elif distrib == 'triang':  # then bounds = width, location of the peak in % of width, Assume lower bound = 0

            
            dists.append(distrib)

            
            bounds.append([Production_mix_dict_distributions_input[param][1][1],
            Production_mix_dict_distributions_input[param][1][3],
            Production_mix_dict_distributions_input[param][1][2]])
                        


            

    names_param = names_param_op + names_param_bio + names_param_geo + names_param_phy + names_param_fish_farm_compound_effect + names_param_prod_mix
    
    

    # Create Empty array for sample
    

    sample_array = np.array([[0]*len(names_param)]*size,dtype="float32")

    def random_draw(dist_param,bounds_param):
        
        if dist_param == "unif":
            value = np.random.uniform(bounds_param[0], bounds_param[1], size=None)
            ##print("value",value)
        if dist_param == "norm":
            value = np.random.normal(bounds_param[0], bounds_param[1], size=None)
            ##print("value",value)
            
        if dist_param == "triang":
            value = np.random.triangular(bounds_param[0], bounds_param[1], bounds_param[2])
            ##print("value",value)
        return(value)    
            


    # Random draw
    for row_index in range(sample_array.shape[0]):

        
        sample_row =[random_draw(dist_param,bounds_param) for dist_param,bounds_param in zip(dists,bounds)]
        sample_array[row_index,:]=sample_row






    return sample_array, names_param, names_param_op, names_param_bio, names_param_geo, names_param_phy, names_param_fish_farm_compound_effect, names_param_prod_mix,   nested_list_techno_op






""" Calculating functions"""





'''Main calculating functions '''

def LCI_one_strain_uniquevalues(Biodict,
                                Physicdict,
                                Tech_opdict,
                                Locationdict,
                                LCIdict_micro,
                                months_suitable_for_cultivation,
                                fishfeed_table_withNP,
                                elemental_contents,
                                biochem_profile_feed_incumbent,
                                electricity_low_voltage_input_per_m3_wastewater,
                                electricity_high_voltage_input_per_m3_wastewater,
                                ingredient_profile_incumbent,
                                N_P_profile_feed_incumbent):

    '''Calculate the LCI for one set of parameters given in input by simulating 
    the cultivation and scaling the values to the FU.
    
    Inputs:
        #Biodict : Dictionnary with biological parameters
        #Physicdict : Dictionnary with physical parameters
        #Tech_opdict : Dictionnary with techno-operational parameters
        #LCIdict_micro : Initialized LCI dictionnary
        #months_suitable_for_cultivation : Months for cultivation ; 
        list of month numbers : [a,b,c]
        
        #fraction_maxyield : Fraction of the maximum yield achieved ; .
        #fishfeed_table : DataFrame with fish feed composition
        #elemental_contents : Table with elemental compositons of macronutrients
        

    Outputs:
        # LCIdict_micro_updated : Dictionnary containing the calculated LCI
        
        Other variables resulting from the simulation  for further investigation 
        and review of the code (Non necessary):
            
        # surfaceyield : Areal yield ; kg.m-2 
        # volumetricyield : Volmetric yield ; kg.m-3 
        # optimization_performance : Fish feed Substitution alrogithm performance
        # needed_dbio_check : Non necessary
        # substitution_check : Non necessary
        # total_production_kg_dw : Total production ; kg dw
        # total_production_harvested_kg_dw : Actual harvested production harvested ; kg dw
        # total_production_loss_kg_dw : Production not harvested  ; kg dw
        # conc_waste_water_nutrient_N, : N Concentration in wastewater ; kg.m-3
        # Same for all elements
        #conc_waste_water_biomass : Biomass Concentration  in wastewater ; kg.m-3
        #bioact_molec_dbio : Molecule Concentration  in the dried biomass ; kg. kg dbio -1
        # min_centrifugation_rate_m3_h : Obsolete
        # max_centrifugation_rate_m3_h : Obsolete
        # totalwatercentrifuged : Total volume of water centrifuged ; m3
        # tubelength : Total tube length over 1m2 ; m
        # facilityvolume : Cultivation volume over 1m2 ; m3
        # exchange_area : Exchange area tube/air over 1m2 ; m
        # totalcooling_thermal : Thermal Cooling needed per m2
        (not via Heat Exchanger) ; kWh
        
    '''

    
        
    
    LCIdict_micro_updated = LCIdict_micro.copy()
    
    # Collecting all values from the dictionnaries and creating local variables

    # Tech_opdict 
    
    height = Tech_opdict['height']
    tubediameter = Tech_opdict['tubediameter']
    gapbetweentubes = Tech_opdict['gapbetweentubes']
    horizontaldistance = Tech_opdict['horizontaldistance']
    length_of_PBRunit = Tech_opdict['length_of_PBRunit']
    width_of_PBR_unit = Tech_opdict['width_of_PBR_unit']
    biomassconcentration = Tech_opdict['biomassconcentration']
    flowrate = Tech_opdict['flowrate']
    centrifugation_efficiency = Tech_opdict['centrifugation_efficiency']
    pumpefficiency = Tech_opdict['pumpefficiency']
    slurry_concentration = Tech_opdict['slurry_concentration']
    water_after_drying = Tech_opdict['water_after_drying']
    recyclingrateaftercentrifuge = Tech_opdict['recyclingrateaftercentrifuge']
    roughness = Tech_opdict['roughness']
    rhosuspension = Tech_opdict['rhosuspension']
    cleaningvolumeVSfacilityvolume = Tech_opdict['cleaningvolumeVSfacilityvolume']
    concentration_hypo = Tech_opdict['concentration_hypo']
    concentration_hydro = Tech_opdict['concentration_hydro']
    boilerefficiency = Tech_opdict['boilerefficiency']
    glass_life_expectancy = Tech_opdict['glass_life_expectancy']
    extraction = Tech_opdict['extraction']
    heat_pump = Tech_opdict['heat_pump']
    COP = Tech_opdict['COP']
    
    heat_input_per_kg_biomass_AD = Tech_opdict['heat_input_per_kg_biomass_AD']
    elec_input_per_kg_biomass_AD = Tech_opdict['elec_input_per_kg_biomass_AD']

    fraction_maxyield = Tech_opdict['fraction_maxyield']
    
    depth_well = Tech_opdict['depth_well']


    ##print("tubediameter",tubediameter)
    ##print("horizontaldistance",horizontaldistance)

    # Physicdict

    Cp = Physicdict['Cp']
    hconv = Physicdict['hconv']
    rhomedium = Physicdict['rhomedium']  
    rhowater = Physicdict['rhowater']  
    Cw = Physicdict['Cw']  
    CH4_LHV = Physicdict['CH4_LHV']
    
    # Biodict

    lipid_af_dw = Biodict['lipid_af_dw']
    rhoalgae = Biodict['rhoalgae']
    MJ_kglip = Biodict['MJ_kglip']
    MJ_kgcarb = Biodict['MJ_kgcarb']
    MJ_kgprot = Biodict['MJ_kgprot']
    PAR = Biodict['PAR']
    losspigmentantenna = Biodict['losspigmentantenna']
    quantumyield = Biodict['quantumyield']
    lossvoltagejump = Biodict['lossvoltagejump']
    losstoATPNADPH = Biodict['losstoATPNADPH']
    losstohexose = Biodict['losstohexose']
    lossrespiration = Biodict['lossrespiration']
    bioact_fraction_molec = Biodict['bioact_fraction_molec']
    random_no3 = Biodict['random_no3']
    Topt = Biodict['Topt']
    T_plateau = Biodict['T_plateau']
    
    ##print("lipid_af_dw",lipid_af_dw)
    #print("rhoalgae",rhoalgae)

    # Conversion to Tmax, Tmin for simpler calculation
    Tmax = Topt+T_plateau/2
    Tmin = Topt-T_plateau/2

    dcell = Biodict['dcell']
    incorporation_rate = Biodict['incorporation_rate']
    ash_dw = Biodict['ash_dw']
    nutrient_utilisation = Biodict['nutrient_utilisation']
    co2_utilisation = Biodict['co2_utilisation']
    phospholipid_fraction = Biodict['phospholipid_fraction']
    
    random_market_subst = Biodict['random_market_subst']

    random_night_monitoring = Biodict['random_night_monitoring']

    random_bioclass = Biodict['random_bioclass']
    
    
    fertiliser_substitutability_AD = Biodict['fertiliser_substitutability_AD']
    
    
    carbon_degradibilidy_AD = Biodict['carbon_degradibilidy_AD']
    
    
    
    

    
    # Locationdict

    lat = Locationdict['lat']
    long = Locationdict['long']
    azimuthfrontal = Locationdict['azimuthfrontal']
    
    #print("dcell",dcell)
    #print("depth_well",depth_well)

    # We assume that the plant will be purposely located next to a water source with Twell<TMin  :
    
    Twell=np.random.uniform(5, Tmin,size=None)   
    
    
        

    # Qualitative parameters are determined based on the probabilities
    # and a new entry is created in the LCI dict




    # Nsource

    if random_no3 < 0.5: # Then the source is Nitrate
        Nsource = 'no3'
    else:
        Nsource = 'nh3'

    LCIdict_micro_updated['Nsource'] = Nsource


    # Biochemical class
    
  
    if random_bioclass < (1/3):
        #print('okC1')
        biochemicalclass='lip'
        
    elif (1/3)<random_bioclass<(2/3)  :
        #print('okC2')
        biochemicalclass='prot'
        
    else:
        #print('okC2')
        biochemicalclass='carb'


    LCIdict_micro_updated['Bio_class'] = biochemicalclass


    # Thermoregulation at night
    if random_night_monitoring < (1/2) :
        
        night_monitoring = 'yes'

    else:
        night_monitoring = 'no'

    LCIdict_micro_updated['night_monitoring'] = night_monitoring


    # Market for substitution

    if random_market_subst < (1/3):
        
        market_for_substitution = 'animal feed'

    elif (1/3)<random_market_subst< (2/3):
        
        market_for_substitution = 'fish feed'  
        
    elif random_market_subst > (2/3):
        
        market_for_substitution = 'anaerobic digestion' 

    #print(market_for_substitution)

    LCIdict_micro_updated['market_for_substitution'] = market_for_substitution




    # Collecting PBR geometry

    geom = functions.PBR_geometry(height,
                                  tubediameter,
                                  gapbetweentubes,
                                  horizontaldistance,
                                  length_of_PBRunit,
                                  width_of_PBR_unit)
    tubelength = geom[1]
    facilityvolume = geom[0]
    exchange_area = geom[-1]

    # LCI values which do not depend on the cultivation simulation

    # Calculating biomass composition at different levels

    biomass_composition = functions.biomasscompo(
        lipid_af_dw,
        ash_dw,
        water_after_drying,
        phospholipid_fraction,
        elemental_contents)


    # ash-free biomass composition (af_dw)
    prot_af_dw = biomass_composition[0]
    carb_af_dw = biomass_composition[1]

    # Including ash (dw)
    lip_dw = biomass_composition[2]
    prot_dw = biomass_composition[3]
    carb_dw = biomass_composition[4]

    # After harvesting and drying  (dbio)
    lip_dbio = biomass_composition[5]
    prot_dbio = biomass_composition[6]
    carb_dbio = biomass_composition[7]
    ash_dbio = biomass_composition[8]

    # Elementary composition

    C_af_dw = biomass_composition[9]
    N_af_dw = biomass_composition[10]
    P_af_dw = biomass_composition[11]
    K_af_dw = biomass_composition[12]
    Mg_af_dw = biomass_composition[13]
    S_af_dw = biomass_composition[14]

    C_dw = biomass_composition[15]
    N_dw = biomass_composition[16]
    P_dw = biomass_composition[17]
    K_dw = biomass_composition[18]
    Mg_dw = biomass_composition[19]
    S_dw = biomass_composition[20]

    # Calculating the absolute bioactive molecule content in the dried biomass
    if biochemicalclass == 'lip':
        bioact_molec_dbio = bioact_fraction_molec * lip_dbio

    elif biochemicalclass == 'carb':
        bioact_molec_dbio = bioact_fraction_molec * carb_dbio

    elif biochemicalclass == 'prot':
        bioact_molec_dbio = bioact_fraction_molec * prot_dbio

    # Nutrients

    # Nitrogen consumption 

    # considering ash content

    N_demand = N_dw * ((1/bioact_molec_dbio) * (1/(1 - water_after_drying)))

    # Recycling part of the nutrients with supernatant
    N_input = (N_demand / nutrient_utilisation + N_demand *
               recyclingrateaftercentrifuge) / (1 + recyclingrateaftercentrifuge)

    N_waste = N_input - N_demand

    #Only the correct source of N is updated

    if Nsource == 'nh3':
        # Already as N in Ecoinvent
        LCIdict_micro_updated['nutrient supply from ammonium sulfate'] = N_input
        
        LCIdict_micro_updated['calcium nitrate production PBR'] = 0

    if Nsource == 'no3':
        LCIdict_micro_updated['nutrient supply from ammonium sulfate'] = 0
        
        # Conversion from N to Calcium nitrate 
        LCIdict_micro_updated['calcium nitrate production PBR'] = N_input/0.15


    # Phosphorus consumption 

    P_demand = P_dw * ((1/bioact_molec_dbio) * (1/(1 - water_after_drying)))

    P2O5_demand = P_demand/0.4366     # Conversion P to P2O5

    P2O5_input = (P2O5_demand / nutrient_utilisation + P2O5_demand *
                  recyclingrateaftercentrifuge) / (1 + recyclingrateaftercentrifuge)  # Recylcing

    P_waste = P2O5_input*0.4366 - P_demand

    LCIdict_micro_updated['P source production PBR'] = P2O5_input
    # C
    
    C_demand = C_dw * ((1 / bioact_molec_dbio) * (1 / (1 - water_after_drying)))
    
    CO2_demand = C_demand * (44/12)  # Conversion C to CO2

    CO2_input = CO2_demand / co2_utilisation
    
    CO2_direct_emission = CO2_input - CO2_demand

    LCIdict_micro_updated['Microalgae CO2 PBR'] = CO2_input
    LCIdict_micro_updated['CO2 direct emissions PBR'] = CO2_direct_emission

    # K

    K_demand = K_dw * ((1/bioact_molec_dbio) * (1/(1 - water_after_drying)))

    K2O5_demand = K_demand*1.2 # Conversion to K2O5

    K2O5_input = (K2O5_demand / nutrient_utilisation + K2O5_demand *
                  recyclingrateaftercentrifuge) / (1 + recyclingrateaftercentrifuge)  # Recycling

    K2O5_waste = K2O5_input - K2O5_demand
    K_waste = K2O5_input/1.2 - K_demand

    LCIdict_micro_updated['K source production PBR'] = K2O5_input

    # Mg
    
    Mg_demand = Mg_dw * ((1 / bioact_molec_dbio)*(1/(1 - water_after_drying)))
    MgSO4_demand = Mg_demand * (120.4/24.3)

    MgSO4_input = (MgSO4_demand / nutrient_utilisation + MgSO4_demand *
                   recyclingrateaftercentrifuge) / (1 + recyclingrateaftercentrifuge)  # Recycling

    Mg_input = MgSO4_input * (24.3/120.4)

    Mg_waste = Mg_input-Mg_demand

    LCIdict_micro_updated['Mg source production PBR'] = MgSO4_input

    # S
    
    S_demand = S_dw * ((1/bioact_molec_dbio) * (1/(1-water_after_drying)))

    S_input = MgSO4_input*(32/120.4)
    
    if Nsource == 'nh3':  # Then ammonium sulfate also brings sulfate
    
        # N input --> (NH4)2SO4 input --> S input
        
        S_input += (N_input/0.21) * 0.24
        
    S_waste = S_input-S_demand




    # Cultivation simulation

    # Intializing variables
    totalcooling_thermal = 0
    totalcooling = 0
    totalheating = 0
    totalproduction = 0
    totalproduction_loss = 0
    totalproduction_harvested = 0
    totalwaterpumpedfromthefacility = 0
    totalwaterpumpedfromthewell = 0
    total_elec_centrifuge = 0
    total_elec_mixing = 0
    totalwatercentrifuged = 0

    meantemp_at_harvest_time_cultivationperiod = 0
    min_centrifugation_rate_list = []
    max_centrifugation_rate_list = []

    # Simulating an average day for each month of the cultivation period
    for month in months_suitable_for_cultivation:

        # Calling the cultivation simulation function
        simulation_averageday = cultsimul.cultivation_simulation_timestep10(hconv,
                                                                            Twell,
                                                                            depth_well,
                                                                            lat,
                                                                            long,
                                                                            azimuthfrontal,  
                                                                            month,  
                                                                            Cp,  
                                                                            height,
                                                                            tubediameter,
                                                                            gapbetweentubes,
                                                                            horizontaldistance,
                                                                            length_of_PBRunit,
                                                                            width_of_PBR_unit,  
                                                                            rhoalgae, 
                                                                            rhomedium,
                                                                            rhosuspension, 
                                                                            dcell, 
                                                                            Tmax,
                                                                            Tmin,
                                                                            Biodict, 
                                                                            ash_dw,
                                                                            Nsource,  
                                                                            fraction_maxyield,  
                                                                            biomassconcentration,
                                                                            flowrate,  
                                                                            centrifugation_efficiency,
                                                                            pumpefficiency,  
                                                                            slurry_concentration,
                                                                            water_after_drying,
                                                                            recyclingrateaftercentrifuge,
                                                                            night_monitoring,
                                                                            elemental_contents)

        # Collecting results and multiplying by 
        # average number of days in a month : 30.4

        monthly_heating_energy = simulation_averageday[1]*30.4 #kWh
      
        monthly_waterpumped_from_the_facility = simulation_averageday[2]*30.4  # L

        monthly_waterpumped_from_the_well = simulation_averageday[3]*30.4  # L

        monthly_production = simulation_averageday[4]*30.4  # g dw

        monthly_production_harvested = simulation_averageday[5]*30.4  # g dw

        monthly_production_loss = simulation_averageday[6]*30.4  # g dw

        monthly_volumetric_yield = simulation_averageday[7]*30.4  # g dw.m-3

        monthly_energy_tocentrifuge = simulation_averageday[8]*30.4  # kWh
        
        collectedtemperatureevolution = simulation_averageday[9]  # °C

        monthly_cooling_energy_thermal = simulation_averageday[17]*30.4  # kWh

        monthly_cooling_energy = simulation_averageday[17]*30.4  # kWh

        # Collection min and max centrifugation rate (Obsolete)

        # list centrifugation rate L.s-1
        list_centrifugation_rate_wholeunit = simulation_averageday[20]

        # list centrifugation rate L.s-1
        water_centrifuged = simulation_averageday[21]*30.4

        
        list_centrifugation_rate_wholeunit_not_0 = [
            i for i in list_centrifugation_rate_wholeunit if i != 0]
        
        min_centrifugation_rate_list.append(
            min(list_centrifugation_rate_wholeunit_not_0))

        max_centrifugation_rate_list.append(
            max(list_centrifugation_rate_wholeunit_not_0))
        
        # Temperature : Some processes have temperature as an input for their 
        # electricity consumption 
        #


        # For Mixing : Mixing is needed at any time of the day and night.
        meantemp_daytotal = (sum(collectedtemperatureevolution) / 
            len(collectedtemperatureevolution))

        monthly_Electricity_mixing_day = functions.mixing_perday(
            rhosuspension,
            tubediameter,
            pumpefficiency,
            flowrate,
            roughness,
            meantemp_daytotal,
            biomassconcentration,
            tubelength)[0] * 30.4  # MJ.m-2.month
        
        # For Drying : Drying requires to heat the slurry to 100 C and the
        # energy will depend on the initial temperature : temperature of the
        # culture at harvest time (9PM)
        
        temp_at_harvest_time = collectedtemperatureevolution[3780]



        # Summing over months in the loop
        totalproduction += monthly_production  # g dw
        
        totalproduction_harvested += monthly_production_harvested  # g dw
        
        totalproduction_loss += monthly_production_loss  # g dw

        totalcooling_thermal += monthly_cooling_energy_thermal # kWh

        totalcooling += monthly_cooling_energy  # kWh
        
        totalheating += monthly_heating_energy  # kWh
        
        total_elec_centrifuge += monthly_energy_tocentrifuge  # kWh
        
        total_elec_mixing += monthly_Electricity_mixing_day/3.6  # conversion to kWh

        totalwaterpumpedfromthefacility += monthly_waterpumped_from_the_facility  # L
        
        totalwaterpumpedfromthewell += monthly_waterpumped_from_the_well  # L

        totalwatercentrifuged += water_centrifuged  # L

        # Collecting the mean temperature over the cultivation period

        # For drying
        meantemp_at_harvest_time_cultivationperiod += temp_at_harvest_time/len(months_suitable_for_cultivation) # °C

    # End of the loop
    # Collecting min and max centrifugation rate during cultivation period (Obsolete)

    min_centrifugation_rate_m3_h = min(
        min_centrifugation_rate_list)*3.6  # m3.h-1
    
    max_centrifugation_rate_m3_h = max(
        max_centrifugation_rate_list)*3.6  # m3.h-1

    # Total production conversion to kg

    total_production_kg_dw = totalproduction/1000  # kg dw

    total_production_harvested_kg_dw = totalproduction_harvested/1000

    total_production_loss_kg_dw = totalproduction_loss/1000

    # Adding the energy for the initial heating of the well water 

    # Water of the well is heaten to Tmin
    if Twell < Tmin:
        initalheating = facilityvolume*Cp*(Tmin-Twell)/3.6  # kWh

    else: # If the water is already warm enough
        initalheating = 0  # KwH

    #Updating LCI with calculated values
    
    # Scaling down to 1 kg of  molecule in dried biomass
    
    
    if heat_pump=="yes": # use COP
    
        LCIdict_micro_updated['Heating kWh PBR'] = ((totalheating + initalheating)/total_production_harvested_kg_dw)*(
            1/bioact_molec_dbio)*(1/(1-water_after_drying))/COP 
        
        LCIdict_micro_updated['Cooling kWh PBR'] = (
            totalcooling_thermal/total_production_harvested_kg_dw)*(1/bioact_molec_dbio)*(1/(1-water_after_drying))/(COP-1)
    
    else:  #electric heater and cooling heat exchanger
        LCIdict_micro_updated['Heating kWh PBR'] = ((totalheating + initalheating)/total_production_harvested_kg_dw)*(
            1/bioact_molec_dbio)*(1/(1-water_after_drying))
        
        LCIdict_micro_updated['Cooling kWh PBR'] = (
            totalcooling/total_production_harvested_kg_dw)*(1/bioact_molec_dbio)*(1/(1-water_after_drying))
        
        
        
        
    LCIdict_micro_updated['Electricity centrifuge kWh PBR'] = (
        total_elec_centrifuge/total_production_harvested_kg_dw) * (1/bioact_molec_dbio)*(1/(1-water_after_drying))

    LCIdict_micro_updated['Electricity mixing kWh PBR'] = (
        total_elec_mixing/total_production_harvested_kg_dw) * (1/bioact_molec_dbio)*(1/(1-water_after_drying))




    # Pumping water from the well and facility
    
    #Calling the function  for depth = well depth
    energy_perm3_fromthewell = functions.pumping_per_m3(
        rhowater, depth_well, pumpefficiency)
    
    # Pumping from the facility
    energy_perm3_fromthefacility = functions.pumping_per_m3(
        rhowater, 1, pumpefficiency)

    # Water pumped from the well
    initialpumping = facilityvolume*energy_perm3_fromthewell

    pumpingforcleaning = (cleaningvolumeVSfacilityvolume
                          * facilityvolume
                          * energy_perm3_fromthewell)

    # Energy consumption for pumping water during the cultivation
    pumping_during_cultiv = (totalwaterpumpedfromthefacility
                             * energy_perm3_fromthefacility
                             + totalwaterpumpedfromthewell
                             * energy_perm3_fromthewell)/1000  # MJ Conversion L to m3

    totalenergypumping = initialpumping+pumpingforcleaning + pumping_during_cultiv

    LCIdict_micro_updated['Electricity pumping kWh PBR'] = ((
        (totalenergypumping/3.6)
        / total_production_harvested_kg_dw)
        * (1/bioact_molec_dbio)
        * (1/(1 - water_after_drying)))  # kWh

    # Glass consumption

    # Assuming a constant wall thickness of 2 mm.
    glass_perm2 = exchange_area * 0.002 # m3 of glass

    glass_volume_perkgmolecule = ((glass_perm2/total_production_harvested_kg_dw) 
                                  * (1/bioact_molec_dbio)
                                  * (1/(1 - water_after_drying))
                                  * 1/(glass_life_expectancy))  # m3

    glass_mass_perkgmolec = glass_volume_perkgmolecule * 2700  # kg # 2700 kg.m-3

    LCIdict_micro_updated['Glass PBR'] = (glass_mass_perkgmolec
                                    *1/(glass_life_expectancy))

    # Drying
    
    water_to_vaporize_perkilo_dbio = ((1/slurry_concentration)
        * (1 - slurry_concentration)
        * (1 - water_after_drying))  # L. kg-1 dbio

    # /1000 for conversion from kJ to MJ and /3.6 from MJ to kWh
    Electricity_drying_perkg = (water_to_vaporize_perkilo_dbio
                                * (Cw + Cp*(100-meantemp_at_harvest_time_cultivationperiod))
                                / (boilerefficiency*1000))/3.6  # kWh.kg dbio-1

    LCIdict_micro_updated['Electricity drying kWh PBR'] = (Electricity_drying_perkg *
        (1/bioact_molec_dbio))  # Scaled up to 1 kg of molecule kWh

    # Water consumption

    initialfilling = facilityvolume  # m3
     
    
    refillingduringcultivation = totalwaterpumpedfromthewell/1000 # m3 

    totalwater = (initialfilling 
                  + refillingduringcultivation 
                  + cleaningvolumeVSfacilityvolume*initialfilling)  # m3
    
    # Water used for cultivation and not for cleaning
    totalwater_cultivation = refillingduringcultivation + initialfilling  # m3

    totalwater_cultivation_perkgmolecule = ((totalwater_cultivation/total_production_harvested_kg_dw) 
                                            * (1/bioact_molec_dbio)
                                            * (1/(1-water_after_drying)))  # m3

    LCIdict_micro_updated['Water(Cultivation) PBR'] = totalwater_cultivation_perkgmolecule # m3

    #Cleaning
    totalwater_cleaning_perkgmolecule = ((cleaningvolumeVSfacilityvolume*initialfilling/total_production_harvested_kg_dw) 
                                            * (1/bioact_molec_dbio)
                                            * (1/(1-water_after_drying)))

    LCIdict_micro_updated['Water Cleaning PBR'] = totalwater_cleaning_perkgmolecule  # m3

    # Wastewater






    # All water used for cultivatiion - what has been vaporized during drying   (scaled per kg molecule)
    
    # / 1000 to convert  water_to_vaporize_perkilo_dbio from L to m3
    totalwater_towaste_perkgmolecule = (totalwater_cultivation_perkgmolecule
                                        - water_to_vaporize_perkilo_dbio
                                        * (1/bioact_molec_dbio) / 1000)  # m3

    
    LCIdict_micro_updated['Wastewater treatment (without electricity) PBR'] =  totalwater_towaste_perkgmolecule
   
    # Electricity separated to adpat the national market
    LCIdict_micro_updated['Electricity wastewater kWh (low voltage) PBR'] = totalwater_towaste_perkgmolecule * electricity_low_voltage_input_per_m3_wastewater
    LCIdict_micro_updated['Electricity wastewater kWh (high voltage) PBR'] = totalwater_towaste_perkgmolecule * electricity_high_voltage_input_per_m3_wastewater
    
    # Not scaled to molecule for easier wastewater concentration calculation
    totalwater_towaste = (totalwater_cultivation 
                          - water_to_vaporize_perkilo_dbio*total_production_harvested_kg_dw 
                          * (1-water_after_drying) / 1000) # m3

    # Average Concentration waste water in biomass
    conc_waste_water_biomass = total_production_loss_kg_dw/totalwater_towaste  # kg.m-3

    # kg.m-3 or g.L-1  Waste nutrient per kg molecule produced/total wastewater per kg molecule produced
    # Includes the elements in the biomass
    conc_waste_water_nutrient_N = ((N_waste/totalwater_towaste_perkgmolecule 
                                    + conc_waste_water_biomass * N_dw)) # kg.m-3 or g.L-1
    

    conc_waste_water_nutrient_P = ((P_waste/totalwater_towaste_perkgmolecule 
                                    + conc_waste_water_biomass * P_dw)) # kg.m-3 or g.L-1
    
    conc_waste_water_nutrient_K = ((K_waste/totalwater_towaste_perkgmolecule 
                                    + conc_waste_water_biomass * K_dw)) # kg.m-3 or g.L-1
    
    conc_waste_water_nutrient_Mg =((Mg_waste/totalwater_towaste_perkgmolecule 
                                    + conc_waste_water_biomass * Mg_dw)) # kg.m-3 or g.L-1
    
    conc_waste_water_nutrient_S = ((S_waste/totalwater_towaste_perkgmolecule 
                                    + conc_waste_water_biomass * S_dw)) # kg.m-3 or g.L-1
    

    # Carbon only in biomass, CO2 is degazed
    conc_waste_water_C = C_dw * conc_waste_water_biomass  # kg.m-3


    # Land

    LCIdict_micro_updated['Land PBR'] = (
        1/total_production_harvested_kg_dw)*(1/bioact_molec_dbio)*(1/(1-water_after_drying))  # m2


    # Cleaning substances

    # Half of the water with 1 substance, half with the other one
    
    #Hypochlorite
    totalhypo = ((cleaningvolumeVSfacilityvolume*facilityvolume)/2) * concentration_hypo  # kg  
    
    #Hydrogen peroxide
    totalhydro = ((cleaningvolumeVSfacilityvolume*facilityvolume)/2) * concentration_hydro  # kg   

    LCIdict_micro_updated['Hydrogen peroxide PBR'] = (
        totalhydro/total_production_harvested_kg_dw)*(1/bioact_molec_dbio)*(1/(1-water_after_drying))  # kg
    
    LCIdict_micro_updated['Hypochlorite PBR'] = (
        totalhypo/total_production_harvested_kg_dw) *(1/bioact_molec_dbio)*(1/(1-water_after_drying))  # kg



    #Extraction and substitution

    if extraction == 'yes':
        # 1 kWh to disrupt 1 kg of microalgal biomass (kg dbio)
        LCIdict_micro_updated['Electricity cell disruption kWh PBR'] = 1 * (1/bioact_molec_dbio)  # kWh.kg-1

        # for anaerobic digestion, we keep the wet biomass 
        bioact_molec_dw = bioact_molec_dbio/(1-water_after_drying)

        # If extraction, then the 
        # remaining biomass composition is changed according to the biochemical class of the extracted molecule

        if biochemicalclass == 'lip':
            
            lip_dbio_after_extract = (lip_dbio - bioact_molec_dbio)/(1-bioact_molec_dbio)
            
            lip_dw_after_extract = (lip_dw - bioact_molec_dw)/(1-bioact_molec_dw)
          
        
            carb_dbio_after_extract = carb_dbio / (1-bioact_molec_dbio)
            
            carb_dw_after_extract = carb_dw /(1-bioact_molec_dw)

            
            prot_dbio_after_extract = prot_dbio / (1-bioact_molec_dbio)
            
            prot_dw_after_extract = prot_dw /(1-bioact_molec_dw)


            ash_dbio_after_extract = ash_dbio / (1-bioact_molec_dbio)
            
            ash_dw_after_extract = ash_dw /(1-bioact_molec_dw)


            water_dbio_after_extract = water_after_drying /(1-bioact_molec_dbio)
            


        elif biochemicalclass == 'carb':
            
            lip_dbio_after_extract = lip_dbio / (1-bioact_molec_dbio)
            
            lip_dw_after_extract = lip_dw/(1-bioact_molec_dw)

            
            carb_dbio_after_extract = (carb_dbio-bioact_molec_dbio) / (1-bioact_molec_dbio)
            
            carb_dw_after_extract = (carb_dw-bioact_molec_dw) /(1-bioact_molec_dw)

            
            prot_dbio_after_extract = prot_dbio / (1-bioact_molec_dbio)
           
            prot_dw_after_extract = prot_dw /(1-bioact_molec_dw)

            
            ash_dbio_after_extract = ash_dbio / (1-bioact_molec_dbio)
            
            ash_dw_after_extract = ash_dw /(1-bioact_molec_dw)
            
            
            water_dbio_after_extract = water_after_drying / (1-bioact_molec_dbio)

        elif biochemicalclass == 'prot':
            
            lip_dbio_after_extract = lip_dbio / (1-bioact_molec_dbio)

            lip_dw_after_extract = lip_dw/(1-bioact_molec_dw)
            
            
            carb_dbio_after_extract = carb_dbio / (1-bioact_molec_dbio)
            
            carb_dw_after_extract = carb_dw /(1-bioact_molec_dw)

            
            prot_dbio_after_extract = (prot_dbio-bioact_molec_dbio) / (1-bioact_molec_dbio)
            
            prot_dw_after_extract = (prot_dw-bioact_molec_dw) /(1-bioact_molec_dw)


            ash_dbio_after_extract = ash_dbio / (1-bioact_molec_dbio)
            
            ash_dw_after_extract = ash_dw /(1-bioact_molec_dw)
           
            water_dbio_after_extract = water_after_drying / (1-bioact_molec_dbio)


        # After extraction, the substitution will occur with the new composition of the biomass
        
        # Call the function which calculates the masses of subsituted fish feed ingredient
        substitution = functions.optimization_for_fishfeed_substitution(fishfeed_table_withNP,
                                                                        lip_dbio_after_extract,
                                                                        prot_dbio_after_extract,
                                                                        carb_dbio_after_extract,
                                                                        water_dbio_after_extract, 
                                                                        ash_dbio_after_extract,
                                                                        incorporation_rate,
                                                                        MJ_kgcarb,
                                                                        MJ_kgprot,
                                                                        MJ_kglip,phospholipid_fraction,elemental_contents)
        # Collect the biochemical profile of the new fish feed after microalgal incorporation
        


    # if no extraction (molecule given to fish directly), 
    #  the biomass composition stays the same (Obsolete)
    else:  
        LCIdict_micro_updated['Electricity cell disruption kWh PBR'] = 0  # kWH.kg-1

        substitution = functions.optimization_for_fishfeed_substitution(fishfeed_table_withNP,
                                                                        lip_dbio,
                                                                        prot_dbio,
                                                                        carb_dbio,
                                                                        water_after_drying,
                                                                        ash_dbio,
                                                                        incorporation_rate,
                                                                        MJ_kgcarb,
                                                                        MJ_kgprot,
                                                                        MJ_kglip,
                                                                        phospholipid_fraction,
                                                                        elemental_contents)




    # Choose the market that the dependent coproducts enter

    if market_for_substitution == 'animal feed':  

        # Without microalgal incorporation, there is no chnage in the biochemical profile
        biochem_profile_feed = biochem_profile_feed_incumbent
        ingredient_profile_feed = ingredient_profile_incumbent
        final_NP_content_feed = N_P_profile_feed_incumbent
        # Model substitution 1  Animal Feed
        # kg #the same subsitution occurs for every kilo
        feedprot_m1 =  substitution[0] * (1/bioact_molec_dbio - 1) # kg
        feedenergy_m1 =  substitution[1] * (1/bioact_molec_dbio - 1)  # MJ

        LCIdict_micro_updated['Feed energy PBR'] = -feedenergy_m1
        LCIdict_micro_updated['Feed protein PBR'] = -feedprot_m1

        optimization_performance = 'No optimization'
        substitution_check = 'No optimization'

    # Model substitution 2 Fish Feed
    
    elif market_for_substitution == 'fish feed':
        
        biochem_profile_feed = substitution[-2]
        ingredient_profile_feed = substitution[2]
        final_NP_content_feed = substitution[-1]

        LCIdict_micro_updated['Feed energy PBR'] = 0  # Do not use Model 1
        LCIdict_micro_updated['Feed protein PBR'] = 0
        

        if extraction == 'yes':  # Then the substituion only takes place with the remaining biomass
        
            # vect_substitution is a list containing the masses of fish feed ingredient substituted by the remaining biomass
            # (1/bioact_molec_dbio-1) = remaining biomass after extraction of the FU : 1kg of molecule
            
            # substitution[5] is the list of masses of fish feed ingredient 
            # replaced by the given biomass composition. in kg. kg dbio-1 (sum =1 kg)

            vect_substitution = substitution[5]*(1/bioact_molec_dbio-1)  # kg

            # evaluation of the optimized recipe
            optimization_performance = substitution[7]
            

            
        else:  # then the molecule incorporated in the biomass takes part in the substutition

            # substitution[5] is the list of masses of fish feed ingredient
            # replaced by the given biomass composition. in kg. kg dbio-1 (sum =1 kg)
            vect_substitution = substitution[5]*(1/bioact_molec_dbio)  # kg
            
            # evaluation of the optimized recipe
            optimization_performance = substitution[7]

        # Adding the ingredients of the fish feed to the LCI for substitution

        substitution_check = 0 # Obsolete
        
        for a in range(0, len(fishfeed_table_withNP['Ingredient'])):
            
            # Ingredients are ranked in the same order in the vector and in the fish feed table
            LCIdict_micro_updated[fishfeed_table_withNP['Ingredient'][a]] = vect_substitution[a]
            
            substitution_check = (substitution_check
                                  + LCIdict_micro_updated[fishfeed_table_withNP['Ingredient'][a]]) #Obsolete


    # biogaz

    elif market_for_substitution == 'anaerobic digestion':
        
        biochem_profile_feed = biochem_profile_feed_incumbent
        ingredient_profile_feed = ingredient_profile_incumbent
        final_NP_content_feed = N_P_profile_feed_incumbent
        
        LCIdict_micro_updated['Feed energy PBR'] = 0  # Do not use Model 1
        LCIdict_micro_updated['Feed protein PBR'] = 0
        LCIdict_micro_updated["Electricity drying kWh PBR"] = 0 # No drying, directly on wet biomass

        remaining_biomass_after_extraction = 1/bioact_molec_dw-1 # kg

        optimization_performance = 'No optimization'
        substitution_check = 'No optimization'




        [N_substituted,
         P_substituted,
            Mg_substituted,
            K_substituted,
            MJ_contained_in_produced_biogas,
            L_produced_biogas,
            kgCO2_biogenic_reemited_total,
            amount_of_upgrading_act]=functions.anaerobic(prot_dw_after_extract,
                                                    lip_dw_after_extract,
                                                    carb_dw_after_extract,
                                                    ash_dw_after_extract,
                                                    fertiliser_substitutability_AD,
                                                    carbon_degradibilidy_AD,
                                                    elemental_contents,
                                                    phospholipid_fraction,
                                                    CH4_LHV)
                                                    
                                                
                        
             
        LCIdict_micro_updated["N fertiliser substitution AD PBR"]   = -remaining_biomass_after_extraction * N_substituted # as N      
        LCIdict_micro_updated["P fertiliser substitution AD PBR"]   = -remaining_biomass_after_extraction * P_substituted/0.4366      
        LCIdict_micro_updated["Mg fertiliser substitution AD PBR"]   = -remaining_biomass_after_extraction * Mg_substituted * (120.4/24.3)
        LCIdict_micro_updated["K fertiliser substitution AD PBR"]   = -remaining_biomass_after_extraction * K_substituted * 1.2   
        
        
        # Inputs to AD

       
        
        LCIdict_micro_updated["Heat market substitution micro"] = -remaining_biomass_after_extraction * MJ_contained_in_produced_biogas
        
        LCIdict_micro_updated["Electricity kWh AD PBR"] = elec_input_per_kg_biomass_AD * remaining_biomass_after_extraction   
        LCIdict_micro_updated["Heat for AD PBR"] = heat_input_per_kg_biomass_AD  * remaining_biomass_after_extraction   


        LCIdict_micro_updated["biogas purification to biomethane by amino washing_micro (no elec)"] = amount_of_upgrading_act * remaining_biomass_after_extraction   
        LCIdict_micro_updated["Electricity biogas upgrading kWh PBR"] = amount_of_upgrading_act * remaining_biomass_after_extraction   
       
        LCIdict_micro_updated["Biogenic_CO2_re-emission micro"] = kgCO2_biogenic_reemited_total * remaining_biomass_after_extraction   






 
    # Yields

    numberofcultivationdays = len(months_suitable_for_cultivation)*30.4 # days
    
    volumetricyield = (total_production_kg_dw 
                       / (facilityvolume*1000*numberofcultivationdays))  # kg.L-1.d-1
    
    surfaceyield = total_production_kg_dw/numberofcultivationdays  # kg.days-1


    # Check mass balance (Obsolete)

    needed_dbio_check = 1/(lipid_af_dw * (1-ash_dw) *
                           (1-water_after_drying) * bioact_molec_dbio)
    
    
    
    
    
    
    

    return [LCIdict_micro_updated,
            surfaceyield,
            volumetricyield,
            optimization_performance,
            needed_dbio_check,
            substitution_check,
            total_production_kg_dw,
            total_production_harvested_kg_dw,
            total_production_loss_kg_dw,
            conc_waste_water_nutrient_N,
            conc_waste_water_nutrient_P,
            conc_waste_water_nutrient_K,
            conc_waste_water_nutrient_Mg,
            conc_waste_water_biomass,
            conc_waste_water_C,
            conc_waste_water_nutrient_S,
            bioact_molec_dbio,
            min_centrifugation_rate_m3_h,
            max_centrifugation_rate_m3_h,
            totalwatercentrifuged,
            tubelength,
            facilityvolume,
            exchange_area,
            totalcooling_thermal,
            Twell,
            biochem_profile_feed,
            ingredient_profile_feed,
            final_NP_content_feed]  # []






@ray.remote
def calculateLCI_1param_parallel(constant_inputs,param_set) :   
        """Function which calculates one LCI and returns the supply vector for
        a set of input parameters.
        The function is built to be called in parallel with ray."""       
    

        print("1LCI")
    
        (Tech_opdict,
         Biodict,
         Locationdict,
         Physicdict,
         Fish_farm_and_compound_effect_dict,
         Production_mix_dict,
         LCIdict_micro,
         months_suitable_for_cultivation,
         fishfeed_table_withNP,
         elemental_contents,
         names_param_op,
         names_param_bio,
         names_param_geo,
         names_param_phy,
         names_param_fish_farm_compound_effect,
         names_param_prod_mix,
         names_param,
         size_nested_sample,
         nested_list_techno_op,
         list_points_grid,
         demand_vector,
         Techno_Matrix_Fish,
         names_values_simu,
         biochem_profile_feed_incumbent,
         index_growth_stages,
         index_feed,
         index_dead_fish,
         index_Nemissions,
         index_Pemissions,
         index_growth_stages_no_filter_laguna,
         index_sludge,
         electricity_low_voltage_input_per_m3_wastewater,
         electricity_high_voltage_input_per_m3_wastewater,
         ingredient_profile_incumbent,
         N_P_profile_feed_incumbent,
         digestibility_list,
         index_roe,
         Dict_FCR_bio,
         ECO_FCR_0,
         ratio_loss_biological_INC_0,
         index_micro_compound,
         index_growing_DK,
         index_300g,
         index_FU,
         list_meth,
         dict_correspondance_techno_growth_stagenames,
         index_growth_stages_to_modif,
         bio_FCR_0,
         index_biogas_updgrade,
         index_N_substitution,
         index_P_substitution,
         index_heat_substitution)=constant_inputs    
    

        # Update the dictionnaries whith the values of the sample

        #print("LCI Calculation")

        # 1 ) Tech_opdict


        for param in Tech_opdict:  
            # We browse the parameters to look for the uncertain ones
            # which need to be updated

            for index in range(len(names_param_op)):
                

                # Looking for the corresponding parameter in the saltelli set
                if names_param[index] == param:
                    
                    # Then it is an uncertain paramater and its value is 
                    # updated with the one from the generated sample
                    
                    Tech_opdict[param] = param_set[index]
                    
        # We will do the same thing for other dictionnaries but there is no need to browse all possible parameters, 
        # just the ones from other dictionnaries which are left.
        
        new_start = len(names_param_op)
        
        #print("NEW STAAAAAAART1",new_start)

        # 2) Biodict

        for param in Biodict:  
            # We browse the parameters to look for the variable ones
            # which need to be updated

            for index in range(new_start, new_start+len(names_param_bio)):

                # Looking for the corresponding parameter in the saltelli set
                if names_param[index] == param:

                    Biodict[param] = param_set[index]

        new_start = new_start+len(names_param_bio)

        # 3) Locationdict

        for param in Locationdict:  
            # We browse the parameters to look for the variable ones
            # that need to be updated
            for index in range(new_start, new_start+len(names_param_geo)):
                #print("LOCATION DICT 3",index)

                # Looking for the corresponding parameter in the saltelli set
                if names_param[index] == param:
                    #print(names_param[index])
                    Locationdict[param] = param_set[index]

        new_start = new_start+len(names_param_geo)
        #print("NEW STAAAAAAART3",new_start)

        # 4) Physicdict
        
        for param in Physicdict: 
            # We browse the parameters to look for the variable ones
            # that need to be updated

            for index in range(new_start, new_start+len(names_param_phy)):

                # Looking for the corresponding parameter in the saltelli set
                if names_param[index] == param:

                    Physicdict[param] = param_set[index]

        new_start = new_start+len(names_param_phy)


        # 5) Fish_farm_and_compound_effect_dict
        
        for param in Fish_farm_and_compound_effect_dict: 
            # We browse the parameters to look for the variable ones
            # that need to be updated

            for index in range(new_start, new_start+len(names_param_fish_farm_compound_effect)):

                # Looking for the corresponding parameter in the saltelli set
                if names_param[index] == param:

                    Fish_farm_and_compound_effect_dict[param] = param_set[index]

        new_start = new_start+len(names_param_fish_farm_compound_effect)
        
         
        # 6) Production_mix_dict
        
        for param in Production_mix_dict: 
            # We browse the parameters to look for the variable ones
            # that need to be updated

            for index in range(new_start, new_start+len(names_param_prod_mix)):

                # Looking for the corresponding parameter in the saltelli set
                if names_param[index] == param:

                    Production_mix_dict[param] = param_set[index]

        new_start = new_start+len(names_param_prod_mix)
               

        """Retrieve value substitution market"""
        
        microalga_incorporation_True_False = (1/3)<Biodict["random_market_subst"]<(2/3) # True if incorporation
        


        """ Nested sampling"""
        
        # For tech_opdict, the values for parameters corresponding to a nested pert distributions  correspond to the mode for a production mix
        
        
        # Sample Prod Mix
        
        
        min_lat= min([coord[0][1] for coord in list_points_grid])  # get minimum latitude
        max_lat= max([coord[0][1] for coord in list_points_grid])  # get maximum latitude
        
    
        

        sample_lat = list(pert(min_lat,Production_mix_dict['mode_lat'], max_lat, size=int(Production_mix_dict['mix_size'])))

        # Create the production mix

        # Get the closest points in the list of pointe
        #listpoints = [[(59,20),"UK"],(20,12,"FR"),...]
        
        prod_mix = []
        for lat in sample_lat:
            

            list_dif = [(point[0][1]-lat)**2 for point in list_points_grid]

            index_selected = list_dif.index(min(list_dif))
            
            point_selected = list_points_grid.pop(index_selected)
            
            
            prod_mix.append(point_selected)
                            


    
        # Initialize array for median results of each loc in the mix    
        res_prod_mix_micro = np.array([[0]*(len(nested_list_techno_op + names_values_simu)+len(LCIdict_micro))]*len(prod_mix),dtype="float32")
        
        res_prod_mix_micro_wastewater = np.array([[0]*7]*len(prod_mix),dtype="float32")

        
        index_row_in_prod_mix = -1
       
        for loc in prod_mix: # For one location in the production mix
        
            Locationdict['lat'] = loc[0][1]
            Locationdict['long'] = loc[0][0]
            

            
            index_row_in_prod_mix+=1
           
            ## CREATE SAMPLE OF UNCERTAIN PERT NESTED TECHNO-OPERATIONAL PARAMETERS

            sample_nested, names_param_nested = sampling_func_nested(Tech_opdict,
                                                          nested_list_techno_op,
                                                          size_nested_sample)

            # Initialize array for results of Techno_operational_set in the location    

            res_LCIs_Loc = np.array([[0]*(len(nested_list_techno_op + names_values_simu)+len(LCIdict_micro))]*size_nested_sample,dtype="float32")
            
            res_LCIs_Loc_wastewater = np.array([[0]*7]*size_nested_sample,dtype="float32")
            

            index_row_in_loc = -1
            
            for tech_op_set in sample_nested : # For one techno-op setup for a location

                
                 index_row_in_loc+=1
                
                 # Update Techo-opdict ( the only one that changes with nested uncertainty) 
                
                 for index in range(len(names_param_nested)):
                    
                     name = names_param_nested[index]
                     Tech_opdict[name] =  tech_op_set[index]

                
                 # Calculate LCI for this location in this techno-operational setup
                 
                 LCI = LCI_one_strain_uniquevalues(Biodict,
                                      Physicdict,
                                      Tech_opdict,
                                      Locationdict,
                                      LCIdict_micro,
                                      months_suitable_for_cultivation,
                                      fishfeed_table_withNP,
                                      elemental_contents,
                                      biochem_profile_feed_incumbent,
                                      electricity_low_voltage_input_per_m3_wastewater,
                                      electricity_high_voltage_input_per_m3_wastewater,
                                      ingredient_profile_incumbent,
                                      N_P_profile_feed_incumbent)
            
            
                # Collect LCI results
            
                
                 LCIdict_micro_collected = LCI[0]


                 
                 surfaceyield = LCI[1]  # kg dw .d-1
    
                 volumetricyield = LCI[2]  # kg dw.L-1.d-1
                        
                 total_production_kg_dw = LCI[6]  # kg dw
        
                 total_production_harvested_kg_dw = LCI[7]  # kg dw
                
                 conc_waste_water_nutrient_N =LCI[9]
                 conc_waste_water_nutrient_P =LCI[10]
                 conc_waste_water_nutrient_K = LCI[11]
                 conc_waste_water_nutrient_Mg = LCI[12]
                 conc_waste_water_biomass = LCI[13]
                 conc_waste_water_C = LCI[14]
                 conc_waste_water_nutrient_S = LCI[15]
        
                 bioact_molec_dbio = LCI[16]  
        
   
                 tubelength = LCI[20]
        
                 facilityvolume = LCI[21]
                
                 totalcooling_thermal = LCI[23]
        
                 Twell = LCI[24]
                 
                 biochem_profile_feed = LCI[25]
                 
                 ingredient_profile_feed = LCI[26]

                 final_NP_content_feed = LCI[27]
                 
                 # List containing simualtions values which are not LCI or LCIA
                 # (same order as their names in names_suppl_info)
                
                 values_simu = [bioact_molec_dbio,
                               surfaceyield,
                               tubelength,
                               facilityvolume,
                               totalcooling_thermal,
                               volumetricyield,
                               total_production_kg_dw,
                               total_production_harvested_kg_dw]
                 
                 
                 # List containing wastewater concentrations
                 
                 values_wastewater = [conc_waste_water_nutrient_N,
                                      conc_waste_water_nutrient_P,
                                      conc_waste_water_nutrient_K,
                                      conc_waste_water_nutrient_Mg,
                                      conc_waste_water_biomass,
                                      conc_waste_water_C,
                                      conc_waste_water_nutrient_S]
                 
                 # (same order as their names in names_suppl_info)
                 values_LCI_micro = [LCIdict_micro_collected[i] for i in LCIdict_micro_collected]
                
                 res_LCIs_Loc[index_row_in_loc] = list(tech_op_set) + list(values_LCI_micro[:-4]) + list(values_simu)   # remove the strings  # Delete the last 4 string values added to the dictionnary
                 
                 res_LCIs_Loc_wastewater[index_row_in_loc] = list(values_wastewater)
                 



            # calculate median of the LCIs for the loc. 
            # This median becomes the the LCi for the loc (nested utecnho-operational uncertainty)


            res_prod_mix_micro[index_row_in_prod_mix] = np.median(res_LCIs_Loc, axis=0)

            res_prod_mix_micro_wastewater[index_row_in_prod_mix] = np.median(res_LCIs_Loc_wastewater, axis=0)
            

        """Calculate the LCI for the fish farms"""
        
        
        # Update values of the parameters for the fish farms and compound effect
        
        

        """ Collect necessary values"""
        
        
        CH4_LHV = Physicdict["CH4_LHV"]
        


        HAIT_FCR_red_ratio_frac = Fish_farm_and_compound_effect_dict["HAIT_FCR_red_ratio_frac"]
        FRIT_FCR_red_ratio_frac = Fish_farm_and_compound_effect_dict["FRIT_FCR_red_ratio_frac"]
        GOIT1_FCR_red_ratio_frac = Fish_farm_and_compound_effect_dict["GOIT1_FCR_red_ratio_frac"]
        GOIT2_FCR_red_ratio_frac = Fish_farm_and_compound_effect_dict["GOIT2_FCR_red_ratio_frac"]
        GOIT1bis_FCR_red_ratio_frac = Fish_farm_and_compound_effect_dict["GOIT1bis_FCR_red_ratio_frac"]
        GODK_FCR_red_ratio_frac = Fish_farm_and_compound_effect_dict["GODK_FCR_red_ratio_frac"]
        SFDK1_FCR_red_ratio_frac = Fish_farm_and_compound_effect_dict["SFDK1_FCR_red_ratio_frac"]
        SFDK2_FCR_red_ratio_frac = Fish_farm_and_compound_effect_dict["SFDK2_FCR_red_ratio_frac"]


        ratio_CO2_CH4_biogas_fish = Fish_farm_and_compound_effect_dict["ratio_CO2_CH4_biogas_fish"] 
        CH4_volume_biogas_fish = Fish_farm_and_compound_effect_dict["CH4_volume_biogas_fish"]
        P_in_fish = Fish_farm_and_compound_effect_dict["P_in_fish"]

        Mg_in_fish =Fish_farm_and_compound_effect_dict["Mg_in_fish"]
        water_in_fish = Fish_farm_and_compound_effect_dict["water_in_fish"]
        gvs_gts_in_fish = Fish_farm_and_compound_effect_dict["gvs_gts_in_fish"]
        
        fraction_non_ingested = Fish_farm_and_compound_effect_dict["fraction_non_ingested"]


        
        # Fish waste management
        [MJ_substituted_per_kilo_deadfish,
         amount_of_upgrading_act, 
         biogenic_CO2_emitted] = functions.natural_gas_subst(
                    ratio_CO2_CH4_biogas_fish,
                    CH4_volume_biogas_fish,
                    CH4_LHV,
                    water_in_fish,
                    gvs_gts_in_fish)
        
                                                             
        # Calculate the minimum FCR with constant digestibility
        
        min_FCR = P_in_fish/(N_P_profile_feed_incumbent[1]*(digestibility_list[-1]-fraction_non_ingested))

        min_HAIT_FCR_red_ratio = min_FCR/Dict_FCR_bio['HAIT_FCR'] 
        min_FRIT_FCR_red_ratio = min_FCR/Dict_FCR_bio['FRIT_FCR'] 
        min_GOIT1bis_FCR_red_ratio = min_FCR/Dict_FCR_bio['GOIT1bis_FCR'] 
        min_GODK_FCR_red_ratio = min_FCR/Dict_FCR_bio['GODK_FCR'] 
        min_SFDK1_FCR_red_ratio = min_FCR/Dict_FCR_bio['SFDK1_FCR'] 
        min_SFDK2_FCR_red_ratio = min_FCR/Dict_FCR_bio['SFDK2_FCR'] 
        
        # Calculate the FCR reduction ratios according to the input parameters and the minimum FCRs


        HAIT_FCR_red_ratio =  (min_HAIT_FCR_red_ratio-1)* HAIT_FCR_red_ratio_frac + 1
        FRIT_FCR_red_ratio =  (min_FRIT_FCR_red_ratio-1)* FRIT_FCR_red_ratio_frac + 1
        GOIT1bis_FCR_red_ratio =  (min_GOIT1bis_FCR_red_ratio-1)* GOIT1bis_FCR_red_ratio_frac + 1
        GODK_FCR_red_ratio =  (min_GODK_FCR_red_ratio-1)* GODK_FCR_red_ratio_frac + 1
        SFDK1_FCR_red_ratio =  (min_SFDK1_FCR_red_ratio-1)* SFDK1_FCR_red_ratio_frac + 1
        SFDK2_FCR_red_ratio =  (min_SFDK2_FCR_red_ratio-1)* SFDK2_FCR_red_ratio_frac + 1
        

        # Store them in a new dictionnary

        dict_param_fish_farms_ready = Fish_farm_and_compound_effect_dict.copy() 
        
        dict_param_fish_farms_ready["HAIT_FCR_red_ratio"] = HAIT_FCR_red_ratio
        dict_param_fish_farms_ready["FRIT_FCR_red_ratio"]  = FRIT_FCR_red_ratio
        dict_param_fish_farms_ready["GOIT1bis_FCR_red_ratio"]  = GOIT1bis_FCR_red_ratio
        dict_param_fish_farms_ready["GODK_FCR_red_ratio"] = GODK_FCR_red_ratio
        dict_param_fish_farms_ready["SFDK1_FCR_red_ratio"]  = SFDK1_FCR_red_ratio
        dict_param_fish_farms_ready["SFDK2_FCR_red_ratio"] = SFDK2_FCR_red_ratio
        

        # Update technosphere with the new parameters

        Tech_num_AH_modif_with_excretion,Tech_num_INC_modif_with_excretion = techno_modif.calculate_fish_technosphere_matrix(Techno_Matrix_Fish,
                                            dict_param_fish_farms_ready,
                                            dict_correspondance_techno_growth_stagenames,
                                            index_growth_stages_to_modif,
                                            index_growth_stages,
                                            index_dead_fish,
                                            index_roe,
                                            MJ_substituted_per_kilo_deadfish,
                                            amount_of_upgrading_act,
                                            biogenic_CO2_emitted,
                                            N_P_profile_feed_incumbent,
                                            biochem_profile_feed_incumbent,
                                            fraction_non_ingested,
                                            digestibility_list,
                                            CH4_LHV,
                                            index_growth_stages_no_filter_laguna,
                                            index_micro_compound,
                                            index_feed,
                                            index_Nemissions,
                                            index_Pemissions,
                                            index_sludge,
                                            index_biogas_updgrade,
                                            index_N_substitution,
                                            index_P_substitution,
                                            index_heat_substitution)


         # Supply vector in the foreground for the technological concept (with microalgae)          

        supply_vector_tot_AH = linalg.solve(Tech_num_AH_modif_with_excretion, demand_vector)
        
        
         # Supply vector in the foreground for the alternative (without microalgae)          
        
        
        supply_vector_tot_INC = linalg.solve(Tech_num_INC_modif_with_excretion, demand_vector)
        
                
        # Activities in supply vectors are ranked as in list_ acitviites -Names

        
        """ Calculate Performance indicators """
        
        [_,
        ECO_FCR_INC_1,
        bio_FCR_INC_1,
        Loss_fraction_INC_1,
        ratio_loss_biological_INC_1,
        _]=calculate_economic_indicators(Tech_num_INC_modif_with_excretion,
                                  index_300g,
                                  index_growing_DK,
                                  index_roe,
                                  index_FU,
                                  index_feed,
                                  index_dead_fish,
                                  index_micro_compound)
                                                                   
        [_,
        ECO_FCR_AH_1,
        bio_FCR_AH_1,
        Loss_fraction_AH_1,
        ratio_loss_biological_AH_1,
        dose_per_total_bio_AH_1]=calculate_economic_indicators(Tech_num_AH_modif_with_excretion,
                                  index_300g,
                                  index_growing_DK,
                                  index_roe,
                                  index_FU,
                                  index_feed,
                                  index_dead_fish,
                                  index_micro_compound)  
                                                      
        # The following lines define several indicators in addition to the ones used in the article                                                       
                                                               
        #    Increase in fraction of losses in the farm 
        #  Before application of the microalgal compound (potential for improvement)    

        Loss_fraction_increase_INC_indic =   ratio_loss_biological_INC_1/ ratio_loss_biological_INC_0                                                             

        #  After application of the microalgal compound                               
        Loss_fraction_increase_AH_indic =   ratio_loss_biological_AH_1/ ratio_loss_biological_INC_0                                                             
        
        #  Total input microalgal compound per FU        
        Input_micro_comp_per_FU_indic = supply_vector_tot_AH[index_micro_compound]  # in kg compound/kg FU
        

        Economic_FCR_red_0_indic = (ECO_FCR_0-ECO_FCR_AH_1)/ECO_FCR_0
        
        Economic_FCR_red_1_indic = (ECO_FCR_INC_1-ECO_FCR_AH_1)/ECO_FCR_INC_1
        
        Biological_FCR_red_0_indic = (bio_FCR_INC_1-bio_FCR_AH_1)/bio_FCR_0
        
        Biological_FCR_red_1_indic = (bio_FCR_INC_1-bio_FCR_AH_1)/bio_FCR_INC_1

        
        # Economic and biological FCR reduction rate with microalgal compound
        if Input_micro_comp_per_FU_indic!=0:
        
            Economic_FCR_red_rate_1 = Economic_FCR_red_1_indic/(Input_micro_comp_per_FU_indic*1000)
            Biological_FCR_red_rate_1 = Biological_FCR_red_1_indic/(Input_micro_comp_per_FU_indic*1000)
            
            Economic_FCR_red_rate_1_tot_bio = Economic_FCR_red_1_indic/(dose_per_total_bio_AH_1*1000)
            Biological_FCR_red_rate_1_tot_bio = Biological_FCR_red_1_indic/(dose_per_total_bio_AH_1*1000)
            
        else:
            Economic_FCR_red_rate_1=-1000 
            Biological_FCR_red_rate_1 = -1000
            Economic_FCR_red_rate_1_tot_bio = -1000
            Biological_FCR_red_rate_1_tot_bio = -1000
            
        # Impact drugs
        list_impacts_drugs=[]
        for meth in list_meth:
            
            name_param = "impact_Drug_prod" + meth[-1]
            list_impacts_drugs.append(Fish_farm_and_compound_effect_dict[name_param])
    
        # Only some indicators were eventually used in the article
        performance_indicators= [Loss_fraction_increase_INC_indic,
                                Loss_fraction_increase_AH_indic,
                                Input_micro_comp_per_FU_indic,
                                Economic_FCR_red_0_indic,
                                Economic_FCR_red_1_indic,
                                Economic_FCR_red_rate_1,
                                ECO_FCR_AH_1,
                                Biological_FCR_red_1_indic,
                                Biological_FCR_red_0_indic,
                                Biological_FCR_red_rate_1,
                                dose_per_total_bio_AH_1,
                                Economic_FCR_red_rate_1_tot_bio,
                                Biological_FCR_red_rate_1_tot_bio]+list_impacts_drugs
        
        
        return res_prod_mix_micro, res_prod_mix_micro_wastewater, supply_vector_tot_AH,supply_vector_tot_INC,prod_mix, param_set, biochem_profile_feed,performance_indicators


















    

def simulations_fish_micro(Tech_opdict_distributions,
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
                           months_suitable_for_cultivation,
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
                           index_heat_substitution):
    
    """Function which simulates all the LCIs and when done, the LCIAs.
    Returns the raw output table with all data points, 
    and the list of table for contribution analysis.
    
    Divides the work into different tasks performed in parallel"""  
  
    

    # Generate microalgae sample
    
    timeA=time()

    
    output_sample_fish_micro = sampling_func_total_montecarlo(Tech_opdict_distributions,
                  Biodict_distributions,
                  Locationdict_distributions, 
                  Physicdict_distributions,
                  Fish_farm_and_compound_effect_dict_distributions,
                  Production_mix_dict_distributions,
                  size)
    

    
    sample_total_input=output_sample_fish_micro[0]
    
    
    names_param_set=output_sample_fish_micro[1]
    
    
    names_param_op=output_sample_fish_micro[2]
    
    names_param_bio=output_sample_fish_micro[3]
    
    names_param_geo=output_sample_fish_micro[4]
    
    names_param_phy=output_sample_fish_micro[5]
    
    names_param_fish_farm_compound_effect=output_sample_fish_micro[6]
    
    names_param_prod_mix=output_sample_fish_micro[7]
        
    nested_list_techno_op=output_sample_fish_micro[8]
    


    names_values_simu = ['bioact_molec_dbio',
           'surfaceyield',
           'tubelength',
           'facilityvolume',
           'totalcooling_thermal',
           'volumetricyield',
           'total_production_kg_dw',
           'total_production_harvested_kg_dw']

    
    time1=time()
    
    """ Calculate LCIs in parallel"""

    # Start Ray.
    ray.shutdown()
    
    ray.init()
    
    # Inputs common to all tasks
    constant_inputs = ray.put([Tech_opdict,
                             Biodict,
                             Locationdict,
                             Physicdict,
                             Fish_farm_and_compound_effect_dict,
                             Production_mix_dict,
                             LCIdict_micro,
                             months_suitable_for_cultivation,
                             fishfeed_table_withNP,
                             elemental_contents,
                             names_param_op,
                             names_param_bio,
                             names_param_geo,
                             names_param_phy,
                             names_param_fish_farm_compound_effect,
                             names_param_prod_mix,
                             names_param_set,
                             size_nested_sample,
                             nested_list_techno_op,
                             list_points_grid,
                             demand_vector,
                             Techno_Matrix_Fish,
                             names_values_simu,
                             biochem_profile_feed_incumbent,
                             index_growth_stages,
                             index_feed,
                             index_dead_fish,
                             index_Nemissions,
                             index_Pemissions,
                             index_growth_stages_no_filter_laguna,
                             index_sludge,
                             electricity_low_voltage_input_per_m3_wastewater,
                             electricity_high_voltage_input_per_m3_wastewater,
                             ingredient_profile_incumbent,
                             N_P_profile_feed_incumbent,
                             digestibility_list,
                             index_roe,
                             Dict_FCR_bio,
                             ECO_FCR_0,
                             ratio_loss_biological_INC_0,
                             index_micro_compound,
                             index_growing_DK,
                             index_300g,
                             index_FU,
                             list_meth,
                             dict_correspondance_techno_growth_stagenames,
                             index_growth_stages_to_modif,
                             bio_FCR_0,
                             index_biogas_updgrade,
                             index_N_substitution,
                             index_P_substitution,
                             index_heat_substitution]) 

  
    # Calculate all the LCIs (foreground supply vectors)
    arrayresult_raw =ray.get([calculateLCI_1param_parallel.remote(constant_inputs,
                                                                  param_set) for param_set in sample_total_input])        

   
    time_B=time()
    time_LCI =time_B-timeA
    print("time_LCI",time_LCI)


    """ Calculate LCIAs in parallel by combining LCis with mc results for the background"""
    
    list_micro_algae_names = [a for a in LCIdict_micro]
    
    # Inputs common to all tasks
    constant_inputs_LCIA = ray.put([list_micro_algae_names,
                                    filter_names_ok_activities_fish,
                                    list_FU_combined_names_mc,
                                    activities_fish_background,
                                    list_processes_electricity_micro,
                                    list_meth,
                                    Dict_incumbent_losses_growth_stages_loss_level,
                                    Dict_incumbent_losses_growth_stages_loss_red,
                                    Dict_incumbent_outputs_growth_stages_loss_red,
                                    Dict_incumbent_outputs_growth_stages_loss_level,
                                    names_param_set,
                                    Fish_farm_and_compound_effect_dict,
                                    nested_list_techno_op,
                                    list_cfs,
                                    drug_inputs_names])
    


    # Calculate all the LCIAs 
    arrayresult_LCIA = ray.get([LCIA_parallel.remote(constant_inputs_LCIA,
                                                                  row_LCI,row_mc) for row_LCI,row_mc in zip(arrayresult_raw,list_array_total_mc_sorted)])
    
    
    ray.shutdown()
    
    print("Done with LCIAS")
    time_LCIA=time()-time_B
    print("time_LCIA",time_LCIA)
    

    # Rebuild a proper dataframe
    
    # Separate the contributions and the absolute results
    
    table_contrib=[pd.DataFrame(np.zeros([len(arrayresult_LCIA),len(activities_fish_background)]),columns=activities_fish_background)for meth in range(len(list_meth))]

    
    list_result_LCIA_without_contrib =[]



    print("Now sorting results")
    
    for row_index in range(len(arrayresult_LCIA)):
        
        row_LCIA =arrayresult_LCIA[row_index]  # Collect the row containing the LCIA
        

        row_performance_indicators = arrayresult_raw[row_index][-1]  # Collect the row containing the performance indicators
        
        list_result_LCIA_without_contrib.append(row_LCIA[0]+row_performance_indicators)


        
        for index_meth in range(len(list_meth)):
            
            table_contrib[index_meth] = row_LCIA[-1][index_meth]            


    # Collect the rest

    
    activities_fish_background_inc = [name+"_inc" for name in activities_fish_background]
    
    list_meth_short = [meth[-1] for meth in list_meth]
    
    names_total_INC = [meth+"INC" for meth in list_meth_short]
    names_total_AH = [meth+"AH" for meth in list_meth_short]
    names_ratio = [meth+"AH/INC" for meth in list_meth_short]
    


    names_nested = [a[0]+"mean" for a in nested_list_techno_op ]
    
    
    # Columns'names for the performance indicators
    list_name_impact_drugs=[]
    
    for meth in list_meth:
        name= "impact_drug_"+meth[-1]
        
        list_name_impact_drugs.append(name)

    name_indicators =["Loss_fraction_increase_INC_indic",
                                "Loss_fraction_increase_AH_indic",
                                "Input_micro_comp_per_FU_indic",
                                "Economic_FCR_red_0_indic",
                                "Economic_FCR_red_1_indic",
                                "Economic_FCR_red_rate_1",
                                "ECO_FCR_AH_1",
                                "Biological_FCR_red_1_indic",
                                "Biological_FCR_red_0_indic",
                                "Biological_FCR_red_rate_1",
                                "dose_per_total_bio_AH_1",
                                "Economic_FCR_red_rate_1_tot_bio",
                                "Biological_FCR_red_rate_1_tot_bio"]+list_name_impact_drugs
 


    names_col_dataframe = names_param_set + names_nested + list_micro_algae_names + names_values_simu  + names_total_AH + names_total_INC + names_ratio + name_indicators




    results_table_df = pd.DataFrame(np.array(list_result_LCIA_without_contrib), columns=names_col_dataframe)    


    return results_table_df, table_contrib
    
    
    
    
@ray.remote
def LCIA_parallel(constant_inputs_LCIA, 
                  row_LCI, row_mc):
    
    """Function which calculates LC impacts associated with one foreground supply vector and background MC sample. 
    To be called in parallel"""  
    
    
    
    [list_micro_algae_names,
     filter_names_ok_activities_fish,
     list_FU_combined_names_mc,
     activities_fish_background,
     list_processes_electricity_micro,
     list_meth,
     Dict_incumbent_losses_growth_stages_loss_level,
     Dict_incumbent_losses_growth_stages_loss_red,
     Dict_incumbent_outputs_growth_stages_loss_red,
     Dict_incumbent_outputs_growth_stages_loss_level,
     names_param_set,
     Fish_farm_and_compound_effect_dict,
     nested_list_techno_op,
     list_cfs,
     drug_inputs_names]= constant_inputs_LCIA
    
    
    print("1LCIA")
    
    # Composition of row_LCI 11.10 :
    # res_prod_mix_micro, res_prod_mix_micro_wastewater, supply_vector_tot_AH,supply_vector_tot_INC,prod_mix, param_set, biochem_profile_feed,performance_indicators

    # Calculate mean LCI and paramters in the production mix for the microalgal compound
    Micro_result = row_LCI[0]
        
    Micro_result_wastewater = row_LCI[1]
    supply_vector_AH = row_LCI[2]
    
    supply_vector_incumbent = row_LCI[3]

    list_points_in_prod_mix = row_LCI[4]
    
    
    param_set = row_LCI[5]
    
    # Collect the biochemical profile of the feed , which depends on if microalgl incrorporation is assumed or not
    biochem_profile_feed = row_LCI[6]
    
    # We will keep the mean values for the simulaito values in the production mix
    mean_results_LCI_micro = np.mean(Micro_result,axis=0)

    # names mean_results_LCI_micro : nested_list_techno_op + list_micro_algae_names + names_values_simu

    # As the production is equally shared in the mix, the emissions from wastewater treatment are averaged for the mix.
        
    mean_results_micro_wastewater = np.mean(Micro_result_wastewater,axis=0)

    
  
    Micro_result_LCI_only = Micro_result[:,len(nested_list_techno_op):(len(nested_list_techno_op)+len(list_micro_algae_names))]
      


    # Keep only the LCI of the  activities directly connected to a background activity
    LCI_fish_AH = list(compress(supply_vector_AH, filter_names_ok_activities_fish))  # Concept (with micro)
    LCI_fish_INC = list(compress(supply_vector_incumbent, filter_names_ok_activities_fish)) # Alternative (without micro)
    
    
    # These activities are in the same order as their names in "activities_fish_background"
    
    
    # Put the mc results for 1 unit of each activity of the foreground
    # in a dictionnary with the name of the activity as the key
    Dict_mc_FU = {}
    
    # Prepare dictionnary such as { "Actname": [impact for this mc with meth1,impact for this mc with meth2],
    #                              "Actname2:"[impact for this mc with meth1,impact for this mc with meth2]]}
        


    for index_name_FU in range(len(row_mc[0])):  
        
        name_FU = list_FU_combined_names_mc[index_name_FU]
        
        
        Dict_mc_FU[name_FU] = [row_mc[meth_index][index_name_FU] for meth_index in range(len(row_mc))]
        
        

    # Put the paramters valus in a dictionnary with paramter name as key
    
    Dict_param_set_values = {}
    
    for index_name_param in range(len(param_set)): 
        Dict_param_set_values[names_param_set[index_name_param]] = param_set[index_name_param]
    
    
    """ Impact of MICROALGAL COMPOUND """
    # Impact of the microalgae compound production mix     
      
    # Initialize with 0
    table_LCIAs_prod_mix = [np.zeros([len(list_points_in_prod_mix), len(list_micro_algae_names)]) for meth_index in range(len(list_meth))]

    """Impact per m3 wastewater in the prodution mix"""

    # Calculate impact due to emissions after wastewater treatment for the microalgal compound production
    conc_waste_water_nutrient_N = mean_results_micro_wastewater[0]
    conc_waste_water_nutrient_P = mean_results_micro_wastewater[1]
    conc_waste_water_nutrient_K = mean_results_micro_wastewater[2]
    conc_waste_water_nutrient_Mg = mean_results_micro_wastewater[3]
    conc_waste_water_biomass = mean_results_micro_wastewater[4]
    conc_waste_water_C = mean_results_micro_wastewater[5]
    conc_waste_water_S = mean_results_micro_wastewater[6]
    
    
    emissions_impacts_list_1m3_waste_water = waste_water_impact_biosphere(conc_waste_water_nutrient_N,
                          conc_waste_water_nutrient_P,
                          conc_waste_water_C,
                          conc_waste_water_nutrient_Mg,
                          conc_waste_water_nutrient_K,
                          conc_waste_water_S,
                          list_meth,
                          list_cfs)

    Dict_mc_FU["Wastewater treatment (without electricity) PBR"] =[a+b for a,b in zip(Dict_mc_FU["Wastewater treatment (without electricity) PBR"],
                                                                   emissions_impacts_list_1m3_waste_water)]
    


    # For each loc in the production mix
    for loc_index in range(len(list_points_in_prod_mix)): 
        
        point = list_points_in_prod_mix[loc_index]
      
        # Regionalize electricity and gas
        for FU in list_processes_electricity_micro: 
            
            Dict_mc_FU[FU] = Dict_mc_FU[point[1]]   # national Electric mic 
            
     

        # Change natural gas supply to national
        
        Dict_mc_FU["Natural gas micro"] = Dict_mc_FU[point[1]+"gas"]

        # Calulate impact associated to all foreground inputs for Microalgal compound production
        for index_meth in range(len(list_meth)):
            
            for_1_meth=[]
            for index_FU in range(Micro_result_LCI_only.shape[1]):
      
                for_1_meth.append(Micro_result_LCI_only[loc_index,index_FU] * Dict_mc_FU[list_micro_algae_names[index_FU]][index_meth])
                
            table_LCIAs_prod_mix[index_meth][loc_index] = for_1_meth


    # As all locations supply equal shares in the production mix, we take the mean impact across locations  
    LCIA_prod_mix_res = [np.mean(table_meth,axis=0) for table_meth in table_LCIAs_prod_mix] 

    # Sum for total impact for 1 kg of microalgal compound
    # Total impact
    impact_micro_compound_res = [np.sum(table_meth) for table_meth in LCIA_prod_mix_res]
        
    # Update dictionnary of LCIA with the impacts for 1 kg of microalgal compound
    Dict_mc_FU["Micro_compound_prod"] = impact_micro_compound_res
    
    
    
    
    """Calculate Total LCA for the fish FU"""
    
  
    # ADD IMPACT DRUGS (uncertain parameters/factors)
    list_impacts_drugs = []
    for meth in list_meth:
        
        name_param = "impact_Drug_prod" + meth[-1]
        list_impacts_drugs.append(Fish_farm_and_compound_effect_dict[name_param])
    
    # assign impacts
    for name_drug_prod in drug_inputs_names:
        
        Dict_mc_FU[name_drug_prod] = list_impacts_drugs
    

    # Impacts for the fish FU
    LCIAs_AH=[]   # Concept (with Microalgae)
    LCIAs_INC=[]  # Concept (without Microalgae)
    for index_meth in range(len(list_meth)):
        
        for_1_meth_AH=[]
        for_1_meth_INC=[]
        
        for index_FU in range(len(LCI_fish_AH)):
        
            for_1_meth_AH.append(LCI_fish_AH[index_FU] * Dict_mc_FU[activities_fish_background[index_FU]][index_meth])
           
            for_1_meth_INC.append(LCI_fish_INC[index_FU] * Dict_mc_FU[activities_fish_background[index_FU]][index_meth])
            
        LCIAs_AH.append(for_1_meth_AH)
        LCIAs_INC.append(for_1_meth_INC)
            
    Total_impacts_per_kg_fish_AH = [np.sum(list_LCIAs_meth) for list_LCIAs_meth in LCIAs_AH]  # [impact meth1, impact meth 2 ---]
    
    Total_impacts_per_kg_fish_INC = [np.sum(list_LCIAs_meth) for list_LCIAs_meth in LCIAs_INC]   
    
    ratios_AH_INC = [AH/INC for AH,INC in zip(Total_impacts_per_kg_fish_AH,Total_impacts_per_kg_fish_INC)]
   
    list_contributions_AH=[]
    list_contributions_INC=[]
    
    # Contributions
    for index_meth in range(len(list_meth)):
        
        cont_for_1_meth_AH = [impact_input/Total_impacts_per_kg_fish_AH[index_meth] for impact_input in LCIAs_AH[index_meth]]
        
        cont_for_1_meth_INC = [impact_input/Total_impacts_per_kg_fish_INC[index_meth] for impact_input in LCIAs_INC[index_meth]]
        
        list_contributions_AH.append(cont_for_1_meth_AH)
        
        list_contributions_INC.append(cont_for_1_meth_INC)
    
 
    row = list(param_set) + list(mean_results_LCI_micro)+ list(Total_impacts_per_kg_fish_AH) + list(Total_impacts_per_kg_fish_INC) + list(ratios_AH_INC) 
   

    return [row, list_contributions_AH]
    


def waste_water_impact_biosphere(my_conc_N,
                      my_conc_P,
                      my_conc_C,
                      my_conc_Mg,
                      my_conc_K,
                      my_conc_S,
                      methods,
                      list_cfs):
    

    '''Function which :
        -Establishes the mass balances which match the input
        concentrations for a given microalgal waste water.
        
        -Calculates the impact of these emissions for each impact category ( method)

        Inputs :
            
            #my_conc_N, my_conc_P, my_conc_C, my_conc_Mg, my_conc_K, my_conc_S:
                concentrations in g.L-1 of different elements in the actual
                microalgal wastewater entering the treatment
                
            #methods: The list of impact assessment methods.    

        Outputs :
            
            #list_sum_impacts_biosphere_waste_water : list containing the total 
            impacts due to biosphere emissions for the treatment of 1 cubic meter of wastewater.
            1 element per Impact category

        '''
    
    # Defining materials from original wastewater activity
    
    
    # All original values of the output flows classified by element
    # ['ID of the exchange, amount of the exchange']
    
    N_original_flows = [['ae70ca6c-807a-482b-9ddc-e449b4893fe3', 0.00049],
                        ['0017271e-7df5-40bc-833a-36110c1fe5d5', 0.000644],
                        ['6dc1b46f-ee89-4495-95c4-b8a637bcd6cb', 0.0001146],
                        ['d068f3e2-b033-417b-a359-ca4f25da9731', 0.00067568],
                        ['b61057a3-a0bc-4158-882e-b819c4797419', 2.393e-05],
                        ['13331e67-6006-48c4-bdb4-340c12010036', 0.011027],
                        ['9990b51b-7023-4700-bca0-1a32ef921f74', 0.00059438],
                        ['7ce56135-2ca5-4fba-ad52-d62a34bfeb35', 0.048295]]
    
    P_original_flows = [['490b267b-f429-4d9a-ac79-224e37fb4d58', 7.2652e-05],
                        ['1727b41d-377e-43cd-bc01-9eaba946eccb', 0.0027476],
                        ['329fc7d8-4011-4327-84e4-34ff76f0e42d', 2.705e-05],
                        ['198ce8e3-f05a-4bec-9f7f-325347453326', 6.2034e-07]]
    
    C_original_flows = [['f65558fb-61a1-4e48-b4f2-60d62f14b085', 0.0072992],
                        ['8734eb08-50cf-4f5a-8d1a-db76d38efe3c', 4.8293e-05],
                        ['725c7923-0ed8-43e5-b485-fad7e34bef08', 4.8293e-05],
                        ['62859da4-f3c5-417b-a575-8b00d8d658b1', 0.012346],
                        ['73ed05cc-9727-4abf-9516-4b5c0fe54a16', 0.17202],
                        ['960c0f37-f34c-4fc1-b77c-22d8b35fd8d5', 0.0075377],
                        ['9afa0173-ecbd-4f2c-9c5c-b3128a032812', 0.0001613],
                        ['baf58fc9-573c-419c-8c16-831ac03203b9', 0.00050213]]
    
    S_original_flows = [['bfc0bf1c-e5e2-4702-a502-08c892031837', 0.0011039],
                        ['d4049741-cef2-4edd-a3af-9728b9e3a568', 0.0010988],
                        ['8c52f40c-69b7-4538-8923-b371523c71f5', 0.000884],
                        ['37d35fd0-7f07-4b9b-92eb-de3c27050172', 0.14465]]
    
    Mg_original_flows = [['ce9fd912-233a-4807-a33e-0323b1e4a7a2', 0.00014782],
                         ['ebfe261d-ab0d-4ade-8743-183c8c6bdcc6', 2.205e-07],
                         ['e8475907-2081-4fd5-9526-bfcef88380db', 0.00039974],
                         ['7bdab722-11d0-4c42-a099-6f9ed510a44a', 0.0051478]]
    
    K_original_flows = [['1653bf60-f682-4088-b02d-6dc44eae2786', 0.0003989]]
    
    Al_original_flows = [['2baa4381-b781-4f5e-90de-508b0fa3fd1f', 0.0010518],
                         ['97e498ec-f323-4ec6-bcc0-d8a4c853bae3', 6.228e-05],
                         ['01056d4b-f9b0-4dfc-b8d9-8407c8376efb', 0.00031181],
                         ['6f0b8b7c-3888-4174-b7e3-916d42d678ee', 6.5822e-07]]
    
    Na_original_flows = [['1fc409bc-b8e7-48b2-92d5-2ced4aa7bae2', 0.002186]]
    
    Ca_original_flows = [['ac066c02-b403-407b-a1f0-b29ad0f8188f', 0.045852],
                         ['ae28c923-a4a3-4f00-b862-1ae6e748efb9', 0.0012412],
                         ['a912f450-5233-489b-a2e9-8c029fab480f', 2.3777e-06],
                         ['f16fa1da-e426-4820-bf9d-71595c22283b', 0.0035605]]
    
    Fe_original_flows = [['9b6d6f07-ebc6-447d-a3c0-f2017d77d852', 0.0017779],
                         ['7c335b9c-a403-47a8-bb6d-2e7d3c3a230e', 0.0036009],
                         ['db364689-e1a3-4629-8835-e6c59d6daf09', 0.009475],
                         ['32cd0492-c0cb-4898-a2b1-675eedc5b688', 1.2671e-07]]
    Cl_original_flows = [['5e050fab-1837-4c42-b597-ed2f376f768f', 0.040484]]
    
    
    
    # The inputs of different elements in the orginal wastewater treament activity.
    # Added to the waste water as treatment.
    
    added_treatment_N = 2.61388*10**-5
    added_treatment_P = 0
    added_treatment_C = 0
    added_treatment_S = 0.003339284
    added_treatment_Al = 0.000497558
    added_treatment_Fe = 0.0098005
    added_treatment_Na = 9.44323*10**-5
    added_treatment_Ca = 4.17657*10**-7
    added_treatment_Cl = 0.010468958
    added_treatment_Mg = 0
    added_treatment_K = 0
    
    # The total outputs (and thus inputs) of elements in the original activity.
    # Including the added treatments.
    
    totalNoutput = 0.020770454172634983
    totalPoutput = 0.0009270353321544698
    totalCoutput = 0.0746397575218131
    totalSoutput = 0.05012543333333333
    totalAloutput = 0.0014265482200000001
    totalFeoutput = 0.01485392671
    totalNaoutput = 0.002186
    totalCaoutput = 0.050656077699999996
    totalCloutput = 0.040484
    totalMgoutput = 0.0056955805
    totalKoutput = 0.0003989
    
    
    # The actual inputs of elements contained in the waste
    # water of the original activity.
    
    # If the value is negative it means that the mass balance in the original
    # activity was not respected and we assume the element was not in the
    # incoming waste water.
    
    absolute_input_N = max(totalNoutput-added_treatment_N, 0)
    absolute_input_C = max(totalCoutput-added_treatment_C, 0)
    absolute_input_P = max(totalPoutput-added_treatment_P, 0)
    absolute_input_S = max(totalSoutput-added_treatment_S, 0)
    absolute_input_Al = max(totalAloutput-added_treatment_Al, 0)
    absolute_input_Fe = max(totalFeoutput-added_treatment_Fe, 0)
    absolute_input_Na = max(totalNaoutput-added_treatment_Na, 0)
    absolute_input_Ca = max(totalCaoutput-added_treatment_Ca, 0)
    absolute_input_Cl = max(totalCloutput-added_treatment_Cl, 0)
    absolute_input_K = max(totalKoutput-added_treatment_K, 0)
    absolute_input_Mg = max(totalMgoutput-added_treatment_Mg, 0)
    
    
    total_flows=(N_original_flows
                 + P_original_flows
                 + C_original_flows
                 + S_original_flows
                 + Mg_original_flows
                 + K_original_flows
                 + Al_original_flows
                 + Na_original_flows
                 + Ca_original_flows
                 + Fe_original_flows
                 + Cl_original_flows)
    
    

    
    # Initialize the dicitonnary that will contain the impacts associated to each substance in the wastewater
    
    dictionnary_original_flows= {flow[0] : [0]*len(methods) for flow in total_flows} 
    
    # Collect the characterization factors for each impact category
    
    #list_cfs= [bw.Method((meth)).load() for meth in methods]



    meth_index = -1
    
    for meth in methods:
        
        meth_index += 1
        
        # #print(meth)
        # #print(meth_index)

        ##print(list_cfs[meth_index])
        
        cfs_dictionnary = { subst[0][1] : subst[1] for subst in list_cfs[meth_index]}
        
        # For all the biosphere flows in the original waste water activity 
        for flow in total_flows: 
            
            if flow in N_original_flows: # If this flow contains nitrogen
                
                # We assume that  the added treatment is
                # equally shared among the flows of a same element.
                
                original_added_treatment = (flow[1] * added_treatment_N/totalNoutput) 
                
                ##print(original_added_treatment)

                if flow[0] in cfs_dictionnary: # if there is a cf for this flow in this method
        
                    ##print("N",flow)
                    # The impact is cf * new value of the flow in the microalgal wastewater
                    
                    # The new value of the flow is : 
                        
                    # (Original flow - part that comes from the treatment)
                    # * (new concentration waste water/original concentration waste water)
                    # + Share of the treatment ending up in this flow.
                    
                    impact_percubic_meter = cfs_dictionnary[flow[0]] * ((flow[1] - original_added_treatment)
                                                 * my_conc_N/absolute_input_N
                                                 + original_added_treatment)
                    
                    # Update the total impact associated to this flow for the right method

                    dictionnary_original_flows[flow[0]][meth_index] = impact_percubic_meter
                  
            elif flow in P_original_flows:
                
                original_added_treatment = (flow[1] * added_treatment_P/totalPoutput) 
                
                if flow[0] in cfs_dictionnary:
                    # #print("P",flow)

                    impact_percubic_meter = cfs_dictionnary[flow[0]]*((flow[1]-original_added_treatment)
                                                 * my_conc_P/absolute_input_P
                                                 + original_added_treatment)
                    
                    dictionnary_original_flows[flow[0]][meth_index] = impact_percubic_meter
           
            elif flow in C_original_flows:
                
                original_added_treatment = (flow[1] * added_treatment_C/totalCoutput) 
                
                if flow[0] in cfs_dictionnary:
                    
                    impact_percubic_meter = cfs_dictionnary[flow[0]]*((flow[1]-original_added_treatment)
                                                 * my_conc_C/absolute_input_C
                                                 + original_added_treatment)
                    
                    dictionnary_original_flows[flow[0]][meth_index] = impact_percubic_meter


    
            elif flow in S_original_flows:
                
                original_added_treatment = (flow[1] * added_treatment_S/totalSoutput) 
                
                if flow[0] in cfs_dictionnary:
                    
                    impact_percubic_meter = cfs_dictionnary[flow[0]]*((flow[1]-original_added_treatment)
                                                 * my_conc_S/absolute_input_S
                                                 + original_added_treatment)
                    
                    dictionnary_original_flows[flow[0]][meth_index] = impact_percubic_meter


            elif flow in Mg_original_flows:
                
                original_added_treatment = (flow[1] * added_treatment_Mg/totalMgoutput) 
                
                if flow[0] in cfs_dictionnary:
                    
                    impact_percubic_meter = cfs_dictionnary[flow[0]]*((flow[1]-original_added_treatment)
                                                 * my_conc_Mg/absolute_input_Mg
                                                 + original_added_treatment)
                    
                    dictionnary_original_flows[flow[0]][meth_index] = impact_percubic_meter

    
            elif flow in K_original_flows:
                
                original_added_treatment = (flow[1] * added_treatment_K/totalKoutput) 
                
                if flow[0] in cfs_dictionnary:
                    
                    impact_percubic_meter = cfs_dictionnary[flow[0]]*((flow[1]-original_added_treatment)
                                                 * my_conc_K/absolute_input_K
                                                 + original_added_treatment)
                    
                    dictionnary_original_flows[flow[0]][meth_index] = impact_percubic_meter

    
            elif flow in Al_original_flows:
                
                original_added_treatment = (flow[1] * added_treatment_Al/totalAloutput) 
                
                if flow[0] in cfs_dictionnary:
                    
                    impact_percubic_meter = cfs_dictionnary[flow[0]]*((flow[1]-original_added_treatment)
                                                 * my_conc_Al/absolute_input_Al
                                                 + original_added_treatment)
                    
                    dictionnary_original_flows[flow[0]][meth_index] = impact_percubic_meter

                
            elif flow in Na_original_flows:
                
                original_added_treatment = (flow[1] * added_treatment_Na/totalNaoutput) 
                
                if flow[0] in cfs_dictionnary:
                    
                    impact_percubic_meter = cfs_dictionnary[flow[0]]*((flow[1]-original_added_treatment)
                                                 * my_conc_Na/absolute_input_Na
                                                 + original_added_treatment)
                    
                    dictionnary_original_flows[flow[0]][meth_index] = impact_percubic_meter

            elif flow in Ca_original_flows:
                
                original_added_treatment = (flow[1] * added_treatment_Ca/totalCaoutput) 
                
                if flow[0] in cfs_dictionnary:
                    
                    impact_percubic_meter = cfs_dictionnary[flow[0]]*((flow[1]-original_added_treatment)
                                                 * my_conc_Ca/absolute_input_Ca
                                                 + original_added_treatment)
                    
                    dictionnary_original_flows[flow[0]][meth_index] = impact_percubic_meter

            elif flow in Fe_original_flows:
                
                original_added_treatment = (flow[1] * added_treatment_Fe/totalFeoutput) 
                
                if flow[0] in cfs_dictionnary:
                    
                    impact_percubic_meter = cfs_dictionnary[flow[0]]*((flow[1]-original_added_treatment)
                                                 * my_conc_Fe/absolute_input_Fe
                                                 + original_added_treatment)
                    
                    dictionnary_original_flows[flow[0]][meth_index] = impact_percubic_meter

    
            elif flow in Cl_original_flows:
                
                original_added_treatment = (flow[1] * added_treatment_Cl/totalCloutput) 
                
                if flow[0] in cfs_dictionnary:
                    
                    impact_percubic_meter = cfs_dictionnary[flow[0]]*((flow[1]-original_added_treatment)
                                                 * my_conc_Cl/absolute_input_Cl
                                                 + original_added_treatment)
                    
                    dictionnary_original_flows[flow[0]][meth_index] = impact_percubic_meter

        
        # Summing the impacts of all flows for 1 cubic meter
        
        list_sum_impacts_biosphere_waste_water =[]
        
        for meth_index_2 in range(len(methods)):
            
            sum_impact = sum([dictionnary_original_flows[flow][meth_index_2] for flow in dictionnary_original_flows ])


            list_sum_impacts_biosphere_waste_water.append(sum_impact)
            
            
            
            
            
        #wastewater_copy.save()
    
    return list_sum_impacts_biosphere_waste_water



    


def calculate_economic_indicators(Tech_Matrix,
                                  index_300g,
                                  index_growing_DK,
                                  index_roe,
                                  index_FU,
                                  index_feed,
                                  index_dead_fish,
                                  index_micro_compound):
    
    """Function which calculates perfomance indicators associated with a tehcnosphere matrix
    Tech_Matrix = Numerical matrix"""
    
    Tech_Matrix_modif=Tech_Matrix.copy()
    # 1) Add the output of 300 g trout as an economic output
    Tech_Matrix_modif[index_FU,index_growing_DK] = Tech_Matrix_modif[index_FU,index_growing_DK] +Tech_Matrix_modif[index_300g,index_growing_DK]
    
    Tech_Matrix_modif[index_300g,index_growing_DK] = 0
    
    # 2) Add the roe output as a an economic output
    
    Tech_Matrix_modif[index_FU,index_FU] = Tech_Matrix_modif[index_FU,index_FU] + Tech_Matrix_modif[index_roe,index_FU]

    
    Tech_Matrix_modif[index_roe,index_FU] = 0
    
    
    demand_vector_unity = [0 for i in range(len(Tech_Matrix_modif))]
    demand_vector_unity[index_FU] = 1

    supply_vector_tot = linalg.solve(Tech_Matrix_modif, demand_vector_unity)
    
    eco_FCR =  supply_vector_tot[index_feed]    # Feed per kilo economic output = Sea raised trout(FU) + 300g trout + Roe
    Loss_ratio = supply_vector_tot[index_dead_fish]
    ratio_loss_biological = Loss_ratio/(Loss_ratio+1)
    
    # FCR biological 
    bio_FCR =  eco_FCR/(1+Loss_ratio) # Feed per kilo biological output = Sea raised trout(FU) + 300g trout + Roe +Dead
    
    # Dose per total biological production (including dead)
    
    dose_per_live_biomass = supply_vector_tot[index_micro_compound]
    
    dose_per_total_bio = dose_per_live_biomass/(1+Loss_ratio)
    
    return Tech_Matrix_modif,eco_FCR,bio_FCR, Loss_ratio, ratio_loss_biological,dose_per_total_bio

    

    