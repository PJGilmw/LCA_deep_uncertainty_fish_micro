# -*- coding: utf-8 -*-
"""
Created on Fri Oct 14 14:24:06 2022

@author: Pierre Jouannais, Department of Planning, DCEA, Aalborg University
pijo@plan.aau.dk

"""
'''
Script containing the functions which modify the growth stages, 
update the technosphere matrix and calculate the excretion.


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

import ray


import pert


import scipy
import numpy as np
from scipy import linalg

from ast import literal_eval



from itertools import compress


import Functions_for_physical_and_biological_calculations_3rd as functions





             
def calculate_A(FCR_red_ratio,Fish_dead_0,Total_output_fish,Fish_input_0):

    num = Fish_dead_0 + Total_output_fish-Fish_input_0 * (1-FCR_red_ratio) 
    A = num/FCR_red_ratio
    
    return A

def calculate_F_dead_1(A,Fish_dead_0,Fish_input_0,loss_lev,loss_red,Total_output_fish):

    Total_bio = Fish_dead_0 + Total_output_fish

    F_dead1 = A * (1-loss_red) * (Fish_dead_0 + loss_lev * Total_output_fish)/Total_bio
    return F_dead1 

def calculate_F_out_x_1(A,Fish_dead_0,Fish_input_0,loss_lev,loss_red,F_out_0,Total_output_fish):
    
    Total_bio = Fish_dead_0 + Total_output_fish

    F_outx1 = (A/Total_bio) *(F_out_0 *( 1-loss_lev + loss_red * loss_lev + loss_red * Fish_dead_0/Total_output_fish))
    
    return F_outx1





 
def calculate_fish_technosphere_matrix(Techno_Matrix_Fish,
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
                                            index_heat_substitution):
    
    
    
    '''
    Function which updates the tehcnosphere matrix according to the sampled parameters. 
    Modification of the fish flows, sludge and N,P excretion

    
      #Outputs : 
          
          Updated technosphere matrix
         
    '''
    
    P_in_fish=dict_param_fish_farms_ready["P_in_fish"]
    N_in_fish=dict_param_fish_farms_ready["N_in_fish"] 
    K_in_fish=dict_param_fish_farms_ready["K_in_fish"]
    Mg_in_fish=dict_param_fish_farms_ready["Mg_in_fish"]
    water_in_fish=dict_param_fish_farms_ready["water_in_fish"]
    gvs_gts_in_fish=dict_param_fish_farms_ready["gvs_gts_in_fish"] 
        

    excretion_N_removal_efficiency =dict_param_fish_farms_ready["excretion_N_removal_efficiency"]
    excretion_P_removal_efficiency =dict_param_fish_farms_ready["excretion_P_removal_efficiency"]
    solid_filter_efficiency =dict_param_fish_farms_ready["solid_filter_efficiency"]
    
   
    fraction_non_ingested = dict_param_fish_farms_ready["fraction_non_ingested"]
    
    CH4_volume_sludge_and_manure = dict_param_fish_farms_ready["CH4_volume_sludge_and_manure"]

    share_fish_sludge_in_substrate = dict_param_fish_farms_ready["share_fish_sludge_in_substrate"]
     
    ratio_CO2_CH4_biogas_sludge = dict_param_fish_farms_ready["ratio_CO2_CH4_biogas_sludge"]
    N_manure = dict_param_fish_farms_ready["N_manure"]
    P_manure= dict_param_fish_farms_ready["P_manure"]
    
    

    
    fertilizer_substi_digest_N = dict_param_fish_farms_ready["fertilizer_substi_digest_N"]
    fertilizer_substi_digest_P = dict_param_fish_farms_ready["fertilizer_substi_digest_P"]
    fertilizer_substi_digest_K = dict_param_fish_farms_ready["fertilizer_substi_digest_K"]
    fertilizer_substi_digest_Mg = dict_param_fish_farms_ready["fertilizer_substi_digest_Mg"]

    fertilizer_substi_manure_field_N = dict_param_fish_farms_ready["fertilizer_substi_manure_field_N"]
    fertilizer_substi_manure_field_P = dict_param_fish_farms_ready["fertilizer_substi_manure_field_P"]
    

    

    
    
    # Initialize the numerical Tehcnosphere
    
    Techno_fish_num_AH = np.zeros((Techno_Matrix_Fish.shape[0],Techno_Matrix_Fish.shape[0]))
    Techno_fish_num_INC = np.zeros((Techno_Matrix_Fish.shape[0],Techno_Matrix_Fish.shape[0]))
    
    
    for i in range(Techno_fish_num_AH.shape[0]):
    
        for j in range(Techno_fish_num_AH.shape[1]):
            
            if type(Techno_Matrix_Fish[i,j])==str: # Then we evaluate the string
            
                Techno_fish_num_AH[i,j] = eval(Techno_Matrix_Fish[i,j])
                Techno_fish_num_INC[i,j] = eval(Techno_Matrix_Fish[i,j])
                
            else:        # Then we just take the float
                Techno_fish_num_AH[i,j] = Techno_Matrix_Fish[i,j]
                Techno_fish_num_INC[i,j] = Techno_Matrix_Fish[i,j]

    
    

    Tech_num_AH_modif = Techno_fish_num_AH.copy()
    
    Tech_num_INC_modif = Techno_fish_num_INC.copy()
    

    # first modif AH
    for column_index in index_growth_stages_to_modif:
        
        # Associate growth stages numbers and parameters
        
        
        # Code
        code_growth_stage = dict_correspondance_techno_growth_stagenames[column_index]
      
        # Collect the names of the parameters corresponding to the growth stage

        name_corresponding_FCR_red_ratio =  code_growth_stage+"_FCR_red_ratio"
        
        name_corresponding_loss_level =  code_growth_stage+"_loss_lev"
        
        name_corresponding_loss_red =  code_growth_stage+"_loss_red"
        
        name_corresponding_micro_dose =  code_growth_stage+"_micro_dose"
        
        # Collect Values parameters with their names
        FCR_red_ratio = dict_param_fish_farms_ready[name_corresponding_FCR_red_ratio]
        
        FCR_red_ratio_INC = 1 # In the alternative, no modification
        
        
        loss_lev = dict_param_fish_farms_ready[name_corresponding_loss_level]
        
        loss_lev_INC = dict_param_fish_farms_ready[name_corresponding_loss_level]

        loss_red = dict_param_fish_farms_ready[name_corresponding_loss_red]
        
        loss_red_INC = 0
        
        # Micro dose 
        
        micro_dose = dict_param_fish_farms_ready[name_corresponding_micro_dose]
        
        
        # Collect incumbent values in the technosphere matrix
        
        # Incumbent dead
        
        Fish_dead_0 = -Tech_num_AH_modif[index_dead_fish,column_index]
        
        # Incumbent fish outputs and inputs to the growth stage

        dict_fish_outputs = {}
        Fish_input_0=0
        
        for index_input in index_growth_stages+[index_roe]: # The produciton outputs have a "+" and the fish inputs from other stages have a "-". So to callcualte the weight gain, we can just add all input values
           

            if Tech_num_AH_modif[index_input,column_index] < 0:  # then it's a fish input from a previous stage
                
                Fish_input_0 = -Tech_num_AH_modif[index_input,column_index]
                #print(Fish_input_0)

            elif Tech_num_AH_modif[index_input,column_index] > 0:  # It's a fish or roe output
                
                dict_fish_outputs[index_input] = Tech_num_AH_modif[index_input,column_index]


        
        
        Total_output_fish_0= sum(dict_fish_outputs.values()) 
        
        """ Modif the matrix"""


        # New total biomass
        
        A =calculate_A(FCR_red_ratio,Fish_dead_0,Total_output_fish_0,Fish_input_0)
        A_INC =calculate_A(FCR_red_ratio_INC,Fish_dead_0,Total_output_fish_0,Fish_input_0) # In the alternative
        
        
    
        # Dead
        dead_1 = calculate_F_dead_1(A,Fish_dead_0,Fish_input_0,loss_lev,loss_red,Total_output_fish_0)
        dead_1_INC = calculate_F_dead_1(A_INC,Fish_dead_0,Fish_input_0,loss_lev_INC,loss_red_INC,Total_output_fish_0)
        
        #print("dead_1",dead_1)
        #print("dead_1_INC",dead_1_INC)

        Tech_num_AH_modif[index_dead_fish,column_index]  =  -dead_1
        Tech_num_INC_modif[index_dead_fish,column_index]  =  -dead_1_INC
        
        # Fish outputs
        summ=0
        summ_INC=0


   

        for index_input_2 in dict_fish_outputs:
            
            Fish_0 = dict_fish_outputs[index_input_2]
            Fish_1 = calculate_F_out_x_1(A,Fish_dead_0,Fish_input_0,loss_lev,loss_red,Fish_0,Total_output_fish_0)
            Fish_1_INC = calculate_F_out_x_1(A_INC,Fish_dead_0,Fish_input_0,loss_lev_INC,loss_red_INC,Fish_0,Total_output_fish_0)
            
            summ+=Fish_1
            summ_INC+=Fish_1_INC
            
            
            Tech_num_AH_modif[index_input_2,column_index]  =  Fish_1
            
            Tech_num_INC_modif[index_input_2,column_index]  =  Fish_1_INC
            
        
        #print("micro_dose",micro_dose)
        Tech_num_AH_modif[index_micro_compound,column_index] = -micro_dose * summ


        
        feed_input = Tech_num_AH_modif[index_feed,column_index]
        new_bio_FCR_AH = -feed_input/(summ+dead_1-Fish_input_0)
        new_bio_FCR_INC = -feed_input/(summ_INC + dead_1_INC-Fish_input_0)
        new_eco_FCR_AH = -feed_input/(summ-Fish_input_0)
        new_eco_FCR_INC = -feed_input/(summ_INC-Fish_input_0)
        

        
        
        # Calculate excretion

    
    Tech_num_AH_modif_with_excretion,Techno_fish_num_INC_with_excretion = excretion(
              Tech_num_AH_modif,
              Tech_num_INC_modif,
              index_growth_stages,
              N_P_profile_feed_incumbent,
              P_in_fish,
              N_in_fish,
              biochem_profile_feed_incumbent,
              index_dead_fish,
              fraction_non_ingested,
              digestibility_list,
              excretion_P_removal_efficiency,
              excretion_N_removal_efficiency,
              solid_filter_efficiency,
              CH4_volume_sludge_and_manure,
              share_fish_sludge_in_substrate,
              CH4_LHV,
              ratio_CO2_CH4_biogas_sludge,
              N_manure,
              P_manure,
              fertilizer_substi_digest_N,
              fertilizer_substi_digest_P,
              fertilizer_substi_manure_field_N,
              fertilizer_substi_manure_field_P,
              index_growth_stages_no_filter_laguna,
              index_feed,
              index_roe,
              index_Nemissions,
              index_Pemissions,
              index_sludge,
              index_biogas_updgrade,
              index_N_substitution,
              index_P_substitution,
              index_heat_substitution) 
        

     
    # Indicators

    
    [_,
        ECO_FCR_INC_1,
        bio_FCR_INC_1,
        Loss_fraction_INC_1,
        ratio_loss_biological_INC_1]=calculate_economic_indicators(Tech_num_INC_modif,
                                  index_300g,
                                  index_growing_DK,
                                  index_roe,
                                  index_FU,
                                  index_feed,
                                  index_dead_fish)
                                                                   
    [_,
        ECO_FCR_AH_1,
        bio_FCR_AH_1,
        Loss_fraction_AH_1,
        ratio_loss_biological_AH_1]=calculate_economic_indicators(Tech_num_AH_modif,
                                  index_300g,
                                  index_growing_DK,
                                  index_roe,
                                  index_FU,
                                  index_feed,
                                  index_dead_fish)  
    

                                                                  
    return Tech_num_AH_modif_with_excretion,Techno_fish_num_INC_with_excretion
    


    










def excretion(Tech_num_AH_modif,
              Tech_num_INC_modif,
              index_growth_stages,
              N_P_profile_feed_incumbent,
              P_in_fish,
              N_in_fish,
              biochem_profile_feed_incumbent,
              index_dead_fish,
              fraction_non_ingested,
              digestibility_list,
              excretion_P_removal_efficiency,
              excretion_N_removal_efficiency,
              solid_filter_efficiency,
              CH4_volume_sludge_and_manure,
              share_fish_sludge_in_substrate,
              CH4_LHV,
              ratio_CO2_CH4_biogas_sludge,
              N_manure,
              P_manure,
              fertilizer_substi_digest_N,
              fertilizer_substi_digest_P,
              fertilizer_substi_manure_field_N,
              fertilizer_substi_manure_field_P,
              index_growth_stages_no_filter_laguna,
              index_feed,
              index_roe,
              index_Nemissions,
              index_Pemissions,
              index_sludge,
              index_biogas_updgrade,
              index_N_substitution,
              index_P_substitution,
              index_heat_substitution):
 
     
     '''
     Function which calculates the exrection according to the new fish flows in the tehcnosphere matrix
    
        
     #Outputs : 
              
     Updated technosphere matrix with excretion
             
     '''

    
    # AH : with micro
    # INC : without

     Tech_num_AH_modif_with_excretion = Tech_num_AH_modif.copy()
     Tech_num_INC_modif_with_excretion = Tech_num_INC_modif.copy()

     for growth_stage_index in index_growth_stages:
         
         
         # Collect feed input 
         total_feed_AH = - Tech_num_AH_modif_with_excretion[index_feed,growth_stage_index]
         total_feed_INC = -Tech_num_INC_modif_with_excretion[index_feed,growth_stage_index]

         total_weight_gain_economic_AH =  0
         total_weight_gain_economic_INC =  0
         
         # Some growth stages have multiple output of fish, so we make sure to collect all outputs
    
         for index_input in index_growth_stages+[index_roe]: # The production outputs have a "+" and the fish inputs from other stages have a "-". So to callcualte the weight gain, we can just add all input values
             
             total_weight_gain_economic_AH += Tech_num_AH_modif_with_excretion[index_input,growth_stage_index]
             total_weight_gain_economic_INC +=Tech_num_INC_modif_with_excretion[index_input,growth_stage_index]
            
         # Biological weight gain includes the losses ( - as there is a - in the matrix)
         total_weight_gain_biological_AH = total_weight_gain_economic_AH - Tech_num_AH_modif_with_excretion[index_dead_fish,growth_stage_index] 
         total_weight_gain_biological_INC = total_weight_gain_economic_INC -Tech_num_INC_modif_with_excretion[index_dead_fish,growth_stage_index] 
         

      
         # Total N and P in feed

         N_in_feed_AH = total_feed_AH * N_P_profile_feed_incumbent[0]
         P_in_feed_AH = total_feed_AH * N_P_profile_feed_incumbent[1]
         N_in_feed_INC = total_feed_INC * N_P_profile_feed_incumbent[0]
         P_in_feed_INC = total_feed_INC * N_P_profile_feed_incumbent[1]
             
         #Total N and P in fish stock (biological = )
         N_stock_fish_AH = total_weight_gain_biological_AH * N_in_fish
         P_stock_fish_AH = total_weight_gain_biological_AH * P_in_fish
         N_stock_fish_INC = total_weight_gain_biological_INC * N_in_fish
         P_stock_fish_INC = total_weight_gain_biological_INC * P_in_fish    
         # Feed non ingested 

         
         non_ingested_feed_AH = total_feed_AH * fraction_non_ingested
         non_ingested_feed_INC = total_feed_INC * fraction_non_ingested
         

     
         """So far we assume the N and P content modification with microalgae is neglectable"""
         N_in_non_ingested_feed_AH = non_ingested_feed_AH * N_P_profile_feed_incumbent[0]
         P_in_non_ingested_feed_AH = non_ingested_feed_AH * N_P_profile_feed_incumbent[1]
         N_in_non_ingested_feed_INC = non_ingested_feed_INC * N_P_profile_feed_incumbent[0]
         P_in_non_ingested_feed_INC = non_ingested_feed_INC * N_P_profile_feed_incumbent[1]
                    
     

         
         #Solid faeces 
         #Total Biochemical class * digestibility
         
         # keep only Lipid, Protein and add phosphorus

    
         biochem_profile_feed_clean_AH = biochem_profile_feed_incumbent[0:4] + [N_P_profile_feed_incumbent[1]]
         biochem_profile_feed_clean_INC = biochem_profile_feed_incumbent[0:4] + [N_P_profile_feed_incumbent[1]]
         

         total_excreted_per_biochemical_class_AH = [biochem * (1-digest) * total_feed_AH for biochem,digest in zip(biochem_profile_feed_clean_AH,digestibility_list)]
         total_excreted_per_biochemical_class_INC = [biochem * (1-digest) * total_feed_INC for biochem,digest in zip(biochem_profile_feed_clean_INC,digestibility_list)]
    
     
         solid_faeces_AH = sum(total_excreted_per_biochemical_class_AH)
         solid_faeces_INC = sum(total_excreted_per_biochemical_class_INC)
         
         #print(total_excreted_per_biochemical_class_AH,"total_excreted_per_biochemical_class_AH")
         
         
         N_in_solid_faeces_AH = total_excreted_per_biochemical_class_AH[1] * 0.16
         N_in_solid_faeces_INC = total_excreted_per_biochemical_class_INC[1] * 0.16
        
         P_in_solid_faeces_AH = total_excreted_per_biochemical_class_AH[4]
         P_in_solid_faeces_INC = total_excreted_per_biochemical_class_INC[4]
         
         # Total solid
         
         Total_solid_AH =solid_faeces_AH + non_ingested_feed_AH
         Total_solid_INC =solid_faeces_INC + non_ingested_feed_INC

         # Liquid excretion
         
         N_liquid_AH = N_in_feed_AH - N_in_non_ingested_feed_AH - N_stock_fish_AH - N_in_solid_faeces_AH
       
         N_liquid_INC = N_in_feed_INC - N_in_non_ingested_feed_INC - N_stock_fish_INC - N_in_solid_faeces_INC
        
         P_liquid_AH = P_in_feed_AH - P_in_non_ingested_feed_AH - P_stock_fish_AH - P_in_solid_faeces_AH
         P_liquid_INC = P_in_feed_INC - P_in_non_ingested_feed_INC - P_stock_fish_INC - P_in_solid_faeces_INC
         
    
         # Filter and lagune efficiency
         N_liquid_after_treatment_AH = N_liquid_AH * (1-excretion_N_removal_efficiency)
         N_liquid_after_treatment_INC = N_liquid_INC * (1-excretion_N_removal_efficiency)
         P_liquid_after_treatment_AH = P_liquid_AH * (1-excretion_P_removal_efficiency)
         P_liquid_after_treatment_INC = P_liquid_INC * (1-excretion_P_removal_efficiency)
         
         Solid_faeces_after_filter_AH = solid_faeces_AH * (1-solid_filter_efficiency)
         Solid_faeces_after_filter_INC = solid_faeces_INC * (1-solid_filter_efficiency)
     
         Solid_faeces_sludge_AH = solid_faeces_AH * solid_filter_efficiency
         Solid_faeces_sludge_INC = solid_faeces_INC * solid_filter_efficiency
         
         total_sludge_AH =  Solid_faeces_sludge_AH + non_ingested_feed_AH*solid_filter_efficiency
         total_sludge_INC =  Solid_faeces_sludge_INC + non_ingested_feed_INC*solid_filter_efficiency
    
         #
         N_total_without_removal_AH = N_liquid_AH + N_in_solid_faeces_AH + N_in_non_ingested_feed_AH 
         N_total_without_removal_INC = N_liquid_INC + N_in_solid_faeces_INC + N_in_non_ingested_feed_INC
         P_total_without_removal_AH = P_liquid_AH + P_in_solid_faeces_AH + P_in_non_ingested_feed_AH 
         P_total_without_removal_INC = P_liquid_INC + P_in_solid_faeces_INC + P_in_non_ingested_feed_INC 
     
         N_solid_and_liquid_after_removal_AH =  (N_in_solid_faeces_AH + N_in_non_ingested_feed_AH) * (1-solid_filter_efficiency)*(1-excretion_N_removal_efficiency) + N_liquid_after_treatment_AH 
         N_solid_and_liquid_after_removal_INC =  (N_in_solid_faeces_INC + N_in_non_ingested_feed_INC) * (1-solid_filter_efficiency)*(1-excretion_N_removal_efficiency) + N_liquid_after_treatment_INC 
         P_solid_and_liquid_after_removal_AH =  (P_in_solid_faeces_AH + P_in_non_ingested_feed_AH) * (1-solid_filter_efficiency)*(1-excretion_P_removal_efficiency) + P_liquid_after_treatment_AH 
         P_solid_and_liquid_after_removal_INC =  (P_in_solid_faeces_INC + P_in_non_ingested_feed_INC) * (1-solid_filter_efficiency)*(1-excretion_P_removal_efficiency) + P_liquid_after_treatment_INC 
     
         # Content N, P per kilo dm sludge
         N_in_sludge_dw = (N_in_solid_faeces_AH +N_in_non_ingested_feed_AH)*solid_filter_efficiency/total_sludge_AH
         P_in_sludge_dw = (P_in_solid_faeces_AH +P_in_non_ingested_feed_AH)*solid_filter_efficiency/total_sludge_AH
 
         
         # Here update sludge treatment
         
         [fertilizer_subst_N_per_kilo_dw_sludge,
          fertilizer_subst_P_per_kilo_dw_sludge,
          MJ_substituted_per_kilo_dw_sludge,
          amount_of_upgrading_act_sludge]=functions.fish_sludge_management(CH4_volume_sludge_and_manure,
                           share_fish_sludge_in_substrate,
                           CH4_LHV,
                           ratio_CO2_CH4_biogas_sludge,
                           N_in_sludge_dw,
                           P_in_sludge_dw,
                           N_manure,
                           P_manure,
                           fertilizer_substi_digest_N,
                           fertilizer_substi_digest_P,
                           fertilizer_substi_manure_field_N,
                           fertilizer_substi_manure_field_P)
                                        
     # Update the slude management activity
         
         # Upgrading biogas
         Tech_num_INC_modif_with_excretion[index_biogas_updgrade,index_sludge] = -amount_of_upgrading_act_sludge
         Tech_num_AH_modif_with_excretion[index_biogas_updgrade,index_sludge] = -amount_of_upgrading_act_sludge
         
         # Fertilizer subst
         Tech_num_INC_modif_with_excretion[index_N_substitution,index_sludge] = fertilizer_subst_N_per_kilo_dw_sludge
         Tech_num_AH_modif_with_excretion[index_N_substitution,index_sludge] = fertilizer_subst_N_per_kilo_dw_sludge

         Tech_num_INC_modif_with_excretion[index_P_substitution,index_sludge] = fertilizer_subst_P_per_kilo_dw_sludge
         Tech_num_AH_modif_with_excretion[index_P_substitution,index_sludge] = fertilizer_subst_P_per_kilo_dw_sludge

         # Heat subst
         Tech_num_INC_modif_with_excretion[index_heat_substitution,index_sludge] = MJ_substituted_per_kilo_dw_sludge
         Tech_num_AH_modif_with_excretion[index_heat_substitution,index_sludge] = MJ_substituted_per_kilo_dw_sludge

                                                        


                                                                 
     # Update the rest of the technosphere matrix
         
         
         if growth_stage_index in index_growth_stages_no_filter_laguna:
             
                 Tech_num_AH_modif_with_excretion[index_Nemissions,growth_stage_index] = -N_total_without_removal_AH
                 Tech_num_AH_modif_with_excretion[index_Pemissions,growth_stage_index] = -P_total_without_removal_AH
                 Tech_num_AH_modif_with_excretion[index_sludge,growth_stage_index] = 0
                 Tech_num_INC_modif_with_excretion[index_Nemissions,growth_stage_index] = -N_total_without_removal_INC
                 Tech_num_INC_modif_with_excretion[index_Pemissions,growth_stage_index] = -P_total_without_removal_INC
                 Tech_num_INC_modif_with_excretion[index_sludge,growth_stage_index] = 0
 

     
         else :
                 Tech_num_AH_modif_with_excretion[index_Nemissions,growth_stage_index] = -N_solid_and_liquid_after_removal_AH
                 Tech_num_AH_modif_with_excretion[index_Pemissions,growth_stage_index] = -P_solid_and_liquid_after_removal_AH
                 Tech_num_AH_modif_with_excretion[index_sludge,growth_stage_index] = -total_sludge_AH
                 Tech_num_INC_modif_with_excretion[index_Nemissions,growth_stage_index] = -N_solid_and_liquid_after_removal_INC
                 Tech_num_INC_modif_with_excretion[index_Pemissions,growth_stage_index] = -P_solid_and_liquid_after_removal_INC
                 Tech_num_INC_modif_with_excretion[index_sludge,growth_stage_index] = -total_sludge_INC

                 
     return Tech_num_AH_modif_with_excretion,Tech_num_INC_modif_with_excretion









    
    