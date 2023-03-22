# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 14:25:06 2022

@author: Pierre Jouannais, Department of Planning, DCEA, Aalborg University
pijo@plan.aau.dk

"""
'''
Script which contains the functions which regionalize the uncertain space, 
calculate the probabilities of success in each region, prepare the inputs for PRIM and apply PRIM.

It also contains the code which calls the functions to produce the final results. 
It needs importing the raw output (data points) and the borgonovo's deltas from their respective files.
Replace the path to these files to your specific locations if needed (line 684 and 750)

Executing the whole script will produce the results and save them in the output folders.

Summary:
    -Accessory functions
    -1) Main processing functions for regionalization
    -2) PRIM functions (l. 418)
    -3) RUN the functions and perform regionalization + PRIM (l.685)
     - 3.1) Regionalize and prepare all the inputs for PRIM  (l.689)
     - 3.2)  Apply PRIM on the regionalized and previously prepared dataframes (l1093)
'''

import numpy as np


import multiprocessing
from functools import reduce
import ray


import itertools

import pandas as pd
import decimal
from random import *
import pstats
from itertools import *
from math import*
import csv
import copy
import seaborn as sns
import random

import time

import math 
from ema_workbench.analysis import prim

from ema_workbench import ema_logging


import Borgonovo_indices_to_du_function_modular as borgtodu

ema_logging.log_to_stderr(ema_logging.INFO)

import matplotlib.pyplot as plt

import datetime

import sys




"""Accessory functions """



def find_index(value,reso,min_seq):
    
    """for a value of a parameter/indicator, finds the corresponding lower boundary 
    for a region  given the resolution for this parameter/indicator"""
    
    return floor((value-min_seq)/reso) * reso + min_seq




def find_corresponding_list_multi_index(row,names_indicators,resolutions_indicators,list_min_indicators):
    
    """finds the corresponding region in the multiindex for a row in the raw result"""
    

    values_indicators = [row[a] for a in names_indicators]

    
    list_multi_index = [trunc(find_index(value,reso,min_seq)*1E8)/1E8 for value,reso,min_seq in zip(values_indicators,resolutions_indicators,list_min_indicators)]

    return tuple(list_multi_index)





def initiate_table_multi_index(lentot,multi_index):
    
    """ Creates a multi index table with cells containing empty lists"""
    
    table_multi_index_counting = pd.DataFrame(np.zeros((1,lentot)),index=[0], columns=multi_index,dtype=object)

    for column in list(table_multi_index_counting.columns):
       #print(column)
       table_multi_index_counting.at[0, column] = []
       
    return table_multi_index_counting




def from_multi_to_mono_df_mono_meth(couple_main_support):
    """Function which converts multi-index dataframes to classic dataframes 
    with index as columns (Preparation to PRIM algorithm)"""
 
    table_mono_success=couple_main_support[0].T.reset_index()
    table_mono_support=couple_main_support[1].T.reset_index()
    
    table_mono_success.columns = list(table_mono_success.columns[:-1])+["Values"]
    table_mono_support.columns = list(table_mono_support.columns[:-1])+["Values"]


    return [table_mono_success,table_mono_support]






"""1) Main processing functions for regionalization """

@ray.remote
def process_table_chunk_mono_meth(constant_inputs,table_slice):
    
     """Function which works on 1 chunck of the total raw output. 
        Creates the multi-index dataframe containing the regions and fills the cells 
        with LCA results corresponding to each region.
        Must be called in parallel with ray"""
        
     [names_indicators,
     resolutions_indicators,
     list_max_indicators,
     list_min_indicators,
     index_AH_INC]=constant_inputs  
     
     
     list_seq_indic = []
     
     # Reset index of the fraction of the raw output
     table_slice = table_slice.reset_index(drop="False")
     
     
     # Creating the list of indexes corresponding to the regions
     
     for indic_index in range(len(names_indicators)):
         
         #print(names_indicators[indic_index])
         #print(names_indicators[indic_index])
         max_indicator = list_max_indicators[indic_index]
         
        # print(max_indicator)
         min_indicator = list_min_indicators[indic_index]
         #print(min_indicator)
         reso = resolutions_indicators[indic_index]
         #print(reso)

         seq_indic_pre =np.arange (min_indicator,max_indicator-0.000000001,reso)

         seq_indic =[trunc(value*1E8)/1E8 for value in seq_indic_pre ]
         
  
         list_seq_indic.append(seq_indic)

     lentot = np.prod([len(seq) for seq in list_seq_indic])


     # Create empty multi index containing all the regions 
     multi_index = pd.MultiIndex.from_product(list_seq_indic, names=names_indicators)
     # Put list of zeros in each of the cells
     table_multi_index_counting = initiate_table_multi_index(lentot,multi_index)


     # Assign LCA results of the raw output to the corresponding region
     for row_index in range(len(table_slice)):
         
        row= table_slice.loc[row_index,:]
        
        # Find the corresponding position for this row in the multiindex table
        
        list_index = find_corresponding_list_multi_index(row,names_indicators,resolutions_indicators,list_min_indicators)
        
        #print("OOOK3")
        #print("OOK4")



        table_multi_index_counting[list_index][0].append([row[index_AH_INC[0]],row[index_AH_INC[1]]])  # Assign at the righ position the duet [concept, alternative]

   
    
     return table_multi_index_counting
 
    
 

def processing_data_mono_meth(tab_res,
                              names_indicators,
                              resolutions_indicators,
                              number_cores,
                              index_AH_INC):
   
    """Function which takes the output table (tab_res), divides into chunk for 
    parallel processing and combines the results into the total multi-index dataframes
    with all regions containing the results of the associated data points"""
            
   # Collect min and max value for all parameters/indicators
    list_max_indicators=[]
    list_min_indicators=[]
    
    for indic_index in range(len(names_indicators)):

        max_indicator = max(tab_res[names_indicators[indic_index]])
        min_indicator = min(tab_res[names_indicators[indic_index]])
        
        list_max_indicators.append(max_indicator)
        list_min_indicators.append(min_indicator)
    

    tab_res_clean = tab_res.copy()
    # Trim max values for the indicators/parameters out of the raw output
    for indic in range(len(names_indicators)):
        tab_res_clean = tab_res_clean[tab_res_clean[names_indicators[indic]]!=list_max_indicators[indic]]

             
    
    ray.shutdown()
    
    ray.init()    
    
    constant_inputs = ray.put([names_indicators,
                               resolutions_indicators,
                               list_max_indicators,
                               list_min_indicators,
                               index_AH_INC])
    
    
    
    # Start Ray.

    # Cut raw output in number of cores
    
    tab_split = np.array_split(tab_res_clean, number_cores) 

    print("Start building multiindex dataframes")
    
    # Collect meta index tables
    collected_meta_index_table = ray.get([process_table_chunk_mono_meth.remote(constant_inputs,table_slice) for table_slice in tab_split])
    

    # Sum all chuinks to obtain final multi-index table
    
    
    combined_multi_index_table=reduce(lambda x, y: x.add(y, fill_value=0), collected_meta_index_table)
       # print("ccc")
    
        
    return collected_meta_index_table,combined_multi_index_table



def apply_resolution(name_indicator,tab_res,resolution_percent_range):
    """Calculate the absolute resolution dU for a parameter/indicator based on its relative resolution """
    
    min_ = min(tab_res[name_indicator])
    max_ = max(tab_res[name_indicator])
    absolute_resolution=(max_- min_) * resolution_percent_range
    
    return absolute_resolution




def process_line_multi_index(cell):
    """Function which, for a cell in the multi-index table (region) which already 
    contains the data points, calculates the proportion of successes and other values"""

    
    
    count_points = 0
    count_sucess = 0
    list_ratios = [] 

    #AH = concept
    # INC = alterative
    if len(cell)!=0:  
        min_AH= cell[0][0] # Initializee min #first value of the first cell
        max_AH= cell[0][0] # Initializee max
        min_INC= cell[0][1] # Initializee min  #second value of the first cell
        max_INC= cell[0][1] # Initializee max
        
        for duet_AH_INC in cell: # For all duets [AH, iNC]
            count_points+=1
            ratio= duet_AH_INC[0]/duet_AH_INC[1]
            list_ratios.append(ratio)
            
            if duet_AH_INC[0]<duet_AH_INC[1]:
                count_sucess+=1
            
            if duet_AH_INC[0] < min_AH:
                
                min_AH =duet_AH_INC[0]
                
            elif duet_AH_INC[0] > max_AH:
            
                max_AH =duet_AH_INC[0]
            
            if duet_AH_INC[1] < min_INC:
                
                min_INC =duet_AH_INC[1]
                
            elif duet_AH_INC[1] > max_INC:
            
                max_INC =duet_AH_INC[1]     
        
        #print("UUU")
        
        prop_sucess = count_sucess/len(cell)

        q25 = np.quantile(list_ratios, 0.25)
        q50 = np.quantile(list_ratios, 0.50)
        q75 = np.quantile(list_ratios, 0.75)
        #print("TTT")
    
        #qthresh= np.quantile(list_ratios, thresh)
        mean= np.mean(list_ratios)
        #print("TVV")
    else :  # there is no point in this cell
        q25 = float("nan")
        q50 = float("nan")
        q75 = float("nan")
        #print("TTT")
    
        #qthresh= float("nan")
        mean= float("nan")
        min_AH= float("nan") 
        max_AH= float("nan") 
        min_INC= float("nan") 
        max_INC= float("nan") 
        prop_sucess= float("nan")
        
    return [q25,q50,q75,mean,count_points,max_INC,min_INC,max_AH,min_AH,prop_sucess]
    
    



def process_success_rates_mono_parallel(combined_multi_index_table):
    """Function which calculates the proportion of successes (and other values) 
    in every cell of the multi-index dataframe
    It calls the next function below in parallel over the different chunks of the whole multi-index dataframe"""

        
    # Divide multi-index dataframes in number of cores
    tab_split = np.array_split(combined_multi_index_table, number_cores,axis=1) 
   
    
    # Start Ray.

    ray.shutdown()
    
    ray.init()    
    

    print("Start building multiindex dataframes")
    
    # Collect processed meta index tables
    collected_tables = ray.get([process_one_chunk.remote(table_slice) for table_slice in tab_split])
    
    # Multi-index dataframe with proportions of successes in each region
    main_result_table = pd.concat([a[0] for a in collected_tables],axis=1)
    
    # Multi-index  with number of data points of successes in each region
    support_result_table = pd.concat([a[1] for a in collected_tables],axis=1)
    
    return [main_result_table,support_result_table]  

@ray.remote
def process_one_chunk(table_slice):
    
    """Function which calculates the proportion of successes (and other values) 
    in every cell of 1 chunk the multi-index dataframe.
    Is called in parallel with ray"""
    
    main_result =copy.deepcopy(table_slice)
    support_result =copy.deepcopy(table_slice)
    
    for index_column in table_slice:
        
        rescell=process_line_multi_index(table_slice[index_column][0])
        
        [q25,q50,q75,mean,count_points,max_INC,min_INC,max_AH,min_AH,prop_sucess] = rescell

    

        main_result[index_column][0] = prop_sucess

        support_result[index_column][0] = count_points
    
    return [main_result,support_result]  






    
"""2) PRIM functions """


def prep_for_prim_1_categ(success_meth, thresh_dec, mode):  # mode = inf or sup to the probability threshold

    """Function which prepares the table for the PRIM algorithm.
    Returns table_parameters (X) and y with True and False depending on the 
    comparison to the decision threshold"""
    
    success_meth_clean = success_meth.copy()
    
        
    # Prepare output table (Values)
    table_output_as_list = list(success_meth_clean["Values"])
    table_output= pd.DataFrame({"Values":table_output_as_list}).reset_index(drop="False")
    
    # Prepare Parameters table
    table_parameters = success_meth_clean[success_meth_clean.columns[:-1]].reset_index(drop="False")
    
    recombine_table = pd.concat([table_parameters,table_output], axis=1)

    if mode == "sup":
        y = table_output["Values"] > thresh_dec  # the cases of interest
        recombine_table["Success"] = recombine_table["Values"]>thresh_dec

    elif mode == "inf": 
        y = table_output["Values"] < thresh_dec  # the cases of interest
        recombine_table["Success"] = recombine_table["Values"]<thresh_dec

    else:
        sys.exit("Invalid mode of success")
    # For manual pairplot
    

    return table_parameters,y, recombine_table







def get_box_data(box_):
    
    """Function which collects the data about a PRIM box"""    
    
    box_data = box_.inspect(style="data")
    box_data_series = box_data[0][0]
    dens = box_data_series["density"] 
    cover = box_data_series["coverage"] 
    mass = box_data_series["mass"]
    mean = box_data_series["mean"]
    
    boundaries_box = box_data[0][1]
    boundaries_box_df = boundaries_box.T.reset_index()

    
    return dens,cover,mass,mean,boundaries_box_df






def go_over_list_PRIM_while_2(tables_for_PRIM,thresh_dec,tresh_density,list_meth_codes, min_mass_param, min_mass_indicator):
    
    """Function which calls PRIM for all impact categories on the regionalized outputs,
    saves and plot discovered boxes"""    

    
    mode="sup"  # Looking for regions with proba > treshold. This could be defined outside and used as a parameter instead


    
    list_result_each_meth = []
    for meth_index in range(len(tables_for_PRIM)): # For all impact categories/methods
        
        if "parameters" in list_meth_codes[meth_index]: # If we are working at high dimensionality / with parameters
            min_mass = min_mass_param  # The minimum mass for PRIM
        else:   # If we are working at low dimensionality / with parameters
            min_mass = min_mass_indicator

        
        
        result_for_1_meth = [] # 
        
        print("meth_type_factor",list_meth_codes[meth_index])
        meth_code = list_meth_codes[meth_index]   # Collect the method code
        table_meth = tables_for_PRIM[meth_index]  # Collect the regionalized output ready for PRIM
        
        print("len(table_meth)",len(table_meth))
        success_table = table_meth[0]
        support_table = table_meth[1]
        
        # Filter and clean spaces
        success_table_clean = success_table.dropna()
        support_table_clean= support_table[support_table['Values']!=0]

        number_dropped = success_table.shape[0]- success_table_clean.shape[0] # how many empty regions
        percentdrop =  number_dropped/ success_table.shape[0]
        print("percentdrop",percentdrop)

        len_list_box = 0 #Initialize with True
        
        thresh_dec_update = thresh_dec
        
        no_box = False
        
        x = datetime.datetime.now()

        month=str(x.month)
        day=str(x.day)
        microsec=str(x.strftime("%f"))
        
        code_boxes = meth_code+"_"+month+"_"+day+"_"+microsec
        
        while len_list_box<=1 and thresh_dec_update > 0.1: # As long as we have not found a single box

            list_box = []

            # Prepare the input for PRIM
            table_parameters,y, recombine_table = prep_for_prim_1_categ(success_table_clean, thresh_dec_update, mode) 
            

            # Call PRIM
            prim_alg = prim.Prim(table_parameters, y, threshold=tresh_density,peel_alpha=0.1, mass_min = min_mass) 
            print("OK")

            # Collect the found boxes
            for i in range(20):
                print(i)
                try:
                    list_box.append(prim_alg.find_box())
                except:
                    print("Error with qp")   # Possible error that must be ignored
            
            list_box = [i for i in list_box if i is not None]
            
            len_list_box = len(list_box)
            
            # Decrease the decsion threshold so that if no box has been found, we try again at a lower treshold
            thresh_dec_update = thresh_dec_update - 0.1
            
            if thresh_dec_update<=0.1:
                no_box = True


        
        if not no_box: # if there are some boxes then study the boxes
            print("Found boxes")
            inital_density = list_box[0].inspect(0,style="data")[0][0]["density"]
            
            result_for_1_meth.append([len(list_box),percentdrop,inital_density])

            cover_tot = 0 #Initialize the coverage for all boxes
            for index_box in range(len(list_box)-1):  # -1 because the last one is empty
                
                box = list_box[index_box]
                print("index_box",index_box )
                [dens,cover,mass,mean,boundaries_box_df] = get_box_data(box)
                cover_tot +=cover
                name_file = "boundary_box_"+str(index_box) + code_boxes
                
                boundaries_box_df.to_csv("../PRIM_process/Boxes/"+str(name_file)+'.csv', sep=';', encoding='utf-8')

                plt.figure()
                plot = box.show_pairs_scatter()  # Plot box
                
                #Quick fix name impact and parameters for plotting
                if meth_code == "GWP100":
                    meth_code = "GW100"
                elif meth_code ==  "TETPinf":
                    meth_code = "TETinf"
                elif meth_code ==  "FEP":
                    meth_code = "FE"    

                dens_trunc = '%.3f'%(dens)    
                inital_density_trunc  = '%.3f'%(inital_density) 
                cover_trunc = '%.3f'%(cover) 
                mass_trunc = '%.3f'%(mass) 

                thresh_dec_trunc ='%.3f'%(thresh_dec_update+0.1) 
                
                name_box = "box_"+str(index_box)+"_"+code_boxes
                
                
                plot.fig.suptitle(str(meth_code)+" "+ "dec_tresh = "+thresh_dec_trunc+" "+"init dens="+inital_density_trunc +"  box dens="+dens_trunc +"  mass="+mass_trunc  +"  cover="+cover_trunc)
                plot.tight_layout()
                
                plot.savefig("../PRIM_process/Box pairplots/"+name_box+".png",dpi=500,bbox_inches="tight")
                
                res_for_1_box =[dens,cover,mass,mean,boundaries_box_df]
                

        
                result_for_1_meth.append(res_for_1_box)
                
                
                # Show peeling trajectory
                
                
                plt.figure()
                plot_peeling = box.show_tradeoff()
                plot_peeling.suptitle(str(meth_code)+" "+ "dec_tresh = "+str(thresh_dec_update+0.1)+" "+"init dens="+inital_density_trunc +"  box dens="+dens_trunc +"  mass="+mass_trunc  +"  cover="+cover_trunc)
                plot_peeling.tight_layout()
                
                plot_peeling.savefig("../PRIM_process/Peeling/"+name_box+".png",dpi=500,bbox_inches="tight")
                
        
        else:
            
            print("No box was found")

            result_for_1_meth ="no_box"
            
            table_parameters_manual,y_manual, recombine_table_manual = prep_for_prim_1_categ(success_table_clean, thresh_dec, mode)
        
            
            #Manual pairplots for the decision threshold
            init_density = sum(recombine_table_manual["Success"])/len(recombine_table_manual["Success"])
            
            plt.figure()

            grid = sns.pairplot(
            data=recombine_table_manual,
           
            hue="Success",
            vars=recombine_table_manual.columns[:-2],
            diag_kind="kde",
            diag_kws={ "common_norm": False, "fill": False},
            plot_kws={'alpha':0.01},
            palette =["red","green"])
            
            
            name_manual = "Manual_"+"_"+code_boxes + "_" +str(thresh_dec_update)
        
            grid.fig.suptitle("Manual_"+str(meth_code)+" "+ "dec_tresh = "+str(thresh_dec)+" "+"init dens="+str(init_density))
            grid.tight_layout()
                
            grid.savefig("../PRIM_process/Manual pairplots/"+name_manual+".png")
                
            
            
            
            
            
            
        list_result_each_meth.append(result_for_1_meth)
   
    

        if not no_box:
            print("OKK")
            table_info=pd.DataFrame(np.array([cover_tot,inital_density,len_list_box])).T
            table_info.columns = ["Total coverage",  "Initial density", " Number boxes"]
            name_info_csv =  name_box = "info_total"+code_boxes
            table_info.to_csv("../PRIM_process/info/"+str(name_info_csv)+'.csv', sep=';', encoding='utf-8')

    return list_result_each_meth








"""3) RUN the functions and perform regionalization + PRIM"""



"""3.1) Regionalize and prepare all the inputs for PRIM"""

# Raw output
tab_res_0 = pd.read_csv("../Outputs/Outputs_recombined/total_output_for_PRIM_1_10_609984.csv",
                                 sep=";",
                                 header=0,
                                 encoding='unicode_escape',
                                 engine='python')


 

# Define list Meth


list_meth_code =  ["GWP100","eutrophication",'FEP','TETPinf','FETPinf']



list_types_factor = ["parameters", "indicators_no_rates"] # High and low dimensionality

    
number_cores = multiprocessing.cpu_count()

time1=time.time()

# Initialize the list which will contain all the tables ready for PRIM (for all impact categories)
tables_for_PRIM = [] 
for meth in list_meth_code:
    
    impact_category = meth
    
    for type_factor in list_types_factor:  
        
        print("meth",meth)

        print("type_factor",type_factor)
        
        if "indicators" in type_factor:
        
            dumax = 1
            dumin = 1/8
        
            
        elif type_factor == "parameters":
        
            dumax = 1
            dumin = 1/4
        
        
        
        
        trim=0.02
        
        
        
        tab_res = tab_res_0.copy()
        
        # Get the columns with the LCA results in the raw output
        list_index_coupled_AH_INC =[[tab_res.columns.get_loc(meth+"AH"),tab_res.columns.get_loc(meth+"INC")] for meth in list_meth_code] 
        
        
        # Here modulated resolution
        
        index_AH_INC = list_index_coupled_AH_INC[0] 
        
        # Load Borgonovo's deltas
        if "indicator" in type_factor: 
            name_csv = "../borgonovo_indices/borgonovo_"+impact_category+"_indicators.csv" 
            
        else:
            name_csv = "../borgonovo_indices/borgonovo_"+impact_category+"_parameters.csv"
            
        borgonovo_indices=pd.read_csv(name_csv,
                                         sep=",",
                                         header=None,
                                         encoding='unicode_escape',
                                         engine='python')
        
        
        # Just to get the column names
        combined_results_head=pd.read_csv("../Outputs/Outputs_recombined/1_9_248054_"+impact_category+".csv",
                                         sep=";",
                                         header=0,
                                         encoding='unicode_escape',
                                         engine='python',
                                         nrows=2)
        
        
        
            
        if "indicator" in type_factor:
            # We will keep only the LCA results and the parameters/indicators at low dimensionality
             columns_to_remove= ['Unnamed: 0',          
                                 'Unnamed: 0.1',
                                 impact_category+'AH',
                                 impact_category+'INC',
                                 impact_category+'AH/INC',
            'HAIT_loss_lev',
             'HAIT_loss_red',
             'FRIT_loss_lev',
             'FRIT_loss_red',
             'GOIT1bis_loss_lev',
             'GOIT1bis_loss_red',
             'GODK_loss_lev',
             'GODK_loss_red',
             'SFDK1_loss_lev',
             'SFDK2_loss_lev',
             'SFDK1_loss_red',
             'FRIT_micro_dose',
             'GODK_micro_dose',
             'GOIT1bis_micro_dose',
             'HAIT_micro_dose',
             'SFDK2_micro_dose',
             'SFDK1_micro_dose',
             'HAIT_FCR_red_ratio_frac',
             'FRIT_FCR_red_ratio_frac',
             'GOIT1bis_FCR_red_ratio_frac',
             'GODK_FCR_red_ratio_frac',
             'SFDK1_FCR_red_ratio_frac',
             'SFDK2_FCR_red_ratio_frac',
             'INC-AH_'+impact_category,
             'impact_Drug_prod'+impact_category]
             
        else:     
             columns_to_remove = ['Unnamed: 0',
                             'Unnamed: 0.1',
                             impact_category+'AH',
         impact_category+'INC',
         impact_category+'AH/INC',
         'Loss_fraction_increase_INC_indic',
         'Loss_fraction_increase_AH_indic',
         'Input_micro_comp_per_FU_indic',
         'Economic_FCR_red_0_indic',
         'Economic_FCR_red_1_indic',
         'Economic_FCR_red_rate_1',
         'ECO_FCR_AH_1',
         'Biological_FCR_red_1_indic',
         'Biological_FCR_red_0_indic',
         'Biological_FCR_red_rate_1',
         'dose_per_total_bio_AH_1',
         'Economic_FCR_red_rate_1_tot_bio',
         'Biological_FCR_red_rate_1_tot_bio',
         'INC-AH_'+impact_category,
         'impact_Drug_prod'+impact_category]
        
        
        
            
         
         
        if type_factor == "parameters":
            
            uncertain_factors = ['HAIT_loss_lev',
             'HAIT_loss_red',
             'FRIT_loss_lev',
             'FRIT_loss_red',
             'GOIT1bis_loss_lev',
             'GOIT1bis_loss_red',
             'GODK_loss_lev',
             'GODK_loss_red',
             'SFDK1_loss_lev',
             'SFDK2_loss_lev',
             'SFDK1_loss_red',
             'FRIT_micro_dose',
             'GODK_micro_dose',
             'GOIT1bis_micro_dose',
             'HAIT_micro_dose',
             'SFDK2_micro_dose',
             'SFDK1_micro_dose',
             'HAIT_FCR_red_ratio_frac',
             'FRIT_FCR_red_ratio_frac',
             'GOIT1bis_FCR_red_ratio_frac',
             'GODK_FCR_red_ratio_frac',
             'SFDK1_FCR_red_ratio_frac',
             'SFDK2_FCR_red_ratio_frac',
             'impact_drug_'+impact_category,
             "bioact_molec_dbio"]
            
           
        # Indicators With rate. Other types of indicators/dimensionalities 
        #not used in the article but that can be used to regionalize
        
        elif type_factor == "indicators_rates":
            uncertain_factors = [
                 'Biological_FCR_red_rate_1',
                 'Loss_fraction_increase_INC_indic',
                 'Economic_FCR_red_rate_1',
                 'impact_drug_'+impact_category,
                 'bioact_molec_dbio']
            
                # Result columns
                
            # Trim results if necessary
            tab_res = tab_res_0[
                            (tab_res_0['Biological_FCR_red_rate_1'] < np.quantile(tab_res_0['Biological_FCR_red_rate_1'], 1-trim))  &
                            (tab_res_0['Biological_FCR_red_rate_1'] > np.quantile(tab_res_0['Biological_FCR_red_rate_1'], trim))  &
                            (tab_res_0['Loss_fraction_increase_INC_indic'] < np.quantile(tab_res_0['Loss_fraction_increase_INC_indic'], 1-trim))  &
                            (tab_res_0['Loss_fraction_increase_INC_indic'] > np.quantile(tab_res_0['Loss_fraction_increase_INC_indic'], trim)) &
                            (tab_res_0['Economic_FCR_red_rate_1'] < np.quantile(tab_res_0['Economic_FCR_red_rate_1'], 1-trim)) &
                            (tab_res_0['Economic_FCR_red_rate_1'] > np.quantile(tab_res_0['Economic_FCR_red_rate_1'], trim)) &
                            (tab_res_0['impact_drug_'+impact_category] < np.quantile(tab_res_0['impact_drug_'+impact_category], 1-trim)) &
                            (tab_res_0['impact_drug_'+impact_category] > np.quantile(tab_res_0['impact_drug_'+impact_category], trim)) &
                            (tab_res_0['bioact_molec_dbio'] < np.quantile(tab_res_0['bioact_molec_dbio'], 1-trim)) &
                            (tab_res_0['bioact_molec_dbio'] > np.quantile(tab_res_0['bioact_molec_dbio'], trim))]
        
           
        
        
        
        elif type_factor =="indicators_rates_bio":
        
            uncertain_factors = [
             'Economic_FCR_red_rate_1_tot_bio',
             'Biological_FCR_red_rate_1_tot_bio',
             'Loss_fraction_increase_INC_indic',
             'impact_drug_'+impact_category,
             'bioact_molec_dbio'
             ]
            
                    
        
            tab_res = tab_res_0[
                            (tab_res_0['Economic_FCR_red_rate_1_tot_bio'] < np.quantile(tab_res_0['Economic_FCR_red_rate_1_tot_bio'], 1-trim))  &
                            (tab_res_0['Economic_FCR_red_rate_1_tot_bio'] > np.quantile(tab_res_0['Economic_FCR_red_rate_1_tot_bio'], trim))  &
                            (tab_res_0['Loss_fraction_increase_INC_indic'] < np.quantile(tab_res_0['Loss_fraction_increase_INC_indic'], 1-trim))  &
                            (tab_res_0['Loss_fraction_increase_INC_indic'] > np.quantile(tab_res_0['Loss_fraction_increase_INC_indic'], trim)) &
                            (tab_res_0['Biological_FCR_red_rate_1_tot_bio'] < np.quantile(tab_res_0['Biological_FCR_red_rate_1_tot_bio'], 1-trim)) &
                            (tab_res_0['Biological_FCR_red_rate_1_tot_bio'] > np.quantile(tab_res_0['Biological_FCR_red_rate_1_tot_bio'], trim)) &
                            (tab_res_0['impact_drug_'+impact_category] < np.quantile(tab_res_0['impact_drug_'+impact_category], 1-trim)) &
                            (tab_res_0['impact_drug_'+impact_category] > np.quantile(tab_res_0['impact_drug_'+impact_category], trim)) &
                            (tab_res_0['bioact_molec_dbio'] < np.quantile(tab_res_0['bioact_molec_dbio'], 1-trim)) &
                            (tab_res_0['bioact_molec_dbio'] > np.quantile(tab_res_0['bioact_molec_dbio'], trim))]
        
           
            
            
        ## Indicators no rates
        
        elif type_factor =="indicators_no_rates":
        
            uncertain_factors = [
              'Economic_FCR_red_1_indic',
              'Biological_FCR_red_1_indic',
             'Loss_fraction_increase_INC_indic',
               'Input_micro_comp_per_FU_indic',
             'impact_drug_'+impact_category,
             'bioact_molec_dbio'
             ]
        
        
            tab_res = tab_res_0[
                            (tab_res_0['Economic_FCR_red_1_indic'] < np.quantile(tab_res_0['Economic_FCR_red_1_indic'], 1-trim))  &
                            (tab_res_0['Economic_FCR_red_1_indic'] > np.quantile(tab_res_0['Economic_FCR_red_1_indic'], trim))  &
                            (tab_res_0['Loss_fraction_increase_INC_indic'] < np.quantile(tab_res_0['Loss_fraction_increase_INC_indic'], 1-trim))  &
                            (tab_res_0['Loss_fraction_increase_INC_indic'] > np.quantile(tab_res_0['Loss_fraction_increase_INC_indic'], trim)) &
                            (tab_res_0['Biological_FCR_red_1_indic'] < np.quantile(tab_res_0['Biological_FCR_red_1_indic'], 1-trim)) &
                            (tab_res_0['Biological_FCR_red_1_indic'] > np.quantile(tab_res_0['Biological_FCR_red_1_indic'], trim)) &
                            (tab_res_0['impact_drug_'+impact_category] < np.quantile(tab_res_0['impact_drug_'+impact_category], 1-trim)) &
                            (tab_res_0['impact_drug_'+impact_category] > np.quantile(tab_res_0['impact_drug_'+impact_category], trim)) &
                            (tab_res_0['bioact_molec_dbio'] < np.quantile(tab_res_0['bioact_molec_dbio'], 1-trim)) &
                            (tab_res_0['bioact_molec_dbio'] > np.quantile(tab_res_0['bioact_molec_dbio'], trim)) &
                            (tab_res_0['Input_micro_comp_per_FU_indic'] < np.quantile(tab_res_0['Input_micro_comp_per_FU_indic'], 1-trim)) &
                            (tab_res_0['Input_micro_comp_per_FU_indic'] > np.quantile(tab_res_0['Input_micro_comp_per_FU_indic'], trim))]
        
           
        
        
        elif type_factor =="indicators_no_rates_bio":
        
            uncertain_factors = [
              'Economic_FCR_red_1_indic',
              'Biological_FCR_red_1_indic',
             'Loss_fraction_increase_INC_indic',
               'dose_per_total_bio_AH_1',
             'impact_drug_'+impact_category,
             'bioact_molec_dbio'
             ]
            
            tab_res = tab_res_0[
                            (tab_res_0['Economic_FCR_red_1_indic'] < np.quantile(tab_res_0['Economic_FCR_red_1_indic'], 1-trim))  &
                            (tab_res_0['Economic_FCR_red_1_indic'] > np.quantile(tab_res_0['Economic_FCR_red_1_indic'], trim))  &
                            (tab_res_0['Loss_fraction_increase_INC_indic'] < np.quantile(tab_res_0['Loss_fraction_increase_INC_indic'], 1-trim))  &
                            (tab_res_0['Loss_fraction_increase_INC_indic'] > np.quantile(tab_res_0['Loss_fraction_increase_INC_indic'], trim)) &
                            (tab_res_0['Biological_FCR_red_1_indic'] < np.quantile(tab_res_0['Biological_FCR_red_1_indic'], 1-trim)) &
                            (tab_res_0['Biological_FCR_red_1_indic'] > np.quantile(tab_res_0['Biological_FCR_red_1_indic'], trim)) &
                            (tab_res_0['impact_drug_'+impact_category] < np.quantile(tab_res_0['impact_drug_'+impact_category], 1-trim)) &
                            (tab_res_0['impact_drug_'+impact_category] > np.quantile(tab_res_0['impact_drug_'+impact_category], trim)) &
                            (tab_res_0['bioact_molec_dbio'] < np.quantile(tab_res_0['bioact_molec_dbio'], 1-trim)) &
                            (tab_res_0['bioact_molec_dbio'] > np.quantile(tab_res_0['bioact_molec_dbio'], trim)) &
                            (tab_res_0['dose_per_total_bio_AH_1'] < np.quantile(tab_res_0['dose_per_total_bio_AH_1'], 1-trim)) &
                            (tab_res_0['dose_per_total_bio_AH_1'] > np.quantile(tab_res_0['dose_per_total_bio_AH_1'], trim))]
        
           
        
        # Get the relative resolutions (and save the associated plots in the folder)
        table_du=borgtodu.get_du_from_borgonovo_indices_modular(borgonovo_indices,
                                          combined_results_head,
                                          impact_category,
                                          type_factor,
                                          dumax,
                                          dumin,
                                          uncertain_factors,
                                          columns_to_remove)
        
        
        
        resolutions_indicators_fraction = list(table_du["dU"])
        print("table_du",table_du)
        
        print("resolutions_indicators_fraction", resolutions_indicators_fraction)
        
        names_indicators = list(table_du["Parameter"])
        
        
        
        # Absolute resolutions for parameters/indicators
        resolutions_indicators = [apply_resolution(name,tab_res,reso_frac) for name,reso_frac in zip(names_indicators,resolutions_indicators_fraction)]
        
        
        
        
        
        
        
        
        # Generate the multi-index tables with the data points in each corresponding region
        collected_meta_index_table,combined_multi_index_table =processing_data_mono_meth(tab_res,
                                     uncertain_factors,
                                     resolutions_indicators,
                                     number_cores,
                                     index_AH_INC)
        
            
           
            
            
        
        
        # Generate multi-index dtaatframes with proportion of successes and 
        #the "mirror"support multindex_dataframe with the number of simulations in each region
        
        couple_main_support=process_success_rates_mono_parallel(combined_multi_index_table)
        
        # Flatten from multi-index dataframe to classic dataframes for PRIM
        table_for_PRIM_1_meth =from_multi_to_mono_df_mono_meth(couple_main_support)
        
        tables_for_PRIM.append(table_for_PRIM_1_meth)
        
        # PRIM input ready
        success_meth_0 = table_for_PRIM_1_meth[0] ### success meth 0
        # Support table ready
        support_meth_0 = table_for_PRIM_1_meth[1] ### success meth 0
        
        x = datetime.datetime.now()
        
        month=str(x.month)
        day=str(x.day)
        microsec=str(x.strftime("%f"))
                     
        
        # Save the dataframe ready for PRIM if needed later
        name_file ="../PRIM_process/intermediate tables for PRIM/success_meth" +impact_category+"_"+type_factor
        
        success_meth_0.to_csv(name_file+'.csv', sep=';', encoding='utf-8')
        
        name_file = "../PRIM_process/intermediate tables for PRIM/support_meth" +impact_category+"_"+type_factor
        
        support_meth_0.to_csv(name_file+'.csv', sep=';', encoding='utf-8')
        
        
        # Check robustness of each region (number of points in each region)
        # Max = maximum coverage
        # Mean covergae when not 0
        # Mean coverage
        sumtot_=0
        max_=0
        count_non_0=0
        for i in list(support_meth_0["Values"]):
            #print(i)
            
            if i!=0:
                count_non_0 +=1
                sumtot_ +=i
                
                if i> max_:
                   max_ = i
        
        # Mean covergae when not 0
        mean_robu_when_not_0 = sumtot_/count_non_0            
        mean_robu = sumtot_/len(list(support_meth_0["Values"]))            
    
        empty_cells = sum(math.isnan(x) for x in list(success_meth_0["Values"]))
        number_regions = len(list(success_meth_0["Values"]))
        # Info robustness
        table_info_robu=pd.DataFrame(np.array([mean_robu_when_not_0,mean_robu,empty_cells,number_regions])).T
        table_info_robu.columns = ["Mean number of points per region if not empty",  "Mean number of points per region", " Number of empty regions" , "Number of regions"]
        name_info_robu_csv =  "info_robu"+impact_category +" " + type_factor +" " + month +"_" + day + "_" + microsec
        table_info_robu.to_csv("../PRIM_process/info/"+str(name_info_robu_csv)+'.csv', sep=';', encoding='utf-8')



"""3.2) Apply PRIM on the regionalized and previously prepared dataframes"""


print(" NOW PRIM ")     




# The function needs to know the method and the dimensionality level to work at.
# This 2 informations are aggregated in a "code" which combines the method and the "type_factor"= which dimensionality level
list_type_factors_adjusted = [[meth + "_" + type_factor for type_factor in list_types_factor] for meth in list_meth_code]
# Flatten to a single list
list_meth_codes_type_factor = list(itertools.chain(*list_type_factors_adjusted))

# Minimum mass for PRIM
min_mass_param = 0.01
min_mass_indicator = 0.01


# Set the axes labels font size
plt.rc('axes', labelsize=10.2)
# Set the font size for x tick labels
plt.rc('xtick', labelsize=10.2)
# Set the font size for y tick labels
plt.rc('ytick', labelsize=10.2)   
        
# Apply PRIM for all categories and save results in the folders.
Results_for_all_categories = go_over_list_PRIM_while_2(tables_for_PRIM,0.85,0.90,list_meth_codes_type_factor, min_mass_param, min_mass_indicator)


time2=time.time()
print("timetot",time2-time1)

print("DONE")