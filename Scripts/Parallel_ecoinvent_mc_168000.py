# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 09:38:46 2022

@author: Pierre Jouannais, Department of Planning, DCEA, Aalborg University
pijo@plan.aau.dk

"""
'''
Script which allows perfoming a very high number of MonteCarlo Iterations for the impacts of 1 unit of each technosphere inputs to the FU in the model.
It saves the result in a pickle object importable into other scripts
The script dvivides a high number of MC into chunks spread in parallel over different CPUS.
Requires instances with multiple CPUS

This version was executed to simulate 168000 MC for the 95 FU of the foreground and all the impact categories. 
See line 647 for "size_diviz" =1000 and number_chunks=168.
It was run on a remote server with 64 CPUs and took some days.


-1) INITIALIZE (l.149)
-2) Function for fast parallel montecarlo iterations in the background (l.487))
-3) RUN functions and generate the MC iterations in the background (l.646))



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
# os.getcwd()
#os.chdir('/home/ubuntu/Work_folder/Code_third_paper/Scripts')
#os.chdir("/home/ubuntu/work_folder/Code_third_paper/Scripts")

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
import multiprocessing


import pert


import scipy
import numpy as np
from scipy import linalg

from ast import literal_eval



from itertools import compress


import functools
import operator




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



"""1) INITIALIZE"""


# Import Technosphere Fish Production

#"C:\Users\GF20PZ\OneDrive - Aalborg Universitet\Dokumenter\AAU\FIsh Farm model\From July 2022\Technosphere matrices"

Techno_Matrix_Fish_loaded = pd.read_excel('../Technosphere_Fish_micro_17_03.xlsx',header=None)


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





                             
                             




# Initialize brightway and calculate LCIA for1 unit of each tehcnosphere inputs




bw.projects.set_current('Fish_and_micro_2103') 



# Loading Ecoinvent
Ecoinvent = bw.Database('ecoinvent 3.8 conseq')


# Loading foreground database

FISHMIC = bw.Database('AH_combi_1')

MICAH = bw.Database('Micro_for_combi')



biosphere=bw.Database('biosphere3')




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
    



# Electricity mixes for europe


list_electricity_CNTR_FU = []

list_electricity_CNTR_names = []



for act in Ecoinvent:
            if 'market for electricity' in act['name']:
                if 'medium voltage' in act['name']:
                    if 'municipal' not in act['name']:
                        if 'label-certified' not in act['name']:
                            if act['location'] in ['ES','SE','IT','PT','GB','FR','GR','CY','IE','DE']:   # 10 countries identified by Skarka
                                print(act)  
                                list_electricity_CNTR_FU.append({act:1})
                                list_electricity_CNTR_names.append(act["location"])




### Additional variables for the regionalization of the wastewater treatment of microalgae

electricity_low_voltage_input_per_m3_wastewater = 0.20571
electricity_high_voltage_input_per_m3_wastewater = 0.0086946



# Choose methods
for meth in list(bw.methods):
    if "IPCC" in meth:
        print(meth)


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


                      
#activities_fish_background

for name in activities_fish_background:

    print(name)
    FISHMIC.search(name)[0]
    
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

    


"""Create a list of the names of the foreground activities which are directly and only connected to an electricity input"""

# Function to modify the electricty impacts according to the national mix

# Necessary
list_processes_electricity_micro = []
for exc in list(Molprod.exchanges()):
#    print(exc)
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
          
 




# For natural_gas regionalization
   
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





# All Fu directly connected to the backkground only
list_FU_combined_mc = list_micro_algae_FU + list_fish_FU + list_electricity_CNTR_FU + list_natural_gas_CNTR_FU
list_FU_combined_names_mc = list_micro_algae_names + activities_fish_background + list_electricity_CNTR_names + list_natural_gas_CNTR_names

        









"""2) Function for fast parallel montecarlo iterations in the background"""


@ray.remote
def parallel_MC_backgound(constant_inputs):
    """Function for fast parallel montecarlo iterations in the background.
    Calculates a number of Montecalro LCA "size" for 1 unit of each technosphere inputs in the foreground .
    For all impact categories.
    
    
    """

    
    [list_micro_algae_FU,
    list_fish_FU,
    list_electricity_CNTR_FU,
    list_natural_gas_CNTR_FU,
    list_micro_algae_names,
    activities_fish_background,
    list_electricity_CNTR_names,
    list_natural_gas_CNTR_names,
    C_matrixes,
    size] = constant_inputs
    
    #Initialize empty list that will contain the results
    
    # For fish
    
    list_array_mc_sample_fish = [np.array([[0]*len(list_fish_FU)  ]*size,dtype="float32") for meth in range(len(list_meth))]
    
    
    
    # For microalgae
    list_array_mc_sample_micro =[np.array([[0]*len(list_micro_algae_FU) ]*size,dtype="float32") for meth in range(len(list_meth))]
    
    
    
    # For electricty
    list_array_mc_sample_elec =[np.array([[0]*len(list_electricity_CNTR_FU) ]*size,dtype="float32") for meth in range(len(list_meth))]
    
    # For natural gas
    
    list_array_mc_sample_natural_gas =[np.array([[0]*len(list_natural_gas_CNTR_FU) ]*size,dtype="float32") for meth in range(len(list_meth))]
    
    
    list_FU_combined_mc = list_micro_algae_FU + list_fish_FU + list_electricity_CNTR_FU + list_natural_gas_CNTR_FU
    list_FU_combined_names_mc = list_micro_algae_names + activities_fish_background + list_electricity_CNTR_names + list_natural_gas_CNTR_names
    
        
        
    # As this function is called in parallel, we need to initiliazie brightway again for each worker
    
    bw.projects.set_current('AH_38_naturalgas_2_0711') 

    
    
# Loading Ecoinvent
    Ecoinvent = bw.Database('ecoinvent 3.8 conseq')


# Loading foreground database

    FISHMIC = bw.Database('AH_combi_1')



    biosphere=bw.Database('biosphere3')



    MICAH = bw.Database('Micro_for_combi')
    
    
    #Initialize a MCLCA object with any FU

    mc=bw.MonteCarloLCA(list_fish_FU[2],('ReCiPe Midpoint (I)', 'climate change', 'GWP20'))
    
    #print("done_initializing")
    
    time1=time()
    for it in range(size):
        
        next(mc) # As the seed is different for each worker, the "next(mc)" is not the same for all workers, which is what we want 
        
        # We iterate over the whole list of functional units but keep track of
        # the position so that results are stored in the right positions in the disctinct list of specific inputs
        # These lists are used differently in the algorithm.
        
        for i in range(0,len(list_FU_combined_mc)): 
            
            if list_FU_combined_names_mc[i] not in list_processes_electricity_micro :
              
                
                
                mc.redo_lcia(list_FU_combined_mc[i])  # redo with new FU
                
                index_array_method=-1
                
                for m in list_meth:
                
        
                    index_array_method+=1
                
                     # Store at the right place
                     
                    if i<len(list_micro_algae_FU):
                        
                        # This calculates the impact for 1 input, for 1 meth    

                        list_array_mc_sample_micro[index_array_method][it,i]=(C_matrixes[m]*mc.inventory).sum()
                        
                    elif len(list_micro_algae_FU) <=i< len(list_micro_algae_FU + list_fish_FU):
        
                        # This calculates the impact for 1 input, for 1 meth    

                        list_array_mc_sample_fish[index_array_method][it,i-len(list_micro_algae_FU)]=(C_matrixes[m]*mc.inventory).sum()
                    
    
                    elif len(list_micro_algae_FU + list_fish_FU) <=i< len(list_micro_algae_FU + list_fish_FU+ list_electricity_CNTR_FU):
                        # This calculates the impact for 1 input, for 1 meth    

                        list_array_mc_sample_elec[index_array_method][it,i-len(list_micro_algae_FU)-len(list_fish_FU)]=(C_matrixes[m]*mc.inventory).sum()
        
                    else:
                        # This calculates the impact for 1 input, for 1 meth    

                        list_array_mc_sample_natural_gas[index_array_method][it,i-len(list_micro_algae_FU)-len(list_fish_FU)-len(list_electricity_CNTR_FU)]=(C_matrixes[m]*mc.inventory).sum()
    
    
    
    # Concatenate /reorganize results
    list_array_total_mc=[]
    for index_meth in range(len(list_meth)):
        concat = np.concatenate([list_array_mc_sample_micro[index_meth],list_array_mc_sample_fish[index_meth],list_array_mc_sample_elec[index_meth],list_array_mc_sample_natural_gas[index_meth]], axis= 1)
        list_array_total_mc.append(concat)
              
        
       
    
    
    list_array_total_mc_sorted =[[[0 for a in range(list_array_total_mc[0].shape[1])] for meth in list_meth] for it in range(size)]
    
    
    for it in range(size):
        row_to_add = [meth[it] for meth in list_array_total_mc]
        #print(row_to_add)
        #print(row_to_add[0])
        #print(len(row_to_add[]))
        #print(len(row_to_add))
        
        list_array_total_mc_sorted[it] = row_to_add
            
    return list_array_total_mc_sorted   
    




"""3) RUN functions and generate the MC iterations in the background"""




number_cores = multiprocessing.cpu_count()



ray.shutdown()

ray.init()
size_divis= 1000   # Each worker performs 1000 MC as one task
number_chunks = 168 # The workers have 168 tasks to share

size= size_divis * number_chunks # Total number of background MC

#The inputs that are constant to all tasks/workers
constant_inputs = ray.put([list_micro_algae_FU,
                           list_fish_FU,
                           list_electricity_CNTR_FU,
                           list_natural_gas_CNTR_FU,
                           list_micro_algae_names,
                           activities_fish_background,
                           list_electricity_CNTR_names,
                           list_natural_gas_CNTR_names,
                           C_matrixes,
                           size_divis])


time3=time()
# Simulate
ray_results = ray.get([parallel_MC_backgound.remote(constant_inputs) for i in range(number_chunks )])

time4=time()

print("timetot_parallel",time4-time3)




list_array_total_mc_sorted = functools.reduce(operator.iconcat, ray_results, [])



x = datetime.datetime.now()

month=str(x.month)
day=str(x.day)
microsec=str(x.strftime("%f"))
             


name_file_background ='montecarlo_background'+"_"+month+"_"+day+"_"+microsec+"size="+str(size)


# Export results as a pickle object 

export_pickle_2(list_array_total_mc_sorted, name_file_background, "Background_mc")
