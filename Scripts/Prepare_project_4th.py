# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 17:23:30 2021

@author: Pierre Jouannais, Department of Planning, DCEA, Aalborg University
pijo@plan.aau.dk
"""

'''
Create a project 'Fish_project_check' with biosphere3, importing Ecoinvent 3.8 consequential and
the foreground databases imported from the json files.

Change the path to the location of your Ecoinvent files and execute the whole script.
This will take some time.
'''

# Change the path to your Ecoinvent files folder

ei38dir = r"C:/Users/GF20PZ/OneDrive - Aalborg Universitet/Dokumenter/AAU/Databases,softwares/Brightway/Ecoinvent/Ecoinvent 3.8/datasets"


import bw2data
import bw2io
from bw2data.parameters import *
import brightway2 as bw
import json
import scipy
import ast

import os
# Set working directory to file location 
# (works only when executing the whole file and not only sections (Run Current cell))

# currentfolder=os.path.dirname(os.path.realpath(__file__))
# os.chdir(currentfolder)


bw.projects.set_current('Fish_and_micro_2103')


bw.databases

bw.bw2setup()




# Import ecoinvent


if 'ecoinvent 3.8 conseq' in bw.databases:
    print("Database has already been imported")
else:
    # Do not change the name
    ei38 = bw.SingleOutputEcospold2Importer(ei38dir, 'ecoinvent 3.8 conseq',use_mp=False) 
    ei38.apply_strategies()
    ei38.statistics()
    ei38.write_database() # This will take some time.


# Import the foreground database


def remap_keys(mapping):
    
    return [{'key':k, 'value': v} for k, v in mapping.iteritems()]




def import_database_data_from_json(jsonfilepath):
      
    with open(jsonfilepath+'.json') as json_file:
        data = json.load(json_file)
      
    database = {ast.literal_eval(k): v for k, v in data.items()}
 
       
    return database   

bw.projects.current

database_dict_1=import_database_data_from_json('../Data/AH_combi_1')
database_dict_2=import_database_data_from_json('../Data/Micro_for_combi')

FISHMIC = bw.Database('AH_combi_1')
MICAH = bw.Database('Micro_for_combi')
                    
MICAH.write(database_dict_2)
                    
FISHMIC.write(database_dict_1)


