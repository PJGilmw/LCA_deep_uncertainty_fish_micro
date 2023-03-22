# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 10:02:21 2022

@author: Pierre Jouannais, Department of Planning, DCEA, Aalborg University
pijo@plan.aau.dk

"""
'''
This script was used to combine different chunks of backrogund MC results to reach certain sizes 
that better match the computing performances of the instances that will use these chunks to simulate the LCA results.
It is not necessariæy needed if the chunks already have adequate sizes.
    
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


import pickle


def importpickle(path):
    with open(path, 'rb') as pickle_load:
        obj = pickle.load(pickle_load)
    return obj    



def export_pickle_2(var, name_var, namefolder_in_root):
    '''Saves a pickle in the working directory and
    saves the object in input across python sessions'''

    path_object = "../"+namefolder_in_root+"/"+name_var+".pkl"
    with open(path_object, 'wb') as pickle_file:
        pickle.dump(var, pickle_file, protocol=pickle.HIGHEST_PROTOCOL)





back_mc_1=importpickle("../Background_mc/500000/montecarlo_background_12_13_311955size=235000.pkl")
back_mc_2=importpickle("../Background_mc/500000/montecarlo_background_12_13_574232size=97800.pkl")
back_mc_3=importpickle("../Background_mc/500000/montecarlo_background_12_14_074170size=168000.pkl")

combined= back_mc_1 + back_mc_2 + back_mc_3 
   

size=len(combined)

x = datetime.datetime.now()

month=str(x.month)
day=str(x.day)
microsec=str(x.strftime("%f"))



name_file='Combined_Background_mc'+"_"+month+"_"+day+"_"+microsec+"_size="+str(size)


export_pickle_2(combined, name_file, "Background_mc")

size1=266800
size2=104832
size3=128320  
           
chunk_1= combined[:size1]
chunk_2= combined[size1:size1+size2]
chunk_3= combined[size1+size2:size1+size2+size3]

name_file_1='chunk_Background_mc'+"_"+month+"_"+day+"_"+microsec+"_size="+str(size1)
name_file_2='chunk_Background_mc'+"_"+month+"_"+day+"_"+microsec+"_size="+str(size2)
name_file_3="chunk_Background_mc"+"_"+month+"_"+day+"_"+microsec+"_size="+str(size3)

export_pickle_2(chunk_1, name_file_1, "Background_mc")
export_pickle_2(chunk_2, name_file_2, "Background_mc")
export_pickle_2(chunk_3, name_file_3, "Background_mc")



#Third chunk was divided in 3 on the 16/12


back_mc_3=importpickle("../Background_mc/chunk_Background_mc_12_15_919097_size=128320.pkl")

len(back_mc_3)/3

sub_size_1=42773
sub_size_2=42773
sub_size_3=42773


sub_chunk_1= back_mc_3[:sub_size_1]
sub_chunk_2= back_mc_3[sub_size_1:sub_size_1+sub_size_2]
sub_chunk_3= back_mc_3[sub_size_1+sub_size_2:sub_size_1+sub_size_2+sub_size_3]


x = datetime.datetime.now()

month=str(x.month)
day=str(x.day)
microsec=str(x.strftime("%f"))


name_sub_file_1='chunk_Background_mc'+"_"+month+"_"+day+"_"+microsec+"_size="+str(sub_size_1) + "_number1" 
name_sub_file_2='chunk_Background_mc'+"_"+month+"_"+day+"_"+microsec+"_size="+str(sub_size_2) + "_number2"
name_sub_file_3="chunk_Background_mc"+"_"+month+"_"+day+"_"+microsec+"_size="+str(sub_size_3)+ "_number3"

export_pickle_2(sub_chunk_1, name_sub_file_1, "Background_mc")
export_pickle_2(sub_chunk_2, name_sub_file_2, "Background_mc")
export_pickle_2(sub_chunk_3, name_sub_file_3, "Background_mc")



sub_size_1+sub_size_2+sub_size_3

size2=104832
size3=128320 