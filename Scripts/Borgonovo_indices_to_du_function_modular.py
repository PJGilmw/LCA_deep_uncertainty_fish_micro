# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 09:46:26 2023


@author: Pierre Jouannais, Department of Planning, DCEA, Aalborg University
pijo@plan.aau.dk

"""
'''
Script containing the function which assigns dU to parameters/indicators based on their ẟ


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

import matplotlib.pyplot as plt







def get_du_from_borgonovo_indices_modular(borgonovo_indices,
                                  combined_results_head,
                                  impact_category,
                                  type_factor,
                                  dumax,
                                  dumin,
                                  uncertain_param,
                                  columns_to_remove):
    """Calculate dUs for all parameters/indicators based on their ẟ (input).
    Plots the dUs in borgonovo_indices"""
    
    # Set the axes labels font size
    plt.rc('axes', labelsize=17)
    # Set the font size for x tick labels
    plt.rc('xtick', labelsize=13)
    # Set the font size for y tick labels
    plt.rc('ytick', labelsize=13)   
        
    
    
    x = datetime.datetime.now()
    
    month=str(x.month)
    day=str(x.day)
    microsec=str(x.strftime("%f"))
                 
    
    # Keep only the parameters/indicators 
    columns_indices = [col for col in list(combined_results_head.columns) if not col in columns_to_remove]
    
    
    borgonovo_indices.columns = columns_indices
    
    # Uncertain parameters

    borgonovo_indices_uncertain = borgonovo_indices[uncertain_param]
    
    
    borgonovo_indices_T = borgonovo_indices.T
    borgonovo_indices_uncertain_T = borgonovo_indices_uncertain.T 
    
    borgonovo_indices_T.reset_index(inplace=True)
    
    borgonovo_indices_uncertain_T.reset_index(inplace=True)
    
    borgonovo_indices_T.columns = ["Parameter","𝛿"]
    borgonovo_indices_uncertain_T.columns = ["Parameter","𝛿"]
    
    borgonovo_indices_T.sort_values(by=["𝛿"],inplace=True)
    borgonovo_indices_uncertain_T.sort_values(by=["𝛿"],inplace=True)
    
    
    
    ax_all = borgonovo_indices_T.plot.bar(x='Parameter', y="𝛿", rot=90)
    
    ax_all.set_xticklabels([])
    ax_all.set_ylabel("Delta")
    ax_all.get_legend().remove()

    fig_all= ax_all.get_figure()

    name_fig_all = "all_indices"+impact_category+"_"+ month+"_"+day+"_"+microsec
   
    fig_all.savefig("../borgonovo_indices/"+name_fig_all,dpi=500,bbox_inches="tight")

 
    
    #du= f(1/indice) = a*(1/indice) + b
    
    
    def du_from_indice(dumax,dumin,min_indice,max_indice, indice):
        """Function which assign the du to a parameter/indicator"""
        
        a = (dumax - dumin ) / ( (1/min_indice)-(1/max_indice))
        
        b = dumax - a*1/min_indice
        du_raw = a*(1/indice) + b
        du_round=1/round(1/du_raw)
    
        
        return [du_raw, du_round]
    
    
    max_indice = max(list(borgonovo_indices_uncertain.iloc[0]))
    min_indice = min(list(borgonovo_indices_uncertain.iloc[0]))

    
    
    list_indice =  list(borgonovo_indices_uncertain.iloc[0])
    list_names = list(borgonovo_indices_uncertain.columns)
    list_du=[]
    
    
    
    for indice in list_indice:
        
        du=du_from_indice(dumax,dumin,min_indice,max_indice, indice)
        #print(du[1])
        list_du.append(du[1])
        
        
    
     
    dataframe_indices = pd.DataFrame(np.array([list_names,list_du]))
    
    dataframe_indices_T = dataframe_indices.T
    
    dataframe_indices_T.columns = ["Parameter", "dU"]
    
    dataframe_indices_T["dU"] = pd.to_numeric(dataframe_indices_T["dU"])
    
  
    
    # Plot combined du and "𝛿" together
    dataframe_indices_combined = pd.DataFrame(np.array([list_names,list_du,list_indice]))

    dataframe_indices_combined_T = dataframe_indices_combined.T

    dataframe_indices_combined_T.columns = ["Parameter", "dU","𝛿"]

    dataframe_indices_combined_T["dU"] = pd.to_numeric(dataframe_indices_combined_T["dU"])
    dataframe_indices_combined_T["𝛿"] = pd.to_numeric(dataframe_indices_combined_T["𝛿"])
    
    dataframe_indices_combined_T.index = dataframe_indices_combined_T["Parameter"]
    
    ax_combined= dataframe_indices_combined_T.plot.bar(rot=90, subplots=True,title=['', ''],legend=None)
  
    ax_combined[0].set_ylabel("dU")
    ax_combined[1].set_ylabel(r'$\delta$')


    
    fig_combined= ax_combined[0].get_figure()



    name_fig_combined = "combined"+impact_category+"_"+ type_factor+"_"+ month+"_"+day+"_"+microsec
   

    fig_combined.savefig("../borgonovo_indices/"+name_fig_combined,dpi=500,bbox_inches="tight")
    
    
    
    
    
    return dataframe_indices_combined_T
    





def get_du_from_borgonovo_indices_modular_order(borgonovo_indices,
                                  combined_results_head,
                                  impact_category,
                                  type_factor,
                                  dumax,
                                  dumin,
                                  uncertain_param,
                                  columns_to_remove):
    """Same as above but plots dUs and "𝛿" in ascending order"""
    
    # Set the axes labels font size
    plt.rc('axes', labelsize=17)
    # Set the font size for x tick labels
    plt.rc('xtick', labelsize=13)
    # Set the font size for y tick labels
    plt.rc('ytick', labelsize=13)   
        
    
    
    x = datetime.datetime.now()
    
    month=str(x.month)
    day=str(x.day)
    microsec=str(x.strftime("%f"))
                 
    



    

    
    columns_indices = [col for col in list(combined_results_head.columns) if not col in columns_to_remove]
    
    
    borgonovo_indices.columns = columns_indices
    

    #Prepare dataframes
    
    borgonovo_indices_uncertain = borgonovo_indices[uncertain_param]
    
    
    borgonovo_indices_T = borgonovo_indices.T
    borgonovo_indices_uncertain_T = borgonovo_indices_uncertain.T 
    
    borgonovo_indices_T.reset_index(inplace=True)
    
    borgonovo_indices_uncertain_T.reset_index(inplace=True)
    
    borgonovo_indices_T.columns = ["Parameter","𝛿"]
    borgonovo_indices_uncertain_T.columns = ["Parameter","𝛿"]
    
    borgonovo_indices_T.sort_values(by=["𝛿"],inplace=True)
    borgonovo_indices_uncertain_T.sort_values(by=["𝛿"],inplace=True)
    
    
    
    ax_all = borgonovo_indices_T.plot.bar(x='Parameter', y="𝛿", rot=90)
    
    ax_all.set_xticklabels([])
    ax_all.set_ylabel("Delta")
    ax_all.get_legend().remove()

    fig_all= ax_all.get_figure()

    name_fig_all = "all_indices"+impact_category+"_"+ month+"_"+day+"_"+microsec
   
    fig_all.savefig("../borgonovo_indices/"+name_fig_all,dpi=500,bbox_inches="tight")


    def du_from_indice(dumax,dumin,min_indice,max_indice, indice):
        
        a = (dumax - dumin ) / ( (1/min_indice)-(1/max_indice))
        
        b = dumax - a*1/min_indice
        du_raw = a*(1/indice) + b
        du_round=1/round(1/du_raw)
    
        
        return [du_raw, du_round]
    
    
    max_indice = max(list(borgonovo_indices_uncertain.iloc[0]))
    min_indice = min(list(borgonovo_indices_uncertain.iloc[0]))

    
    
    list_indice =  list(borgonovo_indices_uncertain.iloc[0])
    list_names = list(borgonovo_indices_uncertain.columns)
    list_du=[]
    
    
    
    for indice in list_indice:
        
        du=du_from_indice(dumax,dumin,min_indice,max_indice, indice)
        #print(du[1])
        list_du.append(du[1])
        
        
    
     
    dataframe_indices = pd.DataFrame(np.array([list_names,list_du]))
    
    dataframe_indices_T = dataframe_indices.T
    
    dataframe_indices_T.columns = ["Parameter", "dU"]
    
    dataframe_indices_T["dU"] = pd.to_numeric(dataframe_indices_T["dU"])
    
  
    
    # Plot combined
    dataframe_indices_combined = pd.DataFrame(np.array([list_names,list_indice,list_du]))

    dataframe_indices_combined_T = dataframe_indices_combined.T

    dataframe_indices_combined_T.columns = ["Parameter", "𝛿","dU"]

    dataframe_indices_combined_T["dU"] = pd.to_numeric(dataframe_indices_combined_T["dU"])
    dataframe_indices_combined_T["𝛿"] = pd.to_numeric(dataframe_indices_combined_T["𝛿"])
    
    dataframe_indices_combined_T.index = dataframe_indices_combined_T["Parameter"]

    ax_combined= dataframe_indices_combined_T.sort_values('𝛿', ascending=False).plot.bar(rot=90, subplots=True,title=['', ''],legend=None)
  
    # ax_du.set_ylabel("dU")"λ"
    ax_combined[0].set_ylabel(r'$\delta$')
    ax_combined[1].set_ylabel("dUre")

    # # ax_du.get_legend().remove()
    #  ax_combined[0].get_legend().remove()
    #  ax_combined[1].get_legend().remove()
    
    fig_combined= ax_combined[0].get_figure()



    name_fig_combined = "combined"+impact_category+"_"+ type_factor+"_"+ month+"_"+day+"_"+microsec
   

    fig_combined.savefig("../borgonovo_indices/"+name_fig_combined,dpi=500,bbox_inches="tight")
    
    
    
    
    
    return dataframe_indices_combined_T
    

    