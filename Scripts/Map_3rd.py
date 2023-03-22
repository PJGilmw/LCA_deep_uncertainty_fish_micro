# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 11:56:25 2022

@author: Pierre Jouannais, Department of Planning, DCEA, Aalborg University
pijo@plan.aau.dk

"""
'''
Script containing the functions generating and handling the geographic files/locations grid/geodataframes.


'''



import datetime
from time import *
import requests
import pickle
import cProfile
from scipy.integrate import odeint

import os

#Set working directory to file location 
#(works only when executing the whole file and not only sections (Run Current cell))

currentfolder=os.path.dirname(os.path.realpath(__file__))
os.chdir(currentfolder)

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
import SALib

import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits import mplot3d

import shapefile as shp
import geopandas as gpd
#import pysal as ps
from shapely.geometry import Polygon, mapping, Point


import Retrieving_solar_and_climatic_data_2nd as solardata



    



def read_shapefile(sf):
    """
    Read a shapefile into a Pandas dataframe with a 'coords' 
    column holding the geometry information. This uses the pyshp
    package
    """
    fields = [x[0] for x in sf.fields][1:]
    records = sf.records()
    shps = [s.points for s in sf.shapes()]
    df = pd.DataFrame(columns=fields, data=records)
    df = df.assign(coords=shps)
    
    return df





def generate_random(number, list_polygon):
    """Function which generates "number" random points in a polygon given as input and 
    returns the shapely objects and the coordinates"""

    points = []
    polygon=Polygon(list_polygon)
    minx, miny, maxx, maxy = polygon.bounds
    while len(points) < number:
        pnt = Point(np.random.uniform(minx, maxx), np.random.uniform(miny, maxy))
        if polygon.contains(pnt):
            points.append(pnt)
            
    coords_points=[]        
    for p in points:
        print(p.x,",",p.y)     
        coords_points.append((p.x,p.y))
    return points,coords_points





def generate_random_latsec_corrected(list_polygon,latsec):
    
    """Function which generates X semi-random points in a polygon given as input and 
    returns the shapely objects and the coordinates. The points are separated on the latitude axis 
    by a same distance "latsec". the number of points generated per country is therefore :
        Country's latitudinal spread//latsec """
    points = []
    polygon=Polygon(list_polygon)
    minx, miny, maxx, maxy = polygon.bounds
    
    # Stop generating points above polar circle
    if maxy>65:
        maxy=65
    #print(minx,maxx)
    print(miny, maxy)
    # Number points 
    number=(maxy-miny-0.5)//latsec
    print('number=',number)
    
    
    
    # first point in the south
    while len(points) < 1:
        ypoint= miny+0.5
        #print(ypoint)
        pnt = Point(np.random.uniform(minx, maxx), ypoint)
        
        # We also check if the coordinates can be found in PVGIS (some points on the coast may be right on the limit)
        if polygon.contains(pnt) and solardata.downloadclimaticdata(round(ypoint,3), round(pnt.x,3), 5, 90, 90) :
            print('ok')
            points.append(pnt)
        else:
            print('outside')
    
      
    ypoint=ypoint+latsec
    
    while len(points) < number:
        
        #print(ypoint)
        a=np.random.uniform(minx, maxx)
        #print(a)
        pnt = Point(a, ypoint)
        
        if polygon.contains(pnt) and solardata.downloadclimaticdata(round(ypoint,3), round(pnt.x,3), 5, 90, 90):
            
            ypoint=ypoint+latsec
            
            #print(ypoint)
            points.append(pnt)
            
        else:
            print('outside') 
            
    coords_points=[]        
    for p in points:
        #print(p.x,",",p.y)     
        coords_points.append((p.x,p.y))
    return points,coords_points






# Put in a function

def geodataframe_initialize(latsec,methods_selected,biosph,MICAH,Ecoinvent ):
    """ Function generating the grid of locations and the associated dataframe. Also collects the electricity mix of each country."""
    
    np.random.seed(1) # The same map is always created 
    
    shp_path = "../Map/Eurostat/NUTS_RG_20M_2021_4326.shp"
    
    sf = shp.Reader(shp_path)
    
    
    # Convert shapefile into dataframe
    
    shapefile_df = read_shapefile(sf)
    
    
    # Now open with geopandas
    
    geodataframe = gpd.read_file(shp_path)   
    
    
    
    
    
    # Keep only the top 10 countries from Skarka 2012 inin the geodataframe
    
    # Top 10 Skarka
    
    list_skarka=['ES','SE','IT','PT','UK','FR','EL','CY','IE','DE']  
    
    #List index to keep
    
    list_index=[a for a in range(geodataframe.shape[0]) if geodataframe.iloc[a]['CNTR_CODE'] in list_skarka]
    
    geodataframe_top=geodataframe.iloc[list_index]
    
    # for shapefile_df
     
    list_index_df=[a for a in range(shapefile_df.shape[0]) if shapefile_df.iloc[a]['CNTR_CODE'] in list_skarka]
    
    shapefile_df_top=shapefile_df.iloc[list_index]
    
    
    
    #Keep only countries level
    
    geodataframe_top_cntr = geodataframe_top[geodataframe_top.LEVL_CODE==0]
    shapefile_df_top_cntr = shapefile_df_top[shapefile_df_top.LEVL_CODE==0]
    
        
    # Reset indexes
    
    # We are going to work with countries only
    
    geodataframe_top_cntr = geodataframe_top_cntr.reset_index(drop=True)
    shapefile_df_top_cntr = shapefile_df_top_cntr.reset_index(drop=True)
    
    
    # Keeping only metropolitan territories
    
    
    geodataframe_metrop = copy.deepcopy(geodataframe_top_cntr)
    shapefile_df_metrop = copy.deepcopy(shapefile_df_top_cntr)
    
    # We make sure that the points will be drawn on the metropolitan areas
    # While the column geometry is used for plotting
    # The sampling is made on the column "coords".
    # We only keep the coordinates of the metropolitan area in "coords".
    
    # We do this for each country :
    
    
    
    # Create a column coords that is updated later
    
    geodataframe_metrop['coords'] = shapefile_df_metrop['coords']
    
    
    # France : keep only metropolitan France
    
    # Coords
    polygon = mapping(geodataframe_metrop['geometry'][0][3])
    
    coord_list = list(polygon['coordinates'][0])
    
    geodataframe_metrop['coords'][0] = coord_list
    
    # We will only plot metropolitan France
    
    geodataframe_metrop['geometry'][0] = geodataframe_metrop['geometry'][0][3]
    
    
    # Germany ok
    
    
    # Greece change only the sampling area to metropolitan
    
    # Coords
    polygon = mapping(geodataframe_metrop['geometry'][3][12])
    
    coord_list=list(polygon['coordinates'][0])
    
    geodataframe_metrop['coords'][3]=coord_list
    
    
    
    # Ireland ok
    
    # Spain no
    
    # Keep only metropolitan Spain
    
    # Coords
    polygon = mapping(geodataframe_metrop['geometry'][5][2])
    
    coord_list=list(polygon['coordinates'][0])
    
    geodataframe_metrop['coords'][5]=coord_list
    
    # We will only plot metropolitan Spain
    
    
    geodataframe_metrop['geometry'][5]=geodataframe_metrop['geometry'][5][2]
    
    
    
    # Portugal 
    
    # keep only metropolitan portugal
    
    # Coords
    polygon = mapping(geodataframe_metrop['geometry'][6][0])
    
    coord_list=list(polygon['coordinates'][0])
    
    geodataframe_metrop['coords'][6]=coord_list
    
    # Plot only metropolitan Portugal
    
    geodataframe_metrop['geometry'][6]=geodataframe_metrop['geometry'][6][0]
    
    
    # Sweden change only drawing coords but plot everything
    
    # Coords
    polygon = mapping(geodataframe_metrop['geometry'][7][1])
    
    coord_list=list(polygon['coordinates'][0])
    
    geodataframe_metrop['coords'][7]=coord_list
    
    
    
    
    #UK change only drawing coords but plot everything
    
    # Coords
    
    polygon = mapping(geodataframe_metrop['geometry'][8][0])
    
    coord_list=list(polygon['coordinates'][0])
    
    geodataframe_metrop['coords'][8]=coord_list
    
    
    # Italy 
    
    # Coords
    polygon = mapping(geodataframe_metrop['geometry'][9][0])
    
    coord_list=list(polygon['coordinates'][0])
    
    geodataframe_metrop['coords'][9]=coord_list
    
    # Keep only metropolitan italy
    
    # We will only plot metropolitan Italy
    
    geodataframe_metrop['geometry'][9]=geodataframe_metrop['geometry'][9][0]
    
    
    
    ### Generate X random points in each region point
    
    
    
    # Necessary to initialize column that will receive the generated points
    
    geodataframe_metrop['random_points']=0
    
    for region_index in range(geodataframe_metrop.shape[0]):
        
        geodataframe_metrop['random_points'].iloc[region_index]=np.array([5,6,7],dtype=object)
        
    
    # Generate random points
    
    # With semi random generation along the longitudinal axis
    
    total_generated_points=[]
    
    for region_index in range(geodataframe_metrop.shape[0]):

        
        generated_points,coords_points= generate_random_latsec_corrected(geodataframe_metrop['coords'].iloc[region_index],latsec) 
        
        
        total_generated_points+=generated_points
    
        #Coordinates for LCA input
        
        geodataframe_metrop['random_points'].iloc[region_index]=coords_points
        
        
        
    #Add these points to a geodataframe
        
    index=[i for i in range(len(total_generated_points))]
    
    gdf_points = gpd.GeoDataFrame(index,geometry=total_generated_points, crs=4326)
    
    
    

    
    
    
    ##### ADDING ECOINVENT ELECTRICITY MIXES
    
    # List all ISO country codes:
        
    listcountries=geodataframe_metrop['CNTR_CODE'].unique()
    
    

    
    list_act=[]
    list_loc=[]
    dict_elec={}
    
    for act in Ecoinvent:
        if 'market for electricity' in act['name'] and 'medium voltage' in act['name'] and 'municipal' not in act['name'] and act['location'] in listcountries and 'label-certified' not in act['name']:
            list_act.append(act)
            list_loc.append(act['location'])
            
    

    
    #Add column with list of impacts for the national electiricity mixes
    
    geodataframe_metrop['Imp_elec']=0
    
    # This column will contains arrays
    for region_index in range(geodataframe_metrop.shape[0]):
        
        geodataframe_metrop['Imp_elec'].iloc[region_index]=np.array([5,6,7],dtype=object)
        
    
    list_FU_elec=[]
    for index_row in range(geodataframe_metrop.shape[0]):
        
        cntr=geodataframe_metrop['CNTR_CODE'][index_row]
        
        #conversion ISO code to ecoinvent code
        
        if cntr=='UK':
            cntr='GB'
        elif cntr=='EL':
            cntr='GR'
            
        for act in Ecoinvent:
            if 'market for electricity' in act['name']:
                if 'medium voltage' in act['name']:
                    if 'municipal' not in act['name']:
                        if 'label-certified' not in act['name']:
                            if act['location']==cntr:
                                print(act)  
                                list_FU_elec.append({act:1})
                             
    
    my_calculation_setup_elec = {'inv': list_FU_elec, 'ia': methods_selected}
       
                             
    bw.calculation_setups['elec_mixes'] = my_calculation_setup_elec
    
    mlca_elec = bw.MultiLCA('elec_mixes')  
    
    res_elec = mlca_elec.results
    
    
    # add results to the column
    # ranked in the same order
    for row_index in range(geodataframe_metrop.shape[0]):
        
        geodataframe_metrop['Imp_elec'].iloc[row_index] = np.array(res_elec[row_index],dtype=object)
        
        #geodataframe_metrop['Imp_elec'].iloc[row_index]=np.array([5,6,7],dtype=object)       
    
    
    
        
        # create columns that will host standard deviation and average of impacts
    for meth_index in range(len(methods_selected)):
        
        name_meth = methods_selected[meth_index][-1]
        
        name_col_std='Std'+'_'+name_meth
        
        name_col_mean='Mean'+'_'+name_meth
        
        name_col_min='Min'+'_'+name_meth
        
        name_col_max='Max'+'_'+name_meth
    
        
        geodataframe_metrop[name_col_std]=0
        
        geodataframe_metrop[name_col_mean]=0
        
        geodataframe_metrop[name_col_min]=0
        
        geodataframe_metrop[name_col_max]=0
        
        
        gdf_points[name_col_max] = 0
        
        gdf_points[name_col_min] = 0
        
        gdf_points[name_col_mean]=0
        
    
    
    return geodataframe_metrop, gdf_points






# Put in a function

def geodataframe_initialize_2(latsec,methods_selected,biosph,MICAH,Ecoinvent ):
    """ Function generating the grid of locations and the associated dataframe. Also collects the electricity mix of each country."""
    
    np.random.seed(1) # The same map is always created 
    
    shp_path = "../Map/Eurostat/NUTS_RG_20M_2021_4326.shp"
    
    sf = shp.Reader(shp_path)
    
    
    # Convert shapefile into dataframe
    
    shapefile_df = read_shapefile(sf)
    
    
    # Now open with geopandas
    
    geodataframe = gpd.read_file(shp_path)   
    
    
    
    
    
    # Keep only the top 10 countries from Skarka 2012 inin the geodataframe
    
    # Top 10 Skarka
    
    list_skarka=['ES','SE','IT','PT','UK','FR','EL','CY','IE','DE']  
    
    #List index to keep
    
    list_index=[a for a in range(geodataframe.shape[0]) if geodataframe.iloc[a]['CNTR_CODE'] in list_skarka]
    
    geodataframe_top=geodataframe.iloc[list_index]
    
    # for shapefile_df
     
    list_index_df=[a for a in range(shapefile_df.shape[0]) if shapefile_df.iloc[a]['CNTR_CODE'] in list_skarka]
    
    shapefile_df_top=shapefile_df.iloc[list_index]
    
    
    
    #Keep only countries level
    
    geodataframe_top_cntr = geodataframe_top[geodataframe_top.LEVL_CODE==0]
    shapefile_df_top_cntr = shapefile_df_top[shapefile_df_top.LEVL_CODE==0]
    
        
    # Reset indexes
    
    # We are going to work with countries only
    
    geodataframe_top_cntr = geodataframe_top_cntr.reset_index(drop=True)
    shapefile_df_top_cntr = shapefile_df_top_cntr.reset_index(drop=True)
    
    
    # Keeping only metropolitan territories
    
    
    geodataframe_metrop = copy.deepcopy(geodataframe_top_cntr)
    shapefile_df_metrop = copy.deepcopy(shapefile_df_top_cntr)
    
    # We make sure that the points will be drawn on the metropolitan areas
    # While the column geometry is used for plotting
    # The sampling is made on the column "coords".
    # We only keep the coordinates of the metropolitan area in "coords".
    
    # We do this for each country :
    
    
    
    # Create a column coords that is updated later
    
    geodataframe_metrop['coords'] = shapefile_df_metrop['coords']
    
    
    # France : keep only metropolitan France
    
    # Coords
    polygon = mapping(geodataframe_metrop['geometry'][0][3])
    
    coord_list = list(polygon['coordinates'][0])
    
    geodataframe_metrop['coords'][0] = coord_list
    
    # We will only plot metropolitan France
    
    geodataframe_metrop['geometry'][0] = geodataframe_metrop['geometry'][0][3]
    
    
    # Germany ok
    
    
    # Greece change only the sampling area to metropolitan
    
    # Coords
    polygon = mapping(geodataframe_metrop['geometry'][3][12])
    
    coord_list=list(polygon['coordinates'][0])
    
    geodataframe_metrop['coords'][3]=coord_list
    
    
    
    # Ireland ok
    
    # Spain no
    
    # Keep only metropolitan Spain
    
    # Coords
    polygon = mapping(geodataframe_metrop['geometry'][5][2])
    
    coord_list=list(polygon['coordinates'][0])
    
    geodataframe_metrop['coords'][5]=coord_list
    
    # We will only plot metropolitan Spain
    
    
    geodataframe_metrop['geometry'][5]=geodataframe_metrop['geometry'][5][2]
    
    
    
    # Portugal 
    
    # keep only metropolitan portugal
    
    # Coords
    polygon = mapping(geodataframe_metrop['geometry'][6][0])
    
    coord_list=list(polygon['coordinates'][0])
    
    geodataframe_metrop['coords'][6]=coord_list
    
    # Plot only metropolitan Portugal
    
    geodataframe_metrop['geometry'][6]=geodataframe_metrop['geometry'][6][0]
    
    
    # Sweden change only drawing coords but plot everything
    
    # Coords
    polygon = mapping(geodataframe_metrop['geometry'][7][1])
    
    coord_list=list(polygon['coordinates'][0])
    
    geodataframe_metrop['coords'][7]=coord_list
    
    
    
    
    #UK change only drawing coords but plot everything
    
    # Coords
    
    polygon = mapping(geodataframe_metrop['geometry'][8][0])
    
    coord_list=list(polygon['coordinates'][0])
    
    geodataframe_metrop['coords'][8]=coord_list
    
    
    # Italy 
    
    # Coords
    polygon = mapping(geodataframe_metrop['geometry'][9][0])
    
    coord_list=list(polygon['coordinates'][0])
    
    geodataframe_metrop['coords'][9]=coord_list
    
    # Keep only metropolitan italy
    
    # We will only plot metropolitan Italy
    
    geodataframe_metrop['geometry'][9]=geodataframe_metrop['geometry'][9][0]
    
    
    
    ### Generate X random points in each region point
    
    
    
    # Necessary to initialize column that will receive the generated points
    
    geodataframe_metrop['random_points']=0
    
    for region_index in range(geodataframe_metrop.shape[0]):
        
        geodataframe_metrop['random_points'].iloc[region_index]=np.array([5,6,7],dtype=object)
        
    
    # Generate random points
    
    # With semi random generation along the longitudinal axis
    
    total_generated_points=[]
    
    for region_index in range(geodataframe_metrop.shape[0]):

        
        generated_points,coords_points= generate_random_latsec_corrected(geodataframe_metrop['coords'].iloc[region_index],latsec) 
        
        
        total_generated_points+=generated_points
    
        #Coordinates for LCA input
        
        geodataframe_metrop['random_points'].iloc[region_index]=coords_points
        
        
        
    #Add these points to a geodataframe
        
    index=[i for i in range(len(total_generated_points))]
    
    gdf_points = gpd.GeoDataFrame(index,geometry=total_generated_points, crs=4326)
    
    
    

    
    
    
    ##### ADDING ECOINVENT ELECTRICITY MIXES
    
    # List all ISO country codes:
        
    listcountries=geodataframe_metrop['CNTR_CODE'].unique()
    
    

    
    list_act=[]
    list_loc=[]
    dict_elec={}
    
    for act in Ecoinvent:
        if 'market for electricity' in act['name'] and 'medium voltage' in act['name'] and 'municipal' not in act['name'] and act['location'] in listcountries and 'label-certified' not in act['name']:
            list_act.append(act)
            list_loc.append(act['location'])
            
    

    
    #Add column with list of impacts for the national electiricity mixes
    
    geodataframe_metrop['Imp_elec']=0
    
    # This column will contains arrays
    for region_index in range(geodataframe_metrop.shape[0]):
        
        geodataframe_metrop['Imp_elec'].iloc[region_index]=np.array([5,6,7],dtype=object)
        
    
    list_FU_elec=[]
    for index_row in range(geodataframe_metrop.shape[0]):
        
        cntr=geodataframe_metrop['CNTR_CODE'][index_row]
        
        #conversion ISO code to ecoinvent code
        
        if cntr=='UK':
            cntr='GB'
        elif cntr=='EL':
            cntr='GR'
            
        for act in Ecoinvent:
            if 'market for electricity' in act['name']:
                if 'medium voltage' in act['name']:
                    if 'municipal' not in act['name']:
                        if 'label-certified' not in act['name']:
                            if act['location']==cntr:
                                print(act)  
                                list_FU_elec.append({act:1})
                             
    

    
    # add results to the column
    # ranked in the same order
    for row_index in range(geodataframe_metrop.shape[0]):
        
        geodataframe_metrop['Imp_elec'].iloc[row_index] = np.array(list_FU_elec[row_index],dtype=object)
        
        #geodataframe_metrop['Imp_elec'].iloc[row_index]=np.array([5,6,7],dtype=object)       
    
            
        
    list_points_grid = []
    for row_index in range(geodataframe_metrop.shape[0]):
        #print(geodataframe_metrop.loc[row_index,["random_points"]][0])
        
        list_points_grid.append([geodataframe_metrop.loc[row_index,["random_points"]][0],geodataframe_metrop.loc[row_index,["CNTR_CODE"]][0]])
    
    
    list_points_grid_disaggregated = []
    
    for group in list_points_grid:
        if group[1] == "UK":   # conversion to ecoinvent codes
            group[1] = "GB"
        elif group[1] == "EL":   # conversion to ecoinvent codes
            group[1] = "GR"

        for point in group[0]:
            
            
            list_points_grid_disaggregated.append([(round(point[0],3),round(point[1],3)),group[1]])
        

    
    return geodataframe_metrop, gdf_points, list_points_grid_disaggregated

