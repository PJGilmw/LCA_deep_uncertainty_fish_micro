# -*- coding: utf-8 -*-
"""
Created on Sat May 22 01:11:59 2021

@author: Pierre Jouannais, Department of Planning, DCEA, Aalborg University
pijo@plan.aau.dk
"""


'''Script containing the functions to retrieve solar climatic data from the API 
of the European Photovoltaic Geographical Information System and estimate
 the solar power received by a given PBR.


# API source : https://ec.europa.eu/jrc/en/PVGIS/docs/noninteractive

'''


import requests
import matplotlib.pyplot as plt
import os

# Set working directory to file location 
# (works only when executing the whole file and not only sections (Run Current cell))

currentfolder=os.path.dirname(os.path.realpath(__file__))
os.chdir(currentfolder)

import Functions_for_physical_and_biological_calculations_3rd as functions
import pandas as pd
import decimal
import random
import numpy as np





# The downloaded file's columns are separated by a tabulation

# The downloaded csv has 4 colums without header : Time with format "XX:XX", Global horizontal irradiciance(direct+diffuse), direct, diffuse, temperature (Celsius)
# takes latitude,longitude, month number, angle and azimuth as inputs

def downloadclimaticdata(lat, long, month, angle, azimuth):
    """API importer. Download a csv file from the European Photovoltaic
    Geographical Information System  with hourly climatic data for
    an average day in the given location and month. The file is downloaded in
    the working directory. Irradiance received on a tilted surface with the
    angle given in input. Data is calculated for the period 2005-2015.

    The downloaded csv has 4 colums without header :
    Time with format "XX:XX", Global horizontal irradiance(direct+diffuse),
    direct, diffuse, temperature (Celsius)

    #Inputs:

        #lat : latitute expressed in format ; XX.XXX or X.XXX
        #long : longitude expressed in format ; XX.XXX or X.XXX
        #month : month of the year ; number of the month
        #angle : Tilt angle of the surface  ; 0= horizontal, 90 = vertical
        #azimuth : azimuth of surface ;  180:-180
        

        

        """

    # No need to download again if already here
    if not os.path.exists("../Climatic_data/dailydataforlat="
                          + str(lat)
                          + "andlong="
                          + str(long)
                          + "formonth"
                          + str(month)
                          + "forangle"
                          + str(angle)
                          + "forazimuth"
                          + str(azimuth)
                          + ".csv"):

        # Create the name of the link from which data will be retrieved
        namelink = ("https://re.jrc.ec.europa.eu/api/DRcalc?lat="
                    + str(lat)
                    + "&lon="
                    + str(long)
                    + "&month="
                    + str(month)
                    + "&global=1&angle="
                    + str(angle)
                    + "&showtemperatures=1&aspect="
                    + str(azimuth)
                    + "&outputformat=basic")

        response = requests.get(namelink)

        # Code 200 means the connection with the API was successful.
        # it can be due to the non existance of data for the zone
        # (for instance : ocean area)
        if response.status_code != 200:
            print('non existing land location or problem with connection')
            return False
        
        else:

            url_content = response.content  # get the content from the link

            # Name of the file in the folder
            nameofthefile = ("dailydataforlat="
                             + str(lat)
                             + "andlong="
                             + str(long)
                             + "formonth"
                             + str(month)
                             + "forangle"
                             + str(angle)
                             + "forazimuth"
                             + str(azimuth)
                             + ".csv")
            # Create a csv file in the folder,
            # ready to be written with the retrieved data
            csv_file = open("../Climatic_data/"+nameofthefile, 'wb')

            # Overwrite the empty csv file with retrieved data
            csv_file.write(url_content)

            csv_file.close()
            
            print("File downloaded")

            return True

    else:
        return True








# Function that takes as inputs : a location, the azimuth of one of the frontal side, the month, and a PBR geometry.

# inputs : latitude (XX.XXX) or (X.XXX),longitude (XX.XXX) or (X.XXX),month [0:12],azimuth of frontal side [-180,+180],PBR geometry.
# Outputs : A list which contains a dataframe in first position [0] with :
# the upperbound,the lowerbound and the average for the averaged Solar power received for each hour of the day
# In second position [1], the temperatures for the day
# In Third position [2], the horizontal irradiance only for yield calculaiton
# This value is also averraged over a month because the climatic data is given for a average day for the whole month.
# Include Scheme for azimut explication


def Qreceived_bym2PBR_month(lat, 
                            long,
                            month,
                            azimuthfrontal,
                            height,
                            tubediameter,
                            gapbetweentubes,
                            horizontaldistance,
                            length_of_PBRunit,
                            width_of_PBR_unit):
    """Functin which estimates the solar power received by 1m2 of a given 
    PBR for each hour of a day.


    #Inputs:

        #lat : latitute expressed in format ; XX.XXX or X.XXX
        #long : longitude expressed in format ; XX.XXX or X.XXX
        #month : month of the year ; number of the month
        #azimuthfrontal : azimuth of the frontal side of the PBR unit ;  180:-180
        #height: Height of the PBR, ; m
        #tubediameter: Tube diameter  ; m
        #gapbetweentubes: Vertical gap between tubes ; m
        #horizontaldistance: Horizontal distance between stacks ; m
        #length_of_PBRunit: Length of the  PBR unit ; m
        #width_of_PBR_unit: width of the PBR unit ; m
        
        
        
    # Outputs :
                    
          
        #daily_received_solarpower : a dataframe with 3 columns :
            -the the upperbound,the lowerbound and the average for 
            Solar power received by 1m2 of PBR for each hour of the day ; W

        #temperatures : series of the air temperature for each hour of the day ; °C

        #Qhorizontal[1] : series of the ground irradiance on a horizontal surface for each hour ; W

        """
        
        

  
  # Azimuth calculation for the 4 sides of the PBR based on the 
  # input azimuth for the frontal side 1. See Azimuth scheme and appendix.

    if azimuthfrontal <= 0:
        azimuthfrontal2 = 180 + azimuthfrontal
        if azimuthfrontal > -90:
            azimuthlateral = azimuthfrontal - 90
            azimuthlateral2 = azimuthfrontal + 90
        elif azimuthfrontal <= - 90:
            azimuthlateral = 270 + azimuthfrontal
            azimuthlateral2 = azimuthfrontal + 90

    elif azimuthfrontal > 0:
        azimuthfrontal2 = azimuthfrontal - 180
        if azimuthfrontal < 90:
            azimuthlateral = azimuthfrontal - 90
            azimuthlateral2 = azimuthfrontal + 90

        elif azimuthfrontal >= 90:
            azimuthlateral = azimuthfrontal-90
            azimuthlateral2 = azimuthfrontal-270

    # Collect data for the 4 sides of the PBR + horizontal one.

    # First we download the files  for each side :

    downloadclimaticdata(lat, long, month, 90, azimuthfrontal)
    downloadclimaticdata(lat, long, month, 90, azimuthfrontal2)
    downloadclimaticdata(lat, long, month, 90, azimuthlateral2)
    downloadclimaticdata(lat, long, month, 90, azimuthlateral)
    downloadclimaticdata(lat, long, month, 0, azimuthlateral)

    # Then, we open the csvs that we just downloaded and create dataframes for
    # each side of the PBR
    
    # Qfrontalside = Solar power received by 1m2 of solid surface oriented 
    # as the frontal side of the PBR. 24 values for 24 hours
    #print('oka')


    #print("THE LAT ISSSSS", lat)
    #print("THE long ISSSSS", long)

    Qfrontalside = pd.read_csv("../Climatic_data/dailydataforlat="
                               + str(lat)
                               + "andlong="
                               + str(long)
                               + "formonth"
                               + str(month)
                               + "forangle"
                               + str(90)
                               + "forazimuth"
                               + str(azimuthfrontal)
                               + ".csv", sep="\t", header=None, encoding='unicode_escape', engine='python')
    #print('okb')
    Qfrontalside2 = pd.read_csv("../Climatic_data/dailydataforlat="
                                + str(lat)
                                + "andlong="
                                + str(long)
                                + "formonth"
                                + str(month)
                                + "forangle"
                                + str(90)
                                + "forazimuth"
                                + str(azimuthfrontal2)
                                + ".csv", sep="\t", header=None, encoding='unicode_escape', engine='python')
    #print('okc')
    Qlateralside = pd.read_csv("../Climatic_data/dailydataforlat="
                               + str(lat)
                               + "andlong="
                               + str(long)
                               + "formonth"
                               + str(month)
                               + "forangle"
                               + str(90)
                               + "forazimuth"
                               + str(azimuthlateral)
                               + ".csv", sep="\t", header=None, encoding='unicode_escape', engine='python')
    #print('okd')
    Qlateralside2 = pd.read_csv("../Climatic_data/dailydataforlat="
                                + str(lat)
                                + "andlong="
                                + str(long)
                                + "formonth"
                                + str(month)
                                + "forangle"
                                + str(90)
                                + "forazimuth"
                                + str(azimuthlateral2)
                                + ".csv", sep="\t", header=None, encoding='unicode_escape', engine='python')
    #print('oke')
    Qhorizontal = pd.read_csv("../Climatic_data/dailydataforlat="
                              + str(lat)
                              + "andlong="
                              + str(long)
                              + "formonth"
                              + str(month)
                              + "forangle"
                              + str(0)
                              + "forazimuth"
                              + str(azimuthlateral)
                              + ".csv", sep="\t", header=None, encoding='unicode_escape', engine='python')  # Note that azimuth is unisignificant for a horizontal surface

    # Collect temperatures
    # The air temperature in the csvs is same regardless of the side of the PBR.
    temperatures = Qhorizontal[4]

    # Collect actual surfaces for each side, as a function of geometry

    frontal_side_surface = functions.PBR_geometry(height, tubediameter,
                                                  gapbetweentubes,
                                                  horizontaldistance,
                                                  length_of_PBRunit,
                                                  width_of_PBR_unit)[2]
    
    lateral_side_surface = functions.PBR_geometry(height, tubediameter,
                                                  gapbetweentubes,
                                                  horizontaldistance,
                                                  length_of_PBRunit,
                                                  width_of_PBR_unit)[3]

    horizontal_surface = functions.PBR_geometry(height, tubediameter,
                                                gapbetweentubes,
                                                horizontaldistance,
                                                length_of_PBRunit,
                                                width_of_PBR_unit)[4]
 
    # Calculate upper and lower bounds for each time of the day

    # Will contain the upper bound, the lower bound,
    # and the average for the solar power received by the PBR for the 24 hours.
    # 1 row = 1 hour

    daily_received_solarpower = pd.DataFrame(np.zeros((24, 3)),
                                             columns=['Upper', 'Lower', 'Average'])

    for time in range(0, 24):  # each hour of the day

        # Upper bound
        
        Qreceived_by_frontal_side_of_the_parallelepiped = (Qfrontalside[1][time] 
                                                           * height
                                                           * length_of_PBRunit)

        Qreceived_by_frontal_side2_of_the_parallelepiped = (Qfrontalside2[1][time]
                                                            * height
                                                            * length_of_PBRunit)

        Qreceived_by_lateral_side_of_the_parallelepiped = (Qlateralside[1][time] 
                                                           * height
                                                           * width_of_PBR_unit)

        Qreceived_by_lateral_side2_of_the_parallelepiped = (Qlateralside2[1][time]
                                                            * height
                                                            * width_of_PBR_unit)

        Qreceived_by_horizontalsurface_of_the_parallelepiped = (Qhorizontal[1][time]
                                                                * width_of_PBR_unit
                                                                * length_of_PBRunit)

        upperbound_totalPBRunit = (Qreceived_by_frontal_side_of_the_parallelepiped
                                   + Qreceived_by_frontal_side2_of_the_parallelepiped
                                   + Qreceived_by_lateral_side_of_the_parallelepiped
                                   + Qreceived_by_lateral_side2_of_the_parallelepiped
                                   + Qreceived_by_horizontalsurface_of_the_parallelepiped)
       
 
        upperbound_per_squaremeter = upperbound_totalPBRunit/(width_of_PBR_unit*length_of_PBRunit)

        # Lower bound
        # The energy received by an actual side of the PBR is :   
        # Energy received by the side of the paralleleliped *(Actual surface / surface of the paralleleliped's side)
        
        Qreceived_by_frontal_side_of_the_actualPBR = (Qreceived_by_frontal_side_of_the_parallelepiped
                                                      * (frontal_side_surface/(height*length_of_PBRunit)))

        Qreceived_by_frontal_side2_of_the_actualPBR = (Qreceived_by_frontal_side2_of_the_parallelepiped
                                                       * (frontal_side_surface/(height*length_of_PBRunit)))

        Qreceived_by_lateral_side_of_the_actualPBR = (Qreceived_by_lateral_side_of_the_parallelepiped
                                                      * (lateral_side_surface/(height*width_of_PBR_unit)))

        Qreceived_by_lateral_side2_of_the_actualPBR = (Qreceived_by_lateral_side2_of_the_parallelepiped
                                                       * (lateral_side_surface/(height*width_of_PBR_unit)))

        Qreceived_by_horizontalsurface_of_the_actualPBR =( Qreceived_by_horizontalsurface_of_the_parallelepiped
                                                          *(horizontal_surface/(length_of_PBRunit*width_of_PBR_unit)))

        lowerbound_totalPBRunit = (Qreceived_by_frontal_side_of_the_actualPBR 
                                   + Qreceived_by_frontal_side2_of_the_actualPBR 
                                   + Qreceived_by_lateral_side_of_the_actualPBR
                                   + Qreceived_by_lateral_side2_of_the_actualPBR
                                   + Qreceived_by_horizontalsurface_of_the_actualPBR)


        lowerbound_per_squaremeter = lowerbound_totalPBRunit /(width_of_PBR_unit*length_of_PBRunit)

        # Updating the table with the calculate value for the correspond time

        daily_received_solarpower['Upper'][time] = upperbound_per_squaremeter

        daily_received_solarpower['Lower'][time] = lowerbound_per_squaremeter

        daily_received_solarpower['Average'][time] = (
            upperbound_per_squaremeter+lowerbound_per_squaremeter)/2

    return [daily_received_solarpower, temperatures, Qhorizontal[1]]





