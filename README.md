# Code documentation for the LCA_deep_uncertainty_fish_micro repository

## Overview of the repository

This repository contains the code used to reproduce the results of the manuscript: *Jouannais.P,Blanco.CF, Pizzol.M, XXX (under review)* 

This is not a generalized platform/package to apply the procedure to any case, but the scripts sepcific to the procedure and not to the case (microalgae and fish production) can be further adapted.

**Cite this repository:**
[![DOI]1)



 
### Overview of folders and files:

**Data**

+ **_elemental_contents.csv_** includes average elemental compositions (N, P, C) of biochemical classes  (Lipid, Phospholipids, proteins, Carbohydrates) 


+ **_Ingredient_composition_withNP.csv_** includes the fish feed composition (wheat, oil etc.) and the calculated composition in term of biochemical classes (Lipid, proteins, carbohydrates, ash, water)+ NP

+ **_Technosphere_Fish_micro_17_03.csv_** is the foreground technosphere matrix  and used for the calculations.

+ **_AH_combi_1.json_**  foreground database importable to Brightway2 (1).


+ **_Micro_for_combi.json_**  foreground database importable to Brightway2 (2).



**Environment**

+ **env_bw_windows.yml** File needed to create the virtual environment on WINDOWS.
+ **env_bw_ubuntu_full.yml** File needed to create the virtual environment on UBUNTU.


**Outputs**

The raw outputs (table with data points without regionalization) of csv and xlsx types are saved in this folder.

**Background_mc**

Folder in which the MonteCarlo iterations for the background are saved as pickle objects.

**Map**

Contains the geographic shapefiles for Europe.

**Climatic_data**

Will contain the climatic data files downoladed from the PVGIS API when the simulations are performed.

**Climatic_data**

Will contain the climatic data files downloaded from the PVGIS API when the simulations are launched.

**geo objects**

Will contain the random grid of locations over Europe for microalgal production as a geopandas pickle object when the simulations are performed.

**Modif ema**

Contains the scripts that must be replaced in the original ema_workbench package.

**borgonovo indices**

Contains the matlab scripts to calculate Borgonovo's GSA deltas.
Will contain the deltas once calculated.

**PRIM_process**

Will contain the results of the regionalization and of PRIM application. Divided into 6 sub-folders containing different types of outputs: plots of the boxes, numerical descriptions of the boxes, peeling trajectories etc.

**Scripts** and **Scripts accessory**

+ fourteen **.py** files: python scripts including the stochastic LCA models itself and the scripts of the procedure. 



Files, scripts, and their functions'interconnections are mapped below.  
<br>  

<img src="Code map_Fish_micro.jpg"
     alt="Markdown Monster icon"
     style="float: left; margin-right: 10px;" />  
<br>  




**Functions_for_physical_and_biological_calculations_3nd**

Contains functions to calculate values related to :

+ Strain biology and metabolism
+ PBR geometry
+ Culture physical properties and centrifugation
+ Optimization for fish feed substitution
+ Anaerobic digestion of microalgae
+ Anaerobic digestion and valorization of dead fish and sludge



**Cultivation_simul_2nd**

Contains functions which simulate the microalgal cultivation with solving of differential equations for temperature evolution and associated thermoregulation.


**Map_3rd**

Contains the functions generating and handling the geographic files/locations grid/geodataframes.

**Retrieving_solar_and_climatic_data_2nd**

Contains functions which download solar and temperature data from the European Photovoltaic Geographical System and estimate the solar power received by a given PBR geometry during the day.


**Technosphere_matrix_modifications_fish_3rd**

Contains functions which modify foreground technospheres  with new mortalities, FCR, excretion etc. 

**Main_functions_Fish_micro**

Script containing the main functions performing the simulations to produce the raw output tables.


**Borgonovo_indices_to_du_function_modular**

Contains functions assigning relative resolutions to parameters/indicators according to their GSA's indices.

**Simulation_fish_micro_104830**

Script which calls the functions to run simulations, LCAs and produce the raw output tables.


 See *Reproducing results from the article.*



**Prepare_project_4th** 

Creates the Brightway2 project and loads the foreground databases in it. Imports your local version of ecoinvent 3.8 consequential in the new project and loads bioshpere3.

See *Reproducing results from the article.*


**Processing_Multiindex_regions_PRIM**

Script which performs regionalization and apply PRIM to discover boxes and final results.


 See *Reproducing results from the article.*


**Parallel_ecoinvent_mc_168000**

Script which performs the MC iterations for the background (1 unit of each technosphere input)


 See *Reproducing results from the article.*


**Combine_background_mc**

Accessory script to combine different outputs of background MC iterations, if different chunks have been generated on different instances.


 See *Reproducing results from the article.*



**Combine_and_clean_results**

Script to clean raw outputs and combine chunks of the raw output, if different chunks have been generated on different instances.

 See *Reproducing results from the article.*
 
 
**Add_bioact_molec_fix_function**

Function which calculates an additional indicator a posteriori on the raw output.





<br>

### Reproducing results of the article

  &#128680;&#128680;&#128680;**WARNING 
Reproducing the results of the article requires a very large amount of data points (500 000 in the article) which is only reasonably achievable in a few days if using multiple instances with multiple cores and recombining results. Most of the functions are written to be called in parallel with the package "ray".

The following instructions are given without assuming any recombination of results from multiple instances.

The regionalization and PRIM application also require substantial computing capacities (large memory/multiple CPUS).


*Requirements*

+ Miniconda or Anaconda
https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html

+ A local version of the ecoinvent 3.8 consequential database

+ A python interface (e.g., Spyder) and a matlab licence


*Step by step procedure:*

1. **Download or clone a local copy of the repository. Keep the folders structure.**

2. **Prepare a conda environment with all needed packages**

+ From terminal or Anaconda/Miniconda terminal access the "Environment" folder. The path depends on where you saved the repository:

```
 cd <yourpathtothefolder/Environment>
```

+ Create the conda environment with all necessary packages using the .yml file corresponding to your OS.

**For Windows:**


```
conda env create --file env_bw_windows.yml
```

+ Activate the newly created environment:

```
conda activate env_bw_windows
```

+ Install the ray package from pip:
```
pip install ray
```

**For Ubuntu:**


```
conda env create --file env_bw_ubuntu.yml
```

+ Activate the newly created environment:

```
conda activate env_bw_ubuntu
```

+ Install the ray package from pip:
```
pip install ray
```

For MACOS, you should try to install the environment from the ubuntu file and, in case of issues, complete the environment by installing the problematic packages manually. 




3. **Set up the Brigtway2 project**

+ In the Scripts directory, open the file **Prepare_project_4th.py** in a text editor or python interface and change the value of the variable ```ei38dir``` by specifying the directory where the ecoinvent files are on your drive: ```ei38dir = <yourpathtoecoinventfiles>```. 

+ From the python interface or from command, execute the whole script to prepare the Brightway2 project (```python Prepare_project_4th.py```) .

4. **Run the simulations** 

4.1 ***Perform MC iterations in the background***

+ In the Scripts directory, open the file **Parallel_ecoinvent_mc_168000**. 

+ The script is parameterized for 168000 iterations. Change the number of simulations if needed and run the simulation by executing the whole script from the python interface or from command line .This will calculate stochastic LCA impacts for 1 unit of each of the 90+ foreground technosphere inputs, for all methods and save results in the folder *Background_mc*. Takes around 2 or 3 days for instances with 64 CPUs.



4.2 ***Calculate LCAs and generate the raw output***

+ In the Scripts directory, open the file **Simulation_fish_micro_104830**. 

+ Change the path to the correspond file of background MC results generated in 4.1. 

+  The script is parameterized for 104830 iterations. Change the number of simulations to match the size of the background MC file. Run the simulation by executing the whole script from the python interface or from command line. 

+ Wait for all the simulations to be finished. The script will export the raw output tables with all data points in "Outputs". Takes around 4 days for instances with 64 CPUs.


4.3 ***Clean (and combine if necessary) raw outputs***

+ In the "Scripts accessory"" directory, open the file **Combine_and_clean_results**. 

+ Change the path to the corresponding raw output files generated in 4.2. 

+ Execute the whole script from the python interface or from command line. 

+ The script will export the clean/combined output tables with all data points in the subfolder "Outputs_recombined"" within "Outputs".


4.4 ***Calculate borgonovo deltas on Matlab***

+ In Matlab, use "Calculate_borgnovo" and the recombined files saved for each method in 4.3 to calculate and export the Borgonovo's deltas for each method/impact category.

+ Export results in the containing folder "borgonovo_indices".




4.5 ***Modify the ema_workbench package to add transparency to the regions***

+ In the conda installation of the ema_workbench package https://emaworkbench.readthedocs.io/en/latest/ (J.H Kwakkel), replace files with their modified counterparts present in the folder "Modif ema". This will add transparency to the dots to see  stacked regions.



4.6 **Apply the regionalization and PRIM to discover boxes**
  
+ In the Scripts directory, open the file **Processing_Multiindex_regions_PRIM**.

+ Replace the paths to the borgonovo indices generated in 4.4.

+ Execute the script. &#128680;&#128680;&#128680;**WARNING . Here again, computing requirements are large and the procedure should be run on instances with large memory and multiple CPUs.

+ This will save box pairplots, boxes numerical descriptions, information about the regionalization step, peeling trajectories into the corresponding subfolders in "PRIM_process".

<br>  

