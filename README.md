# MSFD
In this repo are hosted codes to arrange datasets and perform indicators calculation.


## Script 1 

### Description:
 Starting from medits TA, TB and TC, the code calculate Biomass Index, L95 and some other useful information. Script 1 retain the possibility to calculate the L95 on the total population or on individuals larger than a threshold specified by the user.

### Usage
Inputs are medits_ta, medits_tb medits_tc, gsa_coordinates, Sp_Medits, ASFIS_2017. All files must be stored in the folder specified at line 44 (input_dir).
Directory of output must be specified also, at line 46 (output_dir). Size threshold setting are at line 52 (seize_thresh), here you can put what is the minimum size (in mm) on which to calculate L95. 
