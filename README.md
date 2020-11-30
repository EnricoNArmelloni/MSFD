# MSFD
In this repo are hosted codes to arrange datasets and perform indicators calculation.


## Script 1 

### Description:
 Starting from medits TA, TB and TC, the code calculate Biomass Index, L95 and some other useful information. Script 1 retain the possibility to calculate the L95 on the total population or on individuals larger than a threshold specified by the user.

### Usage
Inputs are medits_ta, medits_tb medits_tc, gsa_coordinates, Sp_Medits, ASFIS_2017. All files must be stored in the folder specified at line 44 (input_dir).
Directory of output must be specified also, at line 46 (output_dir). 
Size threshold setting are at line 52 (seize_thresh), here you can put what is the minimum size (in mm) on which to calculate L95. 
For use the script, please carefully read and set the parameters in the section "Settings", lines 42-70. Then, you can simply run the code.


## Script 2

### Description
This code serve to identify and fit a gam model, and allows to create standardized outputs for data exploration and diagnostic. It also permit to predict the timeseries using a standard value for sampling moth. The code must be run after the script 1, as it take as input data the outputs of script 1. 

### Usage
Input data is the L95 output file from script 1. The input directory should be the same as the output directory of code 1. Settings are at lines 10-18.
Line 16 serve to change the input file type, basing on the setting choosen for code 1. Alternatives are "total", for data withouth any kind of selection, "sizethreshold" for L95 calculated on individuals larger than the threshold and "selection" (not to be used for now). For example, if you run code 1 with the size threshold active, you have to write here "sizethreshold" so the appropriate data will be imported. 
To use the script, read and set the parameters in the section "Settings" and the you have to go line by line to choose the best model for your species and area.


## Script 3

### Description
Starting from code 1 output, the code performs a segmented regression and a Spearman rank correlation analysis. 

### Usage
Input data is the L95 output file from script 1. The input directory should be the same as the output directory of code 1. Settings are at lines 10-25.
Line 17 serve to change the input file type, basing on the setting choosen for code 1. Alternatives are "total", for data withouth any kind of selection, "sizethreshold" for L95 calculated on individuals larger than the threshold and "selection" (not to be used for now). For example, if you run code 1 with the size threshold active, you have to write here "sizethreshold" so the appropriate data will be imported. 
Line 23 serve to inform the code wheter you are using "raw" data or data coming from the gam prediction. If you type "Y", then you could write the month used to predict the dataset on line 25.
For use the script, please carefully read and set the parameters in the section "Settings", lines 10-25. Then, you can simply run the code.
