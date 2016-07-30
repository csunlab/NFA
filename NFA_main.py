'''
Network-based Frequency Analysis (NFA) - Main Script

Derrible S, Ahmad N (2015) Network-Based and Binless Frequency Analyses. PLoS ONE 10(11): e0142108. doi:10.1371/journal.pone.0142108
Available at: http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0142108

version 1.00
'''

#Libraries needed to run the tool

import numpy as np
import csv
import pandas as pd
from NFA_functions import *
import os

#Ask for file name
file_name = raw_input("Name of file:")
file_header = raw_input("File has labels and header (Y):")
sampling = raw_input("Sample data (rows) for faster runtime; enter proportion (e.g., 0.1) or leave blank to keep all:")

#Create a pandas dataframe from the csv file       
if file_header == 'Y':
    data_original = pd.read_csv(file_name + '.csv', header=0, index_col=0, quotechar='"')
else:
    data_original = pd.read_csv(file_name + '.csv', header=None, index_col=None)
    data_original.index = data_original.index.astype(str)

#Randomly selects a sub-sample of the data if desired and reports the number of columns and rows to analyze 
if sampling != '':
    data = data_original.loc[np.random.choice(data_original.index, int(float(sampling)*len(data_original.index)), replace=False)]
else:
    data = data_original

print("{0} rows and {1} columns to analyze".format(len(data.index), len(data.columns.values)))
print("")

#Record all question to save all plots and csv files for all networks generated. Doing this significantly slows down the process but may be useful to invistigate the results
record_all = raw_input("Record results for all networks generated for all columns (Y):")

#Define Zeta Type if not default
zeta_type = raw_input("Zeta coefficient is supplied (S), manual input (M), default (D):")
zeta_list = []

if zeta_type == 'S':
    zeta_file = raw_input("Name of zeta coefficients csv file:")
    zeta_name = pd.read_csv(zeta_file + '.csv', header=None)
    zeta_list = [float(zeta_name.ix[x,0]) for x in arange(len(zeta_name.index))]
elif zeta_type == 'M':
    zeta_list = [float(var) for var in raw_input("Enter coefficients (no comma between numbers):").split()]

zeta_thresh = raw_input("Threshold number to determine optimal zeta coefficient (default if left blank):")


#Collect optional information for graph aesthetics
option = raw_input("Enter optional information for graph scales and labels (Y):")

if option == 'Y':
    x_label = raw_input("Label for x axis:")
    coef_show = raw_input("Show zeta coefficient on figure (Y):")
    column_name = raw_input("Add column name as title (Y):")
    mode_name = raw_input("Show mode name (Y):")
    mode_value = raw_input("Show mode value (Y):")
    x_scale = raw_input("Log (L) scale for x axis (leave blank for normal):")
    y_scale = raw_input("Log (L) scale for y axis (leave blank for normal):")
    x_limits = raw_input("Set boundaries for x axis (Y):")
    if x_limits == 'Y':
        x_limits1 = raw_input("   Lower boundary:")
        x_limits2 = raw_input("   Higher boundary:")
    else:
        x_limits1 = ""
        x_limits2 = ""
    y_limits = raw_input("Set boundaries for y axis (Y):")
    if y_limits == 'Y':
        y_limits1 = raw_input("   Lower boundary:")
        y_limits2 = raw_input("   Higher boundary:")
    else:
        y_limits1 = ""
        y_limits2 = ""
    grid_lines = raw_input("Show gridlines (Y):")

else:
    x_label = ""
    coef_show = ""
    column_name = ""
    x_scale = ""
    y_scale = ""
    mode_name = ""
    mode_value = ""
    x_limits1 = ""
    x_limits2 = ""
    y_limits1 = ""
    y_limits2 = ""
    grid_lines = ""


#Define meta list and properties that will store optimal networks
results = []
results_header = []


#Start looping the columns; e.g., year 2000, then 2001, then 2002
for column in arange(len(data.columns.values)):
    current = data.columns.values[column]
    print("")
    print("Starting {0}".format(current))

    #Find none blank values and assign them in list distribution purely to get the zeta values in default 'D' mode
    if zeta_type != 'S' and zeta_type != 'M':
        from NFA_zeta import zeta_default
        zeta_list = zeta_default(pd.to_numeric(data.ix[:,column]).dropna())

    #Define list to store all the temporarily results for the current column
    temp_results = []

    #Compute, store, and save all networks for list of zeta coefficients
    for i in arange(len(zeta_list)):
        temp_results.append(network_calc(pd.to_numeric(data.ix[:,column]).dropna(), zeta_list[i], column)) 
        if record_all == 'Y':
            temp_zeta = [temp_results[i][1]]
            temp_zeta_current = int(current)
            save_individual_results(temp_results[i], temp_results[i][1], file_name, current)
            plot_final_network(temp_results[i], temp_results[i][1], file_name, current, x_label, coef_show, column_name, mode_value, x_scale, y_scale, x_limits1, x_limits2, y_limits1, y_limits2, grid_lines)

    #Determine "optimal" network
    giant_size = []
    for i in arange(len(temp_results)):
        giant_size.append(temp_results[i][9]) #Store size of giant clusters for all networks created


    list_track = []
    if zeta_thresh != '':
        step = int(zeta_thresh)
    else:
        from NFA_zeta import zeta_threshold
        step = zeta_threshold(temp_results) #See NFA_zeta script to modify default

    for i in arange(len(temp_results)):
        if i <= (len(temp_results) - step):
            track = (1/float(step))*sum(giant_size[i:i+step])/float(giant_size[i])
            list_track.append(track)
        else:
            list_track.append(1)

    #optiNet_id = len(temp_results)-1 #Set temporary optimal network as the last one in case none are found below
    count = 0
    running = 0
    while running == 0:
        if 0.99999 < list_track[count] < 1.00001: #Cannot have = 1 seems python sometimes store 1 as 0.999999999 or 1.00000001
            id_opti = count
            running = 1
        else:
            count += 1
    
    id_name = temp_results[id_opti][0].vs[temp_results[id_opti][7]][name]
    id_value = temp_results[id_opti][0].vs[temp_results[id_opti][7]]['value']

    save_prop(temp_results, file_name, current)
    plot_prop(temp_results, file_name, current, id_opti, temp_results[id_opti][1])
    plot_all_networks(temp_results, file_name, current, x_scale, y_scale, x_limits1, x_limits2, y_limits1, y_limits2, grid_lines)

    save_individual_results(temp_results[id_opti], temp_results[id_opti][1], file_name, current)

    print("Optimal found for coef {0}, with value {1} for {2}".format(temp_results[id_opti][1], id_value, id_name))

    results.append(float(temp_results[id_opti][1]))
    results_header.append(int(current))
    plot_final_network(temp_results[id_opti], temp_results[id_opti][1], file_name, current, x_label, coef_show, column_name, mode_value, x_scale, y_scale, x_limits1, x_limits2, y_limits1, y_limits2, grid_lines)

print("")
print("Done!")
os._exit(0)
