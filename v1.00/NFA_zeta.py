'''
Network-based Frequency Analysis (NFA) - Zeta Script

Derrible S, Ahmad N (2015) Network-Based and Binless Frequency Analyses. PLoS ONE 10(11): e0142108. doi:10.1371/journal.pone.0142108
Available at: http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0142108

version 1.00 - for pandas versions earlier than 0.18.0
'''

#Libraries needed to run the tool
import numpy as np

#Determine the list of zeta coefficients
def zeta_default(distribution):
    zeta_list = []

    median = np.median(distribution)
    std = np.std(distribution)

    #Defining two cases 
    if float(std)/median >= 1:
        zeta_list = np.arange(0.1*median, 1.01*median, 0.1*median).round(2)
    else:
        zeta_list = np.arange(0.05*median, 0.51*median, 0.05*median).round(2)

    return zeta_list
    
#Set threshold for the giant cluster to be identical in size for x consecutive times, where x is 1/3 the total number of networks calculate for the column
def zeta_threshold(temp_results):
    return int(len(temp_results)/3) #Default: 3, change the value of '3' is desired
