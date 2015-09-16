'''
Citation: Derrible, S., and Ahmad, N., "Network-Based and Binless Frequency Analyses", under review
(citation will be updated once the work will be published)
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
