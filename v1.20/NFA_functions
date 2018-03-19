#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Network-based Frequency Analysis (NFA) - Functions Script

Derrible S, Ahmad N (2015) Network-Based and Binless Frequency Analyses. PLoS ONE 10(11): e0142108. doi:10.1371/journal.pone.0142108
Available at: http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0142108

version 1.2 - for Python 3 and pandas 0.20.3 and later
'''

#Libraries needed to run the tool
import numpy as np
from igraph import *
import matplotlib.pyplot as plt
from pylab import *
import csv
import pandas as pd


#Calculate network from individual columns
def network_calc(self, coef, col):
    size = len(self.index)    
    adjacency = np.zeros(size**2).reshape(size, size)
    for i in arange(size):
        temp_a =  self.iloc[i]-coef <= self
        temp_b = self.iloc[i]+coef >= self
        adjacency[i] = temp_a & temp_b
        
    g = Graph.Adjacency((adjacency > 0).tolist(), mode=ADJ_UNDIRECTED)
    g.simplify()

    for i in arange(g.vcount()):
        if i < g.vcount():
            g.vs[i][name] = self.index[i]
            g.vs[i]['value'] = self.iloc[i]

    if max(g.degree()) != 0:
        max_degree = max(g.degree())
    else:
        max_degree = 1
    
    list_max_degree = [g.vs[idx].index for idx, li in enumerate(g.degree()) if li == max_degree]
    index_max_degree = int(np.min(list_max_degree))

    giant = g.clusters(WEAK).giant()

    return [g, coef, g.vcount(), g.ecount(), g.diameter(), g.density(), g.average_path_length(), index_max_degree, giant.density(), float(giant.vcount())/g.vcount()]


#Save general network results for an entire column
def save_prop(self, file_name, current):
    prop_file = open(file_name + '/' + file_name + '_coef_prop_' + str(current) + '.csv', 'w')
    try:
        newf = csv.writer(prop_file)
        newf.writerow(['Coefficient', 'Nodes', 'Links', 'Diameter', 'Density', 'Average Path Length', 'Mode', 'Giant Cluster Density', '% Nodes in Giant'])
        for i in range(len(self)):
            newf.writerow([self[i][1], self[i][2], self[i][3], self[i][4], self[i][5], self[i][6], self[i][0].vs[self[i][7]]['value'] , self[i][8], self[i][9]])
    finally:
        prop_file.close()
        print("Properties saved in csv file for {0}".format(current))


#Save network properties and results for a given coefficient
def save_individual_results(self, coef, file_name, current):
    coef_data_file = open(file_name + '/' + file_name + '_coef_data_' + str(current) + '_' + str(coef) + '.csv', 'w')
    try:
        newf = csv.writer(coef_data_file)
        newf.writerow(['Coefficient', self[1]])
        newf.writerow(['Nodes', self[2]])
        newf.writerow(['Diameter', self[4]])
        newf.writerow(['Giant Cluster Density', self[8]])
        newf.writerow(['% Nodes in Giant Cluster', self[9]])
        newf.writerow(['Mode Value', self[0].vs[self[7]]['value']])
        newf.writerow(['Mode Name', self[0].vs[self[7]][name]])
        newf.writerow([])
        newf.writerow(['Name', 'Value', 'Degree'])
        for i in range(self[2]):
            newf.writerow([self[0].vs[i][name], self[0].vs[i]['value'], self[0].degree(i)])
    finally:
        coef_data_file.close()


#Plot density, diameter and average shortest path lengths of all networks
def plot_prop(self, file_name, current, id_opti, cutoff):
    fig = plt.figure()
    fig.set_size_inches(8,12)
    fig.subplots_adjust(hspace=.5)
    rc('font', family='serif')

    coef = []
    giant_density = []
    giant_node = []
    diameter = []
    avg_path = []
    hindex = []
    mode_value = []
    for i in range(len(self)):
        coef.append(self[i][1])
        giant_density.append(self[i][8])
        giant_node.append(self[i][9])
        diameter.append(self[i][4])
        avg_path.append(self[i][6])
        hindex.append(coef[i] / self[i][0].vs[self[i][7]]['value'])
        mode_value.append(self[i][0].vs[self[i][7]]['value'])

    colors = []
    for i in range(len(self)):
        if i == id_opti:
            colors.append(0.9)
        else:
            colors.append(0.5)
    
    figdensity = fig.add_subplot(3, 2, 1)
    title('Giant Cluster Density')
    plot(coef, giant_density, marker='o', color='0.7')
    plot(coef[id_opti], giant_density[id_opti], marker='o', color='red')

    figavg_path = fig.add_subplot(3, 2, 2)
    title('% Nodes in Giant')
    plot(coef, giant_node, marker='o', color='0.7')
    plot(coef[id_opti], giant_node[id_opti], marker='o', color='red')

    figavg_path = fig.add_subplot(3, 2, 3)
    title('Homogeneity Index')
    plot(coef, hindex, marker='o', color='0.7')
    plot(coef[id_opti], hindex[id_opti], marker='o', color='red')

    figavg_path = fig.add_subplot(3, 2, 4)
    title('Mode Value')
    plot(coef, mode_value, marker='o', color='0.7')
    plot(coef[id_opti], mode_value[id_opti], marker='o', color='red')

    figdiameter = fig.add_subplot(3, 2, 5)
    title('Diameter')
    plot(coef, diameter, marker='o', color='0.7')
    plot(coef[id_opti], diameter[id_opti], marker='o', color='red')
    
    figdiameter = fig.add_subplot(3, 2, 6)
    title('Average Path Length')
    plot(coef, avg_path, marker='o', color='0.7')
    plot(coef[id_opti], avg_path[id_opti], marker='o', color='red')

    savefig(file_name + '/' + file_name + '_FigProp_' + str(current) + '.pdf', dpi=150)
    close()


#Plot one page with networks for all zeta coefficients
def plot_all_networks(self, file_name, current, x_scale, y_scale, x_limits1, x_limits2, y_limits1, y_limits2, grid_lines):

    #To have enough spots on the figure for all sub-figures
    hori_length = int(len(self)**0.5)
    if mod(len(self), hori_length) == 0:
        vert_length = hori_length
    else:
        vert_length = int(hori_length)+1

    check = 0
    while check == 0:
        if hori_length * vert_length < len(self):
            hori_length += 1
        else:
            check = 1

    fig = plt.figure()
    fig.set_size_inches(hori_length*3, vert_length*4) #Each sub-figure will be on drawn on 3in x 4in
    fig.subplots_adjust(hspace=.5)
    rc('font', family='serif')
    
    for i in range(len(self)):
        subfig = fig.add_subplot(vert_length, hori_length, i+1)
        subfig.yaxis.set_visible(False)
        title('Coef: ' + str(self[i][1]))

        #Sort values from smallest to largest to plot graph as a line as opposed to an XY scatter
        to_sort = []
        for j in range(0, self[i][2]):
                       to_sort.append([self[i][0].vs[j]['value'], self[i][0].degree(j)])
        to_sort.sort(key=lambda x:x[0])
        x_val = [x[0] for x in to_sort]
        y_val = [x[1] for x in to_sort]

        plot(x_val, y_val, color='#483D8B')

        if x_scale == 'L':
            plt.xscale('log')
        if y_scale == 'L':
            plt.yscale('log')
        if x_limits1 != "" and x_limits2 != "":
            plt.xlim(x_limits1, x_limits2)
        if y_limits1 != "" and y_limits2 != "":
            plt.ylim(y_limits1, y_limits2)
        if grid_lines == 'Y':
            plt.grid(b=True, which='both', color='0.75', linestyle='--')

    savefig(file_name + '/' + file_name + '_FigAll_' + str(current) + '.pdf', dpi=150)
    close()


#Plot results from individual network for a given coefficient
def plot_final_network(self, coef, file_name, current, x_label, coef_show, column_name, mode_value, x_scale, y_scale, x_limits1, x_limits2, y_limits1, y_limits2, grid_lines):
    id_name = self[0].vs[self[7]][name]
    id_value = self[0].vs[self[7]]['value']

    #Sort values from smallest to largest to plot graph as a line as opposed to an XY scatter
    final_sort = []
    for j in range(self[2]):
        final_sort.append([self[0].vs[j]['value'], self[0].degree(j), self[0].vs[j][name]])

    final_sort.sort(key=lambda x:x[0])
    x_val = [x[0] for x in final_sort]
    y_val = [x[1] for x in final_sort]
    z_val = [x[2] for x in final_sort]

    #Define figure
    text_font1 = {'size':'10', 'color':'0.5', 'weight':'normal','horizontalalignment':'right'}
    text_font2 = {'size':'10', 'color':'0.5', 'weight':'normal','horizontalalignment':'left'}
    title_font = {'size':'12', 'color':'black', 'weight':'bold','horizontalalignment':'center'}
    fig = plt.figure()
    fig.set_size_inches(5,4.5)
    plot(x_val, y_val, color='#4682B4', linewidth=2)
    plt.ylabel('Connections')
    plt.tick_params(axis='both', length=0)

    #Add options based on user inputs
    if x_label != "":
        plt.xlabel(x_label)
    if coef_show == 'Y':
        fig.text(0.9, 0.92, '$\zeta=' + str(coef) + '$', **text_font1)
    if column_name == 'Y':
        title(current, **title_font)
    if mode_value == 'Y':
        fig.text(0.12, 0.92, 'Value {0:.2f} in {1}'.format(id_value, id_name), **text_font2)
    if x_scale == 'L':
        plt.xscale('log')
    if y_scale == 'L':
        plt.yscale('log')
    if x_limits1 != "" and x_limits2 != "":
        plt.xlim(x_limits1, x_limits2)
    if y_limits1 != "" and y_limits2 != "":
        plt.ylim(y_limits1, y_limits2)
    if grid_lines == 'Y':
        plt.grid(b=True, which='both', color='0.75', linestyle='--')
        
    savefig(file_name + '/' + file_name + '_FigFinal_' + str(current) + '_' + str(coef) + '.png', dpi=300)
    close()
