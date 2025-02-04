#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 16:22:50 2024

@author: connor

M. hexanoica grouped by pH, first 5 time points
"""
import os

os.chdir(os.path.dirname(__file__))

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Import growth data

path = "/Users/connor/Documents/Projects/Lactate growth (Roy and Cam)/pH time series/Results/20240130_Megasphaera different pH.xlsx"
path_save = '/Users/connor/Documents/Projects/Lactate growth (Roy and Cam)/Paper/Plots/'
OD_data = pd.read_excel(path, sheet_name = 'MH_growth_graph', index_col = [0,1])

rate_data = pd.read_excel(path, sheet_name = 'MH_rates_graph', index_col = [0,1])

# remove time pointss after 48.5

print(rate_data)
print(OD_data)
OD_remove = [56.50, 75.00, 77.50] # kept time point 50.50 in here
rate_remove = [75.0]
print(rate_data)
OD_data = OD_data.drop(OD_remove, axis = 0, level = 1)
rate_data = rate_data.drop(rate_remove, axis = 0, level = 1)
print(rate_data)

#%%

rate_data_slim = rate_data.drop(['∆Propionate SD', '∆Butyrate SD', '∆Pentanoate SD', '∆Caproate SD'], axis = 1)
rate_SDs = rate_data.drop(['∆Propionate', '∆Butyrate', '∆Pentanoate', '∆Caproate'], axis = 1)
group_list = ['No lactate', 'pH 5.5', 'pH 6.0', 'pH 6.5']

# intitialize x locations for bars and plot grouped bars (same for all groups)

x = rate_data.loc["No lactate"].index
width = np.max(x)/50  # the width of the bars

print(rate_data_slim[['∆Propionate']])
#%%
for i in range(len(group_list)):
    
    multiplier = 0
    name = group_list[i]
    print('plotting', name)
    
    fig, ax = plt.subplots()
    ax.set_ylim([-0.3, 1.2])
    ax.plot(OD_data.loc[name].index, OD_data.loc[name]['OD600'])
    ax.errorbar(OD_data.loc[name].index, OD_data.loc[name]['OD600'], 
                yerr = OD_data.loc[name]['OD600 STD'], fmt = 'none')
    
    # Resize the bars to plot them on right y axis
    
    sizer = np.max(rate_data_slim.loc[name]) / np.max(OD_data.loc[name])
    
    resized_rate_data_slim = rate_data_slim / sizer
    
    resized_rate_SDs = rate_SDs / sizer
    
    for attribute, measurement in resized_rate_data_slim.loc[name].items():
        offset = width * multiplier
        rects = ax.bar(x + offset, measurement, width, label=attribute[1:])
        ax.errorbar(x + offset, measurement, yerr = resized_rate_SDs.loc[name][attribute+' SD'], fmt = 'none')
        # ax.bar_label(rects, padding=3)
        multiplier += 1
    
    def convert(x):
        return x*sizer
    
    def revert(x):
        return x/sizer
    
    secax_y = ax.secondary_yaxis('right', functions = (convert, revert))
    secax_y.set_ylabel('Specific Production Rate (mM/OD/hr)')
    ax.set_xlabel('Time (hrs)')
    ax.set_ylabel('OD600')
    ax.set_title('M. hexanoica '+name)
    ax.legend()
    plt.savefig(path_save+'MH_'+name+'kinetics_products.eps', dpi = 300)
    plt.savefig(path_save+'MH_'+name+'kinetics_products.jpg', dpi = 300)
