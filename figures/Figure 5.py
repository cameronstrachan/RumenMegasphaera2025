#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 11:34:35 2024

@author: connor

Plotting substrate utilization and product formation from batch experiments with
glucose and lactate.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

data_path = '/Users/connor/Documents/Projects/Lactate growth (Roy and Cam)/\
Dual electron donor experiment/August 2024/Dual electron donor compiled data 0824.xlsx'
save_path = '/Users/connor/Documents/Projects/Lactate growth (Roy and Cam)/Paper/Plots/'

ME_data = pd.read_excel(data_path, sheet_name = 'M. elsdenii')
MH_data = pd.read_excel(data_path, sheet_name = 'M. hexanoica')
product_colours = sns.color_palette("viridis", 5)
substrate_colours = sns.color_palette("magma", 4)
linewidth = 1.5
elinewidth = 1
capsize = 1.5



all_products = ['Propionate',
                 'Butyrate',
                 'Pentanoate',
                 'Hexanoate',
                 'Octanoate']

all_substrates = ['Acetate(GC)',
                  'Lactate',
                  'Glucose',
                  'Acetate']

product_colour_dict = {compound : product_colours[i] for i, compound in enumerate(all_products)}
substrate_colour_dict = {compound : substrate_colours[i] for i, compound in enumerate(all_substrates)}

fig, axes = plt.subplots(2, 2, sharex = True, sharey = True, figsize = (9, 6), constrained_layout = True)

# M. elsdenii

### Add to this list when lactate and glucose data comes in ###

ME_products = [
             'Propionate', 
             'Butyrate', 
             'Pentanoate', 
             'Hexanoate', 
             ]

ME_substrates = [
             # 'Acetate (GC)', 
             'Lactate', 
             'Glucose', 
             'Acetate'
             ]

MH_products = [
             'Butyrate', 
             'Hexanoate', 
             'Octanoate'
             ]

MH_substrates = [
             'Lactate', 
             'Glucose', 
             'Acetate'
             ]

time = ME_data['Time (hrs)'].to_list()

ME_OD = ME_data['OD Mean']
ME_OD_std = ME_data['OD Std']
MH_OD = MH_data['OD Mean']
MH_OD_std = MH_data['OD Std']
max_OD = np.max(np.concatenate((np.array(ME_OD), np.array(MH_OD))))
max_conc = max(ME_data[[x+' Mean' for x in ME_products+ME_substrates]].max(axis = None),
               MH_data[[x+' Mean' for x in MH_products+MH_substrates]].max(axis = None))

sizer = max_conc / max_OD

axes[0, 0].errorbar(time,
               sizer*ME_OD,
               sizer*ME_OD_std,
               ecolor = 'black',
               color = 'black', 
               label = 'OD',
               linestyle = 'dashed', 
               elinewidth = elinewidth,
               linewidth = linewidth,
               capsize = capsize)

axes[1, 0].errorbar(time,
               sizer*ME_OD,
               sizer*ME_OD_std,
               ecolor = 'black',
               color = 'black', 
               label = 'OD',
               linestyle = 'dashed',
               elinewidth = elinewidth,
               linewidth = linewidth,
               capsize = capsize)

axes[0, 1].errorbar(time,
               sizer*MH_OD,
               sizer*MH_OD_std,
               ecolor = 'black',
               color = 'black', 
               label = 'OD',
               linestyle = 'dashed', 
               elinewidth = elinewidth,
               linewidth = linewidth,
               capsize = capsize)

axes[1, 1].errorbar(time,
               sizer*MH_OD,
               sizer*MH_OD_std,
               ecolor = 'black',
               color = 'black', 
               label = 'OD',
               linestyle = 'dashed', 
               elinewidth = elinewidth,
               linewidth = linewidth,
               capsize = capsize)

secax_y1 = axes[0, 1].secondary_yaxis('right', functions = (lambda x: x / sizer, 
                                                         lambda x: x*sizer))

secax_y2 = axes[1, 1].secondary_yaxis('right', functions = (lambda x: x / sizer, 
                                                         lambda x: x*sizer))

secax_y1.set_ylabel('OD$_{600}$')
secax_y2.set_ylabel('OD$_{600}$')

# Using HPLC measurements for acetate for the time being
for substrate in ME_substrates:
    conc = ME_data[substrate+' Mean']
    conc_std = ME_data[substrate+' Std']
    axes[0, 0].errorbar(time, conc, conc_std, label = substrate,
                   ecolor = substrate_colour_dict[substrate],
                   color = substrate_colour_dict[substrate],
                   elinewidth = elinewidth,
                   linewidth = linewidth,
                   capsize = capsize)

for product in ME_products:
    conc = ME_data[product+' Mean']
    conc_std = ME_data[product+' Std']
    axes[1, 0].errorbar(time, conc, conc_std, label = product,
                   ecolor = product_colour_dict[product],
                   color = product_colour_dict[product],
                   elinewidth = elinewidth,
                   linewidth = linewidth,
                   capsize = capsize)

# M. hexanoica

# Using HPLC measurements for acetate for the time being
for substrate in MH_substrates:
    conc = MH_data[substrate+' Mean']
    conc_std = MH_data[substrate+' Std']
    axes[0, 1].errorbar(time, conc, conc_std, label = substrate,
                   ecolor = substrate_colour_dict[substrate],
                   color = substrate_colour_dict[substrate],
                   elinewidth = elinewidth,
                   linewidth = linewidth,
                   capsize = capsize)

for product in MH_products:
    conc = MH_data[product+' Mean']
    conc_std = MH_data[product+' Std']
    axes[1, 1].errorbar(time, conc, conc_std, label = product,
                   ecolor = product_colour_dict[product],
                   color = product_colour_dict[product],
                   elinewidth = elinewidth,
                   linewidth = linewidth,
                   capsize = capsize)

axes[0, 0].grid(linestyle = '--')
axes[1, 0].grid(linestyle = '--')

axes[1, 0].set_xlabel('Time (hrs)')
axes[0, 0].set_ylabel('Substrates\n\nConcentration (mM)')
axes[1, 0].set_ylabel('Products\n\nConcentration (mM)')

axes[0, 1].grid(linestyle = '--')
axes[1, 1].grid(linestyle = '--')

axes[1, 1].set_xlabel('Time (hrs)')

axes[0, 0].set_title('M. elsdenii')
axes[0, 1].set_title('M. hexanoica')

lines = [] 
labels = [] 

for ax in fig.axes: 
    Line, Label = ax.get_legend_handles_labels() 
    repeats = []
    for i in range(len(Label)):
        if Label[i] in labels:
            repeats += [i]
    
    new_Label = [name for i, name in enumerate(Label) if i not in repeats]
    new_Line = [name for i, name in enumerate(Line) if i not in repeats]
         
    lines.extend(new_Line) 
    labels.extend(new_Label)

fig.legend(lines, labels, bbox_to_anchor=(1.18, 0.65)) 
fig.suptitle('Megasphaera Dual Electron\nDonor Experiment')
fig.savefig(save_path+'Dual Electron Donor Experiment Results.eps', dpi = 300, bbox_inches = 'tight')

