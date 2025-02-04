#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 14:05:53 2024

@author: connor

Plot initial and final pH from batch experiments
"""

path_save = '/Users/connor/Documents/Projects/Lactate growth (Roy and Cam)/Paper\
/Plots/'

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

### pH dependency experiment November 2023 ###
path_read = '/Users/connor/Documents/Projects/Lactate growth (Roy and Cam)/Paper\
/pH shifts.xlsx'

data = pd.read_excel(path_read, index_col = 0)
print(data)
data = data.set_axis([str(column) for column in data.columns.to_list()], axis = 1)
pH = data[[column for column in data.columns.to_list() if 'std' not in column]]
std = data[[column for column in data.columns.to_list() if 'std' in column]]
print(pH.columns.to_list())

separation = 0.1
width = 0.2

fig, ax = plt.subplots()

counter = 0 

xs = np.arange(len(pH.columns.to_list()))

for group in data.index.to_list():
    
    buffer = (counter - 1)*(width / 2 + 2*separation)
    x = [i + buffer for i in xs]
    y = pH.loc[group]
    yerr = std.loc[group]
    print(y)
    ax.bar(x, y, width = width, yerr = yerr, label = group)
    counter += 1
    
ax.legend()
ax.set_xticks(xs, pH.columns.to_list())
ax.set_ylim(5, 8)
ax.set_ylabel('Final pH')
ax.set_xlabel('Initial pH')
ax.set_title('pH Changes')

fig.savefig(path_save+'pH Change.eps', dpi = 300)

### Dual electron donor experiment August 2024 ###

path_read = '/Users/connor/Documents/Projects/Lactate growth (Roy and Cam)/Paper\
/pH shifts dual electron donor.xlsx'

data = pd.read_excel(path_read, sheet_name = 'Means', index_col = 0)

rows = [row for row in data.index.to_list() if 'neg' not in row]
neg_rows = [row for row in data.index.to_list() if 'neg' in row]

means = data.loc[rows]['Mean']
stds = data.loc[rows]['Std']

neg_means = data.loc[neg_rows]['Mean']
neg_stds = data.loc[neg_rows]['Std']

fig, ax = plt.subplots()

xs = np.arange(2)
sep = 0.2
width = 0.3
colours = ['r', 'b']

group = 'M. elsdenii'
i = 0

x = xs[i]
y = means[group]
y_neg = neg_means[group+' neg']
yerr = stds[group]
yerr_neg = neg_stds[group+' neg']
ax.bar(x - sep, y, yerr = yerr, width = width, color = colours[0])
ax.bar(x + sep, y_neg, yerr = yerr_neg, width = width, color = colours[1])

group = 'M. hexanoica'
i = 1

x = xs[i]
y = means[group]
y_neg = neg_means[group+' neg']
yerr = stds[group]
yerr_neg = neg_stds[group+' neg']
ax.bar(x - sep, y, yerr = yerr, width = width, color = colours[0], label = 'Final pH')
ax.bar(x + sep, y_neg, yerr = yerr_neg, width = width, color = colours[1], label = 'Negative')
    
ax.legend()
ax.set_xticks(xs, means.index.to_list())
ax.set_ylim(5, 7.5)
ax.set_ylabel('pH')
ax.set_title('pH Changes')
fig.savefig(path_save+'Dual Electron Donors pH Change.eps', dpi = 300)
    
