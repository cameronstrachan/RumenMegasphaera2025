#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 11:03:59 2024

@author: connor

July 19th LHR Megasphaera high throughput kinetics
Trying to get growth rates directly from maximum tangent
Getting growth rate from each technical replicate individually
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from math import log, exp, sqrt
from scipy.optimize import curve_fit

#%% User defined functions

def transform(x, OD_min):
    if x != 0.0:
        return log(x / OD_min)

def high_throughput_graphs_replicates2(sheet_name, title, colours, cutoff, drop):
    
    """ Graphing curves of individual replicates """
    
    data = pd.read_excel(path, sheet_name = sheet_name, usecols = lambda x: x != 'Time',
                         index_col = 0)
    
    negatives_list = [c for c in data.columns if data[c]['Conc. (mM)'] == 'Neg']
    negatives = data[negatives_list]

    data_list = [c for c in data.columns if data[c]['Conc. (mM)'] != 'Neg'
                 and data[c]['Conc. (mM)'] != 250 and data[c]['Conc. (mM)'] != 300
                 and data[c]['Conc. (mM)'] != 200]
    data = data[data_list].drop(index = 0)
    
    data = data[data.loc['Conc. (mM)'].sort_values().index.to_list()]
    concentrations = np.sort(list(set(data.loc['Conc. (mM)'].to_list())))
    
    mu_dict = {}

    plt.rcParams.update({'font.size': 5})
    
    fig_OD, axes_OD = plt.subplots(4, 18, figsize = (12, 3), sharex = False, 
                                   sharey = True, tight_layout = True, dpi = 300)
    fig_der, axes_der = plt.subplots(4, 18, figsize = (12, 3), sharex = False, 
                                   sharey = True, tight_layout = True, dpi = 300)
    
    time = data.drop(['Conc. (mM)'], axis = 0).index
    max_log_OD = 0
    
    for i, conc in enumerate(concentrations):
        
        replicate = 1
        mu_list = []
        
        for col in data.columns.to_list():
            
            if data[col].loc['Conc. (mM)'] == conc and replicate < 5:
                
                colour = colours[replicate - 1]
                
                ax_OD = axes_OD[replicate - 1, i]
                ax_der = axes_der[replicate - 1, i]
                
                OD = np.array(data.drop(['Conc. (mM)'], axis = 0)[col].to_list())
                min_OD = OD.min()
                max_ind = np.argmax(OD) + 1
                min_ind = np.argmin(OD)
                
                # Need to mask out all values after the max and all values 
                # before the min
                mask = np.ones(len(OD), dtype = bool)
                mask[max_ind:] = False
                mask[:min_ind] = False
                primetime = time[mask].to_numpy()
                OD = OD[mask]
                
                log_OD = [transform(x, min_OD) for x in OD]
                
                if np.array(log_OD).max() > max_log_OD:
                    max_log_OD = np.array(log_OD).max()
                
                dxdt = [(log_OD[i+1] - log_OD[i]) / (primetime[i+1] - primetime[i]) 
                        for i in range(len(log_OD) - 1)]
                
                mu_mask = np.ones(len(dxdt), dtype = bool)
                mu_mask[:cutoff] = False
                max_mu_ind = np.argmax(np.array(dxdt)[mu_mask]) + int(cutoff)
                mu_max = np.array(dxdt)[max_mu_ind]
                
                ax_OD.plot(primetime, 
                            log_OD,
                            alpha = 0.8,
                            color=colour,
                            linewidth= 0.7)
                
                max_mu_int = log_OD[max_mu_ind] - dxdt[max_mu_ind]*primetime[max_mu_ind]
                
                simtime = np.linspace(0, primetime.max(), 10)
                tangent = [dxdt[max_mu_ind] * t + max_mu_int for t in simtime]
                ax_OD.plot(simtime, tangent, linewidth = 0.5, color = 'black')
                
                ax_OD.yaxis.set_ticks(np.arange(0, max_log_OD, 0.5))
                x_start, x_end = ax_OD.get_xlim()
                ax_OD.xaxis.set_ticks(np.arange(0, x_end, 10))
                ax_OD.set_ylim(0, max_log_OD + 0.2)
                
                ax_der.plot(np.delete(primetime, -1), 
                            dxdt,
                            alpha = 0.8,
                            color=colour,
                            linewidth=0.7)
                
                y_start, y_end = ax_der.get_ylim()
                x_start, x_end = ax_der.get_xlim()
                ax_der.xaxis.set_ticks(np.arange(0, x_end, 10))
                
                mu_list = mu_list + [mu_max]
                
                replicate += 1
            
        
        if conc not in drop:
            mu_mean = np.mean(np.array(mu_list))
            mu_std = np.std(np.array(mu_list))
            mu_dict[conc] = (mu_mean, mu_std)
            axes_OD[replicate - 5, i].set_title(f'{conc} mM\n{round(mu_mean, 3)} 1/h')
            axes_der[replicate - 5, i].set_title(f'{conc} mM\n{round(mu_mean, 3)} 1/h')
            
        else:
            mu_dict[conc] = (None, None)
            axes_OD[replicate - 5, i].set_title(f'{conc} mM\nDropped')
            axes_der[replicate - 5, i].set_title(f'{conc} mM\nDropped')
        
    
    handles, labels = axes_OD[0, 0].get_legend_handles_labels()
    plt.figlegend(handles, labels , bbox_to_anchor=[6.3/8, 0.7/5], 
               loc='center')
    
    plt.rcParams.update({'font.size': 8})
    fig_OD.suptitle(title)
    fig_der.suptitle(title)
    fig_OD.text(0.5, -0.01, 'Time (h)', ha='center')
    fig_OD.text(-0.01, 0.5, 'log(OD/OD min)', va='center', rotation='vertical')
    fig_der.text(0.5, -0.01, 'Time (h)', ha='center')
    fig_der.text(-0.01, 0.5, '(1/OD)dOD/dt', va='center', rotation='vertical')
    
    return fig_OD, fig_der, negatives, time, mu_dict     
        
def growth_rates_graph(dict_list, title, new, ax, colour, label):
    
    mu_dicts = [mu_dict.copy() for mu_dict in dict_list]
    
    growth_rates = np.array([[i[0] for i in list(mu_dicts[j].values())] 
                              for j in range(len(mu_dicts))]) 
    
    errors = np.array([[i[1] for i in list(mu_dicts[j].values())] 
                              for j in range(len(mu_dicts))])
    
    def weighted_mean(growth_rates, errors):
        growth_rates = [x for x in growth_rates if x is not None]
        errors = [x for x in errors if x is not None]
        numerator = 0
        denominator = 0
        for i in range(len(growth_rates)):
            numerator += growth_rates[i] / errors[i]**2
            denominator += 1 / errors[i]**2
        
        return numerator / denominator    

    def weighted_error(errors):
        errors = [x for x in errors if x is not None]
        denominator = 0
        for i in range(len(errors)):
            denominator += 1 / errors[i]**2
        
        return sqrt(1 / denominator)
        
    means = [weighted_mean(growth_rates[:, i], errors[:, i]) 
             for i in range(len(list(mu_dicts[0])))]
    
    weighted_errors = [weighted_error(errors[:, i]) 
             for i in range(len(list(mu_dicts[0])))]
    
    concentrations = list(mu_dicts[0].keys())
    
    if new:
        fig, ax = plt.subplots()
        ax.errorbar(concentrations, means, weighted_errors, fmt='.',
        markersize = 3,
        color=colour,
        ecolor= colour, 
        elinewidth=1, 
        capsize=2)
        ax.set_title(title+' Specific Growth Rates')
        ax.set_ylabel('Specific Growth Rate (1/h)')
        ax.set_xlabel('Lactate Conc. (mM)')
        ax.set_xticks(np.arange(0, np.array(concentrations).max()+10, 10))
        plt.grid(True)
        # ax.set_yscale('log')
    
        return fig, ax
    
    else:
        ax.errorbar(concentrations, means, weighted_errors, fmt='.',
        markersize = 3,
        color=colour,
        ecolor= colour, 
        elinewidth=1, 
        capsize=2,
        label = label)
        ax.set_ylabel('Specific Growth Rate (1/h)')
        ax.set_xlabel('Lactate Conc. (mM)')
        ax.set_xticks(np.arange(0, np.array(concentrations).max()+10, 10))
        plt.grid(True)

def growth_rates_graph_bio_rep(dict_list, strain, title, new, ax, colours):
    
    ### Subtracting growth rate at 0 mM
    
    mu_dicts = [mu_dict.copy() for mu_dict in dict_list]
    markerstyles = ['.', 'v', 's']

    if new:

        fig, ax = plt.subplots()
        
        for i, dic in enumerate(mu_dicts):
            
            growth_rates = np.array([i[0] for i in list(dic.values())]) 
         
            errors = np.array([i[1] for i in list(dic.values())])
            
            conc_mask = np.ones(len(growth_rates), dtype = bool)
            
            none_inds = [i for i, j in enumerate(growth_rates) if j == None]
            
            growth_rates = np.array([i for i in growth_rates if i is not None]) 
         
            errors = np.array([i for i in errors if i is not None])  
            
            for i in none_inds:
                conc_mask[i] = False
            
            concentrations = np.array(list(mu_dicts[0].keys()))[conc_mask]
        
            ax.errorbar(concentrations, growth_rates, errors, fmt= markerstyles[i],
            markersize = 3,
            color=colours[i],
            ecolor= colours[i], 
            elinewidth=2, 
            capsize=1, 
            label = strain+f' {i+1}'
            )
            ax.set_title(title+' Specific Growth Rates')
            ax.set_ylabel('mu (1/h)')
            ax.set_xlabel('Lactate Conc. (mM)')
            ax.set_xticks(np.arange(0, np.array(concentrations).max()+10, 10))
            plt.grid(True)
            
        
        ax.legend()
        return fig
    
    else:
        
        for i, dic in enumerate(mu_dicts):
            
            growth_rates = np.array([i[0] for i in list(dic.values())]) 
         
            errors = np.array([i[1] for i in list(dic.values())])
            
            conc_mask = np.ones(len(growth_rates), dtype = bool)
            
            none_inds = [i for i, j in enumerate(growth_rates) if j == None]
            
            growth_rates = np.array([i for i in growth_rates if i is not None]) 
         
            errors = np.array([i for i in errors if i is not None])  
            
            for i in none_inds:
                conc_mask[i] = False
            
            concentrations = np.array(list(mu_dicts[0].keys()))[conc_mask]
        
            ax.errorbar(concentrations, growth_rates, errors, fmt=markerstyles[i],
            markersize = 3,
            color=colours[i],
            ecolor= colours[i], 
            elinewidth=2, 
            capsize=1, 
            label = strain+f' {i+1}', 
            )
            ax.set_title(title+' Specific Growth Rates')
            ax.set_ylabel('mu (1/h)')
            ax.set_xlabel('Lactate Conc. (mM)')
            ax.set_xticks(np.arange(0, np.array(concentrations).max()+10, 10))
            plt.grid(True)
        
path = 'data/chemical_growth_kinetic/Megapshaera HTK July 19 2024 processed data.xlsx'

#path_save = '/Users/connor/Documents/Projects/Lactate growth\
# (Roy and Cam)/High Throughput Kinetics/LHR Megapshaera HTK July 19th 2024/'
 
colours = ['r', 'b', 'g', 'y']
 
#%% M. elsdenii 1
colours = ['r', 'b', 'g', 'y']
title = 'M. elsdenii 1 0719'
sheet_name = 'ME 1'

fig_OD, fig_der, negatives, time, me_mu_dict_1 = high_throughput_graphs_replicates2(
    sheet_name, title, cutoff = 0, colours = colours, drop = [])
#fig_OD.savefig(f'{path_save+title} logOD replicates.jpg', dpi = 300, bbox_inches = 'tight')
#fig_der.savefig(f'{path_save+title} derivative replicates.jpg', dpi = 300, bbox_inches = 'tight')

#%% M. elsdenii 2

title = 'M. elsdenii 2 0719'
sheet_name = 'ME 2'

fig_OD, fig_der, negatives, time, me_mu_dict_2 = high_throughput_graphs_replicates2(
    sheet_name, title, cutoff = 0, colours = colours, drop = [])
#fig_OD.savefig(f'{path_save+title} logOD replicates.jpg', dpi = 300, bbox_inches = 'tight')
#fig_der.savefig(f'{path_save+title} derivative replicates.jpg', dpi = 300, bbox_inches = 'tight')

#%% M. elsdenii 3

title = 'M. elsdenii 3 0719'
sheet_name = 'ME 3'

fig_OD, fig_der, negatives, time, me_mu_dict_3 = high_throughput_graphs_replicates2(
    sheet_name, title, cutoff = 0, colours = colours, drop = [])
#fig_OD.savefig(f'{path_save+title} logOD replicates.jpg', dpi = 300, bbox_inches = 'tight')
#fig_der.savefig(f'{path_save+title} derivative replicates.jpg', dpi = 300, bbox_inches = 'tight')

#%% Monod plot elsdenii

title = '$\it{M. elsdenii}$'

colours = ['red', 'orange', 'black']
dict_list = [me_mu_dict_1, me_mu_dict_2, me_mu_dict_3]

fig_monod_bio_rep = growth_rates_graph_bio_rep(
    dict_list, 'M. elsdenii', title, new = True, ax = None, colours = colours)
fig_monod, ax_monod = growth_rates_graph(dict_list, title, new = True, ax = None, colour = 'red', label = None)

#fig_monod_bio_rep.savefig(f'{path_save+title} monod biological rep.jpg', dpi = 300, bbox_inches = 'tight')
#fig_monod.savefig(f'{path_save+title} monod.jpg', dpi = 300, bbox_inches = 'tight')

#%% M. hexanoica 1
colours = ['r', 'b', 'g', 'y']
title = 'M. hexanoica 1 0719'
sheet_name = 'MH 1'

fig_OD, fig_der, negatives, time, mh_mu_dict_1 = high_throughput_graphs_replicates2(
    sheet_name, title, cutoff = 10, colours = colours, drop = [])
#fig_OD.savefig(f'{path_save+title} logOD replicates.jpg', dpi = 300, bbox_inches = 'tight')
#fig_der.savefig(f'{path_save+title} derivative replicates.jpg', dpi = 300, bbox_inches = 'tight')

#%% M. hexanoica 2

title = 'M. hexanoica 2 0719'
sheet_name = 'MH 2'

fig_OD, fig_der, negatives, time, mh_mu_dict_2 = high_throughput_graphs_replicates2(
    sheet_name, title, cutoff = 10, colours = colours, drop = [])
#fig_OD.savefig(f'{path_save+title} logOD replicates.jpg', dpi = 300, bbox_inches = 'tight')
#fig_der.savefig(f'{path_save+title} derivative replicates.jpg', dpi = 300, bbox_inches = 'tight')

#%% M. hexanoica 3

title = 'M. hexanoica 3 0719'
sheet_name = 'MH 3'

fig_OD, fig_der, negatives, time, mh_mu_dict_3 = high_throughput_graphs_replicates2(
    sheet_name, title, cutoff = 10, colours = colours, drop = [2.5, 5])
#fig_OD.savefig(f'{path_save+title} logOD replicates.jpg', dpi = 300, bbox_inches = 'tight')
#fig_der.savefig(f'{path_save+title} derivative replicates.jpg', dpi = 300, bbox_inches = 'tight')

#%% Monod plot hexanoica

title = '$\it{M. hexanoica}$'
colours = ['blue', 'purple', 'green']

dict_list = [mh_mu_dict_1, mh_mu_dict_2, mh_mu_dict_3]

fig_monod_bio_rep = growth_rates_graph_bio_rep(
    dict_list, 'M. hexanoica', title, new = True, ax = None, colours = colours)

dict_list = [mh_mu_dict_1, mh_mu_dict_2] # Excluding M hexanoica 3

fig_monod, ax_monod = growth_rates_graph(dict_list, title, new = True, ax = None, colour = 'blue', label = None)

#fig_monod_bio_rep.savefig(f'{path_save+title} monod biological rep.jpg', dpi = 300, bbox_inches = 'tight')
#fig_monod.savefig(f'{path_save+title} monod.jpg', dpi = 300, bbox_inches = 'tight')

#%% Plot both growth rates (individual biological replicates)

title = 'Megapshaera 0719'

me_dict_list = [me_mu_dict_1, me_mu_dict_2, me_mu_dict_3]
mh_dict_list = [mh_mu_dict_1, mh_mu_dict_2, mh_mu_dict_3]
me_colours = ['red', 'orange', 'black']
mh_colours = ['blue', 'purple', 'green']


fig6, ax = plt.subplots()
growth_rates_graph_bio_rep(me_dict_list, 'M. elsdenii', title,
                   new = False, ax = ax, colours = me_colours)
growth_rates_graph_bio_rep(mh_dict_list, 'M. hexanoica', title, 
                   new = False, ax = ax, colours = mh_colours)
fig6.legend(loc = [0.72, 0.6])
#fig6.savefig(f'{path_save+title} monod combined.jpg', dpi = 300, bbox_inches = 'tight')

#%% Plot both growth rates averaged (without M hexanoica 3)

title = '$\it{Megasphaera}$ Specific Growth Rates'

me_dict_list = [me_mu_dict_1, me_mu_dict_2, me_mu_dict_3]
mh_dict_list = [mh_mu_dict_1, mh_mu_dict_2]
me_colour = 'red'
mh_colour = 'blue'

fig6, ax = plt.subplots(dpi = 300)
growth_rates_graph(me_dict_list, title = None, label = '$\it{M. elsdenii}$',
                   new = False, ax = ax, colour = me_colour)
growth_rates_graph(mh_dict_list, title = None, label = '$\it{M. hexanoica}$',
                   new = False, ax = ax, colour = mh_colour)
ax.set_title(title)
fig6.legend(loc = [0.72, 0.78])
#fig6.savefig(f'{path_save+title} monod combined.jpg', dpi = 300, bbox_inches = 'tight')
#fig6.savefig(f'{path_save+title} monod combined.eps', dpi = 300, bbox_inches = 'tight')
#fig6.savefig(f'{path_save+title} monod combined.png', dpi = 300, bbox_inches = 'tight')


num_to_keep = 18; [plt.close(n) for n in plt.get_fignums() if n != num_to_keep]
plt.show()