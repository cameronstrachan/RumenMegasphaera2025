#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 11:45:16 2024

@author: connor

Calculate enzyme cost per ATP for all products. Plot enzyme costs and ATP yields.
Note: Acetate assimilation refers to consumption of one mole of acetate for 
each mole of RBO product, whereas recycling refers to regeneration of one mole
of acetate from acetyl-CoA for each mole of product.
"""

import os
os.chdir(os.path.dirname(__file__))
from equilibrator_api import ComponentContribution, Q_
cc = ComponentContribution()

#%% User defined function for enzyme cost minimization

import pandas as pd
from equilibrator_api import Reaction, R, default_T
from equilibrator_pathway import ThermodynamicModel, EnzymeCostModel, EnzymeCostFunction
import numpy as np
import matplotlib.pyplot as plt

def reactions_from_excel(reactions, compounds, cc):
    
    reaction_dict = {}
    redox_list = []
    proton_list = []
    
    reduction_equivalents = {'fd_ox' : 'ferredoxin',
                             'nad' : 'nadh', 
                             'nadp' : 'nadph',
                             'rub_ox' : 'rubredoxin', 
                             'ferro_c' : 'cytochrome c'}

    reduction_equivalents_names = list(reduction_equivalents.keys())

    oxidation_equivalents = {'fd_red' : 'ferredoxin',
                             'nadh' : 'nadh', 
                             'nadph' : 'nadph',
                             'rub_red' : 'rubredoxin', 
                             'ferri_c' : 'cytochrome c'}

    oxidation_equivalents_names = list(oxidation_equivalents.keys())
    all_equivalents = reduction_equivalents_names + oxidation_equivalents_names
    
    for i in reactions.index:
        
        rxn_name = reactions['Name'].iloc[i]
        reactants = reactions['Reaction'].iloc[i].split(' = ')[0].split(' + ')
        products = reactions['Reaction'].iloc[i].split(' = ')[1].split(' + ')
        dic = {}
        
        for reactant in reactants:
            coeff = reactant.split(' ')[0]
            name = reactant.split(' ')[1]
           
            if name in reduction_equivalents_names:
                redox_list += [(rxn_name, reduction_equivalents[name], 'red', coeff)]
            
            elif name in oxidation_equivalents:
                redox_list += [(rxn_name, oxidation_equivalents[name], 'ox', coeff)]
                
            else:  
                ID = compounds['ID'].loc[name]
                compound = cc.get_compound(ID)
                dic[compound] = -int(coeff)
                
        for product in products:
            coeff = product.split(' ')[0]
            name = product.split(' ')[1]

            if name not in all_equivalents:
                ID = compounds['ID'].loc[name]
                compound = cc.get_compound(ID)
                dic[compound] = int(coeff)
                
        reaction_dict[rxn_name] = Reaction(dic)
        
    for i in reactions.index:
        if reactions['IMF'].fillna('None').iloc[i] != 'None':
            rxn_name = reactions['Name'].iloc[i]
            numb_protons = reactions['IMF'].iloc[i]
            if reactions['IMF'].iloc[i] > 0:
                direction = 'out'
            elif reactions['IMF'].iloc[i] < 0:
                direction = 'in'
        
            proton_list += [(rxn_name, numb_protons, direction)]
    
    return reaction_dict, redox_list, proton_list

def init_thermo_model(reaction_dict, # Dict{Reaction Name : Reaction Object}
                      compound_dict, # Dict{Compound Name : Compound Object}
                      redox_list,
                      proton_list,
                      fluxes,
                      lower_bounds,
                      upper_bounds,
                      E_fe,
                      E_rub,
                      E_cyt,
                      E_nadh,
                      E_nadph,
                      pmf,
                      config_dict,
                      cc):
        
    
    # Free energies of reducing equivalent oxidations
    
    G_fe = 2*96.485*E_fe # kJ/mol
    G_rub = 2*96.485*E_rub
    G_cyt = 1*96.485*E_cyt
    G_nadh = 2*96.485*E_nadh # kJ/mol 
    G_nadph = 2*96.485*E_nadph # kJ/mol
    
    redox_dict = {'ferredoxin' : G_fe,
                  'rubredoxin' : G_rub,
                  'cytochrome c' : G_cyt,
                  'nadh' : G_nadh,
                  'nadph' : G_nadph}
    
    # Free energy change of one mole of protons crossing from inside of the membrane to outside
    
    G_prot = 96.485*pmf # kJ/mol proton, F = 96.485 Coulombs per mole proton, pmf in Volts

    reaction_name_list = list(reaction_dict.keys())
    reaction_list = list(reaction_dict.values())
    compound_name_list = list(compound_dict.keys())
    compound_list = list(compound_dict.values())

    S = cc.create_stoichiometric_matrix_from_reaction_objects(reaction_list)
    S.columns = reaction_name_list
    S.index = [str(compound) for compound in S.index]
    compound_dict_swap = {str(compound) : name for name, compound in compound_dict.items()}
    S.index = [compound_dict_swap[i] for i in S.index]

    # Convert fluxes to dimensionless quantity objects with pint
    
    fluxes = Q_(fluxes, '')

    # Initialize thermo model with explicit standard free energies

 
    thermo_model_explicit = ThermodynamicModel(S = S, 
                                      compound_dict = compound_dict, 
                                      reaction_dict = reaction_dict, 
                                      fluxes = fluxes, 
                                      config_dict = config_dict,
                                      comp_contrib = cc)
    
    thermo_model_explicit.fluxes = fluxes

    # Set standard free energies of reactions to allow for modification of redox reactions

    standard_dg_primes_new, dg_sigma_new = cc.standard_dg_prime_multi(reaction_list, 'fullrank')

    # Add free energies of ferredoxin redox

    for tup in redox_list:
        
        rxn_name = tup[0]
        redox_equivalent = tup[1]
        direction = tup[2]
        coeff = int(tup[3])
        position = [x for x, y in enumerate(reaction_name_list) if y == rxn_name]
        
        if direction == 'ox':
            thermo_model_explicit.standard_dg_primes[position] += coeff*Q_(redox_dict[redox_equivalent], 'kJ/mol')
        
        if direction == 'red':
            thermo_model_explicit.standard_dg_primes[position] -= coeff*Q_(redox_dict[redox_equivalent], 'kJ/mol')
        
    
    for tup in proton_list:
        # Get position of reaction in reaction list
        rxn_name = tup[0]
        position = [x for x, y in enumerate(reaction_name_list) if y == rxn_name]
        proton_number = abs(tup[1])
        
        if tup[2] == 'out':
            thermo_model_explicit.standard_dg_primes[position] += Q_(proton_number*G_prot, 'kJ/mol')
        elif tup[2] == 'in':
            thermo_model_explicit.standard_dg_primes[position] -= Q_(proton_number*G_prot, 'kJ/mol')

    # Set concentration bounds explicitly

    for compound in compound_name_list:
        thermo_model_explicit.set_bounds(compound, lower_bounds[compound], upper_bounds[compound])

    thermo_model_explicit.reaction_list = reaction_list
    thermo_model_explicit.reaction_dict = reaction_dict
    thermo_model_explicit.reaction_names = reaction_name_list
    thermo_model_explicit.compound_names = compound_name_list
    thermo_model_explicit.compound_list = compound_list
    thermo_model_explicit.compound_dict = compound_dict

    return thermo_model_explicit

def MDF_lite(thermo_model):
    
    mdf_sol = thermo_model.mdf_analysis()
    mdf_sol.original_standard_dg_prime = mdf_sol.reaction_df['original_standard_dg_prime'] # these values are without correction for uncertainty
    mdf_sol.standard_dg_prime = mdf_sol.reaction_df['standard_dg_prime']
    mdf_sol.physiological_dg_prime = mdf_sol.reaction_df['physiological_dg_prime']
    mdf_sol.optimized_dg_prime = mdf_sol.reaction_df['optimized_dg_prime']                     
    mdf_sol.concentrations = {mdf_sol.compound_df['compound_id'].iloc[i] : 
                              mdf_sol.compound_df['concentration_in_molar'].iloc[i] 
                              for i in range(len(mdf_sol.compound_df.index))}

    min_force_reactions = [i for i in mdf_sol.reaction_df[mdf_sol.reaction_df.shadow_price != 0].index]
    min_forces = [mdf_sol.reaction_df['optimized_dg_prime'].iloc[i].m for i in min_force_reactions]
    
    return -min_forces[0], mdf_sol

def ECM_sum(ecm_sol, protein_mass_units):
    
    prot = 0
    for i in range(len(ecm_sol.protein_mass_array)):
        factor = np.prod(ecm_sol.costs[i, :])
        prot += factor*ecm_sol.protein_mass_array[i]
        
    if protein_mass_units == 'Da':
        return round(prot/(3600*1000), 3) # g protein / (mM product / h)
    
    elif protein_mass_units == 'kDa':
        return round(prot/(3600), 3) # g protein / (mM product / h)
    
    else:
        print('Protein Mass Units not recognized')
        
def ECM_energies(ecm_sol, thermo_model, title, path, scale):
    
    ecm_sol.original_standard_dg_prime = ecm_sol.reaction_df['original_standard_dg_prime'] # these values are without correction for uncertainty
    ecm_sol.standard_dg_prime = ecm_sol.reaction_df['standard_dg_prime']
    ecm_sol.physiological_dg_prime = ecm_sol.reaction_df['physiological_dg_prime']
    ecm_sol.optimized_dg_prime = ecm_sol.reaction_df['optimized_dg_prime']                
    
    ecm_sol.concentrations = {ecm_sol.compound_df['compound_id'].iloc[i] : 
                              ecm_sol.compound_df['concentration_in_molar'].iloc[i] 
                              for i in range(len(ecm_sol.compound_df.index))}
    
    fig, ax = plt.subplots()
    x = np.arange(len(thermo_model.reaction_list))

    min_force_reactions = []
    
    phys_dg = [i.m for i in list(ecm_sol.reaction_df['physiological_dg_prime'])]
    ax.bar([scale*i for i in x], phys_dg, width = 0.2, 
           label = "∆G$^{0}$\' : Physiological Conc. (1 mM)", color = 'xkcd:grey')

    opt_dg = [i.m for i in list(ecm_sol.reaction_df['optimized_dg_prime'])]
    ax.bar([scale*i for i in x if i not in min_force_reactions], [g for i, 
            g in enumerate(opt_dg) if i not in min_force_reactions], 
           width = 0.4, label = "∆G\' : ECM Optimized Conc.", color = 'xkcd:blue', alpha = 0.7)

    ax.legend()
    ax.set_title(f'{title}')
    ax.set_ylabel('kJ/mol')
    ax.grid(axis = 'x', color = 'xkcd:light grey', linewidth = 0.75)
    ax.set_xlim(scale*x[0]-1*scale, (scale*x[-1]+scale*1))
    ax.plot(np.append((np.append(x[0]-1*scale, x)), (x[-1]+scale*1)) , np.zeros(len(x)+2), color = 'xkcd:light grey', linewidth = 0.75) 
    plt.xticks(scale*x, thermo_model.reaction_names, rotation = 90, fontsize = 8)
    plt.savefig(path+title+'.jpg', dpi = 300, bbox_inches = 'tight')
    
    return ecm_sol

def ECM(thermo_model,
        new_dgs,
        Km_df,
        Km_units,
        crc_array,
        crc_units,
        protein_mass_array,
        protein_mass_units,
        compound_mass_array,
        compound_mass_units,
        title, 
        path, 
        colors):
    
    """ ECM from excel model file. Allows for specification of kinetic constants
    for individual reactions """
    
    ### Create Enzyme Cost Model ###
    
    reaction_list = thermo_model.reaction_list
    reaction_names = thermo_model.reaction_names
    compound_names = thermo_model.compound_names 
    compound_list = thermo_model.compound_list     
    # Parameter Dataframe
    
    # Need to reformate Km dataframe
    
    param_df = pd.DataFrame()
    param_df["QuantityType"] = ['Michaelis constant']*len(Km_df.index)
    param_df['Value'] = Km_df['Michaelis Constant']
    param_df['Compound'] = Km_df['Compound']
    param_df['Reaction'] = Km_df['Reaction']
    param_df['Unit'] = Km_df['Unit']
    
    # Need substrate catalytic rate constant for each reaction
    
    crc_rows = len(reaction_list)
    crc_df = pd.DataFrame([], columns = ["QuantityType", "Value", "Compound", "Reaction", "Unit"])
    crc_df["QuantityType"] = ['substrate catalytic rate constant']*crc_rows
    crc_df["Value"] = crc_array
    crc_df['Unit'] = [crc_units]*crc_rows
    crc_df['Reaction'] = reaction_names
    
    param_df = pd.concat([param_df, crc_df], ignore_index = True, copy = False)
    
    # Need protein and compound molecular mass
    
    prot_mass_rows = len(reaction_list)
    prot_mass_df = pd.DataFrame([], columns = ["QuantityType", "Value", "Compound", "Reaction", "Unit"])
    prot_mass_df["QuantityType"] = ['protein molecular mass']*prot_mass_rows
    prot_mass_df["Value"] = protein_mass_array
    prot_mass_df['Unit'] = [protein_mass_units]*prot_mass_rows
    prot_mass_df['Reaction'] = reaction_names
    
    param_df = pd.concat([param_df, prot_mass_df], ignore_index = True, copy = False)
    
    comp_mass_rows = len(compound_list)
    comp_mass_df = pd.DataFrame([], columns = ["QuantityType", "Value", "Compound", "Reaction", "Unit"])
    comp_mass_df["QuantityType"] = ['molecular mass']*comp_mass_rows
    comp_mass_df["Value"] = compound_mass_array
    comp_mass_df['Unit'] = [compound_mass_units]*comp_mass_rows
    comp_mass_df['Compound'] = compound_names
    
    param_df = pd.concat([param_df, comp_mass_df], ignore_index = True, copy = False)
    
    ### Solve ECM problem ###
    
    model = EnzymeCostModel(thermo_model, param_df)
    
    # Convert new_dgs to pint.Quantity object
    new_dgs_mod = Q_(list([i.m for i in new_dgs]), 'kJ/mol')
    
    # Reshape
    new_dgs_mod_reshaped = np.reshape(new_dgs_mod, (len(model.reaction_ids), 1))
    dir_mat = np.diag(np.sign(model._thermo_model.fluxes.magnitude + 1e-10).flat)
    
    # Set standard_dg_over_rt of ecf object to new values
    model.ecf.standard_dg_over_rt = ((dir_mat @ new_dgs_mod_reshaped) / (R * default_T)).m_as("").flatten()
    np.random.seed(1982)
    ecm_sol = model.optimize_ecm()
    
    ecm_sol.protein_mass_array = protein_mass_array
    ecm_sol.protein_mass_units = protein_mass_units
    ecm_sol.model = model
    ecm_sol.reaction_names = reaction_names
    
    protein_cost = ECM_sum(ecm_sol, protein_mass_units)
    fig1, ax = plt.subplots(1, 1, figsize=(7, 5))
    
    costs = ecm_sol.costs
    base = min(filter(None, costs[:, 0])) / 2.0
    idx_zero = costs[:, 0] == 0
    costs[idx_zero, 0] = base
    costs[idx_zero, 1:] = 1.0

    bottoms = np.hstack([np.ones((costs.shape[0], 1)) * base, np.cumprod(costs, 1)])
    steps = np.diff(bottoms)
    
    top_level = 3
    
    labels = EnzymeCostFunction.ECF_LEVEL_NAMES[0:top_level]

    ind = range(costs.shape[0])  # the x locations for the groups
    width = 0.8
    ax.set_yscale("log")
    for i, label in enumerate(labels):
        ax.bar(
            ind,
            steps[:, i].flat,
            width,
            bottom=bottoms[:, i].flat,
            color=colors[i],
            label=label,
            alpha = 0.9)

    ax.set_xticks(ind)
    ax.set_xticklabels(ecm_sol.reaction_ids, size="medium", rotation=90)
    ax.legend(loc="best", framealpha=0.2)
    ax.set_title(f"{title}\nEnzyme Cost: {round(protein_cost, 3)} g/mmol/h")
    ax.set_ylabel("Enzyme Demand [M]")
    ax.set_ylim(bottom=base)
    
    plt.savefig(path+title+'.jpg', dpi = 300, bbox_inches = 'tight')
    
    ECM_energies(ecm_sol, thermo_model, title+' Energies', path, scale = 0.8)
    
    return protein_cost, ecm_sol, param_df.fillna('None').to_string()

#%% ECM initiation

plt.rcdefaults()

# Directory for saving images and model files

ECM_path_save = '/Users/connor/Documents/Projects/Lactate growth (Roy and Cam)/\
Paper/ECM graphs/'

# Simulation Configuration

config_dict = {'flux_unit': 'Unitless',
               'version': '3',
               'kcat_source': 'fwd',
               'denominator': 'CM',
               'regularization': 'volume', 
               'objective': 'enzyme',
               'standard_concentration': '1 M', 
               'solver': 'CLARABEL', 
               'algorithm': 'ECM', 
               'p_h': '6',
               'ionic_strength': '250 mM',
               'p_mg': '3',
               'dg_confidence' : '0'}

E_fe = -0.400 # Ferredoxin, V

colors = ["xkcd:grey", "xkcd:red", "xkcd:blue"]

def wrap(model_file, title):

    # Standard Concentration Bounds

    standard_lb = 0.001
    standard_ub = 10

    # Specify reduction potentials (V)

    E_rub = -0.075 # Rubredoxin
    E_cyt = 0.254 # Cytochrome C
    E_nadh = -0.280 # NADH
    E_nadph = -0.370 # NADPH

    # Specify chemical potential of proton crossing membrane (in to out)

    compounds = pd.read_excel(model_file, sheet_name = 'Compounds', index_col = 0)
    reactions = pd.read_excel(model_file, sheet_name = 'Reactions')
        
    reaction_dict, redox_list, proton_list = reactions_from_excel(reactions, compounds, cc)
    compound_dict = {name : cc.get_compound(compounds['ID'].loc[name]) for name in compounds.index}
    lower_bounds = compounds['Lower Bound'].fillna(standard_lb).apply(lambda x: Q_(x, 'mM'))
    upper_bounds = compounds['Upper Bound'].fillna(standard_ub).apply(lambda x: Q_(x, 'mM'))
    fluxes = reactions['Flux'].to_list()

    thermo_model = init_thermo_model( reaction_dict, # Dict{Reaction Name : Reaction Object}
                                      compound_dict, # Dict{Compound Name : Compound Object}
                                      redox_list,
                                      proton_list,
                                      fluxes,
                                      lower_bounds,
                                      upper_bounds,
                                      E_fe,
                                      E_rub,
                                      E_cyt,
                                      E_nadh,
                                      E_nadph,
                                      pmf,
                                      config_dict,
                                      cc
                                      )

    mdf, mdf_sol = MDF_lite(thermo_model)
    
    Km_units = 'mM'
    crc_units = '1/s'
    protein_mass_units = 'kDa'
    compound_mass_units = 'Da'
    crc_array = reactions['Forward Catalytic Constant'].to_list()
    protein_mass_array = reactions['Protein Mass'].to_list()
    compound_mass_array = compounds['Compound Mass'].to_list()
    Km_df = pd.read_excel(model_file, sheet_name = 'Michaelis Constants')
    new_dgs = mdf_sol.standard_dg_prime
    
    protein_cost, ecm_sol, params = ECM(thermo_model,
                                    new_dgs,
                                    Km_df,
                                    Km_units,
                                    crc_array,
                                    crc_units,
                                    protein_mass_array,
                                    protein_mass_units,
                                    compound_mass_array,
                                    compound_mass_units,
                                    title, 
                                    ECM_path_save, 
                                    colors
                                    )
    
    return protein_cost

dic = {}

#%% Butyrate (Recycling)

pmf = 0.200 # V

model_file = '/Users/connor/Documents/Projects/Lactate growth (Roy and Cam)/Paper\
/Model Files/Lactate to Butyrate SLP with MWs.xlsx'

title = 'Butyrate (Acetate Recycling)'

cost = wrap(model_file, title)
atp = 0.25
dic[title] = (cost / atp, atp)



#%% Caproate (Recycling)

model_file = '/Users/connor/Documents/Projects/Lactate growth (Roy and Cam)/Paper\
/Model Files/Lactate to Caproate SLP with MWs.xlsx'
 
title = 'Caproate (Acetate Recycling)'

cost = wrap(model_file, title)
atp = 0.33
dic[title] = (cost / atp, atp)

#%% Octanoate (Recycling)

pmf = 0.100

model_file = '/Users/connor/Documents/Projects/Lactate growth (Roy and Cam)/Paper\
/Model Files/Lactate to Octanoate SLP with MWs.xlsx'

title = 'Octanoate (Acetate Recycling)'

cost = wrap(model_file, title)
atp = 0.375
dic[title] = (cost / atp, atp)

#%% Butyrate (Assimilating)

pmf = 0.100 # V

model_file = '/Users/connor/Documents/Projects/Lactate growth (Roy and Cam)/Paper\
/Model Files/Lactate to Butyrate IMF with MWs.xlsx'

title = 'Butyrate (Acetate Assimilation)'

cost = wrap(model_file, title)
atp = 0.5
dic[title] = (cost / atp, atp)

#%% Caproate (Assimilating)

model_file = '/Users/connor/Documents/Projects/Lactate growth (Roy and Cam)/Paper\
/Model Files/Lactate to Caproate IMF with MWs.xlsx'
 
title = 'Caproate (Acetate Assimilation)'

cost = wrap(model_file, title)
atp = 0.5
dic[title] = (cost / atp, atp)

#%% Octanoate (Assimilating)

pmf = 0.100

model_file = '/Users/connor/Documents/Projects/Lactate growth (Roy and Cam)/Paper\
/Model Files/Lactate to Octanoate IMF with MWs.xlsx'

title = 'Octanoate (Acetate Assimilation)'

cost = wrap(model_file, title)
atp = 0.5
dic[title] = (cost / atp, atp)

#%% Propionate (Acrylate Pathway)

pmf = 0.100

model_file = '/Users/connor/Documents/Projects/Lactate growth (Roy and Cam)/Paper\
/Model Files/Lactate to Propionate (Acrylate Pathway).xlsx'

title = 'Propionate'

cost = wrap(model_file, title)
atp = 0.33
dic[title] = (cost / atp, atp)

#%% Bar plot

print(dic)

fig, axes = plt.subplots(2, 1, figsize = (6, 6), constrained_layout = True, sharex = True)
path_save = '/Users/connor/Documents/Projects/Lactate growth (Roy and Cam)/Paper/Plots/'
width = 0.4
color1 = 'xkcd:grey'
color2 = 'xkcd:dark grey'
x = ['Propionate', 'Butyrate', 'Caproate', 'Octanoate']

### Enzyme Cost ###

ax = axes[0]

# Enzyme cost without recycling

groups = ['Propionate', 'Butyrate (Acetate Assimilation)', 'Caproate (Acetate Assimilation)', 
          'Octanoate (Acetate Assimilation)']

costs1 = [dic[group][0] for group in groups]

ax.bar(x, costs1, width = width, color = color1)

# Enzyme cost with recycling
groups = ['Butyrate (Acetate Recycling)', 'Caproate (Acetate Recycling)', 
          'Octanoate (Acetate Recycling)']

bottoms = costs1

costs2 = [dic[group][0] for group in groups]

costs2 = [costs2[i] - costs1[i] for i in range(len(costs2))]

costs2.insert(0, 0)

ax.bar(x, costs2, width = width, bottom = bottoms, color = color2)
ax.set_ylabel('Enzyme Cost per ATP\n(g/mmol/hr)')

### ATP yield ###

ax = axes[1]

# ATP yield without recycling
groups = ['Butyrate (Acetate Recycling)', 'Caproate (Acetate Recycling)', 
          'Octanoate (Acetate Recycling)']


yields1 = [dic[group][1] for group in groups]

yields1.insert(0, 0)

ax.bar(x, yields1, width = width, color = color2, label = 'Acetate\nRecycling')

# ATP yield with recycling

groups = ['Propionate', 'Butyrate (Acetate Assimilation)', 'Caproate (Acetate Assimilation)', 
          'Octanoate (Acetate Assimilation)']

bottoms = yields1

yields2 = [dic[group][1] for group in groups]

yields2 = [yields2[i] - yields1[i] for i in range(len(yields2))]

ax.bar(x, yields2, width = width, bottom = bottoms, color = color1)
ax.set_ylabel('ATP Yield\non Lactate')

h, l = ax.get_legend_handles_labels()
axes[0].legend(h, l)
axes[0].set_title('Enzyme Costs and ATP Yields\nof Lactate Products')
fig.savefig(path_save+'Cost and yield.eps', dpi = 300)
fig.savefig(path_save+'Cost and yield.jpg', dpi = 300)

#%% Bar plot only assimilation

print(dic)

fig, axes = plt.subplots(2, 1, figsize = (6, 6), constrained_layout = True, sharex = True)
path_save = '/Users/connor/Documents/Projects/Lactate growth (Roy and Cam)/Paper/Plots/'
width = 0.4
color1 = 'xkcd:grey'
color2 = 'xkcd:dark grey'
x = ['Propionate', 'Butyrate', 'Caproate', 'Octanoate']

### Enzyme Cost ###

ax = axes[0]

# Enzyme cost without recycling

groups = ['Propionate', 'Butyrate (Acetate Assimilation)', 'Caproate (Acetate Assimilation)', 
          'Octanoate (Acetate Assimilation)']

costs1 = [dic[group][0] for group in groups]

ax.bar(x, costs1, width = width, color = color2)
ax.set_ylabel('Enzyme Cost per ATP\n(g/mmol/hr)')

### ATP yield ###

ax = axes[1]

# ATP yield without recycling
groups = ['Propionate', 'Butyrate (Acetate Assimilation)', 'Caproate (Acetate Assimilation)', 
          'Octanoate (Acetate Assimilation)']

bottoms = yields1

yields2 = [dic[group][1] for group in groups]

# ax.bar(x, yields1, width = width, color = color2, label = 'Acetate\nRecycling')

ax.bar(x, yields2, width = width, color = color2)
ax.set_ylabel('ATP Yield\non Lactate')

axes[0].set_title('Enzyme Costs and ATP Yields\nof Lactate Products')
fig.savefig(path_save+'Cost and yield (only assimilation).eps', dpi = 300)

