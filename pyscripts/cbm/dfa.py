# -*- coding: utf-8 -*-
"""
Created on Wed May 16 19:52:46 2018
Dynamic Flux Analysis adapted for Python

@author: Sriram Chandrasekaran & Scott Campit
"""

from __future__ import print_function
import os, cobra, lxml, gurobipy, json
from os.path import join
from time import time
import escher, cameo

# Let's use RECON3D since that is the latest version
#data_path = r'C:\Users\scott\OneDrive\Desktop\GeneRegulation'
data_path = r'C:\Users\scampit\Desktop\Kinetic Model\Data\Recon'
recon3fil = r'Recon3D.xml'

model = cobra.io.read_sbml_model(join(data_path, recon3fil))
model.solver='gurobi'

# okay let's check to see our model loaded properly:
print("Number of reactions: ", len(model.reactions))
print("Number of metabolites: ", len(model.metabolites))
print("Number of genes: ", len(model.genes))

# Now let's run a simulation
solution = model.optimize(objective_sense='maximize')
print("Here is the solution from maximizing the biomass function: ", model.optimize().objective_value)
print("From maximizing the biomass function, here are the resulting flux values: ", solution.fluxes)
print("Model Objective: ", model.objective)

# Escher visualizations #######################################################
escher.list_available_maps()

# available maps for humans:
# {'organism': 'Homo sapiens','map_name': 'RECON1.Amino acid metabolism (partial)'},
# {'organism': 'Homo sapiens','map_name': 'RECON1.Carbohydrate metabolism'},
# {'organism': 'Homo sapiens','map_name': 'RECON1.Glycolysis TCA PPP'},
# {'organism': 'Homo sapiens','map_name': 'RECON1.Inositol retinol metabolism'},
# {'organism': 'Homo sapiens','map_name': 'RECON1.Tryptophan metabolism'},

# We'll build a glycolysis map: 

# Only works in Jupyter notebook
b = escher.Builder(map_name='RECON1.Glycolysis TCA PPP')
b.display_in_notebook()

b = escher.Builder(map_name='RECON1.Glycolysis TCA PPP', reaction_data=dict(solution.fluxes), # change the default colors
                   reaction_scale=[{'type': 'min', 'color': '#cccccc', 'size': 4},
                                   {'type': 'value', 'value': 0.1, 'color': '#cccccc', 'size': 8},
                                   {'type': 'mean', 'color': '#0000dd', 'size': 20},
                                   {'type': 'max', 'color': '#ff0000', 'size': 40}],
                   # absolute value and no text for data
                   reaction_styles=['size', 'color', 'abs'],
                   # only show the primary metabolites
                   hide_secondary_metabolites=True)
b.display_in_notebook()

# Load Metabolomics data into Escher:
import pandas as pd
fil = r'metabolomics.csv' # not real
metabolomics = pd.read_csv(fil).to_dict()

b = escher.Builder(map_name='RECON1.Glycolysis TCA PPP', reaction_data=dict(solution.fluxes), # change the default colors
                   reaction_scale=[{'type': 'min', 'color': '#cccccc', 'size': 4},
                                   {'type': 'value', 'value': 0.1, 'color': '#cccccc', 'size': 8},
                                   {'type': 'mean', 'color': '#0000dd', 'size': 20},
                                   {'type': 'max', 'color': '#ff0000', 'size': 40}],
                   # absolute value and no text for data
                   reaction_styles=['size', 'color', 'abs'],
                   # only show the primary metabolites
                   hide_secondary_metabolites=True)
b.display_in_notebook()

S = cobra.util.array.create_stoichiometric_matrix(model, array_type="dense")
obj = cobra.util.solver.linear_reaction_coefficients(model)
rhs = model.optimize().fluxes
lb
ub
vtype


shadow = model.optimize().shadow_prices
b = model.optimize().fluxes
c = model.optimize().objective_value

'''
class Dynamic_Flux_Analysis():


def flux_activity_coeff(model, timecourse_metabolomics, sheetname='Sheet1', kappa=None, kappa2=None, genedelflag=None, rxndelflag=None, normquantile=None):
    docstring
    flux_activity_coeff


    # default params
    if kappa is None:
        kappa = 1
    if kappa2 is None:
        kappa2 = 1E-3
    if sheetname is None:
        sheetname = 'Sheet1'
    if genedelflag is None:
        genedelflag = 0
    if rxndelflag is None:
        rxndelflag = 0
    if normquantile is None:
        normquantile = 0
'''
