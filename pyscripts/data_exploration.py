# -*- coding: utf-8 -*-
"""
Created on Fri Aug 17 13:33:10 2018
Preliminary data analysis

Objective: 
@author: scampit
"""

import pandas as pd
import numpy as np
import scipy
from os import walk
import matplotlib.pyplot as plt
import seaborn as sns

# Clean up the sabio RK data to only have Km, S, Vmax, and kcat
kinetic_data_path = r'C:\Users\scampit\Desktop\Kinetic Model\Data\sabiodb'
drop = ['Tissue', 'CellularLocation', 'Enzymename', 'EnzymeType', 
        'parameter.standardDeviation', 'Cofactor', 'Inhibitor', 'ChebiID', 
        'KeggReactionID']
# kinetics data
kin_fil = r'homosapien_clean.csv'
kin = pd.read_csv(kinetic_data_path+'/'+kin_fil).drop(drop, axis=1)
kin_alt = pd.read_csv(kinetic_data_path+'/'+kin_fil)

# Using the KEGG annotation -- may want to look at other annotations later
KEGG_path = r'C:\Users\scampit\Desktop\Kinetic Model\Data\KEGG_Metabolism'
metabolic_pathways=[]
for path, nam, fil in walk(KEGG_path):
    metabolic_pathways.extend(fil)
    break
metabolic_pathways = [s.strip('KEGG') for s in metabolic_pathways]
metabolic_pathways = [s.strip('.txt') for s in metabolic_pathways]
metabolic_pathways = [s.replace('_', ' ') for s in metabolic_pathways]
metabolic_pathways = [s.lower() for s in metabolic_pathways]
metabolic_pathways = [s.title() for s in metabolic_pathways]

# specifically human genes from _ database
#human_gene_path = r'C:\Users\scampit\Desktop\Kinetic Model\Data\human_genes'
#human_gene_fil = r'genage_human.csv'
#hgene = pd.read_csv(human_gene_path+'\\'+human_gene_fil, usecols=['symbol']).drop_duplicates(keep='first')
#kin = pd.merge(kin, hgene, how='inner', left_on='Gene names', right_on='symbol')

###############################################################################
def feature_parser(kin, param_type=''):
    """
    feature_parser will output a csv file containing the kinetic parameter 
    inputted. 
    
    # Idea: have it output a heatmap of the top 10 features, its corresponding
    metabolite or gene, and the pathway it's associated with. 
    
    """

    df = kin.loc[kin['parameter.type'] == param_type]
        
    # construct the pandas dataframe differently depending on which kinetic 
    # parameter we're looking at. 
    if param_type == 'Km':
        df = df.groupby(['Gene names', 'parameter.associatedSpecies', 'Pathway'], 
                    as_index=False)['parameter.startValue'].mean()
        df = df.set_index('Gene names').rename(columns={'parameter.startValue':param_type, 
                     'parameter.associatedSpecies':'Substrate'})
        df[param_type] = np.log(df[param_type])*-1
        df = df.pivot_table(index='Pathway', columns='Substrate', values=param_type)
    elif param_type == 'S/Km':
        df = df.groupby(['Gene names', 'parameter.associatedSpecies', 'Pathway'])
    else:
        df = df.groupby(['Gene names', 'Pathway'], 
                as_index=False)['parameter.startValue'].mean()
        df = df.rename(columns={'Gene names':'Gene symbol', 'parameter.startValue':param_type})
        df[param_type] = np.log(df[param_type])
        df = df.pivot_table(index='Pathway', columns='Gene symbol', values=param_type)
    
    df = df.T
    df = df.fillna(value=0)
    
    #sns.set_style('dark')
    #sns.heatmap(data=df, cmap='RdBu_r', robust=True, annot=True, linewidths=0.5, 
    #            xticklabels=True, yticklabels=True)
    #plt.title("Average metabolite Km values within KEGG pathways")
    #plt.figure(figsize=(15,15))
    #return plt.show()
    
    save_path = r'C:\Users\scampit\Desktop\Kinetic Model\Data\sabiodb\pathway_association'
    if param_type == 'kcat/Km':
        df.to_csv(save_path+'\\'+'substrate_sensor_metabolic.csv')
    else:
        df.to_csv(save_path+'\\'+param_type+'_metabolic.csv')
    return df

# extract and take the average kinetics data by genes 

feature_parser(kin, param_type='Km')
feature_parser(kin, param_type='kcat')
feature_parser(kin, param_type='kcat/Km')

###############################################################################
# for S/Km, I will map it by substrate:


###############################################################################
# for v/Vmax, I need to calculate fluxes, which is going to be a little more challenging

kcat = kin.loc[kin['parameter.type'] == 'kcat']
kcat = kcat.groupby(['Gene names', 'Pathway'], as_index=False)['parameter.startValue'].mean()
kcat = kcat.set_index('Gene names').rename(columns={'parameter.startValue':'kcat'})

vmax = kin.loc[kin['parameter.type'] == 'Vmax']
vmax = vmax.groupby('Gene names', as_index=False)['parameter.startValue'].mean()
vmax = vmax.set_index('Gene names').rename(columns={'parameter.startValue':'vmax'})

merged = pd.merge(km, kcat, how='inner', left_index=True, right_index=True)
km_kcat = pd.DataFrame(merged['kcat'].div(merged['Km']))
merged2 = pd.merge(km_kcat, vmax, how='inner', left_index=True, right_index=True)
merged2 = merged2.rename(columns={0:'kcat/Km'})
merged = pd.merge(merged, vmax, how='inner', left_index=True, right_index=True)

merged.to_csv('kinetics_WT.csv')

###############################################################################

# v/Vmax 

# S/Km

# kcat/Km
kcat_km = kin.loc[kin['parameter.type'] == 'kcat/Km'].set_index('Gene names')

###############################################################################

