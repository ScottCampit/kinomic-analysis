# -*- coding: utf-8 -*-
"""
Created on Wed Aug 22 12:18:05 2018
Kinetic model for glycolysis

@author: scampit
"""

import pandas as pd
import numpy as np
import scipy, os

# objective - map kcat, vmax, and km to glycolytic genes

# Sabio RK Data is in 'Km Model'
sabiork_path = r'C:\Users\scampit\Desktop\Kinetic Model\Data\sabiodb'
kin_fil = r'homosapien_clean.csv'
kin = pd.read_csv(sabio_rk+'/'+kin_fil)

km = kin.loc[kin['parameter.type'] == 'Km'].drop_duplicates(keep='first')
km = km.groupby('Gene names', as_index=False)['parameter.startValue'].mean()
km = km.set_index('Gene names').rename(columns={'parameter.startValue':'Km'})

kcat = kin.loc[kin['parameter.type'] == 'kcat'].drop_duplicates(keep='first')
kcat = kcat.groupby('Gene names', as_index=False)['parameter.startValue'].mean()
kcat = kcat.set_index('Gene names').rename(columns={'parameter.startValue':'kcat'})

vmax = kin.loc[kin['parameter.type'] == 'Vmax'].drop_duplicates(keep='first')
vmax = vmax.groupby('Gene names', as_index=False)['parameter.startValue'].mean()
vmax = vmax.set_index('Gene names').rename(columns={'parameter.startValue':'vmax'})

kcat_km_orig = kin.loc[kin['parameter.type'] == 'kcat/Km'].drop_duplicates(keep='first')
kcat_km_orig = kcat_km_orig.groupby('Gene names', as_index=False)['parameter.startValue'].mean()
kcat_km_orig = kcat_km_orig.set_index('Gene names').rename(columns={'parameter.startValue':'kcat/Km'})

# trying to see how complete the data is for human kinetics:
global_match_part = pd.merge(km, kcat, how='inner', left_index=True, right_index=True).dropna() # ~500 genes
covar = global_match_part.cov()
corr = global_match_part.corr()

# Calculate the kcat/km using provided data from SabioRK
kcat_km_calc = pd.DataFrame()
kcat_km_calc['kcat/Km'] = (global_match_part['kcat']/global_match_part['Km']) # There is a large discrepancy between reported and calculated.


val = pd.merge(kcat_km_calc, kcat_km_orig, how='inner', left_index=True, right_index=True)

#global_match_all = pd.merge(global_match, vmax, how='inner', left_index=True, right_index=True).dropna() # ~260 genes

# mapping glycolysis onto sabioRK data
#glycolysis_biocarta = pd.read_csv('BIOCARTA_Glycolysis.txt', skiprows=1)
#glycolysis_kegg = pd.read_csv('KEGG_GLYCOLYSIS_GLUCONEOGENESIS.txt', skiprows=1)
#match_biokm = pd.merge(how='inner', left=km, right=glycolysis_biocarta, left_index=True, right_on='> Glycolysis Pathway')
#match_keggkcat = pd.merge(how='inner', left=kcat, right=glycolysis_kegg, left_index=True, right_on='> Glycolysis / Gluconeogesis')
# Results: There is no good map of glycolysis on MSIGDB

# Manual curation of glycolysis genes including isoforms
glycolysis = ['HK1', 'HK2', 'HK3', 'GCK', 'GPI', 'PGI', 'PHI', 'PFKM', 'PFKL', 'PFKP', 'ALDOA', 'ALDOB', 'ALDOC', 'GAPDH', 'PGK1', 'PGK2', 'PGAM1', 'PGAM2', 'PGAM4', 'ENO1', 'ENO2', 'ENO3', 'PKLR', 'PKM2', 'ACSS1', 'ACSS2', 'TPI1']
glycolysis = pd.DataFrame({'Glycolysis':glycolysis})
match_km = pd.merge(how='inner', left=km, right=glycolysis, left_index=True, right_on='Glycolysis')
match_km = match_km.set_index('Glycolysis')
match_kcat = pd.merge(how='inner', left=kcat, right=glycolysis, left_index=True, right_on='Glycolysis')
match_kcat = match_kcat.set_index('Glycolysis')
match = pd.merge(match_km, match_kcat, how='inner', left_index=True, right_index=True)
# Results: there are only 5 complete annotations for both the Km and the kcat for glycolysis.

# mapping TCA cycle onto sabioRK data
tca_kegg = pd.read_csv('KEGG_CITRATE_CYCLE_TCA_CYCLE.txt', skiprows=1)
match_tca = pd.merge(how='inner', left=tca_kegg, right=global_match_part, left_on='> Citrate cycle (TCA cycle)', right_index=True)
match_tca_km = pd.merge(how='inner', left=tca_kegg, right=km, left_on='> Citrate cycle (TCA cycle)', right_index=True)
match_tca_kcat = pd.merge(how='inner', left=tca_kegg, right=kcat, left_on='> Citrate cycle (TCA cycle)', right_index=True)

###############################################################################



