# -*- coding: utf-8 -*-
"""
Created on Tue Aug 28 12:14:42 2018
Cleaning up BRENDA data

@author: Scott Campit
"""

import os

# working directory:
path = r"C:\Users\scampit\Desktop\Kinetic Model\Data\brenda"
os.chdir(path)

import pandas as pd

# open maps
nam2symbol = pd.read_excel("enzymename_genesymbol.xlsx")
nam2symbol = nam2symbol[["Old family ID", "New family name"]]
ec2nam = pd.read_excel("brenda_map.xlsx")

# import brenda files
kcat = pd.read_excel("BRENDA.xlsx", sheet_name = "kcat")
km = pd.read_excel("BRENDA.xlsx", sheet_name = "km")
ki = pd.read_excel("BRENDA.xlsx", sheet_name = "ki")

# try merging
df = pd.merge(kcat, km, how='inner', on = ['EC Number', "Comments", 'PMID', "Substrate", "Ligand ID"])

# Some analysis
simple = df.fillna("no comment")
simple = df[~df['Comments'].str.contains("isozyme", na=False)]
simple = simple[~df['Comments'].str.contains("allozyme", na=False)]
simple = simple[~df['Comments'].str.contains("mutant", na=False)]
simple = simple[~df['Comments'].str.contains("isoenzyme", na=False)]
simple = simple[~df['Comments'].str.contains("mutation", na=False)]
simple = simple[~df['Comments'].str.contains("isoform", na=False)]
simple = simple[~df['Comments'].str.contains("variant", na=False)]

drop = simple.pop('Comments')
drop = simple.pop('PMID')

ligandmap = simple[['Ligand ID', 'Substrate']]
ligandmap = ligandmap.drop_duplicates(subset=['Ligand ID'], keep='first')
ligandmap = ligandmap.sort_values(by = 'Ligand ID')

simple = simple.groupby(['EC Number', "Ligand ID"]).mean()
simple = simple.reset_index()
simple = simple.drop_duplicates(keep='first')
simple = pd.merge(ligandmap, simple, how='inner', on = ['Ligand ID'])
simple = simple.drop('Ligand ID', axis=1)
#simple.to_excel(path+'/'+"brenda_clean.xlsx")

# Cleaned up manually
#simple = pd.read_excel('brenda_clean.xlsx')
#simple = pd.merge(simple, ec2nam, how='inner', on=['EC Number'])

import mygene
mg = mygene.MyGeneInfo()
ec_list = simple["EC Number"].drop_duplicates(keep="first").tolist()
genes = mg.querymany(ec_list, scopes='ec', fields='symbol', species='human', as_dataframe=True)
genes = genes.dropna(subset=['symbol'])
genes = genes['symbol']
genes = genes.reset_index()
test = pd.merge(simple, genes, how='inner', left_on='EC Number', right_on='query')
simple = test.drop(columns = ['EC Number', 'query', 'Enzyme Name'])
#simple = simple.to_excel(path+'/'+'brenda_clean_mapped.xlsx')
simple = simple.groupby(['symbol', "Substrate"]).mean()
simple = simple.reset_index()
#simple = simple.to_excel(path+'/'+'brenda_clean_mapped.xlsx')

import numpy as np
import scipy
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.datasets import make_blobs

sns.clustermap(simple, method='centroid')
