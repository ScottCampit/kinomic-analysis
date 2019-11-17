# -*- coding: utf-8 -*-
"""
Data explortation
@author: Scott Campit
"""

import pandas as pd
import numpy as np
import scipy
from os import walk
import matplotlib.pyplot as plt
import seaborn as sns


def getIntersectingGenes(file, modelGeneList):
    """
    getIntersectingGenes outputs a list of the common unique genes between two sets of data
    """
    data = pd.read_csv(file)
    modelGenes = set(modelGeneList)
    dataGenes = set(list(data['Gene names']))
    commonGenes = list(modelGenes.intersection(dataGenes))

    return commonGenes, len(commonGenes)


kineticsData = r'./../Data/sabiodb/homosapien.csv'
model = pd.read_csv(r'./../Data/geneLists/knownEpiFactors.csv')
modelGeneList = list(model.HGNC)
modelGeneList = [gene.strip(' ') for gene in modelGeneList]
commonGenes, numberOfGenes = getIntersectingGenes(kineticsData, modelGeneList)
print(commonGenes)
print(numberOfGenes)


def queryFeatures(kin, param_type=''):
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
        df = df.set_index('Gene names').rename(columns={'parameter.startValue': param_type,
                                                        'parameter.associatedSpecies': 'Substrate'})
        df[param_type] = np.log(df[param_type])*-1
        df = df.pivot_table(
            index='Pathway', columns='Substrate', values=param_type)
    elif param_type == 'S/Km':
        df = df.groupby(
            ['Gene names', 'parameter.associatedSpecies', 'Pathway'])
    else:
        df = df.groupby(['Gene names', 'Pathway'],
                        as_index=False)['parameter.startValue'].mean()
        df = df.rename(columns={'Gene names': 'Gene symbol',
                                'parameter.startValue': param_type})
        df[param_type] = np.log(df[param_type])
        df = df.pivot_table(
            index='Pathway', columns='Gene symbol', values=param_type)

    df = df.T
    df = df.fillna(value=0)

    save_path = r'C:\Users\scampit\Desktop\Kinetic Model\Data\sabiodb\pathway_association'
    if param_type == 'kcat/Km':
        df.to_csv(save_path+'\\'+'substrate_sensor_metabolic.csv')
    else:
        df.to_csv(save_path+'\\'+param_type+'_metabolic.csv')
    return df
