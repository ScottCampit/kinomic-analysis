"""
Get histone reader/writer/eraser kinetomics data sabioRK

@author: Scott Campit
"""

import pandas as pd
import numpy as np


def getGeneSymbols(data, ID_to_transform, output_ID='symbol', species='Human'):
    """
    """
    import mygene

    mg = mygene.MyGeneInfo()
    geneName = list(data.values.flatten())

    geneMap = mg.querymany(geneName,
                           scopes=ID_to_transform,
                           fields=output_ID, species=species, as_dataframe=True)

    MappedSet = pd.merge(data, geneMap,
                         how='inner',
                         left_index=True, right_index=True)
    MappedSet = MappedSet.set_index('symbol')

    return MappedSet


def subsetDF(data, list_of_subset_genes):
    """
    """
    subset = data[data['Gene names'].isin(list_of_subset_genes)]
    return subset


data = pd.read_csv('./../Data/sabiodb/homosapien.csv')
epifactor = pd.read_csv('./../Data/geneLists/knownEpiFactors.csv')
geneList = list(epifactor['HGNC'].values)

filtered = subsetDF(data, geneList)
filtered.to_csv('epigeneticsSabioRK.csv')
