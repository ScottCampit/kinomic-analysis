# -*- coding: utf-8 -*-
"""
`parse-sabiork.py` contains several functions that enables the analysis of the SABIO-RK database.

`construct_db` reads in the text files obtained from the SABIO-RK API and cleans the data up into a single db.

@author: Scott Edward Campit
"""

import os
import mygene
import pandas as pd

dfs = []
def construct_db(fil, spec):
    """
    construct_df will make a dataframe from the .txt files for data visualizations.
    """

    col_to_drop = [
        "EntryID",
        "Tissue",
        "CellularLocation",
        "Enzymename",
        "parameter.standardDeviation",
        "parameter.unit",
        "Product",
        "Cofactor",
        "Inhibitor",
        "ChebiID",
        "KeggReactionID",
        "ReactionEquation",
        "PubMedID",
        "KineticMechanismType"
    ]

    # Read in and clean stuff up
    input = r'/mnt/c/Users/scampit/Desktop/kinomics/Data/sabiodb/'
    df = pd.read_csv(input+fil, sep='\t')
    df = df.drop(col_to_drop, axis=1)
    print(df)

    # Get UniProt IDs -> Gene Symbol
    uniprot = df['UniprotID'].unique()
    mg = mygene.MyGeneInfo()
    query = mg.querymany(uniprot,
        scopes='uniprot',
        fields='symbol',
        species=spec,
        as_dataframe=True
    )
    query = query.loc[query['query', 'symbol']]
    print(query)
    df = pd.merge(df, query, how='inner', left_on='UniprotID', right_on='query')
    df = df.set_index('symbol')

    # Unstack df by substrate
    tmp = df['Substrate'].str.split(';').apply(pd.Series, 1).stack()
    tmp.index = tmp.index.droplevel(-1)
    tmp.name = 'Substrate'
    del df['Substrate']
    df = df.join(tmp)
    df = df.drop_duplicates(keep='first')
    return df

# Compile into database
input = r'/mnt/c/Users/scampit/Desktop/kinomics/Data/sabiodb/'

organisms = pd.read_json(input+'organisms.json', typ='series').to_frame('ID')

for fil in os.listdir(input):
    if fil.endswith('.txt'):
        for idx, spec in organisms['ID'].iteritems():
            if idx == fil:
                nam = spec
                print(nam)
                df = construct_db(fil, nam)
                dfs.append(df)
#db = pd.concat(dfs, axis=1)
#path_to_save = r'/mnt/c/Users/scampit/Desktop/kinomics/Data/sabiodb/'
#db.to_json(path_to_save+'sabiorkdb.json')
