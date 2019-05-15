# -*- coding: utf-8 -*-
"""
Created on Wed Aug 22 14:33:35 2018
Evolutionary comparison of kinetic parameters

@author: scampit
"""

import pandas as pd
import numpy as np
import scipy, os

data_path = r'C:\Users\scampit\Desktop\Kinetic Model\Data\sabiodb\SabioRK_Mammals'

### First: Mammalian -- human, mouse, rat, boar, cow. Need to map UNIPROT ID to gene symbol ###
drop_col = ['Enzymename', 'EnzymeType', 'EntryID', 'Organism', 'Tissue', 'Pathway', 
            'CellularLocation', 'ECNumber', 'Cofactor', 'Inhibitor', 'ChebiID',
            'pH', 'Temperature', 'KeggReactionID', 'ReactionEquation', 
            'KineticMechanismType','PubMedID']

human = pd.read_csv(data_path+'\\'+'homosapien.txt', sep='\t').drop(drop_col, axis=1)
mouse = pd.read_csv(data_path+'\\'+'musmusculus.txt', sep='\t').drop(drop_col, axis=1)
rat = pd.read_csv(data_path+'\\'+'rattusnorvegicus.txt', sep='\t').drop(drop_col, axis=1)
boar = pd.read_csv(data_path+'\\'+'susscrofa.txt', sep='\t').drop(drop_col, axis=1)
cow = pd.read_csv(data_path+'\\'+'bostaurus.txt', sep='\t').drop(drop_col, axis=1)
yeast = pd.read_csv(data_path+'/'+'yeast.txt', sep='\t').drop(drop_col, axis=1)

human_uniprot = human['UniprotID'].drop_duplicates(keep='first') 
mouse_uniprot = mouse['UniprotID'].drop_duplicates(keep='first')
rat_uniprot = rat['UniprotID'].drop_duplicates(keep='first')
boar_uniprot = boar['UniprotID'].drop_duplicates(keep='first')
cow_uniprot = cow['UniprotID'].drop_duplicates(keep='first')
yeast_uniprot = yeast['UniprotID'].drop_duplicates(keep='first') 

###############################################################################
def reshaper(series, column='', sep='', name=''):
    '''
    reshaper turns lists within an element of a pd.Series datatype and returns
    a linear pd.Series
    '''
    s = series.str.split(sep).apply(pd.Series, 1).stack()
    s.index = s.index.droplevel(-1)
    s.name = name
    s = s.drop_duplicates(keep='first')
    return s
###############################################################################

human_uniprot = reshaper(human_uniprot, column='UniprotID', sep=' ', name = 'UniprotID')
mouse_uniprot = reshaper(mouse_uniprot, column='UniprotID', sep=' ', name = 'UniprotID')
rat_uniprot = reshaper(rat_uniprot, column='UniprotID', sep=' ', name = 'UniprotID')
boar_uniprot = reshaper(boar_uniprot, column='UniprotID', sep=' ', name = 'UniprotID')
cow_uniprot = reshaper(cow_uniprot, column='UniprotID', sep=' ', name = 'UniprotID')
yeast_uniprot = reshaper(yeast_uniprot, column='UniprotID', sep=' ', name = 'UniprotID')

# save uniprot ids (may be useful later? idk)
save_path = r'C:\Users\scampit\Desktop\Kinetic Model\Data\sabiodb\uniprot_ids'
human_uniprot.to_csv(save_path+'\\'+'human_uniprot.csv', index=False)
mouse_uniprot.to_csv(save_path+'\\'+'mouse_uniprot.csv', index=False)
rat_uniprot.to_csv(save_path+'\\'+'rat_uniprot.csv', index=False)
boar_uniprot.to_csv(save_path+'\\'+'boar_uniprot.csv', index=False)
cow_uniprot.to_csv(save_path+'\\'+'cow_uniprot.csv', index=False)
yeast_uniprot.to_csv(save_path+'\\'+'yeast_uniprot.csv', index=False)

########## need to work on API method of retrieving gene ids ##################
#def gene_retriever(list_of_genes, df)
#    '''
#    gene_retriever gets the UNIPROT accession number and uses their API to obtain
#    the gene symbol. Then it replaces the uniprot accession ID with the gene
#    symbol within the dataframe of interest, while turning the gene symbol into
#    the dataframe index.
#    '''
    
#    import urllib, urllib.request
#    url = 'https://www.uniprot.org/uploadlists/'
#    params = {
#            'from': 'ACC',
#            'to': 'GENENAME',
#            'format': 'csv',
#            'query': list_of_genes.tolist()
#            }
    
#    data = urllib.parse.urlencode(params)
#    requested = urllib.request.Request(url, data)
#    contact="scampit@umich.edu"
#    request.add_header('User-Agent', 'Python %s' % contact)
#    response = urllib.request.urlopen(requested)
#    page = response.read(200000)
#    return 'All Good!'
###############################################################################

# Now we need a function that will clean up our dataframes:
def uniprot_to_gene_mapper(organism, df):
    gene_dir = r'C:\Users\scampit\Desktop\Kinetic Model\sabiodb\genenames'
    gene_df = pd.read_csv(gene_dir+'\\'+organism+'.txt', sep='\t')
    df = pd.merge(df, gene_df, how='inner', left_on='UniprotID', right_on='From')
    df = df.drop(labels=['From', 'UniprotID'], axis=1)
    df['To'] = df.loc[:, 'To'].apply(lambda x: x.upper())
    df = df.set_index('To')
    save_dir = r'C:\Users\scampit\Desktop\Kinetic Model\sabiodb\clean_df'
    df.to_csv(save_dir+'\\'+organism+'.csv')
    return df

###############################################################################

# clean up dataframe (automatically saved)    
human = uniprot_to_gene_mapper('human', human)
mouse = uniprot_to_gene_mapper('mouse', mouse)
rat = uniprot_to_gene_mapper('rat', rat)
boar = uniprot_to_gene_mapper('boar', boar)
cow = uniprot_to_gene_mapper('cow', cow)
yeast = uniprot_to_gene_mapper('yeast', yeast)

# Okay, we'll now look at each organism and compare how many glycolytic genes each one has
#organism = ['human', 'rat', 'mouse', 'boar', 'cow', 'yeast']
def glycolytic_match(organism, df):
    #for o in organism:
    glycolysis = ['HK1', 'HK2', 'HK3', 'GCK', 'GPI', 'PGI', 'PHI', 'PFKM', 
                  'PFKL', 'PFKP', 'ALDOA', 'ALDOB', 'ALDOC', 'GAPDH', 'PGK1', 
                  'PGK2', 'PGAM1', 'PGAM2', 'PGAM4', 'ENO1', 'ENO2', 'ENO3', 
                  'PKLR', 'PKM2', 'ACSS1', 'ACSS2', 'TPI1']
    glycolysis = pd.DataFrame({'Glycolysis':glycolysis})
    kcat = df.loc[df['parameter.type'] == 'kcat'].drop_duplicates(keep='first')
    km = df.loc[df['parameter.type'] == 'km'].drop_duplicates(keep='first')
    kcat_km = df.loc[df['parameter.type'] == 'kcat/Km'].drop_duplicates(keep='first')
    match = pd.merge(how='inner', left=kcat, right=glycolysis, left_index=True, right_on='Glycolysis')
    match = match.drop(labels = ['parameter.type', 'parameter.associatedSpecies', 
                                 'parameter.endValue', 'parameter.standardDeviation', 
                                 'parameter.unit'], axis=1).rename(columns={'parameter.startValue':'kcat'})
    match = match.groupby(by='Glycolysis')['kcat'].mean()
    #match = match.drop_duplicates(subset='kcat', keep='first').set_index('Glycolysis')
        
    return match

human_gly = pd.DataFrame(glycolytic_match('human', human)).rename(columns={'kcat':'human'})
rat_gly = pd.DataFrame(glycolytic_match('rat', rat)).rename(columns={'kcat':'rat'})
mouse_gly = pd.DataFrame(glycolytic_match('mouse', mouse)).rename(columns={'kcat':'mouse'})
boar_gly = pd.DataFrame(glycolytic_match('boar', boar)).rename(columns={'kcat':'boar'})
cow_gly = pd.DataFrame(glycolytic_match('cow', cow)).rename(columns={'kcat':'cow'})
yeast_gly = pd.DataFrame(glycolytic_match('yeast', yeast)).rename(columns={'kcat':'yeast'})

kcat = pd.concat([human_gly,rat_gly, mouse_gly, boar_gly, cow_gly, yeast_gly], 
                 sort=True, axis=1)
###############################################################################





































