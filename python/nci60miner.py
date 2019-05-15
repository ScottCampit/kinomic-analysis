"""
Objective of this script: to dynamically mine most databases for the NCI-60 cell line data periodically (probably evey month or so).

Eventually, I want it to be able to classify the datasets as well into multiple sub directories, but let's do this one step at a time.
"""

import pandas as pd
import matplotlib.pyplot as pyplot
import numpy as np
import scipy.stats as stats
import seaborn as sns
from matplotlib import rcParams

def miner(path, rna_fil, methyl_fil, mirna_fil, genexp_fil, prot_fil):

    # Import files as dataframes
    rna = pd.read_excel(path+"\\"+"NCI60_Summary.xlsx", sheet_name="RNA_averaged", skiprows=2)
    me = pd.read_excel(path+"\\"+"NCI60_Summary.xlsx", sheet_name="Methylation_averaged", skiprows=2)
    rna = pd.read_excel(path+"\\"+"NCI60_Summary.xlsx", sheet_name="genexp_averaged", skiprows=2)
    me = pd.read_excel(path+"\\"+"NCI60_Summary.xlsx", sheet_name="mirna_averaged", skiprows=2)
    prot = pd.read_excel(path+"\\"+"NCI60_Summary.xlsx", sheet_name="prot_averaged", skiprows=2)

    raw_mat =
