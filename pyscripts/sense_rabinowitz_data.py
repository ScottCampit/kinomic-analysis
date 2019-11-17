""" Rabinowitz et al metabolic data
Need to analyze metabolic sensing in humans

Author: Scott Campit
"""

import pandas as pd

# necessary objects
path = r"C:\Users\scampit\Desktop\Kinetic Model"
fil = r"nchembio.2077-S3.xlsx"

# human dataframe
df = pd.read_excel(path+"\\"+fil)
human = df.loc[df["Organism"] == "Homo sapiens"]

# map EC number to gene symbol
uniq_ecNum = pd.unique(human["EC Number"])
uniq_ecNum = pd.DataFrame(uniq_ecNum)
uniq_ecNum.to_csv("ecNum.csv")

gene_fil = "rabinowitz_results.xlsx"
gene_symbol = pd.read_excel(path+"\\"+gene_fil)

# Map gene symbol back to EC number
human = human.merge(gene_symbol, how="inner", left_on="EC Number", right_on="Input Value")
human = human.drop(columns=["EC Number", "Abbreviation", "Organism", "Input Value", "Input Type"])
human = human.drop_duplicates(keep="first")
human = human[~human["Gene Symbol"].isin(["-"])]

# clean genes
s = human['Gene Symbol'].str.split('; ').apply(pd.Series, 1).stack()
s.index = s.index.droplevel(-1)
s.name = 'Gene Symbol'
del human['Gene Symbol']
# drop duplicates
human = human.join(s)
human = human.drop_duplicates(keep='first')

# s/Km
import matplotlib.pyplot as plt
s_km = human.drop(columns=["Concentration (M)", "Km (M)*", "Metabolite"])
s_km = s_km.set_index("Gene Symbol")
s_km.value_counts().plot(kind='bar')
s_km.plot.hist(ylim=(0,5), bins=500)
s_km.to_excel("S_Km.xlsx")

#from tabula import read_pdf
#pdf = "nchembio.2077-S1.pdf"
#df = read_pdf(path+"\\"+pdf, pages=[2,3])
#header = df.iloc[0]
#df = df.rename(columns=header).drop(df.index[0])
#df = pd.read_csv(path+"\\"+"flux.csv").drop(columns=["Unnamed: 0", "flux", "L.B.",
#                "U.B.", "flux.1", "L.B..1", "U.B..1", "Unnamed: 11"])
df = pd.read_csv(path+"\\"+"flux.csv")
#df["Flux"], df["lb"], df["ub"] = df["flux L.B. U.B."].str.split(' ').str
#df = df.drop(columns="flux L.B. U.B.")
#df.to_csv(path+"\\"+"flux.csv")

# flux values from same paper
net_flux = df.loc[df["Direction"] == "net"]
net_flux = net_flux.drop(columns=["Direction", "Substrates", "Products", "lb", "ub"])
net_flux = net_flux.set_index("Reaction")

# Vmax
kinetic_data_path = r'C:\Users\scampit\Desktop\Kinetic Model\Data\sabiodb'
drop = ['Tissue', 'CellularLocation', 'Enzymename', 'EnzymeType',
        'parameter.standardDeviation', 'Cofactor', 'Inhibitor', 'ChebiID',
        'KeggReactionID']
# kinetics data
kin_fil = r'homosapien_clean.csv'
kin = pd.read_csv(kinetic_data_path+'/'+kin_fil).drop(drop, axis=1)

vmax = kin.loc[kin['parameter.type'] == 'Vmax']
vmax = vmax.groupby('Gene names', as_index=False)['parameter.startValue'].mean()
vmax = vmax.set_index('Gene names').rename(columns={'parameter.startValue':'vmax'})
vmax.to_excel("Vmax_sabiork.xlsx")

#net_flux.to_excel("net_flux.xlsx")
import numpy as np
import seaborn as sns
from scipy.optimize import curve_fit
from sklearn import preprocessing

kinet = pd.read_excel("net_flux.xlsx", sheetname="Sheet2")
kinet = kinet.set_index("Gene")

mms = preprocessing.MinMaxScaler()
x = kinet['Flux'].values.reshape(-1,1)
x = mms.fit_transform(x)
kinet['Norm Flux'] = x

nam = kinet.index
fig = sns.scatterplot(x='S/Km', y='Flux', data=kinet)
plt.title("Substrate sensitivity from Rabinowitz")
plt.xlabel("S/Km")
plt.ylabel("V")
fig = fig.get_figure()
fig.savefig("v_skm.png", dpi=400)
# Number of genes used to sample this:
# 17 mapped to this

# exact match doesn't yield anything
comb = pd.merge(s_km, net_flux, how="inner", left_index=True, right_index=True)
key1 = pd.Series(s_km.index)
key2 = pd.Series(net_flux.index)

# Trying to do fuzzy match with difflib but meh
#from difflib import get_close_matches
#net_flux.index = net_flux.index.map(lambda x: get_close_matches(x, s_km.index)[0])

# So instead, I'll try using fuzzywuzzy -- works better
from fuzzywuzzy import fuzz
from fuzzywuzzy import process

names_array = []
ratio_array = []
def match_names(wrong_names,correct_names):
    for row in wrong_names:
        x=process.extractOne(row, correct_names)
        names_array.append(x[0])
        ratio_array.append(x[1])
    return names_array,ratio_array

name_match, ratio_match = match_names(key2, key1)
results = pd.DataFrame()
results['Name'] = pd.Series(name_match)
results['Ratios'] = pd.Series(ratio_match)
results = results[(results['Ratios'] >= 80)]

match = s_km.index[results.index]



results = results['Name'].unique().tolist()
results = pd.DataFrame(results)
test = s_km.isin(results)