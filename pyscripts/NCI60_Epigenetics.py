# -*- coding: utf-8 -*-
"""
NCI-60 Data Exploration for Epigenetic Modifications

@author: scampit
"""

import pandas as pd
path = r"C:\Users\scampit\Desktop\Kinetic Model\Data\NCI60"
# Import datasets
rna = pd.read_excel(path+"\\"+"NCI60_Summary.xlsx", sheet_name="RNA_averaged",\
                    skiprows=2)
me = pd.read_excel(path+"\\"+"NCI60_Summary.xlsx", sheet_name="Methylation_averaged", \
                   skiprows=2)

# Process datasets
rna = rna.set_index("Row Labels").drop(columns=['Average of ME:MDA-N'], axis=1)
me = me.set_index("Row Labels").drop(columns=['Average of ME:MDA-N'], axis=1)

# Create merged dataframe
df = pd.merge(rna, me, how='inner', left_index=True, right_index=True, \
              suffixes = ('rna', 'me')).dropna()

# I'm choosing three specific cell lines due to the available phosphoproteomics
# data

MCF7rna = pd.DataFrame(df.iloc[:,0])
MCF7me = pd.DataFrame(df.iloc[:,59])

A549rna = pd.DataFrame(df.iloc[:,33])
A549me = pd.DataFrame(df.iloc[:,92])

PC3rna = pd.DataFrame(df.iloc[:,49])
PC3me = pd.DataFrame(df.iloc[:,108])

# Construct raw matricies
MCF7 = pd.concat([MCF7rna, MCF7me], axis=1).round(2)
A549 = pd.concat([A549rna, A549me], axis=1).round(2)
PC3 = pd.concat([PC3rna, PC3me], axis=1).round(2)

# The data is skewed with a tail to the right. I tried to do simple min_max scaling
#from sklearn import preprocessing
#import numpy as np
#min_max_scaler = preprocessing.MinMaxScaler()
#x_scaled = min_max_scaler.fit_transform(x)
#y_scaled = min_max_scaler.fit_transform(y)
#df_scaled1 = pd.DataFrame({"RNA Expression":x_scaled[:,0], "Methylation":y_scaled[:,0]},\
#                         index=df.index.values.tolist(), columns=["RNA Expression",
#                                                     "Methylation"])

#test_scaled = min_max_scaler.fit_transform(test)
#df_scaled2 = pd.DataFrame(test_scaled, index=df.index.values.tolist(),\
#                          columns=["RNA Expression", "Methylation"])

# Unfortunately the data still doesn't have a normal distribution. So I'll
# try a log transform:

from sklearn.preprocessing import FunctionTransformer
import numpy as np
transformer = FunctionTransformer(np.log1p, validate=True)

MCF7_log = transformer.transform(MCF7)
MCF7_log = pd.DataFrame(MCF7_log, index=df.index.values.tolist(),\
                          columns=["RNA Expression", "Methylation"])
A549_log = transformer.transform(A549)
A549_log = pd.DataFrame(A549_log, index=df.index.values.tolist(),\
                          columns=["RNA Expression", "Methylation"])
PC3_log = transformer.transform(PC3)
PC3_log = pd.DataFrame(PC3_log, index=df.index.values.tolist(),\
                          columns=["RNA Expression", "Methylation"])

# Different plots
import matplotlib.pyplot as plt
import seaborn as sns

# Scatter
#with sns.axes_style('white'):
#    sns.scatterplot(data=df_scaled)

# Regression
#with sns.axes_style('white'):
#    sns.jointplot(x="RNA Expression", y="Methylation", data=df_scaled, kind="reg")

# Hex plots to visualize some of the data
with sns.axes_style('white'):
    sns.jointplot(x="RNA Expression", y="Methylation", data=MCF7_log, kind="hex")
    #plt.savefig("MCF7.png", dpi=400)
with sns.axes_style('white'):
    sns.jointplot(x="RNA Expression", y="Methylation", data=A549_log, kind="hex")
    #plt.savefig("A549.png", dpi=400)
with sns.axes_style('white'):
    sns.jointplot(x="RNA Expression", y="Methylation", data=PC3_log, kind="hex")
    #plt.savefig("PC3.png", dpi=400)\

###############################################################################
# Okay, now that we did a global analysis, let's look at specific pathways,
    # in particular Methionine metabolism and TCA
###############################################################################

########################### Methionine metabolism #############################
import pandas as pd

kegg_annot = r"C:\Users\scampit\Desktop\Kinetic Model\Data\KEGG_Metabolism"
met = r"KEGG_CYSTEINE_AND_METHIONINE_METABOLISM.txt"

met_genes = pd.read_csv(kegg_annot+"\\"+met)
met = pd.merge(met_genes, df, how="inner", left_on="GENE", right_index=True)
met = met.set_index("GENE")

# RNA and Methylation patterns for all 59 cell lines
met_rna = pd.merge(met_genes, rna, how="inner", left_on="GENE", \
                   right_index=True).set_index("GENE")


met_methyl = pd.merge(met_genes, me, how="inner", left_on="GENE", \
                   right_index=True).set_index("GENE")

fill_val = met_methyl.mean()
met_methyl = met_methyl.fillna(value=fill_val)
#fill_val = fill_val.mean()


# Get matching set for comparisons
match = met_rna.index.intersection(met_methyl.index)
match = pd.DataFrame(match).set_index("GENE")
met_rna = pd.merge(match, met_rna, how="inner", left_index=True, right_index=True)
met_methyl = pd.merge(match, met_methyl, how="inner", left_index=True, right_index=True)

# Normalize data
from sklearn.preprocessing import FunctionTransformer
import numpy as np

transformer = FunctionTransformer(np.log1p, validate=True)

met_rna_t = transformer.transform(met_rna)
met_rna_t = pd.DataFrame(met_rna_t, index=met_rna.index.values.tolist(),\
                         columns=met_rna.columns.values.tolist())
met_methyl_t = transformer.transform(met_methyl)
met_methyl_t = pd.DataFrame(met_methyl_t, index=met_methyl.index.values.tolist(),\
                         columns=met_methyl.columns.values.tolist())

# Box plots of all NCI60 cell line data for Methionine metabolism
import matplotlib.pyplot as plt
import seaborn as sns

plt.xticks(rotation=90)
sns.boxplot(data=met_rna_t)

plt.xticks(rotation=90)
sns.boxplot(data=met_methyl_t)
plt.title("Methionine metabolism methylation and RNA expression")
plt.savefig("Methionine metabolism expression.png", dpi=400)

# Calculate correlation matrix and output correlation figure
#met_corr = met_rna_t.corrwith(met_methyl_t, axis=1)
#sns.heatmap(met_corr)

# Scatterplot
met_concat = pd.concat([met_rna_t.assign(dataset="RNA Expression"), \
                         met_methyl_t.assign(dataset="Methylation")])

#sns.scatterplot(x=met_rna_t, y =met_methyl_t, data=met_concat, style='dataset')

###############################################################################
# RNA and Methylation patterns for three cell lines with Methionine metabolism
MCF7rna = pd.DataFrame(met_rna.iloc[:,0])
MCF7me = pd.DataFrame(met_methyl.iloc[:,0])

A549rna = pd.DataFrame(met_rna.iloc[:,33])
A549me = pd.DataFrame(met_methyl.iloc[:,33])

PC3rna = pd.DataFrame(met_rna.iloc[:,49])
PC3me = pd.DataFrame(met_methyl.iloc[:,49])

MCF7 = pd.concat([MCF7rna, MCF7me], axis=1).round(2)
A549 = pd.concat([A549rna, A549me], axis=1).round(2)
PC3 = pd.concat([PC3rna, PC3me], axis=1).round(2)

from sklearn.preprocessing import FunctionTransformer
import numpy as np
transformer = FunctionTransformer(np.log1p, validate=True)

MCF7_log = transformer.transform(MCF7)
MCF7_log = pd.DataFrame(MCF7_log, index=met_rna.index.values.tolist(),\
                          columns=["RNA Expression", "Methylation"])
A549_log = transformer.transform(A549)
A549_log = pd.DataFrame(A549_log, index=met_rna.index.values.tolist(),\
                          columns=["RNA Expression", "Methylation"])
PC3_log = transformer.transform(PC3)
PC3_log = pd.DataFrame(PC3_log, index=met_rna.index.values.tolist(),\
                          columns=["RNA Expression", "Methylation"])

# Different plots for gene expression over methylation in Methionine Metabolism
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="white")
# Hex Plots
#with sns.axes_style('white'):
#    sns.jointplot(x="RNA Expression", y="Methylation", data=MCF7_log, kind="hex")
#    plt.savefig("MCF7.png", dpi=400)
#with sns.axes_style('white'):
#    sns.jointplot(x="RNA Expression", y="Methylation", data=A549_log, kind="hex")
#    plt.savefig("A549.png", dpi=400)
#with sns.axes_style('white'):
#    sns.jointplot(x="RNA Expression", y="Methylation", data=PC3_log, kind="hex")
#    plt.savefig("PC3.png", dpi=400)

# Regression plots
sns.jointplot(x="RNA Expression", y="Methylation", data=MCF7_log, kind="reg")
sns.jointplot(x="RNA Expression", y="Methylation", data=A549_log, kind="reg")
sns.jointplot(x="RNA Expression", y="Methylation", data=PC3_log, kind="reg")

################################ TCA Cycle ####################################
import pandas as pd

kegg_annot = r"C:\Users\scampit\Desktop\Kinetic Model\Data\KEGG_Metabolism"
tca = r"KEGG_CITRATE_CYCLE_TCA_CYCLE.txt"

tca_genes = pd.read_csv(kegg_annot+"\\"+tca)
tca = pd.merge(tca_genes, df, how="inner", left_on="GENE", right_index=True)
tca = tca.set_index("GENE")

# RNA and Methylation patterns for all 59 cell lines
tca_rna = pd.merge(tca_genes, rna, how="inner", left_on="GENE", \
                   right_index=True).set_index("GENE")

fill_val = tca_rna.mean()
fill_val = fill_val.mean()
tca_methyl = pd.merge(tca_genes, me, how="inner", left_on="GENE", \
                   right_index=True).set_index("GENE").fillna(value=fill_val)

# Get matching set for comparisons
match = tca_rna.index.intersection(tca_methyl.index)
match = pd.DataFrame(match).set_index("GENE")
tca_rna = pd.merge(match, tca_rna, how="inner", left_index=True, right_index=True)
tca_methyl = pd.merge(match, tca_methyl, how="inner", left_index=True, right_index=True)

# Normalize data
from sklearn.preprocessing import FunctionTransformer
import numpy as np

transformer = FunctionTransformer(np.log1p, validate=True)

tca_rna_t = transformer.transform(tca_rna)
tca_rna_t = pd.DataFrame(tca_rna_t, index=tca_rna.index.values.tolist(),\
                         columns=tca_rna.columns.values.tolist())
tca_methyl_t = transformer.transform(tca_methyl)
tca_methyl_t = pd.DataFrame(tca_methyl_t, index=tca_methyl.index.values.tolist(),\
                         columns=tca_methyl.columns.values.tolist())

# Box plots of all NCI60 cell line data for Methionine metabolism
import matplotlib.pyplot as plt
import seaborn as sns

plt.xticks(rotation=90)
sns.boxplot(data=tca_rna_t)

plt.xticks(rotation=90)
sns.boxplot(data=tca_methyl_t)

# Calculate correlation matrix and output correlation figure
#met_corr = met_rna_t.corrwith(met_methyl_t, axis=1)
#sns.heatmap(met_corr)

# Scatterplot
tca_concat = pd.concat([tca_rna_t.assign(dataset="RNA Expression"), \
                         tca_methyl_t.assign(dataset="Methylation")])

#sns.scatterplot(x=met_rna_t, y =met_methyl_t, data=met_concat, style='dataset')

###############################################################################
# RNA and Methylation patterns for three cell lines with Methionine metabolism
MCF7rna = pd.DataFrame(tca_rna.iloc[:,0])
MCF7me = pd.DataFrame(tca_methyl.iloc[:,0])

A549rna = pd.DataFrame(tca_rna.iloc[:,33])
A549me = pd.DataFrame(tca_methyl.iloc[:,33])

PC3rna = pd.DataFrame(tca_rna.iloc[:,49])
PC3me = pd.DataFrame(tca_methyl.iloc[:,49])

MCF7 = pd.concat([MCF7rna, MCF7me], axis=1).round(2)
A549 = pd.concat([A549rna, A549me], axis=1).round(2)
PC3 = pd.concat([PC3rna, PC3me], axis=1).round(2)

from sklearn.preprocessing import FunctionTransformer
import numpy as np
transformer = FunctionTransformer(np.log1p, validate=True)

MCF7_log = transformer.transform(MCF7)
MCF7_log = pd.DataFrame(MCF7_log, index=tca_rna.index.values.tolist(),\
                          columns=["RNA Expression", "Methylation"])
A549_log = transformer.transform(A549)
A549_log = pd.DataFrame(A549_log, index=tca_rna.index.values.tolist(),\
                          columns=["RNA Expression", "Methylation"])
PC3_log = transformer.transform(PC3)
PC3_log = pd.DataFrame(PC3_log, index=tca_rna.index.values.tolist(),\
                          columns=["RNA Expression", "Methylation"])

# Different plots for gene expression over methylation in Methionine Metabolism
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="white")
# Hex Plots
#with sns.axes_style('white'):
#    sns.jointplot(x="RNA Expression", y="Methylation", data=MCF7_log, kind="hex")
#    plt.savefig("MCF7.png", dpi=400)
#with sns.axes_style('white'):
#    sns.jointplot(x="RNA Expression", y="Methylation", data=A549_log, kind="hex")
#    plt.savefig("A549.png", dpi=400)
#with sns.axes_style('white'):
#    sns.jointplot(x="RNA Expression", y="Methylation", data=PC3_log, kind="hex")
#    plt.savefig("PC3.png", dpi=400)

# Regression plots
sns.jointplot(x="RNA Expression", y="Methylation", data=MCF7_log, kind="reg")
sns.jointplot(x="RNA Expression", y="Methylation", data=A549_log, kind="reg")
sns.jointplot(x="RNA Expression", y="Methylation", data=PC3_log, kind="reg")

###############################################################################
# Kinetics with gene expression in Methionine metabolism and TCA cycle
###############################################################################
import pandas as pd
path = r"C:\Users\scampit\Desktop\Kinetic Model\Data\brenda"
kin = r"kinetic_map.xlsx"

# First the BRENDA database
brenda = pd.read_excel(path+"\\"+kin, sheet_name = "BRENDA")
kcat_b = brenda.groupby(["symbol", "Substrate"], axis=0)["kcat"].max().reset_index().set_index("symbol")
km_b = brenda.groupby(["symbol", "Substrate"], axis=0)["Km"].max().reset_index().set_index("symbol")

kegg_annot = r"C:\Users\scampit\Desktop\Kinetic Model\Data\KEGG_Metabolism"

###############################################################################
# Now Methionine metabolism
met = r"KEGG_CYSTEINE_AND_METHIONINE_METABOLISM.txt"
met_genes = pd.read_csv(kegg_annot+"\\"+met)

met_kcat_b = pd.merge(met_genes, kcat_b, how="inner", left_on="GENE", right_index=True)
met_kcat_b = met_kcat_b.set_index(["GENE", "Substrate"])

met_km_b = pd.merge(met_genes, km_b, how="inner", left_on="GENE", right_index=True)
met_km_b = met_km_b.set_index(["GENE", "Substrate"])

# SAM and SAH consumption
met_concat_b = pd.merge(met_kcat_b, met_km_b, how="inner", left_index=True, right_index=True)
met_concat_b = met_concat_b.reset_index()
met_concat_b.to_excel("Methionine BRENDA data.xlsx")

sah_b = met_concat_b.loc[met_concat_b["Substrate"] == "S-adenosyl-L-homocysteine"]
sam_b = met_concat_b.loc[met_concat_b["Substrate"] == "S-adenosyl-L-methionine"]

###############################################################################
# TCA Cycle
tca = r"KEGG_CITRATE_CYCLE_TCA_CYCLE.txt"
tca_genes = pd.read_csv(kegg_annot+"\\"+tca)

tca_kcat_b = pd.merge(tca_genes, kcat_b, how="inner", left_on="GENE", right_index=True)
tca_kcat_b = tca_kcat_b.set_index(["GENE", "Substrate"])

tca_km_b = pd.merge(tca_genes, km_b, how="inner", left_on="GENE", right_index=True)
tca_km_b = tca_km_b.set_index(["GENE", "Substrate"])

# SAM and SAH consumption
tca_concat_b = pd.merge(tca_kcat_b, tca_km_b, how="inner", left_index=True, right_index=True)
tca_concat_b = tca_concat_b.reset_index()
tca_concat_b.to_excel("TCA Cycle BRENDA data.xlsx")

#sah_b = met_concat_b.loc[met_concat_b["Substrate"] == "S-adenosyl-L-homocysteine"]
#sam_b = met_concat_b.loc[met_concat_b["Substrate"] == "S-adenosyl-L-methionine"]

###############################################################################
# Unfortunately acetyl CoA was not in the KEGG Annotation for TCA Cycle. Let's
# try glycolysis / gluconeogenesis

gly = r"KEGG_GLYCOLYSIS_GLUCONEOGENESIS.txt"
gly_genes = pd.read_csv(kegg_annot+"\\"+gly)

gly_kcat_b = pd.merge(gly_genes, kcat_b, how="inner", left_on="GENE", right_index=True)
gly_kcat_b = gly_kcat_b.set_index(["GENE", "Substrate"])

gly_km_b = pd.merge(gly_genes, km_b, how="inner", left_on="GENE", right_index=True)
gly_km_b = gly_km_b.set_index(["GENE", "Substrate"])

# SAM and SAH consumption
gly_concat_b = pd.merge(gly_kcat_b, gly_km_b, how="inner", left_index=True, right_index=True)
gly_concat_b = gly_concat_b.reset_index()
gly_concat_b.to_excel("Glycolysis Gluconeogenesis BRENDA data.xlsx")

# Unfortunately acetyl CoA was not in the KEGG Annotation for glycolysis or
# gluconeogenesis. I don't think we can study metabolic sensing for these substrates
# with the current data.

###############################################################################



###############################################################################
# Unfortunately, the enzyme that produces SAM is not in the BRENDA dataset.
# Time to look in the SabioRK
###############################################################################

#xl = pd.ExcelFile(path+'\\'+kin)
#xl.sheet_names

sabio = pd.read_excel(path+"\\"+kin, sheet_name = "SABIORK")

s_params = sabio.groupby(["Gene names", "Substrate"], axis=0)["parameter.type",
                      "parameter.associatedSpecies",
                      "parameter.startValue"].max().reset_index().set_index("Gene names")

kcat_s = s_params.loc[s_params["parameter.type"]=="kcat"]
s = kcat_s['Substrate'].str.split(';').apply(pd.Series, 1).stack()
s.index = s.index.droplevel(-1)
s.name = 'Substrate'
del kcat_s['Substrate']
# drop duplicates
kcat_s = kcat_s.join(s)
kcat_s = kcat_s.drop_duplicates(keep='first')

km_s = s_params.loc[s_params["parameter.type"]== "Km"]
s = km_s['Substrate'].str.split(';').apply(pd.Series, 1).stack()
s.index = s.index.droplevel(-1)
s.name = 'Substrate'
del km_s['Substrate']
# drop duplicates
km_s = km_s.join(s)
km_s = km_s.drop_duplicates(keep='first')

met_kcat_s = pd.merge(met_genes, kcat_s, how="inner", left_on="GENE", right_index=True)
met_kcat_s = met_kcat_s.set_index(["GENE", "Substrate"])

met_km_s = pd.merge(met_genes, km_s, how="inner", left_on="GENE", right_index=True)
met_km_s = met_km_s.set_index(["GENE", "Substrate"])

# SAM and SAH consumption
met_concat_s = pd.merge(met_kcat_s, met_km_s, how="inner", left_index=True, right_index=True)
met_concat_s = met_concat_s.reset_index()
met_concat_s.to_excel("Kinetics_Metabolic_Pathway.xlsx", sheet_name="SABIORK")

sah_s = met_concat_s.loc[met_concat_s["Substrate"] == "S-adenosyl-L-homocysteine"]
sam_s = met_concat_s.loc[met_concat_s["Substrate"] == "S-adenosyl-L-methionine"]
