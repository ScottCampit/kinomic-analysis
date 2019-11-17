#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 10 18:01:13 2018
Selected KEGG pathway analysis for BRENDA and Sabio RK kinetic parameters

Path for KEGG Annotations: C\Users\scampit\Desktop\Kinetic Model\Data\KEGG_Metabolism
Path for kinetic data: C/Users/Desktop/Kinetic Model/Data/brenda/human

@author: scampit
"""

###############################################################################
########################### BRENDA Pathway Analysis ###########################
###############################################################################
import pandas as pd
import numpy as np

brenda = pd.read_excel("kinetic_map.xlsx", sheet_name="BRENDA")
sabiork = pd.read_excel("kinetic_map.xlsx", sheet_name="SABIORK")

<<<<<<< Updated upstream
# Clean Sabio-RK Data to have all substrates as a single row
s = sabiork['Substrate'].str.split(';').apply(pd.Series, 1).stack()
s.index = s.index.droplevel(-1)
s.name = 'Substrate'
del sabiork['Substrate']
sabiork = sabiork.join(s)
sabiork = sabiork.drop_duplicates(keep='first')

# There are 3532 total substrates
substrates = [brenda["Substrate"].drop_duplicates(keep='first'), sabiork["Substrate"].drop_duplicates(keep='first')]
substrates = pd.concat(substrates)
substrates = substrates.drop_duplicates(keep='first').dropna()

# Substrates of interest
query = ["ATP", "ADP", "NAD+", "NADH", "NADP+", "NADPH", \
         "S-adenosyl-DL-homocysteine", "S-adenosyl-L-homocysteine", \
         "S-adenosyl-L-methionine", "Acetyl-CoA", "Coenzyme A", "FAD"]
=======
Substrates = [brenda[''].nunique() + sabiork[''].nunique()]

# BRENDA SUBSTRATES
atp = brenda[brenda["Substrate"]=="ATP"]
adp = brenda[brenda["Substrate"]=="ADP"]
nad = brenda[brenda["Substrate"]=="NAD+"]
nadh = brenda[brenda["Substrate"]=="NADH"]
>>>>>>> Stashed changes

###############################################################################
############################# SABORK SUBSTRATES ###############################
###############################################################################
import matplotlib.pyplot as plt
import seaborn as sns

# Separate by pathway
paths = sabiork[~sabiork["Pathway"].isna()]
pathways = sabiork["Pathway"].unique().tolist()

################################## ATP ########################################
atp = paths[paths["Species"] == "ATP"]
atp_km = atp[atp["Type.1"] == "Km"]
atp_kcatkm = atp[atp["Type.1"] == "kcat/Km"]

atp_kcat = paths[paths["Type.1"] == "kcat"]
s = atp_kcat['Substrate'].str.split(';').apply(pd.Series, 1).stack()
s.index = s.index.droplevel(-1)
s.name = 'Substrate'
del atp_kcat['Substrate']
atp_kcat = atp_kcat.join(s)
atp_kcat = atp_kcat.drop_duplicates(keep='first')

atp_kcat = atp_kcat[atp_kcat["Substrate"] == "ATP"]
atp_kcat = atp_kcat[atp_kcat["Type.1"] == "kcat"]

atp_kcatkm = atp_kcatkm[atp_kcatkm["Substrate"] == "ATP"]
atp_kcatkm = atp_kcatkm[atp_kcatkm["Type.1"] == "kcat/Km"]

# Purine Metabolism
purine_atp_km = atp_km[atp_km["Pathway"] == "Purine metabolism"]
purine_atp_km["Start Value"] = -1*np.log(purine_atp_km["Start Value"])
sns.swarmplot(x = "Gene names", y = "Start Value", data=purine_atp_km)
sns.violinplot(x = "Gene names", y = "Start Value", data=purine_atp_km)
plt.xticks(rotation=90)
plt.title("Substrate sensing for ATP in Purine metabolism")
plt.ylabel("- log(Km) ")
plt.xlabel("Gene names")
plt.savefig("Purinemetabolism_ATP_Kmvar")

# Pyrimidine Metabolism

# Km
pyrim_atp_km = atp_km[atp_km["Pathway"] == "Pyrimidine metabolism"]
pyrim_atp_km["Start Value"] = -1*np.log(pyrim_atp_km["Start Value"])
sns.swarmplot(x = "Gene names", y = "Start Value", data=pyrim_atp_km)
sns.violinplot(x = "Gene names", y = "Start Value", data=pyrim_atp_km)
plt.xticks(rotation=90)
plt.ylabel("-log(Km)")
plt.title("Enzyme sensitivity for ATP in Pyrimidine Metabolism")
plt.savefig("pyrimidine_ATP_kmvar.svg", dpi=400)

# kcat
pyrim_atp_kcat = atp_kcat[atp_kcat["Pathway"] == "Pyrimidine metabolism"]
pyrim_atp_kcat["Start Value"] = np.log(pyrim_atp_kcat["Start Value"])
sns.swarmplot(x = "Gene names", y = "Start Value", data=pyrim_atp_kcat)
sns.violinplot(x = "Gene names", y = "Start Value", data=pyrim_atp_kcat)
plt.xticks(rotation=90)
plt.ylabel("log(kcat)")
plt.title("Catalytic efficiency for ATP in Pyrimidine Metabolism")
plt.savefig("pyrimidine_ATP_kcatvar.svg", dpi=400)

# kcat/Km
pyrim_atp_km = atp_km[atp_km["Pathway"] == "Pyrimidine metabolism"]
pyrim_atp_kcat = atp_kcat[atp_kcat["Pathway"] == "Pyrimidine metabolism"]
pyrim_atp_kcat["Start Value"] = np.log(pyrim_atp_kcat["Start Value"])
sns.swarmplot(x = "Gene names", y = "Start Value", data=pyrim_atp_kcat)
sns.violinplot(x = "Gene names", y = "Start Value", data=pyrim_atp_kcat)
plt.xticks(rotation=90)
plt.ylabel("log(kcat)")
plt.title("Catalytic efficiency for ATP in Pyrimidine Metabolism")
plt.savefig("pyrimidine_ATP_kcatvar.svg", dpi=400)

# Glycolysis/Gluconeogenesis Metabolism
gly_atp_km = atp_km[atp_km["Pathway"] == ("Glycolysis/Gluconeogenesis")]
gly_atp_km["Start Value"] = np.log(gly_atp_km["Start Value"])
sns.swarmplot(x = "Gene names", y = "Start Value", data=gly_atp_km)
sns.violinplot(x = "Gene names", y = "Start Value", data=gly_atp_km)
plt.xticks(rotation=90)

# Nicotinate and Nicotinamide metabolism
nic_atp_km = atp_km[atp_km["Pathway"] == ("Nicotinate and Nicotinamide metabolism")]
nic_atp_km["Start Value"] = np.log(nic_atp_km["Start Value"])
sns.swarmplot(x = "Gene names", y = "Start Value", data=nic_atp_km)
sns.violinplot(x = "Gene names", y = "Start Value", data=nic_atp_km)
plt.xticks(rotation=90)

# Sulfur metabolism
sul_atp_km = atp_km[atp_km["Pathway"] == ("Sulfur metabolism")]
sul_atp_km["Start Value"] = np.log(sul_atp_km["Start Value"])
sns.swarmplot(x = "Gene names", y = "Start Value", data=sul_atp_km)
sns.violinplot(x = "Gene names", y = "Start Value", data=sul_atp_km)
plt.xticks(rotation=90)

################################## ADP ########################################
adp_km = paths[paths["Species"] == "ADP"]
adp_km = adp_km[adp_km["Type.1"] == "Km"]

# Purine Metabolism
purine_adp_km = adp_km[adp_km["Pathway"] == "Purine metabolism"]
purine_adp_km["Start Value"] = np.log(purine_adp_km["Start Value"])
sns.swarmplot(x = "Gene names", y = "Start Value", data=purine_adp_km)
sns.violinplot(x = "Gene names", y = "Start Value", data=purine_adp_km)
plt.xticks(rotation=90)

# Glycolysis/Gluconeogenesis Metabolism
gly_adp_km = adp_km[adp_km["Pathway"] == ("Glycolysis/Gluconeogenesis")]
gly_adp_km["Start Value"] = np.log(gly_adp_km["Start Value"])
sns.swarmplot(x = "Gene names", y = "Start Value", data=gly_adp_km)
sns.violinplot(x = "Gene names", y = "Start Value", data=gly_adp_km)
plt.xticks(rotation=90)

# Pyruvate metabolism
pyr_adp_km = adp_km[adp_km["Pathway"] == ("Pyruvate metabolism")]
pyr_adp_km["Start Value"] = np.log(pyr_adp_km["Start Value"])
sns.swarmplot(x = "Gene names", y = "Start Value", data=pyr_adp_km)
sns.violinplot(x = "Gene names", y = "Start Value", data=pyr_adp_km)
plt.xticks(rotation=90)

################################## NAD ########################################
nad_km = paths[paths["Species"] == "NAD+"]
nad_km = nad_km[nad_km["Type.1"] == "Km"]

# Arginine and Proline Metabolism
argpro_nad_km = nad_km[nad_km["Pathway"] == "Arginine and Proline metabolism"]
argpro_nad_km["Start Value"] = np.log(argpro_nad_km["Start Value"])
sns.swarmplot(x = "Gene names", y = "Start Value", data=argpro_nad_km)
sns.violinplot(x = "Gene names", y = "Start Value", data=argpro_nad_km)
plt.xticks(rotation=90)

# Fatty Acid Metabolism
fa_nad_km = nad_km[nad_km["Pathway"] == ("Fatty acid metabolism")]
fa_nad_km["Start Value"] = np.log(fa_nad_km["Start Value"])
sns.swarmplot(x = "Gene names", y = "Start Value", data=fa_nad_km)
sns.violinplot(x = "Gene names", y = "Start Value", data=fa_nad_km)
plt.xticks(rotation=90)

# Pyruvate metabolism
pyr_nad_km = nad_km[nad_km["Pathway"] == ("Pyruvate metabolism")]
pyr_nad_km["Start Value"] = np.log(pyr_nad_km["Start Value"])
sns.swarmplot(x = "Gene names", y = "Start Value", data=pyr_nad_km)
sns.violinplot(x = "Gene names", y = "Start Value", data=pyr_nad_km)
plt.xticks(rotation=90)

###############################################################################
############################# Pathway Analysis ################################
###############################################################################

###############################################################################
# Glycine, Serine, and Threonine metabolism ###################################
glyserthr = sabiork[sabiork["Pathway"]=="Glycine, Serine and Threonine metabolism"]
glyserthr_km = glyserthr[glyserthr["Type.1"] == "Km"]

# Look at individual variance for metabolites:
# ATP
glyserthr_km["Start Value"] = np.log(glyserthr_km["Start Value"])

# Aggregate by gene and species for Km
hm_glyserthr_km = pd.pivot_table(glyserthr_km, values="Start Value", index="Gene names", columns="Species", aggfunc=np.median)
hm_glyserthr_km = np.log(hm_glyserthr_km)
hm_glyserthr_km = hm_glyserthr_km.fillna(0)

# Plot the log(Km)
import matplotlib.pyplot as plt
import seaborn as sns

sns.heatmap(data = hm_glyserthr_km)
plt.title("Gly, Ser, and Thr Metabolism median log(Km)")
plt.xticks(rotation=90)

# Aggregate by gene and species for the kcat
glyserthr_kcat = glyserthr[glyserthr["Type.1"] == "kcat"]
s = glyserthr_kcat['Substrate'].str.split(';').apply(pd.Series, 1).stack()
s.index = s.index.droplevel(-1)
s.name = 'Substrate'
del glyserthr_kcat['Substrate']

# drop duplicates
glyserthr_kcat = glyserthr_kcat.join(s)
glyserthr_kcat = glyserthr_kcat.drop_duplicates(keep='first')

hm_glyserthr_kcat = pd.pivot_table(glyserthr_kcat, values="Start Value", index="Gene names", columns="Substrate", aggfunc=np.median)
hm_glyserthr_kcat = np.log(hm_glyserthr_kcat)
hm_glyserthr_kcat = hm_glyserthr_kcat.fillna(0)

# Plot the log(kcat)
import matplotlib.pyplot as plt
import seaborn as sns

sns.heatmap(data = hm_glyserthr_kcat)
plt.title("Glycine, Serine, and Threonine Metabolism")
plt.xticks(rotation=90)

###############################################################################
# Cysteine and Methionine metabolism ###################################
cysmet = sabiork[sabiork["Pathway"]=="Cysteine and methionine metabolism"]
cysmet_km = cysmet[cysmet["Type.1"] == "Km"]

# Look at individual variance for metabolites:
# ATP
cysmet_km["Start Value"] = -1*np.log(cysmet_km["Start Value"])

# Aggregate by gene and species for Km
hm_cysmet_km = pd.pivot_table(cysmet_km, values="Start Value", index="Gene names", columns="Species", aggfunc=np.median)
hm_cysmet_km = np.log(hm_cysmet_km)
hm_cysmet_km = hm_cysmet_km.fillna(0)

# Plot the log(Km)
import matplotlib.pyplot as plt
import seaborn as sns

sns.heatmap(data = hm_cysmet_km)
plt.title("Cysteine and Methionine metabolism median log(Km)")
plt.xticks(rotation=90)

# Aggregate by gene and species for the kcat
glyserthr_kcat = glyserthr[glyserthr["Type.1"] == "kcat"]
s = glyserthr_kcat['Substrate'].str.split(';').apply(pd.Series, 1).stack()
s.index = s.index.droplevel(-1)
s.name = 'Substrate'
del glyserthr_kcat['Substrate']

# drop duplicates
glyserthr_kcat = glyserthr_kcat.join(s)
glyserthr_kcat = glyserthr_kcat.drop_duplicates(keep='first')

hm_glyserthr_kcat = pd.pivot_table(glyserthr_kcat, values="Start Value", index="Gene names", columns="Substrate", aggfunc=np.median)
hm_glyserthr_kcat = np.log(hm_glyserthr_kcat)
hm_glyserthr_kcat = hm_glyserthr_kcat.fillna(0)

# Plot the log(kcat)
import matplotlib.pyplot as plt
import seaborn as sns

sns.heatmap(data = hm_glyserthr_kcat)
plt.title("Glycine, Serine, and Threonine Metabolism")
plt.xticks(rotation=90)


###############################################################################
########################### KEGG Pathway Annotations ##########################
###############################################################################

# Get gene lists
alaaspglu = r"KEGG_ALANINE_ASPARTATE_AND_GLUTAMATE_METABOLISM.txt"
alaaspglu = pd.read_csv(alaaspglu)
alaaspglu_atp = pd.merge(alaaspglu, atp, how="inner", left_on="GENE", right_on="Gene")
alaaspglu_adp = pd.merge(alaaspglu, adp, how="inner", left_on="GENE", right_on="Gene")
alaaspglu_nad = pd.merge(alaaspglu, nad, how="inner", left_on="GENE", right_on="Gene")
alaaspglu_nadh = pd.merge(alaaspglu, nadh, how="inner", left_on="GENE", right_on="Gene")

aminonuc = r"KEGG_AMINO_SUGAR_AND_NUCLEOTIDE_SUGAR_METABOLISM.txt"
aminonuc = pd.read_csv(aminonuc)
aminonuc_atp = pd.merge(aminonuc, atp, how="inner", left_on="GENE", right_on="Gene")
aminonuc_adp = pd.merge(aminonuc, adp, how="inner", left_on="GENE", right_on="Gene")
aminonuc_nad = pd.merge(aminonuc, nad, how="inner", left_on="GENE", right_on="Gene")
aminonuc_nadh = pd.merge(aminonuc, nadh, how="inner", left_on="GENE", right_on="Gene")

arg_pro = r"KEGG_ARGININE_AND_PROLINE_METABOLISM.txt"
arg_pro = pd.read_csv(arg_pro)
arg_pro_atp = pd.merge(arg_pro, atp, how="inner", left_on="GENE", right_on="Gene")
arg_pro_adp = pd.merge(arg_pro, adp, how="inner", left_on="GENE", right_on="Gene")
arg_pro_nad = pd.merge(arg_pro, nad, how="inner", left_on="GENE", right_on="Gene")
arg_pro_nadh = pd.merge(arg_pro, nadh, how="inner", left_on="GENE", right_on="Gene")

tca = r"KEGG_CITRATE_CYCLE_TCA_CYCLE.txt"
tca = pd.read_csv(tca)
tca_atp = pd.merge(tca, atp, how="inner", left_on="GENE", right_on="Gene")
tca_adp = pd.merge(tca, adp, how="inner", left_on="GENE", right_on="Gene")
tca_nad = pd.merge(tca, nad, how="inner", left_on="GENE", right_on="Gene")
tca_nadh = pd.merge(tca, nadh, how="inner", left_on="GENE", right_on="Gene")

cys_met = r"KEGG_CYSTEINE_AND_METHIONINE_METABOLISM.txt"
cys_met = pd.read_csv(cys_met)
cys_met_atp = pd.merge(cys_met, atp, how="inner", left_on="GENE", right_on="Gene")
cys_met_adp = pd.merge(cys_met, adp, how="inner", left_on="GENE", right_on="Gene")
cys_met_nad = pd.merge(cys_met, nad, how="inner", left_on="GENE", right_on="Gene")
cys_met_nadh = pd.merge(cys_met, nadh, how="inner", left_on="GENE", right_on="Gene")

fatty_acid = r"KEGG_FATTY_ACID_METABOLISM.txt"
fatty_acid = pd.read_csv(fatty_acid)
fatty_acid_atp = pd.merge(fatty_acid, atp, how="inner", left_on="GENE", right_on="Gene")
fatty_acid_adp = pd.merge(fatty_acid, adp, how="inner", left_on="GENE", right_on="Gene")
fatty_acid_nad = pd.merge(fatty_acid, nad, how="inner", left_on="GENE", right_on="Gene")
fatty_acid_nadh = pd.merge(fatty_acid, nadh, how="inner", left_on="GENE", right_on="Gene")

folate = r"KEGG_FOLATE_BIOSYNTHESIS.txt"
folate = pd.read_csv(folate)
folate_atp = pd.merge(folate, atp, how="inner", left_on="GENE", right_on="Gene")
folate_adp = pd.merge(folate, adp, how="inner", left_on="GENE", right_on="Gene")
folate_nad = pd.merge(folate, nad, how="inner", left_on="GENE", right_on="Gene")
folate_nadh = pd.merge(folate, nadh, how="inner", left_on="GENE", right_on="Gene")

glutathione = r"KEGG_GLUTATHIONE_METABOLISM.txt"
glutathione = pd.read_csv(glutathione)
glutathione_atp = pd.merge(glutathione, atp, how="inner", left_on="GENE", right_on="Gene")
glutathione_adp = pd.merge(glutathione, adp, how="inner", left_on="GENE", right_on="Gene")
glutathione_nad = pd.merge(glutathione, nad, how="inner", left_on="GENE", right_on="Gene")
glutathione_nadh = pd.merge(glutathione, nadh, how="inner", left_on="GENE", right_on="Gene")

glyserthr = r"KEGG_GLYCINE_SERINE_AND_THREONINE_METABOLISM.txt"
glyserthr = pd.read_csv(glyserthr)
glyserthr_atp = pd.merge(glyserthr, atp, how="inner", left_on="GENE", right_on="Gene")
glyserthr_adp = pd.merge(glyserthr, adp, how="inner", left_on="GENE", right_on="Gene")
glyserthr_nad = pd.merge(glyserthr, nad, how="inner", left_on="GENE", right_on="Gene")
glyserthr_nadh = pd.merge(glyserthr, nadh, how="inner", left_on="GENE", right_on="Gene")

glycolysis = r"KEGG_GLYCOLYSIS_GLUCONEOGENESIS.txt"
glycolysis = pd.read_csv(glycolysis)
glycolysis_atp = pd.merge(glycolysis, atp, how="inner", left_on="GENE", right_on="Gene")
glycolysis_adp = pd.merge(glycolysis, adp, how="inner", left_on="GENE", right_on="Gene")
glycolysis_nad = pd.merge(glycolysis, nad, how="inner", left_on="GENE", right_on="Gene")
glycolysis_nadh = pd.merge(glycolysis, nadh, how="inner", left_on="GENE", right_on="Gene")

his = r"KEGG_HISTIDINE_METABOLISM.txt"
his = pd.read_csv(his)
his_atp = pd.merge(his, atp, how="inner", left_on="GENE", right_on="Gene")
his_adp = pd.merge(his, adp, how="inner", left_on="GENE", right_on="Gene")
his_nad = pd.merge(his, nad, how="inner", left_on="GENE", right_on="Gene")
his_nadh = pd.merge(his, nadh, how="inner", left_on="GENE", right_on="Gene")

lys_deg = r"KEGG_LYSINE_DEGRADATION.txt"
lys_deg = pd.read_csv(lys_deg)
lys_deg_atp = pd.merge(lys_deg, atp, how="inner", left_on="GENE", right_on="Gene")
lys_deg_adp = pd.merge(lys_deg, adp, how="inner", left_on="GENE", right_on="Gene")
lys_deg_nad = pd.merge(lys_deg, nad, how="inner", left_on="GENE", right_on="Gene")
lys_deg_nadh = pd.merge(lys_deg, nadh, how="inner", left_on="GENE", right_on="Gene")

nitrogen = r"KEGG_NITROGEN_METABOLISM.txt"
nitrogen = pd.read_csv(nitrogen)
nitrogen_atp = pd.merge(nitrogen, atp, how="inner", left_on="GENE", right_on="Gene")
nitrogen_adp = pd.merge(nitrogen, adp, how="inner", left_on="GENE", right_on="Gene")
nitrogen_nad = pd.merge(nitrogen, nad, how="inner", left_on="GENE", right_on="Gene")
nitrogen_nadh = pd.merge(nitrogen, nadh, how="inner", left_on="GENE", right_on="Gene")

oxphos = r"KEGG_OXIDATIVE_PHOSPHORYLATION.txt"
oxphos = pd.read_csv(oxphos)
oxphos_atp = pd.merge(oxphos, atp, how="inner", left_on="GENE", right_on="Gene")
oxphos_adp = pd.merge(oxphos, adp, how="inner", left_on="GENE", right_on="Gene")
oxphos_nad = pd.merge(oxphos, nad, how="inner", left_on="GENE", right_on="Gene")
oxphos_nadh = pd.merge(oxphos, nadh, how="inner", left_on="GENE", right_on="Gene")

onecar = r"KEGG_ONE_CARBON_POOL_BY_FOLATE.txt"
onecar = pd.read_csv(onecar)
onecar_atp = pd.merge(onecar, atp, how="inner", left_on="GENE", right_on="Gene")
onecar_adp = pd.merge(onecar, adp, how="inner", left_on="GENE", right_on="Gene")
onecar_nad = pd.merge(onecar, nad, how="inner", left_on="GENE", right_on="Gene")
onecar_nadh = pd.merge(onecar, nadh, how="inner", left_on="GENE", right_on="Gene")

b12 = r"KEGG_PANTOTHENATE_AND_COA_BIOSYNTHESIS.txt"
b12 = pd.read_csv(b12)
b12_atp = pd.merge(b12, atp, how="inner", left_on="GENE", right_on="Gene")
b12_adp = pd.merge(b12, adp, how="inner", left_on="GENE", right_on="Gene")
b12_nad = pd.merge(b12, nad, how="inner", left_on="GENE", right_on="Gene")
b12_nadh = pd.merge(b12, nadh, how="inner", left_on="GENE", right_on="Gene")

ppp = r"KEGG_PENTOSE_PHOSPHATE_PATHWAY.txt"
ppp = pd.read_csv(ppp)
ppp_atp = pd.merge(ppp, atp, how="inner", left_on="GENE", right_on="Gene")
ppp_adp = pd.merge(ppp, adp, how="inner", left_on="GENE", right_on="Gene")
ppp_nad = pd.merge(ppp, nad, how="inner", left_on="GENE", right_on="Gene")
ppp_nadh = pd.merge(ppp, nadh, how="inner", left_on="GENE", right_on="Gene")

pyruvate = r"KEGG_PYRUVATE_METABOLISM.txt"
pyruvate = pd.read_csv(pyruvate)
pyruvate_atp = pd.merge(pyruvate, atp, how="inner", left_on="GENE", right_on="Gene")
pyruvate_adp = pd.merge(pyruvate, adp, how="inner", left_on="GENE", right_on="Gene")
pyruvate_nad = pd.merge(pyruvate, nad, how="inner", left_on="GENE", right_on="Gene")
pyruvate_nadh = pd.merge(pyruvate, nadh, how="inner", left_on="GENE", right_on="Gene")

pur = r"KEGG_PURINE_METABOLISM.txt"
pur = pd.read_csv(pur)
pur_atp = pd.merge(pur, atp, how="inner", left_on="GENE", right_on="Gene")
pur_adp = pd.merge(pur, adp, how="inner", left_on="GENE", right_on="Gene")
pur_nad = pd.merge(pur, nad, how="inner", left_on="GENE", right_on="Gene")
pur_nadh = pd.merge(pur, nadh, how="inner", left_on="GENE", right_on="Gene")

pyrim = r"KEGG_PYRIMIDINE_METABOLISM.txt"
pyrim = pd.read_csv(pyrim)
pyrim_atp = pd.merge(pyrim, atp, how="inner", left_on="GENE", right_on="Gene")
pyrim_adp = pd.merge(pyrim, adp, how="inner", left_on="GENE", right_on="Gene")
pyrim_nad = pd.merge(pyrim, nad, how="inner", left_on="GENE", right_on="Gene")
pyrim_nadh = pd.merge(pyrim, nadh, how="inner", left_on="GENE", right_on="Gene")

sulfur = r"KEGG_SULFUR_METABOLISM.txt"
sulfur = pd.read_csv(sulfur)
sulfur_atp = pd.merge(sulfur, atp, how="inner", left_on="GENE", right_on="Gene")
sulfur_adp = pd.merge(sulfur, adp, how="inner", left_on="GENE", right_on="Gene")
sulfur_nad = pd.merge(sulfur, nad, how="inner", left_on="GENE", right_on="Gene")
sulfur_nadh = pd.merge(sulfur, nadh, how="inner", left_on="GENE", right_on="Gene")

tryp = r"KEGG_TRYPTOPHAN_METABOLISM.txt"
tryp = pd.read_csv(tryp)
tryp_atp = pd.merge(tryp, atp, how="inner", left_on="GENE", right_on="Gene")
tryp_adp = pd.merge(tryp, adp, how="inner", left_on="GENE", right_on="Gene")
tryp_nad = pd.merge(tryp, nad, how="inner", left_on="GENE", right_on="Gene")
tryp_nadh = pd.merge(tryp, nadh, how="inner", left_on="GENE", right_on="Gene")

valeuiso_syn = r"KEGG_VALINE_LEUCINE_AND_ISOLEUCINE_BIOSYNTHESIS.txt"
valeuiso_syn = pd.read_csv(valeuiso_syn)
valeuiso_syn_atp = pd.merge(valeuiso_syn, atp, how="inner", left_on="GENE", right_on="Gene")
valeuiso_syn_adp = pd.merge(valeuiso_syn, adp, how="inner", left_on="GENE", right_on="Gene")
valeuiso_syn_nad = pd.merge(valeuiso_syn, nad, how="inner", left_on="GENE", right_on="Gene")
valeuiso_syn_nadh = pd.merge(valeuiso_syn, nadh, how="inner", left_on="GENE", right_on="Gene")

valeuiso_deg = r"KEGG_VALINE_LEUCINE_AND_ISOLEUCINE_DEGRADATION.txt"
valeuiso_deg = pd.read_csv(valeuiso_deg)
valeuiso_atp = pd.merge(valeuiso_deg, atp, how="inner", left_on="GENE", right_on="Gene")
valeuiso_adp = pd.merge(valeuiso_deg, adp, how="inner", left_on="GENE", right_on="Gene")
valeuiso_nad = pd.merge(valeuiso_deg, nad, how="inner", left_on="GENE", right_on="Gene")
valeuiso_nadh = pd.merge(valeuiso_deg, nadh, how="inner", left_on="GENE", right_on="Gene")

###############################################################################
############################ SabioRK Pathway Analysis #########################
###############################################################################
import pandas as pd
import numpy as np

sabio = pd.read_excel("kinetic_map.xlsx", sheet_name="SABIORK")
clean_sabio = sabio.drop(columns=["Tissue", "Pathway", "Type", "End Value", \
                                  "Standard Dev", "Unit", "Cofactor", \
                                  "Inhibitor", "ChebiID", "pH", "T", "KEGG",\
                                  "ReactionEquation", "Mechanism", "Product"])

# Organize sabioRK database by separating the substrates
s = clean_sabio['Substrate'].str.split(';').apply(pd.Series, 1).stack()
s.index = s.index.droplevel(-1)
s.name = 'Substrate'
del clean_sabio['Substrate']
# drop duplicates
clean_sabio = clean_sabio.join(s)
clean_sabio = clean_sabio.drop_duplicates(keep='first')

# Reconstruct SabioRK database so genes will map directly to Km and kcat
# To choose a single value, I chose the maximum value reported
km = clean_sabio[clean_sabio["Type.1"] == "Km"].drop(columns=["Type.1", \
                "Species"]).rename(columns={"Start Value":"Km"})
km = km.groupby(["Gene names", "Substrate"], sort=False)["Km"].max().reset_index()
kcat = clean_sabio[clean_sabio["Type.1"] == "kcat"].drop(columns=["Type.1", \
                  "Species"]).rename(columns={"Start Value":"kcat"})
kcat = kcat.groupby(["Gene names", "Substrate"], sort=False)["kcat"].max().reset_index()
comb = pd.merge(kcat, km, how="inner", left_on=["Gene names", "Substrate"],\
                right_on=["Gene names", "Substrate"]).drop_duplicates(keep="first").dropna()

# Get substrates of interest
atp = comb[comb["Substrate"]=="ATP"].sort_values(by="kcat", ascending=False)
adp = comb[comb["Substrate"]=="ADP"].sort_values(by="kcat", ascending=False)
nad = comb[comb["Substrate"]=="NAD+"].sort_values(by="kcat", ascending=False)
nadh = comb[comb["Substrate"]=="NADH"].sort_values(by="kcat", ascending=False)
acoa = comb[comb["Substrate"]=="Acetyl-CoA"].sort_values(by="kcat", ascending=False)
coa = comb[comb["Substrate"]=="Coenzyme A"].sort_values(by="kcat", ascending=False)

# Get gene lists
alaaspglu = r"KEGG_ALANINE_ASPARTATE_AND_GLUTAMATE_METABOLISM.txt"
alaaspglu = pd.read_csv(alaaspglu)
alaaspglu_atp = pd.merge(alaaspglu, atp, how="inner", left_on="GENE", right_on="Gene names")
alaaspglu_adp = pd.merge(alaaspglu, adp, how="inner", left_on="GENE", right_on="Gene names")
alaaspglu_nad = pd.merge(alaaspglu, nad, how="inner", left_on="GENE", right_on="Gene names")
alaaspglu_nadh = pd.merge(alaaspglu, nadh, how="inner", left_on="GENE", right_on="Gene names")

aminonuc = r"KEGG_AMINO_SUGAR_AND_NUCLEOTIDE_SUGAR_METABOLISM.txt"
aminonuc = pd.read_csv(aminonuc)
aminonuc_atp = pd.merge(aminonuc, atp, how="inner", left_on="GENE", right_on="Gene names")
aminonuc_adp = pd.merge(aminonuc, adp, how="inner", left_on="GENE", right_on="Gene names")
aminonuc_nad = pd.merge(aminonuc, nad, how="inner", left_on="GENE", right_on="Gene names")
aminonuc_nadh = pd.merge(aminonuc, nadh, how="inner", left_on="GENE", right_on="Gene names")

arg_pro = r"KEGG_ARGININE_AND_PROLINE_METABOLISM.txt"
arg_pro = pd.read_csv(arg_pro)
arg_pro_atp = pd.merge(arg_pro, atp, how="inner", left_on="GENE", right_on="Gene names")
arg_pro_adp = pd.merge(arg_pro, adp, how="inner", left_on="GENE", right_on="Gene names")
arg_pro_nad = pd.merge(arg_pro, nad, how="inner", left_on="GENE", right_on="Gene names")
arg_pro_nadh = pd.merge(arg_pro, nadh, how="inner", left_on="GENE", right_on="Gene names")

tca = r"KEGG_CITRATE_CYCLE_TCA_CYCLE.txt"
tca = pd.read_csv(tca)
tca_atp = pd.merge(tca, atp, how="inner", left_on="GENE", right_on="Gene names")
tca_adp = pd.merge(tca, adp, how="inner", left_on="GENE", right_on="Gene names")
tca_nad = pd.merge(tca, nad, how="inner", left_on="GENE", right_on="Gene names")
tca_nadh = pd.merge(tca, nadh, how="inner", left_on="GENE", right_on="Gene names")

cys_met = r"KEGG_CYSTEINE_AND_METHIONINE_METABOLISM.txt"
cys_met = pd.read_csv(cys_met)
cys_met_atp = pd.merge(cys_met, atp, how="inner", left_on="GENE", right_on="Gene names")
cys_met_adp = pd.merge(cys_met, adp, how="inner", left_on="GENE", right_on="Gene names")
cys_met_nad = pd.merge(cys_met, nad, how="inner", left_on="GENE", right_on="Gene names")
cys_met_nadh = pd.merge(cys_met, nadh, how="inner", left_on="GENE", right_on="Gene names")

fatty_acid = r"KEGG_FATTY_ACID_METABOLISM.txt"
fatty_acid = pd.read_csv(fatty_acid)
fatty_acid_atp = pd.merge(fatty_acid, atp, how="inner", left_on="GENE", right_on="Gene names")
fatty_acid_adp = pd.merge(fatty_acid, adp, how="inner", left_on="GENE", right_on="Gene names")
fatty_acid_nad = pd.merge(fatty_acid, nad, how="inner", left_on="GENE", right_on="Gene names")
fatty_acid_nadh = pd.merge(fatty_acid, nadh, how="inner", left_on="GENE", right_on="Gene names")

folate = r"KEGG_FOLATE_BIOSYNTHESIS.txt"
folate = pd.read_csv(folate)
folate_atp = pd.merge(folate, atp, how="inner", left_on="GENE", right_on="Gene names")
folate_adp = pd.merge(folate, adp, how="inner", left_on="GENE", right_on="Gene names")
folate_nad = pd.merge(folate, nad, how="inner", left_on="GENE", right_on="Gene names")
folate_nadh = pd.merge(folate, nadh, how="inner", left_on="GENE", right_on="Gene names")

glutathione = r"KEGG_GLUTATHIONE_METABOLISM.txt"
glutathione = pd.read_csv(glutathione)
glutathione_atp = pd.merge(glutathione, atp, how="inner", left_on="GENE", right_on="Gene names")
glutathione_adp = pd.merge(glutathione, adp, how="inner", left_on="GENE", right_on="Gene names")
glutathione_nad = pd.merge(glutathione, nad, how="inner", left_on="GENE", right_on="Gene names")
glutathione_nadh = pd.merge(glutathione, nadh, how="inner", left_on="GENE", right_on="Gene names")

glyserthr = r"KEGG_GLYCINE_SERINE_AND_THREONINE_METABOLISM.txt"
glyserthr = pd.read_csv(glyserthr)
glyserthr_atp = pd.merge(glyserthr, atp, how="inner", left_on="GENE", right_on="Gene names")
glyserthr_adp = pd.merge(glyserthr, adp, how="inner", left_on="GENE", right_on="Gene names")
glyserthr_nad = pd.merge(glyserthr, nad, how="inner", left_on="GENE", right_on="Gene names")
glyserthr_nadh = pd.merge(glyserthr, nadh, how="inner", left_on="GENE", right_on="Gene names")

glycolysis = r"KEGG_GLYCOLYSIS_GLUCONEOGENESIS.txt"
glycolysis = pd.read_csv(glycolysis)
glycolysis_atp = pd.merge(glycolysis, atp, how="inner", left_on="GENE", right_on="Gene names")
glycolysis_adp = pd.merge(glycolysis, adp, how="inner", left_on="GENE", right_on="Gene names")
glycolysis_nad = pd.merge(glycolysis, nad, how="inner", left_on="GENE", right_on="Gene names")
glycolysis_nadh = pd.merge(glycolysis, nadh, how="inner", left_on="GENE", right_on="Gene names")

his = r"KEGG_HISTIDINE_METABOLISM.txt"
his = pd.read_csv(his)
his_atp = pd.merge(his, atp, how="inner", left_on="GENE", right_on="Gene names")
his_adp = pd.merge(his, adp, how="inner", left_on="GENE", right_on="Gene names")
his_nad = pd.merge(his, nad, how="inner", left_on="GENE", right_on="Gene names")
his_nadh = pd.merge(his, nadh, how="inner", left_on="GENE", right_on="Gene names")

lys_deg = r"KEGG_LYSINE_DEGRADATION.txt"
lys_deg = pd.read_csv(lys_deg)
lys_deg_atp = pd.merge(lys_deg, atp, how="inner", left_on="GENE", right_on="Gene names")
lys_deg_adp = pd.merge(lys_deg, adp, how="inner", left_on="GENE", right_on="Gene names")
lys_deg_nad = pd.merge(lys_deg, nad, how="inner", left_on="GENE", right_on="Gene names")
lys_deg_nadh = pd.merge(lys_deg, nadh, how="inner", left_on="GENE", right_on="Gene names")

nitrogen = r"KEGG_NITROGEN_METABOLISM.txt"
nitrogen = pd.read_csv(nitrogen)
nitrogen_atp = pd.merge(nitrogen, atp, how="inner", left_on="GENE", right_on="Gene names")
nitrogen_adp = pd.merge(nitrogen, adp, how="inner", left_on="GENE", right_on="Gene names")
nitrogen_nad = pd.merge(nitrogen, nad, how="inner", left_on="GENE", right_on="Gene names")
nitrogen_nadh = pd.merge(nitrogen, nadh, how="inner", left_on="GENE", right_on="Gene names")

oxphos = r"KEGG_OXIDATIVE_PHOSPHORYLATION.txt"
oxphos = pd.read_csv(oxphos)
oxphos_atp = pd.merge(oxphos, atp, how="inner", left_on="GENE", right_on="Gene names")
oxphos_adp = pd.merge(oxphos, adp, how="inner", left_on="GENE", right_on="Gene names")
oxphos_nad = pd.merge(oxphos, nad, how="inner", left_on="GENE", right_on="Gene names")
oxphos_nadh = pd.merge(oxphos, nadh, how="inner", left_on="GENE", right_on="Gene names")

onecar = r"KEGG_ONE_CARBON_POOL_BY_FOLATE.txt"
onecar = pd.read_csv(onecar)
onecar_atp = pd.merge(onecar, atp, how="inner", left_on="GENE", right_on="Gene names")
onecar_adp = pd.merge(onecar, adp, how="inner", left_on="GENE", right_on="Gene names")
onecar_nad = pd.merge(onecar, nad, how="inner", left_on="GENE", right_on="Gene names")
onecar_nadh = pd.merge(onecar, nadh, how="inner", left_on="GENE", right_on="Gene names")

b12 = r"KEGG_PANTOTHENATE_AND_COA_BIOSYNTHESIS.txt"
b12 = pd.read_csv(b12)
b12_atp = pd.merge(b12, atp, how="inner", left_on="GENE", right_on="Gene names")
b12_adp = pd.merge(b12, adp, how="inner", left_on="GENE", right_on="Gene names")
b12_nad = pd.merge(b12, nad, how="inner", left_on="GENE", right_on="Gene names")
b12_nadh = pd.merge(b12, nadh, how="inner", left_on="GENE", right_on="Gene names")

ppp = r"KEGG_PENTOSE_PHOSPHATE_PATHWAY.txt"
ppp = pd.read_csv(ppp)
ppp_atp = pd.merge(ppp, atp, how="inner", left_on="GENE", right_on="Gene names")
ppp_adp = pd.merge(ppp, adp, how="inner", left_on="GENE", right_on="Gene names")
ppp_nad = pd.merge(ppp, nad, how="inner", left_on="GENE", right_on="Gene names")
ppp_nadh = pd.merge(ppp, nadh, how="inner", left_on="GENE", right_on="Gene names")

pyruvate = r"KEGG_PYRUVATE_METABOLISM.txt"
pyruvate = pd.read_csv(pyruvate)
pyruvate_atp = pd.merge(pyruvate, atp, how="inner", left_on="GENE", right_on="Gene names")
pyruvate_adp = pd.merge(pyruvate, adp, how="inner", left_on="GENE", right_on="Gene names")
pyruvate_nad = pd.merge(pyruvate, nad, how="inner", left_on="GENE", right_on="Gene names")
pyruvate_nadh = pd.merge(pyruvate, nadh, how="inner", left_on="GENE", right_on="Gene names")

pur = r"KEGG_PURINE_METABOLISM.txt"
pur = pd.read_csv(pur)
pur_atp = pd.merge(pur, atp, how="inner", left_on="GENE", right_on="Gene names")
pur_adp = pd.merge(pur, adp, how="inner", left_on="GENE", right_on="Gene names")
pur_nad = pd.merge(pur, nad, how="inner", left_on="GENE", right_on="Gene names")
pur_nadh = pd.merge(pur, nadh, how="inner", left_on="GENE", right_on="Gene names")

pyrim = r"KEGG_PYRIMIDINE_METABOLISM.txt"
pyrim = pd.read_csv(pyrim)
pyrim_atp = pd.merge(pyrim, atp, how="inner", left_on="GENE", right_on="Gene names")
pyrim_adp = pd.merge(pyrim, adp, how="inner", left_on="GENE", right_on="Gene names")
pyrim_nad = pd.merge(pyrim, nad, how="inner", left_on="GENE", right_on="Gene names")
pyrim_nadh = pd.merge(pyrim, nadh, how="inner", left_on="GENE", right_on="Gene names")

sulfur = r"KEGG_SULFUR_METABOLISM.txt"
sulfur = pd.read_csv(sulfur)
sulfur_atp = pd.merge(sulfur, atp, how="inner", left_on="GENE", right_on="Gene names")
sulfur_adp = pd.merge(sulfur, adp, how="inner", left_on="GENE", right_on="Gene names")
sulfur_nad = pd.merge(sulfur, nad, how="inner", left_on="GENE", right_on="Gene names")
sulfur_nadh = pd.merge(sulfur, nadh, how="inner", left_on="GENE", right_on="Gene names")

tryp = r"KEGG_TRYPTOPHAN_METABOLISM.txt"
tryp = pd.read_csv(tryp)
tryp_atp = pd.merge(tryp, atp, how="inner", left_on="GENE", right_on="Gene names")
tryp_adp = pd.merge(tryp, adp, how="inner", left_on="GENE", right_on="Gene names")
tryp_nad = pd.merge(tryp, nad, how="inner", left_on="GENE", right_on="Gene names")
tryp_nadh = pd.merge(tryp, nadh, how="inner", left_on="GENE", right_on="Gene names")

valeuiso_syn = r"KEGG_VALINE_LEUCINE_AND_ISOLEUCINE_BIOSYNTHESIS.txt"
valeuiso_syn = pd.read_csv(valeuiso_syn)
valeuiso_syn_atp = pd.merge(valeuiso_syn, atp, how="inner", left_on="GENE", right_on="Gene names")
valeuiso_syn_adp = pd.merge(valeuiso_syn, adp, how="inner", left_on="GENE", right_on="Gene names")
valeuiso_syn_nad = pd.merge(valeuiso_syn, nad, how="inner", left_on="GENE", right_on="Gene names")
valeuiso_syn_nadh = pd.merge(valeuiso_syn, nadh, how="inner", left_on="GENE", right_on="Gene names")

valeuiso_deg = r"KEGG_VALINE_LEUCINE_AND_ISOLEUCINE_DEGRADATION.txt"
valeuiso_deg = pd.read_csv(valeuiso_deg)
valeuiso_atp = pd.merge(valeuiso_deg, atp, how="inner", left_on="GENE", right_on="Gene names")
valeuiso_adp = pd.merge(valeuiso_deg, adp, how="inner", left_on="GENE", right_on="Gene names")
valeuiso_nad = pd.merge(valeuiso_deg, nad, how="inner", left_on="GENE", right_on="Gene names")
valeuiso_nadh = pd.merge(valeuiso_deg, nadh, how="inner", left_on="GENE", right_on="Gene names")

###############################################################################