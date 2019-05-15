# -*- coding: utf-8 -*-
"""
Created on Sat Nov 10 15:56:08 2018
Kinetic exploratory analysis of the BRENDA and the Sabio RK databases
    * Path is C/Users/Desktop/Kinetic Model/Data/brenda/human

@author: Scott Campit
"""

###############################################################################
# Important functions
import pandas as pd
import numpy as np

def kcat_km_map(df, metabolites = []):
    """
    kcat_km_map will automate the process of creating a dataframe containing
    the metabolites of interest within a list and output the kcat and Kms for
    further processing.

    """
    atp = pd.DataFrame([])
    adp = pd.DataFrame([])
    nad = pd.DataFrame([])
    nadh = pd.DataFrame([])
    acoa = pd.DataFrame([])
    coa = pd.DataFrame([])
    if "ATP" in metabolites:
        atp = df[df["Substrate"] == "ATP"]
    elif "ADP" in metabolites:
        adp = df[df["Substrate"] == "ADP"]
    elif "NAD+" in metabolites:
        nad = df[df["Substrate"] == "NAD+"]
    elif "NADH" in metabolites:
        nadh = df[df["Substrate"] == "NADH"]
    elif "Acetyl-CoA" in metabolites:
        acoa = df[df["Substrate"] == "Acetyl-CoA"]
    elif "Coenzyme A" in metabolites:
        coa = df[df["Substrate"] == "Coenzyme A"]
    else:
        pass
    return atp, adp, nad, nadh, acoa, coa

###############################################################################
########################### BRENDA Pathway Analysis ###########################
###############################################################################
#path = r"C:/Users/scott/OneDrive/Desktop/kinetic-project/Data/brenda/human"

brenda = pd.read_excel("kinetic_map.xlsx", sheet_name="BRENDA")

# Substrates of interest
atp = brenda[brenda["Substrate"]=="ATP"]
adp = brenda[brenda["Substrate"]=="ADP"]
nad = brenda[brenda["Substrate"]=="NAD+"]
nadh = brenda[brenda["Substrate"]=="NADH"]
acoa = brenda[brenda["Substrate"]=="Acetyl-CoA"]
coa = brenda[brenda["Substrate"]=="Coenzyme A"]
akg = brenda[brenda["Substrate"] == "2-oxoglutarate"]
icit = brenda[brenda["Substrate"] == "Isocitrate"]

# Scatterplot of the kcat/Km relationship appears to be nonlinear
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="white")

fig = plt.figure()

ax1 = sns.scatterplot(x="kcat", y="Km", data=atp)
sns.despine()
plt.title("Global kcat/Km associations for ATP corresponding to each gene")
plt.savefig("./figures/brenda_ATP_kcat_km.svg", dpi=400)

ax2 = sns.scatterplot(x="kcat", y="Km", data=adp)
sns.despine()
plt.title("Global kcat/Km associations for ADP corresponding to each gene")
plt.savefig("./figures/brenda_ADP_kcat_km.svg", dpi=400)

ax3 = sns.scatterplot(x="kcat", y="Km", data=nad)
sns.despine()
plt.title("Global kcat/Km associations for NAD+ corresponding to each gene")
plt.savefig("./figures/brenda_NAD_kcat_km.svg", dpi=400)

ax4 = sns.scatterplot(x="kcat", y="Km", data=nadh)
sns.despine()
plt.title("Global kcat/Km associations for NADH corresponding to each gene")
plt.savefig("./figures/brenda_NADH_kcat_km.svg", dpi=400)

ax5 = sns.scatterplot(x="kcat", y="Km", data=akg)
sns.despine()
plt.title("Global kcat/Km associations for 2-oxoglutarate corresponding to each gene")
plt.savefig("./figures/brenda_aKG_kcat_km.svg", dpi=400)

plt.set_xlim(xmin=0, xmax=100)
plt.set_ylim(ymin=0, ymax=5)

###############################################################################
# Barplot for Km values #######################################################

import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="white")

fig = plt.figure()
#fig.subplots_adjust(hspace=0.4, wspace=0.4)

ax1 = fig.add_subplot(2, 3, 1)
ax1 = sns.barplot(x="Gene", y="Km", data=atp)
sns.despine()
plt.xticks(rotation=90)
plt.title("Global gene-Km associations for ATP")
plt.savefig("./figures/brenda_ATP_gene_km.svg", dpi=400)

ax2 = fig.add_subplot(2, 3, 2)
ax2 = sns.barplot(x="Gene", y="Km", data=adp)
sns.despine()
plt.xticks(rotation=90)
plt.title("Global gene-Km associations for ADP")
plt.savefig("./figures/brenda_ADP_gene_km.svg", dpi=400)

ax3 = fig.add_subplot(2, 3, 3)
ax3 = sns.barplot(x="Gene", y="Km", data=nad)
sns.despine()
plt.xticks(rotation=90)
plt.title("Global gene-Km associations for NAD+")
plt.savefig("./figures/brenda_NAD_gene_km.svg", dpi=400)

ax4 = fig.add_subplot(2, 3, 4)
ax4 = sns.barplot(x="Gene", y="Km", data=nadh)
sns.despine()
plt.xticks(rotation=90)
plt.title("Global gene-Km associations for NADH")
plt.savefig("./figures/brenda_NADH_gene_km.svg", dpi=400)

ax5 = fig.add_subplot(2, 3, 4)
ax5 = sns.barplot(x="Gene", y="Km", data=akg)
sns.despine()
plt.xticks(rotation=90)
plt.title("Global gene-Km associations for aKG")
plt.savefig("./figures/brenda_aKG_gene_km.svg", dpi=400)

plt.set_xlim(xmin=0, xmax=1000)
plt.set_ylim(ymin=0, ymax=2)

###############################################################################
# Barplot for kcat values #####################################################

import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="white")
fig = plt.figure()

ax1 = sns.barplot(x="Gene", y="kcat", data=atp)
sns.despine()
plt.xticks(rotation=90)
plt.title("Global gene-kcat associations for ATP")

ax2 = sns.barplot(x="Gene", y="kcat", data=adp)
sns.despine()
plt.xticks(rotation=90)
plt.title("Global gene-kcat associations for ADP")

ax3 = sns.barplot(x="Gene", y="kcat", data=nad)
sns.despine()
plt.xticks(rotation=90)
plt.title("Global gene-kcat associations for NAD+")

ax4 = sns.barplot(x="Gene", y="kcat", data=nadh)
sns.despine()
plt.xticks(rotation=90)
plt.title("Global gene-kcat associations for NADH")

plt.set_xlim(xmin=0, xmax=1000)
plt.set_ylim(ymin=0, ymax=2)

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

###############################################################################
# Scatterplot for global analysis #############################################

import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="white")

fig = plt.figure()

ax1 = sns.scatterplot(x="kcat", y="Km", data=atp)
sns.despine()
plt.title("Global kcat/Km associations for ATP corresponding to each gene")

ax2 = sns.scatterplot(x="kcat", y="Km", data=adp)
sns.despine()
plt.title("Global kcat/Km associations for ADP corresponding to each gene")

ax3 = sns.scatterplot(x="kcat", y="Km", data=nad)
sns.despine()
plt.title("Global kcat/Km associations for NAD+ corresponding to each gene")

ax4 = sns.scatterplot(x="kcat", y="Km", data=nadh)
sns.despine()
plt.title("Global kcat/Km associations for NADH corresponding to each gene")

ax5 = sns.scatterplot(x="kcat", y="Km", data=acoa)
sns.despine()
plt.title("Global kcat/Km associations for Acetyl-CoA corresponding to each gene")

ax6 = sns.scatterplot(x="kcat", y="Km", data=coa)
sns.despine()
plt.title("Global kcat/Km associations for Coenzyme A corresponding to each gene")

plt.set_xlim(xmin=0, xmax=1000)
plt.set_ylim(ymin=0, ymax=2)

###############################################################################
# Barplot for Km values #######################################################

import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="white")

fig = plt.figure()
#fig.subplots_adjust(hspace=0.4, wspace=0.4)

ax1 = fig.add_subplot(2, 3, 1)
ax1 = sns.barplot(x="Gene names", y="Km", data=atp)
sns.despine()
plt.xticks(rotation=90)
plt.title("Global gene-Km associations for ATP")

ax2 = fig.add_subplot(2, 3, 2)
ax2 = sns.barplot(x="Gene names", y="Km", data=adp)
sns.despine()
plt.xticks(rotation=90)
plt.title("Global gene-Km associations for ADP")

ax3 = fig.add_subplot(2, 3, 3)
ax3 = sns.barplot(x="Gene names", y="Km", data=nad)
sns.despine()
plt.xticks(rotation=90)
plt.title("Global gene-Km associations for NAD+")

ax4 = fig.add_subplot(2, 3, 4)
ax4 = sns.barplot(x="Gene names", y="Km", data=nadh)
sns.despine()
plt.xticks(rotation=90)
plt.title("Global gene-Km associations for NADH")

ax5 = fig.add_subplot(2, 3, 5)
ax5 = sns.barplot(x="Gene names", y="Km", data=acoa)
sns.despine()
plt.xticks(rotation=90)
plt.title("Global gene-Km associations for Acetyl-CoA")

ax6 = fig.add_subplot(2, 3, 6)
ax6 = sns.barplot(x="Gene names", y="Km", data=coa)
sns.despine()
plt.xticks(rotation=90)
plt.title("Global gene-Km associations for Coenzyme A")

plt.set_xlim(xmin=0, xmax=1000)
plt.set_ylim(ymin=0, ymax=2)

###############################################################################
# Barplot for kcat values #####################################################

import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="white")

fig = plt.figure()

ax1 = sns.barplot(x="Gene names", y="kcat", data=atp)
sns.despine()
plt.xticks(rotation=90)
plt.title("Global gene-kcat associations for ATP")

ax2 = sns.barplot(x="Gene names", y="kcat", data=adp)
sns.despine()
plt.xticks(rotation=90)
plt.title("Global gene-kcat associations for ADP")

ax3 = sns.barplot(x="Gene names", y="kcat", data=nad)
sns.despine()
plt.xticks(rotation=90)
plt.title("Global gene-kcat associations for NAD+")

ax4 = sns.barplot(x="Gene names", y="kcat", data=nadh)
sns.despine()
plt.xticks(rotation=90)
plt.title("Global gene-kcat associations for NADH")

ax5 = sns.barplot(x="Gene names", y="kcat", data=acoa)
sns.despine()
plt.xticks(rotation=90)
plt.title("Global gene-kcat associations for Acetyl-CoA")

ax6 = sns.barplot(x="Gene names", y="kcat", data=coa)
sns.despine()
plt.xticks(rotation=90)
plt.title("Global gene-kcat associations for Coenzyme A")

plt.set_xlim(xmin=0, xmax=1000)
plt.set_ylim(ymin=0, ymax=2)

###############################################################################