'''
sbml_parse.py function: Takes a directory with sbml files and parses the
following info:
    a) the file name, which is its UNIPROT ID
    b) the gene ID, from a separate .csv file if it's human data
    c) species name with a Km
    d) metabolite with a Km value
    e) the km value
    f) the CHEBI ID associated with the Km value
    g) prints output in .txt file

'''

import xml.etree.cElementTree as ET
import codecs, sys, os, csv
import pandas as pd
from os.path import basename

def sbml_parse(compartment='compartment')
    print("UNIPROT ID"+"\t"+"GENE ID"+'\t'+"NAME"+"\t"+"METABOLITE"+"\t"+"Km Value"+"\t"+"CHEBI ID")
    path1 = 'C:/Users/scott/OneDrive/Desktop/Chandrasekaran/Projects/sabiodb/sabio-sbml'
    path2 = 'C:/Users/scott/OneDrive/Desktop/Chandrasekaran/Projects/uniprot/Human_Uniprot.csv'

    uniprotid = {}
    sbmlist = os.listdir(path1)
    for sbml in sbmlist:
        uniprotid = sbml
        with open(os.path.join(path1, uniprotid)) as sbmlfiles:
            # get rid of .sbml extension in file name for UNIPROT ID
            base=os.path.basename(uniprotid)
            splitbase=os.path.splitext(base)
            uniprotid = os.path.splitext(base)[0]
            # Get Gene IDs from mammaliaGeneID.tab
            df = pd.read_csv(path2, na_values=['.'])
            geneid = []
            for index, row in df.iterrows():
                if uniprotid in row[0]:
                    geneid = row[7]

            # use ET.parse to move through sbml in tree
            root = ET.parse(sbmlfiles)
            tree = root.getroot()
            naam = {
            'xmlns':"http://www.sbml.org/sbml/level2/version4",
            'html':"http://www.w3.org/1999/xhtml",
            'rdf':"http://www.w3.org/1999/02/22-rdf-syntax-ns#",
            'bqbiol':"http://biomodels.net/biology-qualifiers/",
            'xmlns':"http://www.sbml.org/sbml/level3/version1/core",
            'fbc':"http://www.sbml.org/sbml/level3/version1/fbc/version2"
            }

            name_id = []
            metabolite_id = []
            km_value = []

            # print km attributes in local parameters
            for child in tree.findall("./xmlns:model/xmlns:listOfReactions/xmlns:reaction/xmlns:kineticLaw/xmlns:listOfLocalParameters/xmlns:localParameter", naam):
                if ("Km" in child.attrib['id']):
                    try:
                        # prints the km_id
                        name_id = child.attrib["id"].split()
                        metabolite_id = child.attrib["name"].split()
                        km_value = child.attrib["value"].split()
                    except KeyError:
                        pass

                spec_id = []
                species = {}
                # for species, if they match Km name ID, I want to store that name into species
                for spec in tree.findall("./xmlns:model/xmlns:listOfSpecies/xmlns:species",naam):
                    spec_id = spec.attrib['id'].split()
                    for cell_line in spec_id:
                        #print(uniprotid+cell_line)
                        #prints each spec_id per uniprotid
                        for name in name_id:
                            if cell_line not in name:
                                pass
                            else:
                                spec_id = cell_line

                        # prints a specific compartment type
                        for intra in spec_id:
                            if compartment in spec_id:
                                spec_id = spec_id

                                chebids = {}
                                # Find CHEBI IDs and store in chebids dict from spec node
                                for res in spec.findall("./xmlns:annotation/rdf:RDF/rdf:Description/bqbiol:is/rdf:Bag/rdf:li",naam):
                                    if("CHEBI" in res.attrib["{"+naam['rdf']+"}resource"]):
                                        if(spec_id not in species):
                                            species[spec_id] = res.attrib["{"+naam['rdf']+"}resource"].split("/")[-1]
                                        else:
                                            species[spec_id] = species[spec_id]+"|"+res.attrib["{"+naam['rdf']+"}resource"].split("/")[-1]

                                        chebids[res.attrib["{"+naam['rdf']+"}resource"].split("/")[-1]] = 0

    with open('output_'+compartment+'_.txt', 'a') as file:
        file.write((uniprotid + "\t" + str(geneid) + '\t' + ','.join(
            name_id) + "\t" + ','.join(metabolite_id) + "\t" + ','.join(
            km_value) + "\t" + ",".join(chebids))+'\n')
        file.close()
