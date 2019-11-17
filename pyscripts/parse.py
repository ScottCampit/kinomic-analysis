"""
Map metabolic model
@author: Scott Campit
"""
import re

from bs4 import BeautifulSoup as bs
import xml.etree.cElementTree as et

import numpy as np
import pandas as pd


def getGenesFromXML(file):
    """
    """

    dfcols = ['BiGG', 'CCDS', 'Entrez', 'OMIM', 'RefSeq', 'Synonym', 'SBO']
    rows = []

    xmlFile = et.parse(file)
    root = xmlFile.getroot()

    nameSpace = {
        'xmlns': "http://www.sbml.org/sbml/level2/version4",
        'html': "http://www.w3.org/1999/xhtml",
        'rdf': "http://www.w3.org/1999/02/22-rdf-syntax-ns#",
        'bqbiol': "http://biomodels.net/biology-qualifiers/",
        'xmlns': "http://www.sbml.org/sbml/level3/version1/core",
        'group': "http://www.sbml.org/sbml/level3/version1/groups/version1",
        'fbc': "http://www.sbml.org/sbml/level3/version1/fbc/version2"
        }

    """
    # Get Gene ID maps
    for gprs in root.findall(
            './xmlns:model/fbc:listOfGeneProducts/fbc:geneProduct/xmlns:annotation/rdf:RDF/rdf:Description', nameSpace):

        bigg = list(gprs.attrib.values())[0]
        bigg = re.sub('#G_', '', bigg)

        for ids in root.findall(
                    './xmlns:model/fbc:listOfGeneProducts/fbc:geneProduct/xmlns:annotation/rdf:RDF/rdf:Description/bqbiol:is/rdf:Bag/rdf:li', nameSpace):
            html = list(ids.attrib.values())[0]
            idStr = re.sub('http://identifiers.org/', '', html)
            type = idStr.split('/')[0]
            id = idStr.split('/')[1]
    """
    # Reaction subsystem map
    #for subsys in root.findall('./xmlns:model/group:listOfGroups/group:groups', nameSpace):
    #print(subsys)
    #subsystem = list(subsys.attrb.values())

    for reaction in root.findall('./xmlns:model/xmlns:listOfReactions/xmlns:reaction', nameSpace):
        rxn = reaction.attrib
        biggRxn = rxn.get('id', None)
        rxnName = rxn.get('name', None)
        rev = rxn.get('reversible', None)

        for gpa in root.findall('./xmlns:model/xmlns:listOfReactions/xmlns:reaction/fbc:geneProductAssociation/xmlns:fbc/fbc:geneProductRef', nameSpace):
            biggRxn = list(gpa.attrib)
            print(biggRxn)

        #for ele in reaction:
        #print(ele.attrib)
        #biggReaction = reaction.tag['id']
        #print(biggReaction)
        #print(reaction.tags)


file = r'./../Data/models/RECON1.xml'
getGenesFromXML(file)


def getGenesFromJSON(file):
      """
"""
      df = pd.read_json(file, orient='values')
      return df

      #df = getGenesFromJSON(r'./../Data/models/RECON1.json')
      #print(df)
