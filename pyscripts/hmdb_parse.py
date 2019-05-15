# -*- coding: utf-8 -*-
"""
Created on Sun Sep  2 16:15:51 2018
HMDB parser

@author: Scott Campit
"""

import xml.etree.cElementTree as et

path = r'C:\Users\scampit\Desktop\Kinetic Model\Data\hmdb'
fil = r'hmdb_metabolites.xml'

nam = {'xmlns':"http://www.hmdb.ca"}
tree = et.parse(path+'\\'+fil)
root = tree.getroot()

def main():
    '''

    '''
    xmltree = et.parse(path+'\\'+fil)
    col = ['Metabolite', 'Tissues', 'Pathway', 'KEGG ID', 'Concentrations', 'Avg Conc', 'No. of Patients', 'Disease State']        
    df = pd.DataFrame(columns=col)
    
    for node in xmltree.getroot():
        

for child in root.findall 

for child in root.findall('./metabolite/normal_concentrations', nam):
    print(child.tag)
    
main()