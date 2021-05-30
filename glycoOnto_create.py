import os
import sys
import pandas as pd
import pickle as pkl
from owlready2 import *
import re
import pathlib
from pathlib import Path

#Script-specific imports:
from glycoOnto_classObjects import *
from glycoOnto_utils import *

#Load ontology
glycoOnto=get_ontology('./glycoOnto.rdf').load()

# Create class dictionary:
classDict={x.name:x for x in glycoOnto.classes()}

# Read finished glycoEnzymes:
glycoEnzDB_finished=pd.read_csv('./finishedGlycogenes.tsv',sep='\t',header=0)

### Define ontology terms to parse: ###
ontoTerms=['Pathway','Function']
annotationTerms=['Reactant','Product','Constraint']
### Reparse the table to only contain relevant fields ###
glycoEnzDB_finished=glycoEnzDB_finished[['geneName']+ontoTerms+annotationTerms]
#Fix "Pathway" and "Functions" to match with strings:
glycoEnzDB_finished['Pathway']=[re.sub('\ ','_',x) for x in glycoEnzDB_finished['Pathway']]
glycoEnzDB_finished['Function']=[re.sub('\ ','_',x) for x in glycoEnzDB_finished['Function']]

### Check for missing/incorrect classes in ontology:
for i,row in glycoEnzDB_finished.iterrows():
    print_class_errors(i,row,classDict,glycoOnto)

#Helper function to create instances:

# Main iteration:

def procInstance(row,classDict,ontoDict,annotDict,onto):
    # 1. Create processing classes from procDict
    rowDta=row.to_dict()
    print(rowDta['geneName'])
    # Create instance in ontology, use the "Function" class to instantiate it:
    funList=rowDta['Function'].split(',')
    inst=classDict[funList[0]](rowDta['geneName'])
    if len(funList)>1:
        for fun in funList[1:]:
            cls=isa_A(inst,classDict[fun])
            cls.associate()
    # Associate pathways with instance:
    pathList=rowDta['Pathway'].split(',')
    pathList,invalid_pathList=class_process(pathList,onto)
    if len(pathList)>0:
        for pth in pathList:
            cls=inPathway_A(inst,classDict[pth])
            cls.associate()
    #Add annotation information: 
    for k,v in annotDict.items():
        cls=v(inst,rowDta[k])
        cls.associate()

#Blow up all instances, start with clean slate:
[destroy_entity(x) for x in glycoOnto.individuals()]
# Add Glycogene instances to ontology:
with glycoOnto:
    for _,row in glycoEnzDB_finished.iterrows():
        procInstance(row,classDict,ontoDict,annotDict,glycoOnto)

# Save the ontology into new file:
outFile="glycoOnto_individuals_reactions.rdf"
ontoFile_loc=os.path.join(os.getcwd(),outFile)
if outFile not in os.listdir('./'):
    print('Creating \"{0}\" in working dir...'.format(outFile))
    Path(ontoFile_loc).touch()

glycoOnto.save(file=ontoFile_loc,format="rdfxml")
