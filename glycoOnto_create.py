import os
import sys
import pandas as pd
import pickle as pkl
from owlready2 import *
import re
import pathlib
from pathlib import Path

#Load ontology
glycoOnto=get_ontology('./glycoOnto.owl').load()

# Create class dictionary:
classDict={x.name:x for x in glycoOnto.classes()}

# Read finished glycoEnzymes:
glycoEnzDB_finished=pd.read_csv('./finishedGlycogenes.tsv',sep='\t',header=None).iloc[:,[3,6,7]]
glycoEnzDB_finished.columns=['geneName','Pathway','Function']
#Fix "Pathway" and "Functions" to match with strings:
glycoEnzDB_finished['Pathway']=[re.sub('\ ','_',x) for x in glycoEnzDB_finished['Pathway']]

#Helper function to create instances:

def get_inst(geneName,onto):
    # Checks if instance is already present in the ontology
    # Returns object corresponding to instance if present, 
    #   otherwise returns None
    ind=[x for x in onto.individuals() if x.name==geneName]
    if len(ind)==0:
        return(None)
    else:
        return(ind[0])

def check_if_class(cl,onto):
    # Determines if class is valid based on ontology
    return(cl in [x.name for x in onto.classes()])

def class_process(classList,onto):
    # Creates two lists: one with valid classes, one with invalid ones:
    validClasses=[c for c in classList if check_if_class(c,onto)]
    invalidClasses=[c for c in classList if c not in classList]
    return validClasses,invalidClasses

def check_instance_of(geneName,cl,onto):
    # Determines if a gene is an instance of a class:
    inst=get_inst(geneName,onto)
    return(cl in [x.name for x in inst.is_instance_of])

def newClasses(geneName,cl_list,onto):
    # Creates a list of classes that haven't been associated 
    # with the gene being processed:
    return([c for c in cl_list if not check_instance_of(geneName,c,onto)])

def add_classes(geneInst,classDict,classList):
    # Adds classes to existing instance:
    [geneInst.is_a.append(classDict[x]) for x in classList]

def add_instance(geneName,classList,classDict,onto):
    # Adds instance to ontology provided a list of classes:
    #If gene is new instantiate it:
    print(classList[0])
    inst=classDict[classList[0]](geneName)
    # If more than one class, add to ontology
    if len(classList)>1: add_classes(inst,classDict,classList[1:])

# Main iteration:
def create_instance(geneName,classList,classDict,onto):
    # Find valid and invalid classes:
    classList,invalid_classList=class_process(classList,onto)
    if len(invalid_classList)>0:
        print('Found invalid classes associated with {0}: {1}'.format(geneName,invalid_classList))
    if len(classList)==0:
        print('No valid classes for {0}'.format(geneName))
        return(None)
    # Check if instance present in ontology:
    inst=get_inst(geneName,onto) 
    if inst is not None:
        print('{0} already present in ontology'.format(geneName))
        # Filter out already-existing relationships:
        classList=newClasses(geneName,classList,onto)
        if len(classList)==0:
            print('\tNothing to add')
            return(None)
        add_classes(inst,classDict,classList)
    #If gene is new instantiate it:
    add_instance(geneName,classList,classDict,onto)

# END HELPER FUNCTIONS

# Add Glycogene instances to ontology:
for _,row in glycoEnzDB_finished.iterrows():
    # Parse gene, pathway and function information from row:
    gene=row['geneName']
    pth=row['Pathway'].split(',')
    func=row['Function'].split(',')
    #Perform functions below under the "glycoOnto" namespace
    with glycoOnto:
        # Process pathway information:
        print('Processing pathways:')
        create_instance(gene,pth,classDict,glycoOnto)
        # Process function information:
        print('Processing functions:')
        create_instance(gene,func,classDict,glycoOnto)

# Save the ontology into new file:
outFile="glycoOnto_individuals.rdf"
ontoFile_loc=os.path.join(os.getcwd(),outFile)
if outFile not in os.listdir('./'):
    print('Creating \"{0}\" in working dir...'.format(outFile))
    Path(ontoFile_loc).touch()

glycoOnto.save(file=ontoFile_loc,format="rdfxml")
