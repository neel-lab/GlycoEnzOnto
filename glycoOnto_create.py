import os
import sys
import pandas as pd
import pickle as pkl
from owlready2 import *
import re
import pathlib
from pathlib import Path

#Load ontology
glycoOnto=get_ontology('./glycoOnto.rdf').load()

# Create class dictionary:
classDict={x.name:x for x in glycoOnto.classes()}

# Read finished glycoEnzymes:
#glycoEnzDB_finished=pd.read_csv('./finishedGlycogenes.tsv',sep='\t',header=None).iloc[:,[3,6,7]]

glycoEnzDB_finished=pd.read_csv('./finishedGlycogenes.tsv',sep='\t',header=0)

### Define ontology terms to parse: ###
ontoTerms=['Pathway','Function']
annotationTerms=['Reactant','Product','Constraint']
### Reparse the table to only contain relevant fields ###
glycoEnzDB_finished=glycoEnzDB_finished[['geneName']+ontoTerms+annotationTerms]

### Map data table terms to ontology properties:
# Association classes:
class inPathway_A:
    def __init__(self,inst,pth):
        self.inst=inst
        self.pth=pth
        self.datatype='object'
    def check_instance_of(self):
        # Determines if a gene is an instance of a class:
        return(self.pth in self.inst.inPathway)
    def add(self):
        self.inst.inPathway.append(self.pth)
    def associate(self):
        if self.check_instance_of() is False:
            self.add()

class isa_A:
    def __init__(self,inst,obj):
        self.inst=inst
        self.pth=obj
        self.datatype='object'
    def check_instance_of(self):
        # Determines if a gene is an instance of a class:
        return(self.pth in self.inst.is_a)
    def add(self):
        self.inst.is_a.append(self.pth)
    def associate(self):
        if self.check_instance_of() is False:
            self.add()

#Annotation classes
class hasSubstr_A:  
    def __init__(self,inst,substr_string):
        self.inst=inst
        self.string=substr_string
        self.datatype='string'
    def check_instance_of(self):
        # Determines if a gene is an instance of a class:
        return(self.string in self.inst.hasSubstrate)
    def add(self):
        self.inst.hasSubstrate.append(self.string)
    def associate(self):
        if self.check_instance_of() is False:
            self.add()
 
class hasProd_A:  
    def __init__(self,inst,prod_string):
        self.inst=inst
        self.string=prod_string
        self.datatype='string'
    def check_instance_of(self):
        # Determines if a gene is an instance of a class:
        return(self.string in self.inst.hasProduct)
    def add(self):
        self.inst.hasProduct.append(self.string)
    def associate(self):
        if self.check_instance_of() is False:
            self.add()

class hasConstr_A:  
    def __init__(self,inst,constr_string):
        self.inst=inst
        self.string=constr_string
        self.datatype='string'
    def check_instance_of(self):
        # Determines if a gene is an instance of a class:
        return(self.string in self.inst.hasConstraint)
    def add(self):
        self.inst.hasConstraint.append(self.string)
    def associate(self):
        if self.check_instance_of() is False:
            self.add()

# Associate Functions with Dictionary Mappings:
ontoDict={
        'Pathway':inPathway_A,
        'Function':isa_A
}

annotDict={
        'Reactant':hasSubstr_A,
        'Product':hasProd_A,
        'Constraint':hasConstr_A
}

#Fix "Pathway" and "Functions" to match with strings:
glycoEnzDB_finished['Pathway']=[re.sub('\ ','_',x) for x in glycoEnzDB_finished['Pathway']]
glycoEnzDB_finished['Function']=[re.sub('\ ','_',x) for x in glycoEnzDB_finished['Function']]

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

def newClasses(geneName,cl_list,onto):
    # Creates a list of classes that haven't been associated 
    # with the gene being processed:
    return([c for c in cl_list if not check_instance_of(geneName,c,onto)])

def add_classes(geneInst,classDict,classList):
    # Adds classes to existing instance:
    [geneInst.is_a.append(classDict[x]) for x in classList]

#def add_instance(geneName,classList,classDict,onto):
#    # Adds instance to ontology provided a list of classes:
#    #If gene is new instantiate it:
#    print(classList[0])
#    inst=classDict[classList[0]](geneName)
#    # If more than one class, add to ontology
#    if len(classList)>1: add_classes(inst,classDict,classList[1:])

def add_instance(geneName,firstClass):
    # Adds instance to ontology provided a list of classes:
    #If gene is new instantiate it:
    inst=firstClass(geneName)
    return(inst)


# Main iteration:

def procInstance(row,classDict,ontoDict,annotDict,onto):
    # 1. Create processing classes from procDict
    rowDta=row.to_dict()
    print(rowDta['geneName'])
    # Check if instance present in ontology:
    inst=get_inst(rowDta['geneName'],onto) 
    #Add class information:
    for k,v in ontoDict.items():
        #Parse valid classes:
        classList=rowDta[k].split(',')
        classList,invalid_classList=class_process(classList,onto)
        if len(invalid_classList)>0:
            print('Found invalid classes associated with {0}: {1}'.format(geneName,invalid_classList))
        if len(classList)==0:
            print('No valid classes for {0}'.format(rowDta['geneName']))
            continue
        for c in classList:
            print(c)
            obj=classDict[c]
            if inst is None:
                #If gene is new instantiate it:
                inst=add_instance(rowDta['geneName'],obj)
            else:
                cls=v(inst,obj)
                cls.associate()
    #Add annotation information: 
    for k,v in annotDict.items():
        cls=v(inst,rowDta[k])
        cls.associate()

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
