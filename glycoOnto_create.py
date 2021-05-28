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

def add_instance(geneName,cl,classDict):
    # Adds instance to ontology provided a list of classes:
    #If gene is new instantiate it:
    print(cl)
    inst=classDict[cl](geneName)
    return(inst)

#def add_instance(geneName,firstClass):
#    # Adds instance to ontology provided a list of classes:
#    #If gene is new instantiate it:
#    inst=firstClass(geneName)
#    return(inst)

#def add_instance(geneName,firstClass):
#    # Adds instance to ontology provided a list of classes:
#    #If gene is new instantiate it:
#    inst=firstClass(geneName)
#    return(inst)


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
