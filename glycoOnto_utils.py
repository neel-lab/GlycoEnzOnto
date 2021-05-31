import os
import sys
import pandas as pd
import pickle as pkl
from owlready2 import *
import re
import pathlib
from pathlib import Path

def get_inst(geneName,onto):
    # Checks if instance is already present in the ontology
    # Returns object corresponding to instance if present, 
    #   otherwise returns None
    ind=[x for x in onto.individuals() if x.name==geneName]
    if len(ind)==0:
        return(None)
    else:
        return(ind[0])

def get_class(nme,onto):
    #Retrieves class object if name exists:
    classList=[x.name for x in onto.classes()]
    if nme in classList:
        return([x for x in onto.classes() if x.name==nme][0])
    else:
        return(None)

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
