import os
import sys
import pandas as pd
import pickle as pkl
from owlready2 import *
import re
import pathlib
from pathlib import Path

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


