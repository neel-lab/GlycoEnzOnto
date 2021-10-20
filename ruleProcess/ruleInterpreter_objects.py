import pandas as pd
from tokenClasses_lex import *
from ruleInterpreter import *
import traceback
import glypy
from glypy import *
from dataclasses import dataclass

import sys
sys.path.insert(0,'../glycan_rdf')
from glycan_structure_ontology import GlycoCTProcessor
finishedGlycogenes=pd.read_csv('../finishedGlycogenes.tsv',sep='\t',index_col=None)
#Make dict:
objDict=dict()
for _,r in finishedGlycogenes.iterrows():
    if r['Rules']=='no reaction':
        continue
    try:
        objDict[r['geneName']]=reactionRule(lexer(r['Rules']))
    except:
        objDict[r['geneName']]=None

############################
# Glypy's glycan dictionary:
############################
monosaccharide_dict={k:v.serialize() for k,v in glypy.monosaccharides.items()}

#Get the names of unprocessed glycogenes:
unproc_ggenes=[k for k,v in objDict.items() if v==None]
#Filter dictionary with only instantiated reactionRules:
objDict={k:v for k,v in objDict.items() if k not in unproc_ggenes}

############################################
# Typed Dictionaries for Reaction Instances:
############################################

@dataclass
class Entity:
    entity_list: list

    @classmethod
    def factory(cls,*args,**kwargs):
        return cls(*args,**kwargs)

@dataclass
class ReactionEntity(Entity):
    operation: str

    @classmethod
    def factory(cls,*args,**kwargs):
        return cls(*args,**kwargs)

@dataclass
class monoReactionEntity(ReactionEntity):
    mod_entity: list

    @classmethod
    def factory(cls,*args,**kwargs):
        return cls(*args,**kwargs)

@dataclass
class SubstitutionEntity:
    operation: str
    from_entity_list: list
    to_entity_list: list

    @classmethod
    def factory(cls,*args,**kwargs):
        return cls(*args,**kwargs)

@dataclass
class monoSubstitutionEntity(SubstitutionEntity):
    mono_entity: list

    @classmethod
    def factory(cls,*args,**kwargs):
        return cls(*args,**kwargs)



####################################
# Rule String to glycoCT Conversion:
####################################

class glycoCTConvert:

    def __init__(self,reactionRuleObj):
        #Gather rule components:
        self.ruleSets=reactionRuleObj.ruleSets

    def rctTokenProcess(self,token):
        '''
        Returns tuples which contain the reaction strings
        and the operations that were performed on them based
        on the reaction rules:
        '''
        #Get kind of reaction:
        if token.token.__name__ in ['additionToken','subtractionToken']:
            if token.token.__name__=='additionToken':
                monoEntity=token.product()
                opString='addition'
            elif token.token.__name__=='subtractionToken':
                monoEntity=token.substrate()
                opString='subtraction'
            return(ReactionEntity.factory(entity_list=monoEntity,operation=opString))
        elif token.token.__name__ in ['substitutionToken','reversibleToken']:
            monoEntityFrom=token.substrate()
            monoEntityTo=token.product()
            if token.token.__name__=='substitutionToken':
                opString='substitution'
            elif token.token.__name__=='reversibleToken':
                opString='reversible'
            return(SubstitutionEntity.factory(operation=opString,from_entity_list=monoEntityFrom,to_entity_list=monoEntityTo))

    def mono_rctTokenProcess(self,token):
        '''
        Processes monosaccharides which have reaction rules in them.
        '''
        if token.token.reactionToken.token.__name__ in ['additionToken','subtractionToken']:
            if token.token.reactionToken.token.__name__=='additionToken':
                monoEntity=token.substrate()
                modToken=token.token.reactionToken.token.ligand_token[0].product()
                opString='addition'
            elif token.token.reactionToken.token.__name__=='subtractionToken':
                monoEntity=token.product()
                modToken=token.token.reactionToken.token.ligand_token[0].product()
                opString='subtraction'
            return(monoReactionEntity.factory(operation=opString,mod_entity=modToken,entity_list=monoEntity))
        elif token.token.reactionToken.token.__name__ in ['substitutionToken','reversibleToken']:
            monoEntity=token.substrate()
            fromModEntity=token.token.reactionToken.token.from_ligand_token[0].product()
            toModEntity=token.token.reactionToken.token.to_ligand_token[0].product()
            monoEntity=re.sub(fromModEntity,'',monoEntity)
            if token.token.reactionToken.token.__name__=='substitutionToken':
                opString='substitution'
            elif token.token.reactionToken.token.__name__=='reversibleToken':
                opString='reversible'
            return(monoSubstitutionEntity.factory(operation=opString,from_entity=fromModEntity,to_entity=toModEntity,mono_entity=monoEntity))
            
    def structureExtract(self,ruleSet):
        '''
        Reads through each ruleSet and creates a list of tuples
        where the first entity is a tag specifying if an operation
        is to occur with the structure, and the second is the 
        monosaccharide entity.
        '''
        ruleObjList=[]
        for elt in ruleSet:
            if elt.substrate()!=elt.product():
                if elt.__name__=='reactionToken':
                    obj=self.rctTokenProcess(elt)
                elif elt.__name__=='entityToken':
                    obj=self.mono_rctTokenProcess(elt)
            else:
                obj=Entity(entity_list=elt.substrate())
            ruleObjList.append(obj)
        return(ruleObjList)

    def get_structureObj_lists(self):
        '''
        Create nested list of structure objects:
        '''
        return([self.structureExtract(l) for l in self.ruleSets])

    
    
    def buildRuleGraph(self,structList):
        '''
        Builds graph structures out of reaction rules
        '''

    def tokenizeIUPAC(self,res):
        '''
        Takes monosaccharide entities and splits them
        into expected tags required for glycoCT conversion.
        Returns a dictionary for glycoStructOnto instantiation.
        '''
        #First expect monosaccharide:
        monoPattern=re.compile(r'(?P<WildCard>\[?\.\.\.\]?)?(?P<Mono>[A-Za-z]+?)(?P<Modification>\d\,?[SP])?(?P<Linkage>\([ab\?][12\?]\-[\d\?]\))')
        #Pattern for just modifications:
        modPattern=re.compile(r'(?P<Position>\d)(?<Modification>\D)')
        #Linkage Pattern:
        linkPattern=re.compile(r'\((?P<anomerConfig>[ab])(?P<anomerLink>[12\?])\-(?P<linkCarbon>[\d\?])\)')

        #Try matching monosaccharide:
        if re.search(monoPattern,res) is not None:
            #Extract monosaccharide patterns matched:
            monoGroups=re.search(monoPattern,res).groupdict() 
            ### Look up GlycoCT representation of monosaccharide in GlyPy ###
            if monoGroups['Mono'] in monosaccharide_dict.keys():
                monoRep=glypy.monosaccharides[monoGroups['Mono']]
            elif monoGroups['WildCard'] is not None:
                monoRep=monoGroups['WildCard'] 

        monoMatcher=LexMatcher(entityDict['monosaccharides'],
                entityDict['Compartments'],
                entityDict['Modifications'],
                regex="\[?%s(\[%s\])*((?:\d\,|\,\d|\d|\<.+?\>)*)((?:\{\!?\,\d\}|\{\!?\d?\D+?\}|\{.+?\<?\-\>.+?\}))*%s*(\([ab\?][12\?]\-[\d\?]\))*\]?")
        #Monosaccharide matching:
        regex=re.compile()
        #regex="\[?%s(\[%s\])*((?:\d\,|\,\d|\d|\<.+?\>)*)((?:\{\!?\,\d\}|\{\!?\d?\D+?\}|\{.+?\<?\-\>.+?\}))*%s*(\([ab\?][12\?]\-[\d\?]\))*\]?")
        for (_,mono) in ruleTupList:
            pass
