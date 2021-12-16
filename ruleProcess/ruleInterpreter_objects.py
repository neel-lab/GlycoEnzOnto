import pandas as pd
from tokenClasses_lex import *
from ruleInterpreter import *
import traceback
import glypy
from glypy import *
from dataclasses import dataclass
from itertools import tee,product

def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)

def test_wrapper(rule):
    try:
        res=lexer(rule)
    except Exception as exc:
        return(None)
    try:
        gct=reactionRule(res)
    except Exception as exc:
        return(None)
    return(gct)


import sys
sys.path.insert(0,'../glycan_rdf')
from glycan_structure_ontology import GlycoCTProcessor
finishedGlycogenes=pd.read_csv('../finishedGlycogenes.tsv',sep='\t',index_col=None)
import owlready2
from owlready2 import *
glycoStructOnto=get_ontology('../glycan_rdf/glycoStructOnto.rdf').load()
glycoStructOnto.base_iri='http://mzjava.expasy.org/glycan/'

#Make dict:
objDict=dict()
for _,r in finishedGlycogenes.iterrows():
    if r['Rules']=='no reaction' or test_wrapper(r['Rules']) is None:
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

class Entity:
    def __init__(self,entity_list):
        self.entity_list=entity_list
        self.tokenDicts=[self.tokenizer(x) for x in self.entity_list]

    @classmethod
    def factory(cls,*args,**kwargs):
        return cls(*args,**kwargs)

    #Tokenizer:
    def tokenizer(self,ent):
        '''
        Returns a dictionary of Entity features.
        Return dictionary keys:
        - Monosaccharide label
          "One or more monosaccharides (on a branch) is 
          noted if the entity is not a monosaccharide.
        - Anomeric carbon number
        - Linkage carbon number
        - Modifications
          If modifications are detected, a new entity is 
          created
        - Modifications
        '''
        t_dict={
            'monosaccharide':self.mono_label(ent),
            'anomer_carbon':self.anomer_link(ent),
            'anomer_config':self.anomer_config(ent),
            'link_carbon':self.linkage_carbon(ent),
            'mods':self.mod_type(ent)
        }
        return(t_dict)
    
    def parseWrap(fun):
        def _wrap(self,ent):
            #First expect monosaccharide:
            sch=re.compile(r'(?P<WildCard>\[?\.\.\.\]?)?(?P<Mono>[A-Za-z]+?)(?P<Modification>\d(\,\d)*[SP])?(?P<Linkage>\([ab\?][12\?]\-[\d\?]\))')
            sch_res=re.search(sch,ent)
            #Apply "getter" functions to retrieve structural features:
            ft=fun(sch_res) if sch_res is not None else None
            return(ft)
        return(_wrap)

    @parseWrap
    def mono_label(sch_res):
        return(sch_res['Mono'])
    
    @parseWrap
    def mod_type(sch_res):
        return(sch_res['Modification'])

    def linkGetterWrap(fun):
        def _wrap(sch_res):
            #Linkage Pattern:
            linkPattern=re.compile(r'\((?P<anomerConfig>[ab])(?P<anomerLink>[12\?])\-(?P<linkCarbon>[\d\?])\)')
            linkData=re.search(linkPattern,sch_res['Linkage'])
            #Apply getter function:
            ft=fun(linkData)
            return(ft)
        return(_wrap)

    @parseWrap
    @linkGetterWrap
    def anomer_link(sch_res):
        return(sch_res['anomerLink'])

    @parseWrap
    @linkGetterWrap
    def linkage_carbon(sch_res):
        return(sch_res['linkCarbon'])

    @parseWrap
    @linkGetterWrap
    def anomer_config(sch_res):
        return(sch_res['anomerConfig'])


class ReactionEntity(Entity):
    def __init__(self,operation,**kwargs):
        self.operation=operation
        super().__init__(**kwargs)

    @classmethod
    def factory(cls,*args,**kwargs):
        return cls(*args,**kwargs)

class monoReactionEntity(ReactionEntity):
    def __init__(self,mod_entity,**kwargs):
        self.mod_entity=mod_entity
        super().__init__(**kwargs)

    @classmethod
    def factory(cls,*args,**kwargs):
        return cls(*args,**kwargs)

class SubstitutionEntity:
    def __init__(self,operation,from_entity_list,to_entity_list):
        self.operation=operation
        self.from_entity_list=from_entity_list
        self.to_entity_list=to_entity_list

    @classmethod
    def factory(cls,*args,**kwargs):
        return cls(*args,**kwargs)

class monoSubstitutionEntity(SubstitutionEntity):
    def __init__(self,mono_entity,**kwargs):
        self.mono_entity=mono_entity
        super().__init__(**kwargs)

    @classmethod
    def factory(cls,*args,**kwargs):
        return cls(*args,**kwargs)

####################################
# Rule String to glycoCT Conversion:
####################################

class glycoCTConvert:

    def __init__(self,reactionRule_inst,objOnto):
        #Gather rule components:
        self.ruleSets=reactionRule_inst.ruleSets
        self.pairLists=self.entity_pairLists()
        self.entityLists=[list(reversed([self.structureExtract(x) for x in lst])) for lst in self.ruleSets]
        self.objOnto=objOnto

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
            
    def get_branching(self,elt):
        if elt.__name__=='reactionToken':
            eltToken=elt.token.ligand_token[0]
        else:
            eltToken=elt
        if hasattr(eltToken.token,'branching'):
            elt_branch=eltToken.token.branching
        else:
            elt_branch={'leftBracket':False,'rightBracket':False}
        return(elt_branch)

    def structureExtract(self,elt):
        '''
        Reads through each ruleSet and creates a list of tuples
        where the first entity is a tag specifying if an operation
        is to occur with the structure, and the second is the 
        monosaccharide entity.
        '''
        if elt.substrate()!=elt.product():
            if elt.__name__=='reactionToken':
                obj=self.rctTokenProcess(elt)
            elif elt.__name__=='entityToken':
                obj=self.mono_rctTokenProcess(elt)
        else:
            obj=Entity(entity_list=elt.substrate())
        return(obj)

    def connectList(self,ruleSet):
        '''
        Creates connectivity list for elements in a rule set 
        '''
        connectList=[]
        returnInd=None
        ruleSet_reducing=list(reversed(ruleSet))
        for i,((x,p),(y,c)) in enumerate(pairwise(enumerate(ruleSet_reducing))):
            p_branch,c_branch=self.get_branching(p),self.get_branching(c)
            p_elt,c_elt=self.structureExtract(p),self.structureExtract(c)
            if p_branch['leftBracket']: 
                continue
            connectList.append(((x,p_elt),(y,c_elt)))
            if c_branch['leftBracket']:
                if c_branch['rightBracket']:
                    #Make a connection to both the child
                    # and the next following element:
                    lt=self.structureExtract(ruleSet_reducing[y+1])
                    connectList.append(((x,p_elt),(y+1,lt)))
                else:
                    lt_p=self.structureExtract(ruleSet_reducing[returnInd]);lt_c=self.structureExtract(ruleSet_reducing[y+1])
                    connectList.append(((returnInd,lt_p),(y+1,lt_c)))
            elif c_branch['rightBracket']:
                returnInd=y-1
        return(connectList)

    def connectLists(self):
        return([self.connectList(x) for x in self.ruleSets])
    
    def entity_pairLists(self):
        connect_lists=self.connectLists()
        pairLists=[[(x[0][0],x[1][0]) for x in lst] for lst in connect_lists]
        return(pairLists)

    def eltProd(self,entityList):
        return([list(prod(*[x.tokenDicts for x in lst])) for lst in self.entityLists])

    #############################
    # Ontology-Related Functions:
    #############################
    def monoProcess(self,mono_string):
        '''
        Creates an instance of a monosaccharide in the given ontology
        of the class:
        '''
        #IUPAC to GlycoCT:
        mono_glycoCT=[v.serialize() for k,v in glypy.monosaccharides.items() if k==mono_string]
        if not (not mono_glycoCT):
            #Check if there are substituents on it:
            if re.search('LIN',mono_glycoCT[0]) is not None:
                #Split the monosaccharide and substituent parts:
                monoPart,modPart=re.search('(\S+?\-){3}\d\:\d',mono_glycoCT[0]).group(),re.search('2s\:[a-z]+',mono_glycoCT[0]).group()
                monoPart,modPart=''.join(['RES ',monoPart]),''.join(['RES ',re.sub('2','1',modPart)])
                mono_inst=[x for x in self.objOnto.get_children_of(self.objOnto.monosaccharide) if x.label[0]==monoPart]
                mod_inst=[x for x in self.objOnto.get_children_of(self.objOnto.substituent) if x.label[0]==modPart]
                if not (not mono_inst) and not (not mod_inst):
                    return({'mono':mono_inst[0],'mod':mod_inst[0]})
            else:
                monoPart=mono_glycoCT[0]
                mono_inst=[x for x in self.objOnto.get_children_of(self.objOnto.monosaccharide) if x.label[0]==monoPart]
                if not (not mono_inst):
                    return({'mono':mono_inst[0],'mod':None})
        else:
            return(None)

    def mono_mono_connect(self,fromRes,toRes,anomerLink,anomericity,linkageNumber):
        fromRes.is_GlycosidicLinkage.append(toRes)
        #Handle anomer connection:
        if anomerLink==1:
            fromRes.has_anomerCarbon_1.append(toRes)
        elif anomerLink==2:
            fromRes.has_anomerCarbon_2.append(toRes)

        if linkageNumber==1:
            fromRes.has_linkedCarbon_1.append(toRes)
        elif linkageNumber==2:
            fromRes.has_linkedCarbon_2.append(toRes)
        elif linkageNumber==3:
            fromRes.has_linkedCarbon_3.append(toRes)
        elif linkageNumber==4:
            fromRes.has_linkedCarbon_4.append(toRes)
        elif linkageNumber==5:
            fromRes.has_linkedCarbon_5.append(toRes)
        elif linkageNumber==6:
            fromRes.has_linkedCarbon_6.append(toRes)
        elif linkageNumber==8:
            fromRes.has_linkedCarbon_8.append(toRes)

        if anomericity=='a':
            fromRes.has_anomericConnection_alpha.append(toRes)
        elif anomericity=='b':
            fromRes.has_anomericConnection_beta.append(toRes)

    def mono_subst_connect(self,fromRes,toRes,linkageNumber):
        fromRes.is_SubstituentLinkage.append(toRes)
        if linkageNumber==1:
            fromRes.has_linkedCarbon_1.append(toRes)
        elif linkageNumber==2:
            fromRes.has_linkedCarbon_2.append(toRes)
        elif linkageNumber==3:
            fromRes.has_linkedCarbon_3.append(toRes)
        elif linkageNumber==4:
            fromRes.has_linkedCarbon_4.append(toRes)
        elif linkageNumber==5:
            fromRes.has_linkedCarbon_5.append(toRes)
        elif linkageNumber==6:
            fromRes.has_linkedCarbon_6.append(toRes)
        elif linkageNumber==8:
            fromRes.has_linkedCarbon_8.append(toRes)

    def buildRuleGraph(self,onto):
        '''
        Builds graph structures out of reaction rules
        '''
        for eList,pList in zip(self.entityLists,self.pairLists):
            #Create all permutations of glycan graph patterns:
            eList_prods=self.eltProd(eList)
            for elp in eList_prods:
                #Current rule instance:
                rl=onto.reactionRule
                for p in pList:
                    #Parent-Child Relationship
                    parent_ind,child_ind=p[0],p[1]
                
