import re
from tokenMatchers import *


######################
# Reaction Rule Tokens
######################

class reactionToken:

    def __init__(self,inputString):
        '''
        Token that describes any sort of reaction.
        These tokens must start with the '{' character, which
        could mean the beginning of one of the following:
        
        - Addition/Subtraction: { }/{! }
        - Substitution/Reversible reaction: { -> }/{ <-> }
        '''
        self.inputString=inputString
        #Get location of reaction match
        self.match=self.matchFun()
        if self.match is not None:
            #Get reaction type:
            self.token=self.get_reaction_type()
        else:
            self.token=None

    def __repr__(self):
        return(self.get_reaction_type().__repr__())

    def detectFun(self):
        '''
        Reaction Detection Method
        Logic:
           reactionMatcher
        '''
        #Helper function:
        isReaction=reactionMatcher(self.inputString)
        return(isReaction)
    
    def matchFun(self):
        '''
        Matcher function for reaction instances.
        If the reactionMatcher passes, get contents
        in between { }.
        '''
        if self.detectFun():
            return(reactionMatcher(self.inputString,presence=False))
        else:
            return(None)
    
    def get_reaction_type(self):
        '''
        Returns token of the kind of reaction matched.

        Internal methods to these token classes will run
        once instantiated.
        '''
        reactionMatch=self.matchFun()
        if reactionMatch is not None:
            reactionText=reactionMatch.group()
            if additionMatcher(reactionText) and not substitutionMatcher(reactionText):
                return(additionToken(reactionText))
            elif subtractionMatcher(reactionText):
                return(subtractionToken(reactionText))
            elif substitutionMatcher(reactionText):
                return(substitutionToken(reactionText))
            elif reversibleMatcher(reactionText):
                return(reversibleToken(reactionText))


class additionToken(reactionToken):

    def __init__(self,inputString):
        '''
        This class describes addition reactions which
        have the form { }.
        
        The input string to this class is passed
        by the reactionToken __init__ method.

        The input string should be of the form {.+?}
        '''
        self.__name__='additionToken'
        self.inputString=inputString
        ligand_string=self.get_ligand()
        #Call the entityToken class to identify 
        # the class of the ligand:
        self.ligand_token=entityToken(ligand_string)

    def __repr__(self):
        try:
            ligand=self.get_ligand()
        except:
            return('PARSE ERROR')
        return('An addition of %s' %(ligand))

    def get_ligand(self):
        ligand=re.search('\{(.*?)\}',self.inputString).groups()[0]
        return(ligand)

    def substrate(self):
        '''
        Substrate method for token.
        '''
        return('')

    def product(self):
        '''
        Product method for token.
        '''
        return(self.get_ligand())

class subtractionToken(reactionToken):

    def __init__(self,inputString):
        '''
        This class describes subtraction reactions which
        have the form {! }.
        
        The input string to this class is passed
        by the reactionToken __init__ method.

        The input string should be of the form {.+?}
        '''
        self.__name__='subtractionToken'
        self.inputString=inputString
        self.ligand=self.get_ligand()
        #Call the entityToken class to identify 
        # the class of the ligand:
        self.ligand_token=entityToken(ligand)

    def __repr__(self):
        try:
            ligand=self.get_ligand()
        except:
            return('PARSE ERROR')
        return('A removal of %s' %(ligand))

    def get_ligand(self):
        ligand=re.search('\{\!(.+?)\}',self.inputString).groups()[0]
        return(ligand)

    def substrate(self):
        return(self.get_substrate())

    def product(self):
        return('')


class substitutionToken(reactionToken):

    def __init__(self,inputString):
        '''
        This class describes substitution reactions which
        have the form {A -> B}.
        
        The input string to this class is passed
        by the reactionToken __init__ method.
        '''
        self.__name__='substitutionToken'
        self.inputString=inputString
        from_ligand,to_ligand=self.get_from_to()
        #Call the entityToken class to identify 
        # the class of the ligand:
        self.from_ligand_token,self.to_ligand_token=(entityToken(from_ligand),entityToken(to_ligand))

    def __repr__(self):
        try:
            frm,to=self.get_from_to()
        except:
            return('PARSE ERROR')
        return('A conversion from %s to %s' %(frm,to))

    def get_from_to(self):
        frm,to=re.search('\{(.+?)\-\>(.+?)\}',self.inputString).groups()
        return(frm,to)

    def substrate(self):
        frm,_=self.get_from_to()
        return(frm)

    def product(self):
        _,to=self.get_from_to()
        return(to)
    

class reversibleToken(reactionToken):

    def __init__(self,string):
        '''
        This class describes reversible reactions which
        have the form {A <-> B}.
        
        The input string to this class is passed
        by the reactionToken __init__ method.
        '''
        self.__name__='substitutionToken'
        self.inputString=inputString
        from_ligand,to_ligand=self.get_from_to()
        #Call the entityToken class to identify 
        # the class of the ligand:
        self.from_ligand_token,self.to_ligand_token=(entityToken(from_ligand),entityToken(to_ligand))

    def __repr__(self):
        try:
            frm,to=self.get_from_to()
        except:
            return('PARSE ERROR')
        return('A reversible conversion of %s to %s' %(frm,to))

    def get_from_to(self):
        frm,to=re.search('\{(.+?)\<\-\>(.+?)\}',self.inputString).groups()
        return(frm,to)

    def substrate(self):
        frm,_=self.get_from_to()
        return(frm)

    def product(self):
        _,to=self.get_from_to()
        return(to)

###################
# Constraint Tokens
###################
    
class constraintToken:
    
    def __init__(self,inputString):
        '''
        Constraint tokens indicate the beginning of 
        quantitative and structural constraints of 
        rule strings.

        This superclass returns the kind of constraint
        as an attribute as well as its own recursive parsing
        routine which stores a local set of entities which
        are subject to the constraint
        '''
        self.inputString=inputString
        self.match=self.matchFun()
        if self.match is not None:
            self.constr=self.get_constraint_type()
            #Recursively parse entities until 
            # end of string is reached, and return
            # entity list.
            #If there are numerical constraints, return
            # a quant_token, else return none
            idx=self.match.end()
            self.entities,self.end_idx,self.quant_token=self.entity_parse(idx)
        else:
            self.constr=None

    def detectFun(self):
        '''
        Detection method for constraintTokens:
        
        Must be one of:
        - n : Marks beginning of numerical constraint
        - @ : Attachment constraint
        - ! : Logical Negation Constraint
        '''
        isQuantityStart=quantityStartMatcher(self.inputString)
        isAttach=attachRuleMatcher(self.inputString)
        isNegation=negationRuleMatcher(self.inputString)
        return(isQuantityStart or isAttach or isNegation)

    def matchFun(self):
        '''
        Returns string matching object where constraint
        token was detected.
        '''
        if self.detectFun():
            if quantityStartMatcher(self.inputString):
                return(quantityStartMatcher(self.inputString,presence=False))
            elif attachRuleMatcher(self.inputString):
                return(attachRuleMatcher(self.inputString,presence=False))
            elif negationRuleMatcher(self.inputString):
                return(negationRuleMatcher(self.inputString,presence=False))
        else:
            return(None)

    def get_constraint_type(self):
        '''
        Returns string matching object where constraint
        token was detected.
        '''
        if self.detectFun():
            if quantityStartMatcher(self.inputString):
                return(quantityRule_token(self.inputString,presence=False))
            elif attachRuleMatcher(self.inputString):
                return(attachRule_token(self.inputString,presence=False))
            elif negationRuleMatcher(self.inputString):
                return(negationRule_token(self.inputString,presence=False))

    def entity_parse(self,idx,entityList=[]):
        '''
        Use the current index and the input string
        to continue parsing entities until a stopping
        condition is met.

        Encounters a logical tag or the end of the string
        when parsing ends:
        '''
        #Terminate if at the end of inputString
        if idx==len(self.inputString):
            quantToken=None
            return(entityList,idx,quantToken)
        stringFrag=self.inputString[idx:]
        #Instantiate logical, entity, and quantifier tokens
        # given the string fragment
        lg=logicalToken(stringFrag)
        ent=entityToken(stringFrag)
        qtr=quantifierToken(stringFrag)
        #Terminate if the next token is a 
        # logical token, or if quantifier text detected:
        if lg.detectFun():
            quantToken=None
            return(entityList,idx,quantToken)
        elif quantifierToken.detectFun():
            #Instantiates a quantity token
            quantToken=quantifierToken(stringFrag)
            idx+=mtch.end()
            return(entityList,idx,quantToken)
        elif ent.detectFun():
            #Append entity to list, find its matching location,
            # then increment search index by the length of the match:
            entityList.append(ent)
            mtch=ent.matchFun()
            idx+=mtch.end()
            #Run this function again:
            return(self.entity_parse(idx,entityList=entityList))
        else:
            raise Exception('Couldn\'t find valid entity or stopping condition, most likely due to malformed rule.')


class quantityRule_token(constraintToken):

    def __init__(self,string):
        self.__name__='quantityRule_token'
        self.inputString=string

    def __repr__(self):
        return('Quantitiy of')

    def detectFun(self):
        '''
        Detect method for token
        Logic:
           quantityStartMatcher
        '''
        isQuantityStart=quantityStartMatcher(self.inputString,presence=True)
        return(isQuantityStart)

    def matchFun(self):
        '''
        Returns match object.
        If "detectFun" returns True, return a match object
        using a matcherClass.
        '''
        if self.detectFun():
            return(quantityStartMatcher(self.inputString))
        else:
            return(None)

class attachRule_token(constraintToken):

    def __init__(self,string):
        self.__name__='attachRule_token'
        self.inputString=string

    def __repr__(self):
        return('An attachment constraint')
        
    def detectFun(self):
        '''
        Detect method for token
        Logic:
           attachRuleMatcher
        '''
        isAttachConstraint=attachRuleMatcher(self.inputString,presence=True)
        return(isAttachConstraint)

    def matchFun(self):
        '''
        Returns match object.
        If "detectFun" returns True, return a match object
        using a matcherClass.
        '''
        if self.detectFun():
            return(attachRuleMatcher(self.inputString))
        else:
            return(None)

class negationRule_token(constraintToken):

    def __init__(self,string):
        self.__name__='negationRule_token'
        self.inputString=string

    def __repr__(self):
        return('The absence of')

    def detectFun(self):
        '''
        Detect method for token
        Logic:
           negationRuleMatcher
        '''
        isNegationConstraint=negationRuleMatcher(self.inputString,presence=True)
        return(isNegationConstraint)

    def matchFun(self):
        '''
        Returns match object.
        If "detectFun" returns True, return a match object
        using a matcherClass.
        '''
        if self.detectFun():
            return(negationRuleMatcher(self.inputString))
        else:
            return(None)


class quantifierToken(constraintToken):

    def __init__(self,string):
        self.__name__='quantifierToken'
        self.inputString=string
    
    def __repr__(self):
        return("Preceeding pattern matches %s times" %(self.inputString))
    
    def detectFun(self):
        '''
        Detect method for token
        Logic:
           quantityStartMatcher
        '''
        isQuantifier=quantifierMatcher(self.inputString,presence=True)
        return(isQuantifier)
    
    def matchFun(self):
        '''
        Returns match object.
        If "detectFun" returns True, return a match object
        using a matcherClass.
        '''
        if self.detectFun():
            return(quantifierMatcher(self.inputString))
        else:
            return(None)

    def logical_fun(self):
        qt,val=re.search('',self.inputString).groups()
        l_fun=lambda pat,s:re.findall(pat,s)
        if qt=="=":
            return lambda pat,s: len(l_fun(pat,s))==qt
        elif qt==">=":
            return lambda pat,s: len(l_fun(pat,s))>=qt
        elif qt=="<=":
            return lambda pat,s: len(l_fun(pat,s))<=qt
        elif qt==">":
            return lambda pat,s: len(l_fun(pat,s))>qt
        elif qt=="<":
            return lambda pat,s: len(l_fun(pat,s))<qt

###############
# Entity Tokens
###############

class entityToken:

    def __init__(self,inputString):
        '''
        Token that describes any sort of chemical entity.
        
        The following entities will be detected:
        - Monosaccharides
        - Sugar Nucleotides
        - Modifications (Sulfation, Phosphorylation)
        - Compartments ( cis/medial/trans Golgi, ER, Lysosome)
        - Aglyca (Dol-P-P, Asn, Ser/Thr,Cer)
        '''
        #Object String input:
        self.inputString=inputString
        self.entityDict={
                'monosaccharides':['GlcNAc','GlcN','GlcA','Glc','GalNAc','Gal','ManNAc','Man','Neu5Ac','Neu5Gc','Xyl','Fuc','IdoA','Kdn','Ribitol(P5-1)','Ribitol(P5-3)','Ribitol','Fruc'],
                'Nucleotides':['CMP','UDP','GDP'],
                'Modifications':['S','P'],
                'Substrates':['PAPS'],
                'Compartments':['ER','Golgi','Cytoplasm'],
                'Aglyca':['Dol-P-P', 'Asn','Ser/Thr','Ser','Thr','Cer','5-hydroxy-L-lysyl-','Lys','WXXW','PI'],
                'ProteinConstraints':['EGF','Cad']
        }

        #Utilities:
        #Monosaccharide Processing:
        self.monoMatcher=LexMatcher(self.entityDict['monosaccharides'],
                self.entityDict['Compartments'],
                self.entityDict['Modifications'],
                regex="\[?%s(\[%s\])*((?:\d\,|\,\d|\d|\<.+\>)*)((?:\{\!?\,\d\}|\{\!?\d?\D+?\}|\{.+?\<?\-\>.+?\}))*%s*(\([ab\?][12\?]\-[\d\?]\))*\]?")
        #Sugar Nucleotide Matching:
        self.nucleotideSugarMatcher=LexMatcher(self.entityDict['Nucleotides'],self.entityDict['monosaccharides'],regex="%s\-%s")
        #Modification Matcher
        self.modMatcher=LexMatcher(self.entityDict['Modifications'],regex="\d%s",front=True)
        #Aglycon Matcher:
        self.aglyconMatcher=LexMatcher([re.escape(x) for x in self.entityDict['Aglyca']],front=True)
        #Compartment Matcher:
        self.compartmentMatcher=LexMatcher(self.entityDict['Compartments'],regex='\[%s\]',front=True)
        #Transport Matcher:
        self.transportMatcher=matcherClass('\-\>')
        #Protein Constraints:
        self.proteinConstraintMatcher=LexMatcher(self.entityDict['ProteinConstraints'],regex="\[%s\]")

        ### Get token Type: ###
        if self.detectFun() is not None:
            self.match,self.token=self.detectFun()
        else:
            #Try detecting unknown entity:
            self.match=None
            self.token=None

    def __repr__(self):
        _,token=self.detectFun()
        return(token.__repr__())

    def mono_parse(self):
        '''
        Wrapper for the monoMatcher method.

        Returns a dictionary of terms if successfully matched.
        Otherwise returns None
        '''
        if self.monoMatcher(self.inputString):
            #Get the match object into a group:
            monoMatch=self.monoMatcher(self.inputString,presence=False)
            matchDict={
                'mono':monoMatch.groups()[0],
                'compartment':monoMatch.groups()[1],
                'modPos':monoMatch.groups()[3],
                'innerReaction':monoMatch.groups()[4],
                'modType':monoMatch.groups()[5],
                'monoLink':monoMatch.groups()[6]
            }
            return(matchDict)
        else:
            return None

    def detectMono(self):
        '''
        Method for detecting monosaccharide entities.
        Uses the "monoMatcher" method to detect valid monosaccharide
        entities and split them based on their components:

        "monoMatcher"s Regex Groups:
        Group 1 Monosaccharides: Detected from the "entityToken"s allowed monosaccharide list.
        Group 2 Transporter Markers: Detected from the "entityToken"s allowed compartment list.
        Group 3 Carbon positions of modifications "((?:\d\,|\,\d|\d|\<.+\>)*)" RETURNED AS STRING.
        Group 4 Presence of any internal reaction information "((?:\{\!?\,\d\}|\{\!?\d?\D+?\}))*".
        Group 5 Modification name (if present): Detected from "entityToken"'s allowed modification list.
        Group 6 Linkage information: (\([ab\?][12\?]\-[\d\?]\)).

        '''
        #Parse the monosaccharide components using
        # the wrapping function;
        mono_components=self.mono_parse()
        # If not a monosaccharide, return None:
        if mono_components is None:
            return(None)
        #Otherwise, create the monosaccharide token
        # and associate ancestor tokens when necessary
        #Compartment:
        if mono_components['compartment'] is not None:
            compartment_token=compartmentToken(mono_components['compartment'])
        else:
            compartment_token=None
        #Modifications:
        if len(mono_components['modPos'])>0:
            #Create a list with the modification positions
            modPos_lst=re.findall('(\d|\<.+>)',mono_components['modPos'])
            if mono_components['modType'] is None:
                #Infer the mod if it is within
                # the modpos:
                modInfer=[re.search('(\d)(\D)',x) for x in modPos_lst if re.search('(\d)(\D)',x) is not None]
                if len(modInfer)==0:
                    raise Exception("Malformed Monosaccharide: there are modification positions within the monosaccharide with no valid modifications found.  Check the input rule string")
                else:
                    mtch=modInfer[0]
                    mono_components['modType']=mtch.groups()[1]
            modTokens=[]
            for mp in modPos_lst:
                if '<' in mp:
                    modTokens.append(multiToken(mp))
                else:
                    modTokens.append(modToken(mp+mono_components['modType']))
        else:
            modTokens=None
        #Modification Reactions:
        if mono_components['innerReaction'] is not None:
            rct_token=reactionToken(mono_components['innerReaction'])
        else:
            rct_token=None
        #Return All Attributes:
        return(mono_components['mono'],mono_components['monoLink'],modTokens,rct_token,compartment_token)

    def find_unknown_token(self,idx=0):
        '''
        Scans for unknown token contents and returns 
        a match object.
        '''
        idx+=1
        inputString_frag=self.inputString[idx:]
        #Instantiate entity token, runs the whole process again
        ent=entityToken(inputString_frag)
        if ent.token is not None or idx==len(self.inputString):
            #Create match object
            mtch=re.search('.{%d}'%(idx),self.inputString)
            token=unknownToken(mtch.group())

        mtch=self.find_unknown_token(idx=0)
        token=unknownToken(unknownMtch.group())
        return(mtch,token)

    def detectFun(self):
        '''
        Main entity detection wrapper.
        Checks for following patterns in order:
        1. Monosaccharide-linkage: Gal(b1-4)
           - With Compartments
           - With modifications
           - With modification reactions
        2. Nucleotide-Sugar:
           - Must begin with any of the nucleotides in entityDict
        3. Aglyca: 
           - Must match exactly any of the Aglyca entries in entityDict.
        4. Other:
        '''
        
        #Monosaccharide detection:
        if self.detectMono() is not None:
            monoMatch=self.monoMatcher(self.inputString,presence=False)
            (mono,linkage,modTokens,rct_token,compartment_token)=self.detectMono()
            token=monoToken(mono,linkage,modTokens,compartment_token,rct_token)
            return(monoMatch,token)

        #Nucleotide Sugar:
        elif self.nucleotideSugarMatcher(self.inputString):
            mtch=self.nucleotideSugarMatcher(self.inputString,presence=False)
            token=nsToken(mtch.group())
            return(mtch,token)

        #Modifications:
        elif self.modMatcher(self.inputString):
            mtch=self.modMatcher(self.inputString,presence=False)
            token=modToken(mtch.group())
            return(mtch,token)

        #Compartments:
        elif self.compartmentMatcher(self.inputString):
            mtch=self.compartmentMatcher(self.inputString,presence=False)
            token=compartmentToken(mtch.group())
            return(mtch,token)

        #Aglycon Matcher:
        elif self.aglyconMatcher(self.inputString):
            mtch=self.aglyconMatcher(self.inputString,presence=False)
            token=aglycoToken(mtch.group())
            return(mtch,token)
        #Wild Card Matcher:
        elif wildCardMatcher(self.inputString):
            mtch=wildCardMatcher(self.inputString,presence=False)
            token=wildCardToken(mtch.group())
            return(mtch,token)
        #Transport Arrows:
        elif self.transportMatcher(self.inputString):
            mtch=self.transportMatcher(self.inputString,presence=False)
            token=transportToken(mtch.group())
            return(mtch,token)
        #Protein Constraints:
        elif self.proteinConstraintMatcher(self.inputString):
            mtch=self.proteinConstraintMatcher(self.inputString,presence=False)
            token=proteinConstraintToken(mtch.group())
            return(mtch,token)
        #Nothing matched.  Increment within string to find
        # boundaries of unknown object, then return as unknownToken
        else:
            return(None)
        
    def matchFun(self):
        match,token=self.detectFun()
        return(match)

class monoToken(entityToken):

    def __init__(self,mono,linkage=None,modifications=None,
        compartment=None,reaction=None):
        '''
        Models Monosaccharides with no linkage information.

        Monosaccharide properties:
        - Linkage: Describes how the monosaccharide is attached
          to its adjacent moiety on a glycan structure.  If the 
          monosaccharide is not attached to a glycan, this value
          is set to None.
        - Modifications: Describes additions to the monosaccharide
          and which carbon they are attached.
        - Compartments: Describes where the monosaccharide exists.
          Typically used in transport reaction rules
        - Reactions: Reaction tokens can be associated with 
          monosaccharides if some process modifies it
        - Plurality Containers: Describes the optional presence  
          of modifications on monosaccharides.


        Can take arguments for compartment info
        '''
        self.mono=mono
        self.__name__='mono_token'
        self.compartment=compartment
        self.reactionToken=reaction
        self.modTokens=modifications
        self.linkage=linkage

    def __repr__(self):
        res='A %s monosaccharide' %(self.mono)

        if self.linkage is not None:
            res+='\n\twith a %s linkage' %(self.linkage)
        if self.modTokens is not None:
            modString='\n\t\t'.join([x.__repr__() for x in self.modTokens])
            res+='\n\twith the following modifications:\n\t\t%s' %(modString)
        if self.compartment is not None:
            res=res+'\n\tin the %s compartment' %(self.compartment.inputString)
        if self.reactionToken is not None:
            reactString='\n\t\t'+self.reactionToken.__repr__()
            res=res + '\n\thaving the following reaction: %s' %(reactString)
        return(res)

class modToken(entityToken):

    def __init__(self,string):
        self.inputString=string

    def __repr__(self):
        mod,pos=self.get_mod_pos()
        return("A modification of %s on the %s position" %(mod,pos))
    
    def get_mod_pos(self):
        pos,mod=re.search('(\d)(\D)',self.inputString).groups()
        return(mod,pos)

class nsToken(entityToken):

    def __init__(self,string,compartment=None):
        '''
        Models nucleotide sugars.
        (UDP/GDP/CMP)-(sugar)
        '''
        self.__name__='NucleotideSugar_Token'
        self.inputString=string

    def __repr__(self):
        try:
            mono,link=self.get_nt_mono()
        except:
            return('PARSE ERROR')
        return('%s nucleotide attached to %s' %(mono,link))

    def get_nt_mono(self):
        sch=re.search('(.+)\-(.+)',self.inputString)
        nt,mono=sch.groups()
        return(nt,mono)

class compartmentToken(entityToken):

    def __init__(self,inputString):
        '''
        Models instances of compartments.
        Usually are appended to the end of
        monosaccharides when describing transport
        rules.
        '''
        self.inputString=inputString

    def __repr__(self):
        return('The %s compartment' %(self.inputString))

class aglycoToken(entityToken):

    def __init__(self,inputString):
        '''
        Models aglycon instances.
        '''
        self.inputString=inputString

    def __repr__(self):
        return('The %s aglycon' %(self.inputString))

class wildCardToken(entityToken):

    def __init__(self,string):
        self.__name__='wildcardToken'
        self.inputString=string
    
    def __repr__(self):
        if '[' in self.inputString and ']' in self.inputString:
            return('One or more monosaccharide on a distinct branch')
        else:
            return('One or more monosaccharides')

class transportToken(entityToken):

    def __init__(self,string):
        self.__name__='transportToken'
        self.inputString=string

    def __repr__(self):
        return('Transports to')

class proteinConstraintToken(entityToken):

    def __init__(self,string):
        self.__name__='proteinConstraintToken'
        self.inputString=string

    def __repr__(self):
        return('A %s protein constraint' %(self.inputString))

class unknownToken(entityToken):

    def __init__(self,string):
        self.__name__='unknownToken'
        self.inputString=string

    def __repr__(self):
        return('An unknown entity: %s'%(self.inputString))

##############################
# Multi-entity Token Container
##############################

class multiToken:

    def __init__(self,string):
        self.__name__='multiToken'
        self.inputString=string
        self.entity_tokens=self.ligandTokens()

    def __repr__(self):
        try:
            s=self.get_option_list()
        except:
            return("PARSE ERROR")
        if len(s)==1:
            return('Position where %s can be optionally present' %(s[0]))
        else:
            return('Position where following are possible: ' + ' OR '.join(s))
    
    def detectFun(self):
        if multiMatcher(self.inputString):
            return(True)
        else:
            return(False)

    def matchFun(self):
        if self.detectFun():
            return(multiMatcher(self.inputString,presence=False))
        else:
            return(None)

    def get_option_list(self):
        s=re.search('\<(.+)\>',self.inputString).group()
        s=re.sub('(<|>)','',s)
        return(s.split(','))

    def ligandTokens(self):
        ''' Ligand token method for reaction token.
            Returns a list because this is a 
            multi-token container.
        '''
        ligands=self.get_option_list()
        ligand_tokens=[entityToken(l) for l in ligands]
        return(ligand_tokens)

##################
# Separator Tokens
##################

class logicalToken:

    def __init__(self,string):
        self.inputString=string
        self.logicalToken=self.get_type()

    def __repr__(self):
        return(self.get_type().__repr__())

    def detectFun(self):
        if orMatcher(self.inputString) or andMatcher(self.inputString):
            return(True)
        else:
            return(False)
    
    def matchFun(self):
        if self.detectFun():
            if orMatcher(self.inputString):
                return(orMatcher(self.inputString,presence=False))
            elif andMatcher(self.inputString):
                return(andMatcher(self.inputString,presence=False))
        else:
            return(None)
    
    def get_type(self):
        if self.detectFun():
            mtch=self.matchFun()
            if orMatcher(self.inputString):
                return(or_separator(mtch.group()))
            elif andMatcher(self.inputString):
                return(and_separator(mtch.group()))
        else:
            return(None)

class or_separator(logicalToken):

    def __init__(self,string):
        self.inputString=string

    def __repr__(self):
        return('Logical OR')

class and_separator(logicalToken):

    def __init__(self,string):
        self.inputString=string
    
    def __repr__(self):
        return('Logical AND')
