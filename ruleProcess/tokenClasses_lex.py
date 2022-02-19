import re
from tokenMatchers import *
from LexerClass import *
from itertools import chain,product as prod

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
        self.__name__='reactionToken'
        self.inputString=inputString
        #Get location of reaction match
        self.match=self.matchFun()
        if self.match is not None:
            #Get reaction type:
            self.token=self.get_reaction_type()
        else:
            self.token=None

    def __repr__(self):
        '''
        Representation takes on the value of whatever 
        reaction token was instantiated in this object.
        '''
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
            if (additionMatcher(reactionText) or termAdditionMatcher(reactionText)) and not substitutionMatcher(reactionText):
                return(additionToken(reactionText))
            elif subtractionMatcher(reactionText) or termSubtractionMatcher(reactionText):
                return(subtractionToken(reactionText))
            elif substitutionMatcher(reactionText) and not reversibleMatcher(reactionText):
                return(substitutionToken(reactionText))
            elif reversibleMatcher(reactionText):
                return(reversibleToken(reactionText))

    def substrate(self):
        return(self.token.substrate())

    def product(self):
        return(self.token.product())


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
        self.ligand_string=self.get_ligand()
        #Call the rule lexer to identify components
        # within the reaction rule:
        self.ligand_token=lexer(self.ligand_string)

    def __repr__(self):
        if self.ligand_token is not None:
            return('An addition of %s' %(self.ligand_token.__repr__()))
        else:
            return('PARSE ERROR')

    def get_ligand(self):
        ligand=re.search('\{(.*?)\}',self.inputString).groups()[0]
        return(ligand)

    def substrate(self):
        '''
        Return the ligand's representation as a 
        substrate.
        If the addition is on a terminal branch, 
          substrate should be a "[".
        
        '''
        if self.inputString[0]=='[':
            return(['['])
        else:
            return([''])

    def product(self):
        '''
        Return the ligand string
        '''
        prd=self.ligand_token[0].product()
        #Handling possible terminal branch case:
        if self.inputString[0]=='[':
            prd=''.join(['[',prd[0]])
            return(prd)
        else:
            return(prd)

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
        self.ligand_string=self.get_ligand()
        #Call the rule lexer to identify components
        # within the reaction rule:
        self.ligand_token=lexer(self.ligand_string)

    def __repr__(self):
        if self.ligand_token is not None:
            return('A removal of %s' %(self.ligand_token.__repr__()))
        else:
            return('PARSE ERROR')

    def get_ligand(self):
        ligand=re.search('\{\!(.+?)\}',self.inputString).groups()[0]
        return(ligand)

    def substrate(self):
        '''
        Return the ligand string
        '''
        return(self.ligand_token[0].substrate())
        #return(list(chain(*[x.substrate() for x in self.ligand_token])))

    def product(self):
        '''
        Return the ligand's representation as a 
        substrate.
        '''
        return([''])


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
        self.from_ligand_string,self.to_ligand_string=self.get_from_to()
        #Call the entityToken class to identify 
        # the class of the ligand:
        self.from_ligand_token,self.to_ligand_token=(lexer(self.from_ligand_string),lexer(self.to_ligand_string))

    def __repr__(self):
        if self.from_ligand_token is not None and self.to_ligand_token is not None:
            return('A conversion from %s to %s' %(self.from_ligand_token.__repr__(),self.to_ligand_token.__repr__()))
        else:
            return('PARSE ERROR')

    def get_from_to(self):
        frm,to=re.search('\{(.+?)\-\>(.+?)\}',self.inputString).groups()
        return(frm,to)

    def substrate(self):
        '''
        Return the from ligand representation:
        '''
        return(self.from_ligand_token[0].token.product())

    def product(self):
        '''
        Return the to ligand representation:
        '''
        return(self.to_ligand_token[0].token.product())
    

class reversibleToken(reactionToken):

    def __init__(self,inputString):
        '''
        This class describes reversible reactions which
        have the form {A <-> B}.
        
        The input string to this class is passed
        by the reactionToken __init__ method.
        '''
        self.__name__='substitutionToken'
        self.inputString=inputString
        self.from_ligand_string,self.to_ligand_string=self.get_from_to()
        #Call the entityToken class to identify 
        # the class of the ligand:
        self.from_ligand_token,self.to_ligand_token=(lexer(self.from_ligand_string),lexer(self.to_ligand_string))

    def __repr__(self):
        if self.from_ligand_token is not None and self.to_ligand_token is not None:
            return('A reversible conversion from %s to %s' %(self.from_ligand_token.__repr__(),self.to_ligand_token.__repr__()))
        else:
            return('PARSE ERROR')

    def get_from_to(self):
        frm,to=re.search('\{(.+?)\<\-\>(.+?)\}',self.inputString).groups()
        return(frm,to)

    def substrate(self):
        '''
        Return the from ligand representation:
        '''
        return(self.from_ligand_token[0].token.product())

    def product(self):
        '''
        Return the to ligand representation:
        '''
        return(self.to_ligand_token[0].token.product())
    


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
        self.__name__='constraintToken'
        self.inputString=inputString
        self.match=self.matchFun()
        if self.match is not None:
            self.constr=self.get_constraint_type()
        else:
            self.constr=None
        #Equate the substrate and product methods, returns 
        # empty string
        self.product=self.substrate

    def __repr__(self):
        return(self.constr.__repr__())

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
        isQuantity=quantifierMatcher(self.inputString)
        return(isQuantityStart or isAttach or isNegation or isQuantity)

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
            elif quantifierMatcher(self.inputString):
                return(quantifierMatcher(self.inputString,presence=False))
        else:
            return(None)

    def get_constraint_type(self):
        '''
        Returns string matching object where constraint
        token was detected.
        '''
        if self.detectFun():
            if quantityStartMatcher(self.inputString):
                return(quantityRule_token(quantityStartMatcher(self.inputString,presence=False)))
            elif attachRuleMatcher(self.inputString):
                return(attachRule_token(attachRuleMatcher(self.inputString,presence=False)))
            elif negationRuleMatcher(self.inputString):
                return(negationRule_token(negationRuleMatcher(self.inputString,presence=False)))
            elif quantifierMatcher(self.inputString):
                return(quantifierToken(self.inputString))
        else:
            return(None)

    def substrate(self):
        '''
        Returns an empty string, not part of actual string.
        '''
        return('')

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
        self.mtch=quantifierMatcher(self.inputString,presence=False)
    
    def __repr__(self):
        qt,val=self.mtch.groups()
        return("Preceeding pattern matches %s %s times" %(qt,val))

    def get_quantifier_quantity(self):
        '''
        Returns the quantifier as a string
        and the quantity as an integer.
        '''
        qt,val=self.mtch.groups()
        return((str(qt),int(val)))
        
    
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

    def logical_fun(self,mtchs):
        '''
        This function evaluates the length of a set 
        of match objects for a particular glycan constraint.
        Employed in constraint generation functions:
        '''
        qt,val=self.get_quantifier_quantity()

        l_fun=lambda pat,s:re.findall(pat,s)
        if qt=="=":
            return len(mtchs)==qt
        elif qt==">=":
            return len(mtchs)>=qt
        elif qt=="<=":
            return len(mtchs)<=qt
        elif qt==">":
            return len(mtchs)>qt
        elif qt=="<":
            return len(mtchs)<qt

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
        self.__name__='entityToken'
        #Object String input:
        self.inputString=inputString
        ### Get token Type: ###
        if self.detectFun():
            self.match,self.token=self.detectMain()
        else:
            #Try detecting unknown entity:
            self.match=None
            self.token=None

    def __repr__(self):
        try:
            _,token=self.detectMain()
            return(token.__repr__())
        except:
            return('PARSE ERROR')

    def mono_parse(self):
        '''
        Wrapper for the monoMatcher method.

        Returns a dictionary of terms if successfully matched.
        Otherwise returns None
        '''
        if monoMatcher(self.inputString):
            #Get the match object into a group:
            monoMatch=monoMatcher(self.inputString,presence=False)
            if re.search('^\[',monoMatch.group()) is not None:
                if re.search('\]$',monoMatch.group())  is not None:
                    leftBracket=True;rightBracket=True
                else:
                    leftBracket=True;rightBracket=False
            elif re.search('\]$',monoMatch.group()):
                leftBracket=False;rightBracket=True
            else:
                leftBracket=False;rightBracket=False
            matchDict={
                'mono':monoMatch.groups()[0],
                'compartment':monoMatch.groups()[1],
                'modPos':monoMatch.groups()[3],
                'innerReaction':monoMatch.groups()[4],
                'modType':monoMatch.groups()[5],
                'monoLink':monoMatch.groups()[6],
                'branching':{'leftBracket':leftBracket,'rightBracket':rightBracket}
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
        modTokens=[]
        if len(mono_components['modPos'])>0:
            #Create a list with the modification positions
            modPos_lst=re.findall('(\d|\<.+?\>)',mono_components['modPos'])
            if mono_components['modType'] is None:
                #Infer the mod if it is within
                # the modpos:
                modInfer=[re.search('(\d)(\D)',x) for x in modPos_lst if re.search('(\d)(\D)',x) is not None]
                if len(modInfer)==0:
                    raise Exception("Malformed Monosaccharide: there are modification positions within the monosaccharide with no valid modifications found.  Check the input rule string")
                else:
                    mtch=modInfer[0]
                    mono_components['modType']=mtch.groups()[1]
            for mp in modPos_lst:
                if '<' in mp:
                    #If no modification is present within
                    # multi-entity container, associate with
                    # "modType"
                    ents=re.search('\<(.+?)\>',mp).groups()[0].split(',')
                    newEnts=[]
                    for e in ents:
                        if len(e)==0:
                            continue
                        elif mono_components['modType'] not in e:
                            newEnts.append(e+mono_components['modType']) 
                        else:
                            newEnts.append(e)
                    mp=''.join(['<',','.join(newEnts),'>'])
                    modTokens.append(multiToken(mp))
                else:
                    modTokens.append(modToken(mp+mono_components['modType']))

        #The GalN and GlcN exception to modification rules:
        elif mono_components['mono'] in ['GalN','GlcN'] and mono_components['modType'] is not None and mono_components['modPos'] is None:
            modTokens.append(modToken(''.join(['N',mono_components['modType'][0]])))
        else:
            modTokens=None
        #Modification Reactions:
        if mono_components['innerReaction'] is not None:
            #Test different reaction types
            # Code below is meant to detect ambiguous modification additions.
            # Addition/Subtraction entities should only have one 
            # entity within them.  Substitution/Reversible tokens should
            # only have one entity in their "from" and "to" token attributes.
            rct_token=reactionToken(mono_components['innerReaction'])
            if rct_token.token.__name__ in ['additionToken','subtractionToken']:
                if rct_token.token.ligand_token[0].__name__=='unknownToken':
                    #Get Reaction Text:
                    rct_string=rct_token.token.ligand_token[0].product()[0]
                    #If the reaction text contains a comma and a number,
                    # infer the reaction is adding another of the current
                    # modification:
                    #Otherwise keep token the same and handle error later:
                    addMoreMono=re.search('((?<=\,)(\d)|(\d)(?=\,))',rct_string)
                    if addMoreMono is not None:
                        #rct_token.token.ligand_token[0]=entityToken(addMoreMono.groups()[0]+mono_components['modType'])
                        rct_token=reactionToken(''.join(['{',addMoreMono.group(),mono_components['modType'],'}']))
            elif rct_token.token.__name__ in ['substitutionToken','reversibleToken']:
                if rct_token.token.from_ligand_token[0].token.__name__=='unknownToken' or rct_token.token.to_ligand_token[0].token.__name__=='unknownToken':
                    #Get Reaction Text from "from" and "to" ligands:
                    from_rct_string=rctTok.token.from_ligand_token.product()
                    to_rct_string=rctTok.token.to_ligand_token.product()
                    #If the reaction text contains a comma and a number,
                    # infer the reaction is adding another of the current
                    # modification:
                    #Otherwise keep token the same and handle error later:
                    from_addMoreMono,to_addMoreMono=re.search('((?<=\,)(\d)|(\d)(?=\,))',from_rct_string),re.search('((?<=\,)(\d)|(\d)(?=\,))',to_rct_string)
                    if from_addMoreMono is not None:
                        rct_token.token.from_ligand_token[0]=entityToken(from_addMoreMono.groups[0]+mono_components['modType'])
                    if to_addMoreMono is not None:
                        rct_token.token.to_ligand_token[0]=entityToken(to_addMoreMono.groups[0]+mono_components['modType'])
        else:
            rct_token=None
        #Return All Attributes:
        return(mono_components['mono'],mono_components['monoLink'],mono_components['branching'],modTokens,rct_token,compartment_token)

    def detectMain(self):
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
            monoMatch=monoMatcher(self.inputString,presence=False)
            (mono,linkage,branching,modTokens,rct_token,compartment_token)=self.detectMono()
            token=monoToken(mono,linkage,branching,modTokens,compartment_token,rct_token)
            return(monoMatch,token)

        #Nucleotide Sugar:
        elif nucleotideSugarMatcher(self.inputString):
            mtch=nucleotideSugarMatcher(self.inputString,presence=False)
            token=nsToken(mtch.group())
            return(mtch,token)

        #Modifications:
        elif modMatcher(self.inputString):
            mtch=modMatcher(self.inputString,presence=False)
            token=modToken(mtch.group())
            return(mtch,token)

        #Compartments:
        elif compartmentMatcher(self.inputString):
            mtch=compartmentMatcher(self.inputString,presence=False)
            token=compartmentToken(mtch.group())
            return(mtch,token)

        #Aglycon Matcher:
        elif aglyconMatcher(self.inputString):
            mtch=aglyconMatcher(self.inputString,presence=False)
            token=aglycoToken(mtch.group())
            return(mtch,token)
        #Wild Card Matcher:
        elif wildCardMatcher(self.inputString):
            mtch=wildCardMatcher(self.inputString,presence=False)
            token=wildCardToken(mtch.group())
            return(mtch,token)
        #Transport Arrows:
        elif transportMatcher(self.inputString):
            mtch=transportMatcher(self.inputString,presence=False)
            token=transportToken(mtch.group())
            return(mtch,token)
        #Protein Constraints:
        elif proteinConstraintMatcher(self.inputString):
            mtch=proteinConstraintMatcher(self.inputString,presence=False)
            token=proteinConstraintToken(mtch.group())
            return(mtch,token)
        #Substrate:
        elif substrateMatcher(self.inputString):
            mtch=substrateMatcher(self.inputString,presence=False)
            token=substrateToken(mtch.group())
            return(mtch,token)
        #Nothing matched, return None:  
        else:
            return(None)
        
    def detectFun(self):
        try:
            match,token=self.detectMain()
            return(True)
        except:
            return(False)

    def matchFun(self):
        if self.detectFun():
            match,token=self.detectMain()
            return(match)
        else:
            return(None)

    def substrate(self):
        '''
        Returns the entity's substrate representation:
        '''
        return(self.token.substrate())

    def product(self):
        '''
        Returns the entity's product representation:
        '''
        return(self.token.product())

class monoToken(entityToken):

    def __init__(self,mono,linkage=None,branching=None,modifications=None,
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

        This object returns substrate and product representations of this
        monosaccharide using the "substrate" and "product" methods.
        '''
        self.mono=mono
        self.__name__='mono_token'
        self.compartment=compartment
        self.reactionToken=reaction
        self.modTokens=modifications
        self.linkage=linkage
        self.branching=branching
        #Flag for case when a single monosaccharide
        # continues its branching.
        #self.isextendbranch=False

    def __repr__(self):

        #Halde branching representation:
        if self.branching['leftBracket']:
            if self.branching['rightBracket']:
                res='A branched %s monosaccharide' %(self.mono)
            else:
                res='A terminal %s monosaccharide' %(self.mono)
        elif self.branching['rightBracket']:
            res='A %s monosaccharide starting on a branch' %(self.mono)
        else:
            res='A %s monosaccharide' %(self.mono)

        #Linkage Text
        if self.linkage is not None:
            res+=' with a %s linkage' %(self.linkage)
        #Modification Text:
        if self.modTokens is not None:
            modString=' '.join([str([x.__repr__()]) for x in self.modTokens])
            res+=' with the following modifications: %s' %(modString)
        #Compartment Text
        if self.compartment is not None:
            res=res+' in %s' %(self.compartment.__repr__())
        #Reaction Text:
        if self.reactionToken is not None:
            reactString=' '+self.reactionToken.__repr__()
            res=res + ' having the following reaction: %s' %(reactString)
        return(res)

    #Substrate/Product Utilities:
    def modification_perms(self):
        '''
        Creates all possible permutations of 
        a monosaccharide's modification permutations
        Accounts for the possibility of uncertainty markers
        around modification info: <2S> meaning 2S or no 2S
        '''
        #Construct a modification list:
        modTokenStrings=[]
        mod=""
        if self.modTokens is not None:
            #Get modification type:
            for t in self.modTokens:
                if t.__name__=='multiToken':
                    tokStrings=t.product()
                else:
                    tokStrings=[t.product()]
                modTokenStrings.append(tokStrings)
        #Create all permutations of "modTokenStrings", then
        # order based on modification number:
        modTokenCombs=[sorted(x) for x in list(prod(*modTokenStrings))]
        modTokenCombs=[[y for y in x if y!=''] for x in modTokenCombs]
        modTokenCombs=[[re.sub('\d','',y) if i<(len(x)-1) else y for i,y in enumerate(x)] for x in modTokenCombs]
        return(modTokenCombs)

    def react_strings(self):
        '''
        Returns strings for substrate and product elements
        within a monosaccharide element.
        '''
        #Return substrate/product representation of reaction token:
        if self.reactionToken is not None:
            sub,prd=self.reactionToken.token.substrate(),self.reactionToken.token.product()
            return(sub,prd)
        #elif self.isextendbranch:
        #    sub='['+self.mono;prod=self.mono
        #    return(sub,prod)
        else:
            return('','')

    #Creates all permutations of a monosaccharide's representation:
    def rp_stringWrap(fun):
        def _wrap(self):
            #Modifications:
            mod_perms=self.modification_perms()
            ### Decorator Function START: ###
            #Reaction Components, either substrate or product:
            compList=fun(self)
            if compList==[['']]:
                compList=['']
            ### Decorator Function END: ###
            mod_perms=[sorted(x+compList) for x in mod_perms]
            if mod_perms==[[['']]]:
                mod_perms=[['']]
            #Cleanup empties:
            mod_perms=[[y for y in x if y!=''] for x in mod_perms]
            #Generate Linkage possibilities:
            modStrings=list(chain(*[[','.join(x)] for x in mod_perms]))
            modStrings=[','.join([re.sub('\D','',y) if i<(len(x.split(','))-1) else y for i,y in enumerate(x.split(','))]) for x in modStrings]
            #Other Entities:
            #Brackets:
            leftBracket='[' if self.branching['leftBracket'] else ''
            rightBracket=']' if self.branching['rightBracket'] else ''
            #Compartments:
            compartment='' if self.compartment is None else self.compartment.substrate()
            #Linkages:
            linkage='' if self.linkage is None else self.linkage
            #Generate All possible monosaccharide representations:
            #String Order:
            # Left Bracket, Monosaccharide, Compartment, Modification Permutations,Linkage Information, Right Bracket: 
            if self.isextendbranch:
                reactStrings=[''.join(x) for x in prod(*[[leftBracket],[compartment],modStrings,[linkage],[rightBracket]])]
            else:
                reactStrings=[''.join(x) for x in prod(*[[leftBracket],[self.mono],[compartment],modStrings,[linkage],[rightBracket]])]
            return(reactStrings)
        return(_wrap)

    #Returns substrate components:
    @rp_stringWrap
    def substrate(self):
        '''
        Builds monosaccharide "substrate" representation
        Fills in substrate-specific components
        '''
        #Gather monosaccharide modification reaction info, and
        # integrate into the modification permutations:
        (substrs,_)=self.react_strings()
        return([substrs])

    #Returns product components:
    @rp_stringWrap
    def product(self):
        '''
        Builds monosaccharide "product" representation
        Fills in product-specific components
        '''
        #Gather monosaccharide modification reaction info, and
        # integrate into the modification permutations:
        (_,prods)=self.react_strings()
        return([prods])


class modToken(entityToken):

    def __init__(self,string):
        self.__name__='modToken'
        self.inputString=string
        self.product=self.substrate

    def __repr__(self):
        mod,pos=self.get_mod_pos()
        return("A modification of %s on the %s position" %(mod,pos))
    
    def get_mod_pos(self):
        pos,mod=re.search('(\d|N)(\D)',self.inputString).groups()
        return(mod,pos)
    
    def substrate(self):
        return(self.inputString)

class nsToken(entityToken):

    def __init__(self,string,compartment=None):
        '''
        Models nucleotide sugars.
        (UDP/GDP/CMP)-(sugar)
        '''
        self.__name__='NucleotideSugar_Token'
        self.inputString=string
        self.product=self.substrate

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
    
    def substrate(self):
        return([self.inputString])

class compartmentToken(entityToken):

    def __init__(self,inputString):
        '''
        Models instances of compartments.
        Usually are appended to the end of
        monosaccharides when describing transport
        rules.
        '''
        self.__name__='compartmentToken'
        self.inputString=inputString
        self.product=self.substrate

    def __repr__(self):
        return('The %s compartment' %(self.inputString))

    def substrate(self):
        return(self.inputString)

class aglycoToken(entityToken):

    def __init__(self,inputString):
        '''
        Models aglycon instances.
        '''
        self.__name__='aglycoToken'
        self.inputString=inputString
        self.product=self.substrate

    def __repr__(self):
        return('The %s aglycon' %(self.inputString))

    #def substrate(self):
    #    return([self.inputString])

    def substrate(self):
        return(["$"])


class wildCardToken(entityToken):

    def __init__(self,string):
        self.__name__='wildcardToken'
        self.inputString=string
        self.product=self.substrate
    
    def __repr__(self):
        if '[' in self.inputString and ']' in self.inputString:
            return('One or more monosaccharide on a distinct branch')
        else:
            return('One or more monosaccharides')

    def substrate(self):
        return([self.inputString])

class transportToken(entityToken):

    def __init__(self,string):
        self.__name__='transportToken'
        self.inputString=string
        self.product=self.substrate

    def __repr__(self):
        return('Transports to')

    def substrate(self):
        return([self.inputString])

class proteinConstraintToken(entityToken):

    def __init__(self,string):
        self.__name__='proteinConstraintToken'
        self.inputString=string
        self.product=self.substrate

    def __repr__(self):
        return('A %s protein constraint' %(self.inputString))

    def substrate(self):
        return([self.inputString])

class substrateToken(entityToken):

    def __init__(self,string):
        self.__name__='substrateToken'
        self.inputString=string
        self.product=self.substrate

    def __repr__(self):
        if self.inputString=='R':
            return('Some arbitrary substrate')
        else:
            return('A %s substrate' %(self.inputString))
    
    def substrate(self):
        return([self.inputString])

class unknownToken(entityToken):

    def __init__(self,string):
        self.__name__='unknownToken'
        self.inputString=string
        self.product=self.substrate

    def __repr__(self):
        return('An unknown entity: %s'%(self.inputString))

    def substrate(self):
        return([self.inputString])

##############################
# Multi-entity Token Container
##############################

class multiToken:

    def __init__(self,string):
        self.__name__='multiToken'
        self.inputString=string
        self.entity_strings=self.get_option_list()
        self.tokens=list(chain(*[lexer(e) for e in self.entity_strings]))
        self.product=self.substrate

    def __repr__(self):
        if len(self.tokens)>0:
            if len(self.tokens)==1:
                return('Position where %s can be optionally present' %(self.tokens[0]))
            else:
                tokenList=','.join([x.__repr__() for x in self.tokens])
                return(' '.join(['Position where following are possible: ',' OR '.join([str([x.__repr__()]) for x in self.tokens])]))
        else:
            return('PARSE ERROR')
    
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
        s=re.search('\<(.+?)\>',self.inputString).group()
        s=re.sub('(<|>)','',s)
        return(s.split(','))
    
    def substrate(self):
        '''
        Returns substrate/product representations for 
        multiToken-contained entities:
        '''
        if len(self.tokens)==1:
            if type(self.tokens[0].product())==str:
                return(['',self.tokens[0].product()])
            else:
                return(['',self.tokens[0].product()[0]])
        else:
            return([x.product() for x in self.tokens])


##################
# Separator Tokens
##################

class logicalToken:

    def __init__(self,string):
        self.__name__='logicalToken'
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
        self.__name__='or_separator'
        self.inputString=string

    def __repr__(self):
        return('Logical OR')

class and_separator(logicalToken):

    def __init__(self,string):
        self.__name__='and_separator'
        self.inputString=string
    
    def __repr__(self):
        return('Logical AND')

#######################
# Default Lexer:
#######################

lexer=LexerClass(lexicon=[reactionToken,constraintToken,entityToken,multiToken,logicalToken],ukToken=unknownToken)

#########################
# Unexpected Token Errors
#########################

class unexpectedTokenError(Exception):

    def __init__(self,expect_t,t):
        self.msg="Was expecting to find {expect_t} but found {t}"
