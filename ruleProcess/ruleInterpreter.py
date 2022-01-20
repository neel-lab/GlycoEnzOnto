import re
from tokenClasses_lex import *
from LexerClass import *
from itertools import product as prod,chain,tee
from dataclasses import dataclass
import functools
from functools import *
import difflib
from difflib import ndiff
from typing import Any

def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)

#######################
# Default Lexer:
#######################

lexer=LexerClass(lexicon=[reactionToken,constraintToken,entityToken,multiToken,logicalToken],ukToken=unknownToken)

##########################
# Rule String Exceptions
##########################

class reactionRuleErrorException(Exception):
    def __init__(self,basicValidation,hasNoConstraints):
        self.basicValidation=basicValidation
        self.hasNoConstraints=hasNoConstraints
        if not self.basicValidation:
            raise basicValidationError()
        elif not self.hasNoConstraints:
            raise ConflictingRuleException()
        pass

class constraintRuleErrorException(Exception):

    def __init__(self,basicValidation,hasNoReaction):
        self.basicValidation=basicValidation
        self.hasNoReaction=hasNoReaction
        if not self.basicValidationError:
            raise basicValidationError()
        elif not self.hasNoReaction:
            raise ConflictingRuleException()
        else:
            raise malformedConstraintError()

class basicValidationError(Exception):
    def __init__(self):
        self.message="Some tokens are unknown in the input"
        super().__init__(self.message)

class ConflictingRuleException(Exception):
    def __init__(self):
        self.message="Reaction and Constraint component detected in same rule string, this is not allowed"
        super().__init__(self.message)

class malformedConstraintError(Exception):
    def __init__(self):
        self.message="No constraint-specific text detected within rule string."
        super().__init__(self.message)


######################
# Rule Regex Converter
######################

class generatorSet:

    def rule2regex(self,rl):
        '''
        Converts reaction rules into regex string.  Makes
        syntax compatible with regular expressions.
        '''
        # Escape Characters:
        rl=re.sub(r'\(','\(',rl)
        rl=re.sub(r'\)','\)',rl)
        rl=re.sub(r'\[','\[',rl)
        rl=re.sub(r'\]','\]',rl)
        rl=re.sub(r'\-','\-',rl)
        #If a wild card exists in front of pattern,
        # match can happen internally on a glycan structure.
        wild=re.search(r'(?!^)\.\.\.',rl)
        frontwild=re.search(r'^\.\.\.',rl)
        #Order for wild card detection is important:
        if frontwild is not None:
            rl=re.sub('^\.\.\.','',rl)
        else:
            rl=re.sub('^','(?:^|\[)',rl)
        if wild is not None:
            rl=re.sub('\.\.\.','(.+)',rl)
        # Uncertain linkages:
        rl=re.sub(r'\-\?','-[0-9]',rl)
        #If a core linkage is detected, append a "$" to the 
        # string to look at the core:
        if re.search(r'\([ab\?][12\?]\-$',rl) is not None:
            rl=re.sub(r'(\([ab][12]\-$)','\g<1>$')
        return(rl)

##########################
# Rule Superclass
##########################

class Rule:

    def __init__(self,ruleComponents):
        '''
        Interprets parsed rule string lists by:
        - Checking if combination of token types is allowed.
        '''
        self.ruleComponents=ruleComponents
        self.ruleSets,self.logicalSeps=self.getRuleSets(self.ruleComponents)

    #######################
    #Rule Splitter Utility:
    #######################
    @staticmethod
    def getRuleSets(ruleComponents):
        '''
        Separates rule components sequences by 
        logicalToken types.  These sequences are then
        processed for validity.
        '''
        ruleSets=[]
        separators=[]
        c_set=[]
        for e in ruleComponents:
            if e.__name__=='logicalToken':
                ruleSets.append(c_set)
                c_set=[]
                separators.append(e)
                continue
            c_set.append(e)
        ruleSets.append(c_set)
        return(ruleSets,separators)

    ##########################
    # Rule Validation Wrapper:
    ##########################
    def checkWrapper(fun):
        def _check(self,ruleSets):
            '''
            Applies a function that checks if each
            rule set in the rule string is valid.
            The function should return either "True"
            or "False" for a rule set.
            The checker determines if all of these rule
            sets return something that is NOT "False".
            Responses can be "None" or "True"
            '''
            res=all([fun(x) is not False for x in ruleSets])
            return(res)
        return(_check)

    def allTrueWrap(fun):
        def _check_sub(ruleSet):
            '''
            Checks if all tokens within a ruleSet are 
            True.
            '''
            res=all(fun(ruleSet))
            return(res)
        return(_check_sub)

    def possibleTrueWrap(fun):
        def _check_sub(ruleSet):
            '''
            If a ruleSet matches none of the "fun"
            constraints, returns None.

            If only some of the constraints in "fun"
            are matched, returns "False", means something
            is wrong with the rule's formatting.

            Otherwise, returns "True" if all conditions
            in "fun" are satisfied.
            '''
            validQuantityRuleToken=True if ruleSet[0].__name__=='quantityRule_token' else False
            validQuantityToken=True if ruleSet[-1].__name__=='quantifierToken' else False
            #Statements is a list of booleans returned by
            # the "fun":
            statements=fun(ruleSet)
            if sum(statements)==0:
                #None of the tags match, assume the 
                # token type was not evoked.  Return None
                return(None)
            elif sum(statements)<len(statements):
                #Some of the statements were correct about 
                # the token but not all.  Means the token is
                # malformed
                return(False)
            else:
                #All are true, meaning this is a valid quantification rule:
                return(True)
        return(_check_sub)

    #############################
    # Universal Token Validators:
    #############################
    
    @checkWrapper
    @allTrueWrap
    def validTokens(ruleSet):
        return([(x.substrate() is not None and x.product() is not None) for x in ruleSet])

    @checkWrapper
    @allTrueWrap
    def noUnknownTokens(ruleSet):
        return([x.__name__!='unknownToken' for x in ruleSet])

    def basicValidationWrapper(self,ruleSets):
        return(self.validTokens(ruleSets) and self.noUnknownTokens(ruleSets))


##########################
# Reaction Rule Class:
##########################

class GlycanProcessorGenerator(generatorSet):
    '''
    Class which returns from/to strings for 
    processing substrates/products
    '''
    def __init__(self,fromString,toString,reactionType):
        self.fromString=fromString
        self.toString=toString
        print(self.fromString)
        print(self.toString)
        self.reactionType=reactionType

    #Call method for GlycanProcessorGenerator:
    def __call__(self,glycan):
        return(self.makeProducts(glycan))
     
    def constructToString(self,mtch):
        '''
        Calls the difflib "ndiff" routine to construct
        the new "to" string:
        '''
        toString_clean=self.toString
        if len(mtch.groups())>0:
            toString_clean=re.sub('\.\.\.',mtch.groups()[0],toString_clean)
        else:
            #Assume the dots are at the front, indicating a
            # wild card group:
            toString_clean=re.sub('^\.{3}','',toString_clean)
        toString_clean=re.sub('\$$','',toString_clean)
        resList=[]
        #Get Differential Symbol:
        if self.reactionType==additionToken:
            dsym='-'
        elif self.reactionType==subtractionToken:
            dsym='+'
        for i,elt in enumerate(ndiff(toString_clean,mtch.group())):
            df,elt=re.search('(.)\ (.)',elt).groups()
            #Handles replacing ambiguous links here:
            if (df=='-' and elt=='?'):
                continue
            elif (df in [' ',dsym]) or (df=='+' and elt.isdigit()):
                resList.append(elt)
        toString_clean=''.join(resList)
        #Remove do
        #Internal branching error:
        if re.search('^\[.+?\].+?\]$',toString_clean) is not None:
            toString_clean=re.sub('^(\[.+?)\](.+?\])$','\g<1>\g<2>',toString_clean)
        #Beginning branch error:
        elif re.search('^\]',toString_clean) is not None:
            toString_clean=re.sub('^\]','',toString_clean)
        #Fix front branch:
        if mtch.group()[0]=='[' and toString_clean[0]!='[':
            toString_clean=''.join(['[',toString_clean])
        print(toString_clean)
        return(toString_clean)

    def makeToRepString(self,fromWildGrp):
        toRepString=re.sub('(?!^)\.\.\.',fromWildGrp,self.toString)
        return(toRepString)
    
    def getGlycanMatch(self,fromString_regex,glycan):
        '''
        Performs the glycan matching procedure.
        Returns all possible matches
        '''
        return(re.finditer(fromString_regex,glycan))

    def fixGlycanBranching(self,prd):
        '''
        Fixes incorrect branching:
        '''
        newProd=re.sub('(\[.+?)\](.+?\])','\g<1>\g<2>',prd)
        return(newProd)

    def fixGlycanBranching_main(self,prd):
        '''
        Confirms correct branching and removes branching errors:
        '''
        posDta=[(x.start(),x.group()=='[') for x in re.finditer('\[|\]',prd)]
        #If the length of "posDta" is odd,
        # means something is wrong with connectivity.
        resolved=True
        delInds=[]
        for (i,(x_strt,x_isEnd)),(j,(y_strt,y_isEnd)) in pairwise(enumerate(posDta)):
            #Fully resolved or a beginning of branch:
            if (x_isEnd and not y_isEnd):
                continue
            #Nested branch:
            if x_isEnd and y_isEnd:
                #Store index for later:
                delInds.append(x_strt)
                resolved=False
            #Resolving nested branch:
            if not resolved:
                if not y_isEnd:
                    delInds.pop()


    def makeProducts(self,glycan):
        '''
        Method generates a list of strings that look
        like the "to" pattern.
        '''
        fromString_regex=super().rule2regex(self.fromString)
        print(fromString_regex)
        #self.toString=re.sub('^\.\.\.','',self.toString)
        mtchs=self.getGlycanMatch(fromString_regex,glycan)
        #Initialize a product list:
        products=[]
        for m in mtchs:
            to=self.constructToString(m)
            print(to)
            #Re-add a "[" if processing the head of a branch:
            #if re.search('^\[',m.group()) is not None:
            #    to=''.join(['[',to])
            #Constant indicies:
            front_start=0
            front_end=m.start()
            back_start=m.end()
            back_end=len(glycan)
            #Replace text where "frm" was found with the "to" text:
            products.append(''.join([glycan[front_start:front_end],to,glycan[back_start:back_end]]))
            #Fix broken branching here:
            #products=[self.fixGlycanBranching(p) for p in products]
        return(products)



class reactionRule(Rule):

    def __init__(self,ruleComponents):
        super().__init__(ruleComponents)
        #Check Conditions for each rule set:
        # Default constraints:
        basicValidation=self.basicValidationWrapper(self.ruleSets)
        # Reaction-specific constraints:
        hasNoConstraints=self.noConstraints(self.ruleSets)
        #If all are true, returns an instance of class.
        #Otherwise returns none:
        if (not basicValidation) or (not hasNoConstraints):
            raise reactionRuleErrorException(basicValidation,hasNoConstraints)
        #Store reaction type:
        self.reactionType=list(set(chain(*[[type(x.token) for x in rst if isinstance(x,reactionToken)] for rst in self.ruleSets])))[0]
        #If rule is valid, create generator terms:
        #Instantiate forward and referse inference methods:
        self.forward=self.forwardGeneratorMain()
        self.reverse=self.reverseGeneratorMain()

    ############################
    # Reaction Rule Validation:
    ############################

    @Rule.checkWrapper
    @Rule.allTrueWrap
    def noConstraints(ruleSet):
        return([x.__name__!='constraintToken' for x in ruleSet])

    #######################
    # Class Factory Method:
    #######################

    @classmethod
    def fromComponents(cls,ruleComponents):
        '''
        Main method for determining if passed reaction
        rule components are reaction rules.
        Criteria:
        1. All tokens must be error-free instantiations.
        2. All tokens must not be unknown
        3. No tokens should be Constraints.
        '''
        return(cls(ruleComponents))
        
    ######################################
    # Substrate/Product Pair List Builder:
    ######################################
    def pairListBuilder(self,ruleSet):
        '''
        Method that invokes substrate and product 
        building from tokens in a ruleSet:
        Returns pairs of tokens.
        '''
        #Create every permutation of substrate/product strings:
        substrates=[''.join(x) for x in prod(*[y.substrate() for y in ruleSet])]
        products=[''.join(x) for x in prod(*[y.product() for y in ruleSet])]
        #Return substrate/product pairs:
        return([(s,p) for s,p in zip(substrates,products)])

    #def pairListBuilder(self,ruleSet):
    #    '''
    #    Method that invokes substrate and product 
    #    building from tokens in a ruleSet:
    #    Returns pairs of lists containing substrate/product components.
    #    '''
    #    #Create every permutation of substrate/product strings:
    #    substrates=list(prod(*[y.substrate() for y in ruleSet]))
    #    products=list(prod(*[y.product() for y in ruleSet]))
    #    #Return substrate/product pairs:
    #    return([(s,p) for s,p in zip(substrates,products)])

    def pairListGenerator(self):
        return(functools.reduce(lambda x,y: x+y,[self.pairListBuilder(x) for x in self.ruleSets]))

    ###############################################
    # Reaction Search/Replace Generating Functions:
    ###############################################

    def makeGlycanProcessor(self,frm,to):
        '''
        Wrapper to dynamically create instances
        of "proc_" with "frm" "to" pairs.
        '''
        gpg=GlycanProcessorGenerator(frm,to,self.reactionType)
        return(lambda glycan:gpg(glycan))

    def glycanProcAggregator(fun):
        def _wrap(self):
            '''
            Returns a function that converts glycans 
            based on given "frm" "to" strings.
            "fun" reorders the pairList to perform
            "forward" or "reverse" inference.
            '''
            #Generate the pairList:
            pairList=self.pairListGenerator()
            #Reordering the pairList to return either
            # substrate->product or product->substrate.
            frm_list,to_list=fun(pairList)
            #Make a list of replacer functions:
            funList=[self.makeGlycanProcessor(frm,to) for frm,to in zip(frm_list,to_list)]
            #Make conversion function:
            convertMain=functools.reduce(lambda cur,pres: lambda string: cur(string)+pres(string),funList)
            return(convertMain)
        return(_wrap)

    ##############################
    # Inference Generator Wrapper:
    ##############################

    def inferWrapper(fun):
        def _wrap(self):
            pairList=self.pairListGenerator()
            inferFun=fun(pairList)
            return(inferFun)
        return(_wrap)
    
    ######################################
    # Forward/Reverse Inference Functions:
    ######################################

    @glycanProcAggregator
    def forwardGeneratorMain(pairList):
        '''
        Produces a list of all possible glycans where "substrate"
        is taken to "product"
        '''
        return([t[0] for t in pairList],[t[1] for t in pairList])

    @glycanProcAggregator
    def reverseGeneratorMain(pairList):
        '''
        Produces a list of all possible glycans where "product"
        is taken to "substrate"
        '''
        return([t[1] for t in pairList],[t[0] for t in pairList])

##########################
# Constraint Rule Class:
##########################

class ConstraintMethodGenerator(generatorSet):
    '''
    This class creates methods which evaluate 
    substrates/products for their validity.
    This class is instantiated by passing lists of
    rule tokens which are searched for the following attributes:
    1. Negation: The prefix of a constraint string prohibiting a 
       structure to exist.
    2. Numeric: This flag denotes that the structural components
       of the constraint string have some constraint on its 
       frequency.
    3. Attachment: This flag denotes there is an "attachment" 
       constraint in the rule.  Will trigger methods to recognize
       the attachment constraint.
    '''
    def __init__(self,ruleSet,reactionRule=None):
        '''
        Parses constraint rule components within a
        rule set and returns tags used for 
        the constraint constructor.
        Keeps variables tracking if the rule is a negation, attachment,
        or numeric constraint while parsing monosaccharide elements into
        "seqSet"
        '''
        self.negation=False
        self.numeric=None
        self.attachment=None
        self.reactionRule=reactionRule
        self.seqSet=[]
        self.addMono=None
        for i,t in enumerate(ruleSet):
            if t.__name__=='constraintToken':
                if t.constr.__name__=='negationRule_token':
                    self.negation=True
                elif t.constr.__name__=='quantityRule_token':
                    self.numeric=True
                    if ruleSet[-1].constr.__name__!='quantifierToken':
                        raise Exception("Quantity rule detected but no quantifier/quantity provided")
                    else:
                        self.numeric=t.constr
                elif t.constr.__name__=='attachRule_token':
                    self.attachment=True
                    if self.reactionRule is None:
                        raise Exception("Attachment constraint detected but no monosaccharide provided for attachment")
                    else:
                        #Instantiates "addMono" variable within object
                        self.addMono=self.getAddMono()
                        #Append the attachment constraint so that it can be replaced
                        # in the "createSeq" method:
                        self.seqSet.append('@')
            else:
                self.seqSet.append(t.product()[0])

        #Generate the constraint:
        self.constraint=self.constraintGen()

    #Call method for ConstraintMethodGenerator:
    def __call__(self,glycan):
        return(self.constraint(glycan))

    def getAddMono(self):
        '''
        For constraints that have addition constraints,
        find the monosaccharide which is added in a parsed 
        reaction rule.
        '''
        #Search for all addition possibilities in the 
        # reactionRule's ruleSets:
        addMono=list(set(chain(*[[x.product()[0] for x in rst if x.token.__name__=='additionToken'] for rst in self.reactionRule.ruleSets])))[0]
        return(addMono)

    def createSeq(self):
        '''
        Generates monosaccharide sequence from
        "seqSet"
        '''
        seq=''.join([x for x in self.seqSet])
        if self.attachment:
            if self.addMono is None:
                #This exception shouldn't trigger as this is 
                # processed in the __init__ method of this class.
                #Keeping just in case for debugging reasons:
                raise Exception("Attachment constraint detected but no monosaccharide provided for attachment")
            else:
                #Replace the @ symbol with the added monosaccharide
                seq=re.sub('\@',self.addMono,seq)
        return(seq)
    
    def constraintGen(self):
        '''
        Generates a constraint method which evaluates if a 
        generated product is valid.
        '''
        #Base search function:
        stringSearch=self.createSeq()
        stringSearch_regex=super().rule2regex(stringSearch)
        funOut_search=lambda glycan: re.finditer(stringSearch_regex,glycan)
        #Negation:
        if self.negation:
            funOut=lambda glycan: len(list(funOut_search(glycan)))==0
        else:
            funOut=lambda glycan: len(list(funOut_search(glycan)))>0
        #Numerical Constraint:
        if self.numeric is not None:
            #The quantifier object contains a method
            # Which automatically evaluates if the quantity
            # of patterns matched fulfills the quantifier.
            funOut=lambda glycan: self.numeric.logical_fun(funOut_search(glycan))
        return(funOut)




class constraintRule(Rule):

    def __init__(self,ruleComponents,reactionRule=None):
        #Creates ruleSets and logical seps:
        super().__init__(ruleComponents)
        #Check Conditions for each rule set:
        # Default constraints:
        basicValidation=self.basicValidationWrapper(self.ruleSets)
        # Reaction-specific constraints:
        hasNoReactions=self.noReactions(self.ruleSets)
        # Valid constraint formatting:
        hasValidNegation=self.validNegation(self.ruleSets);hasValidQuantityConstraint=self.validNumeric(self.ruleSets)
        #There is an issue if the passed string has NO valid 
        # forms of constraints:
        if all([not x for x in [hasValidNegation,hasValidQuantityConstraint]]):
            raise constraintRuleErrorException(basicValidation,hasNoReaction)
        #If all checks pass, create the constraint method:
        self.reactionRule=reactionRule
        self.constraint=self.ConstraintGenerator_Aggregator()

    def __call__(self,glycString):
        return(self.constraint(glycString))

    #############################
    # Constraint Rule Validation:
    #############################

    @Rule.checkWrapper
    @Rule.allTrueWrap
    def noReactions(ruleSet):
        return([x.__name__!='reactionToken' for x in ruleSet])

    @Rule.checkWrapper
    @Rule.possibleTrueWrap
    def validNumeric(ruleSet):
        '''
        Numeric constraints must have the quantity rule prefix "n"
        as well as the quantity constraint as the suffix.
        '''
        validQuantityRuleToken=True if ruleSet[0].__name__=='quantityRule_token' else False
        validQuantityToken=True if ruleSet[-1].__name__=='quantifierToken' else False
        return([validQuantityRuleToken,validQuantityToken])

    @Rule.checkWrapper
    @Rule.possibleTrueWrap
    def validNegation(ruleSet):
        '''
        Checks for negation rule presence and validity  
        Negation rules must come at the FRONT of any constraint.
        Thus, only the first reaction rule token can be the negation
        rule token.
        '''
        validNegation=True if ruleSet[0]=='negationRule_token' else False
        return([validNegation])

    #######################
    # Class Factory Method:
    #######################

    @classmethod
    def fromComponents(cls,ruleComponents):
        '''
        Main method for determining if passed reaction
        rule components are constraint rules.
        Criteria:
        1. All tokens must be error-free instantiations.
        2. All tokens must not be unknown
        3. No tokens should be Reaction Rules.
        4. Constraints should be formated correctly:
          4a. Correctly formatted negation if present.
          4b. Correctly formatted numerical constraints if present.
        '''
        #Check Conditions for each rule set:
        # Default constraints:
        basicValidation=self.basicValidationWrapper(self.ruleSets)
        # Reaction-specific constraints:
        hasNoReactions=self.noReactions(self.ruleSets)
        # Valid constraint formatting:
        hasValidNegation=self.validNegation();hasValidQuantityConstraint=self.validNumeric()
        #If all are true, returns an instance of class.
        #Otherwise returns none:
        if basicValidation and all([hasNoReactions,hasValidNegation,hasValidQuantityConstraint]):
            return(cls(ruleComponents))
        else:
            return(None)
    
    def ConstraintGenerator_Aggregator(self):
        '''
        Takes a set of ConstraintGenerator classes and 
        concatenates these into one function.
        This allows the substrates/products to be evaluated 
        with only one function
        '''
        #Make Constraint objects:
        ConstraintFuns=[ConstraintMethodGenerator(r,self.reactionRule) for r in self.ruleSets]
        #ConstraintFuns=[ConstraintMethodGenerator.fromComponents(r).constraintGen() for r in self.ruleSets]
        #Constraint_Classes=[ConstraintMethodGenerator.fromComponents(r) for r in self.ruleSets]
        #Merge into one function:
        def mergeWrapper(acc,_zipInfo):
            sep,cur=_zipInfo
            if sep=='&':
                res=lambda within: acc(within) and cur(within)
            elif sep=='|':
                res=lambda within: acc(within) or cur(within)
            else:
                res=lambda within: acc(within)
            return(res)
        #Merge the functions together
        # "initialFunction" is the first constraint in the "Constraint_Classes" list.
        initialFunction=lambda within:ConstraintFuns[0](within)
        pairedIterator=zip(ConstraintFuns[1:],self.logicalSeps)
        constraintMain=functools.reduce(lambda acc,x:mergeWrapper(acc,x),pairedIterator,initialFunction)
        return(constraintMain)
