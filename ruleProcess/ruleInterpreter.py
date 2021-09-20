import re
from tokenClasses_lex import *
from LexerClass import *
from itertools import product as prod
from dataclasses import dataclass
import functools

#######################
# Default Lexer:
#######################

lexer=LexerClass(lexicon=[reactionToken,constraintToken,entityToken,multiToken,logicalToken],ukToken=unknownToken)

##########################
# Rule String Exceptions
##########################

class ConflictingRuleException(Exception):

    def __init__(self):
        self.message="Reaction and Constraint component detected in same rule string, this is not allowed"
        super().__init__(self.message)

class ConflictingRuleException(Exception):

    def __init__(self):
        self.message="Reaction and Constraint component detected in same rule string, this is not allowed"
        super().__init__(self.message)

######################
# Rule Regex Converter
######################

class regexMethodSet:

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
        if frontwild is not None:
            rl=re.sub('^\.\.\.','',rl)
        elif frontwild is None:
            rl=re.sub('^','(?:^|\[)',rl)
        elif frontwild is None and wild is not None:
            rl=re.sub('\.\.\.','(.+?)',rl)
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
        self.ruleSets,self.logicalSeps=self.getRuleSets()

    #######################
    #Rule Splitter Utility:
    #######################
    def getRuleSets(self,ruleComponents):
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
        def _check(self,ruleSet):
            '''
            Checks if all tokens within a ruleSet are 
            True.
            '''
            res=all(fun(ruleSet))
            return(res)
        return(_check)

    def possibleTrueWrap(fun):
        def _check(self,ruleSet):
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
        return(_check)



    #############################
    # Universal Token Validators:
    #############################
    
    @checkWrapper
    @allTrueWrap
    def validTokens(self,ruleSet):
        return([(x.substrate() is not None and x.product() is not None) for x in ruleSet])

    @checkWrapper
    @allTrueWrap
    def noUnknownTokens(self,ruleSet):
        return([x.__name__!='unknownToken' for x in ruleSet])

    def basicValidationWrapper(self,ruleSets):
        return(validTokens(ruleSets) and noUnknownTokens(ruleSets))


##########################
# Reaction Rule Class:
##########################

@dataclass
class GlycanProcessorGenerator(regexMethodSet):
    '''
    Class which returns from/to strings for 
    processing substrates/products
    '''
    fromString: str
    toString: str

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
        if frontwild is not None:
            rl=re.sub('^\.\.\.','',rl)
        elif frontwild is None:
            rl=re.sub('^','(?:^|\[)',rl)
        elif frontwild is None and wild is not None:
            rl=re.sub('\.\.\.','(.+?)',rl)
        # Uncertain linkages:
        rl=re.sub(r'\-\?','-[0-9]',rl)
        #If a core linkage is detected, append a "$" to the 
        # string to look at the core:
        if re.search(r'\([ab\?][12\?]\-$',rl) is not None:
            rl=re.sub(r'(\([ab][12]\-$)','\g<1>$')
        return(rl)

    def makeToRepString(self,fromWildGrp):
        toRepString=re.sub('(?!^)\.\.\.',fromWildGrp,self.toString)
        return(toRepString)
    
    def getGlycanMatch(self,fromString_regex,glycan):
        '''
        Performs the glycan matching procedure.
        Returns all possible matches
        '''
        return(re.finditer(fromString_regex,glycan))
    
    def makeProducts(self,glycan):
        '''
        Method generates a list of strings that look
        like the "to" pattern.
        '''
        fromString_regex=super().rule2regex(self.fromString)
        self.toString=re.sub('^\.\.\.','',self.toString)
        mtchs=self.getGlycanMatch(fromString_regex,glycan)
        #Initialize a product list:
        products=[]
        for m in mtchs:
            #Matches with groups are assumed to have 
            # wild card regions within the match string.
            #Take the matched wild card text and replace it
            # in the to string using the "makeToRepString"
            # method:
            if len(m.groups())!=0:
                #Assume 1 group:
                fromWildGroup=m.groups()[0]
                to=self.makeToRepString(fromWildGroup)
            #Otherwise, there are no wild card token in the 
            # substrate/product strings, and just return the
            # string contents.
            else:
                to=self.toString
            #Add a "[" if the head of a branch:
            if re.search('^\[',m.group()) is not None:
                to=''.join(['[',to])
            #Constant indicies:
            front_start=0
            front_end=m.start()
            back_start=m.end()
            back_end=len(glycan)
            #Replace text where "frm" was found with the "to" text:
            products.append(''.join([glycan[front_start:front_end],to,glycan[back_start:back_end]]))
        return(products)

    def __call__(self,glycan):
        return(self.makeProducts(glycan))


class reactionRule(Rule):

    def __init__(self):
        super().__init__(self,ruleComponents)
        #Instantiate forward and referse inference methods:
        self.forward=self.forwardGeneratorMain()
        self.reverse=self.reverseGeneratorMain()


    ############################
    # Reaction Rule Validation:
    ############################
    @Rule.checkWrapper
    @Rule.allTrueWrap
    def noConstraints(self,ruleSet):
        return([x.__name__!='constraintToken' for x in ruleSet])

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
        #Split rule:
        ruleSets,logicalSeps=self.getRuleSets(ruleComponents)
        #Check Conditions for each rule set:
        # Default constraints:
        basicValidation=self.basicValidationWrapper(self.ruleSets)
        # Reaction-specific constraints:
        hasNoConstraints=self.noConstraints(ruleSets)
        #If all are true, returns an instance of class.
        #Otherwise returns none:
        if basicValidation and hasNoConstraints:
            return(cls(ruleComponents))
        else:
            return(None)
        
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

    def pairListGenerator(self):
        return(functools.reduce(sum,[pairListBuilder(x) for x in self.ruleSets]))

    ###############################################
    # Reaction Search/Replace Generating Functions:
    ###############################################

    def makeGlycanProcessor(self,frm,to):
        '''
        Wrapper to dynamically create instances
        of "proc_" with "frm" "to" pairs.
        '''
        gpg=GlycanProcessGenerator(frm,to)
        return lambda glycan:gpg(glycan)
        #return lambda glycan: proc_(glycan,frm,to)

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
            funList=[makeGlycanProcessor(frm,to) for frm,to in zip(frm_list,to_list)]
            #Make conversion function:
            convertMain=reduct(lambda cur,pres: lambda string: cur(string)+pres(string),funList)
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
    def forwardGeneratorMain(self,pairList):
        '''
        Produces a list of all possible glycans where "substrate"
        is taken to "product"
        '''
        return([t[0] for t in pairList],[t[1] for t in pairList])

    @glycanProcAggregator
    def reverseGeneratorMain(self,pairList):
        '''
        Produces a list of all possible glycans where "product"
        is taken to "substrate"
        '''
        return([t[1] for t in pairList],[t[0] for t in pairList])

##########################
# Constraint Rule Class:
##########################

@dataclass
class ConstraintMethodGenerator(regexMethodSet):
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
    negation: bool
    numeric: quantifierToken = None
    attachment: bool
    seqSet: list
    addMono: str = None

    def createSeq(self):
        '''
        Generates monosaccharide sequence from
        "seqSet"
        '''
        seq=''.join([x.product() for x in self.seqSet])
        if self.attachment:
            if addMono is None:
                raise Exception("Attachment constraint detected but no monosaccharide provided for attachment")
            else:
                seq=re.sub('\@',self.addMono)
        return(seq)
    
    # "searchDecorator" takes arguments "isNegation" and "numericObj"
    # which are passthrough variables for "negation" and "numeric"
    # for this object.  Additional function modifiers can be made
    # within this block to expand constraint processing:
    def searchDecorator(isNegation=negation,numericObj=numeric):
        def _searchDecWrap(self,fun):
            '''
            Decorates the main search function "fun"
            with additional logical/numeric constraints.
            '''
            #Base search function:
            funOut=fun()
            #Negation
            if isNegation:
                funOut=lambda glycan: funOut(glycan) is None
            else:
                funOut=lambda glycan: funOut(glycan) is not None
            #Numerical Constraint
            if numericObj is not None:
                #The quantifier object contains a method
                # Which automatically evaluates if the quantity
                # of patterns matched fulfills the quantifier.
                funOut=lambda glycan: numericObj.logical_fun(funOut(glycan))
            return(funOut)
        return(_searchDecWrap)
                
    @searchDecorator(negation,numeric)
    def constraintGen(self):
        '''
        Converts constraint sequence into a 
        regular expression function
        '''
        stringSearch=self.createSeq()
        stringSearch_regex=super().rule2regex(stringSearch)
        return lambda glycan: re.finditer(stringSearch_regex,glycan)

    @classmethod
    def fromComponents(cls,ruleSet):
        '''
        Parses constraint rule components within a
        rule set and returns tags used for 
        the constraint constructor.
        '''
        negation=False
        numeric=quantifierToken
        attachment=None
        seqSet=[]
        for i,t in enumerate(ruleSet):
            if t.__name__=='constraintToken':
                if t.constr.__name__=='negationRule_token':
                    negation=True
                elif t.constr.__name__=='quantityRule_token':
                    numeric=True
                    if ruleSet[-1].__name__!='quantifierToken':
                        raise Exception("Quantity rule detected but no quantifier/quantity provided")
                    else:
                        numeric=t.constr
                elif t.constr.__name__=='attachRule_token':
                    attachment=True
            else:
                seqSet.append(t)
        return(cls(negation,numeric,attachment,seqSet))


class constraintRule(Rule):

    def __init__(self):
        super().__init__(ruleComponents)

    ############################
    # Constraint Rule Validation:
    ############################

    @Rule.checkWrapper
    @Rule.allTrueWrap
    def noReactions(self,ruleSet):
        return([x.__name__!='reactionToken' for x in ruleSet])

    @Rule.checkWrapper
    @Rule.possibleTrueWrap
    def validNumeric(self,ruleSet):
        '''
        Numeric constraints must have the quantity rule prefix "n"
        as well as the quantity constraint as the suffix.
        '''
        validQuantityRuleToken=True if ruleSet[0].__name__=='quantityRule_token' else False
        validQuantityToken=True if ruleSet[-1].__name__=='quantifierToken' else False
        return([validQuantityRuleToken,validQuantityToken])

    @Rule.checkWrapper
    @Rule.possibleTrueWrap
    def validNegation(self,ruleSet):
        '''
        Checks for negation rule presence and validity  
        '''
        validNegation=True if ruleSet[0]=='negationRule_token' else False
        return([validNegation])

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
        #Split rule:
        ruleSets,logicalSeps=self.getRuleSets()
        #Check Conditions for each rule set:
        # Default constraints:
        basicValidation=self.basicValidationWrapper(ruleSets)
        # Reaction-specific constraints:
        hasNoReactions=self.noReactions(ruleSets)
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
        #Split rule:
        ruleSets,logicalSeps=self.getRuleSets()
        #Make Constraint objects:
        Constraint_Classes=[ConstraintMethodGenerator.fromComponents(r) for r in ruleSets]
        #Merge into one function:
        def mergeWrapper(acc,cur,sep):
            if sep=='&':
                return lambda in: acc(in) and cur(in)
            elif sep=='|':
                return lambda in: acc(in) or cur(in)
        #Merge the functions together
        constraintMain=functools.reduce(lambda acc,(cur,sep):lambda in: acc(in) and x[0](in) if x[1]=='&' else acc(in) or x[0](in),lambda in:Constraint_Classes[0](in))
        return constraintMain
