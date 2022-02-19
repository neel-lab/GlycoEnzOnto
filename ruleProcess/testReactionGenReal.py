import re
import functools
from dataclasses import dataclass 

################
# Test pair list
################

pairList=[('Gal(b1-4)GlcNAc(b1-?)','Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-?)')]
def pairListGenerator():
    return([('Gal(b1-4)GlcNAc(b1-?)','Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-?)')])

###############################################
# Reaction Search/Replace Generating Functions:
###############################################

@dataclass
class GlycanProcessGenerator:
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
            rl=re.sub('^','(^|\[)',rl)
        elif frontwild is None and wild is not None:
            rl=re.sub('\.\.\.','.+?',rl)
        # Uncertain linkages:
        rl=re.sub(r'\-\?','-[0-9]',rl)
        #If a core linkage is detected, append a "$" to the 
        # string to look at the core:
        if re.search(r'\([ab\?][12\?]\-$',rl) is not None:
            rl=re.sub(r'(\([ab][12]\-$)','\g<1>$')
        return(rl)

    def makeToRepString(self,fromWildGrp):
        toRepString=re.sub('\(\.\+\?\)',fromWildGrp,self.toString)
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
        fromString_regex=self.rule2regex(self.fromString)
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


def makeGlycanProcessor(frm,to):
    '''
    Wrapper to dynamically create instances
    of "GlycanProcessGenerator" classes with "frm" "to" pairs.
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
        pairList=pairListGenerator()
        #Reordering the pairList to return either
        # substrate->product or product->substrate.
        frm_list,to_list=fun(pairList)
        #Make a list of replacer functions:
        funList=[makeGlycanProcessor(frm,to) for frm,to in zip(frm_list,to_list)]
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


#############################
# Inferrence Generator Tests:
#############################

forwardInferelator=forwardGeneratorMain(pairList)
reverseInferelator=reverseGeneratorMain(pairList)

testGlycan='Gal(b1-4)GlcNAc(b1-2)[Gal(b1-4)GlcNAc(b1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-?)'

# These inferelators should output a list of sequences
# which modify the input by only one recognition sequence.
print(forwardInferelator(testGlycan))
