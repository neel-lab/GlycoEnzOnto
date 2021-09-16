import re
import functools
from functools import reduce as reduct

################
# Test pair list
################

pairList=[('a','A'),('b','B'),('c','C'),('d','D'),('e','E')]

###############################################
# Reaction Search/Replace Generating Functions:
###############################################

def makeGlycanProcessor(frm,to):
    '''
    Wrapper to dynamically create instances
    of "proc_" with "frm" "to" pairs.
    '''
    return lambda glycan: proc_(glycan,frm,to)

def proc_(glycan,frm,to):
    '''
    Finds fragments in "glycan" where "frm" is located
    and replaces each with text in "to"
    '''
    mtchs=re.finditer(frm,glycan)
    #Initialize a product list:
    products=[]
    for m in mtchs:
        #Constant indicies:
        front_start=0
        front_end=m.start()
        back_start=m.end()
        back_end=len(glycan)
        #Get Match Fragment:
        p=m.group()
        #Replace text where "frm" was found with the "to" text:
        products.append(''.join([glycan[front_start:front_end],to,glycan[back_start:back_end]]))
    return(products)

def glycanProcAggregator(fun):
    def _wrap(pairList):
        '''
        Returns a function that converts glycans 
        based on given "frm" "to" strings.
        "fun" reorders the pairList to perform
        "forward" or "reverse" inference.
        '''
        #Reordering the pairList to return either
        # substrate->product or product->substrate.
        frm_list,to_list=fun(pairList)
        #Make a list of replacer functions:
        funList=[makeGlycanProcessor(frm,to) for frm,to in zip(frm_list,to_list)]
        #Make conversion function:
        convertMain=reduct(lambda cur,pres: lambda string: cur(string)+pres(string),funList)
        return(convertMain)
    return(_wrap)

######################################
# Forward/Reverse Inference Functions:
######################################

@glycanProcAggregator
def forwardInferenceGenerator(pairList):
    '''
    Produces a list of all possible glycans where "substrate"
    is taken to "product"
    '''
    return([t[0] for t in pairList],[t[1] for t in pairList])

@glycanProcAggregator
def reverseInferenceGenerator(pairList):
    '''
    Produces a list of all possible glycans where "product"
    is taken to "substrate"
    '''
    return([t[1] for t in pairList],[t[0] for t in pairList])


########################
# Constraint Generators:
########################

def makeConstraintProcessor(



#############################
# Inferrence Generator Tests:
#############################

forwardInferelator=forwardInferenceGenerator(pairList)
reverseInferelator=reverseInferenceGenerator(pairList)

forwardTestSeq="""Hello, my name is Inigo Montoya, you killed my father, prepare to die."""
reverseTestSeq="""Agile Beavers Conjure Dastardly Evils"""

# These inferelators should output a list of sequences
# which modify the input by only one recognition sequence.
print(forwardInferelator(forwardTestSeq))
print(reverseInferelator(reverseTestSeq))

#############################
# Constraint Generator Tests:
#############################
