import pandas as pd
from tokenClasses_lex import *
from ruleInterpreter import *
import traceback

#finisehdGlycogenes.tsv:
finishedGlycogenes=pd.read_csv('../finishedGlycogenes.tsv',sep='\t',index_col=None)
#Create glycogene objects:
#  Methods:
#    Forward and Reverse reaction inference methods
#    ConstraintMethod

class glycogene:

    def __init__(self,name,ruleString,constraintString):
        '''
        Glycogene Class:
            Takes reaction rule strings and constraint
            strings as input and creates
            the following methods:
        - forward: takes a glycan substrate as an argument and 
                   predicts a set of products, returned as a list.
        - reverse: takes a glycan product and infers substrate where
                   it was derived, returned as a list.
        - constraint: Internal method which checks reactant and product
                      strings to validate the reaction takes place
        '''
        #Glycogene Name:
        self.name=name
        try:
            self.reactionRule=reactionRule(lexer(ruleString))
            #Define the "forward" and "reverse" methods:
            self.forward,self.reverse=self.reactionRule.forward,self.reactionRule.reverse
        except:
            raise Exception(f"Could not create reaction rule method for {self.name}")
        try:
            if constraintString=="None":
                self.constraint=lambda x:True
            else:
                #Define "constraint" method:
                self.constraintRule=constraintRule(lexer(constraintString),self.reactionRule)
                self.constraint=self.constraintRule.constraint
        except:
            raise Exception(f"Could not create constraint rule method for {self.name}")

###########################
# Create Glycogene Objects:
###########################

#Keep track of unprocessed glycogenes:
# No reactions:
noRuleGlycogenes=[]
# Unprocessed:
notProcessed=[]
# Processed: 
processed=[]
 
ggenes=dict()
for i,r in finishedGlycogenes.iterrows():
    ruleString=r['Rules']
    constraintString=r['Constraints']
    geneName=r['geneName']
    if r['Rules']=='no reaction':
        noRuleGlycogenes.append(geneName)
        continue
    try:
        ggenes[geneName]=glycogene(geneName,ruleString,constraintString)
        #ggene=glycogene(geneName,ruleString,constraintString)
        processed.append(geneName)
    except :
        notProcessed.append(geneName)

####################
# Processing Report:
####################

print('Processing Summary:')
print(f'{len(noRuleGlycogenes)} do not have reaction rules and were not processed.')
print(f'Of the {len(processed)+len(notProcessed)} glycogenes:\n{len(processed)} were successfully processed.')
print(f'Unprocessed genes: {notProcessed}')
