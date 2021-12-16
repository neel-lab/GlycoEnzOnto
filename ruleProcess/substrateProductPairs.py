import pandas as pd
from tokenClasses_lex import *
from ruleInterpreter import *
import sys
from sys import argv
sys.path.insert(0,'../glycan_rdf')
finishedGlycogenes=pd.read_csv('../finishedGlycogenes.tsv',sep='\t',index_col=None)

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

#Make dict:
objDict=dict()
for _,r in finishedGlycogenes.iterrows():
    if r['Rules']=='no reaction' or test_wrapper(r['Rules']) is None:
        continue
    try:
        objDict[r['geneName']]=reactionRule(lexer(r['Rules']))
    except:
        objDict[r['geneName']]=None

import argparse
parser=argparse.ArgumentParser(description='Provide glycogene name, get table of substrate/product pairs')
parser.add_argument('--ggene',type=str,help='HGNC symbol of glycogene goes here')
parser.add_argument('--print-all',help='Prints substrate/product pairs for all glycogenes which have been processed into reaction rules.  Records are delimited by \"===\" markers.',type=bool)

def rulePairPrinter(ggene_obj):
    #Get rule Pair
    rulePairs=[ggene_obj.pairListBuilder(x) for x in ggene_obj.ruleSets]
    for r in rulePairs:
        for p in r:
            print('\t'.join([p[0],p[1]]))

def printAllRules(objDict):
    for k,v in objDict.items():
        #Print Gene Name Header:
        print(k)
        #Print rule pair 
        rulePairPrinter(v)
        print('===')

if __name__=='__main__':
    args=parser.parse_args(sys.argv[1:])
    if args.print_all:
        try:
            printAllRules(objDict)
        except:
            raise Exception('Something is wrong with one of the glycogene entries')
    elif args.ggene:
        try:
            rulePairPrinter(objDict[args.ggene])
        except:
            Exception('The glycogene entered does not have a valid rule yet.  Choose another.')
    else:
        print('Please pass a glycogene name.')
