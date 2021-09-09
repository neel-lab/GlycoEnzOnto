import pandas as pd
from tokenClasses_lex import *


def test_wrapper(rule,ruleName):
    '''
    Prints header information and logs
    outpu from "ruleParser"
    '''
    print('\n************TEST: %s************'%(ruleName))
    print(f'Rule String: {rule}')
    try:
        res=lexer(rule)
        print('Parser Results (reading left to right):')
        for r in res:
            print(r)
    except Exception as exc:
        print('ERROR MESSAGE:')
        print(exc)
    print('\n************END************')

finishedGlycogenes=pd.read_csv('../finishedGlycogenes.tsv',sep='\t',index_col=None)
for i,r in finishedGlycogenes.iterrows():
    if r['NewRule']=='no reaction':
        continue
    test_wrapper(r['NewRule'],r['geneName']+' Rule')
