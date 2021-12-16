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
    except Exception as exc:
        print('PARSING ERROR\nERROR MESSAGE:')
        print(exc)
    if any([type(x)==unknownToken for x in res]):
        print('PARSING FAILURE')
    print('Parser Results (reading left to right):')
    for r in res:
        print(r)
    print('\n************END************')

finishedGlycogenes=pd.read_csv('../finishedGlycogenes.tsv',sep='\t',index_col=None)
for i,r in finishedGlycogenes.iterrows():
    if r['Rules']=='no reaction':
        continue
    test_wrapper(r['Rules'],r['geneName']+' Rule')
