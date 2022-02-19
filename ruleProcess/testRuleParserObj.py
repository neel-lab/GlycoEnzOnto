import pandas as pd
from tokenClasses_lex import *
from ruleInterpreter_objects import *
import traceback

def test_wrapper(rule,ruleName):
    '''
    Prints header information and logs
    outpu from 
    '''
    print('\n************TEST: %s************'%(ruleName))
    print(f'Rule String: {rule}')
    try:
        res=lexer(rule)
    except Exception as exc:
        print('PARSING ERROR\nERROR MESSAGE:')
        print(exc)
    try:
        gct=reactionRule(res)
    except Exception as exc:
        print('INSTANTIATION ERROR')
        print(exc)
    print('Parser Results (reading left to right):')
    for r in res:
        print(r)
    print('\n************END************')

finishedGlycogenes=pd.read_csv('../finishedGlycogenes.tsv',sep='\t',index_col=None)
for i,r in finishedGlycogenes.iterrows():
    if r['Rules']=='no reaction':
        continue
    test_wrapper(r['Rules'],r['geneName']+' Rule')
