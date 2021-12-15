import pandas as pd
from tokenClasses_lex import *
from ruleInterpreter import *
import traceback

def test_wrapper(constraint,rule,c_name):
    '''
    Prints header information and logs
    outpu from 
    '''
    print('\n************TEST: %s************'%(c_name))
    print(f'Constraint string: {constraint}')
    try:
        res=lexer(constraint)
    except Exception as exc:
        print('PARSING ERROR\nERROR MESSAGE:')
        print(exc)
    try:
        rct=reactionRule(lexer(rule))
        gct=constraintRule(res,rct)
    except Exception as exc:
        print('INSTANTIATION ERROR')
        print(exc)
    print('Parser Results (reading left to right):')
    for r in res:
        print(r)
    print('\n************END************')

finishedGlycogenes=pd.read_csv('../finishedGlycogenes.tsv',sep='\t',index_col=None)
for i,r in finishedGlycogenes.iterrows():
    if r['Constraints']=='None':
        continue
    test_wrapper(r['Constraints'],r['Rules'],r['geneName']+' Constraint')
