from ruleParser import *
import logging


def test_wrapper(rule,ruleName):
    '''
    Prints header information and logs
    outpu from "ruleParser"
    '''
    print('\n************TEST: %s************'%(ruleName))
    print(f'Rule String: {rule}')
    #Rule Parser Instance:
    rp=RuleParser(rule)
    try:
        res=rp.parseMain()
        print('Parser Results (reading left to right):')
        for r in res:
            print(r)
    except Exception as exc:
        print('ERROR MESSAGE:')
        print(exc)
    print('\n************END************')


#1. Addition Reaction: 
B3GNT6_rule='{GlcNAc(b1-3)}GalNAc(a1-?)'
test_wrapper(B3GNT6_rule,'Addition Examle (B3GNT6)')

#2. Addition Reaction with multiple Reactions: 
B3GNT5_rule='{GlcNAc(b1-3)}Gal(b1-4)Glc(b1-?)|{GlcNAc(b1-3)}Gal(b1-4)'
test_wrapper(B3GNT6_rule,'Addition Examle with Multiple Reactions: (B3GNT5)')

#3. Subtraction Reaction: 
FUCA1_rule='...{![Fuc(a1-6)]}GlcNAc(b1-?)'
test_wrapper(FUCA1_rule,'Subtraction Example with uncertain reducing end (FUCA1)')

#4. Substitution Reaction:
HS6ST2_rule='IdoA(a1-4){GlcNS(a1-4)->GlcNS,6OS(a1-4)}|GlcA(a1-4){GlcNS(a1-4)->GlcNS,6OS(a1-4)}'
test_wrapper(HS6ST2_rule,'Substitution Example, modifications within glycans, multiple rules (HS6ST2)')

#5. Constraint Test:
B3GALT1_constr='!Gal(b1-4)GlcNAc(b1-?)&!@GlcNAc(b1-4)...Man(b1-4)'
test_wrapper(B3GALT1_constr,'Constraint Example, with logical NOT,AND, attachment constraint, and uncertain linkage. (B3GALT1)')

#6. Multi-entity Test:
multi_test='<Gal(b1-3),GalNAc(b1-3)>GlcNAc(b1-?)'
test_wrapper(multi_test,'Multiple substrates at position test (arbitrary rule)')

#6. Parsing Error
broken_rule_1='Gal(aaa)GlcNAc(a1-4)'
test_wrapper(broken_rule_1,'Typographical error: monosaccharide-linkage formatting')

#7. Parsing Error
broken_rule_2='Gal(b1-4)GlcNAc(b1-3)[Fuc(a1-3)]}...'
test_wrapper(broken_rule_2,'Typographical error: Misformatted rule')
