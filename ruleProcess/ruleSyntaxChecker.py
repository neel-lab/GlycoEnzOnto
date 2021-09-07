#!/usr/bin/python3
import os
import sys
import re
import argparse

def argParser():
    '''
    Parses command line arguments for syntax checker.
    '''
    parser=argparse.ArgumentParser(description="Glycogene reaction rule syntax checker")
    parser.add_argument('--rule',type=str,help='Reaction rule string to check')
    parser.add_argument('--explain',help='Prints verbal explanation of reaction rule')
    return(parser)

def split_rule(ruleString):
    '''
    Split rules by their major components and return a list:
    The regular expression matches these following patterns in
    order:
     1. Any uncertain monosaccharide/linkage information (\<.+\>)
     2. Any reaction rule contents: (\{.+?\})
     3. Any uncertain linkage information \[?\.\.\.\]?
     4. Any Monosaccharide-linkage string (\[?.+?\([ab][12]\-\d\)\]?)
      4a. Any quantifier information found in a constraint ((?:(?:\=|\>\=|\<\=|\>\=)\d)?)
    These pieces are returned as a list and are used in 
    substrate-product matching or for validating rule strings
    '''
    return(re.split(r'(?:\<.+\>|\{.+?\}|\[?\.\.\.\]?|\[?.+?\([ab][12]\-\d\)\]?|(?:\=|\>\=|\<\=|\>\=)\d)',ruleString))

def rule2regex(rl):
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
    wild=re.search(r'\.\.\.',rl)
    frontwild=re.search(r'^\.\.\.',rl)
    if frontwild is not None:
        rl=re.sub('^\.\.\.','',rl)
    elif frontwild is None:
        rl=re.sub('^','(?:^|\[)',rl)
    elif frontwild is None and wild is not None:
        rl=re.sub('\.\.\.','.*?',rl)
    # Uncertain linkages:
    rl=re.sub(r'\-\?','-[0-9]',rl)
    #If a core linkage is detected, append a "$" to the 
    # string to look at the core:
    if re.search(r'\([ab][12]\-$',rl) is not None:
        rl=re.sub(r'(\([ab][12]\-$)','\g<1>$')
    return(rl)

def ruleProcess(rl):
    '''
    Transforms a rule string into substrate and product string.
    '''
    #Detect Change:
    chg=re.compile(r'\{.*?\}')
    chgMatch=re.search(chg,rl);chgFrag=getFragment(chgMatch)
    if re.search('\{\!',rl) is not None:
        #Is a deletion reaction:
        substr=re.sub('\{\!|\}','',rl)
        prod=re.sub(chg,'',rl)
    elif re.search('\-\>',chgFrag) is not None:
        #Is a substitution reaction:
        from_to_search=re.compile('\{(.*)\-\>(.*)\}')
        fromMono=re.sub(from_to_search,'\g<1>',rl)
        toMono=re.sub(from_to_search,'\g<2>',rl)
        substr=re.sub(chg,fromMono,rl)
        prod=re.sub(chg,toMono,rl)
    else:
        #Is an addition reaction:
        if (re.search('^\[',chgFrag) is not None) and (re.search('\]$',chgFrag) is None):
            #The monosaccharide was added onto an existing branch.
            #Substrate string needs a '[' appended to where the
            # insertion ends
            appends='['
        else:
            appends=''
        substr=''.join([rl[0:chgMatch.start()],appends,rl[chgMatch.end():len(rl)]])
        prod=re.sub('\{|\}','',rl)
    return(substr,prod)

def make_rule_regex(rl):
    '''
    Wraps the rule processing and the regular expression
    processing in one function:
    '''
    substr,prod=ruleProcess(rl)
    substr=rule2regex(substr)
    prod=rule2regex(prod)
    return substr,prod

