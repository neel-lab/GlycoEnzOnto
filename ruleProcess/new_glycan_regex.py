import re

######################
# Rule Regex Converter
######################

def rule2regex(rl):
    '''
    Converts reaction rules into regex string.  Makes
    syntax compatible with regular expressions.
    '''
    # Uncertain linkages:
    uncertainty=re.search('\(([\D\?])([\d\?])\-([\d\?])\)',rl)
    if uncertainty is not None:
        if uncertainty.groups()[0]=='?':
            pos1='[ab]'
        else:
            pos1=uncertainty.groups()[0]
        if uncertainty.groups()[1]=='?':
            pos2='[12]'
        else:
            pos2=uncertainty.groups()[1]
        if uncertainty.groups()[2]=='?':
            pos3='[0-9]'
        else:
            pos3=uncertainty.groups()[2]
    #Reassembl
    stereo_u=re.search('\(\?[12]\-[0-9]\)',rl)
    number_u=re.search('\([ab]\?\-[0-9]\)|\([ab][0-9]\-\?\)')
    if stereo_u is not None:
        rl=re.sub(
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
    else:
        rl=re.sub('^','(?:^|\[)',rl)
    if frontwild is None and wild is not None:
        rl=re.sub('(?!^)\.\.\.','.+?',rl)
    rl=re.sub(r'\-\?','-[0-9]',rl)
    return(rl)
