import re

######################
# Rule Regex Converter
######################

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

