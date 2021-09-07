import re

#########################################
# Rule and Attribute Class Matching Class
#########################################
class matcherClass:

    def __init__(self,regex,front=True):
        '''
        The matcherClass stores regular expression
        information used to lexicographically 
        parse rule strings and extract rule tokens.
        
        The matcherClass will use the "regex" variable
        to find matches within the rule string, starting
        at the beginning of the string.

        '''
        self.matchString=regex
        #Append "^" for stepping into the rule string.
        if front:
            self.matchString='^'+regex
        else:
            self.matchString=regex

    def search(self,string):
        return(re.search(self.matchString,string))

    def isPresent(self,string):
        searchRes=self.search(string)
        return(True if searchRes is not None else False)

    def __call__(self,string,presence=True):
        if presence:
            return(self.isPresent(string))
        else:
            return(self.search(string))

############################
# Lexion-based Matcher Class
############################

class LexMatcher:

    def __init__(self,*args,regex=None,front=True):
        '''
        The LexMatcher class creates a regular expression
        function which searches for expected values 
        provided in an input list.

        An optional "regex" variable can be passed to a
        LexMatcher where the lexicon can be used with other 
        regular expressions. Use string replacement syntax as
        shown below.

        Example:
        INPUT
        Lexicon: [Gal,Glc,GlcNAc]
        Regex: %s\[.+?\]
        front=True
        OUTPUT
        Search String: "^(Gal|Glc|GlcNAc)\[.+?\]"
        '''
        self.lexicons=args
        self.regex=regex
        self.lexRegex=self.make_regex()
        #Append "^" for stepping into the rule string.
        if front:
            self.lexRegex='^'+self.lexRegex
    
    def make_regex(self):
        '''
        Creates regular expression string to
        match all values.

        If "regex" parameter passed to class,
        processes the pattern to include the 
        lexicon in the right place.
        '''
        #Create lexicon patterns for every kwargs input:
        lexPatterns=['('+'|'.join([re.escape(i) for i in x])+')' for x in self.lexicons]
        #lexPattern='('+'|'.join([re.escape(x) for x in self.lexicon])+')'
        if self.regex is not None:
            #Replace "%s"s with regex patterns in "lexPatterns":
            regex=self.regex %(tuple(lexPatterns))
            return(regex)
        else:
            if len(lexPatterns)>1:
                raise Exception('Multiple lexicons passed, must provide a regex string')
            else:
                lexPattern=lexPatterns[0]
            return(lexPattern)

    def search(self,string):
        return(re.search(self.lexRegex,string))

    def isPresent(self,string):
        searchRes=self.search(string)
        return(True if searchRes is not None else False)

    def __call__(self,string,presence=False):
        if presence:
            return(self.isPresent(string))
        else:
            return(self.search(string))

####################
# Reaction Matchers:
####################
reactionMatcher=matcherClass('\{.+\}')
#innerReactionMatcher used to find reactions within
# monosaccharide entities (for modification reactions)
innerReactionMatcher=matcherClass('\{',front=False)
additionMatcher=matcherClass('\{(?!\!).*?\}')
subtractionMatcher=matcherClass('\{(?=\!).*?\}')
substitutionMatcher=matcherClass('\{.*?\-\>.*?\}')
reversibleMatcher=matcherClass('\{.*?\<\-\>.*?\}')
###################
# Entity Detection:
###################
multiMatcher=matcherClass('\<.+\>')
wildCardMatcher=matcherClass('\[?\.\.\.\]?')
#Monosaccharide-linkage matcher:
#Match general monosaccharide pattern
# ensure no { } in front of the string:
monoLinkMatcher=matcherClass('(?!\{)\[?([A-Za-z0-9\,\{\}]+?)(\([ab\?][12\?]\-[\d\?]\))\]?')
# Excludes any special characters assoicated
# with reaction and constraint string (n,@,{},<>)
#monoLinkMatcher=matcherClass('(?!n|\@|\{|\}|\<|\>)\[?(.+?)(\([ab\?][12\?]\-[\d\?]\))\]?')
#######################
# Constraint Detection:
#######################
quantityStartMatcher=matcherClass('n')
attachRuleMatcher=matcherClass('\@')
negationRuleMatcher=matcherClass('\!')
quantifierMatcher=matcherClass('(\=|\>\=|\<\=|\>\=)(\d)')
#######################
# Separator Detection
#######################
orMatcher=matcherClass('\|')
andMatcher=matcherClass('\&')
