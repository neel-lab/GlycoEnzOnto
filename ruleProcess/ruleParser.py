import os
import sys
import pandas as pd
#tokenClasses imports tokenMatchers
from tokenClasses import *


#########################
# Custom Error Classes
#########################

class MultipleTokensError(Exception):

    def __init__(self,idx,toks):
        self.index=idx
        self.tokenList=toks
        self.message="""Multiple tokens matched at position %d in ruleString.\nTokens: %s""" %(self.index,self.tokenList)
        super().__init__(self.message)

class ruleStringParsingError(Exception):

    def __init__(self,ruleString,uk_startList,uk_endList):
        self.ruleString=ruleString
        self.uk_start=uk_startList
        self.uk_end=uk_endList
        #Build message using interal method:
        self.message=self.error_message()
        super().__init__(self.message)

    def error_message(self):
        msg_top="Couldn\'t understand underlined parts of ruleString:"
        msg_ruleString=self.ruleString
        msg_error=self.error_underline()
        msg='\n'.join([msg_top,'\n',msg_ruleString,msg_error])
        return(msg)

    def error_underline(self):
        errorStringList=[' ']*len(self.ruleString)
        for s,e in zip(self.uk_start,self.uk_end):
            for i in range(s,e):
                errorStringList[i]='~'
        errorString=''.join(errorStringList)
        return(errorString)

class RuleParser:

    def __init__(self,ruleString):
        '''
        Steps through rule strings searching for 
        high-level rule string components.
        
        Searches for the following high-level tokens
        in order:
        - Reaction Token
        - Constraint Token
        - Entity Token
        - Plurality Token
        - Logical AND/OR Tokens
        '''
        self.ruleString=ruleString
        #RuleParser token library:
        #FULL TOKEN LIST:
        #self.token_lib=[additionToken,mono_additionToken,mod_additionToken,subtractionToken,substitutionToken,reversibleToken,quantityRule_token,attachRule_token,negationRule_token,quantifierToken,mono_token,wild_token,branch_wild_token,multiToken,or_separator,and_separator]
        #DEV TOKEN LIST
        self.token_lib=[reactionToken,constraintToken,entityToken,multiToken,logicalToken]
        
    def search_cur_pos(self,idx):
        '''
        Finds matches at current position of ruleString
        using token library
        '''
        #Start reading rule string from current index:
        ruleString_frag=self.ruleString[idx:]
        #Search the current token libraries:
        goodTokens=[]
        for tok_c in self.token_lib:
            try:
                tok_i=tok_c(ruleString_frag)
            except:
                continue
            if tok_i.detectFun():
                goodTokens.append(tok_i)
        return(goodTokens)

    def parseMain(self):
        '''
        Main method to parse rule string into a list of 
        tokens.

        Uses ruleString input and increments along the
        ruleString to identify tokens.

        If no tokens are identified, the current index
        is incremented by 1 until a match is found.

        Portions of the rule string that aren't recognized
        are returned as "unknownTokens".  Warnings are 
        returned if this is the case.
        '''
        #Current index initialized to 0
        idx=0
        #Resolution flag: are we recognizing the string 
        # currently?
        resolved=True
        #Unknown location lists for error logging:
        uk_start=[]
        uk_end=[]
        #Output Token List:
        tokens=[]
        while idx<len(self.ruleString):
            #Call the token search function:
            tokList=self.search_cur_pos(idx)
            # Confirm the token list has one 
            # element.
            # - Cases where two tokens are matched
            #   will result in an error.
            # - An empty list will trigger an 
            #   increment of +1 to the current index:
            if len(tokList)>1:
                raise MultipleTokensError(idx,tokList)
            elif len(tokList)==0:
                #Catalog the position where no
                # match was found:
                # Increment by 1, then search again
                resolved=False
                uk_start.append(idx)
                idx+=1
            else:
                # If a match was found after
                # not being resolved, set
                # resolved to True and append
                # the current index to the uk_end 
                # list.
                if resolved is False:
                    resolved=True
                    uk_end.append(idx)
                #Get matched token:
                token=tokList[0]
                #Get matched position using token matchFun method:
                sch=token.matchFun()
                sch_end=sch.end()
                #Append the found token to the token list:
                tokens.append(token)
                #Increment starting position by
                # distnace to the end of match
                idx+=sch_end
            
        #If uk_start has an index and uk_end
        # is missing a paired index, means the 
        # end of the string was reached.  End 
        # index is also unknown:
        if len(uk_start)>0 and len(uk_end)==0:
            uk_end.append(len(self.ruleString))
        #If there were any unknown regions 
        # in the ruleString, raise a
        # ruleStringParsingError.
        #Otherwise, return the list of tokens
        if len(uk_start)>0:
            raise ruleStringParsingError(self.ruleString,uk_start,uk_end)
        else:
            return(tokens)
