#########################
# Custom Error Classes
#########################

class LexerClass:

    def __init__(self,lexicon,ukToken):
        '''
        Steps through rule strings searching for 
        high-level rule string components.
        
        Searches for the following high-level tokens
        in order:
        - Reaction Token
        - Constraint Token
        - Entity Token
        - Plurality Token (uncertainty)
        - Logical AND/OR Tokens
        '''
        self.lexicon=lexicon
        self.ukToken=ukToken
        
    def search_cur_pos(self,inputString,idx):
        '''
        Finds matches at current position of inputString
        using token library
        '''
        #Start reading rule string from current index:
        inputString_frag=inputString[idx:]
        #Search the current token libraries:
        goodTokens=[]
        for tok_c in self.lexicon:
            try:
                tok_i=tok_c(inputString_frag)
            except:
                continue
            if tok_i.detectFun():
                goodTokens.append(tok_i)
        return(goodTokens)

    def parseMain(self,inputString):
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
        m_start=[]
        m_end=[]
        #Output Token List:
        tokens=[]
        while idx<len(inputString):
            #Call the token search function:
            tokList=self.search_cur_pos(inputString,idx)
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
                if resolved==True:
                    m_start.append(idx)
                resolved=False
                idx+=1
            else:
                # If a match was found after
                # not being resolved, set
                # resolved to True and append
                # the current index to the uk_end 
                # list.
                if resolved is False:
                    resolved=True
                    m_end.append(idx)
                    if self.ukToken is not None:
                        tok=self.ukToken(inputString[m_start[-1]:m_end[-1]])
                        tokens.append(tok)
                    else:
                        tokens.append(None)
                #Get matched token:
                tok=tokList[0]
                #Mark the match start:
                m_start.append(idx)
                #Get matched position using token matchFun method:
                sch=tok.matchFun()
                sch_end=sch.end()
                ### Branching check for 2nd stub:
                #  If the last added mono has left square bracket
                #  and current mono has right bracket, set the mono
                #  "isextendbranch" to "True"
                #if len(tokens)>0:
                #    if tokens[-1].__name__=='reactionToken':
                #        if tok.token.__name__=='monoToken':
                #            if tokens[-1].token.ligand_token[0].token.branching['leftBracket'] and tok.token.branching['rightBracket']:
                #                tok.token.isextendbranch=True
                #Append the found token to the token list:
                tokens.append(tok)
                #Mark the match end:
                m_end.append(idx+sch_end)
                #Increment starting position by
                # distnace to the end of match
                idx+=sch_end
            
        #If uk_start has an index and uk_end
        # is missing a paired index, means the 
        # end of the string was reached.  End 
        # index is also unknown:
        if len(m_start)>len(m_end):
            m_end.append(idx)
            tok=self.ukToken(inputString[m_start[-1]:m_end[-1]])
            tokens.append(tok)
        #If there were any unknown regions 
        # in the ruleString, raise a
        # ruleStringParsingError.
        #Otherwise, return the list of tokens
        if None in tokens:
            #Find where the None token is
            # then report error:
            missingLocs=[i for x,i in enumerate(tokens) if x is None]
            uk_start,uk_end=[[m_start[i],m_end[i]] for i in missingLocs]
            raise ruleStringParsingError(self.ruleString,uk_start,uk_end)
        else:
            return(tokens)

    def __call__(self,inputString):
        return(self.parseMain(inputString))

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

