import re

##########################
# Rule String Exceptions
##########################

class ConflictingRuleException(Exception):

    def __init__(self):
        self.message="Reaction and Constraint component detected in same rule string, this is not allowed"
        super().__init__(self.message)

class ConflictingRuleException(Exception):

    def __init__(self):
        self.message="Reaction and Constraint component detected in same rule string, this is not allowed"
        super().__init__(self.message)

##########################
# Rule Interpreter Class
##########################

class ruleInterpreter:

    def __init__(self,reactionComponents,constraintComponents):
        '''
        Interprets parsed rule string lists by:
        - Checking if combination of token types is allowed.
        - If passed Reaction rule, creates methods for recognizing
          glycan substrates and their conversion into products.
        - If passed a constraint, creates methods to logically evaluate
          substrates/products to determine validity of reaction
        '''
