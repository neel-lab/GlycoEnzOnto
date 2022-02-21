# Rule Processing Repository:

This repository contains code which parses the GlycoEnzOnto reaction rules and constraints.  The tokenization relies on a series of pattern matching functions defined in the ```tokenMatchers.py``` file.  There are five types of tokens recognized by the parser, and each has its own set of subtoken classes to match particular kinds of strings.


* **Reaction Tokens:** Addition,Subtraction,Substitution,Reversible,Transport
* **Constraint Token** Numerical Constraint,Logical Negation,Addition Location Constraint.
* **Entity Token** Any of the below entities will be detected:
  - Monosaccharides
  - Sugar Nucleotides
  - Modifications (Sulfation, Phosphorylation)
  - Compartments ( cis/medial/trans Golgi, ER, Lysosome)
  - Aglyca (Dol-P-P, Asn, Ser/Thr,Cer)
* **Plurality Token** Any uncertainty pattern in rule strings:
  - One or more monosaccharides: "..."
  - One or more monosaccharides on a branch: "[...]"
  - Multiple Token Matching : <  >
    + If the Multiple match token contains one entity, then the entity is optionally
      present.  Otherwise if there is more than one entity, then the entity could be
      any of the entities list, but not none.
* **Logical AND/OR Tokens** Either the "|" or "&" strings matched in a rule.

The ```LexerClass.py``` file contains a class which iteratively searches through a rule string and creates a list of all tokens detected.  The tokens are python classes, all of which are defined in ```tokenClasses_lex.py```.  The subtoken classes inherit methods from their respective superclasses such that token representations and substrate/product string generation is consistent across tokens.  

Each major token class has a common set of methods.  These are described below:

*```__repr__```: This method prints a plain English description of what the token represents.  This description can be viewed if a token is ```print```ed in a python console.
* ```detectFun```: This method searches for strings at the current lexing point of a rule. It utilizes a matching function created in the ```tokenMatchers.py``` file to achieve this.  These methods return True if a string is detected at the front of the rule string.
* ```matchFun```: This method returns a matched pattern corresponding to the token type.  This uses the same matching function as the ```detectFun``` method, but instead returns a string by setting the ```presence``` parameter equal to False.
* ```substrate/product```: These methods print the substrate or product representation of a token.  If the token is a reaction token type, additional methods will be called to detect what kind of entity exists inside the curly braces.  This entity token will be instantiated within the reaction token.  
  * If an addition reaction, nothing is printed in the substrate case, and the entity string is printed in the product case.
  * If a subtraction reaction, the entity representation is printed in substrate case, and nothing in the product case.
  * If a substitution,reversible, or transport reaction, the entity on the left hand side of the "->" is printed in the substrate case and the entity on the right hand side in the product case.

Subtokens of the five major token types each have these methods, and additionally have their own methods for parsing string types which they are supposed to inherit. 


## Running the Rule Parser:

Running the rule parser requires sourcing the ```tokenClasses_lex.py``` file, which instantiates all token classes, and builds a ```lexer``` function which reads rule strings and returns token lists.

### Rule Parser Example:

```
from tokenClasses_lex import *

ruleString="{Gal(b1-4)}GlcNAc(b1-?)...Man(b1-4)GlcNAc(b1-4GlcNAc(b1-?)Asn"
#Run the lexer function:
lexer(ruleString)
```

Running this script results in the following list:

```
[An addition of [A Gal monosaccharide with a (b1-4) linkage],
 A GlcNAc monosaccharide with a (b1-?) linkage,
 One or more monosaccharides,
 A Man monosaccharide with a (b1-4) linkage,
 A GlcNAc monosaccharide,
 An unknown entity: (b1-4,
 A GlcNAc monosaccharide with a (b1-?) linkage,
 The Asn aglycon]
```

## Creating Reaction Rule Methods:

Once glycoenzyme reaction rules and constraints are tokenized, reaction rule and constraint processing objects can be created with classes defined in the `ruleInterpreter.py` file.  Reaction rule tokens are processed in the `reactionRule` object, and constraint tokens are processed in the `constraintRule` object.  Each of these classes is a subclass of the `Rule` superclass.  Below each class is described in more detail:

* `Rule`: This is the superclass from which `reactionRule` and `constraintRule` inherit.  The main function of `Rule` is to separate tokens into groups if logical separators are present in the rule string.  Thus, if an or "|" operator is present, two lists of tokens would be generated.  Two class variables,`ruleSets`  and `logicalSeps` store the separate rule tokens and logical tokens separately.
* `reactionRule`: This class processes sets of tokens extracted from reaction rules.
  + Firstly, the `Rule` method to extract sets of reaction rule tokens and logical separators is run to get the `ruleSets` and `logicalSeps`.
  + Next, a series of token checking processes are run to see if all tokens being passed to the reaction rule object are valid.  This involves:
    i. checking if all tokens have no parsing errors and include known entities (`basicValidationWrapper`), and
    ii. makes sure there are no constraint tokens in the reaction rule (`noConstraints`).
  + Next, the type of reaction is parsed from reaction rule tokens and is stored in a variable called `reactionType`.
  + Finally, two methods for inferring products from substrates (`forward`) and substrates from products (`reverse`) are generated using the `forwardGeneratorMain` and `reverseGeneratorMain` methods, respectively.
* `constraintRule`: Similarly to the `reactionRule` class, this class processes constraint strings.
  * Firstly, the `Rule` method to extract sets of reaction rule tokens and logical separators is run to get the `ruleSets` and `logicalSeps`.
  * Next, a series of token checking processes are run to see if all tokens being passed to the reaction rule object are valid.  This involves:
    i. checking if all tokens have no parsing errors and include known entities (`basicValidationWrapper`),
    ii. makes sure there are no reaction rule tokens in the reaction rule (`noConstraints`), 
    iii. checks if negation and quantity constraints are instantiated correctly.
  * Finally, if the constraint involves a reaction site constraint ("@"), then a corresponding `reactionRule` object is passed to the `constraintRule` class to extract the monosaccharide entity being added.  The `ConstraintGenerator_Aggregator` method is then executed to dynamically create a method which validates glycan structures as valid substrates. 
  
## The `glycoenzyme` class:

This class uses the tokenizer, reaction rule, and constraint rule objects described above and encapsulates them in an object which models an instance of a glycoenzyme.  These classes allow one to process rule strings into an object and use it for substrate/product inference. 

### Reaction Rule and Constraint Rule Example:

```
#Importing glycoenzyme class dependency:
from glycogeneObjs import *

#Reaction rule and constraint for MGAT4A from GlycoOnto:
mgat4a_reactionRule_string="GlcNAc(b1-2){[GlcNAc(b1-4)]}Man(a1-3)"
mgat4a_constraintRule_string="!Gal(b1-?)&!GlcNAc(b1-4)[…]Man(b1-4)"
#Creating MGAT4A object:
mgat4a=glycoenzyme(ruleString=mgat4a_reactionRule_string,constraintString=mgat4a_constraintRule_string)
# According to constraints, any N-linked glycan with
# terminal galactose should not be a substrate of 
# MGAT4A:
validSubstrate='GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-1)'
invalidSubstrate='Gal(b1-4)GlcNAc(b1-2)[GlcNAc(b1-4)]Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-1)'
```
