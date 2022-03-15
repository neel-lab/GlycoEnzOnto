#GlycoEnzOnto

GlycoEnzOnto is an ontology which describes the molecular functions and metabolic pathway involvement of 386 glycosylation-related genes, called glycogenes.  There are four main kinds of glycosylation-related pathways modeled in modelled in GlycoEnzOnto:

1. Metabolism of glycosylation donor substrates: Donor substrates are classified as any molecule that is added onto a glycan structure.  These molecules include:
  * Nucleotide Sugars, and
  * Monosaccharide substituents, which includes sulfate, phosphate, and ethanolamine for example
2. Transport of donor substrates into glycosylation compartments such as the endoplasmic reticulum and Golgi
3. Glycan biosynthesis: these pathways describe the sequences of metabolic reactions whereby glycans are synthesized.  There are three phases of glycan biosynthesis modeled in GlycoEnzOnto:
  * Core structure biosynthesis: these pathways are responsible for creating the core structures for various glycan types (*N*-linked, *O*-linked, Glycosphingolipids, Glycosaminoglycans, *etc.*)  Enzymes involved in these pathways are unique to the kinds of glycans they synthesize.
  * Elongation and Branching pathways: these pathways increase the size and degree of branching of glycan structures, and can take place on a variety of glycan types. These pathways include synthesis of LacNAc structures, as well as the glycosaminoglycan polymerization reactions.   
  * Capping pathways: these pathways transfer monosaccharides and substituents onto glycans such that, after being added, precludes any further attachment of monosaccharides onto a glycan structure.

4. Degradation of glycans:  These pathways describe how the various glycan types are degraded.  These classes are defined in terms of i. the type of glycan being degraded, and ii. the compartment in which the glycan is degraded.

Auxilliary processes which serve to modulate the glycosylation process or substrates for glycosylation are also modelled in GlycoEnzOnto.  The mechanistic functions of each enzyme are encoded as a single string in the form of reaction rules and reaction constraints.  Included alongside our glycoenzyme pathway classification are all biological process (BP), molecular function (MF), and cellular compartment (CC) annotations from the gene ontology (GO) to provide interoperability in function classification, as well as provide cross references to other resources.  GlycoEnzOnto pathways which were deemed to be encompassed by a GO biological process term were related with the "occurrent part of" relationship in the Relation Ontology (RO).
