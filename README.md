# GlycoEnzOnto

**Citation**: Theodore Groth, Alexander D Diehl, Rudiyanto Gunawan and Sriram Neelamegham, "GlycoEnzOnto: A GlycoEnzyme Pathway and Molecular Function Ontology", Bioinformatics. 2022 Oct 25;btac704. ([doi: 10.1093/bioinformatics/btac704](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btac704/6772808))

GlycoEnzOnto is an ontology which describes the molecular functions and metabolic pathway involvement of 403 glycosylation-related genes or "glycogenes". These genes are involved in regulating four different processes, each of which are modeled in GlycoEnzOnto:

1. Glycosylation donor substrate biosynthesis: These metabolic reactions regulate the biosynthsis of "donors". Here, donors are any molecule that provide compounds for glycan biosynthesis.  These donors include:
  * Nucleotide-sugars like CMP-Neu5Ac, GDP-Fucose *etc.*, and
  * Molecules that provide monosaccharide substituents, like PAPS for for sulfate

2. Transport of donor substrates into glycosylation compartments such as the endoplasmic reticulum and Golgi

3. Glycan biosynthesis: these pathways describe the sequence of metabolic reactions whereby glycans are synthesized.  There are three phases of glycan biosynthesis modeled in GlycoEnzOnto:
  * Core structure biosynthesis: these pathways are responsible for creating the core structures for various glycan types (*N*-linked, *O*-linked, Glycosphingolipids, Glycosaminoglycans, *etc.*). Enzymes involved in these pathways are unique to the kinds of glycans they synthesize, and thus they tend to be quite specific in their action.
  * Elongation and Branching pathways: these pathways increase the size and degree of branching of glycan structures, and can take place on a variety of glycan types. These pathways include synthesis of LacNAc structures, as well as the glycosaminoglycan polymerization reactions.   
  * Capping pathways: these pathways transfer monosaccharides and substituents onto glycans such that, after being added, they preclude any further extension of the  glycan.

4. Degradation of glycans:  These pathways describe how the various glycan types are degraded.  These classes are defined in terms of i. the type of glycan being degraded, and ii. the compartment in which the glycan is degraded.

Auxilliary processes which serve to modulate the glycosylation process or substrates for glycosylation are also modeled in GlycoEnzOnto.  The mechanistic functions of each enzyme are encoded as a single string in the form of "reaction rules" and "reaction constraints".  Included alongside our glycoenzyme pathway classification are biological process (BP), molecular function (MF), and cellular component (CC) annotations from the Gene Ontology (GO). This incusion aims to provide interoperability in function classification, as well to provide cross-references to other resources.  GlycoEnzOnto pathways which were deemed to be encompassed by a GO biological process term were related with the "occurrent part of" relationship in the Relation Ontology (RO).

## Enrichment Analysis Tutorial:

Instructions on how to perform GSEA analysis using `GlycoEnzOnto.gmt` are provided [here](https://github.com/neel-lab/GlycoEnzOnto/tree/main/enrichment_analysis).

## Suggesting Additions and Changes to GlycoEnzOnto:

If you would like to contribute to GlycoEnzOnto, such as adding new glycosylation pathway terms, suggesting new glycoenzyme membership to pathways, or any other changes, please submit an Issue on this Github page.  Please describe your suggested changes in the Issue you create, and we will contact you to discuss these new additions.
