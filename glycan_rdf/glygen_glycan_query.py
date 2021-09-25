import pandas as pd
from SPARQLWrapper import SPARQLWrapper, JSON, CSV, TSV
import requests
import json
import urllib
import numpy as np
import glypy
#The database we are querying is GlyGen:
db="http://sparql.glygen.org:8880/sparql"

sparql = SPARQLWrapper(db)
query="""
    PREFIX glycan: <http://purl.jp/bio/12/glyco/glycan#>
    PREFIX gly: <https://sparql.glygen.org/ontology/>
    SELECT DISTINCT ?glycan,?glycanSeq
    WHERE {
        ?glycan glycan:has_glycosequence ?seqAcc .
        ?glycan glycan:is_from_source ?source .
        ?source glycan:has_taxon <http://purl.uniprot.org/core/taxonomy/9606> .
        ?seqAcc glycan:has_sequence ?glycanSeq .
        ?seqAcc glycan:in_carbohydrate_format <http://purl.jp/bio/12/glyco/glycan#carbohydrate_format_glycoct> .
        BIND(STRBEFORE(?glycanSeq,"LIN") as ?Residues)
        BIND(STRAFTER(?glycanSeq,"LIN") as ?Linkages)
        FILTER( !REGEX(?Linkages,"\\\\d\\\\|\\\\d") && REGEX(?glycanSeq,"RES") && REGEX(?glycanSeq,"LIN") && !REGEX(?Residues,"\\\\-x|x\\\\-|\\\\:x|x\\\\:|x\\\\||\\\\|x") )
    }
"""
sparql.setQuery(query)

sparql.setReturnFormat(JSON)

#Returns a list of objects to parse:
results=sparql.query().convert()['results']['bindings']
#Paired tuples of glycan data:
glygen_glycans=[(x['glycan']['value'],x['glycanSeq']['value']) for x in results]
