import pandas as pd
from SPARQLWrapper import SPARQLWrapper, JSON, CSV, TSV
import requests
import json
import urllib
import numpy as np
#The database we are querying is GlyGen:
db="http://sparql.glygen.org:8880/sparql"

sparql = SPARQLWrapper(db)
#GlycoMotif Query:
# Retrieves glycan motifs from glycoMotif database
# curated by the folks at GlyGen.

query1="""
    PREFIX glycan: <http://purl.jp/bio/12/glyco/glycan#>
    PREFIX gly: <https://sparql.glygen.org/ontology/>
    SELECT DISTINCT ?glycan,?glycanSeq
    WHERE {
        ?glycan glycan:has_glycosequence ?seqAcc .
        ?glycan glycan:is_from_source ?source .
        ?source glycan:has_taxon <http://purl.uniprot.org/core/taxonomy/9606> .
        ?seqAcc glycan:has_sequence ?glycanSeq .
        ?seqAcc glycan:in_carbohydrate_format <http://purl.jp/bio/12/glyco/glycan#carbohydrate_format_wurcs> .
    }
"""
query2=''.join([query1,'LIMIT 10000 OFFSET 10000'])

sparql.setQuery(query1)

def wurcs2iupacCondensed(wurcs):
	'''
	Calls to a glycan structure converter API from GlyCosmos.
    input: url-escaped WURCS string
    output: IUPAC condensed string
	'''
	baseURL="https://api.glycosmos.org/glycanformatconverter/2.4.1/wurcs2iupaccondensed/"
	wurcs=urllib.parse.quote(wurcs,safe='')
	totalURL=baseURL+wurcs
	res=json.loads(requests.get(totalURL).content)
	return(res['IUPACcondensed'])

sparql.setReturnFormat(JSON)

if __name__=='__main__':
    #results = sparql.query().convert().decode('utf-8').replace('\"','')
    #Returns a list of objects to parse:
    results=sparql.query().convert()['results']['bindings']
    if len(results)==10000:
        sparql.setQuery(query2)
        results=results+(sparql.query().convert()['results']['bindings'])
    #Read data as pandas dataframe
    df=pd.DataFrame([[x[k]['value'] for k in x.keys()] for x in results],columns=['glycan','WURCS'])
    iupac_strings=[]
    for i,w in enumerate(df['WURCS'].to_list()):
        try:
            newString=wurcs2iupacCondensed(w)
        except:
            iupac_strings.append(np.nan)
            continue
        iupac_strings.append(newString)
    df['IUPAC']=iupac_strings;df.dropna(axis=1,how='any')
    df.to_csv('glygen_glycans_iupac.csv')
