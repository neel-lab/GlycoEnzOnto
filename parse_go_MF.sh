#!/bin/bash

source /home/ted/Documents/Data/glycoNetworkOnto/glycoNetworkOnto_venv/bin/activate

jq '.[]["GO_accs"][]' ggene_json_aggData.json | sed -e 's/\"//g' | sort | uniq > mf_low_terms.txt


#High level ontology : molecular function:
MF=GO:0003674

#Call ROBOT:

robot extract --prefix 'GO: http://purl.obolibrary.org/obo/GO_' \
	--method MIREOT \
	--input-iri http://purl.obolibrary.org/obo/go.owl \
	--lower-terms mf_low_terms.txt  \
	--upper-term ${MF} \
	--copy-ontology-annotations true \
	--annotate-with-source true \
	remove --term GO:0008150 --term GO:0005575 --select "self descendants" --signature true \
	--output ./GO_glycogene_MF.owl -v 
#Cleanup:
rm mf_low_terms.txt
