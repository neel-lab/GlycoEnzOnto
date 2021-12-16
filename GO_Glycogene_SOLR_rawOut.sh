#!/bin/bash

#SOLR Query Parameters:
FL_STRING="bioentity_label%2Cannotation_class%2Cannotation_class_label%2Cdescription%2Cevidence_type%2Contology_class"

if [ -e ggene_json_dta_new.json ]
then
	rm ggene_json_dta_new.json
fi

glycogenes=./finishedGlycogenes.tsv
for GGENE in $(awk 'BEGIN{FS="\t"}FNR==1{next}{print $4}' $glycogenes)
do
	echo "Parsing ${GGENE}"
	GOLR_QUERY=$(echo "http://golr-aux.geneontology.io/solr/select?qt=standard&?fq=document_category:%22annotation%22&q=*:*&fq=bioentity_label:%22${GGENE}%22&fq:taxon=%22NCBITaxon%3A9606%22&&&wt=json&rows=10000&start=0&indent=on&fl=${FL_STRING}")
	echo ${GOLR_QUERY}
	#Store GOlr glycogene query in JSON string
	#curl -X GET -s $(echo "${GOLR_QUERY}") | jq .response.docs | jq '. = (. | unique)' | jq .[] >> ggene_json_dta_new.json
	curl -X GET -s $(echo "${GOLR_QUERY}") | jq 'def IN(s): first((s == .) // empty) // false; ["EXP","IDA","IPI","IMP","IGI","IEP"] as $codes | .response.docs | map(select(.evidence_type | IN($codes[]))) | .[]' >> ggene_json_dta_new.json
done

jq -s 'group_by(.bioentity_label) | map({geneName:.[0].bioentity_label,GO_accs:map(.annotation_class)})' ggene_json_dta_new.json > ggene_json_aggData.json
