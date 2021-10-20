import re
from owlready2 import * 
from glygen_glycan_query import * # Gets glygen_glycans list with GlycoCT
import types
# From Frederique Lisacek's RDF Glycan Triplestore paper:
# They have a prefix in their sparql queries that designate predicates (Object Properties)
# I'll probably end up ignoring this, as I think it's more clear to just use 
# the namespace where things are defined rather than create nested prefixes.
# prefix glycan:  <http://mzjava.expasy.org/glycan/> (I will be using "http://purl.jp/12/glycan#")
# prefix predicate:  <http://mzjava.expasy.org/predicate/>
# prefix component:  <http://mzjava.expasy.org/component/>
# prefix structureConnection:  <http://mzjava.expasy.org/structureConnection/>

glycoStructOnto=get_ontology('http://mzjava.expasy.org/glycan/glycoStructOnto.owl')
glycoStructOnto.base_iri='http://mzjava.expasy.org/glycan/'
#gso_namespace=glycoStructOnto.get_namespace('http://mzjava.expasy.org/glycan/')
#This loads the SKOS ontology
skos=get_ontology('http://www.w3.org/2004/02/skos/core').load()
#GlycoRDF:
glycoRDF=get_ontology('https://raw.githubusercontent.com/ReneRanzinger/GlycoRDF/master/ontology/glycan.owl').load()


################################
# Components possibly for later:
################################
#allBaseComponents=set(list(chain(*[[re.sub('^\d+?b\:','',r) for r in g[1].split('LIN')[0].split(' ') if re.search('^\d+?b\:',r) is not None] for g in glygen_glycans])))


class GlycoCTProcessor:

    def __init__(self,accNum,glycoCT):
        '''
        This class reads in a raw glycoCT string and splits it up
        into its residue and linkage components
        '''
        self.accNum=accNum
        self.glycoCT=glycoCT
        self.residues=re.sub('RES','',self.glycoCT.split('LIN')[0]).split(' ')
        self.residues=[x for x in self.residues if x!='']
        self.residues=[''.join(['RES ',x]) for x in self.residues]
        self.linkages=self.glycoCT.split('LIN')[1].split(' ')
        self.linkages=[x for x in self.linkages if x!='']

    def match_residue(self,res):
        '''
        Matches Components of Residue strings:
        '''
        front_pat=re.compile(r'(?P<ResidueNumber>\d+?)(?P<ResidueType>[bsox])\:')
        front_match=re.search(front_pat,res)
        if front_match is not None:
            #What kind of residue is it?
            if front_match.groupdict()['ResidueType']=='b':
                #Is a base type:
                pat=re.compile(r'(?P<ResidueNumber>\d+?)(?P<ResidueType>[bsox])\:(?P<Anomericity>[abxo])\-(?P<CarbClasses>\D+)\-(?P<RingStartCarbon>\d)\:(?P<RingEndCarbon>\d)(?P<Modifications>(?:\|\d+?\:(?:d|a|keto|sp|sp2|keto|geminal|adli|en|enx))+)?')
            elif front_match.groupdict()['ResidueType']=='s':
                #Is a substituent:
                pat=re.compile(r'(?P<ResidueNumber>\d+?)(?P<ResidueType>[bsox])\:(?P<SubstituentName>.+)')
            mtch=re.search(pat,res).groupdict()
            for k,v in mtch.items():
                if v is None:
                    mtch[k]=''
            return((mtch['ResidueNumber'],mtch))
        else:
            raise Exception('Residue doesnt parse correctly')
    
    def match_link(self,link):
        pat=re.compile(r'(?P<LinkNumber>\d+?)\:(?P<ParentIndex>\d+?)(?P<ParentModType>[odhnxrs])\((?P<ParentCarbonLinks>(?:\-?\d+?(?:\|\d+?)*))\+(?P<ChildCarbonLink>\-?\d+?)\)(?P<ChildIndex>\d+?)(?P<ChildModType>[odhnxrs])')
        mtch=re.search(pat,link).groupdict()
        for k,v in mtch.items():
            if v is None:
                mtch[k]=''
        return((mtch['LinkNumber'],mtch))
    
    def make_mono_residue(self,baseResidue,onto):
        onto_monos=onto.get_children_of(onto.monosaccharide)
        onto_monos_labels=[x.label[0] for x in onto_monos]
        onto_monos_names=[x.name for x in onto_monos]
        onto_monos={ct:v for ct,v in zip(onto_monos_labels,onto_monos_names)}
        try:
            with onto:
                mono_res=onto[onto_monos[baseResidue]]()
        except:
            raise Exception(f'passed residue {baseResidue} not found in ontology monosaccharides')
        return(mono_res)

    def make_substituent_residue(self,substituent,onto):
        onto_subs=onto.get_children_of(onto.substituent)
        onto_subs_labels=[x.label[0] for x in onto_subs]
        onto_subs_names=[x.name for x in onto_subs]
        onto_subs={sb:v for sb,v in zip(onto_subs_labels,onto_subs_names)}
        try:
            with onto:
                subst_res=onto[onto_subs[substituent]]()
        except:
            raise Exception(f'passed residue {substituent} not found in ontology substituents')
        return(subst_res)

    def residue_creator(self,residueData,onto):
        #For Residue Objects:
        residueObjects=dict()
        with onto:
            #Residues:
            for k,v in residueData.items():
                if v['ResidueType']=='b':
                    #Process as a base residue type:
                    #What is the base residue?:
                    baseResidue=''.join(['RES 1b:',v['Anomericity'],'-',v['CarbClasses'],'-',v['RingStartCarbon'],':',v['RingEndCarbon'],v['Modifications']])
                    mono_res_class=self.make_mono_residue(baseResidue,onto)
                    residueObjects[int(v['ResidueNumber'])]=mono_res_class
                elif v['ResidueType']=='s':
                    #Process as a substituent:
                    substResidue=''.join(['RES 1s:',v['SubstituentName']])
                    subst_res_class=self.make_substituent_residue(substResidue,onto)
                    residueObjects[int(v['ResidueNumber'])]=subst_res_class
        return(residueObjects) 

    def mono_mono_connect(self,fromRes,toRes,anomerLink,anomericity,linkageNumber):
        fromRes.is_GlycosidicLinkage.append(toRes)
        #Handle anomer connection:
        if anomerLink==1:
            fromRes.has_anomerCarbon_1.append(toRes)
        elif anomerLink==2:
            fromRes.has_anomerCarbon_2.append(toRes)

        if linkageNumber==1:
            fromRes.has_linkedCarbon_1.append(toRes)
        elif linkageNumber==2:
            fromRes.has_linkedCarbon_2.append(toRes)
        elif linkageNumber==3:
            fromRes.has_linkedCarbon_3.append(toRes)
        elif linkageNumber==4:
            fromRes.has_linkedCarbon_4.append(toRes)
        elif linkageNumber==5:
            fromRes.has_linkedCarbon_5.append(toRes)
        elif linkageNumber==6:
            fromRes.has_linkedCarbon_6.append(toRes)
        elif linkageNumber==8:
            fromRes.has_linkedCarbon_8.append(toRes)

        if anomericity=='a':
            fromRes.has_anomericConnection_alpha.append(toRes)
        elif anomericity=='b':
            fromRes.has_anomericConnection_beta.append(toRes)

    def mono_subst_connect(self,fromRes,toRes,linkageNumber):
        fromRes.is_SubstituentLinkage.append(toRes)
        if linkageNumber==1:
            fromRes.has_linkedCarbon_1.append(toRes)
        elif linkageNumber==2:
            fromRes.has_linkedCarbon_2.append(toRes)
        elif linkageNumber==3:
            fromRes.has_linkedCarbon_3.append(toRes)
        elif linkageNumber==4:
            fromRes.has_linkedCarbon_4.append(toRes)
        elif linkageNumber==5:
            fromRes.has_linkedCarbon_5.append(toRes)
        elif linkageNumber==6:
            fromRes.has_linkedCarbon_6.append(toRes)
        elif linkageNumber==8:
            fromRes.has_linkedCarbon_8.append(toRes)
    
    def glyco_object(self,onto):
        '''
        Processes residue and linkage information: 
        '''
        with onto:
            #Residue String Data:
            residueData=dict([self.match_residue(res) for res in self.residues])
            #Linkage String Data:
            linkData=dict([self.match_link(link) for link in self.linkages])
            #Process the residue data:
            residueObjects=self.residue_creator(residueData,onto)
            #Link residues with a glycan:
            glyc=onto.glycan()
            #Set the glycan iri to be defined by GlyGen:
            glyc.iri=self.accNum
            for r,e in residueObjects.items():
                glyc.has_residue.append(e)
            #Make links between residues in Glycan:
            for l in linkData.values():
                #From/To Information:
                fromRes=residueObjects[int(l['ParentIndex'])]
                toRes=residueObjects[int(l['ChildIndex'])]
                #Get the anomer carbon number
                anomerLink=int(l['ChildCarbonLink'])
                #Get the linkage carbon number:
                linkageNumber=[int(x) for x in l['ParentCarbonLinks'].split('|')][0]
                #Ontology Stuff:
                if onto.monosaccharide in toRes.is_instance_of[0].ancestors():
                    #Get the anomericity of the child monosaccharide and
                    # pass it to the mono_mono_connect function:
                    anomericity=residueData[l['ChildIndex']]['Anomericity']
                    #Mono-Mono connection:
                    self.mono_mono_connect(fromRes,toRes,anomerLink,anomericity,linkageNumber)
                elif onto.substituent in toRes.is_instance_of[0].ancestors():
                    #Mono-substituent connection:
                    self.mono_subst_connect(fromRes,toRes,linkageNumber)

##########################
# Residue Types in Glygen:
##########################
#Remove any glycans with repeat groups in their GlycoCT:
glygen_glycans=[x for x in glygen_glycans if 'REP' not in x[1]]
gob_dict={x[0]:GlycoCTProcessor(x[0],x[1]) for x in glygen_glycans}
unique_residues=list(set(list(chain(*[[re.sub('\d+?([bs])\:(.+)','1\g<1>:\g<2>',x) for x in oj.residues] for oj in gob_dict.values()]))))
unique_monosaccharide_residues=[x for x in unique_residues if re.search('^RES\ 1b',x) is not None]
unique_substituent_residues=[x for x in unique_residues if re.search('^RES\ 1s',x) is not None]

#Tally:
unique_monosaccharide_tally={x:sum([re.sub('RES\ 1b\:','',x) in g[1] for g in glygen_glycans]) for x in unique_monosaccharide_residues}
lowFreqMonos=[k for k,v in unique_monosaccharide_tally.items() if v<=1]
highFreqMonos=[k for k,v in unique_monosaccharide_tally.items() if v>1]

########################################
# Glypy's Unique Monosaccharide Residues
########################################
monoRes=dict()
for (m,v) in glypy.monosaccharides.items():
    mr=v.serialize()
    if re.search('LIN',mr) is None and mr not in monoRes.values():
        monoRes[''.join(['monores',str(len(monoRes)+1)])]=mr
    elif re.search('LIN',mr) is not None:
        mr=re.split('\ 2s\:',mr)[0]
        if mr not in monoRes.values():
            monoRes[''.join(['monores',str(len(monoRes)+1)])]=mr

###############################################################
# Monosaccharide residues in Glygen & NOT represented in Glypy:
###############################################################
#Find the non-represented monosaccharide residues:
noRep=[x for x in unique_monosaccharide_residues if x not in monoRes]
# If these residues happen in at least 10 monosaccharides or more,
# add them to the monoRep list:
for nr in noRep:
    if unique_monosaccharide_tally[nr]>=10:
        if nr not in monoRes.values():
            monoRes[''.join(['monores',str(len(monoRes)+1)])]=nr

#Filter out glycans that do not have residues in the monoRes
# list:
new_gob_dict=dict()
for k,g in gob_dict.items():
    mono_residues=[x for x in g.residues if re.search('\d+?b\:',x) is not None]
    mono_residues=[re.sub('RES\ \d+?b','RES 1b',x) for x in mono_residues]
    if all([m in monoRes.values() for m in mono_residues]):
        new_gob_dict[k]=g

print('Instantiating %i glycans into structure ontology' %(len(new_gob_dict)))

#############################
# glycoStructOnto Definition:
#############################

with glycoStructOnto:
    ##########
    # Classes:
    ##########
    class glycan(Thing):
        label=['Glycan']
    class residue(Thing):
        label=['Residue']
    class monosaccharide(residue):
        label=['Monosaccharide']
    #Dynamically create subclasses of "monosaccharide"
    # using monosaccharides defined in glypy:
    for k,v in monoRes.items():
        nc=types.new_class(k,(monosaccharide,))
        nc.label=v

    class substituent(residue):
        label=['Substituent']
    #Gather substituents in glygen's set of unambiguously
    # defined glycans:
    substituents=dict()
    for g in glygen_glycans:
        residues=g[1].split('LIN')[0].split(' ')
        for r in residues:
                if re.search('^\d+?s',r) is not None:
                    rs=re.sub('^(\d+?)s\:','RES 1s:',r)
                    if rs not in substituents.values():
                        substituents[''.join(['subsres',str(len(substituents)+1)])]=rs
    for s,v in substituents.items():
        subst_in=types.new_class(s,(substituent,))
        subst_in.label=v
    ####################
    # Object Properties:
    ####################
    #Associates residues with a glycan definition:
    class has_residue(ObjectProperty):
        domain=[glycan]
        range=[residue]
    ### Linkages: 
    class links_To(ObjectProperty):
        pass
    class is_SubstituentLinkage(links_To):
        domain=[monosaccharide]
        range=[substituent]
    class is_GlycosidicLinkage(links_To):
        domain=[monosaccharide]
        range=[monosaccharide]
    ### Anomeric carbon linkage positions:
    class has_anomerCarbon_1(ObjectProperty):
        domain=[monosaccharide]
        range=[monosaccharide]
    class has_anomerCarbon_2(ObjectProperty):
        domain=[monosaccharide]
        range=[monosaccharide]
    ### Linkage positions:
    class has_linkedCarbon_1(ObjectProperty):
        domain=[monosaccharide]
        range=[monosaccharide]
    class has_linkedCarbon_2(ObjectProperty):
        domain=[monosaccharide]
        range=[monosaccharide]
    class has_linkedCarbon_3(ObjectProperty):
        domain=[monosaccharide]
        range=[monosaccharide]
    class has_linkedCarbon_4(ObjectProperty):
        domain=[monosaccharide]
        range=[monosaccharide]
    class has_linkedCarbon_5(ObjectProperty):
        domain=[monosaccharide]
        range=[monosaccharide]
    class has_linkedCarbon_6(ObjectProperty):
        domain=[monosaccharide]
        range=[monosaccharide]
    class has_linkedCarbon_8(ObjectProperty):
        domain=[monosaccharide]
        range=[monosaccharide]
    class has_anomericConnection_alpha(ObjectProperty):
        domain=[monosaccharide]
        range=[monosaccharide]
    class has_anomericConnection_beta(ObjectProperty):
        domain=[monosaccharide]
        range=[monosaccharide]

if __name__=='__main__':
    #Ontology-defined monosaccharides:
    onto_monos=[x.name for x in glycoStructOnto.get_children_of(glycoStructOnto['monosaccharide'])]
    onto_monos={x:y for x,y in glypy.monosaccharides.items() if x in onto_monos}
    #Ontology-defined substituents:
    onto_subst=[x.name for x in glycoStructOnto.get_children_of(glycoStructOnto['substituent'])]

    #Testing:
    couldntDos=0
    for g,v in new_gob_dict.items():
        try:
            v.glyco_object(glycoStructOnto)
        except Exception as e: print(f'Couldn\'t do {g}');print(e);couldntDos+=1

    print("Total unprocessed: {0}".format(couldntDos))
    glycoStructOnto.save(file='glycoStructOnto.rdf',format='rdfxml')
    print('Saved \"glycoStructOnto.rdf\"')
