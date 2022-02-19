from glycogeneObjs import *
import itertools

#N-linked glycosylation:

s1='GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-1)'
P1='GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-1)'
P2='GlcNAc(b1-2)[GlcNAc(b1-4)]Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-1)'
P3='Gal(b1-4)GlcNAc(b1-2)[GlcNAc(b1-4)]Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-1)'
P4='Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-1)'
P5='Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-1)'
P6='Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-1)'

#Glycan substrate list:
glycan_substrates=[s1,P1,P2,P3,P4,P5,P5,P6]

#Glycogenes:
mgat2=ggenes['MGAT2']
mgat4=ggenes['MGAT4A']
b4galt=ggenes['B4GALT1']
st6gal1=ggenes['ST6GAL1']
st3gal1=ggenes['ST3GAL1']
#Glycogene list:
ggene_list={'MGAT2':mgat2,'MGAT4':mgat4,'B4GALT':b4galt,'ST6GAL1':st6gal1,'ST3GAL1':st3gal1}

#Loop through all substrates and predict connectivity:
for (gg,gg_obj) in ggene_list.items():
    substrate_prod_dict={g:gg_obj.forward(g) for g in glycan_substrates}
    if all([l is None or len(l)==0 for l in substrate_prod_dict.values()]):
        print(f"Glycogene {gg} did not produce any products")
        next
    print(f"Reactions catalyzed by glycogene {gg}")
    for sub,prods in substrate_prod_dict.items():
        if prods is not None:
            if len(prods)>0:
                for p in prods:
                    print(f"Substrate: {sub} --> Product: {p}")
            else:
                continue
        else:
            continue
    print('\n')
