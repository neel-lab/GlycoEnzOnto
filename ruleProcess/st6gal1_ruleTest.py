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
#Glycogene list:
ggene_list={'MGAT2':mgat2,'MGAT4':mgat4,'B4GALT1':b4galt1,'B4GALT2':b4galt2,'ST6GAL1':st6gal1}

#Loop through all substrates and predict connectivity:
for (g,(gg,gg_obj)) in itertools.product(glycan_substrates,ggene_list.items()):
    prods=gg_obj.forward(g)
    if len(prods)>0:
        print(f"Glycogene {gg} catalyzed the following reactions:")
        for p in prods:
            print(f"Substrate: {g}, Product {p}")
