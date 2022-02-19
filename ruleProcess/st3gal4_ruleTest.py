from glycogeneObjs import *

#ST6GAL1
st6gal1=ggenes['ST6GAL1']
#ST6GAL1 substrate:
substrate=['Gal(b1-4)GlcNAc(b1-2)[Gal(b1-4)GlcNAc(b1-4)]Man(a1-3)[Gal(b1-4)GlcNAc(b1-2)[Gal(b1-4)GlcNAc(b1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-1)']

#Predict Products:
st6gal1.forward(substrate[0])

