import os
import sys
# import math as m
# import numpy as np
# import MyLib as ML
import matplotlib.pyplot as plt

os.chdir("../")
# print "DIR=", os.getcwd()

fName= "minerals.restab"
fName= "species_dominant.restab"

with open(fName,'r') as  f:
  lines = f.readlines()

Data=  []
Names= []
for i,line in enumerate(lines):
  if line.strip()=="": 
    continue
  DataX= []
  ww= line.split()
  for j,w in enumerate(ww):
    if j>0:
      if not w in Names:
        Names.append(w)
      DataX.append(float(Names.index(w)))
  if i>1:
    Data.append(DataX)

for i,N in enumerate(Names): print i,N

# sys.exit()

#plt.figure(1)
#plt.pcolor(Data)

plt.figure()
# plt.imshow(Data, interpolation='none')
plt.imshow(Data, aspect='auto', interpolation='none', origin='lower')
           # extent=extents(x) + extents(y), 
plt.savefig("path_XY.png")

# Extent defines the images max and min of the horizontal
# and vertical values.
# It takes four values like so: 
# extent=[horizontal_min,horizontal_max,vertical_min,vertical_max]

plt.show()

sys.exit()
