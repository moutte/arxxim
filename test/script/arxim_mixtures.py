import os, sys
import pylab as plt

#-----------------------------------------------test mixture mean comp'n
def mix_is_same(v1,v2,tol):
  #return abs(v1[0]-v2[0])+abs(v1[1]-v2[1])+abs(v1[2]-v2[2])<tol
  #return sum(abs(v1[:]-v2[:]))<tol*len(v1)
  return max(abs(v1[:]-v2[:]))<tol

os.chdir("../")

tMix= []
with open("mixture.tab",'r') as f:
  lines= f.readlines()
for il,line in enumerate(lines):
  ww= line.split()
  if il>0:
    x= ww[0]
    y= ww[1]
    v= plt.zeros(len(ww)-2,'float')
    for i,w in enumerate(ww):
      if i>1: v[i-2]= float(w)
    tMix.append((x,y,v))

f= open("mix_mean.tab",'w')

v0n= plt.zeros(len(ww)-2,'float')
v1n= plt.zeros(len(ww)-2,'float')
v2n= plt.zeros(len(ww)-2,'float')

# 3 poles = 3 paires, 1 triples
# 1-2 1-3
# 2-3
# 4 poles = 6 paires, 3 triples
# 1-2 1-3 1-4
# 2-3 2-4
# 3-4
# 5 poles = 10 paires, 10 triples
# 1-2 1-3 1-4 1-5
# 2-3 2-4 2-5
# 3-4 3-5
# 4-5
# 1-2-3 1-2-4 1-2-5 1-3-4 1-3-5 1-4-5
# 2-3-4 2-3-5 2-4-5
# 3-4-5


tol= 0.05
for x,y,v in tMix:
  nmix= 3
  p0= v[0]  ;  v0= v[1:4]   ;  print p0,v0 #
  p1= v[4]  ;  v1= v[5:8]   ;  print p1,v1 #
  p2= v[8]  ;  v2= v[9:12]  ;  print p2,v2 #
  v0n= v0
  v1n= v1
  v2n= v2
  if mix_is_same(v0,v1,tol): #--------------v0=v1
    nmix= nmix-1
    v0n= (p0*v0+p1*v1)/(p0+p1)
    p0= p0+p1
    p1= 0.
    if mix_is_same(v1,v2,tol): #------------v0=v1=v2
      nmix= nmix-1
      v0n= (p0*v0n+p2*v2)/(p0+p2)
      p0= p0+p2
      p2= 0.
  else:
    if mix_is_same(v1,v2,tol): #------------v1=v2!=v0
      nmix= nmix-1
      v1n= (p1*v1+p2*v2)/(p1+p2)
      p1= p1+p2
      p2= 0.
    else:
      if mix_is_same(v2,v0,tol): #----------v0=v2!=v1
        nmix= nmix-1
        v0n= (p0*v0+p2*v2)/(p0+p2)
        p0= p0+p2
        p2= 0.
  #
  print "n=",nmix
  #if nmix>2: raw_input()
  f.write("%s\t%s\t" % (x,y))
  f.write("%d\t" % nmix)
  f.write("%.4g\t" % p0)
  for x in v0n: f.write("%.4g\t" % x)
  if nmix>1:
    f.write("%.4g\t" % p1)
    for x in v1n: f.write("%.4g\t" % x)
    if nmix>2:
      f.write("%.4g\t" % p2)
      for x in v2: f.write("%.4g\t" % x)
    
  f.write('\n')

f.close()

sys.exit()
#---------------------------------------------//test mixture mean comp'n
