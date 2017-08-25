#!/usr/bin/python 

import glob #for listing files, Langtangen, p118
import os

os.chdir("../")
print "DIR=", os.getcwd()

filelist= glob.glob('inn/map2b*.inn')

filelist.sort()

for i in range(len(filelist)):
  print filelist[i]

sProg= "../bin/arxim"

Debug= "3"
i0= 0

i= 0
for sFile in filelist:
  print '\n=========processing '+ sFile + '==========================\n'
  os.system("%s %s %s" % (sProg,sFile,Debug))
  
  print '\n=========done '+ sFile + '==========================\n\n'
  i += 1
  
  if(i==i0):
    raw_input()
    i= 0
