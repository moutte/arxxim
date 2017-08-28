#!/usr/bin/python 

import glob #for listing files, Langtangen, p118
import os

os.chdir("../")
#print "DIR=", os.getcwd()

#----------------------------------------------------------INPUT FILE(S)
files= glob.glob('inn/map1a*.inn')
files.sort()
for f in files: print f
#raw_input()
#---------------------------------------------------------------------//

#---------------------------------------------------------------EXE FILE
sProg= "../bin/arxim"
Debug= "3"
#---------------------------------------------------------------------//

i0= 0

i= 0
for sFile in files:
  print '\n=========processing '+ sFile + '==========================\n'
  os.system("%s %s %s" % (sProg,sFile,Debug))
  
  print '\n=========done '+ sFile + '==========================\n\n'
  i += 1
  
  if(i==i0):
    raw_input()
    i= 0
