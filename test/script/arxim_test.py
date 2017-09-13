#!/usr/bin/python 

import glob #for listing files, Langtangen, p118
import os
import sys

s=sys.platform
print s
if sys.platform.startswith("linux"):
  print "linux"
sys.exit()


os.chdir("../")
#print "DIR=", os.getcwd()

#----------------------------------------------------------INPUT FILE(S)
files= glob.glob("tmp/f1e*.inn")
files= glob.glob("inn/map3d_tp.inn")
files= glob.glob("tmp/map2a*.inn")
files.sort()
for f in files: print f
#raw_input()
#---------------------------------------------------------------------//

#---------------------------------------------------------------EXE FILE
if sys.platform.startswith("win"):
#windows
  sExe= "arxim.exe"
  sExe= os.path.join("..","bin",sExe)

if sys.platform.startswith("linux"):
#linux
  sExe= "arx-basis"
  sExe= os.path.join("..","..","arx-basis","bin",sExe)

  sExe= "a.out"
  sExe= os.path.join("..","bin",sExe)

Debug= "3"
#---------------------------------------------------------------------//

i0= 0

i= 0
for sFile in files:
  print '\n=========processing '+ sFile + '==========================\n'
  os.system("%s %s %s" % (sExe,sFile,Debug))
  
  print '\n=========done '+ sFile + '==========================\n\n'
  i += 1
  
  if(i==i0):
    raw_input()
    i= 0
