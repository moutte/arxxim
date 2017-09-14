#!/usr/bin/python 
import glob
import os
import sys

os.chdir("../valid/")
#print "DIR=", os.getcwd()

#---------------------------------------------------------------EXE FILE
#windows
if sys.platform.startswith("win"): sExe= "arxim.exe"
#linux
if sys.platform.startswith("linux"): sExe= "arx_optim"

sExe= os.path.join("..","..","bin",sExe)

Debug= "3"
#---------------------------------------------------------------------//

#----------------------------------------------------------INPUT FILE(S)
files= glob.glob("inn/d2*.inn")
files.sort()
for f in files: print f
#raw_input()
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
