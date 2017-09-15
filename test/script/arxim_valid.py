#!/usr/bin/python 
import glob
import os
import sys

os.chdir("../valid/")
#print "DIR=", os.getcwd()

#---------------------------------------------------------------EXE FILE
if sys.platform.startswith("win"): sExe= "arxim.exe"   #windows
if sys.platform.startswith("linux"): sExe= "a.out" #linux
sExe= os.path.join("..","..","bin",sExe)

sExe= "../../../arx-basis/bin/a.out"
sExe= "../../../mybin_debug/arx-ifp"


Debug= "2"
#---------------------------------------------------------------------//

#----------------------------------------------------------INPUT FILE(S)
files= glob.glob("inn/b1b*.inn")
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
