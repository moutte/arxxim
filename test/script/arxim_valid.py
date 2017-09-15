#!/usr/bin/python 
import glob
import os
import sys

os.chdir("../valid/")
#print "DIR=", os.getcwd()

#---------------------------------------------------------------EXE FILE
if sys.platform.startswith("win"): sExe= "arxim.exe"   #windows
if sys.platform.startswith("linux"): sExe= "a.out" #linux
sDir= os.path.join("..","..","bin")
sExe= os.path.join(sDir,sExe)

if 0:
  sExe= "../../../arx-basis/bin/a.out"
  sExe= "../../../mybin_debug/arx-ifp"

Debug= "2"
#---------------------------------------------------------------------//

#----------------------------------------------------------INPUT FILE(S)
files= glob.glob("inn/a2*.inn")
files.sort()
for f in files: print f
#raw_input()
#---------------------------------------------------------------------//

i0= 0

i= 0
for sFile in files:
  print '\n=========processing '+ sFile + '==========================\n'
  
  if Debug=="1":
    for f in glob.glob("*.*"):
      if os.path.isfile(f): os.remove(f)

  os.system("%s %s %s" % (sExe,sFile,Debug))
  
  print '\n=========done '+ sFile + '==========================\n\n'
  i += 1
  
  if(i==i0):
    raw_input()
    i= 0
