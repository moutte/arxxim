#!/usr/bin/python 

import glob #for listing files, Langtangen, p118
import os

os.chdir("../")
#print "DIR=", os.getcwd()

#----------------------------------------------------------INPUT FILE(S)
files= glob.glob('inn/map3b*.inn')
files.sort()
for f in files: print f
#raw_input()
#---------------------------------------------------------------------//

#---------------------------------------------------------------EXE FILE
sExe= "arxim.exe"    #windows
sExe= os.path.join("..","bin",sExe)

sExe= "arx-cell"
sExe= "arx-win"
sExe= os.path.join("..","..","mybin_debug",sExe)

sExe= "arxxim"        #linux
sExe= "arx_bis"
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
