#!/usr/bin/python 
import glob
import os
import sys

os.chdir("../valid/")
#print "DIR=", os.getcwd()

#--scan the inn file
def innfile_scan_list(sName,sBlock,sKey):
  list_=[]
  with open(sName,'r') as f: lines= f.readlines()
  Ok= False
  for ll in lines:
    ww=ll.split()
    if len(ww)<1: continue
    if ww[0].upper()==sBlock: Ok=True
    if Ok:
      if ww[0].upper()==sKey: list_.append(ww[1].upper())
      if ww[0].upper()=="END": Ok=False
  return list_
  
def innfile_scan_word(sName,sBlock,sKey):
  word=""
  with open(sName,'r') as f: lines= f.readlines()
  Ok= False
  for ll in lines:
    ww=ll.split()
    if len(ww)<1: continue
    if ww[0].upper()==sBlock: Ok=True
    if Ok:
      if ww[0].upper()==sKey: word= ww[1]
      if ww[0].upper()=="END": Ok=False
  return word
#--//
#---------------------------------------------------------------EXE FILE
if sys.platform.startswith("win"):   #windows
  sExe= "arxim.exe"
if sys.platform.startswith("linux"): #linux
  sExe= "a.out"
sDir= os.path.join("..","..","bin")
sExe= os.path.join(sDir,sExe)

if 0:
  sExe= "../../../arx-basis/bin/a.out"
  sExe= "../../../mybin_debug/arx-ifp"

sDebug= "2"
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
  
  sCommand= innfile_scan_list(sFile,"TEST","COMPUTE")
  sDirout=  innfile_scan_word(sFile,"CONDITIONS","OUTPUT")
  sDirout= sDirout.replace('\\','/')
  print sCommand, sDirout
  raw_input()
  
  if sDebug=="1":
    for f in glob.glob("*.*"):
      if os.path.isfile(f): os.remove(f)

  os.system("%s %s %s" % (sExe,sFile,sDebug))
  
  if "SPCPATH" in sCommand:
    s=sDirout+"_molal.restab"
    with open(s,'r') as f: lines= f.readlines()
  
  print '\n=========done '+ sFile + '==========================\n\n'
  i += 1
  
  if(i==i0):
    raw_input()
    i= 0
