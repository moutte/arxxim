#!/usr/bin/python 

import glob,os,sys
import MyLib as ML
import pylab as plt

VALID= False

os.chdir("../")
if VALID: os.chdir("valid/")
#print "DIR=", os.getcwd()

#---------------------------------------------------------------EXE FILE
if sys.platform.startswith("win"):    sExe= "arxim.exe"  #windows
if sys.platform.startswith("linux"):  sExe= "a.out"      #linux
if VALID: sDir= os.path.join("..","..","bin")
else:     sDir= os.path.join("..","bin")
sExe= os.path.join(sDir,sExe)

Debug= "2"
#---------------------------------------------------------------------//

#----------------------------------------------------------INPUT FILE(S)
files= glob.glob("inn/map_dom_fe_ox.inn")
files.sort()
for f in files: print f
#raw_input()
#---------------------------------------------------------------------//

#-------------------------------------------------------------check_done
def check_done():
  Ok= False
  if os.path.isfile("error.log"):
    res= open("error.log",'r').read()
    if res.strip()=="PERFECT": Ok= True
    else: print "error.log="+ll ; raw_input()
  else:
    print "error.log NOT FOUND"
  return Ok
#---------------------------------------------------------------------//

i0= 0
i= 0
for sFile in files:
  if os.path.isfile("error.log"): os.remove("error.log")

  print '\n=========processing '+ sFile + '==========================\n'
  os.system("%s %s %s" % (sExe,sFile,Debug))
  
  if not check_done(): continue
  
  sCommand= ML.inn_scan_list(sFile,"TEST","COMPUTE")
  sDirout=  ML.inn_scan_word(sFile,"CONDITIONS","OUTPUT").replace('\\','/')
  
  if "SPCPATH" in sCommand:
    #-------------------------------------------------------------------
    s=sDirout+                                           "_molal.restab"
    lines= open(s,'r').readlines()
    labels,tData= ML.lines2table(lines)
    #
    if "pH" in labels:
      iX= labels.index("pH")
      labX= "pH"
    else:
      iX= 0
      labX= "x"
    dataX= (labX,tData[:,iX])
    #
    dataY= ML.table_select(1,"H2O",labels,tData)
    #
    fig= plt.subplot()
    ML.plot(fig,dataX,dataY,False,True)
    plt.show()
    #-----------------------------------------------------------------//
    
    #-------------------------------------------------------------------
    s=sDirout+                                           "_activ.restab"
    lines= open(s,'r').readlines()
    labels,tData= ML.lines2table(lines)
    #
    if "pH" in labels:
      iX= labels.index("pH")
      labX= "pH"
    else:
      iX= 0
      labX= "x"
    dataX= (labX,tData[:,iX])
    #
    dataY= ML.table_select(1,"H2O",labels,tData)
    #
    fig= plt.subplot()
    ML.plot(fig,dataX,dataY,False,False)
    plt.show()
    #-----------------------------------------------------------------//
    
  if "DYN" in sCommand:
    #-------------------------------------------------------------------
    if 0:
      s=sDirout+"_activ.restab"
      lines= open(s,'r').readlines()
      labels,tData= lines2table(lines)
      #
      if "TIME/YEAR" in labels:
        iX=   labels.index("TIME/YEAR")
        labX= labels[iX]
      else:
        iX= 0
        labX= "x"
      dataX= (labX,tData[:,iX])
      dataY= ML.table_select(1,"H2O",labels,tData)
      #
      fig= plt.subplot()
      ML.plot(fig,dataX,dataY,True,False)
      plt.show()
    #-----------------------------------------------------------------//

    #-------------------------------------------------------------------
    s=sDirout+"_minmol.restab"
    lines= open(s,'r').readlines()
    labels,tData= ML.lines2table(lines)
    #
    if "Time/YEAR" in labels:
      iX=   labels.index("Time/YEAR")
      labX= labels[iX]
    else:
      iX= 0
      labX= "x"
    dataX= (labX,tData[:,iX])
    dataY= ML.table_select(2,"PhiM_",labels,tData)
    #
    fig= plt.subplot()
    ML.plot(fig,dataX,dataY,False,False)
    plt.show()
    fig= plt.subplot()
    ML.plot(fig,dataX,dataY,True,True)
    plt.show()
    #-----------------------------------------------------------------//
  
  print '\n=========done '+ sFile + '==========================\n\n'
  i += 1
  
  if(i==i0):
    raw_input()
    i= 0
