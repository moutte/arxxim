#!/usr/bin/python 

import glob,os,sys
import pylab as plt

os.chdir("../")
#print "DIR=", os.getcwd()

#---------------------------------------------------------------EXE FILE
if sys.platform.startswith("win"):
#windows
  sExe= "arxim.exe"
  sExe= os.path.join("..","bin",sExe)

if sys.platform.startswith("linux"):
#linux
  sExe= "arx-basis"
  sExe= os.path.join("..","..","arx-basis","bin",sExe)

  sExe= "arx_debug"
  sExe= "a.out"
  sExe= os.path.join("..","bin",sExe)

Debug= "3"
#---------------------------------------------------------------------//

#----------------------------------------------------------INPUT FILE(S)
files= glob.glob("tmp/test*.inn")
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
    
  Xlabel= "xAb"
  Ylabel= "G_mix"
  figName= "test_mixmodel"
  
  with open("test_mixmodel.tab",'r') as f: lines= f.readlines()
  #--------------------------------------------------------plot XY diagram
  plt.rcParams['figure.figsize']= 8.,6.   #ratio 4/3
  plt.rcParams['figure.figsize']= 8.,8.   #ratio 4/4
  plt.rcParams.update({'font.size': 9})

  fig= plt.subplot(1,1,1)
  symbols=['bo','go','ro','cs','mD','yd']
  symbols=['b','g','r','c','m','y']
  fig.grid(color='r', linestyle='-', linewidth=0.2)
  fig.grid(True)
  
  #fig.set_xlim(Xmin,Xmax)
  #fig.set_ylim(Ymin,Ymax)
  fig.set_xlim(0.,1.)

  for i,line in enumerate(lines):
    if i==0: continue #skip first line
    vx= []
    vy= []
    ww= line.split()
    x= 0.
    for w in ww:
      x= x+0.01
      y= float(w)
      vx.append(x)
      vy.append(y)
    fig.plot(vx, vy, symbols[i%len(symbols)], linestyle='-', linewidth=2.0)
  
  plt.xlabel(Xlabel)
  plt.ylabel(Ylabel)
    
  plt.savefig(figName+".png")
  plt.show()
  #------------------------------------------------------//plot XY diagram
  
  if(i==i0):
    raw_input()
    i= 0
