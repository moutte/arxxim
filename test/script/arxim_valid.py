#!/usr/bin/python 
import glob,os,sys
import pylab as plt
import MyLib as ML

os.chdir("../valid/")
#print "DIR=", os.getcwd()

#--scan the inn file -> list of words
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
#--//  
#--scan the inn file -> word
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

def lines2table(lines):
  lines= [x for x in lines if len(x.strip())>0] #remove empty lines
  labels= lines[0].split()
  nC= len(labels)
  nL= len(lines)-1
  TT= plt.zeros((nL,nC),'float')
  for i,line in enumerate(lines):
    if i==0: continue #skip first line
    ww= line.split()
    for j in range(nC): TT[i-1,j]= ML.num(ww[j])
  return labels,TT
  
def tableSort(labels,TT):
  nlin,ncol= TT.shape
  T= []
  i0= labels.index("H2O")
  for i in range(ncol):
    if i>i0:
      t= (labels[i],TT[:,i],max(TT[:,i]))
      T.append(t)
  return T

#-------------------------------------------------------------------plot
def plot(labels,data,iX,iY0,xLog=False,yLog=False):
  nlin,ncol= data.shape
  
  fig= plt.subplot()
  fig.grid(color='r', linestyle='-', linewidth=0.2)
  fig.grid(True)
  
  fig.set_xlim(1.e1,1.e7)
  fig.set_ylim(1.e-4,1.0)
  
  symbols=['bo','go','ro','cs','mD','yd']
  colors = ['cyan', 'lightblue', 'lightgreen', 'tan', 'pink','red', 'blue']

  for i in range(ncol):
    if "PhiM_" in labels[i]:
    #if i>iY0:
      sy= symbols[(i-i0-1)%len(symbols)]
      lb= labels[i].replace("PhiM_","")
      vx= data[:,iX]
      vy= data[:,i]
      if xLog and yLog:
        fig.loglog(vx, vy, sy, linestyle='-',linewidth=1.0,label=lb)
      else:
        if xLog:
          fig.semilogx(vx, vy, sy, linestyle='-',linewidth=1.0,label=lb)
        elif yLog:
          fig.semilogy(vx, vy, sy, linestyle='-', linewidth=1.0,label=lb)
        else:
          fig.plot(vx, vy, sy, linestyle='-', linewidth=1.0,label=lb)
      
  #-legend= fig.legend(loc='upper left')
  #-legend= fig.legend(bbox_to_anchor=(1.1, 1.05))
  legend= fig.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1),
          ncol=3, fancybox=True) #, shadow=True)
  for lb in legend.get_texts(): lb.set_fontsize('small')
  plt.show()
#-----------------------------------------------------------------//plot

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
files= glob.glob("inn/d2b*.inn")
files.sort()
for f in files: print f
#raw_input()
#---------------------------------------------------------------------//

i0= 0

i= 0
for sFile in files:
  print '\n=========processing '+ sFile + '==========================\n'
  
  sCommand= innfile_scan_list(sFile,"TEST","COMPUTE")
  sDirout=  innfile_scan_word(sFile,"CONDITIONS","OUTPUT").replace('\\','/')
  #print sCommand, sDirout
  #raw_input()
  
  if sDebug=="1":
    for f in glob.glob("*.*"):
      if os.path.isfile(f): os.remove(f)

  os.system("%s %s %s" % (sExe,sFile,sDebug))
  
  if "SPCPATH" in sCommand:
    s=sDirout+"_molal.restab"
    lines= open(s,'r').readlines()
    labels,tData= lines2table(lines)
    iX= 0
    if "pH" in labels: iX= labels.index("pH")
    iY0= labels.index("H2O")
    plot(labels,tData,iX,iY0,False,True)
  
  if "DYN" in sCommand:
    s=sDirout+"_activ.restab"
    s=sDirout+"_minmol.restab"
    lines= open(s,'r').readlines()
    labels,tData= lines2table(lines)
    iX= 0
    if "Time/YEAR" in labels: iX= labels.index("Time/YEAR")
    if "PhiFluid" in labels: iY0= labels.index("PhiFluid")
    if "H2O" in labels: iY0= labels.index("H2O")
    plot(labels,tData,iX,iY0,False,False)
    plot(labels,tData,iX,iY0,True,True)
  
  print '\n=========done '+ sFile + '==========================\n\n'
  i += 1
  
  if(i==i0):
    raw_input()
    i= 0


