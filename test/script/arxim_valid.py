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
def plot(labels,data): #,xLog,yLog):
  nlin,ncol= data.shape
  
  fig= plt.subplot()
  fig.grid(color='r', linestyle='-', linewidth=0.2)
  fig.grid(True)
  
  symbols=['bo','go','ro','cs','mD','yd']
  colors = ['cyan', 'lightblue', 'lightgreen', 'tan', 'pink','red', 'blue']
  lenSym= len(symbols)

  # fig.set_xlabel(titX)
  # fig.set_ylabel(titY)
  # fig.set_title(titre) #, fontsize=fontsize)
  # fig.text(min(vx),max(vy),textt,
  #  fontsize=16,ha = 'left', va = 'top')
  
  # fig.set_xlim(0,1.)
  fig.set_xlim(1e-4,1.0)
  #fig.semilogx(A[iMin:iMax,iX],A[iMin:iMax,iY],  'bo')
  
  xLog= False
  yLog= False # True # 
  
  i0= labels.index("H2O")
  if "pH" in labels:
    iX= labels.index("pH")
  else:
    iX= 0
  for i in range(ncol):
    if i>i0:
      symb= symbols[(i-i0-1)%len(symbols)]
      vx= data[:,iX]
      vy= data[:,i]
      if xLog and yLog:
        fig.loglog(vx, vy, symb, 
          linestyle='-', linewidth=1.0, label=labels[i])
      else:
        if xLog:
          fig.semilogx(vx, vy, symb, 
            linestyle='-', linewidth=1.0, label=labels[i])
        elif yLog:
          fig.semilogy(vx, vy, symb, 
            linestyle='-', linewidth=1.0, label=labels[i])
        else:
          fig.plot(vx, vy, symb, 
            linestyle='-', linewidth=1.0, label=labels[i])
      
  fig.legend(loc='upper left')
  plt.show()
#-----------------------------------------------------------------//plot

#-------------------------------------------------------------------plot
def plot_wrk(labels,data,xLog=False,yLog=False):
  nlin,ncol= data.shape
  
  fig= plt.subplot()
  fig.grid(color='r', linestyle='-', linewidth=0.2)
  fig.grid(True)
  
  fig.set_ylim(1e-4,1.0)
  
  symbols=['bo','go','ro','cs','mD','yd']
  colors = ['cyan', 'lightblue', 'lightgreen', 'tan', 'pink','red', 'blue']

  #i0= labels.index("H2O")
  i0= labels.index("PhiFluid")
  iX= 0
  if "pH" in labels: iX= labels.index("pH")
  if "Time/YEAR" in labels: iX= labels.index("Time/YEAR")
  for i in range(ncol):
    #if i>i0:
    if "PhiM_" in labels[i]:
      symb= symbols[(i-i0-1)%len(symbols)]
      labl= labels[i].replace("PhiM_","")
      vx= data[:,iX]
      vy= data[:,i]
      if xLog and yLog:
        fig.loglog(vx, vy, symb, linestyle='-',linewidth=1.0,label=labl)
      else:
        if xLog:
          fig.semilogx(vx, vy, symb, linestyle='-',linewidth=1.0,label=labl)
        elif yLog:
          fig.semilogy(vx, vy, symb, linestyle='-', linewidth=1.0,label=labl)
        else:
          fig.plot(vx, vy, sym, linestyle='-', linewidth=1.0,label=labl)
      
  legend= fig.legend(loc='upper left')
  for lb in legend.get_texts(): lb.set_fontsize('small')
  for lbl in legend.get_lines(): lbl.set_lineheight(0.8)
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
files= glob.glob("inn/d2*.inn")
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
    #plot(labels,tData) #,False,True)
    plot_wrk(labels,tData,False,True)
  
  if "DYN" in sCommand:
    s=sDirout+"_minmol.restab"
    lines= open(s,'r').readlines()
    labels,tData= lines2table(lines)
    #plot(labels,tData) #,False,True)
    plot_wrk(labels,tData,False,True)
  
  print '\n=========done '+ sFile + '==========================\n\n'
  i += 1
  
  if(i==i0):
    raw_input()
    i= 0


