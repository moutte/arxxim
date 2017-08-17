import os, glob, sys, time
import pylab as plt
import numpy as np

import exceptions
def mynum(s):
  try:
    return float(s)
  except exceptions.ValueError:
    return 0.
    
sExe= "../bin/arxxim"
sDebug= "1"
sCmd= "SPC"

DEBUG= False

#-------------------------------------------------------------------INIT
os.chdir("../")
#----------------------------------------------------cleaning tmp_ files
for l in glob.glob("tmp_*"): os.remove(l)
if os.path.isfile("error.log"): os.remove("error.log")
#sys.exit()  
#---------------------------------------------------/cleaning tmp_ files

#------------------------------------------------------------input files
fInn= "inn/a5a3_fe_ox.inn"

iKeyword= 0
iValue=   3

lisX= ["H"]
Xmin,Xmax,Xdelta= 1., 13., 1.
tol_x= 0.02

lisY= ["OX"]
Ymin,Ymax,Ydelta= 4., 40., 2.
tol_y= 0.02

#----------------------------------------------------------//input files

fInclude= fInn.replace(".inn",".include")

with open(fInclude,'r') as f:
  lines = f.readlines()
f=open(fInclude,'w')
for line in lines:
  ll= line.strip()
  if ll=='':     continue
  if ll[0]=='!': continue
  f.write(line)
f.close()

sArximCommand= sExe,fInn,sDebug,sCmd
fResult= "tmp_species_lact.dat"
fResult= "tmp_species_mole.dat"

fMoles= open("species_moles.restab",'w')
fLngam= open("species_lngam.restab",'w')
fLnact= open("species_lnact.restab",'w')
spc_names= []
spc_types= []
spc_names_fe= []
lis_spc=[]

DEBUG= False
#------------------------------------------------------------------/INIT


#------------------------------------------------initialize the x,y grid
Xdim= int(abs(Xmax-Xmin)/Xdelta)+1
Xser= np.linspace(Xmin,Xmax,num=Xdim)

Ydim= int(abs(Ymax-Ymin)/Ydelta)+1
Yser= np.linspace(Ymin,Ymax,num=Ydim)

if False:
  print Xser
  print Yser
  sys.exit()
#----------------------------------------------//initialize the x,y grid

#------------------------------------------------------saveChemicalSpace
def readChemicalSpace():
  #--------------------------------------------------read chemical space
  with open("tmp_species.tab") as f:
    s_names= f.readline().split()
    s_types= f.readline().split()
  for name in s_names: print name
  if DEBUG: raw_input()
  #------------------------------------------------//read chemical space
  return s_names,s_types

def saveChemicalSpace(s_names):
  #--------------------------------------------------------write headers
  fMoles.write("%s\t" % "Step")
  fLnact.write("%s\t" % "Step")
  for w in s_names:
    fMoles.write("%s\t" % w)
    fLnact.write("%s\t" % w)
  fMoles.write('\n')
  fLnact.write('\n')
  #------------------------------------------------------//write headers
#----------------------------------------------------//saveChemicalSpace

#--------------------------------------------------------------arxim run
def arxim_execute(sCommand):
  OK= True
  #
  os.system("%s %s %s %s" % (sCommand)) #---execute arxim
  #
  if os.path.isfile("error.log"):
     with open("error.log",'r') as f:
       ll= f.readline().strip()
     if ll!="PERFECT":
       print "error.log="+ll
       raw_input()
       OK= False
  else:
    OK= False
  return OK
#------------------------------------------------------------//arxim run

#----------------------------------------------------------arxim results
def arxim_result(s_include):
  with open(fResult,'r') as f:
    ww= f.readline().split()
  #
  max_idx = -1
  max_val = float('-inf')
  for i,w in enumerate(ww):
    #if "FE" in spc_names[i]:
    if s_include[i]:
      val= mynum(w)
      if val > max_val:
        max_val = val
        max_idx = i
    else:
      continue
  #
  return max_idx
#----------------------------------------------------------arxim results

#------------------------------------------------modify the include file
def include_modify(lis,x):
# global variables used: fInclude, iKeyword, iValue
  with open(fInclude,'r') as f:
    lines = f.readlines()
  f=open(fInclude,'w')
  for line in lines:
    ww= line.split()
    if len(ww)>iValue and ww[iKeyword] in lis:
      for i in range(iValue):
        f.write("  %s" % (ww[i]))
      f.write("  %.4g\n" % (x))
    else:
      f.write(line)
  f.close()
#----------------------------------------------//modify the include file

#-------------------------------------------------------refining routine
def refine_xy(lis,F0,F1,x0,x1,tolerance):
  OK= True
  #
  while abs(x0-x1)>tolerance:
    x= (x0+x1)/2.
    include_modify(lis,x)
    OK= arxim_execute(sArximCommand)
    if OK:
      F= arxim_result(spc_include)
      if F in lis_spc:
        # if the species at x is the same as the species at x0,
        # then the species change occurs between x and x1, 
        # then x becomes the new x0
        if   F==F0: x0= x
        elif F==F1: x1= x
        else:
          print "F,F0,F1=",F,F0,F1
          OK= False
          break
      else:
        lis_spc.append(F)
        print "NEW:",F
        OK= False
        break
    else:
      print "error in refine"
      break
    # print x0,x1
    # raw_input()
  x= (x0+x1)/2.
  #
  return x,OK
#-----------------------------------------------------//refining routine

#-----------------------------------------------map the dominant species
#----------------------------------------------------on the initial grid
tab_spc=plt.zeros((Xdim,Ydim),'int')
tab_y=plt.zeros((Xdim,Ydim),'float')
tab_x=plt.zeros((Xdim,Ydim),'float')
first= True

#--vars for refining function
lis_reac=  []
lis_xy=    []
lis_false_x= []
lis_false_y= []
#--
START= time.time()
#-----------------------------------------------------------------y-loop
for iY,y in enumerate(Yser):
  include_modify(lisY,y) #-------------------modify the include file for y
  #---------------------------------------------------------------x-loop
  for iX,x in enumerate(Xser):
    #print iY,iX
    #
    include_modify(lisX,x) #-----------------modify the include file for x
    OK= arxim_execute(sArximCommand) #---------------------execute arxim
    #
    #-------------------------------------------------for first run only
    if first:
      with open("tmp_species.tab") as f:
        spc_names= f.readline().split()
        spc_types= f.readline().split()
      for name in spc_names: print name
      if DEBUG: raw_input()
      spc_include=np.zeros(len(spc_names), dtype=bool)
      for i in range(len(spc_names)):
        if "FE" in spc_names[i]:
          spc_include[i]= True
      for (i,name) in enumerate(spc_names):
        if spc_include[i]: print spc_names[i]
      if DEBUG: raw_input()
      first= False
    #-----------------------------------------------//for first run only
    #
    if OK:
      spc= arxim_result(spc_include)
      if DEBUG: print spc, spc_names[spc]
      #raw_input()
      if not spc in lis_spc:
        lis_spc.append(spc)
      tab_spc[iX,iY]= spc
      tab_x[iX,iY]= x
      tab_y[iX,iY]= y
      #--refining the x-coordinate of a species change
      if iX>0:
        if tab_spc[iX-1,iY] != tab_spc[iX,iY]:
          if DEBUG: print iX,iY
          F0= tab_spc[iX-1,iY]
          F1= tab_spc[iX,  iY]
          print F0,F1
          x0= tab_x[iX-1,iY]
          x1= tab_x[iX,  iY]
          x,OK= refine_xy(lisX,F0,F1,x0,x1,tol_x)
          if OK:
            if F0>F1: s= str(F1)+'='+str(F0)
            else:     s= str(F0)+'='+str(F1)
            if not s in lis_reac: lis_reac.append(s)
            val= s,x,y
            #lis_xy_x.append(val)
            lis_xy.append(val)
          else:
            val= iX,iY
            lis_false_x.append(val)
      #--refining the y-coordinate of a species change        
      if iY>0:
        if tab_spc[iX,iY-1] != tab_spc[iX,iY]:
          F0= tab_spc[iX,iY-1]
          F1= tab_spc[iX,iY]
          print F0,F1
          y0= tab_y[iX,iY-1]
          y1= tab_y[iX,iY]
          y,OK= refine_xy(lisY,F0,F1,y0,y1,tol_y)
          if OK:
            if F0>F1: s= str(F1)+'='+str(F0)
            else:     s= str(F0)+'='+str(F1)
            if not s in lis_reac: lis_reac.append(s)
            val= s,x,y
            #lis_xy_y.append(val)
            lis_xy.append(val)
          else:
            val= iX,iY
            lis_false_y.append(val)
    else:
      spc= -1
      tab_spc[iX,iY]= spc
      tab_x[iX,iY]= x
      tab_y[iX,iY]= y
    
  #raw_input()
  #-------------------------------------------------------------//x-loop
#---------------------------------------------------------------//y-loop
if DEBUG:
  print tab_spc
  with open("res.tab",'w') as f:
    for j in range(Ydim):
      for i in range(Xdim):
        f.write("%s\t" % str(tab_spc[i,j]))
      f.write('\n')
#--------------------------------------------//map the stable assemblage

#--from the lists of points, build the lines for each "reaction"
lines= []
for reac in lis_reac:
  points= []
  for val in lis_xy:
    if val[0]==reac:
      point= val[1],val[2]
      points.append(point)
  points= sorted(points,key=lambda x: x[0])
  lines.append(points)

for species in lis_spc:
  print species
for reac in lis_reac:
  print reac
#raw_input()

END= time.time()
print "TIME=", END-START
if DEBUG: raw_input("type ENTER")

#sys.exit()
#-------------------------------------------------------------//refining

#--------------------------------------------------------plot XY diagram
plt.rcParams['figure.figsize']= 8,6
fig= plt.subplot(1,1,1)
symbols=['bo','go','ro','cs','mD','yd','bo','go','ro','cs','mD','yd']
fig.grid(color='r', linestyle='-', linewidth=0.2)
fig.grid(True)
  
for i,points in enumerate(lines):
  vx= []
  vy= []
  for x,y in points:
    vx.append(x)
    vy.append(y)
  fig.plot(vx, vy, symbols[i], linestyle='-', linewidth=1.0)
  
plt.savefig("000_map_arxim_2"+".png")
plt.show()
#------------------------------------------------------//plot XY diagram

