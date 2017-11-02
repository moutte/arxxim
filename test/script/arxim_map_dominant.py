import os, glob, sys, time
import pylab as plt

import exceptions
def mynum(s):
  try:
    return float(s)
  except exceptions.ValueError:
    return 0.

#-------------------------------------------------------------------INIT
os.chdir("../")
if sys.platform.startswith("win"):   #windows
  sExe= "arxim.exe"
if sys.platform.startswith("linux"): #linux
  sExe= "a.out"
sDir= os.path.join("..","bin")
sExe= os.path.join(sDir,sExe)

sDebug= "1"
sCmd= "SPC"

DEBUG= False

#----------------------------------------------------cleaning tmp_ files
for l in glob.glob("tmp_*"): os.remove(l)
if os.path.isfile("error.log"): os.remove("error.log")
#sys.exit()  
#---------------------------------------------------/cleaning tmp_ files

#-----------------------------------------------------------user-defined
fInn= "inn/map_dom_fe_ox.inn"
#Elements=["FE"]
Selection="FE"

fInn= "inn/map_dom_cr_ox.inn"
#Elements=["CR"]
Selection="CR"

Xlis= ["H"]
Ylis= ["OX"]

Xmin,Xmax,Xdelta,Xtol=  1., 13., 1., 0.02
Ymin,Ymax,Ydelta,Ytol= 40., 10., 1., 0.01

Xlabel= "pH"
Ylabel= "colog(f_O2(g))"

#---------------------------------------------------------//user-defined

sBlock= "SYSTEM"
iKeyword= 0
iValue=   3
nValue=   1

#--clean the inn file
with open(fInn,'r') as f:
  lines = f.readlines()
f=open(fInn,'w')
for line in lines:
  ll= line.strip()
  if ll=='':     continue
  if ll[0]=='!': continue
  f.write(line)
f.close()
#--/
#-----read the inn file and store line numbers and 'headers' for X and Y
Xhead= "  "
Yhead= "  "
Xindex= -1
Yindex= -1
with open(fInn,'r') as f:
  lines = f.readlines()
OK= False
for i,line in enumerate(lines):
  ww= line.split()
  if not OK:
    if ww[0].upper()==sBlock: OK=True
  else:
    if ww[0].upper()=="END": OK= False
  if OK:
    if len(ww)>iValue:
      if ww[iKeyword].upper() in Xlis:
        Xindex= i
        for j in range(iValue):
          Xhead= Xhead + ww[j] + "  "
      if ww[iKeyword].upper() in Ylis:
        Yindex= i
        for j in range(iValue):
          Yhead= Yhead + ww[j] + "  "
if 1:
  print Xindex, Xhead
  print Yindex, Yhead
  raw_input()
  #--sys.exit()
Xlis= Xindex, Xhead
Ylis= Yindex, Yhead
#--//

sArximCommand= sExe,fInn,sDebug,sCmd
fResult= "tmp_species_lact.dat"
fResult= "tmp_species_mole.dat"

if DEBUG: fMoles= open("species_moles.restab",'w')
if DEBUG: fLngam= open("species_lngam.restab",'w')
if DEBUG: fLnact= open("species_lnact.restab",'w')
spc_names= []
spc_types= []
spc_names_fe= []
lis_spc=[]
#------------------------------------------------------------------/INIT


#------------------------------------------------initialize the x,y grid
Xdim= int(round(abs(Xmax-Xmin)/Xdelta))+1
Xser= plt.linspace(Xmin,Xmax,num=Xdim)

Ydim= int(round(abs(Ymax-Ymin)/Ydelta))+1
Yser= plt.linspace(Ymin,Ymax,num=Ydim)

if 0:
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
  #-raw_input()
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
def arxim_ok(sCommand):
  Ok= False
  #
  os.system("%s %s %s %s" % (sCommand)) #---execute arxim
  #
  if os.path.isfile("error.log"):
    res= open("error.log",'r').read()
    if res.strip()=="PERFECT": Ok= True
    else: print "error.log="+ll ; raw_input()
  else:
    print "error.log NOT FOUND"
  return Ok
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
def input_modify(lis,x):
  idx,headd= lis
  with open(fInn,'r') as f:
    lines = f.readlines()
  f=open(fInn,'w')
  for i,line in enumerate(lines):
    if i==idx:
      f.write("%s" % (headd))
      for j in range(nValue): f.write("%.4g  " % (x))
      f.write('\n')
    else:
      f.write(line)
  f.close()
#----------------------------------------------//modify the include file

#-------------------------------------------------------refining routine
def refine_xy(lis,F0,F1,x0,x1,tolerance):
  F= -1
  OK= True
  #
  while abs(x0-x1)>tolerance:
    x= (x0+x1)/2.
    input_modify(lis,x)
    OK= arxim_ok(sArximCommand)
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
      OK= False
      break
  #
  if OK: x= (x0+x1)/2.
  #
  return x,F,OK
#-----------------------------------------------------//refining routine

#-----------------------------------------------map the dominant species
#----------------------------------------------------on the initial grid
tab_spc=plt.zeros((Xdim,Ydim),'int')
tab_y=plt.zeros((Xdim,Ydim),'float')
tab_x=plt.zeros((Xdim,Ydim),'float')

START= time.time()
#-----------------------------------------------------------------y-loop
for iY,y in enumerate(Yser):
  input_modify(Ylis,y) #-----------------modify the include file for y
  #---------------------------------------------------------------x-loop
  for iX,x in enumerate(Xser):
    print x,y
    #
    input_modify(Xlis,x) #---------------modify the include file for x
    OK= arxim_ok(sArximCommand) #---------------------execute arxim
    #
    #--------------------------------------------------done on first run
    if (iX,iY)==(0,0):
      with open("tmp_species.tab") as f:
        spc_names= f.readline().split()
        spc_types= f.readline().split()
      for i,name in enumerate(spc_names):
        print spc_types[i]," ",name
      if DEBUG: raw_input()
      spc_include=plt.zeros(len(spc_names), dtype=bool)
      for i in range(len(spc_names)):
        if spc_types[i]=="AQU":
          if Selection in spc_names[i]:
            spc_include[i]= True
      for (i,name) in enumerate(spc_names):
        if spc_include[i]: print spc_names[i]
      if DEBUG: raw_input()
    #------------------------------------------------//done on first run
    #
    if OK:
      spc= arxim_result(spc_include)
      if DEBUG: print spc, spc_names[spc]
      #raw_input()
      if not spc in lis_spc:
        lis_spc.append(spc)
    else:
      spc= -1
    tab_spc[iX,iY]= spc
    tab_x[iX,iY]= x
    tab_y[iX,iY]= y
    
  #raw_input()
  #-------------------------------------------------------------//x-loop
#---------------------------------------------------------------//y-loop
    
if 1:
  plt.figure()
  plt.imshow(
    tab_spc,
    aspect='auto', 
    interpolation='none')
  plt.show()

centroids=[ (0.,0.,0) for i in range(len(spc_names))]
for iY,y in enumerate(Yser):
  for iX,x in enumerate(Xser):
    i= tab_spc[iX,iY]
    x_,y_,n= centroids[i]
    x_= (n*x_+ x)/(n+1)
    y_= (n*y_+ y)/(n+1)
    n +=1
    centroids[i]= x_,y_,n
    
print tab_spc
with open("res.tab",'w') as f:
  for j in range(Ydim):
    for i in range(Xdim):
      f.write("%s\t" % str(tab_spc[i,j]))
    f.write('\n')
#--------------------------------------------//map the stable assemblage

if False:
  plt.figure()
  # plt.imshow(Data, interpolation='none')
  plt.imshow(tab_spc, aspect='auto', interpolation='none', origin='lower')
             # extent=extents(x) + extents(y), 
  plt.savefig("path_XY.png")

  # Extent defines the images max and min of the horizontal
  # and vertical values.
  # It takes four values like so: 
  # extent=[horizontal_min,horizontal_max,vertical_min,vertical_max]

  plt.show()
  sys.exit()
  
  for species in lis_spc:
    print species,spc_names[species]

  raw_input()

#---------------------------------------------------------------refining
lis_limits=  []
lis_xy_x=  []
lis_xy_y=  []
lis_false_x= []
lis_false_y= []

#--refining the x-coordinate of the species change, at fixed y
for iY in range(Ydim):
  y= tab_y[0,iY]
  input_modify(Ylis,y)
  for iX in range(1,Xdim):
    if tab_spc[iX-1,iY] != tab_spc[iX,iY]:
      print tab_x[iX,iY],tab_y[iX,iY]
      F0= tab_spc[iX-1,iY]
      F1= tab_spc[iX,  iY]
      x0= tab_x[iX-1,iY]
      x1= tab_x[iX,  iY]
      x,F,OK= refine_xy(Xlis,F0,F1,x0,x1,Xtol)
      if OK:
        if F1>F0: s= (F0,F1)
        else:     s= (F1,F0)
        if not s in lis_limits: lis_limits.append(s)
        val= s,x,y
        lis_xy_x.append(val)
      else:
        if F>0:
        #-there is a species Fm on [x0,x1] that is different from F0 & F1
        #--we need two second stage refinements
          xm=x
          Fm=F
          #--------------------second stage refinement between x0 and xm
          x,F,OK= refine_xy(Xlis,F0,Fm,x0,xm,Xtol)
          if OK:
            if F0>Fm: s= (Fm,F0)
            else:     s= (F0,Fm)
            if not s in lis_limits: lis_limits.append(s)
            val= s,x,y
            lis_xy_x.append(val)
          else:
            val= iX,iY
            lis_false_x.append(val)
          #--------------------second stage refinement between xm and x1
          x,F,OK= refine_xy(Xlis,Fm,F1,xm,x1,Xtol)
          if OK:
            if F1>Fm: s= (Fm,F1)
            else:     s= (F1,Fm)
            if not s in lis_limits: lis_limits.append(s)
            val= s,x,y
            lis_xy_x.append(val)
          else:
            val= iX,iY
            lis_false_x.append(val)
        else:
          val= iX,iY
          lis_false_x.append(val)
          continue
        
#--refining the y-coordinate of the species change, at fixed x
for iX in range(Xdim):
  x= tab_x[iX,0]
  input_modify(Xlis,x)
  for iY in range(1,Ydim):
    if tab_spc[iX,iY-1] != tab_spc[iX,iY]:
      F0= tab_spc[iX,iY-1]
      F1= tab_spc[iX,iY]
      y0= tab_y[iX,iY-1]
      y1= tab_y[iX,iY]
      y,F,OK= refine_xy(Ylis,F0,F1,y0,y1,Ytol)
      if OK:
        if F1>F0: s= (F0,F1)
        else:     s= (F1,F0)
        if not s in lis_limits: lis_limits.append(s)
        val= s,x,y
        lis_xy_y.append(val)
      else:
        if F>0:
        #-there is a paragen' Fm on [y0,y1] that is different from F0 & F1
        #--we need two second stage refinements
          ym=y
          Fm=F
          #----------------------second stage refinement between y0 and ym
          y,F,OK= refine_xy(Ylis,F0,Fm,y0,ym,Ytol)
          if OK:
            if F0>Fm: s= (Fm,F0)
            else:     s= (F0,Fm)
            if not s in lis_limits: lis_limits.append(s)
            val= s,x,y
            lis_xy_y.append(val)
          else:
            val= iX,iY
            lis_false_y.append(val)
          #----------------------second stage refinement between ym and y1
          y,F,OK= refine_xy(Xlis,Fm,F1,ym,y1,Ytol)
          if OK:
            if F1>Fm: s= (Fm,F1)
            else:     s= (F1,Fm)
            if not s in lis_limits: lis_limits.append(s)
            val= s,x,y
            lis_xy_y.append(val)
          else:
            val= iX,iY
            lis_false_y.append(val)
        else:
          val= iX,iY
          lis_false_y.append(val)
          continue

#--from the lists of points, build the lines for each "reaction"
lines= []
for reac in lis_limits:
  points= []
  for val in lis_xy_x:
    if val[0]==reac:
      point= val[1],val[2]
      points.append(point)
  for val in lis_xy_y:
    if val[0]==reac:
      point= val[1],val[2]
      points.append(point)
  points= sorted(points,key=lambda x: x[0])
  lines.append(points)

print "LIST OF SPECIES :"
for i,species in enumerate(lis_spc):
  print i, species, spc_names[species]

for reac in lis_limits:
  print reac
#raw_input()

END= time.time()
print "TIME=", END - START
#raw_input("type ENTER")

#sys.exit()
#-------------------------------------------------------------//refining

#-----------------------------------------------file name for the figure
head,tail= os.path.split(fInn)
if '.' in tail:
  figName= tail.split('.')[0]
else:
  figName= tail
if os.path.isdir("png"):
  pass
else:
  os.mkdir("png")
figName= "png/"+figName
#---------------------------------------------------------------------//
  
#--------------------------------------------------------plot XY diagram
plt.rcParams['figure.figsize']= 8,6
plt.rcParams['figure.figsize']= 8.,6.   #ratio 4/3
plt.rcParams['figure.figsize']= 5.,3.75 #ratio 4/3
plt.rcParams['figure.figsize']= 6.,6.   #ratio 4/4
plt.rcParams.update({'font.size': 12})

fig= plt.subplot(1,1,1)
symbols=['bo','go','ro','cs','mD','yd','bo','go','ro','cs','mD','yd']
symbols=['b','g','r','c','m','y']
len_symb= len(symbols)
fig.grid(color='r', linestyle='-', linewidth=0.1)
fig.grid(True)

fig.set_xlim(Xmin,Xmax)
fig.set_ylim(Ymin,Ymax)

for i,points in enumerate(lines):
  vx= []
  vy= []
  for x,y in points:
    vx.append(x)
    vy.append(y)
  fig.plot(vx, vy, symbols[i%len_symb], linestyle='-', linewidth=2.0)
  
for i,centroid in enumerate(centroids):
  x,y,n= centroid
  textstr= spc_names[i]
  fig.text(x,y,
    textstr,
    verticalalignment='center',
    horizontalalignment='center')
  
plt.xlabel(Xlabel)
plt.ylabel(Ylabel)
  
plt.savefig(figName+".png")
plt.show()
#------------------------------------------------------//plot XY diagram
