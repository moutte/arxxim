import os
import glob
import pylab as plt
import numpy as np
import sys

sExe= "..\\bin\\arxim.exe"  #windows
sExe= "../bin/arxim"        #linux
sDebug= "1"
sCmd= "GEM"


#-------------------------------------------------------------------INIT
os.chdir("../")
#----------------------------------------------------cleaning tmp_ files
for l in glob.glob("tmp_*"): os.remove(l)
if os.path.isfile("error.log"): os.remove("error.log")
#sys.exit()  
#---------------------------------------------------/cleaning tmp_ files

#------------------------------------------------------------input files
fInn= "inn/map2a_activ.inn"
fInn= "inn/map2b_activ.inn"

'''
the include block to be modified:
SPECIES
MIN	VIRTUAL	SIO2_BUFF	SI(1)O(2)	2700.	3.5
MIN	VIRTUAL	KOH_BUFF	K(1)O(1)H(1)	2700.	-1.5
END
Keyword is SIO2_BUFF,  its index is 2  -> iKeyword= 2
Value is the last one, its index is 5  -> iValue=   5
'''
iKeyword= 2
iValue=   5
nValue=   8

Xlis= ["SIO2_BUFF"]
Ylis= ["KOH_BUFF"]

Xlabel= "-log[SiO2]"
Ylabel= "-log[K+]/[H+]"

Xmin,Xmax,Xdelta,Xtol= 3.,  2., 0.05,  0.01
Ymin,Ymax,Ydelta,Ytol= -2., -4., 0.1,  0.01

Xmin,Xmax,Xdelta,Xtol= 5.,  2., 0.25,  0.02
Ymin,Ymax,Ydelta,Ytol= 0., -8., 0.25,  0.02

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
fResult= "tmp_gem.tab"
#------------------------------------------------------------------/INIT


#------------------------------------------------initialize the x,y grid
Ydim= int(abs(Ymax-Ymin)/Ydelta)+1
Yser= np.linspace(Ymin,Ymax,num=Ydim)

Xdim= int(abs(Xmax-Xmin)/Xdelta)+1
Xser= np.linspace(Xmin,Xmax,num=Xdim)

if 0:
  print Xser
  print Yser
  sys.exit()
#----------------------------------------------//initialize the x,y grid

#--------------------------------------------------------------arxim run
def arxim_ok(sCommand):
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
    print "error.log: NOT FOUND"
    OK= False
  #
  return OK
#------------------------------------------------------------//arxim run
  
#---------------------------------------------------------read results
def arxim_result(fResult): #(sCommand):
  with open(fResult,'r') as f:
    lines= f.readlines()
  lis=""
  for i,line in enumerate(lines):
    mineral= line.split()[0]
    if i>0:
      lis= lis + mineral + "="
  return lis 
#--------------------------------------------------------/read results

#------------------------------------------------modify the include file
def input_modify(lis,x):
  with open(fInclude,'r') as f:
    lines = f.readlines()
  f=open(fInclude,'w')
  for line in lines:
    ww= line.split()
    if len(ww)>iValue and ww[iKeyword] in lis:
      for i in range(iValue): f.write("  %s"   % (ww[i]))
      for i in range(nValue): f.write("  %.4g" % (x))
      f.write('\n')
    else:
      f.write(line)
  f.close()
#----------------------------------------------//modify the include file

#----------------------------------------------map the stable assemblage
#----------------------------------------------------on the initial grid
tab_paragen=plt.zeros((Xdim,Ydim),'int')
tab_y=plt.zeros((Xdim,Ydim),'float')
tab_x=plt.zeros((Xdim,Ydim),'float')
lis_paragen=[]
#-----------------------------------------------------------------y-loop
for iY,y in enumerate(Yser):
  input_modify(Ylis,y) #-------------------modify the include file for y
  #---------------------------------------------------------------x-loop
  for iX,x in enumerate(Xser):
    print iY,iX
    #
    input_modify(Xlis,x) #-----------------modify the include file for x
    OK= arxim_ok(sArximCommand) #--------------------------execute arxim
    if OK:
      paragen= arxim_result(fResult) #-------------------read arxim result
      print paragen
      #raw_input()
      if not paragen in lis_paragen:
        lis_paragen.append(paragen)
      j= lis_paragen.index(paragen)
      tab_paragen[iX,iY]= j
      tab_x[iX,iY]= x
      tab_y[iX,iY]= y
  #-------------------------------------------------------------//x-loop
#---------------------------------------------------------------//y-loop
    
#---------------------------------compute the positions for field labels
centroids=[ (0.,0.,0) for i in range(len(lis_paragen))]
for iY,y in enumerate(Yser):
  for iX,x in enumerate(Xser):
    i= tab_paragen[iX,iY]
    x_,y_,n= centroids[i]
    x_= (n*x_+ x)/(n+1)
    y_= (n*y_+ y)/(n+1)
    n +=1
    centroids[i]= x_,y_,n
if True:
  for centroid in centroids:
    print centroid
#-------------------------------//compute the positions for field labels
    
print tab_paragen
#--------------------------------------------//map the stable assemblage
#sys.exit()

#-------------------------------------------------------refining routine
def refine_xy(lis,F0,F1,x0,x1,tolerance):
  OK= True
  #
  while abs(x0-x1)>tolerance:
    x= (x0+x1) /2.
    input_modify(lis,x)
    OK= arxim_ok(sArximCommand)
    if OK:
      paragen= arxim_result(fResult)
      if paragen in lis_paragen:
        F= lis_paragen.index(paragen)
        # if the paragen at x is the same as the paragen at x0,
        # then the paragen change occurs between x and x1, 
        # then x becomes the new x0
        if   F==F0: x0= x
        elif F==F1: x1= x
        else:
          print "F,F0,F1=",F,F0,F1
          OK= False
          break
      else:
        lis_paragen.append(paragen)
        F= lis_paragen.index(paragen)
        #print "NEW:",paragen
        OK= False
        break
      # print x0,x1
      # raw_input()
  if OK: x= (x0+x1)/2.
  #
  return x,F,OK
#-----------------------------------------------------//refining routine

#---------------------------------------------------------------refining
lis_reac=  []
lis_xy_x=  []
lis_xy_y=  []
lis_false_x= []
lis_false_y= []

#--refining the x-coordinate of the paragenesis change, at fixed y
for iY in range(Ydim):
  y= tab_y[0,iY]
  input_modify(Ylis,y)
  for iX in range(1,Xdim):
    if tab_paragen[iX-1,iY] != tab_paragen[iX,iY]:
      F0= tab_paragen[iX-1,iY]
      F1= tab_paragen[iX,  iY]
      x0= tab_x[iX-1,iY]
      x1= tab_x[iX,  iY]
      x,F,OK= refine_xy(Xlis,F0,F1,x0,x1,Xtol)
      if OK:
        if F0>F1: s= str(F1)+'='+str(F0)
        else:     s= str(F0)+'='+str(F1)
        if not s in lis_reac: lis_reac.append(s)
        val= s,x,y
        lis_xy_x.append(val)
      else:
      #-there is a paragen' Fm on [x0,x1] that is different from F0 & F1
      #--we need two second stage refinements
        xm=x
        Fm=F
        #----------------------second stage refinement between x0 and xm
        x,F,OK= refine_xy(Xlis,F0,Fm,x0,xm,Xtol)
        if OK:
          if F0>Fm: s= str(Fm)+'='+str(F0)
          else:     s= str(F0)+'='+str(Fm)
          if not s in lis_reac: lis_reac.append(s)
          val= s,x,y
          lis_xy_x.append(val)
        else:
          val= iX,iY
          lis_false_x.append(val)
        #----------------------second stage refinement between xm and x1
        x,F,OK= refine_xy(Xlis,Fm,F1,xm,x1,Xtol)
        if OK:
          if F1>Fm: s= str(Fm)+'='+str(F1)
          else:     s= str(F1)+'='+str(Fm)
          if not s in lis_reac: lis_reac.append(s)
          val= s,x,y
          lis_xy_x.append(val)
        else:
          val= iX,iY
          lis_false_x.append(val)
        #
        
#--refining the y-coordinate of the paragenesis change, at fixed x
for iX in range(Xdim):
  x= tab_x[iX,0]
  input_modify(Xlis,x)
  for iY in range(1,Ydim):
    if tab_paragen[iX,iY-1] != tab_paragen[iX,iY]:
      F0= tab_paragen[iX,iY-1]
      F1= tab_paragen[iX,iY]
      y0= tab_y[iX,iY-1]
      y1= tab_y[iX,iY]
      y,F,OK= refine_xy(Ylis,F0,F1,y0,y1,Ytol)
      if OK:
        if F0>F1: s= str(F1)+'='+str(F0)
        else:     s= str(F0)+'='+str(F1)
        if not s in lis_reac: lis_reac.append(s)
        val= s,x,y
        lis_xy_y.append(val)
      else:
      #-there is a paragen' Fm on [y0,y1] that is different from F0 & F1
      #--we need two second stage refinements
        ym=y
        Fm=F
        #----------------------second stage refinement between y0 and ym
        y,F,OK= refine_xy(Ylis,F0,Fm,y0,ym,Ytol)
        if OK:
          if F0>Fm: s= str(Fm)+'='+str(F0)
          else:     s= str(F0)+'='+str(Fm)
          if not s in lis_reac: lis_reac.append(s)
          val= s,x,y
          lis_xy_y.append(val)
        else:
          val= iX,iY
          lis_false_y.append(val)
        #----------------------second stage refinement between ym and y1
        y,F,OK= refine_xy(Xlis,Fm,F1,ym,y1,Ytol)
        if OK:
          if F1>Fm: s= str(Fm)+'='+str(F1)
          else:     s= str(F1)+'='+str(Fm)
          if not s in lis_reac: lis_reac.append(s)
          val= s,x,y
          lis_xy_y.append(val)
        else:
          val= iX,iY
          lis_false_y.append(val)

#--
lines= []
for reac in lis_reac:
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

for paragen in lis_paragen:
  print paragen
for reac in lis_reac:
  print reac
#raw_input()

#sys.exit()
#-------------------------------------------------------------//refining

#-----------------------------------------processing the paragenese list
#----------i.e. find the phases that are common to all parageneses found
lis_phase= []
for paragen in lis_paragen:
  ww= paragen.split('=')
  for w in ww:
    if not w in lis_phase:
      lis_phase.append(w)
isPartout=[True for i in range(len(lis_phase))]
for paragen in lis_paragen:
  ww= paragen.split('=')
  for i,phase in enumerate(lis_phase):
    if not phase in ww:
      isPartout[i]= False
phase_Partout=[]
for i,phase in enumerate(lis_phase):
  if isPartout[i]:
    phase_Partout.append(phase)
for i,paragen in enumerate(lis_paragen):
  ww= paragen.split('=')
  newp= ""
  for w in ww:
    if not w in phase_Partout: 
      newp= newp + w + '='
  lis_paragen[i]= newp
#--print simplified paragesis
print "SIMPLIFIED PARAGENESIS="
for i,paragen in enumerate(lis_paragen):
  print i, paragen
print "PHASES FOUND IN ALL PARAGENESIS="
for phase in phase_Partout:
  print phase
#--
#---------------------------------------//processing the paragenese list

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

#-----------------------------------------------------------plot diagram
plt.rcParams['figure.figsize']= 10,8
fig= plt.subplot(1,1,1)
symbols=['bo','go','ro','cs','mD','yd','bo','go','ro','cs','mD','yd']
fig.grid(color='r', linestyle='-', linewidth=0.2)
fig.grid(True)
  
fig.set_xlim(Xmin,Xmax)
fig.set_ylim(Ymin,Ymax)

for i,points in enumerate(lines):
  vx= []
  vy= []
  for x,y in points:
    vx.append(x)
    vy.append(y)
  fig.plot(vx, vy, symbols[i], linestyle='-', linewidth=1.0)
  
for i,centroid in enumerate(centroids):
  x,y,n= centroid
  textstr= lis_paragen[i].replace('=','\n')
  fig.text(x,y,
    textstr,
    horizontalalignment='center',
    verticalalignment='center')

plt.xlabel(Xlabel)
plt.ylabel(Ylabel)
  
plt.savefig(figName+".png")
plt.show()
#---------------------------------------------------------//plot diagram
