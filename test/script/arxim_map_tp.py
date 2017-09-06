import os, glob, sys
import pylab as plt

#windows
sExe= "arxim.exe"
sExe= os.path.join("..","bin",sExe)

#linux
sExe= "arxxim"
sExe= "arx-bis"
sExe= os.path.join("..","bin",sExe)

sExe= "arx-basis"
sExe= os.path.join("..","..","arx-basis","bin",sExe)

sExe= "arx-win"
sExe= os.path.join("..","..","arx-win","bin",sExe)

sDebug= "2"
sCmd=  "GEM"

#-------------------------------------------------------------------INIT
os.chdir("../")
#----------------------------------------------------cleaning tmp_ files
for l in glob.glob("tmp_*"): os.remove(l)
if os.path.isfile("error.log"): os.remove("error.log")
#sys.exit()  
#---------------------------------------------------/cleaning tmp_ files

#-----------------------------------------------------------user defined
fInn= "inn/map3c_tp.inn"

Xlis= ["TDGC"]  # must be uppercase ...
Ylis= ["PBAR"]  # must be uppercase ...

Xlabel= "T /DGC"
Ylabel= "P /BAR"

Xmin,Xmax,Xdelta,Xtol= 300., 900.,   50.,  5.
Ymin,Ymax,Ydelta,Ytol= 500., 6500., 500., 10.

Xmin,Xmax,Xdelta,Xtol=  500.,   600., 10.,   1.
Ymin,Ymax,Ydelta,Ytol= 4000.,  5000., 50.,   5.

Xmin,Xmax,Xdelta,Xtol=  400.,   700.,  20.,   5.
Ymin,Ymax,Ydelta,Ytol= 3000.,  6000., 200.,  20.

#---------------------------------------------------------//user defined

'''
the include block to be modified:
SYSTEM.GEM
  TDGC  700
  PBAR  5638
  SIO2   SiO2  1.0
  AL2O3  Al2O3 1.0
END
Keyword is TDGC / PBAR,  its index is 0  -> iKeyword= 0
Value is 1200 / 400,     its index is 1  -> iValue=   1
'''

iKeyword= 0
iValue=   1

fInclude= fInn.replace(".inn",".include")
#--clean the include file
with open(fInclude,'r') as f:
  lines = f.readlines()
#--clean the include file
f=open(fInclude,'w')
for line in lines:
  ll= line.strip()
  if ll=='':     continue
  if ll[0]=='!': continue
  f.write(line)
f.close()

#--store line numbers and header for X and Y
Xhead= "  "
Yhead= "  "
with open(fInclude,'r') as f:
  lines = f.readlines()
for i,line in enumerate(lines):
  ww= line.split()
  if len(ww)>iValue:
    if ww[iKeyword].upper() in Xlis:
      Xindex= i
      for j in range(iValue):
        Xhead= Xhead + ww[j] + "  "
    if ww[iKeyword].upper() in Ylis:
      Yindex= i
      for j in range(iValue):
        Yhead= Yhead + ww[j] + "  "
if 0:
  print Xindex, Xhead
  print Yindex, Yhead
  raw_input()
Xlis= Xindex, Xhead
Ylis= Yindex, Yhead
#--/

#--write fInclude to the global variable IncLines -> don't work ...
with open(fInclude,'r') as f:
  IncLines = f.readlines()
#--/

sArximCommand= sExe,fInn,sDebug,sCmd
# the arxim executable (sExe) accepts 3 arguments, successively:
# > fInn: the path of the arxim script
# > sDebug: an integer that controls the screen output at runtime
# for the present case, sDebug=1, which restricts the screen output
# > sCmd: the calcualtion command, e.g. SPC, GEM, EQU, ...
fResult= "tmp_gem.tab"
#------------------------------------------------------------------/INIT

#DEBUG= True
#fLog= open("log",'w',0)
fParagen= open("paragen.txt",'w',0)

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

#------------------------------------------------modify the include file
def include_modify(lis,x):
  idx,headd= lis
  with open(fInclude,'r') as f:
    lines = f.readlines()
  f=open(fInclude,'w')
  for i,line in enumerate(lines):
    if i==idx: f.write("%s  %.4g\n" % (headd,x))
    else:      f.write(line)
    #old version
    if 0:
      ww= line.split()
      if len(ww)>iValue:
        if ww[iKeyword].upper() in lis:
          for i in range(iValue): f.write("  %s" % (ww[i]))
          f.write("  %.4g\n" % (x))
        else: f.write(line)
      else:
        f.write(line)
    #--/
  f.close()
#----------------------------------------------//modify the include file

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
  
#-----------------------------------------------------------read results
def arxim_result(fResult): #(sCommand):
  with open(fResult,'r') as f:
    lines= f.readlines()
  lis=""
  for i,line in enumerate(lines):
    mineral= line.split()[0]
    if i>0:
      lis= lis + mineral + "="
  return lis 
#----------------------------------------------------------/read results

#-------------------------------------------------------refining routine
def refine_xy(lis,F0,F1,x0,x1,tolerance):
  OK= True
  F= -1
  #
  while abs(x0-x1)>tolerance:
    x= (x0+x1) /2.
    #if DEBUG: fLog.write("%s  %.4g\n" % (lis[0],x))
    include_modify(lis,x)
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
          print "F,F0,F1=", F,F0,F1
          OK= False
          break
      else:
        lis_paragen.append(paragen)
        fParagen.write("%s\n" % paragen)
        F= lis_paragen.index(paragen)
        print "F,F0,F1=", F,F0,F1
        OK= False
        break
    # print x0,x1
    # raw_input()
  if OK: x= (x0+x1)/2.
  #
  return x,F,OK
#-----------------------------------------------------//refining routine

#----------------------------------------------map the stable assemblage
#----------------------------------------------------on the initial grid
tab_paragen=plt.zeros((Xdim,Ydim),'int')
tab_y=plt.zeros((Xdim,Ydim),'float')
tab_x=plt.zeros((Xdim,Ydim),'float')
lis_paragen=[]

#--for refining:
lis_reac=  []
lis_xy_x=  []
lis_xy_y=  []
lis_false_x= []
lis_false_y= []
#--
#-----------------------------------------------------------------y-loop
for iY,y in enumerate(Yser):
  include_modify(Ylis,y) #-----------------modify the include file for y
  #---------------------------------------------------------------x-loop
  for iX,x in enumerate(Xser):
    print x,y
    #
    include_modify(Xlis,x) #---------------modify the include file for x
    OK= arxim_ok(sArximCommand) #--------------------------execute arxim
    if OK:
      paragen= arxim_result(fResult) #-----------------read arxim result
      if not paragen in lis_paragen:
        print paragen
        #raw_input()
        lis_paragen.append(paragen)
        fParagen.write("%s\n" % paragen)
        #centroids.append((x,y,0))
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
if False:
  for centroid in centroids:
    print centroid
#-------------------------------//compute the positions for field labels
    
print tab_paragen
#--------------------------------------------//map the stable assemblage
#sys.exit()

#---------------------------------------------------------------refining
#--refining the x-coordinate of the paragenesis change, at fixed y
for iY in range(Ydim):
  y= tab_y[0,iY]
  include_modify(Ylis,y)
  for iX in range(1,Xdim):
    if tab_paragen[iX-1,iY] != tab_paragen[iX,iY]:
      F0= tab_paragen[iX-1,iY]
      F1= tab_paragen[iX,iY]
      x0= tab_x[iX-1,iY]
      x1= tab_x[iX,iY]
      x,F,OK= refine_xy(Xlis,F0,F1,x0,x1,Xtol)
      if OK:
        if F0>F1 : s= str(F1)+'='+str(F0)
        else     : s= str(F0)+'='+str(F1)
        if not s in lis_reac: lis_reac.append(s)
        val= s,x,y
        lis_xy_x.append(val)
      else:
        if F>0:
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
        else: #F<0
          val= iX,iY
          lis_false_y.append(val)
        
#--refining the y-coordinate of the paragenesis change, at fixed x
for iX in range(Xdim):
  x= tab_x[iX,0]
  include_modify(Xlis,x)
  for iY in range(1,Ydim):
    if tab_paragen[iX,iY-1] != tab_paragen[iX,iY]:
      F0= tab_paragen[iX,iY-1]
      F1= tab_paragen[iX,iY]
      y0= tab_y[iX,iY-1]
      y1= tab_y[iX,iY]
      y,F,OK= refine_xy(Ylis,F0,F1,y0,y1,Ytol)
      if OK:
        if F0>F1 : s= str(F1)+'='+str(F0)
        else     : s= str(F0)+'='+str(F1)
        if not s in lis_reac: lis_reac.append(s)
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
            if F0>Fm: s= str(Fm)+'='+str(F0)
            else:     s= str(F0)+'='+str(Fm)
            if not s in lis_reac: lis_reac.append(s)
            val= s,x,y
            lis_xy_y.append(val)
          else:
            val= iX,iY
            lis_false_y.append(val)
          #----------------------second stage refinement between ym and y1
          y,F,OK= refine_xy(Ylis,Fm,F1,ym,y1,Ytol)
          if OK:
            if F1>Fm: s= str(Fm)+'='+str(F1)
            else:     s= str(F1)+'='+str(Fm)
            if not s in lis_reac: lis_reac.append(s)
            val= s,x,y
            lis_xy_y.append(val)
          else:
            val= iX,iY
            lis_false_y.append(val)
        else: #F<0
          val= iX,iY
          lis_false_y.append(val)
          
#--build the lines for each reaction from the lists of points
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

for i,paragen in enumerate(lis_paragen):
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
      if len(w)>7: w=w[0:7]
      newp= newp + w + '='
  lis_paragen[i]= newp 
#--print simplified paragesis
if len(phase_Partout):
  print "PHASES FOUND IN ALL PARAGENESIS="
  for phase in phase_Partout:
    print phase
print "(SIMPLIFIED) PARAGENESIS="
for i,paragen in enumerate(lis_paragen):
  print i, paragen
#---------------------------------------//processing the paragenese list

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
plt.rcParams['figure.figsize']= 8.,8.   #ratio 4/4
plt.rcParams.update({'font.size': 9})

fig= plt.subplot(1,1,1)
symbols=['bo','go','ro','cs','mD','yd']
symbols=['b','g','r','c','m','y']
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
  j= i%6
  fig.plot(vx, vy, symbols[j], linestyle='-', linewidth=2.0)

for i,centroid in enumerate(centroids):
  x,y,n= centroid
  textstr= lis_paragen[i].replace('=','\n')
  fig.text(x,y,
    textstr,
    verticalalignment='center',
    horizontalalignment='center')
  
plt.xlabel(Xlabel)
plt.ylabel(Ylabel)
  
plt.savefig(figName+".png")
plt.show()
#------------------------------------------------------//plot XY diagram
