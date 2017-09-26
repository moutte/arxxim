import os, glob, sys, time
import pylab as plt

os.chdir("../")

#-------------------------------------------------------------------INIT
if sys.platform.startswith("win"):
#windows
  sExe= "arxim.exe"
  sExe= os.path.join("..","bin",sExe)

if sys.platform.startswith("linux"):
#linux
  sExe= "a.out"
  sExe= os.path.join("..","..","arx-git","bin",sExe)
  
  sExe= "arx_optim"
  sExe= "arx_debug"
  sExe= "a.out"
  sExe= os.path.join("..","bin",sExe)

sDebug= "1"
sCmd=  "GEM"

#----------------------------------------------------cleaning tmp_ files
for l in glob.glob("tmp_*"): os.remove(l)
if os.path.isfile("error.log"): os.remove("error.log")
#sys.exit()  
#---------------------------------------------------/cleaning tmp_ files

#-----------------------------------------------------------user-defined
fInn= "inn/map_xt_fsp_01.inn"

mixture_names=[]

Xlis= ["OR_","AB_"]         # must be uppercase ...
Ylis= ["TDGC","PBAR"]       # must be uppercase ...
Xlabel= "AB"
Ylabel= "T/DGC"

Xmin,Xmax,Xdelta,Xtol=  2.,  98., 4.,  1.
Ymin,Ymax,Ydelta,Ytol=  400.,   700.,  20.,   5.

def Xfunc(x): return 100. - x
def Yfunc(y): return 1000. #+ (y-400.)*10.
#---------------------------------------------------------//user-defined

'''
the include block to be modified:
SYSTEM.GEM
  TDGC  1200
  PBAR  400
  SIO2   SiO2  1.0
  AL2O3  Al2O3 1.0
END
for T,P:
Keyword is TDGC / PBAR,  its index is 0  -> iKeyword= 0
Value is 1200 / 400,     its index is 1  -> iValue=   1
for composition:
iKeyword= 0
iValue=   2
'''
# keywords must be uppercase
Xlis=[x.upper() for x in Xlis]
Ylis=[y.upper() for y in Ylis]

sBlock= "SYSTEM.GEM"

XKeyword= 0
XValue=   2

YKeyword= 0
YValue=   1

#--clean the inn file
with open(fInn,'r') as f:
  lines = f.readlines()
#--clean the include file
f=open(fInn,'w')
for line in lines:
  ll= line.strip()
  if ll=='':     continue
  if ll[0]=='!': continue
  f.write(line)
f.close()

#--store line numbers and header for X and Y
Xindex= [-1 for s in Xlis] 
Yindex= [-1 for s in Ylis]

Xhead= ["  " for s in Xlis]
Yhead= ["  " for s in Ylis]

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
    if len(ww)>XValue:
      s= ww[XKeyword].upper()
      if s in Xlis:
        j= Xlis.index(s)
        Xindex[j]= i
        for k in range(XValue):
          Xhead[j]= Xhead[j] + ww[k] + "  " 
    if len(ww)>YValue:
      s= ww[YKeyword].upper()
      if s in Ylis:
        j= Ylis.index(s)
        Yindex[j]= i
        for k in range(YValue):
          Yhead[j]= Yhead[j] + ww[k] + "  "
if 1:
  print Xindex, Xhead
  print Yindex, Yhead
  #sys.exit()

Xlis= Xindex, Xhead
Ylis= Yindex, Yhead
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

#-----------------------------------------------------read_mixture_names
def read_mixture_names():
  with open("tmp_mixmodels_names.tab") as f:
    s= f.readline()
  s_names= s.split()
  return s_names
#--//
def save_mixture_moles(x,y):
  with open("tmp_mixmodels_moles.tab") as f:
    s= f.readline()
  #fMixtur.write(s)
  tMixtur.append((x,y,s))
 
#------------------------------------------------initialize the x,y grid
Xdim= int(round(abs(Xmax-Xmin)/Xdelta))+1
Xser= plt.linspace(Xmin,Xmax,num=Xdim)

Ydim= int(round(abs(Ymax-Ymin)/Ydelta))+1
Yser= plt.linspace(Ymin,Ymax,num=Ydim)

if 1:
  Xser2= [Xfunc(x) for x in Xser]
  Yser2= [Yfunc(y) for y in Yser]
  print Xser
  print Xser2
  print Yser
  print Yser2
  #sys.exit()
#----------------------------------------------//initialize the x,y grid

if 0:
  def test_func_in_func(x,f):
    y= f(x)
    return y
  x= 0.95
  y= test_func_in_func(x,Xfunc)
  print x,y
  sys.exit()
#------------------------------------------------modify the include file
def include_modify(lis,func,x):
  idx,headd= lis
  with open(fInn,'r') as f:
    lines = f.readlines()
  f=open(fInn,'w')
  for i,line in enumerate(lines):
    if i==idx[0]:
      f.write("%s  %.4g\n" % (headd[0],x))
    elif i==idx[1]:
      y= func(x)
      f.write("%s  %.4g\n" % (headd[1],y))
    else:
      f.write(line)
  f.close()
#----------------------------------------------//modify the include file

if 0:
  x= 0.45
  include_modify(Xlis,Xfunc,x)
  sys.exit()
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
def refine_xy(lis,func,F0,F1,x0,x1,tolerance):
  OK= True
  F= -1
  #
  while abs(x0-x1)>tolerance:
    x= (x0+x1) /2.
    #if DEBUG: fLog.write("%s  %.4g\n" % (lis[0],x))
    include_modify(lis,func,x)
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

tMixtur= [] #table of mixture compositions

#--for refining:
lis_reac=  []
lis_xy_x=  []
lis_xy_y=  []
lis_false_x= []
lis_false_y= []
#--
START= time.time()
#-----------------------------------------------------------------y-loop
for iY,y in enumerate(Yser):
  include_modify(Ylis,Yfunc,y) #-----------modify the include file for y
  #---------------------------------------------------------------x-loop
  for iX,x in enumerate(Xser):
    print x,y
    #
    include_modify(Xlis,Xfunc,x) #---------modify the include file for x
    OK= arxim_ok(sArximCommand) #--------------------------execute arxim
    if OK:
      #
      if iY==0 and iX==0:
        mixture_names= read_mixture_names()
        for s in mixture_names: print s
        raw_input("type ENTER")
      save_mixture_moles(x,y)
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
#------refining the x-coordinate of the paragenesis change, at fixed y--
for iY in range(Ydim):
  y= tab_y[0,iY]
  include_modify(Ylis,Yfunc,y)
  for iX in range(1,Xdim):
    if tab_paragen[iX-1,iY] != tab_paragen[iX,iY]:
      F0= tab_paragen[iX-1,iY]
      F1= tab_paragen[iX,iY]
      x0= tab_x[iX-1,iY]
      x1= tab_x[iX,iY]
      x,F,OK= refine_xy(Xlis,Xfunc,F0,F1,x0,x1,Xtol)
      if OK:
        if F0>F1 : s= str(F1)+'='+str(F0)
        else     : s= str(F0)+'='+str(F1)
        if not s in lis_reac: lis_reac.append(s)
        val= s,x,y
        lis_xy_x.append(val)
        save_mixture_moles(x,y)
      else:
        if F>0:
        #-there is a paragen' Fm on [x0,x1] that is different from F0 & F1
        #--we need two second stage refinements
          xm=x
          Fm=F
          #--------------------second stage refinement between x0 and xm
          x,F,OK= refine_xy(Xlis,Xfunc,F0,Fm,x0,xm,Xtol)
          if OK:
            if F0>Fm: s= str(Fm)+'='+str(F0)
            else:     s= str(F0)+'='+str(Fm)
            if not s in lis_reac: lis_reac.append(s)
            val= s,x,y
            lis_xy_x.append(val)
            save_mixture_moles(x,y)
          else:
            val= iX,iY
            lis_false_x.append(val)
          #--------------------second stage refinement between xm and x1
          x,F,OK= refine_xy(Xlis,Xfunc,Fm,F1,xm,x1,Xtol)
          if OK:
            if F1>Fm: s= str(Fm)+'='+str(F1)
            else:     s= str(F1)+'='+str(Fm)
            if not s in lis_reac: lis_reac.append(s)
            val= s,x,y
            lis_xy_x.append(val)
            save_mixture_moles(x,y)
          else:
            val= iX,iY
            lis_false_x.append(val)
        else: #F<0
          val= iX,iY
          lis_false_y.append(val)
#--//      
#------refining the y-coordinate of the paragenesis change, at fixed x--
for iX in range(Xdim):
  x= tab_x[iX,0]
  include_modify(Xlis,Xfunc,x)
  for iY in range(1,Ydim):
    if tab_paragen[iX,iY-1] != tab_paragen[iX,iY]:
      F0= tab_paragen[iX,iY-1]
      F1= tab_paragen[iX,iY]
      y0= tab_y[iX,iY-1]
      y1= tab_y[iX,iY]
      y,F,OK= refine_xy(Ylis,Yfunc,F0,F1,y0,y1,Ytol)
      if OK:
        if F0>F1 : s= str(F1)+'='+str(F0)
        else     : s= str(F0)+'='+str(F1)
        if not s in lis_reac: lis_reac.append(s)
        val= s,x,y
        lis_xy_y.append(val)
        save_mixture_moles(x,y)
      else:
        if F>0:
        #-there is a paragen' Fm on [y0,y1] that is different from F0 & F1
        #--we need two second stage refinements
          ym=y
          Fm=F
          #--------------------second stage refinement between y0 and ym
          y,F,OK= refine_xy(Ylis,Yfunc,F0,Fm,y0,ym,Ytol)
          if OK:
            if F0>Fm: s= str(Fm)+'='+str(F0)
            else:     s= str(F0)+'='+str(Fm)
            if not s in lis_reac: lis_reac.append(s)
            val= s,x,y
            lis_xy_y.append(val)
            save_mixture_moles(x,y)
          else:
            val= iX,iY
            lis_false_y.append(val)
          #--------------------second stage refinement between ym and y1
          y,F,OK= refine_xy(Ylis,Yfunc,Fm,F1,ym,y1,Ytol)
          if OK:
            if F1>Fm: s= str(Fm)+'='+str(F1)
            else:     s= str(F1)+'='+str(Fm)
            if not s in lis_reac: lis_reac.append(s)
            val= s,x,y
            lis_xy_y.append(val)
            save_mixture_moles(x,y)
          else:
            val= iX,iY
            lis_false_y.append(val)
        else: #F<0
          val= iX,iY
          lis_false_y.append(val)
#--//          
#-------------------------------------------------------------//refining
#-----------build the lines for each reaction from the lists of points--
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
#--//
END= time.time()
print "TIME=", END - START
raw_input("type ENTER")

for i,paragen in enumerate(lis_paragen):
  print paragen
  
for reac in lis_reac:
  print reac
#raw_input()

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
      if len(w)>3: w=w[0:3]
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

#---------------------------------------------------mixture compositions
dim= len(s)
tMix= []
for x,y,s in tMixtur:
# convert tMixtur (table of x,y,string)
# to tMix (table of x,y,vector)
  v=[]
  for s in s.split():
    v.append(float(s))
  tMix.append((x,y,v))
tMix= sorted(tMix,key=lambda x: y)

f= open("mixture.tab",'w',0)
f.write("%s\t%s\t" % (Xlabel,Ylabel))
for s in mixture_names: f.write("%s\t" % s)
f.write('\n')
for x,y,mix in tMix:
  f.write("%.4g\t%.4g\t" % (x,y))
  for m in mix: f.write("%.4g\t" % (m))
  f.write('\n')
f.close()
#-------------------------------------------------//mixture compositions

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
symbols=['b','g','r','c','m','y']
symbols=['bo','go','ro','cs','mD','yd']
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
  fig.plot(vx, vy, symbols[i%len(symbols)], linestyle='-', linewidth=2.0)

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

sys.exit()

#-----------------------------------------------test mixture mean comp'n
def mix_is_same(v1,v2,tol):
  #return abs(v1[0]-v2[0])+abs(v1[1]-v2[1])+abs(v1[2]-v2[2])<tol
  return sum(abs(v1[:]-v2[:]))<tol*len(v1)

tMix= []
with open("mixture.tab",'r') as f:
  lines= f.readlines()
for il,line in enumerate(lines):
  ww= line.split()
  if il>0:
    x= ww[0]
    y= ww[1]
    v= plt.zeros(len(ww)-2,'float')
    for i,w in enumerate(ww):
      if i>1: v[i-2]= float(w)
    tMix.append((x,y,v))

f= open("mix_mean.tab",'w')

v0n= plt.zeros(len(ww)-2,'float')
v1n= plt.zeros(len(ww)-2,'float')
v2n= plt.zeros(len(ww)-2,'float')

# 3 poles = 3 paires, 1 triples
# 1-2 1-3
# 2-3
# 4 poles = 6 paires, 3 triples
# 1-2 1-3 1-4
# 2-3 2-4
# 3-4
# 5 poles = 10 paires, 10 triples
# 1-2 1-3 1-4 1-5
# 2-3 2-4 2-5
# 3-4 3-5
# 4-5
# 1-2-3 1-2-4 1-2-5 1-3-4 1-3-5 1-4-5
# 2-3-4 2-3-5 2-4-5
# 3-4-5


tol= 0.05
for x,y,v in tMix:
  nmix= 3
  p0= v[0]  ;  v0= v[1:4]   ;  print p0,v0 #
  p1= v[4]  ;  v1= v[5:8]   ;  print p1,v1 #
  p2= v[8]  ;  v2= v[9:12]  ;  print p2,v2 #
  v0n= v0
  v1n= v1
  v2n= v2
  if mix_is_same(v0,v1,tol): #--------------v0=v1
    nmix= nmix-1
    v0n= (p0*v0+p1*v1)/(p0+p1)
    p0= p0+p1
    p1= 0.
    if mix_is_same(v1,v2,tol): #------------v0=v1=v2
      nmix= nmix-1
      v0n= (p0*v0n+p2*v2)/(p0+p2)
      p0= p0+p2
      p2= 0.
  else:
    if mix_is_same(v1,v2,tol): #------------v1=v2!=v0
      nmix= nmix-1
      v1n= (p1*v1+p2*v2)/(p1+p2)
      p1= p1+p2
      p2= 0.
    else:
      if mix_is_same(v2,v0,tol): #----------v0=v2!=v1
        nmix= nmix-1
        v0n= (p0*v0+p2*v2)/(p0+p2)
        p0= p0+p2
        p2= 0.
  #
  print "n=",nmix
  if nmix>2: raw_input()
  f.write("%s\t%s\t" % (x,y))
  f.write("%d\t" % nmix)
  f.write('\n')

f.close()

sys.exit()
#---------------------------------------------//test mixture mean comp'n
