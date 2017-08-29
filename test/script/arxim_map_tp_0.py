import MyThermo as MT
import os, glob, sys
import pylab as plt
import numpy as np

sExe= "arxim.exe"    #windows
sExe= "arxim"        #linux
sExe= os.path.join("..","bin",sExe)
sDebug= "1"
sCmd=  "GEM"

#-------------------------------------------------------------------INIT
os.chdir("../")
#----------------------------------------------------cleaning tmp_ files
for l in glob.glob("tmp_*"): os.remove(l)
if os.path.isfile("error.log"): os.remove("error.log")
#sys.exit()  
#---------------------------------------------------/cleaning tmp_ files

#------------------------------------------------------------input files
fInn= "inn/f1r_sialh.inn"
fInn= "inn/f1q_sial.inn"

'''
the include block to be modified:
SYSTEM.GEM
  TDGC  1200
  PBAR  400
  SIO2   SiO2  1.0
  AL2O3  Al2O3 1.0
END
Keyword is TDGC / PBAR,  its index is 0  -> iKeyword= 0
Value is 1200 / 400,     its index is 1  -> iValue=   1
'''

iKeyword= 0
iValue=   1

lisX= ["TDGC"]
lisY= ["PBAR"]

Xmin,Xmax,Xdelta= 540., 560., 5.
tol_x= 1.
Ymin,Ymax,Ydelta= 4000., 5000.,50.
tol_y= 5.

Xmin,Xmax,Xdelta= 300., 900., 50.
tol_x= 5.
Ymin,Ymax,Ydelta= 500., 6500., 500.
tol_y= 10.

#----------------------------------------------------------//input files

fInclude= fInn.replace(".inn",".include")
#--clean the include file
with open(fInclude,'r') as f: lines = f.readlines()
f=open(fInclude,'w')
for line in lines:
  ll= line.strip()
  if ll=='':     continue
  if ll[0]=='!': continue
  f.write(line)
f.close()
#--

sArximCommand= sExe,fInn,sDebug,sCmd
# the arxim executable (sExe) accepts 3 arguments, successively:
# > fInn: the path of the arxim script
# > sDebug: an integer that controls the screen output at runtime
# for the present case, sDebug=1, which restricts the screen output
# > sCmd: the calcualtion command, e.g. SPC, GEM, EQU, ...
fResult= "tmp_gem.tab"
#------------------------------------------------------------------/INIT


#------------------------------------------------initialize the x,y grid
Xdim= int(abs(Xmax-Xmin)/Xdelta)+1
Xser= np.linspace(Xmin,Xmax,num=Xdim)

Ydim= int(abs(Ymax-Ymin)/Ydelta)+1
Yser= np.linspace(Ymin,Ymax,num=Ydim)

if 0:
  print Xser
  print Yser
  sys.exit()
#----------------------------------------------//initialize the x,y grid

#-------------------------------------------------------refining routine
def refine_xy(lis,F0,F1,x0,x1,tolerance):
  OK= True
  #
  while abs(x0-x1)>tolerance:
    x= (x0+x1) /2.
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
        print "NEW:",paragen
        OK= False
        break
    # print x0,x1
    # raw_input()
  x= (x0+x1)/2.
  #
  return x,OK
#-----------------------------------------------------//refining routine

#----------------------------------------------map the stable assemblage
#----------------------------------------------------on the initial grid
tab_phase=plt.zeros((Xdim,Ydim),'int')
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
  include_modify(lisY,y) #-----------------modify the include file for y
  #---------------------------------------------------------------x-loop
  for iX,x in enumerate(Xser):
    print iY,iX
    #
    TK=   x+273.15
    Pbar= y
    G_and= MT.compute_gibbs("and",TK,Pbar)
    G_sil= MT.compute_gibbs("sil",TK,Pbar)
    G_kya= MT.compute_gibbs("kya",TK,Pbar)
    G= [G_and,G_sil,G_kya]
    
    Gmin= min(G)
    # index_Gmin= min(xrange(len(G)), key=G.__getitem__)
    index_Gmin= G.index(min(G))
      paragen= arxim_result(fResult) #-------------------read arxim result
      print paragen
      # raw_input()
      if not paragen in lis_paragen:
        lis_paragen.append(paragen)
        #centroids.append((x,y,0))
      j= lis_paragen.index(paragen)
      tab_phase[iX,iY]= j
      tab_x[iX,iY]= x
      tab_y[iX,iY]= y
  #-------------------------------------------------------------//x-loop
#---------------------------------------------------------------//y-loop

#---------------------------------compute the positions for field labels
centroids=[ (0.,0.,0) for i in range(len(lis_paragen))]
for iY,y in enumerate(Yser):
  for iX,x in enumerate(Xser):
    i= tab_phase[iX,iY]
    x_,y_,n= centroids[i]
    x_= (n*x_+ x)/(n+1)
    y_= (n*y_+ y)/(n+1)
    n +=1
    centroids[i]= x_,y_,n
if False:
  for centroid in centroids:
    print centroid
#-------------------------------//compute the positions for field labels
    
print tab_phase
#--------------------------------------------//map the stable assemblage
#sys.exit()

#---------------------------------------------------------------refining
#--refining the x-coordinate of the paragen change, at fixed y
for iY in range(Ydim):
  y= tab_y[0,iY]
  include_modify(lisY,y)
  for iX in range(1,Xdim):
    if tab_phase[iX-1,iY] != tab_phase[iX,iY]:
      F0= tab_phase[iX-1,iY]
      F1= tab_phase[iX,iY]
      x0= tab_x[iX-1,iY]
      x1= tab_x[iX,iY]
      x,OK= refine_xy(lisX,F0,F1,x0,x1,tol_x)
      if OK:
        if F0>F1 : s= str(F1)+'='+str(F0)
        else     : s= str(F0)+'='+str(F1)
        if not s in lis_reac: lis_reac.append(s)
        val= s,x,y
        lis_xy_x.append(val)
      else:
        val= iX,iY
        lis_false_x.append(val)
        
#--refining the y-coordinate of the paragen change, at fixed x
for iX in range(Xdim):
  x= tab_x[iX,0]
  include_modify(lisX,x)
  for iY in range(1,Ydim):
    if tab_phase[iX,iY-1] != tab_phase[iX,iY]:
      F0= tab_phase[iX,iY-1]
      F1= tab_phase[iX,iY]
      y0= tab_y[iX,iY-1]
      y1= tab_y[iX,iY]
      y,OK= refine_xy(lisY,F0,F1,y0,y1,tol_y)
      if OK:
        if F0>F1 : s= str(F1)+'='+str(F0)
        else     : s= str(F0)+'='+str(F1)
        if not s in lis_reac: lis_reac.append(s)
        val= s,x,y
        lis_xy_y.append(val)
      else:
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
#---------------------------------------//processing the paragenese list

#--------------------------------------------------------plot XY diagram
plt.rcParams['figure.figsize']= 8,6
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
  fig.text(x,y,textstr,horizontalalignment='center')
  
plt.savefig("0_arxim_map_tp"+".png")
plt.show()
#------------------------------------------------------//plot XY diagram
