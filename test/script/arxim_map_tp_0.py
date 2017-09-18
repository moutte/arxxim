import MyThermo as MT
import os, glob, sys
import pylab as plt
import numpy as np

#-------------------------------------------------------------------INIT
os.chdir("../")
#-----------------------------------------------------------user defined
Xmin,Xmax,Xdelta= 540., 560., 5.
Xtol= 1.
Ymin,Ymax,Ydelta= 4000., 5000.,50.
Ytol= 5.

Xmin,Xmax,Xdelta= 300., 900., 50.
Xtol= 5.
Ymin,Ymax,Ydelta= 500., 6500., 500.
Ytol= 10.

phases=["and","sil","kya"]

#---------------------------------------------------------//user defined
#------------------------------------------------------------------/INIT


#------------------------------------------------initialize the x,y grid
Xdim= int(round(abs(Xmax-Xmin)/Xdelta))+1
Xser= np.linspace(Xmin,Xmax,num=Xdim)

Ydim= int(round(abs(Ymax-Ymin)/Ydelta))+1
Yser= np.linspace(Ymin,Ymax,num=Ydim)

if 0:
  print Xser
  print Yser
  sys.exit()
#----------------------------------------------//initialize the x,y grid

#--------------------------------------------------------------------GEM
def phase_with_min_gibbs(TC,Pbar):
  TK= TC+273.15
  #
  G= []
  for fs in phases:
    x= MT.compute_gibbs(fs,TK,Pbar)
    G.append(x)
  #
  # index_Gmin= min(xrange(len(G)), key=G.__getitem__)
  #phase= phases[G.index(min(G))]
  return G.index(min(G))
  
#------------------------------------------------------------------//GEM

#-------------------------------------------------------refining routine
def refine_xy(selec,F0,F1,x0,x1,y,tolerance):
  OK= True
  F= -1
  #
  while abs(x0-x1)>tolerance:
    x= (x0+x1)/2.
    if selec=="x": F= phase_with_min_gibbs(x,y)
    if selec=="y": F= phase_with_min_gibbs(y,x)
    # if the paragen at x is the same as the paragen at x0,
    # then the paragen change occurs between x and x1, 
    # then x becomes the new x0
    if   F==F0: x0= x
    elif F==F1: x1= x
    else:
      print "F,F0,F1=", F,F0,F1
      OK= False
      break
  x= (x0+x1)/2.
  #
  return x,F,OK
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
  #---------------------------------------------------------------x-loop
  for iX,x in enumerate(Xser):
    print iY,iX
    j= phase_with_min_gibbs(x,y)
    if not j in lis_paragen: lis_paragen.append(j)
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
Xlis= "x"
Ylis= "y"
#---------------------------------------------------------------refining
#--refining the x-coordinate of the paragen change, at fixed y
for iY in range(Ydim):
  y= tab_y[0,iY]
  for iX in range(1,Xdim):
    if tab_phase[iX-1,iY] != tab_phase[iX,iY]:
      F0= tab_phase[iX-1,iY]
      F1= tab_phase[iX,iY]
      x0= tab_x[iX-1,iY]
      x1= tab_x[iX,iY]
      x,F,OK= refine_xy(Xlis,F0,F1,x0,x1,y,Xtol)
      if OK:
        if F0>F1 : s= str(F1)+'='+str(F0)
        else     : s= str(F0)+'='+str(F1)
        if not s in lis_reac: lis_reac.append(s)
        val= s,x,y
        lis_xy_x.append(val)
      else:
      #-there is a paragen' Fm on [x0,x1] that is different from F0 & F1
      #--we need two second stage refinements
        xm=x
        Fm=F
        #----------------------second stage refinement between x0 and xm
        x,F,OK= refine_xy(Xlis,F0,Fm,x0,xm,y,Xtol)
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
        x,F,OK= refine_xy(Xlis,Fm,F1,xm,x1,y,Xtol)
        if OK:
          if F1>Fm: s= str(Fm)+'='+str(F1)
          else:     s= str(F1)+'='+str(Fm)
          if not s in lis_reac: lis_reac.append(s)
          val= s,x,y
          lis_xy_x.append(val)
        else:
          val= iX,iY
          lis_false_x.append(val)
        
#--refining the y-coordinate of the paragen change, at fixed x
for iX in range(Xdim):
  x= tab_x[iX,0]
  for iY in range(1,Ydim):
    if tab_phase[iX,iY-1] != tab_phase[iX,iY]:
      F0= tab_phase[iX,iY-1]
      F1= tab_phase[iX,iY]
      y0= tab_y[iX,iY-1]
      y1= tab_y[iX,iY]
      y,F,OK= refine_xy(Ylis,F0,F1,y0,y1,x,Ytol)
      if OK:
        if F0>F1 : s= str(F1)+'='+str(F0)
        else     : s= str(F0)+'='+str(F1)
        if not s in lis_reac: lis_reac.append(s)
        val= s,x,y
        lis_xy_y.append(val)
      else:
      #-there is a paragen' Fm on [x0,x1] that is different from F0 & F1
      #--we need two second stage refinements
        ym=y
        Fm=F
        #----------------------second stage refinement between x0 and xm
        y,F,OK= refine_xy(Ylis,F0,Fm,y0,ym,x,Xtol)
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
        x,F,OK= refine_xy(Xlis,Fm,F1,ym,y1,x,Xtol)
        if OK:
          if F1>Fm: s= str(Fm)+'='+str(F1)
          else:     s= str(F1)+'='+str(Fm)
          if not s in lis_reac: lis_reac.append(s)
          val= s,x,y
          lis_xy_x.append(val)
        else:
          val= iX,iY
          lis_false_x.append(val)

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
  textstr= str(lis_paragen[i])
  textstr= phases[lis_paragen[i]]
  fig.text(x,y,textstr,horizontalalignment='center')
  
plt.savefig("0_arxim_map_tp_0"+".png")
plt.show()
#------------------------------------------------------//plot XY diagram
