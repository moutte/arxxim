import os, glob, sys, time
import pylab as plt

import exceptions
def mynum(s):
  try:
    return float(s)
  except exceptions.ValueError:
    return 0.
    
#windows
sExe= "arxim.exe"
sExe= os.path.join("..","bin",sExe)

#linux
sExe= "arx_debug"
sExe= "arx_optim"
sExe= "a.out"
sExe= os.path.join("..","bin",sExe)

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

#-----------------------------------------------------------user-defined
fInn= "inn/map1a_fe_ox.inn"
Selection="FE"

Xlis= ["H"]
Ylis= ["OX"]

Xlabel= "pH"
Ylabel= "colog(f_O2(g))"

Xmin,Xmax,Xdelta,Xtol=  1., 13., 2., 0.02
Ymin,Ymax,Ydelta,Ytol= 40.,  4., 2., 0.1

#---------------------------------------------------------//user-defined

iKeyword= 0
iValue=   3

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
Xdim= int(round(abs(Xmax-Xmin)/Xdelta))+1
Xser= plt.linspace(Xmin,Xmax,num=Xdim)

Ydim= int(round(abs(Ymax-Ymin)/Ydelta))+1
Yser= plt.linspace(Ymin,Ymax,num=Ydim)

if False:
  print Xser
  print Yser
  sys.exit()
#----------------------------------------------//initialize the x,y grid

#--------------------------------------------------read chemical space
def readChemicalSpace():
  with open("tmp_species.tab") as f:
    s_names= f.readline().split()
    s_types= f.readline().split()
  for i,name in enumerate(s_names):
    print s_types[i]," ",name
  if DEBUG: raw_input()
  return s_names,s_types
#--------------------------------------------------//read chemical space

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
      F= -1
      break
    # print x0,x1
    # raw_input()
  if OK: x= (x0+x1)/2.
  #
  return x,F,OK
#-----------------------------------------------------//refining routine

#-----------------------------------------------map the dominant species
#----------------------------------------------------on the initial grid
tab_spc=plt.zeros((Xdim,Ydim),'int')
tab_y=plt.zeros((Xdim,Ydim),'float')
tab_x=plt.zeros((Xdim,Ydim),'float')

#--vars for refining function
lis_limits=  []
lis_xy=    []
lis_false_x= []
lis_false_y= []
#--
START= time.time()
#-----------------------------------------------------------------y-loop
for iY,y in enumerate(Yser):
  include_modify(Ylis,y) #-----------------modify the include file for y
  #---------------------------------------------------------------x-loop
  for iX,x in enumerate(Xser):
    #print iY,iX
    #
    include_modify(Xlis,x) #---------------modify the include file for x
    OK= arxim_execute(sArximCommand) #---------------------execute arxim
    #
    #-------------------------------------------------for first run only
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
          if "FE" in spc_names[i]:
            spc_include[i]= True
      for (i,name) in enumerate(spc_names):
        if spc_include[i]: print spc_names[i]
      if DEBUG: raw_input()
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
          x,F,OK= refine_xy(Xlis,F0,F1,x0,x1,Xtol)
          if OK:
            if F1>F0: s= (F0,F1)
            else:     s= (F1,F0)
            if not s in lis_limits: lis_limits.append(s)
            val= s,x,y
            #lis_xy_x.append(val)
            lis_xy.append(val)
          else:
            if F>0:
            #-there is a species Fm on [x0,x1] that is different from F0 & F1
            #--we need two second stage refinements
              xm=x
              Fm=F
              #------------------second stage refinement between x0 and xm
              x,F,OK= refine_xy(Xlis,F0,Fm,x0,xm,Xtol)
              if OK:
                if F1>Fm: s= (Fm,F1)
                else:     s= (F1,Fm)
                if not s in lis_limits: lis_limits.append(s)
                val= s,x,y
                lis_xy.append(val)
              else:
                val= iX,iY
                lis_false_x.append(val)
              #------------------second stage refinement between xm and x1
              x,F,OK= refine_xy(Xlis,Fm,F1,xm,x1,Xtol)
              if OK:
                if F1>Fm: s= (Fm,F1)
                else:     s= (F1,Fm)
                if not s in lis_limits: lis_limits.append(s)
                val= s,x,y
                lis_xy.append(val)
              else:
                val= iX,iY
                lis_false_x.append(val)
            else:
              val= iX,iY
              lis_false_x.append(val)
              continue
      
      #--refining the y-coordinate of a species change        
      if iY>0:
        if tab_spc[iX,iY-1] != tab_spc[iX,iY]:
          F0= tab_spc[iX,iY-1]
          F1= tab_spc[iX,iY]
          print F0,F1
          y0= tab_y[iX,iY-1]
          y1= tab_y[iX,iY]
          y,F,OK= refine_xy(Ylis,F0,F1,y0,y1,Ytol)
          if OK:
            if F1>F0: s= (F0,F1)
            else:     s= (F1,F1)
            if not s in lis_limits: lis_limits.append(s)
            val= s,x,y
            #lis_xy_y.append(val)
            lis_xy.append(val)
          else:
            if F>0:
            #-there is a paragen' Fm on [y0,y1] that is different from F0 & F1
            #--we need two second stage refinements
              ym=y
              Fm=F
              #------------------second stage refinement between y0 and ym
              y,F,OK= refine_xy(Ylis,F0,Fm,y0,ym,Ytol)
              if OK:
                if F1>Fm: s= (Fm,F1)
                else:     s= (F1,Fm)
                if not s in lis_limits: lis_limits.append(s)
                val= s,x,y
                lis_xy.append(val)
              else:
                val= iX,iY
                lis_false_y.append(val)
              #------------------second stage refinement between ym and y1
              y,F,OK= refine_xy(Xlis,Fm,F1,ym,y1,Ytol)
              if OK:
                #if F1>Fm: s= str(Fm)+'='+str(F1)
                #else:     s= str(F1)+'='+str(Fm)
                if F1>Fm: s= (Fm,F1)
                else:     s= (F1,Fm)
                if not s in lis_limits: lis_limits.append(s)
                val= s,x,y
                lis_xy_y.append(val)
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
for limit in lis_limits:
  points= []
  for val in lis_xy:
    if val[0]==limit:
      point= val[1],val[2]
      points.append(point)
  points= sorted(points,key=lambda x: x[0])
  lines.append(points)

print "LIST OF SPECIES :"
for i,species in enumerate(lis_spc):
  print i, species, spc_names[species]
print "LIST OF LIMITS :"
for limit in lis_limits:
  print limit
#raw_input()

END= time.time()
print "TIME=", END-START
if DEBUG: raw_input("type ENTER")

#sys.exit()
#-------------------------------------------------------------//refining

centroids=[ (0.,0.,0) for i in range(len(spc_names))]
for iY,y in enumerate(Yser):
  for iX,x in enumerate(Xser):
    i= tab_spc[iX,iY]
    x_,y_,n= centroids[i]
    x_= (n*x_+ x)/(n+1)
    y_= (n*y_+ y)/(n+1)
    n +=1
    centroids[i]= x_,y_,n
    
#---------------------------------------------the file name for the plot
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
plt.rcParams['figure.figsize']= 8.,6.   #ratio 4/3
plt.rcParams['figure.figsize']= 5.,3.75 #ratio 4/3
plt.rcParams['figure.figsize']= 6.,6.   #ratio 4/4
plt.rcParams.update({'font.size': 12})

fig= plt.subplot(1,1,1)
symbols=['bo','go','ro','cs','mD','yd','bo','go','ro','cs','mD','yd']
symbols=['b','g','r','c','m','y']
len_symb= len(symbols)
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
  fig.plot(vx, vy, symbols[i%len_symb], linestyle='-', linewidth=2.0)
  
for i,centroid in enumerate(centroids):
  x,y,n= centroid
  textstr= spc_names[i]
  fig.text(x,y,textstr,horizontalalignment='center')
  
plt.xlabel(Xlabel)
plt.ylabel(Ylabel)
  
plt.savefig(figName+".png")
plt.show()
#------------------------------------------------------//plot XY diagram

