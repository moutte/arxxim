#----------------------------------------------------------------myMkDir
import os
def myMkDir(newdir):
  if os.path.isdir(newdir):
    pass
  elif os.path.isfile(newdir):
    raise OSError("a file with the same name '%s' already exists." % newdir)
  else:
    head, tail = os.path.split(newdir)
    if head and not os.path.isdir(head): myMkDir(head)
    if tail: os.mkdir(newdir)
#--------------------------------------------------------------//myMkDir

def extractFileName(s):
  head, tail = os.path.split(s)
  if '.' in tail:
    name_= tail.split('.')[0]
  else:
    name_= tail
  return name_
  
#-------------------------------------------------------------------data
import exceptions
def num(s):
  try:
    return float(s)
  except exceptions.ValueError:
    return 0.
    
def TextToLines(InName,vL):
  f= open(InName,'r')
  for line in f:
    if line.strip() =="": break
    if line[0] != '!': #skip comment lines
      vL.append(line)
  f.close()
  #return vL

def vLines2Table(vL,nL,nC,TT):
  # conversion of a list of lines, with numeric values only,
  # to a table of float,
  # skipping the first line, which contains the labels
  for i,line in enumerate(vL):
    if i==0: continue #skip first line
    ww= line.split('\t')
    for j in range(nC):
      print ww[j]
      TT[i-1,j]= num(ww[j])
      #raw_input()
  return TT

import pylab as plt
def lines2table(lines):
  lines= [x for x in lines if len(x.strip())>0] #remove empty lines
  labels= lines[0].split()
  TT= plt.zeros((len(lines)-1,len(labels)),'float')
  for i,line in enumerate(lines):
    if i==0: continue #skip first line
    ww= line.split()
    for j in range(len(ww)): TT[i-1,j]= num(ww[j])
  return labels,TT
  
def table_sort(labels,TT):
  nlin,ncol= TT.shape
  T= []
  i0= labels.index("H2O")
  for i in range(ncol):
    if i>i0:
      t= (labels[i],TT[:,i],max(TT[:,i]))
      T.append(t)
  return T
  # return a list of tuples (label, array of floats, max)

def table_select(iMode,sKey,labels,TT):
  nlin,ncol= TT.shape
  T= []
  if iMode==1:
    i0= labels.index(sKey)
    for i in range(ncol):
      if i>i0:
        t= (labels[i],TT[:,i],max(TT[:,i]))
        T.append(t)
  if iMode==2:
    for i in range(ncol):
      if sKey in labels[i]:
        lb= labels[i].replace(sKey,"")
        t= (lb,TT[:,i],max(TT[:,i]))
        T.append(t)
  return T
  # return a list of tuples (label, array of floats, max)
#-------------------------------------------------------------------data

#-------------------------------------------------------------------plot
def plot(fig,titr,dataX,dataY,xLog,yLog):
  
  fig.grid(color='r', linestyle='-', linewidth=0.2)
  fig.grid(True)
  
  
  symbols=['bo','go','ro','co','mo','yo',\
           'bd','gd','rd','cd','md','yd',\
           'bD','gD','rD','cD','mD','yD']
  colors= ['cyan', 'lightblue', 'lightgreen', 'tan', 'pink','red', 'blue']

  labX,vX=dataX 
  for i,TT in enumerate(dataY):
    lb,vY,maxx= TT
    sy= symbols[(i)%len(symbols)]
    if xLog and yLog:
      fig.loglog(vX, vY, sy, linestyle='-',linewidth=1.0,label=lb)
      xmin,xmax= fig.get_xlim()
      ymin,ymax= fig.get_ylim()
      fig.set_xlim(xmax/1.e4,xmax)
      fig.set_ylim(ymax/1.e4,ymax)
    else:
      if xLog:
        fig.semilogx(vX, vY, sy, linestyle='-',linewidth=1.0,label=lb)
      elif yLog:
        fig.semilogy(vX, vY, sy, linestyle='-', linewidth=1.0,label=lb)
      else:
        fig.plot(vX, vY, sy, linestyle='-', linewidth=1.0,label=lb)
        #ymin,ymax= fig.get_ylim()
        #fig.set_ylim(ymax/1.e4,ymax)
        #fig.set_ylim(-15.,-1.)
      
  fig.set_title(titr) #, fontsize=fontsize)
  #-legend= fig.legend(loc='upper left')
  #-legend= fig.legend(bbox_to_anchor=(1.1, 1.05))
  legend= fig.legend(loc='lower center', bbox_to_anchor=(0.5,0.),
          ncol=3, fancybox=True) #, shadow=True)
  for lb in legend.get_texts(): lb.set_fontsize('small')
#-----------------------------------------------------------------//plot

#--scan the inn file -> list of words
def inn_scan_list(sName,sBlock,sKey):
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
def inn_scan_word(sName,sBlock,sKey):
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

#----------------------------------------------------------plot_multiple
def plot_multiple(fig,xLog,yLog,data,labels,xMin,xMax,indexX,indexY,titre):
  # fig= subplot(1,1,1)
  fig.grid(color='r', linestyle='-', linewidth=0.2)
  fig.grid(True)
  
  symbols=['bo','go','ro','co','mo','yo','bd','gd','rd','cd','md','yd']
  # colors = ['cyan', 'lightblue', 'lightgreen', 'tan', 'pink','red', 'blue']

  # fig.set_xlabel(titX)
  # fig.set_ylabel(titY)
  fig.set_title(titre) #, fontsize=fontsize)
  # fig.text(min(vx),max(vy),textt,
  #  fontsize=16,ha = 'left', va = 'top')
  
  # fig.set_xlim(0,1.)
  # fig.set_xlim(1e-3,1e+3)
  # fig.semilogx(A[iMin:iMax,iX],A[iMin:iMax,iY],  'bo')
  for i in range(len(indexY)):
    if i<len(symbols):
      # vx=A[:,indexX]
      # vy=A[:,indexY[i]]
      vx= []
      vy= []
      for j in range(len(data)):
        x= data[j,indexX]
        y= data[j,indexY[i]]
        if x>xMin and x<xMax:
          vx.append(x)
          vy.append(y)
      if xLog and yLog:
        fig.loglog(vx, vy, symbols[i], 
          linestyle='-', linewidth=1.0, label=labels[indexY[i]])
      else:
        if xLog:
          fig.semilogx(vx, vy, symbols[i], 
            linestyle='-', linewidth=1.0, label=labels[indexY[i]])
        elif yLog:
          fig.semilogy(vx, vy, symbols[i], 
            linestyle='-', linewidth=1.0, label=labels[indexY[i]])
        else:
          fig.plot(vx, vy, symbols[i], 
            linestyle='-', linewidth=1.0, label=labels[indexY[i]])
      
    fig.legend(loc='upper left')

#--------------------------------------------------------//plot_multiple

#------------------------------------------------------------plot_single
def plot_single_1(fig,xLog,yLog,vx,vy,vlabs,titX,titY,titre):
  # fig= subplot(1,1,1)
  fig.grid(color='r', linestyle='-', linewidth=0.2)
  fig.grid(True)

  fig.set_xlabel(titX)
  fig.set_ylabel(titY)
  fig.set_title(titre) #, fontsize=fontsize)
  # fig.text(min(vx),max(vy),textt,
  #  fontsize=16,ha = 'left', va = 'top')
  
  # fig.set_xlim(0,1.)
  # fig.set_xlim(1e-3,1e+3)
  #fig.semilogx(A[iMin:iMax,iX],A[iMin:iMax,iY],  'bo')
  
  if xLog and yLog:
    fig.loglog(vx, vy, 'bo')
  else:
    if xLog:
      fig.semilogx(vx, vy, 'bo')
    elif yLog:
      fig.semilogy(vx, vy, 'bo')
    else:
      fig.plot(vx, vy, 'bo')
  
  if len(vlabs)>0:
    for label, x, y in zip(vlabs, vx, vy):
      plt.annotate(
        label, 
        xy = (x, y),
        xytext = (-2,2),
        textcoords = 'offset points', 
        ha = 'right', 
        va = 'bottom')
        #bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
        #arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
      
#----------------------------------------------------------//plot_single

if __name__ == '__main__':
  # raw_input("in progress ... type return ...")
  #
  s="MyLib.restab"
  lines= open(s,'r').readlines()
  labels,tData= lines2table(lines)
  #
  if "Time/YEAR" in labels:
    iX=   labels.index("Time/YEAR")
    labX= "Time/YEAR"
  else:
    iX= 0
    labX= "x"
  dataX= (labX,tData[:,iX])
  dataY= table_select(2,"PhiM_",labels,tData)
  #
  fig= plt.subplot()
  plot(fig,dataX,dataY,False,False)
  plt.show()
  fig= plt.subplot()
  plot(fig,dataX,dataY,True,True)
  plt.show()
  
