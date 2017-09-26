#--------------------------------------------------------------mymakedir
import os
def myMkDir(newdir):
  """works the way a good mkdir should :)
    - already exists, silently complete
    - regular file in the way, raise an exception
    - parent directory(ies) does not exist, make them as well
  """
  if os.path.isdir(newdir):
    pass
  elif os.path.isfile(newdir):
    raise OSError("a file with the same name as the desired " \
                  "dir, '%s', already exists." % newdir)
  else:
    head, tail = os.path.split(newdir)
    if head and not os.path.isdir(head) :
      mymakedir(head)
    #print "_mkdir %s" % repr(newdir)
    if tail :
      os.mkdir(newdir)
#------------------------------------------------------------//mymakedir

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
#-------------------------------------------------------------------data

#----------------------------------------------------------plot_multiple
def plot_multiple(fig,xLog,yLog,data,labels,xMin,xMax,indexX,indexY,titre):
  # fig= subplot(1,1,1)
  fig.grid(color='r', linestyle='-', linewidth=0.2)
  fig.grid(True)
  
  symbols=['bo','go','ro','cs','mD','yd']
  # colors = ['cyan', 'lightblue', 'lightgreen', 'tan', 'pink','red', 'blue']

  # fig.set_xlabel(titX)
  # fig.set_ylabel(titY)
  fig.set_title(titre) #, fontsize=fontsize)
  # fig.text(min(vx),max(vy),textt,
  #  fontsize=16,ha = 'left', va = 'top')
  
  # fig.set_xlim(0,1.)
  # fig.set_xlim(1e-3,1e+3)
  #fig.semilogx(A[iMin:iMax,iX],A[iMin:iMax,iY],  'bo')
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
  raw_input("in progress ... type return ...")
  
