import sys
import pylab as plt
import math as m

fData= "logk_cemdata07_phreeqc"
f= open(fData+".new",'r')
#lines = [ l.strip() for l in instream if (l.strip() and l.strip()[0]!='#') ]
#print lines

def logK(fit,T):
  # PHREEQC
  x= fit[0] + fit[1]*T + fit[2]/T + fit[3]*m.log(T,10) + fit[4]/T/T
  if len(fit)>5: x= x + fit[5]*T*T
  return x

#----------------------------------------------------------ChargeReplace
def ChargeReplace(w):
  w= w.replace(' + ', ' +')
  w= w.replace('++++++','+6')
  w= w.replace('+++++', '+5')
  w= w.replace('++++',  '+4')
  w= w.replace('+++',   '+3')
  w= w.replace('++',    '+2')
  
  w= w.replace(' - ', ' -')
  w= w.replace('------','-6')
  w= w.replace('-----', '-5')
  w= w.replace('----',  '-4')
  w= w.replace('---',   '-3')
  w= w.replace('--',    '-2')
  return w
#--------------------------------------------------------//ChargeReplace

#---------------------------------------------------read PRIMARY species
listPrim= []
for line in f:
  line= ChargeReplace(line)
  if not line.strip(): continue
  if line.strip()[0]=='#': continue
  ww= line.split()
  if ww[0]=="SECONDARY": break
  if "=" in ww: listPrim.append(ww[2])
    
for prim in listPrim: print prim
raw_input()    
#-------------------------------------------------//read PRIMARY species

#---------------------------------------------read SECONDARY aqu'species
listGamm= [] # open(s+".gam",'w')
listAnly= [] # open(s+".anl",'w')
listReac= []
listLogK= []
listSeco= []

def stoikio(s):
  ww= s.split()
  coeffs=[]
  specis=[]
  for i,w in enumerate(ww):
    if i%2==0: coeffs.append(w)
    else:      specis.append(w)
  return coeffs,specis

fi= open("logk_cemdata07_phreeqc_orig.new",'r')
for line in fi:
  line= line.strip()
  if not line.strip(): continue
  ww= line.split()
  if ww[0]=="PHASES": break
  if '=' in line: nam= line.split('=')[1]
  if ww[0]== "-analytical_expression": listAnly.append((nam,ww[1:]))

Xmin,Xmax,Xdim= 0.,150.,7
Tser= plt.linspace(Xmin,Xmax,num=Xdim)
fLogk= open(fData+".logk",'w')
for rec in listAnly:
  coefs=[]
  nam, vec= rec
  for w in vec:
    coefs.append(float(w))
  print nam,len(coefs),coefs
  logKs=[]
  for T in Tser:
    TK= T + 273.15
    logKs.append(-logK(coefs,TK)) # MINUS SIGN !!!
  fLogk.write("%s\t" % nam)
  for x in logKs: fLogk.write("%.4g\t" % x)
  fLogk.write('\n')
raw_input()

for line in f:
  line= line.strip()
  line= ChargeReplace(line)
  if not line.strip(): continue
  if line.strip()[0]=='#': continue
  ww= line.split()
  if ww[0]=="PHASES": break
  if '#' in ww:
    i= ww.index('#')
    ww= ww[0:i+1]
  if   ww[0]== "-gamma":                 listGamm.append(ww[1:])
  elif ww[0]== "-analytical_expression": continue #listAnly.append(ww[1:])
  elif ww[0].lower()== "log_k":          listLogK.append(ww[1:])
  elif "=" in ww:
    ww= line.split('=')
    listReac.append(ww[0])
    listSeco.append(ww[1])
    coeffs,specis= stoikio(ww[0])
    print ww[1], specis
raw_input()

for reac in listReac:
  print reac
  #raw_input()
  
fo= open(fData+"_aqu.tab",'w')
for prim in listPrim:
  fo.write("%s\n" % prim)
for i,reac in enumerate(listReac):
  fo.write("%s\t%s\t" % (listSeco[i],reac))
  # for w in listAnly[i]: fo.write("%s\t" % w)
  fo.write("\n")
fo.close()
#-------------------------------------------//read SECONDARY aqu'species

#------------------------------------------------------------read PHASES
listGamm= [] # open(s+".gam",'w')
listAnly= [] # open(s+".anl",'w')
listReac= []
listLogK= []
listName= []
listForm= []

for line in fi:
  line= line.strip()
  if not line.strip(): continue
  ww= line.split()
  #if ww[0]=="PHASES": break
  if '=' in line: nam= line.split('=')[0]
  if ww[0]== "-analytical_expression": listAnly.append((nam,ww[1:]))

for line in f:
  line= line.strip()
  if not line: continue
  if line[0]=='#': continue
  line= ChargeReplace(line)
  ww= line.split()
  if '#' in ww:
    i= ww.index('#')
    ww= ww[0:i+1]
  if   ww[0]== "-gamma":
    listGamm.append(ww[1:-1])
  elif ww[0]== "-analytical_expression": continue #listAnly.append(ww[1:-1])
  elif ww[0].lower()== "log_k":
    listLogK.append(ww[1:-1])
  elif "=" in ww:
    ww= line.split('=')
    listReac.append(ww[1])
    listForm.append(ww[0])
  else:
    listName.append(line.strip())
  #if ww[0]=="PHASES": break

for reac in listReac:
  print reac
raw_input()
  
fo= open(fData+"_min.tab",'w')
for i,reac in enumerate(listReac):
  fo.write("%s\t" % listName[i])
  fo.write("%s\t%s\t" % (listForm[i],reac))
  # for w in listAnly[i]: fo.write("%s\t" % w)
  fo.write("\n")
fo.close()

sys.exit()

if 0:
  fi= open("logk_cemdata07_phreeqc.orig",'r')
  fo= open("logk_cemdata07_phreeqc_orig.new",'w')
  for line in fi:
    if not line.strip(): continue
    if line.strip()[0]=='#': continue
    if ';' in line:
      ww= line.split(';')
      for w in ww: fo.write(w.strip()+'\n')
    else:
      fo.write(line.strip()+'\n')
    if 0:
      if ';' in line:
        ww= line.split(';')
        for w in ww:
          if "-analytical_expression" in w:
            xx= w.split()
            print len(xx),xx[1:] 

  sys.exit()





