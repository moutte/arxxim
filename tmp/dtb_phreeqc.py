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

Xmin,Xmax,Xdim= 0.,150.,7
Tser= plt.linspace(Xmin,Xmax,num=Xdim)

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
raw_input("ENT")    
#-------------------------------------------------//read PRIMARY species

#---------------------------------------------read SECONDARY aqu'species
# listGamm= [] # open(s+".gam",'w')
listAnly= [] # open(s+".anl",'w')
listReac= []
listLogK= []  ;  listNamK= []
listSeco= []
listName= []
listForm= []
listStok= []

def species_scan(s):
  ww= s.split()
  coeffs=[]
  specis=[]
  for i,w in enumerate(ww):
    if i%2==0: coeffs.append(float(w))
    else:      specis.append(w)
  return coeffs,specis

#------------read analytical_expression of sec'species and compute logks
fi= open("logk_cemdata07_phreeqc_orig.new",'r')
for line in fi:
  line= line.strip()
  if not line.strip(): continue
  ww= line.split()
  if ww[0]=="PHASES": break
  if '=' in line: nam= line.split('=')[1]
  if ww[0]== "-analytical_expression":
    listAnly.append((nam,ww[1:]))
    coefs=[float(w) for w in ww[1:]]
    logKs=[]
    for T in Tser:
      TK= T + 273.15
      logKs.append(-logK(coefs,TK)) # MINUS SIGN !!!
    listLogK.append(logKs)
    listNamK.append(nam)
#---------------------------------------------------------------------//

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
  if   ww[0]== "-gamma":                 continue #listGamm.append(ww[1:])
  elif ww[0]== "-analytical_expression": continue #listAnly.append(ww[1:])
  elif ww[0].lower()== "log_k":          continue #listLogK.append(ww[1:])
  elif "=" in ww:
    ww= line.split('=')
    ww= [w.strip() for w in ww]
    listReac.append(ww[0])
    listSeco.append(ww[1])
    listName.append(ww[1])
    listForm.append(ww[1])
#raw_input()

if 0:
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
#------------read analytical_expression of sec'species and compute logks
for line in fi:
  line= line.strip()
  if not line.strip(): continue
  ww= line.split()
  #if ww[0]=="PHASES": break
  if '=' in line: nam= line.split('=')[0]
  if ww[0]== "-analytical_expression":
    listAnly.append((nam,ww[1:]))
    coefs=[]
    for w in ww[1:]: coefs.append(float(w))
    logKs=[]
    for T in Tser:
      TK= T + 273.15
      logKs.append(-logK(coefs,TK)) # MINUS SIGN !!!
    listLogK.append(logKs)
    listNamK.append(nam)
#---------------------------------------------------------------------//

for line in f:
  line= line.strip()
  if not line: continue
  if line[0]=='#': continue
  line= ChargeReplace(line)
  ww= line.split()
  if '#' in ww:
    i= ww.index('#')
    ww= ww[0:i+1]
  if   ww[0]== "-gamma": continue # listGamm.append(ww[1:-1])
  elif ww[0]== "-analytical_expression": continue # listAnly.append(ww[1:-1])
  elif ww[0].lower()== "log_k": continue # listLogK.append(ww[1:-1])
  elif "=" in ww:
    ww= line.split('=')
    listReac.append(ww[1])
    listForm.append(ww[0])
  else:
    listName.append(line.strip())
  #if ww[0]=="PHASES": break

#---------------stoikio of sec'species with respect to prim'species only
print "==sec'species stoikio=="
for i,reac in enumerate(listReac):
  coeffs,specis= species_scan(reac)
  print listName[i],coeffs, specis
  stoikio= [0. for sp in listPrim]
  for j,sp in enumerate(specis):
    if sp in listPrim:
      stoikio[listPrim.index(sp)]= coeffs[j]
  for j,sp in enumerate(specis):
    if sp in listSeco:
      isec= listSeco.index(sp)
      if isec<i:
        for k in range(len(stoikio)):
          kk= listStok[isec][k]
          if kk!=0:
            stoikio[k]= stoikio[k] + listStok[isec][k]*coeffs[j]
            for T in range(len(Tser)):
              listLogK[i][T]= listLogK[i][T] + kk*listLogK[isec][T]
      else:
        print sp, "=species without stoik"
    elif not sp in listPrim:
      print sp, "=species not found"
  print stoikio
  listStok.append(stoikio)
  # raw_input()
raw_input("ENT")
#---------------------------------------------------------------------//

fLogk= open(fData+".logk",'w')
for i,logk in enumerate(listLogK):
  fLogk.write("%s\t" % listNamK[i])
  for x in logk: fLogk.write("%.4g\t" % x)
  fLogk.write('\n')
raw_input("ENT")

#-----------------------------------------------------write stokio table
fo= open(fData+"_stoikio.tab",'w')
fo.write(".formula\t")
for prim in listPrim: fo.write("%s\t" % prim)
fo.write("\n")
for i,frm in enumerate(listForm):
  fo.write("%s\t" % frm)
  for x in listStok[i]: fo.write("%.4g\t" % x)
  fo.write("\n")
fo.close()
#---------------------------------------------------------------------//
  
if 0:
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

if 1:
  fi= open("logk_cemdata07_phreeqc.orig",'r')
  fo= open("logk_cemdata07_phreeqc_orig.new",'w')
  for line in fi:
    line= line.strip()
    if not line: continue
    if line[0]=='#': continue
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





