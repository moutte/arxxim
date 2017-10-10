s= "logk_cemdata07_phreeqc"
f= open(s+".new",'r')
#lines = [ l.strip() for l in instream if (l.strip() and l.strip()[0]!='#') ]
#print lines

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
  if "=" in ww: listPrim.append(ww[2])
  if ww[0]=="SECONDARY":
    break
    
for prim in listPrim: print prim
raw_input()    
#-------------------------------------------------//read PRIMARY species

#---------------------------------------------read SECONDARY aqu'species
listGamm= [] # open(s+".gam",'w')
listAnly= [] # open(s+".anl",'w')
listReac= []
listLogK= []
listSeco= []


for line in f:
  line= ChargeReplace(line)
  if not line.strip(): continue
  if line.strip()[0]=='#': continue
  ww= line.split()
  if '#' in ww:
    i= ww.index('#')
    ww= ww[0:i+1]
  if   ww[0]== "-gamma":                 listGamm.append(ww[1:-1])
  elif ww[0]== "-analytical_expression": listAnly.append(ww[1:-1])
  elif ww[0].lower()== "log_k":          listLogK.append(ww[1:-1])
  elif "=" in ww:
    j= ww.index("=")
    for p in listPrim:
      if p in line: line=line.replace(p, ' '+p)
    for p in listSeco:
      if p in line: line=line.replace(p, ' '+p)
    listReac.append(line.strip())
    listSeco.append(ww[j+1])
  if ww[0]=="PHASES":
    break

for reac in listReac:
  print reac
  #raw_input()
  
fo= open(s+"_aqu.tab",'w')
for prim in listPrim:
  fo.write("%s\n" % prim)
for i,reac in enumerate(listReac):
  ww= reac.split("=")
  fo.write("%s\t%s\t" % (ww[1],ww[0]))
  for w in listAnly[i]: fo.write("%s\t" % w)
  fo.write("\n")
fo.close()
#-------------------------------------------//read SECONDARY aqu'species

#------------------------------------------------------------read PHASES
listGamm= [] # open(s+".gam",'w')
listAnly= [] # open(s+".anl",'w')
listReac= []
listLogK= []
listName= []

for line in f:
  line= ChargeReplace(line)
  if not line.strip(): continue
  if line.strip()[0]=='#': continue
  
  ww= line.split()
  if '#' in ww:
    i= ww.index('#')
    ww= ww[0:i+1]
  if   ww[0]== "-gamma":                 listGamm.append(ww[1:-1])
  elif ww[0]== "-analytical_expression": listAnly.append(ww[1:-1])
  elif ww[0].lower()== "log_k":          listLogK.append(ww[1:-1])
  elif "=" in ww:
    listReac.append(line.strip())
  else:
    listName.append(line.strip())
  if ww[0]=="PHASES":
    break

for reac in listReac:
  print reac
  #raw_input()
  
fo= open(s+"_min.tab",'w')
for i,reac in enumerate(listReac):
  ww= reac.split("=")
  line= ww[1]
  for p in listPrim:
    if p in line: line=line.replace(p, ' '+p)
  for p in listSeco:
    if p in line: line=line.replace(p, ' '+p)
  fo.write("%s\t" % listName[i])
  fo.write("%s\t%s\t" % (line,ww[0]))
  for w in listAnly[i]: fo.write("%s\t" % w)
  fo.write("\n")
fo.close()

    
