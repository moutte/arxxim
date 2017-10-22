import sys

lisTdg= [0.,25.,60.,100.,150.,200.,250.,300.]
lisMaj=['H2O','H+','O2(aq)','H4(SiO4)','Al+3','Fe+2','Mg+2','Ca+2','Na+','K+']
lisAni=['CO3-2','Cl-','F-','H2(PO4)-','SO4-2','B(OH)4-','NO3-']
lisTrc=['Mn+2','CrO4-2','Ni+2','Cu+2','Zn+2','MoO4-2','AsO4-3','Sr+2','Ba+2','Pb+2']

fname= "THC_crunch_v9b0"
sSource="THCv9b0"

#---------------------------------------------------read primary species
lisPrimary= []
lisSize= []
fi=open(fname+".primary","r")
if 0: fo=open("lisprimary","w")
lisPrimary=[]
for ll in fi:
  if ll.strip():
    w= ll.split()[0].replace("'","")
    lisPrimary.append(w)
    lisSize.append(float(ll.split()[1]))
    if 0: fo.write('%s\n' % w)
fi.close()
if 0: fo.close()
  
if 0: sys.exit()
#-------------------------------------------------//read primary species
lisSelect= lisMaj + lisAni + lisTrc
nPrim= len(lisPrimary)

#Coeffs= ["0.0" for x in range(nPrim)]
#print Coeffs
#sys.exit()
#-----------------------------------------------------------species_scan
def species_scan(n,ww):
  coeffs=[]
  specis=[]
  for i,w in enumerate(ww):
    if i%2==0: coeffs.append(float(w))
    else:      specis.append(w.replace("'",""))
  return coeffs,specis
#---------------------------------------------------------//species_scan

#-----------------------------------------------select secondary species
if 0:
  fi= open(fname+".second","r")
  fs= open(fname+".second.select","w")
  fr= open(fname+".second.reject","w")
  lisSecond= []
  for ll in fi:
    if not ll.strip(): continue
    ww= ll.split()
    namSp= ww[0].replace("'","")  ;  print namSp
    nReac= int(ww[1]) # number of prim'species in formation reaction
    sReac= ww[nReac:nReac+2*nReac]
    isValid= True
    for i in range(nReac):
      sp= sReac[2*i+1].replace("'","")
      if not sp in lisSelect + lisSecond: isValid= False
    if isValid:
      lisSecond.append(namSp)
      fs.write(ll)
    else:
      fr.write(ll)
  fi.close()
  fs.close()
  fr.close()
  sys.exit()
#---------------------------------------------//select secondary species

lisSize= [lisSize[lisPrimary.index(sp)] for sp in lisSelect]
lisPrimary= [sp for sp in lisPrimary if sp in lisSelect]
lisPrimary= lisSelect

#-----------------------------------------------------species to element
fi= open(fname+".primary.stoikio",'r')
ll= fi.readline()
ww= ll.split()
lisElem= ww[1:]
print lisElem

lisTmpStoikio= []
lisTmpPrimary= []
for ll in fi:
  if ll.strip():
    ww= ll.split()
    lisTmpPrimary.append(ww[0])
    stoik= [float(x) for x in ww[1:]]
    lisTmpStoikio.append(stoik)
lisStoikioPrim= []
for sp in lisPrimary:
  stoik= lisTmpStoikio[lisTmpPrimary.index(sp)]
  lisStoikioPrim.append(stoik)
  print sp, stoik
  
fi.close()
#sys.exit()
#---------------------------------------------------//species to element

#---------------------------------------------------build ECFORM formula
def buildEcform(stoikio):
  ecform= [0. for el in lisElem]
  for i,x in enumerate(stoikio):
    if x!=0.:
      for j in range(len(lisElem)):
        ecform[j]= ecform[j] + x*lisStoikioPrim[i][j]
  needDiv= False
  for x in ecform:
     if round(x)!=x: needDiv= True ;  break
  if needDiv: ecform= [x*100. for x in ecform]
  s= ""
  for i,x in enumerate(ecform):
    if abs(x)>0.001:
      s= s+ lisElem[i] + '(' + str(x) + ')'
  if needDiv: s= s + "Div(100)"
  s=s.replace("Chg(-","-(")
  s=s.replace("Chg(","+(")
  sEcform= s.replace(".0)",")")
  print sEcform
  return sEcform
#-------------------------------------------------//build ECFORM formula
  
#-------------------------------------------------read secondary species
fi= open(fname+".second.select","r")
ft= open("logk_"+fname+"_aqu_stoik.tab","w")
fo= open("logk_"+fname+"_aqu.tab","w")
fo25= open("logk_"+fname+"_aqu_25c.tab","w")
fw= open("tmp","w")

fo.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t" % 
        ("TYPE","SOURCE","INDEX","NAME","ECFORM","SIZE","PARAMETERS"))
#for t in lisTdg: fo.write("LOGK\t")
fo.write('\n')
fo25.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t" % 
        ("TYPE","SOURCE","INDEX","NAME","ECFORM","SIZE","PARAMETERS"))
#for t in lisTdg: fo.write("LOGK\t")
fo25.write('\n')

ft.write("%s\t%s\t" % ("TYP","SPECIES"))
for sp in lisPrimary: ft.write("%s\t" % sp)
ft.write('\n')

lisSecond= []
lisStoik=  []
lisLogks=  []

sTyp= "AQU"
iIdx= 0

for i,sNam in enumerate(lisPrimary):
  iIdx+=1
  sIdx='_'+str(iIdx).zfill(4) 
  fo.write('%s\t%s\t%s\t%s\t' % (sTyp,sSource,sIdx,sNam))
  stoikio= [0. for p in lisPrimary]  ;  stoikio[i]=1
  ss= buildEcform(stoikio)
  #
  fo.write('%s\t' % ss)
  fo.write('%.4g\t' % lisSize[i])
  for x in lisTdg: fo.write('0.0\t')
  fo.write('\n')
  #
  fo25.write('%s\t' % ss)
  fo25.write('%.4g\t' % lisSize[i])
  fo25.write('0.0\t')
  fo25.write('\n')
  
for ll in fi:
  if not ll.strip(): continue
  ww= ll.split()
  
  namSp= ww[0].replace("'","")  ;  print namSp
  lisSecond.append(namSp)
  
  nReac= int(ww[1]) # number of prim'species in formation reaction
  sReac= ww[2:2+2*nReac]
  sLogk= ww[2+2*nReac:2+2*nReac+len(lisTdg)]
  size= float(ww[2+2*nReac+len(lisTdg)])
  # print namSp, sReac
  # raw_input()
  
  stoikio= [0. for s in lisPrimary]
  logks= [float(w) for w in sLogk]
  coeffs,specis= species_scan(nReac,sReac)
  
  for j,sp in enumerate(specis):
    if sp in lisPrimary:
      stoikio[lisPrimary.index(sp)]= coeffs[j]
      
  for j,sp in enumerate(specis):
    if sp in lisSecond:
      fw.write("%s %s\n" % (namSp,sp))
      isec= lisSecond.index(sp)
      for k in range(len(stoikio)):
        stoikio[k]= stoikio[k] + lisStoik[isec][k]*coeffs[j]
    elif not sp in lisPrimary:
      print sp, "=species not found"
      raw_input()
      
  for j,sp in enumerate(specis):
    if sp in lisSecond:
      isec= lisSecond.index(sp)
      for T in range(len(lisTdg)):
        if lisLogks[isec][T]!=500. and logks[T]!=500:
          logks[T]= logks[T] + lisLogks[isec][T]*coeffs[j]

  lisStoik.append(stoikio)
  lisLogks.append(logks)
  
  sEcform= buildEcform(stoikio)
        
  # write stoikio's
  ft.write('%s\t%s\t' % (sTyp,namSp))
  for x in stoikio: ft.write('%.4g\t' % x)
  ft.write('\n')
  # write logk's
  iIdx+=1
  sIdx='_'+str(iIdx).zfill(4)
  #
  fo.write('%s\t%s\t%s\t%s\t' % (sTyp,sSource,sIdx,namSp))
  fo.write('%s\t' % sEcform)
  fo.write('%.4g\t' % size)
  for x in logks: fo.write('%.4g\t' % -x)
  fo.write('\n')
  #
  fo25.write('%s\t%s\t%s\t%s\t' % (sTyp,sSource,sIdx,namSp))
  fo25.write('%s\t' % sEcform)
  fo25.write('%.4g\t' % size)
  fo25.write('%.4g\t' % -logks[1])
  fo25.write('\n')

fi.close()
fo.close()
fo25.close()
fw.close()

#sys.exit()
#---------------------------------------------end read secondary species

#--------------------------------------------------------select minerals
if 0:
  fi= open(fname+".minerals","r")
  fs= open(fname+".minerals.select","w")
  fr= open(fname+".minerals.reject","w")
  for ll in fi:
    if not ll.strip(): continue
    if ll.strip()[0]=='!': continue
    ww= ll.split()
    nReac= int(ww[2]) # number of prim'species in formation reaction
    sReac= ww[3:3+2*nReac]
    coeffs,specis= species_scan(nReac,sReac)
    isValid= True
    for sp in specis:
      if not sp in lisSelect + lisSecond: isValid= False
    if isValid:
      print specis
      fs.write(ll)
    else:
      fr.write(ll)
  fi.close()
  fs.close()
  fr.close()
  sys.exit()
#------------------------------------------------------//select minerals

#----------------------------------------------------------read minerals
fi= open(fname+".minerals.select","r")
fo= open("logk_"+fname+"_min.tab","w")
fo25= open("logk_"+fname+"_min_25c.tab","w")
ft= open("logk_"+fname+"_min_stoik.tab","w")
fw= open("tmp","w")

fo.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t" % 
        ("TYPE","SOURCE","INDEX","NAME","ECFORM","SIZE","PARAMETERS"))
fo.write('\n')

fo25.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t" % 
        ("TYPE","SOURCE","INDEX","NAME","ECFORM","SIZE","PARAMETERS"))
fo25.write('\n')

ft.write("%s\t%s\t" % ("TYP","SPECIES"))
for sp in lisPrimary: ft.write("%s\t" % sp)
ft.write('\n')

lisMin=   []

for ll in fi:
  if not ll.strip(): continue
  if ll.strip()[0]=='!': continue
  ww= ll.split()
  
  namSp= ww[0].replace("'","")  ;  print namSp
  lisMin.append(namSp)
  
  if "(g)" in namSp: sTyp= "GAS"
  else:              sTyp= "MIN"
  
  nReac= int(ww[2]) # number of prim'species in formation reaction
  sReac= ww[3:3+2*nReac]
  sLogk= ww[3+2*nReac:3+2*nReac+len(lisTdg)]
  size= float(ww[1])
  
  stoikio= [0. for s in lisPrimary]
  logks= [float(w) for w in sLogk]
  coeffs,specis= species_scan(nReac,sReac)
  print coeffs
  print specis
  
  for j,sp in enumerate(specis):
    if sp in lisPrimary:
      stoikio[lisPrimary.index(sp)]= coeffs[j]
      
  for j,sp in enumerate(specis):
    if sp in lisSecond:
      fw.write("%s %s\n" % (namSp,sp))
      isec= lisSecond.index(sp)
      for k in range(len(stoikio)):
        stoikio[k]= stoikio[k] + lisStoik[isec][k]*coeffs[j]
    elif not sp in lisPrimary:
      print sp, "=species not found"
      raw_input()
      
  for j,sp in enumerate(specis):
    if sp in lisSecond:
      isec= lisSecond.index(sp)
      for T in range(len(lisTdg)):
        if lisLogks[isec][T]!=500. and logks[T]!=500:
          logks[T]= logks[T] + lisLogks[isec][T]*coeffs[j]
  
  lisStoik.append(stoikio)
  lisLogks.append(logks)
  sEcform= buildEcform(stoikio)
  
  # write stoikio's
  ft.write('%s\t%s\t' % (sTyp,namSp))
  for x in stoikio: ft.write('%.4g\t' % x)
  ft.write('\n')
  
  # write logk's
  iIdx+=1
  sIdx='_'+str(iIdx).zfill(4)
  #if only25:
  if '.' in sEcform: fo.write('!')
  fo25.write('%s\t%s\t%s\t%s\t' % (sTyp,sSource,sIdx,namSp))
  fo25.write('%s\t' % sEcform)
  fo25.write('%.4g\t' % size)
  fo25.write('%.4g\t' % -logks[1])
  fo25.write('\n')
  #else:
  if '.' in sEcform: fo25.write('!')
  fo.write('%s\t%s\t%s\t%s\t' % (sTyp,sSource,sIdx,namSp))
  fo.write('%s\t' % sEcform)
  fo.write('%.4g\t' % size)
  for x in logks: fo.write('%.4g\t' % -x)
  fo.write('\n')

fi.close()
fo.close()
fo25.close()
fw.close()
#----------------------------------------------------------read minerals
sys.exit()

#----------------------------------------------------------ChargeReplace
def ChargeReplace(w):
  #w= w.replace(' + ', ' +')
  w= w.replace('++++++','+6')
  w= w.replace('+++++', '+5')
  w= w.replace('++++',  '+4')
  w= w.replace('+++',   '+3')
  w= w.replace('++',    '+2')
  #
  #w= w.replace(' - ',   ' -')
  w= w.replace('------','-6')
  w= w.replace('-----', '-5')
  w= w.replace('----',  '-4')
  w= w.replace('---',   '-3')
  w= w.replace('--',    '-2')
  return w
#--------------------------------------------------------//ChargeReplace
#-----------------------------------------------------------FormulaClean
def FormulaClean(w):
  if ':' in w:
    w= w.replace(':H2O',   '(H2O)'    )
    w= w.replace(':2H2O',  '(H2O)2'   )
    w= w.replace(':3H2O',  '(H2O)3'   )
    w= w.replace(':4H2O',  '(H2O)4'   )
    w= w.replace(':5H2O',  '(H2O)5'   )
    w= w.replace(':6H2O',  '(H2O)6'   )
    w= w.replace(':7H2O',  '(H2O)7'   )
    w= w.replace(':8H2O',  '(H2O)8'   )
    w= w.replace(':10H2O', '(H2O)10'  )
    w= w.replace(':26H2O', '(H2O)26'  )
    w= w.replace(':5.5H2O','(H2O)11/2')
  return w
#---------------------------------------------------------------------//

fi= open(namInn+".dbs",'r')
fo= open(namInn+".dtb",'w')
for line in fi:
  line= ChargeReplace(line)
  line= FormulaClean(line)
  fo.write(line)
fi.close()
fo.close()

sys.exit()


