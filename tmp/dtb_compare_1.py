import sys, os

#-----------------------------------------------------------function num
import exceptions
def num(s):
  try:
    return float(s)
  except exceptions.ValueError:
    return 0.
#------------------------------------------------------------------//num

sdir= '../test/dtb/'

s1= "hsvthr_thcalc_db55_min.tab"
s2= "hsvthr_twq_dec06_min.dtb"

s1= sdir+s1
s2= sdir+s2

fi_1= open(s1,'r')
fi_2= open(s2,'r')

fCmp= open('compare_HSV.tab','w')

ok= 'z'

vNam_1=[]
vNam_2=[]
vDat_1=[]
vDat_2=[]

for line in fi_1:
  ww= line.split()
  if(len(ww)>11):
    nam= ww[1].lower()
    dat= (num(ww[8]),num(ww[9]),num(ww[10]))
    vNam_1.append(nam)
    vDat_1.append(dat)
    print nam
    #fCmp.write('%s\t%s\n' %(ww[0].lower(),ww[4].lower()))
raw_input()

for line in fi_2:
  ww= line.split()
  if(len(ww)>11):
    nam= ww[1].lower()
    dat= (num(ww[8]),num(ww[9]),num(ww[10]))
    vNam_2.append(nam)
    vDat_2.append(dat)
    if nam in vNam_1:
      j= vNam_1.index(nam)
      H1,S1,V1= vDat_1[vNam_1.index(nam)]
      H2,S2,V2= dat
      fCmp.write('%s\t'  % nam)
      if H1+H2 != 0.: delta= (H1-H2)/(H1+H2)*100.
      else:        delta= 0.
      fCmp.write('%.4g\t%.4g\t%.4g\t' % (H1,H2,delta) )
      if S1+S2 != 0.: delta= (S1-S2)/(S1+S2)*100.
      else:        delta= 0.
      fCmp.write('%.4g\t%.4g\t%.4g\t' % (S1,S2,delta) )
      if V1+V2>0.: delta= (V1-V2)/(V1+V2)*100.
      else:        delta= 0.
      fCmp.write('%.4g\t%.4g\t%.4g\t' % (V1,V2,delta) )
      fCmp.write('\n')

fi_1.close()
fi_2.close()
fCmp.flush()
fCmp.close()

sys.exit()
#=======================================================================

fi_1= open(s1,'r')
fi_2= open(s2,'r')

fCmp= open('compare_S.tab','w')

ok= 'z'

vNam_1=[]
vNam_2=[]
vDat_1=[]
vDat_2=[]

for line in fi_1:
  ww= line.split()
  if(len(ww)>0):
    w1A= ww[0].lower()
    w1B= ww[5].lower()
    vNam_1.append(w1A)
    vDat_1.append(w1B)
    print w1A, w1B
    #fCmp.write('%s\t%s\n' %(ww[0].lower(),ww[4].lower()))

for line in fi_2:
  ww= line.split()
  if(len(ww)>0):
    w2A= ww[1].lower()
    w2B= ww[9].lower()
    if w2A in vNam_1:
      j= vNam_1.index(w2A)
      print 
      fCmp.write( '%s\t%s\t%s\n' % (w2A, w2B, vDat_1[j]) )

fi_1.close()
fi_2.close()
fCmp.flush()
fCmp.close()

#=======================================================================

fi_1= open(s1,'r')
fi_2= open(s2,'r')

fCmp= open('compare_V.tab','w')

ok= 'z'

vNam_1=[]
vNam_2=[]
vDat_1=[]
vDat_2=[]

for line in fi_1:
  ww= line.split()
  if(len(ww)>0):
    w1A= ww[0].lower()
    w1B= ww[6].lower()
    vNam_1.append(w1A)
    vDat_1.append(w1B)
    print w1A, w1B
    #fCmp.write('%s\t%s\n' %(ww[0].lower(),ww[4].lower()))

for line in fi_2:
  ww= line.split()
  if(len(ww)>0):
    w2A= ww[1].lower()
    w2B= ww[10].lower()
    if w2A in vNam_1:
      j= vNam_1.index(w2A)
      print 
      fCmp.write( '%s\t%s\t%s\n' % (w2A, w2B, vDat_1[j]) )

fi_1.close()
fi_2.close()
fCmp.flush()
fCmp.close()

#=======================================================================

fi_1= open(s1,'r')
fi_2= open(s2,'r')

fCmp= open('compare_Cp.tab','w')

ok= 'z'

vNam_1=[]
vNam_2=[]
vDat_1=[]
vDat_2=[]

for line in fi_1:
  ww= line.split()
  if(len(ww)>0):
    w1A= ww[0].lower()
    w1B= ww[8].lower()
    vNam_1.append(w1A)
    vDat_1.append(w1B)
    print w1A, w1B
    #fCmp.write('%s\t%s\n' %(ww[0].lower(),ww[4].lower()))

for line in fi_2:
  ww= line.split()
  if(len(ww)>0):
    w2A= ww[1].lower()
    w2B= ww[11].lower()
    if w2A in vNam_1:
      j= vNam_1.index(w2A)
      print 
      fCmp.write( '%s\t%s\t%s\n' % (w2A, w2B, vDat_1[j]) )

fi_1.close()
fi_2.close()
fCmp.flush()
fCmp.close()
