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
idxNam,idxH,idxS,idxV= 1,8,9,10
HSV= True
LOG= False

s1= sdir+s1
s2= sdir+s2

fi_1= open(s1,'r')
fi_2= open(s2,'r')

if HSV: fCmp= open("compare_HSV.tab",'w')
if LOG: fCmp= open("compare_logk.tab",'w')

ok= 'z'

vNam_1=[]
vNam_2=[]
vDat_1=[]
vDat_2=[]

#--------------------------------------------------------------------HSV
if HSV:
  for line in fi_1:
    ww= line.split()
    if(len(ww)>11):
      nam= ww[idxNam].lower()
      dat= (num(ww[idxH]),num(idxS),num(idxV))
      vNam_1.append(nam)
      vDat_1.append(dat)
      print nam
      #fCmp.write('%s\t%s\n' %(ww[0].lower(),ww[4].lower()))
  raw_input()

  for line in fi_2:
    ww= line.split()
    if(len(ww)>11):
      nam= ww[idxNam].lower()
      dat= (num(ww[idxH]),num(ww[idxS]),num(ww[idxV]))
      vNam_2.append(nam)
      vDat_2.append(dat)
      if nam in vNam_1:
        print nam
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
#------------------------------------------------------------------//HSV

fi_1.close()
fi_2.close()
fCmp.flush()
fCmp.close()

sys.exit()
#=======================================================================

