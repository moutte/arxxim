import sys
import math as m
import pylab as plt

def formula_read(s):
  if s[-1]!=')':
    s= s+"(1)"
  list_elements= []
  list_coeffs= []
  w1= s.split(')')
  for x in w1:
    if '(' in x:
      w2= x.split('(')
      y= w2[0] # ;  print y
      z= w2[1] # ;  print z
      list_elements.append(y.capitalize())
      list_coeffs.append(z)
  return list_elements, list_coeffs
  
#--------------------------------------------------compute analytic logK
#--subroutine DtbLogKAnl_Calc
def logK(fit,T):
  # PHREEQC
  x= fit[0] + fit[1]*T + fit[2]/T + fit[3]*m.log(T,10) + fit[4]/T/T
  if len(fit)>5: x= x + fit[5]*T*T
  return x

s= "dtb_anl_test.tab"
instream = open(s,'r')
coeffs = [ line.strip() for line in instream if line.strip() ]
for line in coeffs:
  ww= line.split()
  Min= ww[0]
  fit= []
  for i,w in enumerate(ww):
    if i>1: fit.append(float(w))
  print Min,fit
  #
  Xmin,Xmax,Xdim= 0.,150.,7
  Tser= plt.linspace(Xmin,Xmax,num=Xdim)
  logKs=[]
  for T in Tser:
    TK= T + 273.15
    x= -logK(fit,TK) # MINUS SIGN !!!
    print T,x
    logKs.append(x)
  #raw_input()
  #
  fig= plt.subplot(1,1,1)
  fig.grid(color='r', linestyle='-', linewidth=0.2)
  fig.grid(True)
  fig.plot(Tser, logKs, 'o',linestyle='-',linewidth=1.0)
  fig.set_title(Min)
  plt.show()
  #
sys.exit()
#------------------------------------------------//compute analytic logK

#---------------------------------------------sorting species by formula
s= "list-species.txt"
# formulas= open(s,'r').readlines()
instream = open(s,'r')
formulas = [ line.strip() for line in instream if line.strip() ]

select_maj= ['Al','Fe','Mg','Ca','Na','K','Si']
select_trc= ['Sr','Ba'] + ['Mn','V','Cr','Ni','Cu','Zn'] + ['Pb','Mo','As']
select_ani= ['C','S','N','P','F','Cl'] + ['O','H'] + ['+','-']

el_select_maj= select_maj + select_ani
el_select_trc= select_trc + select_maj + select_ani

list_select= open("list-select.txt",'w')

for f in formulas:
  elements,coeffs= formula_read(f)
  #print list_elements
  #raw_input()
  isvalid= True
  isvalid_trc= False
  for el in elements:
    if not el in el_select_trc:
      isvalid= False
      break
    else:
      if el in select_trc: isvalid_trc= True
  if isvalid_trc: cod= 'TRC'
  elif isvalid:   cod= 'MAJ'
  else:           cod= 'Z'
  ll=''
  if isvalid:
    for el in el_select_trc:
      if el in elements:
        i= elements.index(el)
        if coeffs[i]!='0': 
          ll= ll + el + '('+coeffs[i]+')' 
  else:
    for i,el in enumerate(elements):
      if coeffs[i]!='0': 
        ll= ll + el + '('+coeffs[i]+')'
  list_select.write("%s\t%s\n" % (ll,cod))
  
list_select.close()

sys.exit()
#-------------------------------------------//sorting species by formula

#-------------------------------------------------------------test log_P
import math as m

r= m.sqrt(10.)
r= m.sqrt(r)
r= m.sqrt(r)
print r
x= 1.
#r= 1.333521
for i in range(50):
  x= x*r
  if x>5000.: break
  print i,x,m.log(x,10)

sys.exit()
#-------------------------------------------------------------test log_P

