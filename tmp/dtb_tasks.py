import sys

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

