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
select_trc= ['Sr','Ba'] + ['V','Cr','Ni','Cu','Zn'] + ['Pb','Mo','As']
select_ani= ['C','S','N','P','F','Cl'] + ['O','H'] + ['+','-']

el_select= select_maj + select_ani
el_select= select_trc + select_maj + select_ani

list_select= open("list-select.txt",'w')

for f in formulas:
  elements,coeffs= formula_read(f)
  #print list_elements
  #raw_input()
  isvalid= True
  for el in elements:
    if not el in el_select:
      isvalid= False
      break
  if isvalid:
    print f, elements
    cod= 'A'
  else:
    cod= 'Z'
  ll=''
  if cod=='A':
    for el in el_select:
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

