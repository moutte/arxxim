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

el_select= ['Al','Fe','Mg','Ca','Na','K','Si']
el_select= el_select + ['C','S','P','F','Cl']
el_select= el_select + ['O','H']
el_select= el_select + ['+','-']

list_select= open("list_select.txt",'w')

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
    if 0:
      if 'C' in elements:
        coef= coeffs[elements.index('C')]
        cod= 'C'
      else:
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

