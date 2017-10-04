import sys

class Element(object):
  def __init__(self,name,weight,entropy):
    self.name=    name
    self.weight=  weight
    self.entropy= entropy

class Species(object):
  def __init__(self,name,typ,formula):
    self.name=    name
    self.typ=     typ
    self.formula= formula
    self.list_elem_names,self.list_coeffs= formula_read(formula)

def formula_read(s):
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

class PhaseModel(object):
  def __init__(self,name,eos):
    self.name=    name
    self.eos=     eos
    self.species= []
  def add(self, x):
    self.species.append(x)

class Phase(object):
  def __init__(self,name,model,composition):
    self.name= name
    self.model= model
    self.composition= composition
  
class ChemSpace(object):
  def __init__(self,name,elements,species,models):
    self.name= name
    self.elements= elements
    self.species= species
    self.models= models



#----------------------------------------------------------read elements
def elements_readfile(path):
  f=open(path,'r')
  elements= []
  for line in f:
    if line.strip()=='': continue
    if line[0]=="%":     continue
    w= line.split()
    if len(w)>1:
      x= Element(name=w[0].capitalize(),entropy=0.,weight=w[1])
      elements.append(x)
  f.close()
  return elements
#--------------------------------------------------------//read elements

#-----------------------------------------------------------read species
def species_readfile(path):
  f=open(path,'r')
  species= []
  for line in f:
    if line.strip()=='': continue
    if line[0]=="%": continue
    w= line.split()
    if len(w)>2:
      x= Species(typ=w[1],name=w[1],formula=w[2])
      species.append(x)
  f.close()
  return species
#---------------------------------------------------------//read species

#-----------------------------------------------------------reading test
elements= elements_readfile("data/atomic_masses.dat")
for x in elements:
  print x.name, x.weight

species= species_readfile("data/species.dat")
for x in species:
  print x.typ, x.name, x.formula, x.list_elem_names, x.list_coeffs

elements_all= elements
list_elem_all= [x.name for x in elements_all]
elements= []
elem_names= []
for spc in species:
  for ele in spc.list_elem_names:
    if ele in list_elem_all:
      if not ele in elem_names:
        elem_names.append(ele)
        j= list_elem_all.index(ele)
        elements.append(elements_all[j])
print elem_names
for x in elements:
  print x.name, x.weight

water_model= PhaseModel("WATER","DILUTE")
for spc in species:
  if spc.typ=="AQU":
    water_model.species.append(spc)
  
def formula_table_compute(species,elem_names):
  formula_table= []
  for spc in species:
    v= [0 for j in range(len(elem_names))]
    i=0
    for ele in spc.list_elem_names:
      #print ele
      j= elem_names.index(ele)
      v[j]= spc.list_coeffs[i]
      i+=1
    formula_table.append(v)
  return formula_table

formula_table= formula_table_compute(species,elem_names)

print formula_table
#---------------------------------------------------------//reading test

models= []
model= PhaseModel(name="name_1",eos="ideal")
for spc in species:
  model.add(spc)
models.append(model)

model= PhaseModel(name="name_2",eos="ideal")
for spc in species:
  model.add(spc)
models.append(model)

myspace= ChemSpace(
  name="myspace",
  elements= elements,
  species= species,
  models= models)

sys.exit()


