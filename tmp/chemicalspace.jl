type Element
  name::   String
  weight:: String  #Float32
end

type Species
  name::        String
  typ::         String
  formula::     String
  list_elem::   Vector{String}
  list_coeffs:: Vector{Int8}
  #list_elem::   String[]
  #list_coeffs:: String[]
  #x::Array{Float64,1}
  #dim::         Integer
  #elements::    Array{Element} #[dim]
  #stoichio::    Array{Integer} #[dim]
end

type PhaseModel
  name::    String
  species:: Vector{Species}
  eos::     String
end

type Phase
  name::        String
  model::       PhaseModel
  composition:: Vector{Float64}
end
  
type ChemicalSpace
  name::        String
  elements::    Vector{Element}
  species::     Vector{Species}
  models::      Vector{PhaseModel}
end

function formula_read(s)
  list_elements= String[]
  list_coeffs=   Int8[]
  if contains(s,")")
    w1= split(s,")")
    for x in w1
      if contains(x,"(")
        w2= split(x,"(")
        y= ucfirst(lowercase(w2[1])) # ;  print y
        z= parse(Int8,w2[2]) # ;  print z
        push!(list_elements, y)    # list_elements.append(y.capitalize())
        push!(list_coeffs,   z)    # list_coeffs.append(z)
      end
      #println(list_elements)
    end
  end
  return list_elements, list_coeffs
end

if false
s= "AL(1)O(1)H(1)+(2)"
list_elements, list_coeffs= formula_read(s)
println(list_elements)
end

#-----------------------------------------------------------read species
function species_readfile(path)
  f=open(path)
  lines = readlines(f)
  close(f)

  species= Species[]
  for line in lines
    if strip(line)==""        ; continue ; end
    if startswith(line,"%")   ; continue ; end
    w= split(line)
    if length(w)>2
      list_elements, list_coeffs= formula_read(w[3])
      x= Species(w[2],w[1],w[3],list_elements, list_coeffs)
      push!(species, x)
    end
  end
  
  return species
end
#---------------------------------------------------------//read species

#----------------------------------------------------------read elements
function elements_readfile(path)
  f=open(path)
  lines = readlines(f)
  close(f)
  
  elements= Element[]
  for line in lines
    if strip(line)==""        ; continue ; end
    if startswith(line,"%")   ; continue ; end
    w= split(line)
    if length(w)>1
      x= Element(ucfirst(lowercase(w[1])),w[2])
      push!(elements, x)
    end
  end
  
  return elements
end
#--------------------------------------------------------//read elements

#-----------------------------------------------------------reading test
elements= elements_readfile("data/atomic_masses.dat")
if false
  for x in elements
    println(x.name," ",x.weight)
  end
end

species= species_readfile("data/species.dat")
if false
  for x in species
    println( x.typ," ", x.name," ", x.formula ) #, x.list_elem_names, x.list_coeffs )
  end
end
if false
  for x in species
    println(x.list_coeffs )
  end
end

elements_all= elements
list_elem_all= [x.name for x in elements_all]
elements= Element[]
elem_names= String[]
for spc in species
  for ele in spc.list_elem
    #println(ele)
    if ele in list_elem_all
      if !(ele in elem_names)
        push!(elem_names, ele)
        j= findfirst(list_elem_all,ele)
        push!(elements, elements_all[j])
      end
    end
  end
end

#println(elem_names)
for x in elements
  println(x.name," ",x.weight)
end

formula_table= [0 for i in 1:length(species), c in 1:length(elements)]
#formula_table= zeros(length(species),length(elements))
i=1
for spc in species
  j=1
  for ele in spc.list_elem
    k= findfirst(elem_names,ele)
    formula_table[i,k]= spc.list_coeffs[j]
    j+=1
  end
  i+=1
end

for i in 1:length(species)
  println(formula_table[i,:]," ",species[i].name)
end

#println(formula_table)
#---------------------------------------------------------//reading test

mixmodels= PhaseModel[]
model= PhaseModel("name_1",species,"ideal")
push!(mixmodels,model)
model= PhaseModel("name_2",species,"ideal")
push!(mixmodels,model)

myspace= ChemicalSpace("myspace", elements,species,mixmodels)

quit()


