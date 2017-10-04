
# it is standard practice in Julia to capitalize type names
type Element
  name::        String
  weight::      String #Float32
end

type Species
  typ::         String
  name::        String
  formula::     String
  # dim::         Integer
  # elements::    Element[] #(dim)
  # stoichio::    Integer[] #(dim)
end

function elements_readfile(path) #-------------------------read elements
  f=open(path)
  lines = readlines(f)
  close(f)
  elements=Element[] # Vector{Element}[]
  for line in lines
    if line[1] in ['%','!']   ; continue ; end
    w= split(line)
    #println(w[1]," ",w[2])
    if length(w)>1
      element= Element(ucfirst(lowercase(w[1])),w[2])
      #println(element.name," ",element.weight)
      push!(elements,element)
    end
  end
  return elements
end #----------------------------------------------------//read elements

function species_readfile(path) #---------------------------read species
  f=open(path)
  lines = readlines(f)
  close(f)
  
  species= Species[]
  for line in lines
    if line[1] in ['%','!']   ; continue ; end
    w= split(line)
    #println(w[1]," ",w[2])
    if length(w)>1
      x= Species(w[1],w[2],w[3])
      #println(element.name," ",element.weight)
      push!(species,x)
    end
  end
  return species
end #-----------------------------------------------------//read species

function formula_read(s) #----------------------------------formula_read
  list_elements= String[]
  list_stoikio=   Int8[]
  w1= split(s,')')
  for x in w1
    if '(' in x
      w2= split(x,'(')
      y= w2[1] # ;  print y
      z= w2[2] # ;  print z
      push!(list_elements,ucfirst(lowercase(y)))
      push!(list_stoikio,parse(Int,z))
    end
  end
  return list_elements, list_stoikio
end #-----------------------------------------------------//formula_read

v_element_all= elements_readfile("data/atomic_masses.dat")
for x in v_element_all
  println(x.name," ",x.weight)
end
list_elements_all= [x.name for x in v_element_all]

v_species= species_readfile("data/species.dat")
v_element= Element[]
list_elements=String[]

for spc in v_species
  println(spc.typ," ",spc.name," ",spc.formula)
  (elements, coeffs)= formula_read(spc.formula)
  for s in elements
    if ! (s in list_elements)
      push!(list_elements,s)
      if s in list_elements_all
        i= findfirst(list_elements_all,s)
        push!(v_element,v_element_all[i])
      end
    end
  end
end
# println(list_elements)
list_elements= [x.name for x in v_element]
println(list_elements)

#---------------------------------------------------------t_stoikio_read 
function t_stoikio_read(v_element,v_species)
  t_stoikio= Int8[length(v_element),length(v_species)]
  t_stoikio= zeros(Int8,length(v_element),length(v_species))
  for j in eachindex(v_species)
    (elements, coeffs)= formula_read(v_species[j].formula)
    for k in eachindex(elements)
      i= findfirst(list_elements,elements[k])
      println(j,k,i)
      if i>0
       t_stoikio[i,j]= coeffs[k]
      end
    end
  end
  return t_stoikio
end
#-------------------------------------------------------//t_stoikio_read

t_stoikio= t_stoikio_read(v_element,v_species)

for j in eachindex(v_species)
  println(v_species[j].name)
  println(t_stoikio[:,j])
end

type PhaseModel
  name::        String
  species::     Vector{Species}
  eos::         String
end

type Phase
  name::        String
  model::       String
  composition:: Vector{Float32}
end
  
type ChemicalSpace
  name::        String
  elements::    Vector{Element}
  species::     Vector{Species}
  models::      Vector{PhaseModel}
end

function parab3p(sigma0,sigma1, lambdac,lambdam, ff0,ffc,ffm)
# Apply three-point safeguarded parabolic model for a line search.
# C. T. Kelley, April 1, 2003
# input:
#   sigma0,sigma1= parameters, safeguarding bounds for the linesearch
#   lambdac=  current steplength
#   lambdam=  previous steplength
#   ff0=  |F(x_c)|^2
#   ffc=  |F(x_c + lambdac)|^2
#   ffm=  |F(x_c + lambdam)|^2
# output:
#   lambdap=  new value of lambda given by parabolic model
#
#-------------------- Compute coefficients of interpolation polynomial
# p(lambda)=  ff0 + (-c1 lambda + c2 lambda^2)/d1
# d1=  (lambdac - lambdam)*lambdac*lambdam < 0
# so, if c2 > 0 we have negative curvature and default to
# lambdap=  sigma1 * lambda.
#
  c2= lambdam*(ffc-ff0)-lambdac*(ffm-ff0)
  if(c2 >=0.0)
    lambdap=  sigma1*lambdac
  else
    c1= lambdam *lambdam *(ffc-ff0) - lambdac *lambdac *(ffm-ff0)
    lambdap= c1 *0.5 /c2
    lambdap= max(lambdap, sigma0*lambdac)
    lambdap= min(lambdap, sigma1*lambdac)
  end
  
  return lambdap
end

"""
function Newton_Walker()
  SigmaMax= 0.5
  SigmaMin= 0.1
  Tau=      1.0E-4
  iArmMax= 12
  INTEGER :: iArm
  REAL(dp):: Norm_vF0,Norm_vF,delta,lambda,Sigma
  REAL(dp):: X
  !
  ForcePositive= COUNT(vIsPlus(:))>0
  !
  vTolF(:)= TolF
  !
  iErr=-1
  !
  for Its=1:MaxIts
  
    vX0(:)= vX(:)
    vFunc0(:)= Residual(vX0)
    Norm_vF0= SQRT(SUM(vFunc0(:)*vFunc0(:)))
    
    if(bFinDIF)
      CALL Jacobian_Numeric(Residual,vX0,vFunc0,tJac)
    else
      CALL Jacobian(vX0,tJac)
    end
    
    #----------------------solve linear equations using LU decomposition
    CALL LU_Decomp(tJac, vIndex, D, bSingul)
    if(bSingul)
      iErr=-2
      break
    end
    vDX(:)= -vFunc0(:)
    CALL LU_BakSub(tJac,vIndex,vDX)
    #---------------------/solve linear equations using LU decomposition
    
    if(ForcePositive)
      X= MAXVAL(-vDX(:)/vX0(:)*2.0D0,MASK=vIsPlus(:))
      if(iDebug>2 .AND. X>1.0D0) write(69,*) "Walker-UR, ", X
      if(X>1.0D0) 
        vDX(:)= vDX(:) /X
      end
    end
    !
    vX(:)= vX0(:) + vDX(:)
    !
    vFunc(:)= Residual(vX)
    Norm_vF= SQRT(SUM(vFunc(:)*vFunc(:)))
    !
    lambda= 1.0d0
    iArm=   0
    !
    !---------------------------Test the step and backtrack as necessary
    DO WHILE ( Norm_vF > (1.0d0-tau*lambda)*Norm_vF0 )
      iArm= iArm +1
      
      if(iArm > iArmMax)  #error('Maximum number of backtracks reached.')
        iErr= -4
        break
      end
      !
      delta= (Norm_vF /Norm_vF0)**2 -1.0d0 +2.0d0*lambda
      
      if (delta > Zero)
        Sigma=  lambda/delta
        Sigma=  MIN(Sigma,SigmaMax)
        Sigma=  MAX(Sigma,SigmaMin)
      else
        Sigma=  SigmaMax
      end
      
      vDX(:)=  Sigma*vDX(:)
      lambda=  Sigma*lambda
      vX(:)=   vX0(:) + vDX(:)
      vFunc(:)= Residual(vX)
      Norm_vF= SQRT(SUM(vFunc(:)*vFunc(:)))
      
    end
    #--------------------------/Test the step and backtrack as necessary
    Error_F= MAXVAL(ABS(vFunc(:)))
    Delta_X= MAXVAL(ABS(vDX(:)))

    #-- convergence --
    if(Converge(vFunc,vTolF))
      iErr= 0
      break
    end
    !
    #-- stationary --
    if(Delta_X < TolX)
      iErr=-5
      break
    ENDIF
    !
  end
  
  nIts=Its
  
end
"""




