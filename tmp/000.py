import re
import sys
import numpy as np
from fractions import Fraction
import pkgutil

def read_masses(): 
    """
    A simple function to read a file with a two column list of 
    elements and their masses into a dictionary
    """
    datastream = pkgutil.get_data('burnman', 'data/input_masses/atomic_masses.dat')
    datalines = [ line.strip() for line in datastream.split('\n') if line.strip() ]
    lookup=dict()
    for line in datalines:
        data="%".join(line.split("%")[:1]).split()
        if data != []:
            lookup[data[0]]=float(data[1])
    return lookup

def read_masses_mine(): 
    """
    A simple function to read a file with a two column list of 
    elements and their masses into a dictionary
    """
    datastream = open('data/atomic_masses.dat')
    datalines = [ line.strip() for line in datastream if line.strip() ]
    lookup=dict()
    for line in datalines:
        data="%".join(line.split("%")[:1]).split()
        if data != []:
            lookup[data[0]]=float(data[1])
    return lookup

def dictionarize_formula(formula):
    """
    A function to read a chemical formula string and 
    convert it into a dictionary
    """
    f=dict()
    elements=re.findall('[A-Z][^A-Z]*',formula)
    for element in elements:
        element_name=re.split('[0-9][^A-Z]*',element)[0]
        element_atoms=re.findall('[0-9][^A-Z]*',element)
        if len(element_atoms) == 0:
            element_atoms=Fraction(1.0)
        else:
            element_atoms=Fraction(element_atoms[0])
        f[element_name]=f.get(element_name, 0.0) + element_atoms

    return f

def formula_mass(formula, atomic_masses):
    """
    A function to take chemical formula and atomic mass
    dictionaries and 
    """
    mass=sum(formula[element]*atomic_masses[element] for element in formula)
    return mass

s= 'NaAlSi3O8'
formula= dictionarize_formula(s)
print formula

atomic_masses= read_masses_mine()
print atomic_masses

for mass in atomic_masses:
  print mass
  
x= formula_mass(formula, atomic_masses)
print x

sys.exit()


