print '\n'
print 'tests sur les classes'
print '\n'

class Mineral:
  def __init__(self, str):
    self.nom= str
    self.G=   0.
    self.H=   0.
    self.S=   0.
  def printIt(self):
    print "nom=", self.nom
    print "G=  ", self.G
    print "H=  ", self.H
  def readFromLine(self,str):
    words= str.split()
    self.nom= words[0]
    self.G=   float(words[1])

Min_1= Mineral("")
str= "Andalusite 0.100 0.200 0.300"
Min_1.readFromLine(str)

Min_2= Mineral("Sillimanite")

MinList= [Min_1,Min_2]

Min_1.printIt()

#class C:
  #def __init__(self, val): self.val = val
  #def f(self): print "hello, my value is:", self.val

#a = C(27)
#a.f()

#a.val= 33
#a.f()
