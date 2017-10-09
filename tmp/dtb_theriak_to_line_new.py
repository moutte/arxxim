import os

sSource= 'DEC06'
s= "thr_twq_dec06_min"

os.chdir("../test/dtb")

fInn= open(s+".dtb",'r')
fOut= open(s+"_line.dtb",'w')

# DEFAULT FIELD LIST in dtb_read_dtbminthr :
# "TYPE NAME ECFORM SKIP SOURCE SKIP PARAMETERS"

Ok= False

for line in fInn:
  words= line.split()
  
  if len(words)>0:
    w0= words[0]
    if w0[0]=='!': continue
    if w0=='MINERAL':
      Ok= True
      sTyp= 'MIN'
      continue
    if w0=='GAS':
      Ok= True
      sTyp= 'GAS'
      continue
    if w0=='END': Ok= False
    if w0=='ENDMINERAL': Ok= False
    if w0=='ENDGAS': Ok= False
    
    if Ok:
      print line
      while len(words)<7: words.append('_')
      if words[0]=='&':
        for j in range(6):
          fOut.write('%s\t' % (words[j+1]))
      else:
        fOut.write('\n')
        #for j in range(5):
        #  fOut.write('%s\t' % (words[j]))
        fOut.write('%s\t' % sTyp)
        fOut.write('%s\t' % words[0])
        fOut.write('%s\t' % words[1])
        fOut.write('%s\t' % words[2])
        fOut.write('%s\t' % sSource)
        fOut.write('%s\t' % words[3])

fInn.close()
fOut.close()
