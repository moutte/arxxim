s= 'dtb_thr/thr_2010jun92_aqu'
sTyp= 'AQU'

sSource= 'tcdb55c2p'
s= 'dtb_thr/thr_2010tcdb55c2p_min'
sTyp= 'MIN'
s= 'dtb_thr/thr_2010tcdb55c2p_gas'
sTyp= 'GAS'

sSource= '2010_JUN92'
s= 'dtb_thr/thr_2010jun92_min'
sTyp= 'MIN'
s= 'dtb_thr/thr_2010jun92_gas'
sTyp= 'GAS'

fInn= open(s+'.dtb','r')
fOut= open(s+'_line.tab','w')

Ok= 0

for line in fInn:
  words= line.split()
  
  if len(words)>0:
    #if words[0][0] != ' ':
    if(line[0]=='*'): continue
    if(line[0]!=' '):
      #fOut.write('%s\t %s\t %s\t' % (words[0],words[1],words[2]) )
      fOut.write('\n')
      fOut.write('%s\t %s\t' % (sTyp,sSource))
      fOut.write('%s\t %s\t' % (words[0],words[1]))
    else:
      while len(words)<7: words.append('0.0')
      for j in range(6):
        fOut.write('%s\t' % (words[j]))
      
'''
      if words[0]=='MINERAL':
        Ok= 1
        sTyp= 'MIN'
        continue
      if words[0]=='GAS':
        Ok= 1
        sTyp= 'GAS'
        continue
      if words[0]=='ENDMINERAL': Ok= 0
      if words[0]=='ENDGAS': Ok= 0
      
      if Ok==1:
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
            fOut.write('%s\t' % '_')
'''

fInn.close()
fOut.close()
