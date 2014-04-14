# -*- coding: utf-8 -*-
"""
Created on Wed Mar 12 11:02:07 2014

@author: jahad
"""

from sequence import Sequence
seq = Sequence('QHGQLWFPEGFKVSEASKKKRREPLGEDSAGLKPLKNASDGALMDDNQNEWGDEDLETKKFRFEEPVVRGDLDDQTDHRQWTQQHLDAADLSMSAMAPTPPQGEVDADCMDVNVRGPDGF')
'''
print('N:\t ' + `seq.len`)
print('f+:\t' + `seq.Fplus()`)
print('f:\t' + `seq.Fminus()`)
print('NCPR:\t' + `seq.NCPR()`)
print('FCR:\t' + `seq.FCR()`)
print('D:\t' + `seq.delta()`)
print('Dmax:\t' + `seq.deltaMax()`)
print('Kappa:\t' + `seq.kappa()`)
print('<H>:\t' + `seq.meanHydropathy()`)
'''
#print(seq.toFileString())
'''
import random as rand
import time as t
import numpy as np
rng = rand.Random()
rng.seed(t.time())

seq = Sequence('QHGQLWFPEGFKVSEASKKKRREPLGEDSAGLKPLKNASDGALMDDNQNEWGDEDLETKKFRFEEPVVRGDLDDQTDHRQWTQQHLDAADLSMSAMAPTPPQGEVDADCMDVNVRGPDGF')
start = t.time()
for i in range(0,10000):
    #swapPair = rng.sample(np.arange(0,seq.len),2)
    #posInd = np.where(seq.chargePattern>0)[0]
    #negInd = np.where(seq.chargePattern<0)[0]
    #neutInd = np.where(seq.chargePattern==0)[0]
    #test = seq.swapRes(swapPair[0],swapPair[1])
    test = seq.swapRandChargeRes().delta()

    #charge1 = seq.chargePattern[swapPair[0]]
    #charge2 = seq.chargePattern[swapPair[1]]
    #tempChargeSeq = seq.chargePattern
    #tempChargeSeq[swapPair[0]] = charge2
    #tempChargeSeq[swapPair[1]] = charge1

    #seq.delta()
end = t.time()
runtime = end-start
print(runtime)
'''
'''
binmin = .002
nbins = 10
binsz = .01
startTime = t.time()
for testNum in np.arange(0,10000):
    k = rng.uniform(.002,.452)
    idx_new = int((k-binmin)/binsz)
endTime = t.time()
print(endTime-startTime)

startTime = t.time()
for testNum in np.arange(0,10000):
    k = rng.uniform(.002,.452)
    for i in np.arange(0,nbins):
        if(k-binmin-binsz*i<binsz):
            idx_new = i
            break
    else:
        idx_new = nbins+1
endTime = t.time()
print(endTime-startTime)
'''

import wanglandau as wl
import numpy as np
seq = Sequence('EEEEEEEEEEEEEEEEEEEEEEEEEKKKKKKKKKKKKKKKKKKKKKKKKK')
#seq = Sequence('QHGQLWFPEGFKVSEASKKKRREPLGEDSAGLKPLKNASDGALMDDNQNEWGDEDLETKKFRFEEPVVRGDLDDQTDHRQWTQQHLDAADLSMSAMAPTPPQGEVDADCMDVNVRGPDGF')
outdir = 'test/'
f = open(outdir + "seqparam.txt",'w')
f.write(seq.toFileString())
print(seq.toFileString())
f.close()
#print(seq.seq)
#print(seq.chargePattern)
#testseq = seq.swapRandChargeRes().swapRandChargeRes().swapRandChargeRes()
#print(testseq.seq)
#print(testseq.chargePattern)
#print(testseq.kappa())
#print(Sequence(testseq.seq).kappa())
#print(seq.chargePattern-testseq.chargePattern)

#g = wl.WangLandau(seq,outdir,nbins=10,binmin=.02,binmax=.452)
g = wl.WangLandau(seq,nbins=10,binmin=.000,binmax=1,nflatchk=10000,convergence=np.exp(.000001),genPermutants=False)

#assert(seq.meanHydropathy() == seq.cumMeanHydropathy()[-1])

#import random as rm
#rand= rm.Random()
#rand.shuffle('0123456789')