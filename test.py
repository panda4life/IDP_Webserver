# -*- coding: utf-8 -*-
"""
Created on Wed Mar 12 11:02:07 2014

@author: jahad
"""

from sequence import Sequence
seq = Sequence('EEEEEEEEEEKKKKK')
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
import time
rng = rand.Random()
rng.seed(time.time())
start = time.time()
for i in range(0,100000):
    swappair = rng.sample(range(0,seq.len),2)
    test = seq.swapRandChargeRes().delta()
end = time.time()
runtime = end-start
print(runtime)
'''
import wanglandau as wl
import numpy as np
seq = Sequence('QHGQLWFPEGFKVSEASKKKRREPLGEDSAGLKPLKNASDGALMDDNQNEWGDEDLETKKFRFEEPVVRGDLDDQTDHRQWTQQHLDAADLSMSAMAPTPPQGEVDADCMDVNVRGPDGF')
print(seq.toFileString())
g = wl.WangLandau(seq,nbins=10,binmin=.002,binmax=.452,nflatchk=10000,convergence=np.exp(.000001),genPermutants=False)
#assert(seq.meanHydropathy() == seq.cumMeanHydropathy()[-1])

#import random as rm
#rand= rm.Random()
#rand.shuffle('0123456789')