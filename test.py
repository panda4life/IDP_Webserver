# -*- coding: utf-8 -*-
"""
Created on Wed Mar 12 11:02:07 2014

@author: jahad
"""

from sequence import Sequence
seq = Sequence('EEEEEEEEEEEKKKKKKAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA')
print('N:\t ' + `seq.len`)
print('f+:\t' + `seq.Fplus()`)
print('f:\t' + `seq.Fminus()`)
print('NCPR:\t' + `seq.NCPR()`)
print('FCR:\t' + `seq.FCR()`)
print('D:\t' + `seq.delta()`)
print('Dmax:\t' + `seq.deltaMax()`)
print('Kappa:\t' + `seq.kappa()`)
print('<H>:\t' + `seq.meanHydropathy()`)
#assert(seq.meanHydropathy() == seq.cumMeanHydropathy()[-1])

#import random as rm
#rand= rm.Random()
#rand.shuffle('0123456789')