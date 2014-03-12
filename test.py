# -*- coding: utf-8 -*-
"""
Created on Wed Mar 12 11:02:07 2014

@author: jahad
"""

from residues import resTable
rs = resTable('residueData.csv')

from sequence import Sequence
seq = Sequence('0123456789')
seq = seq.swapRes(2,5)
print(seq.seq)