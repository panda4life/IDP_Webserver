# -*- coding: utf-8 -*-
"""
Created on Wed Mar 12 11:10:55 2014

@author: jahad
"""

import numpy as np

class Sequence:
    def __init__(self,seq):
        self.seq = seq
        self.len = len(seq)
    """    
    def FCR(self):
        
    def NCPR(self):
        
    def Fminus(self):
        
    def Fplus(self):
        
    def meanHydropathy(self):
        
    def cumMeanHydropathy(self):
        
    def sigma(self):
        
    def delta(self):
        
    def kappa(self):
        """
    def swapRes(self,index1,index2):
        if(index1 == index2):
            return Sequence(self.seq)
        elif(index2<index1):
            temp = index2
            index2 = index1
            index1 = temp
        else:
            pass
        
        tempseq = self.seq[:index1] + self.seq[index2] + self.seq[(index1+1):(index2)]+ self.seq[index1] + self.seq[(index2+1):]
        return Sequence(tempseq)