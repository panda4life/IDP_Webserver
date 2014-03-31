# -*- coding: utf-8 -*-
"""
Created on Wed Mar 12 11:10:55 2014

@author: jahad
"""

import numpy as np
from residues import resTable

lkupTab = resTable('residueData.csv')

class Sequence:
    def __init__(self,seq):
        self.seq = seq.upper()
        self.len = len(seq)
        self.dmax = -1 #initializing to prevent extra computational time
        self.N_ITERS = 5
        self.N_STEPS = 5000

    def countPos(self):
        ans = 0
        for i in np.arange(0,self.len):
            if(lkupTab.lookUpCharge(self.seq[i])>0):
                ans += 1
        return ans

    def countNeg(self):
        ans = 0
        for i in np.arange(0,self.len):
            if(lkupTab.lookUpCharge(self.seq[i])<0):
                ans += 1
        return ans

    def countNeut(self):
        ans = 0
        for i in np.arange(0,self.len):
            if(lkupTab.lookUpCharge(self.seq[i])==0):
                ans += 1
        return ans

    def Fplus(self):
        return self.countPos()/(self.len+0.0)

    def Fminus(self):
        return self.countNeg()/(self.len+0.0)

    def FCR(self):
        return self.Fplus() + self.Fminus()

    def NCPR(self):
        return self.Fplus() - self.Fminus()

    def meanHydropathy(self):
        ans = 0
        for i in np.arange(0,self.len):
            ans += lkupTab.lookUpHydropathy(self.seq[i])/self.len
        return ans

    def cumMeanHydropathy(self):
        ans = [lkupTab.lookUpHydropathy(self.seq[0])]
        for i in np.arange(1,self.len):
            ans.append(ans[i-1]+lkupTab.lookUpHydropathy(self.seq[i]))
        ans /= (np.arange(0,self.len)+1)
        return ans

    def sigma(self):
        if(self.FCR() > 0):
            return self.NCPR()**2/self.FCR()
        else:
            return 0

    def deltaForm(self,bloblen):
        sigma = self.sigma()
        nblobs = self.len-bloblen+1
        ans = 0
        for i in np.arange(0,nblobs):
            blob = Sequence(self.seq[i:(i+bloblen)])
            ans += ((sigma - blob.sigma())**2)/nblobs
        return ans

    def delta(self):
        return (self.deltaForm(5)+self.deltaForm(6))/2

    def deltaMax(self):
        #if this has been computed already, then return it
        if(self.dmax != -1):
            return self.dmax
        elif(self.FCR() == 0):
            self.dmax = 0
        #first computational trick (Maximum Charge Separation)
        elif(self.countNeut() == 0):
            setupSequence = ''
            for i in np.arange(0,self.countPos()):
                setupSequence += '+'
            for i in np.arange(0,self.countNeg()):
                setupSequence += '-'
            assert(self.len == len(setupSequence))
            self.dmax = Sequence(setupSequence).delta()
        #second computational trick (Maximization of # of Charged Blobs)
        elif(self.countNeut() >= 18):
            #DEBUG
            #maxSeq = ''
            nneuts = self.countNeut()
            posBlock = ''
            negBlock = ''
            for i in np.arange(0,self.countPos()):
                posBlock += '+'
            for i in np.arange(0,self.countNeg()):
                negBlock += '-'
            for startNeuts in np.arange(0,7):
                for endNeuts in np.arange(0,7):
                    setupSequence = ''
                    midBlock = ''
                    endBlock = ''
                    for i in np.arange(startNeuts):
                        setupSequence += '0'
                    for i in np.arange(0,endNeuts):
                        endBlock += '0'
                    for i in np.arange(0,nneuts-startNeuts-endNeuts):
                        midBlock += '0'
                    setupSequence += posBlock
                    setupSequence += midBlock
                    setupSequence += negBlock
                    setupSequence += endBlock
                    #DEBUG
                    #print(setupSequence)
                    #assert(len(setupSequence) == self.len)
                    nseq = Sequence(setupSequence)
                    if(nseq.delta()>self.dmax):
                        self.dmax = nseq.delta()
                        #DEBUG
                        #maxSeq = nseq.seq
        #third computational trick (Search through set of sequences that fit DMAX pattern)
        else:
            #DEBUG
            #maxSeq = ''
            posBlock = ''
            negBlock = ''
            nneuts = self.countNeut()
            for i in np.arange(0,self.countPos()):
                posBlock += '+'
            for i in np.arange(0,self.countNeg()):
                negBlock += '-'
            for midNeuts in np.arange(0,nneuts+1):
                midBlock = ''
                for i in np.arange(0,midNeuts):
                    midBlock += '0'
                for startNeuts in np.arange(0,nneuts-midNeuts+1):
                    setupSequence = ''
                    for i in np.arange(0,startNeuts):
                        setupSequence += '0'
                    setupSequence += posBlock
                    setupSequence += midBlock
                    setupSequence += negBlock
                    for i in np.arange(0,nneuts-startNeuts-midNeuts):
                        setupSequence += '0'
                    #DEBUG
                    #print(setupSequence)
                    #assert(len(setupSequence) == self.len)
                    nseq = Sequence(setupSequence)
                    if(nseq.delta()>self.dmax):
                        self.dmax = nseq.delta()
                        #DEBUG
                        #maxSeq = nseq.seq
            #DEBUG
            #print('\n' + maxSeq)
        return self.dmax

    def kappa(self):
        return self.delta()/self.deltaMax()

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

    def toString(self):
        s = "%i\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d" % (self.len,self.fminus(),self.fplus(),self.FCR(),self.NCPR(),self.sigma(),self.delta(),self.deltaMax(),self.kappa(),self.meanHydropathy())
        return s

    def toFileString(self):
        s =  "N:\t%i\n" % (self.len)
        s +=  "f-:\t%d\n" % (self.fminus())
        s += "f+:\t%d\n" % (self.fplus())
        s += "FCR:\t%d\n" % (self.FCR())
        s += "NCPR:\t%d\n" % (self.NCPR())
        s += "Sigma:\t%d\n" % (self.sigma())
        s += "Delta:\t%d\n" % (self.delta())
        s += "Max Delta:\t%d\n" % (self.deltaMax())
        s += "Kappa:\t%d\n" % (self.kappa())
        s += "<H>:\t%d\n" % (self.meanHydropathy())
        s += "Phase Plot Region:" % (self.phasePlotRegion())
        s += "Phase Plot Annotation:" % (self.phasePlotAnnotation())