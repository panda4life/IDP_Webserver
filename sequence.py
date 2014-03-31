# -*- coding: utf-8 -*-
"""
Created on Wed Mar 12 11:10:55 2014

@author: jahad
"""

import numpy as np
import random as rng
import time as t
from residues import resTable

lkupTab = resTable('residueData.csv')

class Sequence:
    def __init__(self,seq,dmax = -1):
        self.seq = seq.upper()
        self.len = len(seq)
        chargePattern = []
        for i in np.arange(0,self.len):
            if(lkupTab.lookUpCharge(self.seq[i])>0):
                chargePattern = np.append(chargePattern,1)
            elif(lkupTab.lookUpCharge(self.seq[i])<0):
                chargePattern = np.append(chargePattern,-1)
            else:
                chargePattern = np.append(chargePattern,0)
        self.chargePattern = chargePattern
        self.dmax = dmax #initializing to prevent extra computational time
        self.N_ITERS = 5
        self.N_STEPS = 5000

    def countPos(self):
        return len(np.where(self.chargePattern>0)[0])

    def countNeg(self):
        return len(np.where(self.chargePattern<0)[0])

    def countNeut(self):
        return len(np.where(self.chargePattern==0)[0])

    def Fplus(self):
        return self.countPos()/(self.len+0.0)

    def Fminus(self):
        return self.countNeg()/(self.len+0.0)

    def FCR(self):
        return (self.countPos() + self.countNeg())/(self.len+0.0)

    def NCPR(self):
        return (self.countPos() - self.countNeg())/(self.len+0.0)
    
    def phasePlotRegion(self):
        fcr = self.FCR()
        ncpr = self.NCPR()
        if(fcr < .25 and ncpr<.25):
            return 1
        elif(fcr >= .25 and fcr <= .35 and ncpr <= .35):
            return 2
        elif(fcr > .35 and ncpr <= .35):
            return 3
        elif(fcr > .35 and ncpr > .35):
            if(self.Fplus>.35):
                return 4
            elif(self.Fminus>.35):
                return 5
            else: #This case is impossible but here for completeness
                return None
        else: #This case is impossible but here for completeness
            return None
                
        
    def phasePlotAnnotation(self):
        region = self.phasePlotRegion()
        if(region == 1):
            return 'Globule/Tadpole'
        elif(region == 2):
            return 'Boundary Region'
        elif(region == 3):
            return 'Coils,Hairpins and Chimeras'
        elif(region == 4):
            return 'Negatively Charged Swollen Coils'
        elif(region == 5):
            return 'Positively Charged Swollen Coils'
        else :
            return 'ERROR, NOT A REAL REGION'
        
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
        if(self.countNeut() == self.len):
            return 0
        else:
            return self.NCPR()**2/self.FCR()

    def deltaForm(self,bloblen):
        sigma = self.sigma()
        nblobs = self.len-bloblen+1
        ans = 0
        for i in np.arange(0,nblobs):
            blob = self.chargePattern[i:(i+bloblen)]
            bpos = len(np.where(blob>0)[0])
            bneg = len(np.where(blob<0)[0])
            bncpr = (bpos-bneg)/(bloblen+0.0)
            bfcr = (bpos+bneg)/(bloblen+0.0)
            if(bfcr == 0):
                bsig = 0
            else:
                bsig = bncpr**2/bfcr
            ans += (sigma - bsig)**2/nblobs
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

    def swapRandChargeRes(self):
        rand = rng.Random()
        rand.seed(t.time())
        pickableRes = np.arange(0,self.len)
        swapPair = rand.sample(pickableRes,2)
        while(self.chargePattern[swapPair[0]] == self.chargePattern[swapPair[1]]):
            swapPair = rand.sample(pickableRes,2)
        return self.swapRes(swapPair[0],swapPair[1])

    def toString(self):
        s = "%i\t%3.5f\t%3.5f\t%3.5f\t%3.5f\t%3.5f\t%3.5f\t%3.5f\t%3.5f\t%3.5f" % (self.len,self.Fminus(),self.Fplus(),self.FCR(),self.NCPR(),self.sigma(),self.delta(),self.deltaMax(),self.kappa(),self.meanHydropathy())
        return s

    def toFileString(self):
        s =  "N:\t%i\n" % (self.len)
        s += "f-:\t%3.5f\n" % (self.Fminus())
        s += "f+:\t%3.5f\n" % (self.Fplus())
        s += "FCR:\t%3.5f\n" % (self.FCR())
        s += "NCPR:\t%3.5f\n" % (self.NCPR())
        s += "Sigma:\t%3.5f\n" % (self.sigma())
        s += "Delta:\t%3.5f\n" % (self.delta())
        s += "Max Delta:\t%3.5f\n" % (self.deltaMax())
        s += "Kappa:\t%3.5f\n" % (self.kappa())
        s += "<H>:\t%3.5f\n" % (self.meanHydropathy())
        s += "Phase Plot Region: %i\n" % (self.phasePlotRegion())
        s += "Phase Plot Annotation: %s\n" % (self.phasePlotAnnotation())
        return s