# -*- coding: utf-8 -*-
"""
Created on Mon Mar 31 13:45:04 2014

@author: jahad
"""
import numpy as np
import random as rng
import time as t
from sequence import Sequence
def WangLandau(seq,nbins=10,binmin=0.0,binmax=1.0,nflatchk=10000,convergence=np.exp(.000001),genPermutants=False):
    # WL histogram information
    binsz = (binmax-binmin)/nbins
    bincts = binsz/2+binsz*np.arange(0,nbins)+binmin
    g = [0]*nbins
    H = [0]*nbins
    f = np.exp(1)

    #RNG
    rand = rng.Random()
    rand.seed(t.time())

    #original sequence information #OBSOLETE, Swapping residues stores this info to next sequence
    #gDmax = seq.deltaMax()

    #initial starting conditions
    oseq = seq
    kold = oseq.kappa()
    idx_old = np.argmin(abs(bincts-kold))
    nstep = 0
    niter = 0
    globalStartTime = t.time()
    startTime = t.time()
    while(f > convergence):
        nseq = oseq.swapRandChargeRes()
        knew = nseq.kappa()
        idx_new = np.argmin(abs(bincts-knew))
        #WL Acceptance Criterion
        #http://www.pages.drexel.edu/~cfa22/msim/node53.html
        if(idx_new < nbins):
            acceptProb = min([1,np.exp(g[idx_old]-g[idx_new])])
            # if new sequence kappa is less visited than old sequence kappa, visit it
            if(rand.random() < acceptProb):
                g[idx_new] = g[idx_new] + np.log(f)
                H[idx_new] = H[idx_new] + 1
                oseq = nseq
                kold = knew
                idx_old = np.argmin(abs(bincts-kold))
                nseq = None
                idx_new = 0
            else:
                g[idx_new] = g[idx_new] + np.log(f)
                H[idx_new] = H[idx_new] + 1
                nseq = None
                idx_new = 0
            nstep = nstep + 1
            if(nstep == nflatchk):
                print(H)
                nstep = 0
                if len(np.where(H/np.mean(H) >= .7)) == nbins:
                    f = f**.5
                    H = [0]*nbins
                    niter = niter + 1
                endTime = t.time()
                print(endTime-startTime)
                startTime = t.time()
    globalEndTime = t.time()
    print("Total Run Time in Seconds")
    print(globalEndTime-globalStartTime)
    return np.vstack(bincts,g)
