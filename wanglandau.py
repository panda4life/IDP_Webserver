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
    bincts = binsz/2+np.arange(0,nbins)
    g = [0]*nbins
    H = [0]*nbins
    f = np.exp(1)
    
    #RNG
    rand = rng.Random()
    rand.seed(t.time())
    
    #original sequence information
    gDmax = seq.deltaMax()
    
    #initial starting conditions
    oseq = seq
    idx_old = np.where(abs(bincts-oseq.kappa())<binsz)[0]
    nstep = 0
    niter = 0
    
    while(f > convergence):
        nseq = oseq.swapRandChargeRes()
        nseq.dmax = gDmax
        k = nseq.kappa()
        idx_new = np.where(abs(bincts-k)<binsz)[0]
        if(idx_new < nbins):
            # if new sequence kappa is less visited than old sequence kappa, visit it
            if(np.exp(g[idx_new])<np.exp(g[idx_old])):
                g[idx_new] = g[idx_new] + np.log(f)
                H[idx_new] = H[idx_new] + 1
                oseq = nseq
                nseq = None
                idx_old = idx_new
                idx_new = 0
            # if new sequence kappa is more visited than old sequence kappa, visit it
            # with a given probability (exp(oldvisits - newvisits))
            elif(rand.random()<np.exp(g[idx_old]-g[idx_new])):
                g[idx_new] = g[idx_new] + np.log(f)
                H[idx_new] = H[idx_new] + 1
                oseq = nseq
                nseq = None
                idx_old = idx_new
                idx_new = 0
            # if we decide not to visit at all, we say we saw the place, but never stayed
            else:
                g[idx_new] = g[idx_new] + np.log(f)
                H[idx_new] = H[idx_new] + 1
                nseq = None
                idx_new = 0
            nstep = nstep + 1
            if(np.mod(nstep,100)):
                print(nstep)
            if(nstep == nflatchk):
                print(H)
                nstep = 0
                if len(np.where(H/np.mean(H) >= .7)) == nbins:
                    f = f**.5
                    H = [0]*nbins
                    niter = niter + 1
    return np.vstack(bincts,g)
                