# -*- coding: utf-8 -*-
"""
Created on Fri Mar 28 01:09:46 2014

@author: James Ahad
"""

import sequence as s
import numpy as np
import datetime as datetime
import os as os

class Client:
    def __init__(self, directory, clientId = None, email = None):
        self.directory = directory
        self.clientId = clientId
        self.email = email

        self.mkDirs()
        self.sequences = []
        self.tstamp = self.mkTimestamp(0)

    def mkDirs(self):
        if not os.path.exists(self.directory):
            os.makedirs(self.directory)
        os.makedirs(self.directory + 'Sequence_Parameters')
        os.makedirs(self.directory + 'WangLandau_Results')
        os.makedirs(self.directory + 'PDB_Library')

    def mkTimestamp(self, event):
        currTime = datetime.datetime.now()
        fileName = self.directory + "tstamp.dat"
        if(event == 0): #Client just connected
            f = open(fileName, 'w')
        else:
            f = open(fileName, 'a')
        f.write('%s\t%i\n'%(str(currTime),event))
        f.close()
        return currTime

    def addSequence(self, seq):
        self.sequences = np.append(self.sequences,s.Sequence(seq))

    def genSeqIDFile(self):
        fileName = self.directory + "sequences.dat"
        f = open(fileName, 'w')
        f.write('SeqID\tResidue Sequence\n')
        for i in np.arange(0,len(self.sequences)):
            f.write('%i\t%s\n'%(i,self.sequences[i].seq))
        f.close()

    def genGlobalSeqParameterFile(self):
        fileName = self.directory + "sequences.dat"
        f = open(fileName, 'w')
        f.write('SeqID\tResidue Sequence\n')
        for i in np.arange(0,len(self.sequences)):
            f.write('%i\t%s\n'%(i,self.sequences[i].seq))
        f.close()

    def genLocalSeqParameterFiles(self):
        blah

    def WL(self):
        blah


