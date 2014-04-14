# -*- coding: utf-8 -*-
"""
Created on Fri Mar 28 01:09:46 2014

@author: James Ahad
"""

import sequence as s
import numpy as np
import datetime as datetime
import os as os
import wanglandau as wl
import subprocess as sb
import serverVariables as sv

class Client:
    def __init__(self, clientId, email = None):
        self.directory = sv.BASEDIR + ("%d/" % clientId)
        self.clientId = clientId
        self.email = email
        self.sequences = []
        self.wlPID = None
        self.heteroPID = None
        self.chkDirs()

    def chkDirs(self):
        if not os.path.exists(self.directory):
            os.makedirs(self.directory)
            os.makedirs(self.directory + 'Sequence_Parameters')
            os.makedirs(self.directory + 'WangLandau_Results')
            os.makedirs(self.directory + 'PDB_Library')
            self.tstamp = [[0,self.mkTimestamp(0)]]
        else:
            '''
            FIXME: add handling for pre-existing client data
            '''


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
        fileName = self.directory + "allSeqParams.dat"
        f = open(fileName, 'w')
        f.write('N\tF-\tF+\tFCR\tNCPR\tsigma\tdelta\tdelta max\tkappa\t<H>\n')
        for seq in self.sequences:
            f.write(seq.toString()+"\n")
        f.close()

    def genLocalSeqParameterFiles(self):
        baseDir = self.directory + "Sequence_Parameters/"
        for i in range(len(self.sequences)):
            filename = baseDir + ("seq%d.dat" % i)
            f = open(filename, 'w')
            f.write(self.sequences[i].toFileString())
            f.close()

    def WL(self,seqid):
        if(seqid>=len(self.sequences)):
            print("Invalid sequence number. Exiting")
            exit
        pid = os.fork()
        if pid > 0:
            self.wlPID = pid
            self.tstamp.append([1,self.mkTimestamp(1)])
        else:
            wl.WangLandau(self.sequences(seqid))
            os._exit(0)

    def heteroLib(self,seqid):
        if(seqid>=len(self.sequences)):
            print("Invalid sequence number. Exiting")
            exit
        pid = os.fork()
        if pid > 0:
            self.heteroPID = pid
            self.tstamp.append([2,self.mkTimestamp(2)])
        else:
            os.chdir(self.directory + 'PDB_Library')
            sb.call([sv.CAMPARI,'-k',sv.HETEROKEY])
            os._exit(0)


