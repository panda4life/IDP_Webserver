# -*- coding: utf-8 -*-
"""
Created on Wed Mar 12 09:58:32 2014

@author: jahad
"""
import csv
import numpy as np

class Residue:
    def __init__(self,name,letterCode3,letterCode1,hydropathy,charge):
        self.name = name
        self.letterCode3 = letterCode3
        self.letterCode1 = letterCode1
        self.hydropathy = hydropathy
        self.charge = charge

class resTable:
    def __init__(self,filename):
        self.filename = filename
        self.resTab = []
        with open(self.filename, 'rb') as f:
            reader = csv.reader(f)
            for row in reader:
                res = Residue(row[0],row[1],row[2],float(row[3]),int(row[4]))
                self.resTab = np.append(self.resTab,res)
    
    def lookForRes(self,resCode,codeType=1):    
        if(codeType == 1):
            return next((x for x in self.resTab if x.letterCode1 == resCode), None)
        elif(codeType == 3):
            return next((x for x in self.resTab if x.letterCode3 == resCode), None)
        else:
            return None
    
    def lookUpHydropathy(self,resCode,codeType=1):
        res = self.lookForRes(resCode,codeType)
        if(res == None):
            print('Illegal residue code or codeType\nLegal code types are 1 and 3\nResidue codes must be the corresponding 1 letter or 3 letter code for a given residue\n')
            return None
        return self.lookForRes(resCode,codeType).hydropathy        
    
    def lookUpCharge(self,resCode,codeType=1):
        #first check for reduced residue types (aka charge +/-/0)
        if(resCode == '+'):
            return 1
        elif(resCode == '-'):
            return -1
        elif(resCode == '0'):
            return 0
        else:    
            res = self.lookForRes(resCode,codeType)
            if(res == None):
                print('Illegal residue code or codeType\nLegal code types are 1 and 3\nResidue codes must be the corresponding 1 letter or 3 letter code for a given residue\n')
                return None
            return self.lookForRes(resCode,codeType).charge


    
    
    
    
