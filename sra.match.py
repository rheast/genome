#!/usr/bin/python3
import os, sys
sys.dont_write_bytecode = True
from rheast import rheast

class retrovirus:
    def __init__(self):
        for i in ['run', 'seq', 'sra']:
            p = rheast.locatePath()+'\\'+i+'\\'
            rheast.createPath(p)
            setattr(self,i,p)

        self.fastq = 'SRR.......'
        self.fasta = 'NC_001802'
        self.hunt = 32

        self.runSraCompress()
        self.runSraToNc()
        self.runSraRepeat()
        return
    #Complete

    def runSraRepeat(self):
        sequence = rheast.getPathList(self.run,self.fastq+'.fasta')[0]
        sequence = rheast.getFileSequence(sequence)
        sequence = rheast.getSeqTwin(sequence,12)
        rheast.writeData(self.run+self.fastq+'.log',sequence,'w')

        robot = rheast.getPathList(self.seq,self.fasta+'.fasta')[0]
        robot = rheast.getFileSequence(robot)
        robot = rheast.getSeqTwin(robot,12)
        rheast.writeData(self.run+self.fasta+'.log',robot,'w')
        return
    #Complete

    def runSraToNc(self):
        path = self.run+self.fastq
        if os.path.exists(path+'.fasta'):
            return

        rheast.writeData(path+'.dat',[],'w')
        robot = rheast.getPathList(self.seq,self.fasta+'.fasta')[0]
        title = rheast.getFile(robot,True)[0]

        for i, e in enumerate(title):
            if e == ' ':
                title = title[i:len(title)]
                break

        robot = rheast.getFileSequence(robot)
        sequence = rheast.getPathList(self.run,self.fastq+'.tmp')[0]
        sequence = rheast.getFile(sequence)

        bidden = rheast.getSeqMould(self.hunt,sequence,robot)
        bidden = rheast.getSeqSplice(self.hunt,sequence,bidden,path+'.dat')
        bidden = rheast.getSeqArrange(self.hunt,bidden,robot)
        bidden = rheast.getArrayPackage(bidden,70)
        rheast.writeData(path+'.fasta',['>'+self.fastq+title]+bidden,'w')
        return
    #Complete

    def runSraCompress(self):
        path = self.run+self.fastq
        if os.path.exists(path+'.tmp'):
            return
        
        sequence = rheast.getPathList(self.seq,self.fastq+'.fastq')[0]
        sequence = rheast.getFileSequence(sequence)
        sequence = rheast.getSeqCompress(sequence)
        sequence = rheast.getSeqSplit(sequence,self.hunt)
        rheast.writeData(path+'.tmp',sequence,'w')
        return
    #Complete
#Terminate

retrovirus()