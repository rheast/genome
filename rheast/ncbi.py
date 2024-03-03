#!/usr/bin/python3
from .rheast import rheast
from urllib.request import urlopen, urlretrieve
import re
import os


class NCBI:
    def __init__(self) -> None:
        self.ncbi = 'https://www.ncbi.nlm.nih.gov/'
        self.viewer = self.ncbi+'sviewer/viewer.cgi'
        self.sralite = self.ncbi.replace('www', 'sra-downloadb.be-md')
        self.sralite += 'sos5/sra-pub-zq-14/SRR005/513/'
        self.report = {
            'fasta': 'fasta',
            'genbank': 'gb'
        }
        return

    def getSralite(self, path='', name=''):
        p = f'{path}{name}.sra'
        url = self.sralite+name+'.sralite.1'
        if not os.path.exists(p):
            rheast.downloadFile(p, url)
        return

    def getNuccore(self, path='', name=''):
        uid = ''
        for i in self.report:
            p = f'{path}{name}.{self.report[i]}'
            if os.path.exists(p):
                continue
            if not uid:
                uid = self.getNuccoreUidlist(name)
            url = f'{self.viewer}?save=file&report={i}&id={uid}'
            rheast.downloadFile(p, url)
        return

    def getNuccoreUidlist(self, name=''):
        uid = urlopen(f'{self.ncbi}nuccore/{name}')
        uid = uid.read().decode('utf-8')
        uid = uid.split('ncbi_uidlist')[-1].split('/>')[0]
        uid = re.search(r'\d+', uid)[0]
        return uid


ncbi = NCBI()
