#!/usr/bin/python3
#pip install beautifulsoup4
import json,math,re,os
from concurrent.futures import ProcessPoolExecutor
from urllib.request import urlopen, urlretrieve
from bs4 import BeautifulSoup

path = os.path.split(os.path.realpath(__file__))[0]+'\\seq\\'
if not os.path.exists(path):
    os.makedirs(path)

data = ['NC_001436', 'NC_001488', 'NC_000858', 'NC_001815', 'NC_001802', 'NC_001722', 'NC_001549', 'NC_001482']
ncbi = 'https://www.ncbi.nlm.nih.gov'
report = {
    'fasta': 'fasta',
    'genbank': 'gb'
}

def getNuccore(name):
    html = urlopen(ncbi+'/nuccore/'+name)
    html = BeautifulSoup(html.read(),'html.parser')
    id = html.find('meta',{'name':'ncbi_uidlist'})['content']

    for rep in report:
        a = ncbi+'/sviewer/viewer.cgi?save=file&report='+rep+'&id='+id
        b = path+name+'.'+report[rep]
        urlretrieve(a,b)
        print(name,a,b)

if __name__ == '__main__':
    for i in data:
        ProcessPoolExecutor().submit(getNuccore,i)