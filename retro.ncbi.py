#!/usr/bin/python3
import time, sys
sys.dont_write_bytecode = True
from concurrent.futures import ProcessPoolExecutor
from rheast import rheast
from spider import spider

if __name__ == '__main__':
    path = rheast.locatePath()
    ncbi = 'https://www.ncbi.nlm.nih.gov/'
    data = ['NC_001436', 'NC_001488', 'NC_000858', 'NC_001815', 'NC_001802', 'NC_001722', 'NC_001549', 'NC_001482']
    path += '\\seq\\'
    rheast.createPath(path)
    
    pro = ProcessPoolExecutor()
    for i in data:
        pro.submit(spider.ncbiNuccore,path,ncbi,i)
    pro.shutdown()