#!/usr/bin/python3
import time, sys
sys.dont_write_bytecode = True
from concurrent.futures import ProcessPoolExecutor
from rheast import rheast
from spider import spider

if __name__ == '__main__':
    path = rheast.locatePath()
    ncbi = 'https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos5/sra-pub-zq-14/SRR005/513/'
    sufx = '.sralite.1'
    data = spider.getSraExtend('SRR5513579','SRR5513615')
    
    pro = ProcessPoolExecutor()
    for i in data:
        pro.submit(spider.ncbiSra,path,ncbi,i,sufx)
    pro.shutdown()
