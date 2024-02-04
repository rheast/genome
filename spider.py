#pip install beautifulsoup4
from urllib.request import urlopen, urlretrieve
from bs4 import BeautifulSoup
import os

class Spider:
    def ncbiNuccore(self,path,ncbi,name):
        report = {
            'fasta': 'fasta',
            'genbank': 'gb'
        }

        html = urlopen(ncbi+'nuccore/'+name)
        html = BeautifulSoup(html.read(),'html.parser')
        id = html.find('meta',{'name':'ncbi_uidlist'})['content']

        for rep in report:
            a = ncbi+'sviewer/viewer.cgi?save=file&report='+rep+'&id='+id
            b = path+name+'.'+report[rep]

            if os.path.exists(b):
                print(True)
                continue

            try:
                urlretrieve(a,b+'.tmp')

            except Exception as error:
                print(name,error)

            else:
                os.rename(b+'.tmp',b)
                print(name,a,b)
        return
    #Complete

    def ncbiSra(self,path,ncbi,name,suffix):
        a = ncbi+name+suffix
        b = path

        for i in ['\\sra\\',name+'\\']:
            b += i

            if not os.path.exists(b):
                os.makedirs(b)
        
        b += name
        
        if os.path.exists(b+'.sra'):
            return
        
        try:
            urlretrieve(a,b+'.tmp')

        except Exception as error:
            print(name,error)
            
        else:
            os.rename(b+'.tmp',b+'.sra')
            print(name,True)
        return
    #Complete
    
    def ncbiFastqDump(self,tool,item,output):
        if os.path.exists(item.replace('sra','seq')+'.fastq'):
            return

        os.system(tool+' '+item+' --outdir '+output)
        return
    #Complete

    def getSraExtend(self,*array):
        robot = []

        for i in range(100):
            num = int(array[0].replace('SRR','')) + i
            num = 'SRR'+str(num)
            robot.append(num)

            if num == array[1]:
                break

        print(len(robot),robot)
        return robot
    #Complete

spider = Spider()