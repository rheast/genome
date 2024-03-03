#!/usr/bin/python3
from urllib.request import urlopen, urlretrieve
from operator import itemgetter
import concurrent.futures
import atexit
import shutil
import time
import os
import re


class RHEast:
    def __init__(self) -> None:
        self.time = time.time()
        self.setPathFolder()
        return

    def getPath(self, path=''):
        if not path:
            path = __file__
        path = os.path.abspath(path)
        return path

    def getPathFolder(self, path=''):
        if not path:
            path = self.getPath()
        path = os.path.dirname(path)
        path += re.findall(r'[\\/]', path)[0]
        return path

    def getPathFiles(self, path=''):
        data = os.listdir(path)
        data = [i for i in data if os.path.isfile(os.path.join(path, i))]
        return data

    def setPathCreate(self, path=''):
        if not os.path.exists(path):
            os.makedirs(path)
        return

    def setPathFolder(self, path=''):
        if not path or not os.path.isdir(path):
            path = self.getPathFolder(path)
        self.path = path
        atexit.register(self.setPathFolderRemove, path+'__pycache__')
        return

    def setPathFolderRemove(self, path=''):
        shutil.rmtree(path, ignore_errors=True)
        return

    def getTime(self):
        data = round(time.time()-self.time, 3)
        print(f'Program running time: {data} s')
        return data

    def runTime(self):
        atexit.register(self.getTime)
        return

    def writeData(self, path='', data=[], style='w'):
        if type(data) != type([]):
            data = [data]
        with open(path, style, encoding='utf-8') as file:
            for i in data:
                file.write(str(i)+'\n')
        return

    def getArrayClean(self, data=[]):
        data = [i for i in data if i]
        return data

    def getArrayPackage(self, data=[], length=1, delete=False):
        data = [data[i:i+length] for i in range(0, len(data), length)]
        if delete and data and len(data[-1]) < length:
            del data[-1]
        return data

    def getArraySorted(self, data=[], style='n'):
        if re.findall('^'+style, 'number'):
            data = sorted(data, key=itemgetter(0))
        if re.findall('^'+style, 'length'):
            data = sorted(data, key=lambda index: len(index[-1]))
        return data

    def getFile(self, path=''):
        data = ''
        with open(path, 'r', encoding='utf-8') as file:
            data = file.read()
        return data

    def downloadFile(self, path='', url=''):
        name = ''
        try:
            name = os.path.basename(path)
            urlretrieve(url, path+'.tmp')
        except Exception as error:
            print('Error:', name, error)
        else:
            os.rename(path+'.tmp', path)
            print('Completed:', name, url)
        return

    def runExecutor(self, function={}, data=[]):
        with concurrent.futures.ProcessPoolExecutor() as executor:
            for i in data:
                executor.submit(function, i)
        return

    def getExecutor(self, function={}, data=[]):
        with concurrent.futures.ProcessPoolExecutor() as executor:
            results = list(executor.map(function, data))
        return results

    def getSeqFasta(self, data='', single=False):
        info, titie = [[]], ''
        data = data.split('\n')
        for i in data:
            i = i.strip()
            if not i:
                titie = ''
                info.append([])
            elif i[0] == '>':
                titie = i
                info.append([titie, ''])
            elif i.isalpha():
                info[-1][-1] += i
        data = self.getArrayClean(info)
        if single:
            data = data[0][-1]
        return data

    def getSeqFastq(self, data='', single=False):
        data = data.split('\n')
        data = self.getArrayPackage(data, 4, True)
        if single:
            data = [i[1] for i in data]
        return data

    def getGenBank(self, data='', style=False):
        data = self.getGenBankGroup(data)
        for i in data:
            if len(data[i]) == 1:
                data[i] = data[i][0]
                continue
            match i:
                case 'ORIGIN':
                    data[i] = self.getGenBankOrigin(data[i])
                case _:
                    data[i] = self.getGenBankInfo(data[i], i)
        if style and re.findall('^'+style, 'translation'):
            data = [e['CDS'] for e in data['FEATURES'] if 'CDS' in e]
            data = [e for e in data if 'translation' in e]
            data = [[e['join'], e['translation']] for e in data]
        return data

    def getGenBankDetail(self, data=''):
        data = data.split(' /')
        info = {}
        for e in data:
            if e.replace('.', '').isdigit():
                e = 'join='+e
            e = e.replace('"', '').replace(')', '')
            e = e.replace('(', '=').split('=')
            if e[0] == 'translation':
                e[1] = e[1].replace(' ', '')
            if len(e) == 1:
                info[e[0]] = e[0]
            else:
                info[e[0]] = e[1]
        return info

    def getGenBankOrigin(self, data=[]):
        info = ''
        for e in data:
            if type(e) == type(''):
                continue
            if e[1]:
                info += e[1]
            else:
                break
        info = info.replace(' ', '')
        return info

    def getGenBankInfo(self, data=[], name=''):
        info, text = [], ''
        for e in data:
            if type(e) == type(''):
                if e.strip():
                    info.append({})
                    text = name
                    info[-1][text] = e
            elif e[0]:
                if name == 'FEATURES':
                    info[-1][text] = self.getGenBankDetail(info[-1][text])
                    info.append({})
                text = e[0]
                info[-1][text] = e[1]
            else:
                info[-1][text] += ' '+e[1]
        if len(info) == 1:
            info = info[0]
        return info

    def getGenBankGroup(self, data=''):
        data = data.split('\n')
        info, title, indent = {}, '', 0
        for i in data:
            j = re.findall(r'(\s*\S+\s+).*', i)
            if j and j[0][0] != ' ':
                text = j[0].strip()
                if text == title:
                    info[title].append(i[indent:])
                else:
                    title = text
                    indent = len(j[0])
                    info[title] = [i[indent:]]
            else:
                info[title].append([i[:indent].strip(), i[indent:]])
        return info

    def getSeqAnticodon(self, data=''):
        data = data.replace('U', 'T')
        info = ''
        for i in data:
            n = '.'
            for e in [['A', 'T'], ['C', 'G']]:
                if i in e:
                    for a in e:
                        if i != a:
                            n = a
            info += n
        return info

    def getCodon(self, data='', single=False):
        data = self.getArrayPackage(data, 3, True)
        codon = [
            ['A', 'GCT', 'GCC', 'GCA', 'GCG'],
            ['R', 'CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
            ['N', 'AAT', 'AAC'],
            ['D', 'GAT', 'GAC'],
            # ['B','AAT','AAC','GAT','GAC'],
            ['C', 'TGT', 'TGC'],
            ['Q', 'CAA', 'CAG'],
            ['E', 'GAA', 'GAG'],
            # ['Z','CAA','CAG','GAA','GAG'],
            ['G', 'GGT', 'GGC', 'GGA', 'GGG'],
            ['H', 'CAT', 'CAC'],
            ['I', 'ATT', 'ATC', 'ATA'],
            ['L', 'CTT', 'CTC', 'CTA', 'CTG', 'TTA', 'TTG'],
            ['K', 'AAA', 'AAG'],
            ['M', 'ATG'],  # True
            ['F', 'TTT', 'TTC'],
            ['P', 'CCT', 'CCC', 'CCA', 'CCG'],
            ['S', 'TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
            ['T', 'ACT', 'ACC', 'ACA', 'ACG'],
            ['W', 'TGG'],
            ['Y', 'TAT', 'TAC'],
            ['V', 'GTT', 'GTC', 'GTA', 'GTG'],
            ['#', 'TAA', 'TGA', 'TAG'],  # False
        ]
        info = ''
        for i in data:
            n = '#'
            for e in codon:
                if i in e:
                    n = e[0]
            if single and n == '#':
                break
            info += n
        return info

    def getCodonPossibility(self, data='', minimum=1):
        info = []
        for i in range(len(data)):
            seq = data[i:]
            if len(seq) < 3:
                break
            if seq[0:3] != 'ATG':
                continue
            seq = self.getCodon(seq, True)
            seq = seq.split('#')[0]
            if len(seq) < minimum:
                continue
            for e in info:
                if re.search(seq, e[1]):
                    seq = None
            if seq:
                seq = [i+1, seq]
                info.append(seq)
        return info


rheast = RHEast()

if __name__ == '__main__':
    pass
