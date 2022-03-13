#!/usr/bin/python3
import json,math,re,os
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import cpu_count
from operator import itemgetter, attrgetter

class RHEast:
    """ ARRAY """
    def getArrayClean(self,array,*number): # -> "change type {) to []":
        if number and str(number[0]).isdigit():
            number = int(number[0])
        else:
            number = 0
        
        robot = []
        for arr in array:
            if type(arr) == type(1) or len(arr) > number:
                robot.append(arr)
        
        return robot
    #Complete

    def getArrayPackage(self,array,number,*delete): # -> "delete: False or True|whatever, default: False":
        array = [array[i:i+number] for i in range(0,len(array),number)]
        if delete:
            if len(array[-1]) < number and len(array) > 1:
                del array[-1]
        return array
    #Complete

    def getArraySorted(self,array,*style): # -> "style: number[big>small] or length[small>big], default: number":
        if style:
            style = style[0].lower()
        else:
            style = 'n'

        if len(re.findall('^'+style,'number')) > 0:
            array = sorted(array,key=itemgetter(0))[::-1]

        if len(re.findall('^'+style,'length')) > 0:
            array = sorted(array,key=lambda index:len(index[-1]))
        
        return array
    #Complete

    none = {
        ' ': ' ',
        '.': '.',
        '#': '#',
    }
    
    for n in none:
        for i in range(12):
            none[n] += none[n]

    """ PATH """
    def locatePath(self):
        path = os.path.split(os.path.realpath(__file__))[0]
        return path
    #Complete

    def createPath(self,path):
        if not os.path.exists(path):
            os.makedirs(path)
        return
    #Complete

    def getPathName(self,path):
        robot = self.getArrayClean(path.split('\\'))[-1]
        return robot
    #Complete

    def getPathList(self,path,*name): # -> "name: matching path, default: False":
        if not path:
            path = self.locatePath()

        if path[-1] != '\\':
            path += '\\'
        
        if name and type(name[0]) == type([]):
            name = name[0]

        robot = []
        array = os.listdir(path)

        for arr in array:
            if name:
                for n in name:
                    n = re.findall(n,arr)
                    if len(n) > 0:
                        robot.append(path+arr)
            else:
                robot.append(path+arr)
        
        if len(robot) > 1:
            robot = sorted(self.getArrayClean(set(robot)))
        
        return robot
    #Complete

    """ FILE """
    def getFile(self,path,*keep): # -> "keep: False or True|whatever, default: False":
        robot = open(path,'r').read()

        if keep:
            keep = keep[0]
        
        if not keep:
            for s in ['[',']','{','}','(',')',' ','"',"'"]:
                robot = robot.replace(s,'')

        robot = robot.split('\n')

        if not keep:
            for r, rob in enumerate(robot):
                rob = rob.split(',')

                for i, e in enumerate(rob):
                    if e.isdigit():
                        rob[i] = int(e)
                
                robot[r] = self.getArrayClean(rob)

        robot = self.getArrayClean(robot)
        return robot
    #Complete

    def getFileSize(self,path):
        size = os.path.getsize(path)
        print(self.getPathName(path),size)
        return size
    #Complete

    def getFileSequence(self,path,*name):
        style = [
            ['fasta', self.getFasta],
            ['fastq', self.getFastq],
            ['gb', self.getGenBank],
        ]

        if not name:
            style = style[0:2]
        else:
            if type(name[0]) == type([]):
                name = name[0]
        
            for i, arr in enumerate(style):
                if not arr[0] in name:
                    style[i] = []
            
        if type(path) != type([]):
            path = [path]

        style = self.getArrayClean(style)
        sequence = []

        for p in path:
            if sequence:
                break

            for s in style:
                if p and re.findall(s[0],p):
                    sequence = s[1](p)
                    break
        
        return sequence
    #Complete

    def getFasta(self,path,*part):
        sequence = open(path,'r').read()
        sequence = sequence.split('\n')

        if len(part) == 1:
            part = part[0]
        
        if type(part) == type('N'):
            part = part.split('.')
            part = self.getArrayClean(part)

        for i, p in enumerate(part):
            p = p.replace(',','')
            if p.isdigit():
                part[i] = int(int(p)/70)
            else:
                part[i] = ''
        
        part = self.getArrayClean(part)
        if len(part) != 2:
            part = False

        robot = ''
        for s, seq in enumerate(sequence):
            if seq.isalpha():
                if not part or part[0]-2 < s < part[1]+2:
                    robot += seq
            
            if robot:
                if len(seq) == 0 or seq[0] == '>':
                    break
        
        return robot
    #Complete

    def getFastq(self,path):
        sequence = open(path,'r').read()
        sequence = sequence.split('\n')
        sequence = self.getArrayClean(sequence)
        sequence = self.getArrayPackage(sequence,4,True)
        return sequence
    #Complete

    def getGenBank(self,path):
        sequence = open(path,'r').read()
        if sequence[0:5] != 'LOCUS':
            return []
        
        for i, name in enumerate(['FEATURES', 'ORIGIN']):
            sequence = sequence.split(name)

            if len(sequence) < 2:
                sequence = []
                break
            else:
                sequence = sequence[i-1]
        
        if not sequence:
            return []
        
        sequence = sequence.replace(self.none[' '][0:21],'')
        sequence = sequence.split('\n'+self.none[' '][0:5])
        
        for s, seq in enumerate(sequence):
            seq = seq.replace('\n/',' /')
            seq = seq.replace('\n','')

            if seq[0] == ' ':
                sequence[s] = ()
                continue

            array = seq.split(' /')

            for a, arr in enumerate(array):
                if a == 0:
                    arr = arr.split(' ')
                else:
                    arr = arr.split('=')
                
                array[a] = self.getArrayClean(arr)
            
            sequence[s] = array
        
        sequence = self.getArrayClean(sequence)

        for s, seq in enumerate(sequence):
            array = [[0, 0], '']

            for i, arr in enumerate(seq):
                if i == 0:
                    arr = arr[-1].split('(')[-1].split(')')[0]
                    arr = arr.replace('..',',')

                    for t in ['<','>']:
                        arr = arr.replace(t,'')
                    
                    arr = arr.split(',')

                    for n, num in enumerate(arr):
                        if num.isdigit():
                            arr[n] = int(num)
                        else:
                            arr[n] = ''

                    array[0] = self.getArrayClean(arr)
                    continue

                if len(arr) != 2:
                    continue

                if arr[0] == 'translation':
                    array[1] = arr[1].replace('"','')
            
            if not array[1]:
                array = []
            
            sequence[s] = array
        
        sequence = self.getArrayClean(sequence)
        return sequence
    #Complete

    def writeData(self,path,text,*style): # -> "style: 'a' or 'w', default('a')":
        if style:
            style = style[0]
        
        if not style:
            style = 'a'

        if type(text) in [type(''),type(1)]:
            text = [text]

        path = open(path,style,encoding='utf-8')
        
        for t in text:
            path.write(str(t)+'\n')

        path.close()
        return
    #Complete

    def printPercent(self,index,number,percent,*text):
        if text:
            text = text[0]
        
        if not text:
            text = ''
        
        centum = int(100*index/number)
        if centum > percent:
            percent = centum
            print(text,str(percent)+'%')
        
        return percent
    #Complete

    def runProcess(self,array,program):
        process = ProcessPoolExecutor()
        for arr in array:
            run = process.submit(program,arr)
        process.shutdown(wait=True)
        return
    #Complete

    def getProcess(self,program,parameter):
        project = []
        process = ProcessPoolExecutor(cpu_count())

        for i in range(cpu_count()):
            parameter['process'] = i
            run = process.submit(program,**parameter)
            project.append(run)
            #print(i,run)
        
        process.shutdown(wait=True)

        array = []
        for pro in project:
            array += pro.result()
        
        print(len(array))
        return array
    #Complete

    """ SEQUENCE """
    def getSeqCompress(self,sequence): # -> "return: [1|number,N|string]":
        if len(sequence) < 2:
            return sequence
        
        #Get Repeat Times
        robot = {}
        for seq in sequence:
            txt = ()
            num = 1

            if type(seq) == type([]):
                if str(seq[0]).isdigit():
                    num = int(seq[0])
                txt = seq[1]
            
            elif seq.isalpha():
                txt = seq
            
            if txt and len(set(txt)) < 2 < len(txt):
                txt = ()
            
            if txt:
                if txt in robot:
                    robot[txt] += num
                else:
                    robot[txt] = num

        #Json To Array
        robot = str(robot)
        for s in ['{','}',"'",'"',' ']:
            robot = robot.replace(s,'')
            
        robot = robot.split(',')
        robot = self.getArrayClean(robot)

        for r, rob in enumerate(robot):
            rob = rob.split(':')[::-1]
            rob[0] = int(rob[0])
            robot[r] = rob
        
        robot = self.getArraySorted(robot)

        for r, rob in enumerate(robot):
            if rob[0] < 2:
                robot = robot[0:r]
                break
        
        return robot
    #Complete

    def getSeqScope(self,sequence,*number): # -> "number: 1-N, default(2)":
        if number:
            number = number[0]
        
        if type(number) != type(1) or number < 2:
            number = 2
        
        if not sequence:
            return []
        
        for s, seq in enumerate(sequence):
            if seq[0] < number:
                sequence = sequence[0:s]
                break
        
        return sequence
    #Complete

    def getSeqAnticodon(self,sequence,*inverted): # -> "use True or [::-1] to change anticodon to inverted":
        anticodon = ''
        sequence = sequence.replace('U','T')

        for s in sequence:
            if s in ['A','T','C','G']:
                for arr in [['A','T'],['C','G']]:
                    if s in arr:
                        for a in arr:
                            if s != a:
                                anticodon += a
                                break
                        break
            else:
                a = '.'
                if s in ['N','^','$','*','+','?']:
                    a = s
                anticodon += a

        if inverted:
            anticodon = anticodon[::-1]
        
        return anticodon
    #Complete

    def getSeqPosition(self,txt,sequence):
        if '.' in set(txt):
            arr = re.findall(txt,sequence)
            if arr:
                txt = arr[0]

        array = sequence.split(txt)

        for a, arr in enumerate(array):
            if a == 0:
                array[a] = len(arr)
            else:
                array[a] = array[a-1] + len(txt+arr)
        
        del array[-1]
        return array
    #Complete

    def getSeqRepeat(self,hunt,*array):
        array = self.getArrayClean(array)

        if len(array) == 1 and type(array[0]) == type([]):
            array = array[0]

        for a, arr in enumerate(array):
            array[a] = arr.replace(' ','')
        
        array = self.getArrayClean(array)
        if len(array) != 2:
            return
        
        if len(array[0]) > len(array[1]):
            array = array[::-1]
        
        index = 0
        robot = []

        for i in range(len(array[0])):
            if i < index:
                continue
            
            for j in range(len(array[0])-i):
                arr = array[0][i:i+j+hunt]

                if len(arr) < hunt:
                    break

                if re.findall(arr,array[1]):
                    if not arr in robot:
                        robot.append(arr)
                        index = i+len(arr)-hunt
                else:
                    break

        bidden = str(robot)
        for r, rob in enumerate(robot):
            if len(re.findall(rob,bidden)) > 1:
                robot[r] = ''
        
        robot = self.getArrayClean(robot)
        return robot
    #Complete

    def getSeqSimilar(self,hunt,*array):
        if len(array) == 1 and type(array[0]) == type([]):
            array = array[0]

        for a, arr in enumerate(array):
            array[a] = arr.replace(' ','')
        
        array = self.getArrayClean(array)
        if len(array) != 2:
            return
        
        if len(array[0]) > len(array[1]):
            array = array[::-1]
        
        robot = []

        for i in range(len(array[0])):
            reg = array[0][i:i+hunt]

            if i+hunt > len(array[0]) or len(reg) < hunt:
                break
            
            slide = array[1]
            rob = []

            for j in range(len(slide)):
                if j+hunt > len(slide) or len(slide) < hunt:
                    break

                sld = slide[i:i+hunt]
                txt = ''
                num = 0

                for n in range(hunt):
                    if sld[n] == reg[n]:
                        num += 1
                        txt += reg[n]
                    else:
                        txt += '.'
                
                if num >= int(hunt/3):
                    rob.append(txt)
            
            for r in set(rob):
                if robot and r[0:hunt-1] == robot[-1][::-1][0:hunt-1][::-1]:
                    robot[-1] += r[-1]
                else:
                    robot.append(r)
        
        robot = self.getArrayClean(robot)
        return robot
    #Complete

    def getSeqQuality(self,array,*quality):
        if len(array) != 4 or len(array[1]) != len(array[-1]):
            return []
        
        if quality:
            quality = quality[0]
        
        if type(quality) == type(''):
            quality = ord(quality)
        
        if type(quality) != type(1):
            quality = ord('#')
        
        sequence = ''
        for i, e in enumerate(array[1]):
            if ord(array[-1][i]) > quality:
                sequence += e
            else:
                sequence += 'N'
        
        number = len(sequence.replace('N',''))/len(sequence)
        if number < 0.25:
            return []

        array[1] = sequence
        return array
    #Complete

    def getSeqExplore(self,sequence,*number):
        if number and str(number[0]).isdigit():
            number = int(number[0])
        else:
            number = 4
        
        array = self.getArrayPackage(sequence,number,True)
        robot = []

        for i in range(len(array)):
            arr = array[i:i+number]
            if len(arr) < number:
                break

            txt = ''
            for a in arr:
                txt += a

            robot.append(txt)
        
        robot = self.getArrayClean(robot)
        return robot
    #Complete

    def getSeqSplit(self,sequence,hunt):
        prog = []
        for i in range(hunt):
            prog.append(10000*i)

        robot = []
        for s, seq in enumerate(sequence):
            arr = []
            txt = seq[1]
            inv = self.getSeqAnticodon(txt,True)

            if inv == txt:
                arr = self.getArrayPackage(txt,int(len(txt)/2))
                
            else:
                rep = self.getSeqRepeat(hunt,inv,txt)
                for r in rep:
                    txt = txt.replace(r,'#'+r+'#')
                
                if rep:
                    arr = txt.split('#')
                    arr = self.getArrayClean(arr,hunt)
                
            for a in arr:
                robot.append([seq[0],a])
            
            if s in prog:
                print(s,arr)
        
        robot = self.getSeqCompress(robot)
        return robot
    #Complete

    def getSeqMould(self,hunt,array,mould):
        robot = ''
        for a, arr in enumerate(array):
            if type(arr) == type([]):
                arr = arr[-1]
            
            rep = self.getSeqRepeat(hunt,arr,mould)
            for r in rep:
                if len(r) > len(robot):
                    robot = r
            
            if robot:
                break
        
        return robot
    #Complete

    def getSeqMajor(self,sequence):
        major = ''
        number = 0

        for seq in sequence:
            number += seq[0]

        sequence = self.getArraySorted(sequence,'l')

        for i in range(len(sequence)):
            if not sequence:
                break
            
            text = sequence[0][-1]
            match = []
            length = 0
            overage = 0

            for seq in sequence:
                overage += seq[0]
                if re.findall(text,seq[1]):
                    length += seq[0]
                    match.append(seq)
            
            if overage < number/2:
                break
            
            if length > number/2:
                major = text
                sequence = match
            
            sequence = sequence[1:len(sequence)]
        
        return major
    #Complete

    def getSeqOverlay(self,sequence,side):
        overlay = ''
        number = 0

        for seq in sequence:
            number += seq[0]

        for i in range(64):
            robot = []

            for seq in sequence:
                if len(seq[1]) > i:
                    robot.append([seq[0],seq[1][::side*2-1][i]])
            
            robot = self.getSeqCompress(robot)
            
            if robot and robot[0][0] > number/2:
                overlay += robot[0][1]
            else:
                break

        overlay = overlay[::side*2-1]
        return overlay
    #Complete

    def getSeqSplice(self,hunt,array,mould,*data):
        if data:
            data = data[0]
        
        sequence = mould

        for side in range(2):
            for i in range(128):
                mould = sequence[::1-side*2][0:hunt][::1-side*2]
                extend = self.getSeqExtend(hunt,array,mould,side,data)

                if extend and re.findall(extend,sequence):
                    mould = ()
                    print(i,'repeat',extend)

                if extend:
                    sequence = [sequence,extend][::side*2-1]
                    sequence = sequence[0] + sequence[1]
                    print(i,extend)
                else:
                    print(i,False)
                
                if not extend or not mould:
                    break
                
            if data:
                self.writeData(data,[sequence,'\n'])

        return sequence
    #Complete

    def getSeqExtend(self,hunt,array,mould,side,*data):
        if data:
            data = data[0]
        
        if not side in [0,1]:
            print('Side:',False)
            return

        major = ''
        sequence = []

        for arr in array:
            remain = arr[1].split(mould)[side-1]
            remain = arr[1].replace(remain,'')
            remain = remain[::side*2-1][hunt:len(remain)+hunt][::side*2-1]

            if len(remain) > hunt/2:
                sequence.append([arr[0],remain])
            
            if len(sequence) > 10 and arr[0] < 3:
                break
            
            if arr[0] < 2:
                break
        
        sequence = self.getSeqCompress(sequence)

        if sequence and sequence[0][0] > 10:
            major = self.getSeqMajor(sequence)
        else:
            major = self.getSeqOverlay(sequence,side)
                    
        if data:
            self.writeData(data,sequence[0:100]+['',major,'\n'])

        return major
    #Complete

    def getSeqArrange(self,hunt,sequence,mould):
        for i in range(4):
            remain = sequence[0:int(hunt/(i+1))]

            if re.findall(remain,sequence[::-1][0:100][::-1]):
                sequence = sequence.replace(remain,'#'+remain).split('#')
                sequence = self.getArrayClean(sequence)[0]
                sequence += sequence
                break

        for side in range(2):
            for i in range(4):
                remain = mould[::1-side*2][0:int(hunt/(i+1))][::1-side*2]

                if re.findall(remain,sequence):
                    sequence = sequence.replace(remain,'#'+remain+'#').split('#')
                    sequence = self.getArrayClean(sequence)
                    del sequence[-side]

                    robot = ''
                    for seq in sequence:
                        robot += seq
                    sequence = robot
                    break

        print(len(sequence))
        return sequence
    #Complete

    codon = [
        ['A','GCT','GCC','GCA','GCG'],
        ['R','CGT','CGC','CGA','CGG','AGA','AGG',],
        ['N','AAT','AAC',],
        ['D','GAT','GAC',],
        #['B','AAT','AAC','GAT','GAC'],
        ['C','TGT','TGC'],
        ['Q','CAA','CAG',],
        ['E','GAA','GAG'],
        #['Z','CAA','CAG','GAA','GAG'],
        ['G','GGT','GGC','GGA','GGG'],
        ['H','CAT','CAC'],
        ['I','ATT','ATC','ATA',],
        ['L','CTT','CTC','CTA','CTG','TTA','TTG',],
        ['K','AAA','AAG'],
        ['M','ATG'], #True
        ['F','TTT','TTC',],
        ['P','CCT','CCC','CCA','CCG',],
        ['S','TCT','TCC','TCA','TCG','AGT','AGC',],
        ['T','ACT','ACC','ACA','ACG'],
        ['W','TGG',],
        ['Y','TAT','TAC'],
        ['V','GTT','GTC','GTA','GTG'],
        ['#','TAA','TGA','TAG',], #False
    ]

    def getCodon(self,sequence,*reverse):
        if reverse:
            reverse = reverse[0]
        
        robot = '#'
        for c in self.codon:
            if sequence in c:
                if not reverse:
                    robot = c[0]
                else:
                    robot = c[1:len(c)]
                break

        return robot
    #Complete

    def getCodonTranslate(self,array,*end):
        if end:
            end = end[0]

        if type(array) == type('N'):
            array = self.getArrayPackage(array,3,True)

        robot = ''
        for arr in array:
            arr = self.getCodon(arr)
            robot += arr

            if end and arr == '#':
                break

        return robot
    #Complete

    def getCodonPotential(self,sequence,*number):
        if number:
            number = number[0]
        
        if not number:
            number = 0
        
        robot = []
        
        for i in range(len(sequence)-3):
            seq = sequence[i:i+3]

            if seq != 'ATG':
                continue
            
            seq = sequence[i:len(sequence)]
            txt = self.getArrayPackage(seq,3,True)
            txt = self.getCodonTranslate(txt,True)

            if len(txt) < 4:
                continue

            if txt[-1] == '#':
                txt = txt[0:len(txt)-1]

            index = 0
            for rob in robot:
                if type(rob) == type([]):
                    for r in rob:
                        if type(r) == type('') and r[0] == 'M':
                            rob = r
                            break
                
                if re.findall(txt,rob):
                    index += 1
                    break
            
            if index:
                continue

            num = [i+1,i+1+len(txt)*3]
            seq = sequence[num[0]-1:num[1]-1]
            
            if len(txt) > number:
                robot.append([num,txt,seq])

        return robot
    #Complete

    def getCodonExtract(self,hunt,match,protein,sequence,data):
        text = match
        robot = []

        if data['seq']:
            text = self.none['.'][0:3*hunt-len(data['seq'])] + data['seq']
            text = re.findall(text,match)[0]
        
        for i in range(len(text)-3*hunt+1):
            txt = text[i:i+3*hunt]
            i = self.getSeqPosition(txt,match)[0]
            print(i,txt)

            for x in [1,-1]:
                for y in [1,-1]:
                    #if x != -1 or y != -1 or i != 1:
                        #continue

                    seq = self.getArrayPackage(txt,3)[::x]
                    for h in range(hunt):
                        seq[h] = seq[h][::y]

                    seq = self.getCodonTranslate(seq)
                    if re.findall('#',seq):
                        continue

                    if data['codon'] and not re.findall(data['codon'],seq):
                        continue

                    for pro in protein:
                        codon = self.getSeqPosition(seq,pro[1])

                        for c in codon:
                            array = {
                                'y': (i+1)*y,
                                'x': (pro[0][0]+(c+hunt)*3)*x,
                                'type': data['type'],
                                'codon': {
                                    'num': [c+1, c+hunt],
                                    'light': seq,
                                    'extend': re.findall(
                                        str(self.none['.'][0:4]+pro[1][c-1:c+hunt+1]+self.none['.'][0:4]),
                                        str(self.none[' '][0:5]+pro[1]+self.none[' '][0:5])
                                    )[0],
                                },
                                'sequence': {
                                    'num': [pro[0][0]+(c-5)*3, pro[0][0]+(c+hunt+5)*3-1],
                                    'light': sequence[pro[0][0]+c*3-1:pro[0][0]+(c+hunt)*3-1],
                                    'extend': sequence[pro[0][0]+(c-5)*3-1:pro[0][0]+(c+hunt+5)*3-1],
                                },
                                'bind': {
                                    'num': [data['bind'][0]+i,data['bind'][0]+i+hunt*3][::y],
                                    'light': txt[::y],
                                    'extend': txt[::y],
                                },
                                'match': {
                                    'num': [data['bind'][0]-6,data['bind'][1]+6],
                                    'light': txt,
                                    'extend': re.findall(
                                        self.none['.'][0:6]+match+self.none['.'][0:6],
                                        sequence
                                    )[0],
                                },
                            }

                            if not array in robot:
                                robot.append(array)
        
        return robot
    #Complete

    def getCodonAxis(self,data):
        index = 0
        robot = {
            'x': [[],[]],
            'y': [[],[]],
        }
        
        for a in ['x','y']:
            i = int(data[a][1]/data[a][0])
            if i > index:
                index = i

        for a in ['x','y']:
            for i in range(20):
                num = data[a][0]*(i+1)

                if num > data[a][1]:
                    continue
                
                for n in range(2):
                    robot[a][n].append(num*(n*2-1))
            
            robot[a] = robot[a][0] + robot[a][1]

        return robot
    #Complete

    def getNum1(self,number):
        if number > 0:
            number = 1
        else:
            number = -1
        return number
    #Complete

#Terminate