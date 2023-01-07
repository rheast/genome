#!/usr/bin/python3
#pip install matplotlib
import matplotlib.pyplot as pyplot
from rheast import rheast

class retrovirus:
    def __init__(self):
        for i in ['run', 'seq']:
            p = rheast.locatePath()+'\\'+i+'\\'
            rheast.createPath(p)
            setattr(self,i,p)
        
        self.runCodonMatch()
        return
    #Complete

    def runCodonMatch(self):
        rheast.writeData(self.run+'retro.figure.dat',[],'w')
        pyplot.rcParams['font.sans-serif'] = 'Times New Roman'
        figure = pyplot.figure(figsize=(32,6))

        hunt = 3
        image = {
            'color': ['#ee6677', '#eccd00', '#3399ff', '#00c800'],
            'shape': ['s', 'o'],
        }

        data = {
            'HTLV': {
                'NC_001436': {'type': 'HTLV1', 'bind': [406, 423]},
                'NC_001488': {'type': 'HTLV2', 'bind': [766, 783]},
                'NC_000858': {'type': 'STLV1', 'bind': [758, 775]},
                'NC_001815': {'type': 'STLV2', 'bind': [715, 732]},
            },
            'HIV': {
                'NC_001802': {'type': 'HIV1', 'bind': [182, 199]},
                'NC_001722': {'type': 'HIV2', 'bind': [859, 876]},
                'NC_001549': {'type': 'SIV',  'bind': [689, 706]},
                'NC_001482': {'type': 'FIV',  'bind': [358, 375]},
            },
        }

        for n, d in enumerate(data):
            title = []
            number = {}
            bidden = []

            for nc in data[d]:
                sequence = rheast.getPathList(self.seq,nc)
                sequence = rheast.getFileSequence(sequence)

                protein = rheast.getPathList(self.seq,nc)
                protein = rheast.getFileSequence(protein,'gb')

                if not protein:
                    protein = rheast.getCodonPotential(sequence,8)
                
                data[d][nc]['type'] = nc + '   ' + data[d][nc]['type']
                bind = data[d][nc]['bind']
                match = sequence[bind[0]-1:bind[1]]
                
                robot = rheast.getCodonExtract(hunt,match,protein,sequence,data[d][nc])
                rheast.writeData(self.run+'retro.figure.dat',['>'+nc,match,'\n\n'])
                
                title.append(nc + rheast.none[' '][0:3] + data[d][nc]['type'])
                bidden.append(robot)

                if not number:
                    for i in range(len(match)):
                        for y in [1,-1]:
                            number[i*y] = {}
                            for x in [1,-1]:
                                number[i*y][x] = []

                for rob in robot:
                    x = rheast.getOne(rob['x'])
                    y = rob['y']
                    number[y][x].append(nc)

            for y in number:
                for x in number[y]:
                    number[y][x] = list(set(number[y][x]))
                    if len(number[y][x]) < len(data[d]):
                        number[y][x] = False
            
            rheast.writeData(self.run+'retro.figure.dat',['',number,'\n'])

            for i, robot in enumerate(bidden):
                array = {
                    'y': [],
                    'x': [],
                }

                for r, rob in enumerate(robot):
                    x = rheast.getOne(rob['x'])
                    y = rob['y']

                    if number[y][x]:
                        array['x'].append(rob['x'])
                        array['y'].append(rob['y'])
                    else:
                        robot[r] = {}

                bidden[i] = rheast.getArrayClean(robot)
                rheast.writeData(self.run+'retro.figure.dat',[title[i].split(' ')[0]]+bidden[i]+['\n'])

                pyplot.scatter(
                    array['x'], array['y'],
                    marker = image['shape'][n],
                    color = image['color'][i],
                    linewidths = (len(bidden)-i-1)*2,
                )
            
                pyplot.scatter(
                    [], [],
                    marker = image['shape'][n],
                    color = image['color'][i],
                    linewidths = 4,
                    label = title[i],
                )

        array = {
            'x': [500, 10000],
            'y': [2, 10],
        }
        array = rheast.getCodonAxis(array)
        gca = pyplot.gca()

        pyplot.xlim(array['x'][0],array['x'][-1])
        pyplot.xticks(array['x'])
        pyplot.ylim(array['y'][0],array['y'][-1])
        pyplot.yticks(array['y'])
        pyplot.legend(loc='upper right')

        gca.spines['top'].set_color('none')
        gca.spines['right'].set_color('none')
        gca.xaxis.set_ticks_position('bottom')
        gca.yaxis.set_ticks_position('left')
        gca.spines['bottom'].set_position(('data', 0))
        gca.spines['left'].set_position(('data', 0))
        figure.savefig(self.run+'retro.figure.svg')
        return
    #Complete
#Terminate

retrovirus()