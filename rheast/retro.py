#!/usr/bin/python3
from turtle import title
from .rheast import rheast
import numpy
import re


class Retro:
    def __init__(self) -> None:
        self.direction = [[a, b] for a in [1, -1] for b in [1, -1]]
        self.image = ()
        self.pyplot = ()
        self.patches = ()
        return

    def getRetroMatch(self, sequence='', protein=[], info={}, length=1):
        left, right = info['bind']
        bind = sequence[left-1:right]
        data = []
        for i in range(len(bind)):
            text = bind[i:i+length*3]
            if len(text) < length*3:
                break
            for d in self.direction:
                span = rheast.getArrayPackage(text, 3)
                span = [rheast.getCodon(e[::d[1]]) for e in span]
                span = ''.join(span)[::d[0]]
                for pro in protein:
                    c, codon = pro
                    s = re.search(span, codon)
                    if not s:
                        continue
                    array = {
                        'direction': d,
                        'index': {'x': 0, 'y': i+1},
                        'array': [],
                        'info': info,
                    }
                    for n in range(4):
                        index, light, normal, extend = [], '', '', ''
                        match n:
                            case 0:
                                index = [left, right]
                                light = text
                                normal = bind
                                extend = '.'*15+bind+'.'*30
                                extend = re.search(extend, sequence)[0]
                                extend = extend[i+length*3-1:]
                                extend = extend[:30+length*3]
                            case 1:
                                index = [left+i, left+i+length*3-1]
                                index = index[::d[1]]
                                light = text[::d[1]]
                            case 2:
                                index = list(s.span())
                                light = span
                                extend = codon[index[0]-5:index[1]+5]
                                index[0] += 1
                            case 3:
                                index = c+s.span()[0]*3
                                index = [index, index+length*3-1]
                                light = sequence[index[0]-1:index[1]]
                                extend = sequence[index[0]-15-1:index[1]+15]
                                array['index']['x'] = index[0]
                        arr = {
                            'index': index,
                            'light': light,
                            'normal': normal,
                            'extend': extend,
                        }
                        array['array'].append(arr)
                    data.append(array)
        return data

    def getRetroGroup(self, data=[]):
        info = {}
        for d in self.direction:
            for i, array in enumerate(data):
                for arr in array:
                    if arr['direction'] == d:
                        x, y = d
                        y = y*arr['index']['y']
                        n = f'{x}:{y}'
                        if not n in info:
                            info[n] = {}
                        if not i in info[n]:
                            info[n][i] = []
                        info[n][i].append(arr)
        for n in info:
            if len(info[n]) < len(data):
                info[n] = ()
        data = {}
        for n in info:
            for i in info[n]:
                array = info[n][i]
                if not n in data:
                    data[n] = []
                data[n].append(array)
        return data

    def getRetroCross(self, data=[]):
        info = {}
        for e in data:
            for i in e:
                if not i in info:
                    info[i] = 0
                info[i] += 1
        for d, e in enumerate(data):
            arr, num, seq = {}, {}, {}
            for i in e:
                if info[i] == len(data):
                    if not i in data:
                        arr[i], num[i], seq[i] = [], [], {}
                    for m in e[i]:
                        for n in m:
                            arr[i].append(n)
                            num[i].append(n['index']['x'])
                            s = n['array'][2]
                            s = re.findall('.'+s['light']+'.', s['extend'])[0]
                            if not s in seq[i]:
                                seq[i][s] = []
                            seq[i][s].append(n['index']['x'])
            array = []
            for i in num:
                for s in seq[i]:
                    if len(seq[i][s]) == len(e[i]):
                        num[i] = seq[i][s]
                        break
                while len(num[i]) > len(e[i]):
                    a = numpy.array(num[i])
                    m = int(numpy.mean(a))
                    b = [abs(x-m) for x in a]
                    c = b.index(max(b))
                    del num[i][c]
                for a in arr[i]:
                    x = a['index']['x']
                    if x in num[i]:
                        array.append(a)
            data[d] = array
        return data

    def getRetroRange(self, data=[]):
        x, y = [], []
        for i in data:
            for e in i:
                x.append(e['index']['x']*e['direction'][0])
                y.append(e['index']['y']*e['direction'][1])
        data = [sorted(i) for i in [x, y]]
        return data

    def getRetroMap(self, path='', data=[], info=[]):
        pyplot = self.pyplot
        patches = self.patches
        pyplot.rcParams['font.sans-serif'] = 'Times New Roman'
        figure = pyplot.figure(figsize=self.image['size'])
        for i, array in enumerate(data):
            marker = self.image['marker'][i]
            head = []
            for a in array:
                for n, arr in enumerate(array[a]):
                    color = self.image['color'][n]
                    title = ''
                    x, y = [], []
                    for e in arr:
                        x.append(e['index']['x']*e['direction'][0])
                        y.append(e['index']['y']*e['direction'][1])
                        if not title:
                            j, k = e['info']['name'], e['info']['type']
                            title = j+' '*3+k
                    pyplot.scatter(
                        x, y,
                        marker=marker,
                        color=color,
                        linewidths=len(array[a])-n-1,
                    )
                    if not title in head:
                        head.append(title)
                        pyplot.scatter(
                            [], [],
                            marker=marker,
                            color=color,
                            linewidths=3,
                            label=title,
                        )
        point = self.image['point']
        for i in point:
            point[i] = [a for a in point[i] if a]
        pyplot.xlim(point['x'][0], point['x'][-1])
        pyplot.xticks(point['x'])
        pyplot.ylim(point['y'][0], point['y'][-1])
        pyplot.yticks(point['y'])
        pyplot.legend(loc='upper right', framealpha=1)
        pyplot.gca().spines['top'].set_color('none')
        pyplot.gca().spines['right'].set_color('none')
        pyplot.gca().xaxis.set_ticks_position('bottom')
        pyplot.gca().yaxis.set_ticks_position('left')
        pyplot.gca().spines['bottom'].set_position(('data', 0))
        pyplot.gca().spines['left'].set_position(('data', 0))
        rect = self.image['rect']
        w, h, l = rect['w'], rect['h'], rect['l']
        x, y = info
        rect = patches.Rectangle(
            (x[0]-w, y[0]-h),
            abs(abs(x[0]-w)-abs(x[-1]+w)),
            abs(abs(y[0]-h)-abs(y[-1]+h)),
            linewidth=l,
            edgecolor=self.image['color'][0],
            facecolor='none'
        )
        pyplot.gca().add_patch(rect)
        figure.savefig(path, format='svg', bbox_inches='tight')
        return pyplot


retro = Retro()
