#!/usr/bin/python3
# pip install matplotlib
from rheast.rheast import rheast
from rheast.retro import retro
from rheast.ncbi import ncbi


class Analyze:
    def __init__(self) -> None:
        self.nc, self.sra = (), ()
        self.setPath()
        return

    def operate(self):
        match 5:
            case 0:
                # Download NC from NCBI
                self.setDownload('NC')
            case 1:
                # Download SRA from NCBI
                self.setDownload('SRA')
                # Use fasterq-dump
            case 5:
                # Compare sequences
                self.getMatch()
            case 9:
                # Calculate mutation probability
                self.getProbability()
        return

    def setPath(self):
        rheast.setPathFolder(__file__)
        p = rheast.path
        self.run = f'{p}run{p[-1]}'
        self.seq = f'{p}seq{p[-1]}'
        for i in [self.run, self.seq]:
            rheast.setPathCreate(i)
        return

    def setDownload(self, style=''):
        data = []
        match style:
            case 'NC':
                data = [arr['name'] for array in self.nc for arr in array]
                data = [(ncbi.getNuccore, self.seq, i) for i in data]
            case 'SRA':
                data = [(ncbi.getSralite, self.seq, i) for i in self.sra]
        if data:
            rheast.runExecutor(self.runDownload, data)
        return

    def runDownload(self, args):
        f, a, b = args
        f(a, b)
        return

    def getMatch(self):
        data = []
        for array in self.nc:
            array = rheast.getExecutor(self.getMatrix, array)
            array = retro.getRetroGroup(array)
            data.append(array)
        array = retro.getRetroCross(list(data))
        rheast.writeData(self.run+'retro.js', 'data='+str(array))
        import matplotlib.pyplot as pyplot
        import matplotlib.patches as patches
        retro.pyplot = pyplot
        retro.patches = patches
        array = retro.getRetroRange(array)
        retro.getRetroMap(self.run+'retro.svg', data, array)
        return

    def getMatrix(self, data={}):
        path = self.seq+data['name']
        sequence = rheast.getFile(path+'.fasta')
        sequence = rheast.getSeqFasta(sequence, 1)
        protein = rheast.getCodonPossibility(sequence, 256)
        data = retro.getRetroMatch(sequence, protein, data, 3)
        return data

    def getProbability(self):
        S = (0.592, 0.893, 0.61, 0.511, 0.541, 0.491)
        M = (11/189, 11/189, 11/189, 13/189, 13/189, 13/189)
        P = 1
        for i in range(len(S)):
            P = P * (M[i] ** (3 * (1-S[i])))
        P = '{:.64}'.format(P)
        print(P)
        return P


if __name__ == '__main__':
    rheast.runTime()
    analyze = Analyze()
    analyze.nc = [
        [
            {'name': 'NC_001436', 'type': 'HTLV-1', 'bind': [406, 423]},
            {'name': 'NC_001488', 'type': 'HTLV-2', 'bind': [766, 783]},
            {'name': 'NC_000858', 'type': 'STLV-1', 'bind': [758, 775]},
            {'name': 'NC_001815', 'type': 'STLV-2', 'bind': [715, 732]},
        ],
        [
            {'name': 'NC_001802', 'type': 'HIV-1', 'bind': [182, 199]},
            {'name': 'NC_001722', 'type': 'HIV-2', 'bind': [859, 876]},
            {'name': 'NC_001549', 'type': 'SIV',   'bind': [689, 706]},
            {'name': 'NC_001482', 'type': 'FIV',   'bind': [358, 375]},
        ],
    ]
    analyze.sra = [f'SRR{i}' for i in list(range(5513579, 5513615))]
    retro.image = {
        'color': ['#ee6677', '#eccd00', '#3399ff', '#00c800'],
        'marker': ['s', 'o'],
        'point': {
            'x': list(range(-9000, 9000+1, 1000)),
            'y': list(range(-10, 10+1, 2)),
        },
        'size': [16, 5],
        'rect': {
            'w': 400,
            'h': 1.5,
            'l': 1.5,
        }
    }
    analyze.operate()
