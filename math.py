H = {'Gag','Pol','Env'}
I = {'ENSE00003838363','ENSE00001316091'}
A = {'androgen','androgen receptor'}
def f(x):
    if x == A:
        return I | H
print(H <= f(A)) #True
print(H == f(A)) #False
