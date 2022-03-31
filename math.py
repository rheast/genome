H = {'gag','pol','env'}
I = {'ENSE00003838363','ENSE00001316091'}
I.update(H)
A = {'androgen','androgen receptor'}
def f(x):
    if x == A:
        return I
print(H <= f(A)) #True
print(H == f(A)) #False