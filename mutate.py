S = (0.592, 0.893, 0.61, 0.511, 0.541, 0.491)
M = (11/189, 11/189, 11/189, 13/189, 13/189, 13/189)
P = 1
for i in range(len(S)):
    P = P * (M[i] ** (3 * (1-S[i])))
print('{:.100f}'.format(P))