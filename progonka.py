import math
import numpy as np
import matplotlib.pyplot as plt
from math import gamma
from Changing_parameters import N,T,rho
from Tabular_values import a_0, Na, E_h
from Dirak_functions import integral_1_2, integral_3_2, integral_minus_1_2

from Atom_parameters import Atom_weight,z
from Cell import z_0, r_0, volume, theta
from State_functions import eta, rho_e



# ВСЁ ДЛЯ ЛИНЕАРИЗАЦИИ И ПРОГОНКИ С ИТЕРАЦИЯМИ

s = 30
N=50000

PHI_S = [0] * (s+1)
for i in range(s+1):
    PHI_S[i] = [0] * (N+1)


RESULT = [0] * s
for i in range(s):
    RESULT[i] = [0] * (N+1)

X = [0] * (N + 1)
A = [0] * (N)
B = [0] * (N)
C = [0] * (N)
D = [0] * (N)

ALPHA = [0] * (N + 1)
BETA = [0] * (N + 1)
Y = [0] * (N + 1)

# константа а в уравнении
const = 4 * (2 * theta(T)) ** (1/2) / math.pi * (r_0(rho)) ** 2


# СЕТКА

for i in range(N + 1):
    X[i] = (i / N) ** 2
for i in range(N + 1):

    if i == 0:
        PHI_S[0][0] = z / (theta(T) * r_0(rho))
    else:
        PHI_S[0][i] = z / (theta(T) * r_0(rho)) * (1 - 3 / 2 * X[i] + 1 / 2 * (X[i])**3) - eta(T,rho) * X[i]



# ШАГ СЕТКИ

def h(N):
    return 1 - X[N - 1]


# КОЭФФИЦИЕНТЫ

for i in range(1, N):
    A[i] = 1 / (X[i] - X[i - 1])



for i in range(1, N):
    C[i] = 1 / (X[i + 1] - X[i])


def progonka(s):

    s_current = 0

    while s_current < s:
        for i in range(1, N):
            B[i] = -A[i] - C[i] - const / 4 * (X[i + 1] - X[i - 1]) * integral_minus_1_2(PHI_S[s_current][i] / X[i])

            D[i] = const * (X[i + 1] - X[i - 1]) * (
                            (X[i] / 2 * integral_1_2(PHI_S[s_current][i]/X[i])) - PHI_S[s_current][i] / 4 * integral_minus_1_2(
                        PHI_S[s_current][i] / X[i]))



        for i in range(N - 1, -1, -1):
            if i == N - 1:
                ALPHA[i] = 1 / (1 - h(N) + h(N) ** 2 / 4 * const * integral_minus_1_2(PHI_S[s_current][N]))

                BETA[i] = -h(N) ** 2 / 2 * const * (
                        integral_1_2(PHI_S[s_current][N]) - PHI_S[s_current][i] / 2 * integral_minus_1_2(
                    PHI_S[s_current][N])) * ALPHA[N - 1]
            else:
                ALPHA[i] = -A[i + 1] / (B[i + 1] + C[i + 1] * ALPHA[i + 1])

                BETA[i] = (D[i + 1] - C[i + 1] * BETA[i + 1]) / (B[i + 1] + C[i + 1] * ALPHA[i + 1])



        # ИСКОМОЕ УРАВНЕНИЕ
        for i in range(N + 1):
            if i == 0:
                Y[i] = z / (theta(T) * r_0(rho))
            else:
                Y[i] = ALPHA[i - 1] * Y[i - 1] + BETA[i - 1]

        for i in range(N + 1):
            RESULT[s_current][i] = PHI_S[s_current][i] * theta(T)

        for i in range(N + 1):
            PHI_S[s_current+1][i] = Y[i]



        s_current+=1


        #PHI = [0]*(N+1)
        #for i in range(N+1):
        #    PHI[i] = RESULT4[i]/theta(T)


progonka(s)
fig = plt.figure()
for s_current in range(s):
    graph1 = plt.plot(X, RESULT[s_current])

plt.grid(True)
plt.show()


#K = [0]*(N)
#proverka =[0]*(N)
#for i in range(1,N):
#    K[i] = A[i]*PHI_S[0][i-1]+B[i]*PHI_S[0][i]+C[i]*PHI_S[0][i+1]
#    proverka[i] = D[i] - K[i]
##plt.plot(proverka)
##plt.plot(D)
#plt.plot(K)
#plt.title('коэффициент через уравнение 13')
#plt.grid(True)
#plt.show()

PHI = [0]*(N+1)
PHI1 = [0]*(N+1)

for i in range(N+1):
    PHI1[i] = RESULT[10][i]/theta(T)

for i in range(N+1):
    PHI[i] = RESULT[29][i]/theta(T)

print(PHI[0:100])

