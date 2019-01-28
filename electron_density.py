
from Cell import theta
import math
from Dirak_functions import integral_1_2
import scipy
from scipy import integrate
import sympy
from sympy import *
import matplotlib.pyplot as plt
from progonka import PHI
import numpy as np
T = 0.001
rho =1
N=50000
A = 196.96657
r0 = 1.388*(A/rho)**(1/3)


X=[0]*(N+1)

const = ((2)**(0.5))/((math.pi)**2)
const_e = ((2*theta(T))**1.5)/(2*((math.pi)**2))
def Chem_potential_e():
    return theta(T)*PHI[N]
print(Chem_potential_e())

def x():
    for i in range(1,N + 1):
        X[i] = (i / N) ** 2
x()
X[0] = X[1]/2
H =[0]*(N+1)
def setka():
    for i in range(1, N + 1):
        H[i] = (i/N)
setka()
def int1_2(x):
    if x>10**7:
        return (x**1.5)*(2/3)
    else:
        return integral_1_2(x)

def Vr(i):
    return ((PHI[i])*theta(T)*r0 - Chem_potential_e()*X[i]*r0)
def V(i):
    return Vr(i)/(X[i]*r0)

#def V(i):
#    return (theta(T)*(PHI[i])/X[i] - Chem_potential_e())

betta = 1/T

def electron_density(i):
    return const_e*int1_2((V(i)+Chem_potential_e())/theta(T))

def n_e_full(i):
    return (const / (betta ** 1.5))*integral_1_2(betta*(Chem_potential_e() - V(i)))

n_e_full1 = [0]*(N+1)
for i in range(N):
    n_e_full1[i] = electron_density(i)

print(n_e_full1)
Z_calc1=scipy.integrate.simps(n_e_full1, X)

n_e = [0]*(N+1)
for i in range(N):
    n_e[i] = electron_density(i)

error = [0]*N


j = [0]*N
for i in range(N):
    j[i] += i

G = [0]*(N+1)
for i in range(N):
    G[i] = X[i]*r0
integrand = [0]*(N+1)
for i in range(N):
    integrand[i] = n_e[i]
integrand1 = [0]*(N+1)
for i in range(N):
    integrand1[i] = n_e[i]*4*math.pi*((X[i]*r0)**2)
Z_calc = scipy.integrate.simps(integrand1, G)
print('Заряд =', Z_calc)

with open("integrate.txt","w") as out:
   for i in n_e:
      print(i,file=out)

print(n_e)



h = [0]*(N+1)
def intsimps(x,y):
    I = 0
    I1 = 0
    for i in range(N):
        h[i] = (X[i+1] - X[i])*r0
        I = (h[i]/2)*(y[i+1]+y[i])*(4*math.pi*((X[i]*r0)**2))
        I1 += I
    return I1
intsimps(H, integrand)
Z_calc1 = intsimps(H, integrand)
print('Заряд =', Z_calc1)
print('Заряд =', Z_calc)