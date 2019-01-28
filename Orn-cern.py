#from scratch import RESULT3
#from scratch import PHI_S
from Cell import theta
import math
from Dirak_functions import integral_1_2
import scipy
from scipy import integrate
import sympy
from sympy import *
from scratch import h
import matplotlib.pyplot as plt
from progonka import PHI

T=0.001
Z = 79
R = 3.63
N=100
n = 100
r = 0.5
rho = 1
A = 196
r0 = 1.388*(A/rho)**(1/3)
k = 10000
h = 1/N
X=[0]*(N+1)
theta1 = 36.75*T
def x():
    for i in range(1,N + 1):
        X[i] = (i / N) ** 2
x()
X[0] = X[1]/2

n_e_full = [0]*(N+1)
#def Chem_potential_e():
#    return (RESULT3[N])
def Chem_potential_e():
    return -1.33*(10**(-2))
def n_e_pribl():
    return (((2*theta(T))**1.5)/(2*(math.pi)**2))*integral_1_2((Chem_potential_e())/theta(T))
u = n_e_pribl()
Z0 = u * (4/3) * math.pi * ((r0)**3)
print(((-Chem_potential_e())/theta(T)))
print(r0)
print(theta(T))
print(Z0)
#print(RESULT3[0])
#print(Chem_potential_e())
#Chem_potential_e =9.780787802584433
#print(Chem_potential_e())
Chem_potential_e_preved = -Chem_potential_e()/theta(T)
def Vr(i):
    return (PHI[i] + Chem_potential_e_preved*X[i]*r0)*theta(T)
def V(i):
    return Vr(i)/(X[i]*r0)
def g(r):
    return 1/(1+(math.e**(k*(r-R))))

const = ((2)**(0.5))/((math.pi)**2)
betta = 1/T
def n_e_full1(i):
    return (((2*theta(T))**1.5)/(2*(math.pi)**2))*integral_1_2((V(i)+Chem_potential_e())/theta(T))
print(n_e_full1(1))
for i in range(N):
    n_e_full[i] = n_e_full1(i)*4*math.pi*((X[i]*r0)**2)
#print(n_e_full)
#print(PHI_S[3][1])
#print(V(1))
#print(Vr(1))
#print(betta*(Chem_potential_e()-V(i)))
#print(integral_1_2(betta*(Chem_potential_e()-V(1))))
Z_calc = scipy.integrate.simps(n_e_full, X)
def mean(numbers):
    return float(sum(numbers)) / max(len(numbers), 1)

print('Заряд:',Z_calc)
c1 = (((2*theta(T))**1.5)/(2*(math.pi)**2))
rho = [0]*(N+1)
def rho_e(i):
    return c1*integral_1_2((V(i)+Chem_potential_e())/theta(T))
for i in range(N):
    rho[i] = rho_e(i)
plt.plot(X,rho)
plt.show()
def V_Ne_eff(r):
    f = lambda x: ((n_e_full1 - n_e_0() * g(x)) / (math.fabs(r - x-0.00001)))
    y, err = scipy.integrate.quad(f,0,R)
    return (-Z/r + y)
def n_e_full(r):
    return (const/(betta**1.5))*integral_1_2(betta*(Chem_potential_e()-V_Ne_eff(r)))
def n_e_0():
    return (const/(betta**1.5))*integral_1_2(betta*(Chem_potential_e()))
    #n_e_0 = limit(n_e_full(r),r,0)
F = lambda x: n_e_full(r=x)
Z_calc = scipy.integrate.quad(F,0,3.63)
print(n_e_full(r=1))
G = lambda x: n_e_full(r=x)
#Z,err1 = scipy.integrate.quad(G,0,3.63)
#V_Ne_eff(1)
#n_e_full1 = n_e_full(0.5)
#V_Ne_eff(1)
#n_e_full(1)
#Z,err1 = scipy.integrate.quad(G,0,3.63)
#print(Z)
#print(V_Ne_eff(1), n_e_full(1))
#print(Chem_potential_e())
#print(1)
#while Z_calc != Z:
#   V_Ne_eff(r)
#    n_e_full(r)
#    Z_calc = scipy.integrate(n_e_full(r), 0, R)




