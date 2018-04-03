%matplotlib inline
from __future__ import division
import numpy as np
from numpy.random import rand
import matplotlib.pyplot as plt

def initialstate(N):   
    ''' generates a random configuration of -1 and +1 spins for initial condition'''
    state = 2*np.random.randint(2, size=(N,N))-1
    return state


def mcmove(config, beta):
    '''Monte Carlo move using Metropolis algorithm '''
    for i in range(N):
        for j in range(N):
                a = np.random.randint(0, N) # picks a random x-coordinate on the lattice
                b = np.random.randint(0, N) # picks a random y-coordinate on the lattice
                s =  config[a, b]           # centers on that random point
                nb = config[(a+1)%N,b] + config[a,(b+1)%N] + config[(a-1)%N,b] + config[a,(b-1)%N] + config[(a+2)%N,b] + config[a,(b+2)%N] + config[(a-2)%N,b] + config[a,(b-2)%N]
                cost = s*nb                 # calculates deltaE for the change of spin by multiplying value of chosen spin by sum of all neighbouring spins
                if cost < 0:                # if cost is less than 0 then spin switches
                    s *= -1
                elif rand() < np.exp(-cost*beta):   #if cost greater than 0 then set probability of spin switching
                    s *= -1
                config[a, b] = s                           # sets new spin value whether switched or not
    return config


def calcEnergy(config):
    '''Total Energy of a given configuration'''
    energy = 0
    for i in range(len(config)):                           # calculates energy for each spin then sums over all spins in lattice
        for j in range(len(config)):
            S = config[i,j]
            nb = config[(i+1)%N, j] + config[i,(j+1)%N] + config[(i-1)%N, j] + config[i,(j-1)%N] + config[(i+2)%N, j] + config[i,(j+2)%N] + config[(i-2)%N, j] + config[i,(j-2)%N]
            energy += -nb*S/2.0
    return energy


def calcMag(config):
    '''Total Magnetization of a given configuration'''
    mag = np.sum(config)
    return mag
    
    nt      = 2**8        # number of temperature points
N       = 2**4        # size of the lattice, N x N
eqSteps = 2**10       # number of MC sweeps for equilibration
mcSteps = 2**10       # number of MC sweeps for calculation

n1, n2  = 1.0/(mcSteps*N*N), 1.0/(mcSteps*mcSteps*N*N)
tc = 2.269            #critical temperature
T = np.random.normal(tc, .64, nt)
T  = T[(T>1.2) & (T<3.8)]
nt = np.size(T)                                                                                                                                                                                                                                                                                                                                   

Energy         = np.zeros(nt)
Magnetization  = np.zeros(nt)
SpecificHeat   = np.zeros(nt)

for m in range(len(T)):
    E1 = M1 = E2 = M2 = 0
    config = initialstate(N)
    iT = 1.0 / T[m]; iT2 = iT * iT;
    
    for i in range(eqSteps):         # equilibrate the system
        mcmove(config, iT)           # Monte Carlo moves

    for i in range(mcSteps):
        mcmove(config, iT)           
        Ene = calcEnergy(config)     # calculate the energy
        Mag = calcMag(config)        # calculate the magnetisation

        E1 = E1 + Ene
        M1 = M1 + Mag
        M2 = M2 + Mag * Mag 
        E2 = E2 + Ene * Ene

        Energy[m]         = n1 * E1
        Magnetization[m]  = n1 * M1
        SpecificHeat[m]   = (n1 * E2 - n2 * E1 * E1) * iT2
       
f = plt.figure(figsize=(18, 10)); # plot the calculated values    

sp =  f.add_subplot(2, 2, 1 );
plt.plot(T, Energy, 'o', color = "#A60628");
plt.xlabel("Temperature (T)", fontsize = 20);
plt.ylabel("Energy ", fontsize = 20);

sp =  f.add_subplot(2, 2, 2 );
plt.plot(T, abs(Magnetization), 'o', color = "#348ABD");
plt.xlabel("Temperature (T)", fontsize = 20);
plt.ylabel("Magnetization ", fontsize = 20);

sp =  f.add_subplot(2, 2, 3 );
plt.plot(T, SpecificHeat, 'o', color = "#A60628");
plt.xlabel("Temperature (T)", fontsize = 20);
plt.ylabel("Specific Heat ", fontsize = 20);

