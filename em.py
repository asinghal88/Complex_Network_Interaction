#!/usr/bin/env python

from __future__ import print_function,division

from sys import stderr
from numpy import loadtxt,concatenate,max,arange
from random import random
import pandas as pd

n = 113      # Number of nodes
N = 25          # Number of measurements

# Read in the observations
data = loadtxt("out!.txt",int)
m = len(data)   # Number of edges

print(m)

# Set up the initial values
alpha = 0.5*(random()+1)
beta = 0.5*random()
rho = 2*m/(n*(n-1))
nch2 = n*(n-1)//2
Q = dict()

for r in range(100):

    # Calculate Q_ij
    Qmissing = rho*(1-alpha)**N/(rho*(1-alpha)**N + (1-rho)*(1-beta)**N)
    sumE = sumEQ = 0.0
    sumQ = nch2*Qmissing
    for i in range(m):
        uv = (data[i,0],data[i,1])
        E = data[i,2]
        Q[uv] = rho*(alpha**E)*(1-alpha)**(N-E)/ \
                (rho*(alpha**E)*(1-alpha)**(N-E)  \
                + (1-rho)*(beta**E)*(1-beta)**(N-E))
        sumE += E 
        sumEQ += E*Q[uv]
        sumQ += Q[uv] - Qmissing

    # Calculate updated values for parameters
    alpha = sumEQ/(N*sumQ)
    beta = (sumE-sumEQ)/(N*(nch2-sumQ))
    rho = sumQ/nch2
    #print(r,alpha,beta,rho,file=stderr)

# Print out the final network
print(alpha,beta,rho)       
for uv in Q:
      uv[0],uv[1],Q[uv]
k1=pd.DataFrame(list(Q.items()), columns=['node_pair', 'prob'])

k1.to_csv('out2!.csv')
