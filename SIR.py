#MDM1 Rep 3 - SIR 
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import csv

'''
Euler Method for solving (numerical method, high inaccuracy),
Susceptible doesnt work properly so needs fixing
'''

f = open('UK_Covid_Data.csv')
csv_reader = csv.reader(f)

R_values = []
#create a list of all the current values, ignoring the first line
for line_no, line in enumerate(csv_reader, 1):
    if line_no == 1:
        pass
    else:
        R_values.append(float(line[1]))

# The three parameters of the SIR model
# Using them as global variables so as not to clutter up the notation
Npeople = 67330000 # population of the UK

gamma = 1/14 #14 days to recover

#r = 2.39 # R value

#beta =  gamma*r

## Computes all the right-hand side functions
def func_SIR(t,S,I,R,gamma,r):
    gS = - ((gamma*r[n+1])*S*I)/Npeople
    gI = (gamma*r[n+1]*S*I/Npeople) - gamma*I
    gR = gamma*I
    # We return the g functions
    print(gI)
    return gS,gI,gR

## Choose a domain of integration
tmax = 365  # days
## Choose a time step (try different values)
dt = 1

## Number of time steps:
N = int(tmax/dt)

## Allocate the arrays
t = np.zeros((N,))
S = np.zeros((N,))
I = np.zeros((N,))
R = np.zeros((N,))

## Initial condition 1: A handful of infected
I[0] = 24

## and nobody has recovered and is therefore immune
R[0] = 68

## Everybody else is healthy but susceptible
S[0] = Npeople - R[0] - I[0]


## Implementing the Euler algorithm
for n in range(0,N-1):
    # Compute all three right hand side functions at t[n] 
    gS,gI,gR = func_SIR(t[n],S[n],I[n],R[n],gamma,R_values[n])
    # Update time, S, I and R according to the Euler rule
    t[n+1] = t[n] + dt
    S[n+1] = S[n] + gS*dt
    I[n+1] = I[n] + gI*dt
    R[n+1] = R[n] + gR*dt
    
# Plot of the results
plt.figure(figsize=(8,6))
plt.plot(t,S,'-b',lw=2,label='Susceptible')
plt.plot(t,I,'-r',lw=2,label='Infected')
plt.plot(t,R,'-g',lw=2,label='Recovered')
plt.xlabel('time (days)')
plt.xlim(0,tmax)
plt.ylim(0,Npeople)
plt.ylabel('Number of people (millions)')
plt.title('SIR model for ' + r'$\beta =$' + str(beta)+ r', $\gamma =$' + str(gamma))
plt.show()
