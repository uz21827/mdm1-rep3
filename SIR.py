import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import csv

f = open('UK Covid Data.csv')
csv_reader = csv.reader(f)

# read in R values from CSV file
R_values_raw = []
with open('UK Covid Data.csv', 'r') as f:
    csv_reader = csv.reader(f)
    for i in range(0,244):
        next(csv_reader)  # skip header row
    R_values_raw = [float(row[3]) for row in csv_reader]

Npeople = 67330000.0 # population of the UK
gamma = 1/14.0 #14 days to recover
tmax = 820  # days
dt = 1 #time step - 1 day
N = int(tmax/dt) # umber of time steps


def func_SIR(t, S, I, gamma, R_current):
    gS_local = - ((gamma * R_current) * S * I) / Npeople
    gI_local = ((gamma * R_current * S * I) / Npeople) - gamma * I
    gR_local = gamma * I
    return gS_local, gI_local, gR_local

# Allocate the arrays
t = np.linspace(0, tmax, N)
S = np.zeros((N,))
I = np.zeros((N,))
R = np.zeros((N,))

# Initial condition 1: A handful of infected
I[0] = 20349

# and nobody has recovered and is therefore immune
R[0] = 451718.0

# Everybody else is healthy but susceptible
S[0] = Npeople - R[0] - I[0]

# Implementing the Euler algorithm
for n in range(0, N-1):
    # Compute all three right hand side functions at t[n] 
    gS, gI, gR = func_SIR(t[n], S[n], I[n], gamma, R_values_raw[n])
    # Update time, S, I and R according to the Euler rule
    S[n+1] = S[n] + gS * dt
    I[n+1] = I[n] + gI * dt
    R[n+1] = R[n] + gR * dt
    t[n+1] = t[n] + dt
    print("R: ", R_values_raw[n],"    Infected:", I[n],"    Day:", t[n])
   
    
# Plot of the results
plt.figure(figsize=(8,6))
plt.semilogy(t,S,'-b',lw=2,label='Susceptible')
plt.semilogy(t,I,'-r',lw=2,label='Infected')
plt.semilogy(t,R,'-g',lw=2,label='Recovered')
plt.xlabel('time (days)')
plt.xlim(0,tmax)
plt.ylim(0,Npeople)
plt.ylabel('Number of peope')
plt.title('SIR model for UK Covid Data: 01/10/2020 - 02/01/2023')
plt.show()
