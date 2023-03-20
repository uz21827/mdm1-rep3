#MDM1 Rep 3 - SIR 
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

'''
Euler Method for solving (numerical method, high inaccuracy),
Susceptible doesnt work properly so needs fixing
'''
# The three parameters of the SIR model
# Using them as global variables so as not to clutter up the notation
Npeople = 6600000 # population of the UK

gamma = 1/10 # estimating about 10 days to recovery
# We'll specify r and deduce beta from beta = gamma*r

# r1 - Unmitigated COVID: r = 3 with no social distancing
# r = 3 
# r2 - Suppressed transmission due to lockdown

r = 2

beta =  gamma*r

## Computes all the right-hand side functions
def func_SIR(t,S,I,R):
    gS = - (beta*S*I)/Npeople
    gI = (beta*S*I/Npeople) - gamma*I
    gR = gamma*I
    # We return the g functions
    return gS,gI,gR

## Choose a domain of integration
tmax = 365  # days
## Choose a time step (try different values)
dt = 1  # i.e.,  a couple of hours
## Number of time steps:
N = int(tmax/dt)

## Allocate the arrays
t = np.zeros((N,))
S = np.zeros((N,))
I = np.zeros((N,))
R = np.zeros((N,))

## Initial condition 1: A handful of infected
## I[0] = 100

## Initial condition 2: 10% of people infected
I[0] = 0.01*Npeople

## Everybody else is healthy but susceptible
S[0] = Npeople - I[0]

## and nobody has recovered and is therefore immune
R[0] = 0

## Implementing the Euler algorithm
for n in range(0,N-1):
    # Compute all three right hand side functions at t[n] 
    gS,gI,gR = func_SIR(t[n],S[n],I[n],R[n])
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
plt.ylabel('Number of people')
plt.title('SIR model for ' + r'$\beta =$' + str(beta)+ r', $\gamma =$' + str(gamma))
plt.show()




# Total population, N.
N = 66000000
# Initial number of infected and recovered individuals, I0 and R0.
I0, R0 = 1000, 0
# Everyone else, S0, is susceptible to infection initially.
S0 = N - I0 - R0
# Contact rate, beta, and mean recovery rate, gamma, (in 1/days).
beta, gamma = 0.2, 1/10
# A grid of time points (in days)
t = np.linspace(0, 365, 365)

# The SIR model differential equations.
def deriv(y, t, N, beta, gamma):
    S, I, R = y
    dSdt = -beta * S * I / N
    dIdt = beta * S * I / N - gamma * I
    dRdt = gamma * I
    return dSdt, dIdt, dRdt

# Initial conditions vector
y0 = S0, I0, R0
# Integrate the SIR equations over the time grid, t.
ret = odeint(deriv, y0, t, args=(N, beta, gamma))
S, I, R = ret.T

# Plot the data on three separate curves for S(t), I(t) and R(t)

plt.figure(figsize=(8,6))
plt.plot(t,S,'b',lw=2,label='Susceptible')
plt.plot(t,I,'r',lw=2,label='Infected')
plt.plot(t,R,'g',lw=2,label='Recovered')
plt.xlabel('time (days)')
plt.ylabel('Number of people')
plt.xlim(0,365)
plt.ylim(0,N)
plt.title('SIR model for ' + r'$\beta =$' + str(beta)+ r', $\gamma =$' + str(gamma))
plt.show()

