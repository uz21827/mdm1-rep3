#MDM1 Rep 3 - SIR 
import numpy as np
import matplotlib.pyplot as plt

# The three parameters of the SIR model
# Using them as global variables so as not to clutter up the notation
Npeople = 67330000 # population of the UK

gamma = 2.0/10 # estimating about 10 days to recovery
# We'll specify r and deduce beta from beta = gamma*r

# r1 - Unmitigated COVID: r = 3 with no social distancing
# r = 3 
# r2 - Suppressed transmission due to lockdown

r = 1.6

beta =  gamma*r

# Computes all the right-hand side functions
def func_SIR(t,S,I,R):
    gS = - beta*S*I/Npeople
    gI = beta*S*I/Npeople - gamma*I
    gR = gamma*I
    # We return the g functions
    return gS,gI,gR

# Choose a domain of integration
tmax = 365  # days
# Choose a time step (try different values)
dt = 0.02  # i.e.,  a couple of hours
# Number of time steps:
N = int(tmax/dt)

# Allocate the arrays
t = np.zeros((N,))
S = np.zeros((N,))
I = np.zeros((N,))
R = np.zeros((N,))

# Initial condition 1: A handful of infected
# I[0] = 100

# Initial condition 2: 10% of people infected
I[0] = 0.01*Npeople

# Everybody else is healthy but susceptible
S[0] = Npeople - I[0]
# and nobody has recovered and is therefore immune
R[0] = 0

# Implementing the Euler algorithm

'''
The Euler method is a first-order method, which means that the local error (error per step) 
is proportional to the square of the step size, 
and the global error (error at a given time) is proportional to the step size. 
The Euler method often serves as the basis to construct more complex method
'''

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
plt.semilogy(t,S,'-b',lw=2,label='Susceptible')
plt.semilogy(t,I,'-r',lw=2,label='Infected')
plt.semilogy(t,R,'-g',lw=2,label='Recovered')
plt.xlabel('time (days)')
plt.ylabel('Number of people')
plt.legend()
plt.title('SIR model for ' + r'$\beta =$' + str(beta)+ r', $\gamma =$' + str(gamma))
plt.show()

