#MDM1 Rep 3 - SIR 
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# The three parameters of the SIR model
# Using them as global variables so as not to clutter up the notation
Npeople = 67330000 # population of the UK

gamma = 1.0/14 # estimating about 10 days to recovery
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
# Choose a time step
dt = 0.02  
# Number of time steps:
N = int(tmax/dt)

# Allocate the arrays
t = np.zeros((N,))
S = np.zeros((N,))
I = np.zeros((N,))
R = np.zeros((N,))

# Initial condition 1: A handful of infected
#I[0] = 100

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

for n in range(0,n-1):
    # compute all three right hand side functions at t[n] 
    gs,gi,gr = func_sir(t[n],s[n],i[n],r[n])

    # update time, s, i and r according to the euler rule
    t[n+1] = t[n] + dt
    s[n+1] = s[n] + gs*dt
    i[n+1] = i[n] + gi*dt
    r[n+1] = r[n] + gr*dt
    
# plot of the results
plt.figure(figsize=(8,6))
plt.semilogy(t,s,'-b',lw=2,label='susceptible')
plt.semilogy(t,i,'-r',lw=2,label='infected')
plt.semilogy(t,r,'-g',lw=2,label='recovered')
plt.xlabel('time (days)')
plt.ylabel('number of people')
plt.legend()
plt.title('sir model for ' + r'$\beta =$' + str(beta)+ r', $\gamma =$' + str(gamma))
plt.show()




'''
Notes:
Suseptible doesnt seem to work very well, not sure why, Euler method not great
Look into using odeint from scripy.integrate instead as is more an "exact" result
Keep both methods though to show progress etc etc
'''
