#!/usr/bin/env python
# coding: utf-8

# # Understanding Chowdhary and Assisi, 2019

# #### To-do list
# 
# 1. Two reciprocally coupled **excitatory** oscillators
# 	- (to observe synchrony)
# 2. Two reciprocally coupled **inhibitory** oscillators
# 	- (to observe asynchrony)
# 3. Bipartite networks
# 	1. No within-group connections
# 	2. Complementary edges are excitatory connections
# 	3. Ratio of cumulative strength of excitation/inhibition is varied
# 4. Balanced sudoku network

# ## Two reciprocally coupled oscillators

# #### Imports/Init

# In[1]:


from brian2 import *


# In[2]:


get_ipython().run_line_magic('matplotlib', 'inline')


# In[3]:


# Constants
tau = 4*ms
I = 1.1/ms
gamma = 1/ms

# Non-dimensionalized constants
tp = tau/ms
curr = I*ms
g = gamma*ms
out = curr/g

# equation (1)
eqs = '''
dv/dt = I - gamma*v : 1
'''

# equation (2) but with t replaced by (phi*tau/gamma) [para equation]
# converted to D.E. by differentiation
# dv/dt = (I/gamma) - v is the actual equation, but tau is introduced to keep the equation dimensionally consistent.
eqs2 = '''
dv/dt = ((I/gamma) - v)/tau : 1
'''


# In[4]:


# Sample first spikes

start_scope()

G1 = NeuronGroup(1, eqs2, threshold='v>1', reset='v = 0', method='exact')
G2 = NeuronGroup(1, eqs2, threshold='v>1', reset='v = 0', method='exact')

G2.v = [0.8]

M1 = StateMonitor(G1, 'v', record=True)
M2 = StateMonitor(G2, 'v', record=True)

run(10*ms)

plot(M1.t/ms, M1.v[0], label='Neuron 1')
plot(M2.t/ms, M2.v[0], label='Neuron 2')
xlabel('Time (ms)')
ylabel('v');


# In[5]:


# plot(M1.t/ms, I/gamma(1-exp(-M1.t/tau)), 'C1--',label='Analytic')
# Plotting the analytic solution

start_scope()

G1 = NeuronGroup(1, eqs2, threshold='v>1', reset='v = 0', method='exact')

M1 = StateMonitor(G1, 'v', record=True)
S1 = SpikeMonitor(G1)

run(20*ms)

# solution to du/dt = (I/gamma - U)/tau is:
# I/gamma * (1 - exp(-t/tau))

plot(M1.t/ms, M1.v[0], label='Neuron 1')
plot(M1.t/ms, (I/gamma)*(1-exp(-M1.t/tau)), 'C1--',label='Analytic')
xlabel('Time (ms)')
ylabel('v');


# In[6]:


# Unconnected neurons

start_scope()

G1 = NeuronGroup(1, eqs2, threshold='v>1', reset='v = 0', method='exact')
G2 = NeuronGroup(1, eqs2, threshold='v>1', reset='v = 0', method='exact')

G2.v = [0.8]

# S1 = Synapses(G1, G2, on_pre='v_post += 0.2')
# S2 = Synapses(G2, G1, on_pre='v_post += 0.2')

# S1.connect()
# S2.connect()

M1 = StateMonitor(G1, 'v', record=True)
M2 = StateMonitor(G2, 'v', record=True)
Sp1 = SpikeMonitor(G1)
Sp2 = SpikeMonitor(G2)

run(50*ms)

plot(M1.t/ms, M1.v[0], label='Neuron 1')
plot(M2.t/ms, M2.v[0], label='Neuron 2')
xlabel('Time (ms)')
ylabel('v');


# In[7]:


# Plotting phase separation between the two neurons over time
# No connections, so stay constant over time

spikes = min(len(Sp2.t[:]), len(Sp1.t[:]))
Vals = Sp1.t[:spikes] - Sp2.t[:spikes]

plot(Sp1.t/ms, Vals/ms)
xlabel('Time')
ylabel('Interspike Intervals')
ylim(-20, 20);


# #### Two reciprocally coupled excitatory oscillators

# In[8]:


# Now with excitatory connections
# Updates with a constant 0.1 addition upon spike.

start_scope()

G1 = NeuronGroup(1, eqs2, threshold='v>1', reset='v = 0', method='exact')
G2 = NeuronGroup(1, eqs2, threshold='v>1', reset='v = 0', method='exact')

G2.v = [0.8]

S1 = Synapses(G1, G2, on_pre='v_post += 0.1')
S2 = Synapses(G2, G1, on_pre='v_post += 0.1')

S1.connect()
S2.connect()

M1 = StateMonitor(G1, 'v', record=True)
M2 = StateMonitor(G2, 'v', record=True)
Sp1 = SpikeMonitor(G1)
Sp2 = SpikeMonitor(G2)

run(100*ms)

plot(M1.t/ms, M1.v[0], label='Neuron 1')
plot(M2.t/ms, M2.v[0], label='Neuron 2')
xlabel('Time (ms)')
ylabel('v');


# In[43]:


# Plotting separations again

spikes = min(len(Sp2.t[:]), len(Sp1.t[:]))
Vals = Sp1.t[:spikes] - Sp2.t[:spikes]

plot(Sp1.t/ms, Vals/ms)
plt.axhline(y = 0, color = 'r', linestyle = '--')
xlabel('Time')
ylabel('Interspike Intervals')
ylim(-10, 10);


# In[10]:


# Analytic solution of eqs2 that we're using: U = (I/gamma)*(1-exp(-M1.t/tau))
# function for this (computes analytic U for given value of t (phi)) => U_forward
# === can i make brian compute this? ===
# need to make a function that gives us the inverse, that is :
# we know U at time point t_-, we update it with an excitatory U(t_-) + e, we need to find the 'timepoint at which U' equals this update.
# i.e. U^-1(U(t_-) + e)
# on_pre='v_post = U_forward(min(U_inverse(v_post + e)), 1))'
# where U_inverse computes the time for U(t_-) + e and returns that time point. 
# inverse of a*(1 - e^(-bx)) is given by:  (1/b) * ln(1/(1 - x/a))
# inverse of (I/gamma)*(1 - e^(-1/tau * t)) will be (tau) * ln(1/(1 - gamma*t/I))


# In[135]:


def U_inverse(v):
        inv = tp * log(1/(1 - g*v/curr))
        return inv


# In[136]:


def U_forward(phi):
        forw = (out)*(1-exp(-phi/tp))
        return forw


# In[54]:


print(U_inverse(M1.v[0]))


# In[133]:


# Now with excitatory connections, according to update rule

start_scope()

exc = 0.1

G1 = NeuronGroup(1, eqs2, threshold='v>1', reset='v = 0', method='exact')
G2 = NeuronGroup(1, eqs2, threshold='v>1', reset='v = 0', method='exact')

G2.v = [0.8]

# S1 = Synapses(G1, G2, on_pre='v_post = U_forward(U_inverse(v_post + exc))')
# S2 = Synapses(G2, G1, on_pre='v_post = U_forward(U_inverse(v_post + exc))')
# Can't figure out how to call these functions to make the equations. Have instead just substituted the formula:
# (i/g)*(1-exp(-(tp * log(1/(1 - g*(v_post + exc)/i)))/tp))
# i/g has been defined as the 'out' variable

S1 = Synapses(G1, G2, on_pre='v_post = out * (1-exp(-(tp * log(1/(1 - g*(v_post + exc)/curr)))/tp))')
S2 = Synapses(G2, G1, on_pre='v_post = out * (1-exp(-(tp * log(1/(1 - g*(v_post + exc)/curr)))/tp))')

S1.connect()
S2.connect()

M1 = StateMonitor(G1, 'v', record=True)
M2 = StateMonitor(G2, 'v', record=True)
Sp1 = SpikeMonitor(G1)
Sp2 = SpikeMonitor(G2)

run(100*ms)

plot(M1.t/ms, M1.v[0], label='Neuron 1')
plot(M2.t/ms, M2.v[0], label='Neuron 2')
xlabel('Time (ms)')
ylabel('v');


# In[134]:


# Plotting separations for excitatory coupling

spikes = min(len(Sp2.t[:]), len(Sp1.t[:]))
Vals = Sp1.t[:spikes] - Sp2.t[:spikes]

plot(Sp1.t/ms, Vals/ms)
plt.axhline(y = 0, color = 'r', linestyle = '--')
xlabel('Time')
ylabel('Interspike Intervals')
ylim(-10, 10);


# #### Two reciprocally coupled inhibitory oscillators

# In[142]:


# Now with inhibitory connections, according to update rule

start_scope()

exc = 0.1

G1 = NeuronGroup(1, eqs2, threshold='v>1', reset='v = 0', method='exact')
G2 = NeuronGroup(1, eqs2, threshold='v>1', reset='v = 0', method='exact')

G2.v = [0.05]

# replace (v_post + exc) with (v_post - exc)

S1 = Synapses(G1, G2, on_pre='v_post = out * (1-exp(-(tp * log(1/(1 - g*(v_post - exc)/curr)))/tp))')
S2 = Synapses(G2, G1, on_pre='v_post = out * (1-exp(-(tp * log(1/(1 - g*(v_post - exc)/curr)))/tp))')

S1.connect()
S2.connect()

M1 = StateMonitor(G1, 'v', record=True)
M2 = StateMonitor(G2, 'v', record=True)
Sp1 = SpikeMonitor(G1)
Sp2 = SpikeMonitor(G2)

run(100*ms)

plot(M1.t/ms, M1.v[0], label='Neuron 1')
plot(M2.t/ms, M2.v[0], label='Neuron 2')
xlabel('Time (ms)')
ylabel('v');


# In[146]:


# Plotting separations for inhibitory coupling

spikes = min(len(Sp2.t[:]), len(Sp1.t[:]))
Vals = Sp1.t[:spikes] - Sp2.t[:spikes]

plot(Sp1.t/ms, Vals/ms)
plt.axhline(y = 0, color = 'r', linestyle = '--')
plt.axhline(y = 11/2, color = 'g', linestyle = '--')
xlabel('Time')
ylabel('Interspike Intervals')
ylim(-10, 10);

