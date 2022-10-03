#!/usr/bin/env python
# coding: utf-8

# ## Bipartite networks from Chowdhary and Assisi, 2019

# #### To-do list
# 3. Bipartite networks
# 	1. No within-group connections
# 		- Only inhibitory
# 		- Only excitatory
# 	2. Complementary edges are excitatory connections
# 		- (Mix of inhibitory and excitatory)
# 	3. Ratio of cumulative strength of excitation/inhibition is varied
# 		- Varied by changing p (changing number of connections of each type)
# 		- Varied by changing strengths of each type
# 
# 
# - Bipartite network
# 	- Is there synchrony between neurons of the same group, even without excitatory connections?
# 	> Our simulations showed that while connected neurons did
# 	> not fire together, nodes from the same partition did not
# 	> fire synchronously either.
# - Sudoku network
# - Vogels' model
# 	- What kind of connections between inhibitory and excitatory groups (E-I, I-I, I-E)

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
out = curr/g  #should be > the threshold of firing

# equation (1)
eqs = '''
dv/dt = I - gamma*v : 1
'''

# equation (2) 
eqs2 = '''
dv/dt = ((I/gamma) - v)/tau : 1
'''


# ## No within-group connections

# #### Inhibitory connections only

# In[4]:


# From starter.ipynb
# Init all 100 neurons with random starting V, the groups G1 and G2 are connected with probability of 0.6 as in the paper
# (only inhibitory connections across the groups, no other connections)

start_scope()

exc = 0.001
n = 100

np.random.seed(10)

G1 = NeuronGroup(n, eqs2, threshold='v>1', reset='v = 0', method='euler')
G2 = NeuronGroup(n, eqs2, threshold='v>1', reset='v = 0', method='euler')

group1init = np.random.randint(0, 9, size=n)/10
group2init = np.random.randint(0, 9, size=n)/10

G1.v = group1init
G2.v = group2init

# G1.v = 'rand()'
# G2.v = 0.


# In[5]:


S1 = Synapses(G1, G2, on_pre='v = out * (1-exp(-(tp * log(1/(1 - g*(v_post - exc)/curr)))/tp))')
S2 = Synapses(G2, G1, on_pre='v = out * (1-exp(-(tp * log(1/(1 - g*(v_post - exc)/curr)))/tp))')

# S1 = Synapses(G1, G2, on_pre='v_post = v_post - exc')
# S2 = Synapses(G2, G1, on_pre='v_post = v_post - exc')

p = 0.6

S1.connect(p = p)
S2.connect(p = p)

M1 = StateMonitor(G1, 'v', record=True)
M2 = StateMonitor(G2, 'v', record=True)
Sp1 = SpikeMonitor(G1)
Sp2 = SpikeMonitor(G2)

run(500*ms)


# In[113]:


plot(Sp1.t/ms, Sp1.i, '.b', markersize=5)
plot(Sp2.t/ms, Sp2.i, '.r', markersize=5)
xlabel('Time (ms)')
ylabel('Neuron index')
text(150, 110, 'e = {}   |   n = {}   |   p = {}'.format(exc, n, p))
show()

# !neurons of the same indices across the two groups are not necessarily connected to each other!


# Phase asynchrony is observed.

# #### Excitatory connections only

# In[85]:


# Init all 100 neurons with random starting V, the groups G1 and G2 are connected with probability of 0.6 as in the paper
# (only excitatory connections across the groups, no other connections)

start_scope()

exc = 0.001
n = 100

np.random.seed(10)

G1 = NeuronGroup(n, eqs2, threshold='v>1', reset='v = 0', method='euler')
G2 = NeuronGroup(n, eqs2, threshold='v>1', reset='v = 0', method='euler')

group1init = np.random.randint(0, 9, size=n)/10
group2init = np.random.randint(0, 9, size=n)/10

G1.v = group1init
G2.v = group2init


# In[86]:


S1 = Synapses(G1, G2, on_pre='v = out * (1-exp(-(tp * log(1/(1 - g*(v_post + exc)/curr)))/tp))')
S2 = Synapses(G2, G1, on_pre='v = out * (1-exp(-(tp * log(1/(1 - g*(v_post + exc)/curr)))/tp))')

# S1 = Synapses(G1, G2, on_pre='v_post = v_post - exc')
# S2 = Synapses(G2, G1, on_pre='v_post = v_post - exc')

p = 0.6

S1.connect(p = p)
S2.connect(p = p)

M1 = StateMonitor(G1, 'v', record=True)
M2 = StateMonitor(G2, 'v', record=True)
Sp1 = SpikeMonitor(G1)
Sp2 = SpikeMonitor(G2)

run(500*ms)


# In[87]:


plot(Sp1.t/ms, Sp1.i, '.b', markersize=5)
plot(Sp2.t/ms, Sp2.i, '.r', markersize=5)
xlabel('Time (ms)')
ylabel('Neuron index')
text(150, 110, 'e = {}   |   n = {}   |   p = {}'.format(exc, n, p))
show()

# !neurons of the same indices across the two groups are not necessarily connected to each other!


# Phase synchrony is observed.

# #### Phase Space plots

# In[136]:


polared = lambda t, M : M.v[:,t] * 2*pi
dummy = ones(n)*5
alph = 0.2
polar(polared(0, M1), dummy, 'bo', alpha = alph)
polar(polared(0, M2), dummy, 'ro', alpha = alph)


# In[137]:


polar(polared(10, M1), dummy, 'bo', alpha = alph)
polar(polared(10, M2), dummy, 'ro', alpha = alph)


# In[138]:


polar(polared(100, M1), dummy, 'bo', alpha = alph)
polar(polared(100, M2), dummy, 'ro', alpha = alph)


# In[139]:


polar(polared(210, M1), dummy, 'bo', alpha = alph)
polar(polared(210, M2), dummy, 'ro', alpha = alph)


# In[140]:


polar(polared(400, M1), dummy, 'bo', alpha = alph)
polar(polared(400, M2), dummy, 'ro', alpha = alph)


# In[141]:


polar(polared(450, M1), dummy, 'bo', alpha = alph)
polar(polared(450, M2), dummy, 'ro', alpha = alph)


# In[142]:


polar(polared(430, M1), dummy, 'bo', alpha = alph)
polar(polared(430, M2), dummy, 'ro', alpha = alph)


# In[143]:


polar(polared(490, M1), dummy, 'bo', alpha = alph)
polar(polared(490, M2), dummy, 'ro', alpha = alph)


# In[144]:


# should plot phase space at the start vs towards the end
# starts distributed across phase space
# becomes clustered


# In[145]:


# ISI Plotter (needs to be modded for this version with multiple neurons)
# spikes = min(len(Sp2.t[:]), len(Sp1.t[:]))
# Vals = Sp1.t[:spikes] - Sp2.t[:spikes]

# plot(Sp1.t/ms, Vals/ms)
# plt.axhline(y = 0, color = 'r', linestyle = '--')
# plt.axhline(y = 11/2, color = 'g', linestyle = '--')
# xlabel('Time')
# ylabel('Interspike Intervals')
# ylim(-50, 500);


# ## Complementary edges are excitatory connections

# #### 1. Pairs of neurons which don't have inhibitory connections are made excitatory

# In[24]:


# Init all 100 neurons with random starting V
# Make a matrix representing random connections, the groups G1 and G2 are connected with probability of 0.6
# (Every neuron from group G1 has 0.6 probability of connecting to each neuron of group G2)
# 
# C = [pre, post]


start_scope()

inh = 0.001
exc = 0.001
n = 100

seed(11)

G1 = NeuronGroup(n, eqs2, threshold='v>1', reset='v = 0', method='euler')
G2 = NeuronGroup(n, eqs2, threshold='v>1', reset='v = 0', method='euler')

group1init = randint(0, 9, size=n)/10
group2init = randint(0, 9, size=n)/10

G1.v = group1init
G2.v = group2init


# In[25]:


p = 0.8

I12 = choice(2, (n, n), p=[1-p, p])
sI1, tI1 = I12.nonzero()
I21 = choice(2, (n, n), p=[1-p, p])
sI2, tI2 = I21.nonzero()

E12 = 1 - I12
sE1, tE1 = E12.nonzero()
E21 = 1 - I21
sE2, tE2 = E21.nonzero()


imshow(I12);


# In[26]:


# Inhibitory synapses
I1 = Synapses(G1, G2, on_pre='v = out * (1-exp(-(tp * log(1/(1 - g*(v_post - inh)/curr)))/tp))')
I2 = Synapses(G2, G1, on_pre='v = out * (1-exp(-(tp * log(1/(1 - g*(v_post - inh)/curr)))/tp))')

# Excitatory synapses
E1 = Synapses(G1, G2, on_pre='v = out * (1-exp(-(tp * log(1/(1 - g*(v_post + exc)/curr)))/tp))')
E2 = Synapses(G2, G1, on_pre='v = out * (1-exp(-(tp * log(1/(1 - g*(v_post + exc)/curr)))/tp))')

# S1 = Synapses(G1, G2, on_pre='v_post = v_post - exc')
# S2 = Synapses(G2, G1, on_pre='v_post = v_post - exc')

I1.connect(i=sI1, j=tI1)
I2.connect(i=sI2, j=tI2)

E1.connect(i=sE1, j=tE1)
E2.connect(i=sE2, j=tE2)


M1 = StateMonitor(G1, 'v', record=True)
M2 = StateMonitor(G2, 'v', record=True)
Sp1 = SpikeMonitor(G1)
Sp2 = SpikeMonitor(G2)

run(500*ms)


# In[27]:


plot(Sp1.t/ms, Sp1.i, '.b', markersize=5)
plot(Sp2.t/ms, Sp2.i, '.r', markersize=5)
xlabel('Time (ms)')
ylabel('Neuron index')
text(150, 110, 'e = {}   |   n = {}   |   p = {}'.format(exc, n, p))
show()

# !neurons of the same indices across the two groups are not necessarily connected to each other!


# #### 2. Above, plus excitatory connections between neurons of the same group as well.

# In[4]:


# Init all 100 neurons with random starting V
# Make a matrix representing random connections, the groups G1 and G2 are connected with probability of 0.6
# (Every neuron from group G1 has 0.6 probability of connecting to each neuron of group G2)
# 
# C = [pre, post]


start_scope()

inh = 0.001
exc = 0.001
n = 100

seed(11)

G1 = NeuronGroup(n, eqs2, threshold='v>1', reset='v = 0', method='exact')
G2 = NeuronGroup(n, eqs2, threshold='v>1', reset='v = 0', method='exact')

group1init = randint(0, 9, size=n)/10
group2init = randint(0, 9, size=n)/10

G1.v = group1init
G2.v = group2init


# In[5]:


p = 0.8

I12 = choice(2, (n, n), p=[1-p, p])
sI1, tI1 = I12.nonzero()
I21 = choice(2, (n, n), p=[1-p, p])
sI2, tI2 = I21.nonzero()

E12 = 1 - I12
sE1, tE1 = E12.nonzero()
E21 = 1 - I21
sE2, tE2 = E21.nonzero()

imshow(I12);


# In[6]:


# Inhibitory synapses
I1 = Synapses(G1, G2, on_pre='v = out * (1-exp(-(tp * log(1/(1 - g*(v_post - inh)/curr)))/tp))')
I2 = Synapses(G2, G1, on_pre='v = out * (1-exp(-(tp * log(1/(1 - g*(v_post - inh)/curr)))/tp))')

# Excitatory synapses
E1 = Synapses(G1, G2, on_pre='v = out * (1-exp(-(tp * log(1/(1 - g*(v_post + exc)/curr)))/tp))')
E2 = Synapses(G2, G1, on_pre='v = out * (1-exp(-(tp * log(1/(1 - g*(v_post + exc)/curr)))/tp))')

# Self-excitatory synapses
S1 = Synapses(G1, G1, on_pre='v = out * (1-exp(-(tp * log(1/(1 - g*(v_post + exc)/curr)))/tp))')
S2 = Synapses(G2, G2, on_pre='v = out * (1-exp(-(tp * log(1/(1 - g*(v_post + exc)/curr)))/tp))')

# S1 = Synapses(G1, G2, on_pre='v_post = v_post - exc')
# S2 = Synapses(G2, G1, on_pre='v_post = v_post - exc')

I1.connect(i=sI1, j=tI1)
I2.connect(i=sI2, j=tI2)

E1.connect(i=sE1, j=tE1)
E2.connect(i=sE2, j=tE2)

S1.connect(condition='i != j')
S2.connect(condition='i != j')

M1 = StateMonitor(G1, 'v', record=True)
M2 = StateMonitor(G2, 'v', record=True)
Sp1 = SpikeMonitor(G1)
Sp2 = SpikeMonitor(G2)

run(500*ms)


# In[7]:


plot(Sp1.t/ms, Sp1.i, '.b', markersize=5)
plot(Sp2.t/ms, Sp2.i, '.r', markersize=5)
xlabel('Time (ms)')
ylabel('Neuron index')
text(150, 110, 'e = {}   |   n = {}   |   p = {}'.format(exc, n, p))
show()

# !neurons of the same indices across the two groups are not necessarily connected to each other!


# In[7]:


polared = lambda t, M : (M.v[:,t] * 2*pi, median(M.v[:,t] * 2*pi))
dummy = ones(n)*5
alph = 0.2

timer = 450

pld1 = polared(timer, M1)

pld2 = polared(timer, M2)

M1_avg = pld1[1]
M2_avg = pld2[1]

polar(pld1[0], dummy, 'bo', alpha = alph)
polar(pld2[0], dummy, 'ro', alpha = alph)

polar(M1_avg, 3, 'bo')
polar(M2_avg, 3, 'ro')


# In[34]:


times = range(0, 501, 10)
M1_avgs = []
M2_avgs = []

for t in times:
        pd1 = polared(t, M1)
        pd2 = polared(t, M2)
        M1_avgs.append(pd1[1])
        M2_avgs.append(pd2[1])

plot(times, M1_avgs, color='red')
plot(times, M2_avgs, color='blue')
plot(times, subtract(M2_avgs, M1_avgs), color='green')

