# Metropolis-Hastings: A Simple Example 

from __future__ import division
import numpy as np

import matplotlib.pyplot as plt
plt.style.use('bmh')

from scipy.stats import rv_discrete

# specify the posterior distribution 
pi1 = 0.2
thet1, thet2 = 0, 1

tau = rv_discrete(values=[(thet1, thet2), (pi1, 1-pi1)])

# specify the proposal distribution
def prop(q):
	"""
	Returns a dictionary with keys tau_1 and tau_2 and values discrete RV
	with probabilities determined by $q$ as in the matrix above.
	"""
	proposal = {}
	proposal[tau.xk[0]] = rv_discrete(values = [(tau.xk[0], tau.xk[1]), 
	                      (q, 1-q)])
	proposal[tau.xk[1]] = rv_discrete(values = [(tau.xk[0], tau.xk[1]), 
						  (1-q, q)])
	return proposal 
	
# calculate the probability of accepting proposal \vartheta given the current 
# \theta value

def alpha(vartheta, theta, q=0.2):
	"""
	Computes acceptance probability.
	"""
	post_ratio = tau.pmf(vartheta)/tau.pmf(theta)
	prop_ratio = (prop(q)[vartheta].pmf(theta)/prop(q)[theta].pmf(vartheta))
	return min(post_ratio/prop_ratio, 1)
	
from itertools import product

for vartheta, theta in product((thet1, thet2), (thet1, thet2)):
	alpha_str = (r'$\alpha(\vartheta=\tau_{}|\theta=\tau_{})={}$'
                 .format(vartheta+1,theta+1,alpha(vartheta,theta)))
	print alpha_str
                 
def transition_matrix(q=0.2):
	k11 = ( alpha(thet1,thet1,q=q)*prop(q)[thet1].pmf(thet1) 
           + (1-alpha(thet2,thet1,q=q))*prop(q)[thet1].pmf(thet2))

	k12 = alpha(thet2,thet1,q=q)*prop(q)[thet1].pmf(thet2)
	k21 = alpha(thet1,thet2,q=q)*prop(q)[thet2].pmf(thet1)

	k22 = ( alpha(thet2,thet2,q=q)*prop(q)[thet2].pmf(thet2) + 
           + (1-alpha(thet1,thet2,q=q))*prop(q)[thet2].pmf(thet1))

	return np.array([[k11,k12],[k21,k22]])
	
K = transition_matrix()
print 'K = ', K 

eig, eigv = np.linalg.eig(K)
print eig

qvec = np.array([0.000001, 0.2, 0.5, 0.99])

for q in qvec: 
	K = transition_matrix(q)
	eig, eigv = np.linalg.eig(K)
	plt.plot(eig[0]**np.arange(10), linewidth = 4)
	
_ = plt.legend([r'$q=0.00$', '$q=0.20$', '$q=0.50$', '$q=0.99$'], fontsize = 14, loc = 'center right')

plt.show()

























