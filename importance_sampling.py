# let's say this is in test_branch

from __future__ import division

import numpy as np
from scipy.stats import norm, uniform 

import pandas as p
import matplotlib.pyplot as plt

plt.style.use('bmh')
np.random.seed(1234)


def f(x, phibar = -0.5, Vbar = 0.25):
	"""
	Posterior distribution of the set identified parameter theta. 
	See Equation (3.40) in the text. 
	"""
	phi1 = norm.cdf( (x-phibar)/np.sqrt(Vbar) )
	phi2 = norm.cdf( (x - 1.0 - phibar)/np.sqrt(Vbar) )
	return phi1 - phi2 
	
# Vbar controls the curvature around the mode of the posterior. 
x = np.linspace(2, -2, 100)
f1 = f(x, Vbar = 1/20)
f2 = f(x)
f3 = f(x, Vbar = 1/100) 

fig, ax = plt.subplots() 
ax.plot(x, f3, linewidth = 3)
ax.plot(x, f1, linestyle = 'dashed', linewidth = 3)
ax.plot(x, f2, linestyle = '-.', linewidth = 3)
ax.set_ylim(0, 1.1)
_ = ax.legend(map(r'$\bar V_{{\phi}} = {}$'.format,['1/100', '1/20', '1/4']))
# visualize the posterior: 
# plt.show()

# now we implement the importance sampling 
gc = norm(loc = 0.5, scale = np.sqrt(0.125))
gd = norm(loc = 0.5, scale = np.sqrt(0.500))
gc1 = gc.pdf(x)
gd1 = gd.pdf(x)

fig, ax = plt.subplots()
ax.plot(x, f3, linewidth = 3)
ax.plot(x, gc1, linewidth = 3, linestyle = 'dashed')
ax.plot(x, gd1, linewidth = 3, linestyle = '-.')
ax.set_ylim(0, 1.2)
_ = ax.legend([r'$\pi(\theta)$', r'$g_{c}(\theta)$', r'$g_{d}(\theta)$'] )
# plt.show()

# we'll compare the performance of importance sampling estimators to the (typical infeasible) performance of iid sampling. We assess the estimation for two choices of h.

def direct_sample(nsim, phibar = -0.5, Vbar = 0.25):
	""" 
	Directly sample from the set-identified posterior.
	"""
	phisim = phibar + norm.rvs(size = (nsim, ))*np.sqrt(Vbar)
	return phisim + uniform.rvs(size = (nsim,))
	
thsim = direct_sample(1e7, Vbar=1/100)
Eth_direct = np.mean(thsim)
Eth2_direct = np.mean(thsim**2)
Vth_direct = np.std(thsim)**2
Vth2_direct = np.std(thsim**2)**2

print "h(theta) = theta:\tE[h] = %s, V[h] = %s" % (Eth_direct, Vth_direct)
print "h(theta) = theta^2:\tE[h] = %s, V[h] = %s" % (Eth2_direct, Vth2_direct)

from itertools import product

mc = np.array([1e2, 5e2, 1e3, 5e3, 75e2, 1e4])
def monte_carlo_experiment(g = gc, mc = mc, nrun = 1000):
	"""
	Runs a monte carlo experiment
	For every nsim in mc, do nrep replications of importance, save important 
	estimates.	 	       
	"""
	resdf = []
	for nsim, i in product(mc, range(nrun)):
		# draw theta ~ g
		thsim = g.rvs(size = (nsim, ))
		# compute wt = f/g
		wtsim = f(thsim, Vbar = 1/100)/g.pdf(thsim)
		# normalize
		wtsim = wtsim/np.mean(wtsim)
		# compute
		Eth = np.mean(thsim * wtsim)/np.mean(wtsim)
		Eth2 = np.mean(thsim**2 * wtsim)/np.mean(wtsim)
		
		Oth = np.std(wtsim * (thsim-Eth_direct))**2
		Oth2 = np.std(wtsim * (thsim**2 - Eth2_direct))**2
		
		Ineff = 1 + np.std(wtsim)**2
		resdict = {'nsim': nsim, 'rep': i, 'Eth': Eth, 'Eth2': Eth2, 'Oth': Oth, 
					'Oth2': Oth2, 'Ineff': Ineff}
		resdf.append(resdict)
		
	IS = p.DataFrame(resdf)
	IS_mean = IS.groupby('nsim').mean()
	IS_var = IS.groupby('nsim').var()
	
	return IS_mean, IS_var	
	
	
IS_mean, IS_var = monte_carlo_experiment()

print IS_mean

print IS_var
		
IS_var['InEff_h'] = mc * IS_var['Eth'] / Vth_direct
IS_var['InEff_h2'] = mc * IS_var['Eth2'] / Vth2_direct
	
print IS_var
fig, ax = plt.subplots()
(mc * IS_var['Eth'] / Vth_direct).plot(ax=ax, marker='^', markersize=10)
(mc * IS_var['Eth2'] / Vth2_direct).plot(ax=ax, marker='o', markersize=10)
(IS_mean['Oth']/Vth_direct).plot(ax=ax, linestyle='dashed', marker='^', markersize=10)
(IS_mean['Oth2']/Vth2_direct).plot(ax=ax, linestyle='dashed', marker='o', markersize=10)
_ = IS_mean['Ineff'].plot(ax=ax)


IS_mean, IS_var = monte_carlo_experiment(g=gd)
fig, ax = plt.subplots()
(mc * IS_var['Eth'] / Vth_direct).plot(ax=ax, marker='^', markersize=10)
(mc * IS_var['Eth2'] / Vth2_direct).plot(ax=ax, marker='o', markersize=10)
(IS_mean['Oth']/Vth_direct).plot(ax=ax, linestyle='dashed', marker='^', markersize=10)
(IS_mean['Oth2']/Vth2_direct).plot(ax=ax, linestyle='dashed', marker='o', markersize=10)
_ = IS_mean['Ineff'].plot(ax=ax)
		
# plt.show()		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		






















































