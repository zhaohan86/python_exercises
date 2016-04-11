"""
Implement MCMC Metropolis-Hasting algorithm for a simple gaussian mixture. Example taken from:
http://www.nehalemlabs.net/prototype/blog/2014/02/24/an-introduction-to-the-metropolis-method-with-python/
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab 

def q(x, y):
	g1 = mlab.bivariate_normal(x, y, 1.0, 1.0, -1, -1, -0.8)
	g2 = mlab.bivariate_normal(x, y, 1.5, 0.8, 1, 2, 0.6)
	return  0.6*g1 + 28.4*g2  / (0.6 + 28.4) 
	
''' Metropolis Hastings'''

N = 100000
s = 10
r = np.zeros(2)
p = q(r[0], r[1])
# p = 0.00528913714392,  r = [ 0.  0.]

samples = []
for i in xrange(N):
	rn = r + np.random.normal(size = 2)
	pn = q(rn[0], rn[1])
	if pn >= p: # accept the draw with probability 1, move the chain 
		p = pn
		r = rn
	else: # accept the raw with probability pn/p
		u = np.random.rand()
		if u < pn/p:
			p = pn
			r = rn
	if i % s == 0: # pick every 10th draw in the markov chain
		samples.append(r)

samples = np.array(samples)
# print samples.shape
# plot draws
plt.scatter(samples[:, 0], samples[:, 1], alpha=0.5, s=1)
# plot target
dx = 0.01
x = np.arange(np.min(samples), np.max(samples), dx)
y = np.arange(np.min(samples), np.max(samples), dx)
X, Y = np.meshgrid(x, y)
Z = q(X,Y)
CS = plt.contour(X, Y, Z, 10)
plt.clabel(CS, inline = 1, frontsize = 10)
plt.show()