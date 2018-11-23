from astropy.table import Table
import matplotlib.pyplot as plt
import sys
import numpy as np

from scipy.optimize import curve_fit

def cum_dist(x, xmin, xmax, p):
    return (x**(1-p)-xmin**(1-p))/(xmax**(1-p)-xmin**(1-p))


elem='a'
names=np.genfromtxt('names', dtype=str)

for ii in range(0, 4999):
	elem_list=np.concatenate([Table.read(name, '/{0}'.format(ii))[elem] for name in names])

	xx=np.sort(elem_list)
	yy=(np.array(range(len(xx)))).astype(float)/float(len(xx))
	popt,pcov=curve_fit(cum_dist, xx, yy, [1., 2., 2.])
	print popt

