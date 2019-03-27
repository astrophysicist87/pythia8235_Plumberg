import numpy as np
import matplotlib as mpl, matplotlib.pyplot as plt
import sys

mpl.rcParams['pdf.fonttype'] = 42

data = np.loadtxt(sys.argv[1])

# get rid of weird outlier sub-collisions
data = data[np.where((np.abs(data[:,0])<=15.0) & (np.abs(data[:,1])<=15.0))]

#plt.scatter(data[:,0], data[:,1], s=10.0, c='blue', alpha=0.1, edgecolors='none')
#plt.hist2d(data[:,0], data[:,1], bins=50)
nbins=40
H, xedges, yedges = np.histogram2d(data[:,0], data[:,1], bins=nbins)

xcenters = (xedges[:-1] + xedges[1:]) / 2
ycenters = (yedges[:-1] + yedges[1:]) / 2
X, Y = np.meshgrid(xcenters, ycenters)

#plt.contourf(X, Y, H)
cm = plt.cm.gnuplot
#plt.imshow(H, cmap=cm, extent=[xedges.min(), xedges.max(), yedges.min(), yedges.max()], origin='lower', interpolation='bicubic')
plt.imshow(H, cmap=cm, extent=[-15.0, 15.0, -15.0, 15.0], origin='lower', interpolation='bicubic')
#plt.pcolormesh(xi, yi, H.reshape(xi.shape))

plt.xlabel(r'$x$ (fm)')
plt.ylabel(r'$y$ (fm)')
#plt.xlim(-15.0, 15.0)
#plt.ylim(-15.0, 15.0)
plt.axes().set_aspect('equal')

plt.savefig("scatter.pdf")




'''
# libraries
#import matplotlib.pyplot as plt
#import numpy as np
from scipy.stats import kde
 
# create data
print 'Loading data...'
data = np.loadtxt(sys.argv[1])
print 'done.'
x = data[:,0]
y = data[:,1]
 
# Evaluate a gaussian kde on a regular grid of nbins x nbins over data extents
nbins=300
k = kde.gaussian_kde([x,y])
xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
zi = k(np.vstack([xi.flatten(), yi.flatten()]))
 
# Make the plot
plt.pcolormesh(xi, yi, zi.reshape(xi.shape))
plt.savefig('scatter.pdf')
'''
