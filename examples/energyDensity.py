import numpy as np
import matplotlib as mpl, matplotlib.pyplot as plt
import sys

mpl.rcParams['pdf.fonttype'] = 42

data = np.loadtxt(sys.argv[1])
(nRows, nColumns) = data.shape	#nColumns is 11, at the moment
data = data.reshape([nRows / 2, 2, nColumns])
indices = data[:,:,0:3]	#event ID, particle ID, and particle status
pVecs = data[:,:,3:7]
xVecs = data[:,:,7:11]


# for convenience
minkowskiNorm = lambda p: np.sqrt(p[0]**2-p[1]**2-p[2]**2-p[3]**2)


sqrt_sInvs=minkowskiNorm((pVecs[:,0]+pVecs[:,1]).T)
xCens = (0.5*(xVecs[:,0]+xVecs[:,1]))[:,[1,2]]	# only need average x and y positions for colliding pair of partons


#data = np.c_[ xCens, sqrt_sInvs ]
data = (np.c_[ xCens, sqrt_sInvs ])[np.where((np.abs(xCens[:,0])>0.5) | (np.abs(xCens[:,1])>0.5))]



#generate density plot as usual (note weights option for center-of-mass energies)
nbins=40
H, xedges, yedges = np.histogram2d(data[:,0], data[:,1], weights=data[:,2], bins=nbins)

xcenters = (xedges[:-1] + xedges[1:]) / 2
ycenters = (yedges[:-1] + yedges[1:]) / 2

X, Y = np.meshgrid(xcenters, ycenters)

cm = plt.cm.gnuplot
plt.imshow(H, cmap=cm, extent=[xedges.min(), xedges.max(), yedges.min(), yedges.max()], origin='lower', interpolation='bicubic')
plt.colorbar()

plt.xlabel(r'$x$ (fm)')
plt.ylabel(r'$y$ (fm)')

plt.axes().set_aspect('equal')

plt.savefig("energyDensity.pdf")

# End of file
