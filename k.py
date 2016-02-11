import math, numpy as np
#Small program to determine differential rotation for a given set of parameters
pars = np.array(input('Specify equatorial period, spot period, and spot latitude: '))
pars = pars.astype(np.float)

lat = pars[2]*math.pi/180
k = (1/(math.sin(lat))**2)*(1-pars[0]/pars[1])
print k
