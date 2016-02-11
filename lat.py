import math, numpy as np
#Small program to determine spot latitude for a given set of parameters
pars = np.array(input('Specify equatorial period, spot period, and k: '))
pars = pars.astype(np.float)

Beta = (180/math.pi)*math.asin(math.sqrt((1-pars[0]/pars[1])/pars[2]))

print Beta
