import numpy as np
#Short program designed to determine longitude for a given epoch
#and then the final longitude for said given epoch after some time
pars = np.array(input('Specify epoch, period, time (all in days): '))
pars = pars.astype(np.float)
start = -360*pars[0]/pars[1]
final = start + ((pars[2] % pars[1])/pars[1])*360
print 'starting longitude: ', start
print 'ending longitude: ', final
