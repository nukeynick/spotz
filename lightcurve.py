#program plots data file against model file and displays residuals
import matplotlib.pyplot as plt
import pandas as pd
import os,numpy as np
from scipy.interpolate import interp1d

#read in segment filename
segment = raw_input('Specify file name (format: kplr##.##.##.k##.dat): ')
#split file name into constituent parts (kplr###, range1, range2, k##)
name = str.split(segment,'.')

#Create title, target name and time range
tit_var = name[0] + ' spot positions over days ' + name[1] + '-' + name[2]

#determine length of time for rotation
time = float(int(name[2])-int(name[1]))

#column headers (needed by pandas)
col_names=('n','i','epoch1','per1','lat1','rad1','epoch2','per2','lat2','rad2','epoch3','per3','lat3','rad3','chi2')

#locate file
#stat_file = name[0] + '/' + segment
#I typically store data in a subdirectory
stat_file = segment

#read in file
stats = pd.read_table(stat_file, names=col_names, header=None, sep=' ')

#locate index of best chi2 result (only records first instance)
best = stats['chi2'].idxmin()

#determine k from filename
tmp=name[3]
k=float(tmp[-2:])*.01

#determine pre-defined stellar parameters
def params(name):
    if name == 'kplr004831454_45':
        name = 'kplr004831454'
    star = np.array(['kplr003128488','kplr003743511','kplr004831454','kplr005524879','kplr006188162','kplr006607150','kplr011304230','kplr011560431'])
    s = np.where(star == name)
    R = np.array([.688,.99,.92,1.36,.83,1.05,.96,.85])
    P = np.array([5.99,5.57,4.97,11.17,7.78,3.80,7.01,3.02])
    F = np.array([1.06107,1.02089,1.00877,1.00318,1.0208,1.01118,1.016,1.02209])
    return R[s],P[s],F[s]

R, P, F = params(name[0])

#gives the best fit model parameters
print ''
print 'Minimum chi2   ',stats.chi2[best]
print 'best epochs    ',stats.epoch1[best],stats.epoch2[best],stats.epoch3[best]
print 'best periods   ',stats.per1[best],stats.per2[best],stats.per3[best]
print 'best latitudes ',stats.lat1[best],stats.lat2[best],stats.lat3[best]
print 'best radii     ',stats.rad1[best],stats.rad2[best],stats.rad3[best]
#needed to convert from epoch to longitude
lon = np.array([stats.epoch1[best],stats.epoch2[best],stats.epoch3[best]])
pers = np.array([stats.per1[best],stats.per2[best],stats.per3[best]])
#keeps longitude greater than -360 degrees
for j in range(len(lon)):
    while lon[j] > P:
        lon[j] -= pers[j]

sl = -360*(lon/pers)
print 'starting longitudes:',sl
print ''

#open original data
star_loc = name[0] + '.' + name[1] + '.' + name[2] + '.dat'
data_names=('x','y','err')
data = pd.read_table(star_loc, names=data_names, header=None, sep=' ')

print 'Pertinent model information for',name[0]
print 'Inclination:',stats.i[0]
print 'Equatorial Period:',P
print 'Stellar Radius:',R
print 'Max Flux:',F
print ''

#create file needed by generate.cpp to produce a model
temp_out = [F,stats.i[0],R,k,P,stats.epoch1[best],stats.per1[best],stats.lat1[best],stats.rad1[best],stats.epoch2[best],stats.per2[best],stats.lat2[best],stats.rad2[best],stats.epoch3[best],stats.per3[best],stats.lat3[best],stats.rad3[best]]
with open('temp.dat', "w") as tfile:
    for j in temp_out:
        tfile.write('%.5f\n' % j)

#create model and read in the model data
os.system('./gen')
model = pd.read_table('model.data.dat', names=data_names, header=None, sep=' ')

#interpolate so model 'x' points match data 'x' points
intfunc = interp1d(model.x,model.y,kind='linear')
yi = intfunc(data.x)

#plot model on top of data with residuals
plt.plot(data.x,data.y,'.',label='Kepler Data')
plt.plot(model.x,model.y,'.',label='Model Data')
plt.plot(data.x,data.y-yi+.98,'-',label='Residuals')
plt.plot([0,max(data.x)],[.98,.98],'-')
plt.ylim(min(data.y)-.05,max(data.y)+.02)
plt.xlim(min(data.x),max(data.x))
plt.xlabel('Time (days)')
plt.ylabel('Relative Flux')
legend = plt.legend(loc = 'lower center')

#option to save plot as a pdf
save = raw_input('Create a .pdf? (y/n) ')
if save=='y' or save=='yes':
    fname = name[0] + '.' + name[1] + '.' + name[2] + '.pdf'
    plt.savefig(fname)
    print 'Plot saved as',fname

if save=='n' or save=='no':
    plt.show()
