#For a given segment, determine the five best models for a k value
#output is either a .mean file or a kfit.pdf

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
import os, numpy as np

#read in segment filename
segment = raw_input('Specify file name (format: kplr##.##.##): ')
#split file name into constituent parts (kplr###, range1, range2)
dir = str.split(segment,'.')
#create list of all corresponding segments
os.system('ls ' + dir[0] + '/' + segment + '* > temp.lis')
names=[]
with open('temp.lis','r') as f:
    for line in f:
        cols = line.strip().split()
        names.append(cols[0])

#select options for program
opt = raw_input('Save? ((f)ile/(p)lot/(n)one): ')
#set up output file for use
if opt == 'f' or opt == 'file':
    out_file = dir[0] + '.mean.dat'
    span = dir[1] + ':' + dir[2]

#arrays and other important things
limit = len(names)
ktmp=[0]*limit
chib=[0]*limit
ks=[0]*5*limit
chis=[0]*5*limit
chim=np.ndarray((2,limit)) #mean and std
col_names=('n','i','epoch1','per1','lat1','rad1','epoch2','per2','lat2','rad2','epoch3','per3','lat3','rad3','chi2')
#work on each model result
for j in range(0,limit):
#first determine k from filename
    tmp=str.split(names[j],'.')
    tmp=tmp[3]
    ktmp[j]= float(tmp[-2:])*.01
#read in data using pandas
    info = pd.read_table(names[j], names=col_names, header=None, sep=' ')
#sort chi2 and create plotting arrays for chi2 and k
    s_val = np.sort(info.chi2)
    for jj in range(0,5):
        chis[j*5+jj] = s_val[jj]
        ks[j*5+jj] = ktmp[j]

#find best, mean, and std
    chib[j]=s_val[0]
    chim[0][j]=np.mean(s_val[0:5])
    chim[1][j]=np.std(s_val[0:5])

    print 'For k = ',ktmp[j]
    print 'Mean chi2 = ',chim[0,j],' +- ',chim[1,j]
    print 'Best chi2 = ',chib[j]
    print ''

#appending .mean file
    if opt == 'f' or opt =='file':
        with open(out_file, "a") as ofile:
            ofile.write('%s %.2f %.5f %.5f %.5f \n' % (span,ktmp[j],chim[0,j],chim[1,j],chib[j]))

#plot final results
plt.plot(ks,chis,'*')
plt.xlabel('k values')
plt.ylabel(r'$\chi^2$')
plt.xlim([min(ks)-.01,max(ks)+.01])
plt.errorbar(ktmp,chim[0,:],yerr=chim[1,:],fmt='x',markersize=20)

if opt=='n' or opt=='none':
    plt.show()

if opt=='p' or opt=='plot':
    fname = segment + '.kfit.pdf'
    plt.savefig(fname)
    print 'Plot saved as',fname


#kplr006607150.425.435
