#Program calculates and displays initial spot longitudes for a segment
#Then rotates and displays final spot longitudes for a segment
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd

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
#I typically store data in subdirectories
#stats = name[0] + '/' + segment
stats = segment

#read in file (PANDAS IS EASY!)
info = pd.read_table(stats, names=col_names, header=None, sep=' ')
#info['xxx'] gives the column labeled 'xxx'

#locate index of best chi2 result (only records first instance)
best = info['chi2'].idxmin()

#determine k from filename
tmp=name[3]
ktmp=float(tmp[-2:])*.01

#best value arrays
EA = [info.epoch1[best],info.epoch2[best],info.epoch3[best]]
PA = [info.per1[best],info.per2[best],info.per3[best]]
LA = [info.lat1[best],info.lat2[best],info.lat3[best]]
RA = [info.rad1[best],info.rad2[best],info.rad3[best]]

print ''
print 'Best values:'
print 'chi2 =',info.chi2[best]
print ''

for j in range(0,3):
    while EA[j] > PA[j]:
        EA[j] = EA[j]-PA[j]
    print 'Parameters for spot:', j+1
    print 'epoch =',EA[j]
    EA[j] = -(EA[j]/PA[j])*360
    print 'Longitude =',EA[j]
    print 'period =',PA[j]
    print 'latitude =',LA[j]
    print 'radius =',RA[j]
    print ''
    EA[j] += 360
    #using pyplot to create correctly sized circular spots
    spts = plt.Circle((EA[j],LA[j]),radius=RA[j],fc='g')
    plt.gca().add_patch(spts)

#rotate = raw_input('Hit enter to rotate star...')

#Determine spot longitudes after time has elapsed
#(time range is specified in input file name)
for j in range(0,3):
    EA[j] = EA[j] + ((time % PA[j])/PA[j])*360
    while EA[j] > 360: EA[j] -= 360
    print 'Final longitude for spot',j+1,'=',EA[j]
    fpts = plt.Circle((EA[j],LA[j]),radius=RA[j],fc='r')
    plt.gca().add_patch(fpts)

#plot parameters
plt.axis('scaled')
plt.xlim([0,360])
plt.ylim([0,90])
plt.xlabel('Longitude')
plt.ylabel('Latitude')
#legend parameters
plt.legend([spts,fpts],['Start Positions','Final Positions'],bbox_to_anchor=(0., 1.02, 1., 1.02), loc=3, ncol=2, mode="expand", borderaxespad=0.)
plt.title(tit_var,y=1.2)
plt.show()

#kplr006607150.425.435.k08.dat
