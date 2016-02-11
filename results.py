#Results.py is the plot making program for the final results
#Is normally run after the analysis routines have been completed

import pylab, math
import numpy as np

star=['kplr003128488','kplr003743511','kplr004831454','kplr005524879','kplr006188162','kplr006607150','kplr011304230','kplr011560431']
k = [0.08,0.09,0.14,0.14,0.09,0.085,0.09,0.09]
err = [0.005,0.01,0.03,0.03,0.01,0.01,0.01,0.01]
peq = [5.99,5.57,4.97,11.17,7.78,3.80,7.01,3.02]
T = [4713,5746,5615,6200,4968,5775,5608,5094]
vsini = [5.6,7.6,6.6,5.6,5.06,12.7,6.74,13.36]
dOmega = [0.08,0.10,0.18,0.08,0.07,0.14,0.08,0.19]
#Oerr = ((k+err)*2*math.pi/peq)-dOmega
Oerr = [0.00916039,0.01280405,0.03491781,0.01562592,0.01076074,0.01707963,0.00963174,0.01805249]


Sci_star = ['Sun','EpsEri','KIC5110407','KIC7985370','KIC7765135','KIC8429280','CoRoT-2','CoRoT-6','KCeti']#,'CoRoT-7','LOPeg','ABDor','RXJ','PZTel','LQHya','V374Peq','HD197890']
S_k = [0.19,0.11,0.12,0.08,0.07,0.05,0.11,0.12,0.09]#,0.06,0.0494,0.005,0.006,0.015,0.0004,0.002,0.002]
S_err = [0.0,0.02,0.04,0.01,0.01,0.01,0.02,0.02,0.01]#,0.00,0.01,0.001,0.001,0.001,0.00003,0.0005,.0005]
S_peq = [25.05,11.2,3.47,2.856,2.41,1.16,4.5,6.35,8.77]#,23.6,1.6,0.515,0.31,0.95,0.445654,0.38,.423]
S_T = [5778,5084,5200,5815,5835,5238,5625,6090,5560]#,5275,5000,5000,5750,5070,3200,4750]
S_vsini = [2.0,2.4,31.0,18.0,21.6,37.0,10.1,7.5,4.5]#,3.5,28.0,91.0,115.0,68.0,36.5,132.0,69] #31 may be wrong
S_dOmega = [0.048,0.06,0.22,0.18,0.18,0.194,0.153,0.119,0.06]#,0.016,0.27,0.0618,0.12,0.101,0.0063,0.03,0.035]
S_Oerr = [0.00,0.01,0.07,0.005,0.005,0.03,0.02,0.02,0.01]#,0.00]

#sources
#EpsEri: Croll 2006b
#KIC5110407: Roettenbacher 2013, bad values?
#KIC7985370: Frohlich 2012
#KIC7765135: Frohlich 2012
#KIC8429280: Frasca 2011
#CoRoT-2: Frohlich 2009, Silva-Valio 2011 vsini=10.1 assuming i=87
#CoRoT-6: Lanza 2011
#KCeti: Walker 2007
#unused sources
#CoRoT-7: Lanza 2010, k is a lower limit
#LOPeg: Barnes 2005b
#ABDor: Cameron 2001, very young, polar spots, using average dOmega
#RXJ: Donati 2000, not a MS star, polar spots
#PZTel: Barnes 2000, does not contain the values here
#LQHya: Kovari 2005, alternate dOmega = 0.022
#V374Peg: Morin 2008, M dwarf, strong magnetic fields
#HD197890: Barnes2005a, strong magnetic fields (Speedy Mic)
#Hill 2014 AE Aqr: tidally locked star
#Hussain 2006 V471 Tau: tidally locked star, no dOmega
#Petit 2004 HR1099: subgiant

axes=[k,dOmega,peq,T,vsini,err,Oerr]
Saxes=[S_k,S_dOmega,S_peq,S_T,S_vsini,S_err,S_Oerr]
titles = ['k','Differential Angular Velocity (dOmega) [rad/day]','Period [days]','Temperature [K]','vsini [km/s]']

j=0
opt='c'

while opt=='c' or opt=='s' or opt=='continue' or opt=='save':
#Selection Guide, appears every 5 iterations
    if j % 5 == 0:
        print ''
        print('Axis selection guide (x,y)')
        print('k (1), dOmega (2), Period (3), Temperature (4), vsini (5)')
        print('')
#Select type of plot to make
    xin = int(raw_input('Select x-axis (3, 4, or 5: '))
    yin = int(raw_input('Select y-axis (1 or 2: '))
    while yin != 1 and yin != 2:
        print 'ERROR: INVALID Y-AXIS PLEASE SELECT EITHER k (1) OR dOmega (2)!'
        yin = int(raw_input('Select y-axis: '))

    xin -= 1
    yin -= 1

    pylab.plot(axes[xin],axes[yin],'.',label='Science Targets')
    pylab.plot(Saxes[xin],Saxes[yin],'*',label='Literature Targets')
    pylab.xlabel(titles[xin])
    pylab.ylabel(titles[yin])
    if xin != 2 or yin != 1:
        legend = pylab.legend(loc = 'lower right')
    elif xin == 2 or yin == 1:
        legend = pylab.legend(loc = 'center right')

    if opt=='c' or opt=='continue':
        pylab.show()
    if opt=='s' or opt=='save':
        fname = raw_input('Specify file name: ')
        pylab.savefig(fname)
        print 'Plot saved as',fname

    if j % 5 == 0:
        print 'Options: (c)ontinue, (e)nd, (s)ave plot'
        print ''
    opt = raw_input('Option Choice? (c,e,s): ')
    if opt=='e' or opt=='end':
        print 'Program Complete'

    j += 1
