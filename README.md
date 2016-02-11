# spotz
Star spot model and analysis
The contents of /spotz represents the bulk of the coding done for my thesis work.
There are other programs that I wrote over the past 7 years, but those here represent the more important parts.
The qspotz program had a lengthy history before I began working on it and what I have here can be considered the next step in its evolution.
My thesis also includes a lot of science that previous astronomers hadn't worked with but that isn't too relevant here.
The main change I made to the program is the genetic algorithm along with a closer examination of the statistical significance of the model fit that wasn't always as prevalent in earlier work.
All of the python programs were written by me, originally in IDL, and are designed for different reasons.

I have also included some example data from KIC 6607150, there are two types.
  kplr006607150.425.435.dat contains 10 days worth of Kepler observations on this star, this is the data that is used in qspotz.
  kplr006607150.425.435.k##.dat are the 6 different qspotz model results, using 6 different values of k.
k is a differential rotation parameter that I was trying to solve for, it is disscused in greater detail in my thesis summary: http://nukeynick.github.io/ .

Each code has some comments to help clarify what they do.
I'll first talk about the little programs, lon.py, k.py, and lat.py.
These are small programs that I wrote up to help visualize the behavior of star spots given a few basic parameters.
k and lat are actually the same function, just written differently (one returns k given latitude the other returns latitidue given k).
lon.py was a little more useful and it makes an appearance in kchart.py and Cchart.py.  Given a star at some origin Epoch (another way of describing longitude) and given its rotational period and a time interval what will its final longitude be?
I actually did a lot more work with this when it was in IDL and some of the results appear in my final thesis.  However, the coding isn't too exciting and the science is only marginally useful.

The poorly named kchart.py (I have some history with programs named chart) uses lon.py to create a visual representation of how the spots move around a star given their different rotational periods.
This program was really just meant for me to help visualize the surface map of the stars.  I found the results to be useful though the math and science are not exactly advanced.
The main programming feature (for what it's worth) is that I plotted the spots with their actual sizes.

Cchart.py was the main tool I used in analyzing the model results.
A file like kplr006607150.425.435.k16.dat contains the top 10% of n runs of qspotz on a lightcurve segment of kplr006607150 covering days 425-435 for a differential rotation value (k) = 0.16.
Cchart.py takes all the k values for a segment (input notation: kplr006607150.425.435) and returns the five best chi2 values for each value of k.
I did a bunch of additional statistical analysis on these before settling on the current state of the program.
By comparing the results for one segment to the dozens of other modeled segements a picture emerges of what the true value of k is.
I like this program because it gives the user options on what type of output they want, including a text file which I did a bunch of work with but won't go into here.
There is also a bunch of basic string manipulation because I keep the data files located in different places (I took some of that stuff out for githib though).

lightcurve.py is a program that I would used to examine how effective my model reproduced the Kepler data.
It produces a plot with the Kepler data and model data on top of each other, as well as the residuals located just below.
This program also calls to generate.cpp to produce the model lightcurve since the model results are only in the form of parameter values (lightcurve.py also passes these parameters to generate.cpp).
As I'll descibe later the model takes a lightcurve and produces a set of parameters which can reproduce it.
generate.cpp is simply a stripped down version of the model working backwards, it takes parameters and produces a lightcurve.
This process is also a lot quicker since the model isn't trying to find a best solution to the model.

results.py may be my favorite program due to the various options available to the viewer.
Very little science is being done in this program as its main purpose is to produce the plots for my dissertation presentation.
However those plots contain a great deal of accumulated scientific knowledge.
The first part of the program is a little bulky, and based on what I know of Python it would be easy to clean up.
This is how it needed to be done in IDL though, highlighting one that languages shortcomings.
There are also a bunch of notes on the data I collected that have nothing to do with the functionality of the code.
The main code is essentially one big loop which allows the user to examine the data in different ways while maintaining a certain degree of usability.
Since the plots are all slightly different from each other there are also some conditional statements in there to keep each plot looking correct.

The real workhorse of this project was qspotz.cpp which is what took the vast amounts of data from Kepler, ran billions of individual models on it all and produced the most accurate spot model solutions possible (there may be non-spot related phenomena which also affect the lightcurve).
Everytime it is run there are a few parameters that need to be updated, the file name (and output name), radius, peak flux, inclination, equatorial period, and k.
These are the parameters which I had to solve for before modeling could begin (and that was a lot of work too, just not coding specific).
As I sit here typing this it also occurs to me that I could have automated the process to update these parameters in the code, so this part isn't as efficient as it could have been.
I then had the model produce a set of initial parameters for the lightcurve and the genetic algorithm takes over.
A genetic algorithm has NPOP number of members in a population, each member has all the parameters needed to produce a lightcurve (9 depending on certain factors).
The model, plus a bunch of additional checks, performs a linear least-squares regression between the model and the raw data and returns a chi2 for each memeber of NPOP.
Lower chi2 means the fit is good, the model then saves the top 10% of NPOP and creates a new generation of models.
The model then runs through NGEN number of generations and outputs the top 10% into a stat file (kplr006607150.425.435.k08.dat is one of these).
I typically ran the model with NPOP = 250 and NGEN = 200 to ensure the model reduced to good solutions everytime it ran.
qspotz was then run 50 times for each model and each value of k.
Due to some degeneracys in the physical nature of the stars it was not possible to directly solve for k with the accuracy I wanted, so I ran the models with different values of k and used Cchart.py to find which was best.
This model was then tested extensively to ensure it could perfom correctly, those tests are scientific in nature so I won't go into them here.
