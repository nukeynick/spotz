#include <cmath>
#include <stdarg.h>	
#include <stdio.h>
#include <cstdlib>
#include <time.h>
#include "spotz.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <omp.h>

using namespace std;

//these values need to be global parameters
//Spot model constants
const double m_Timeinc = .02;  //model will have five times as many data points
const double m_TimeMin = 0.0;
const double m_TimeMax = 40;
const double m_Period = 0;
const double m_Limit = 1e5;
const double m_LimitValue = 1e5;
const int m_NumSpots = 3;
const int m_FitSpot = 3;
//int m_RepeatOverlap = 0; //may need to assign a starting value for these two
//int m_SpotsOverlap = 0;  //otherwise they are defined in the header file
const int m_Iter = 0;
const int m_IterFitPeriod = 0;
const double m_IterFitTime_Min = 0;
const double m_IterFitTime_Max = 9999;
const int m_StarWant = 0;
const int m_NeedToStop = 0;
const int m_NeedToStopTISMO = 0;
const int m_StripWant = 0;
const int m_gMagnitude = 0;
const int m_Magnitude = 0;//boolean
const int m_BuddingBool = 1;//boolean
//const int m_Degrees_of_Freedom = 1;  //define below
const int m_Number_Fitted_Params = 6;
const int m_Number_Of_Steps_Run_Fit = 0;
const double m_ONE_PI = 3.1415926535897932384626433832795;
//AGA constants
const int NPOP = 250;  //number of memebers in the AGA population
const int NGEN = 200;  //number of generations to explore
const int EFAC = 10;  //elitism factor (save best 10%)
const double randmax = RAND_MAX;

/*
paramslocal is the initialization file, most of this information will be selected by the Genetic Algorithm, some of the values can be set to constants.
paramslocal[#]
0 = USI unspotted magnitude of the star.  This will be turned off.
1 = U unspotted intensity of the star/system.  Set to 1.
2 = L1 unspotted intensity of the spotted star.  Set to 1.
3 = i inclination of the stellar rotation axis to the line of sight.  57.3
11 = numspots number of spots.  GA!
12 = RStar the radius of the star in solar radii.  Use exisiting models
13 = PHASEmin the eopch of minimum phase (HJD).  
14 = NUMFit = the number of fitted parameters (used in reduced chi2).
30 = k differential rotation parameter. Solve using the model
31 = Peq equatorital period (days).  Need this, get from Tom
#00 = kw# the flux ratio between spot # and the unspotted photosphere (0-1) GA!
#01 = u# linear limb darkening coefficient (0-1).  Contact Kepler people?
#02 = E# spot epoch in HJD.  GA!
#03 = p# spot period in days.  GA!
#04 = l# spot longitude (0-360), lambda.  GA! (epoch covers this)
#05 = b# spot latitude (-90-90), beta.  GA!
#06 = g# spot angular radius (0-90), gamma.  GA!
*/

/*
Changes to be made for GA:
kw# needs to be constrained, check literature to get some ideas.
u# will need to come from current limb darkening tables (depends on log(g) and T)
Current plan for numspots is to find the best fits for 1,2,3,... spots and
compare afterward.
Rstar should come from current models, use Kepler's temperatures for now, of
course these values are wrong at present and will need to be updated to reflect
changes. (6/20/11)
PHASEmin is somewhat arbitrary, good idea to start at a flux max or min?
If left undefined the program will just start with the lowest time epoch.
NUMFit should be the number of parameters being fitted by the GA.  Which are
inclination, numspots, and the 4 spot parameters for each spot.
At present define the P_eq using Period04, use gaussian values for spot models
7nd solve for the differential rotation k parameter when a fit is found, so
set params_local[30] = -10000
*/

/*
printinfo_int should be a boolean operator, 1=print 0=no printing
paramslocal is defined above
Ntop defines the max number of data points in model/data
jd_array is the time of each data point in the original light curve
flux_array is the data from the original light curve, counterpart of jd_array
merr_array is flux error, used in the chi square analysis
model_array is the model flux after the spots have been considered
Ntop_model_return is the new Ntop after we increase number of model points
jd_model_array same as above except for time
model_plot_array same as model_array except with more points
residual_array is the difference between model and data
flat_array is used to make a nice line in the plot of residuals = 0
Timemin and TimeMax are pulled out of the maxmin function
Timeinc is the increment of the time from the origial light curve
Degrees_of_Freedom is used in chi square analysis, definition uncertain
ChiCare_return is the result of chi square analysis*/

main() {    
  ofstream outfile ("kplr004831454_45.865.895.k18.dat",ofstream::binary);
  if(!outfile) {
    cout << "Cannot open file.";
    return 1;
  }
/* 
  ofstream outfile2 ("trend.dat",ofstream::binary);
  if(!outfile2) {
    cout << "Cannot open file.";
    return 1;
  }
*/
  srand ( time(NULL) ); 
  time1 = time(NULL);
  FILE *fp;  //input file
  //char infile[100];
  //  printf("Enter light-curve: ");
  //  scanf("%s", &infile);
  char infile[100] = "kplr004831454.865.895.dat";  //filename in the program
  char new_line[100];  //array for one line of the input file

  //stuff to go into runmodel_special
  printinfo_int = 0;

  i = 0;
  if ((fp = fopen(infile,"r")) == NULL)  {
    cout << "Could not open file for reading: " << infile;
    abort();
  }
  while (fgets(new_line, sizeof(new_line), fp))  {
    if (sscanf(new_line,"%lf %lf %lf",&jd,&mag,&err) == 3)
      {
	m_jd[i] = jd;
	m_flux[i] = mag;
	m_merr[i] = err;
	i++;
      }
  }
  m_Ntop = i;
  fclose(fp);
  m_Degrees_of_Freedom = m_Ntop-m_Number_Fitted_Params;

  //Moved from runmodel_special, only need to do once
  returnminmax(m_Ntop,m_jd,&jdmin,&jdmax);
  TimeMin = (jdmin);
  TimeMax = (jdmax)+1;

  //paramslocal_array will come from GA portion, eventual GA for kw and u
  for(m=0;m<NPOP;m++) {
    paramslocal_array[m][0] = 1.00877;
    paramslocal_array[m][1] = 1.00877;
    paramslocal_array[m][2] = 1.00877;
    paramslocal_array[m][11] = 3;//Numspots        ADJUSTABLE!!
    paramslocal_array[m][12] = .922;//Rstar_solar   UPDATE!!
    paramslocal_array[m][13] = TimeMin;//PhaseMin
    paramslocal_array[m][30] = 0.18;//k
    paramslocal_array[m][31] = 4.96623;//              UPDATE!!
    paramslocal_array[m][100] = .3;//kw1
    paramslocal_array[m][101] = .684;//u1
    paramslocal_array[m][104] = 0;
    paramslocal_array[m][200] = .3;
    paramslocal_array[m][201] = .684;
    paramslocal_array[m][204] = 0;
    paramslocal_array[m][300] = .3;
    paramslocal_array[m][301] = .684;
    paramslocal_array[m][304] = 0;
  }
 
  int numspots = int(paramslocal_array[0][11]);

  //AGA TIME!!!
  //    outfile.open("garbage");
  // nvar == number of parameters (inclination + 4 for every spot) maybe do small scale fitting for kw?
  nvar = 7;
  paramslocal_array[0][14] = nvar;//NumFit

  //Run for ~1000 times to examine success rate (Monte-Carlo)

  for(nn=0;nn<50;nn++)  {

  spread[0] = 1.;
  spread[1] = .8;
  spread[2] = .15;
  spread[3] = .8;
  spread[4] = .15;
  spread[5] = .8;
  spread[6] = 1.1;
  spread[7] = .15;

  // Populate initial generation based on each parameters limits
#pragma omp parallel for
  for(m=0;m<NPOP;m++) {

   fit_params[m][0] = 45.0;//inclination
   fit_params[m][1] = m_jd[0]+GRAND((paramslocal_array[m][31]),spread[1]);
   while(fit_params[m][1] < m_jd[0]) fit_params[m][1] = m_jd[0]+GRAND((paramslocal_array[m][31]),spread[1]);//Epoch
   fit_params[m][2] = ERAND(20,1.5);
   fit_params[m][3] = m_jd[0]+GRAND(paramslocal_array[m][31],spread[3]);
   fit_params[m][4] = ERAND(20,1.5);
   fit_params[m][5] = m_jd[0]+GRAND(paramslocal_array[m][31],spread[5]);
   fit_params[m][6] = ERAND(1.08*paramslocal_array[m][31],paramslocal_array[m][31]);//period
   
   fit_params[m][7] = ERAND(20,1.5);

  }

  // Begin looping over generations and members        
  for(n=0;n<NGEN;n++)  //n is generation identifier
    {
      for(m=0;m<NPOP;m++)  //m is population member
	{
	  paramslocal_array[m][3] = fit_params[m][0];
          paramslocal_array[m][102] = fit_params[m][1];
	  paramslocal_array[m][103] = GRAND(5.15204,.005);
	  paramslocal_array[m][105] = asin(sqrt((1-paramslocal_array[m][31]/paramslocal_array[m][103])/paramslocal_array[m][30]))*180./m_ONE_PI;
          paramslocal_array[m][106] = fit_params[m][2];
          paramslocal_array[m][202] = fit_params[m][3];
          paramslocal_array[m][203] = GRAND(5.26376,.005);
	  paramslocal_array[m][205] = asin(sqrt((1-paramslocal_array[m][31]/paramslocal_array[m][203])/paramslocal_array[m][30]))*180./m_ONE_PI;
          paramslocal_array[m][206] = fit_params[m][4];
          paramslocal_array[m][302] = fit_params[m][5];
          paramslocal_array[m][303] = fit_params[m][6];
          paramslocal_array[m][305] = asin(sqrt((1-paramslocal_array[m][31]/fit_params[m][6])/paramslocal_array[m][30]))*180./m_ONE_PI;
          paramslocal_array[m][306] = fit_params[m][7];

          runmodel_special(printinfo_int,paramslocal_array[m],m_Ntop,m_jd,m_flux,m_merr,m_model,&m_Ntop_model,m_jd_model,m_model_plot,m_resid,m_flat,m_Timeinc,m_Degrees_of_Freedom,numspots,TimeMin,TimeMax,&chicare);
	  stdevvals[m] = chicare;
	}

      // Breed solutions
      help[0] = 9.9E9;
      for(m=0;m<NPOP;m++)  // Find lowest chi2
	if(stdevvals[m] < help[0])
	  help[0] = stdevvals[m];  // help[0] is lowest chi2

      for(m=1;m<NPOP/EFAC;m++)
	{
	  help[m] = 9.9E9;
	  for(i=0;i<NPOP;i++)
	    if(stdevvals[i] < help[m] && stdevvals[i] > help[m-1])
	      help[m] = stdevvals[i];  // Makes help[m] the next lowest stdev after help[m-1]
	}

      for(m=0;m<NPOP/EFAC;m++)  // Elitism - Save top 10%
	for(i=0;i<NPOP;i++)
	  if(stdevvals[i] == help[m])
	    for(j=1;(j<nvar+1);j++)  {
              offspring[m][j] = fit_params[i][j];
          //recalculate respective latitudes
              paramslocal_array[m][305] = asin(sqrt((1-paramslocal_array[m][31]/offspring[m][6])/paramslocal_array[m][30]))*180./m_ONE_PI;
              if(paramslocal_array[m][305] > 75.0)
                help[m] += 10;

            }
      if(n == NGEN-1) {
	for(m=0;m<NPOP/EFAC;m++) {
	  //output file
	  outfile << nn << " " << fit_params[m][0] << " " << offspring[m][1] << " " << paramslocal_array[m][103] << " " << paramslocal_array[m][105] << " " << offspring[m][2] << " " << offspring[m][3] << " " << paramslocal_array[m][203] << " " << paramslocal_array[m][205] << " " << offspring[m][4] << " " << offspring[m][5] << " " << offspring[m][6] << " " << paramslocal_array[m][305] << " " << offspring[m][7] << " " << help[m] << endl;
	}
      }
//      outfile2 << n << " " << help[0] << endl;

      wsum = 0.0;  // Create weighting scheme
      for(m=0;m<NPOP;m++)
      wsum += 1.0/pow(stdevvals[m],2);  // Weight by 1 over the chi^2  higher value of exponent more selective over whole parameter space. 2 is pretty good

      wvals[0] = 0.0;
      for(m=0;m<NPOP;m++)
	wvals[m+1] = wvals[m] + 1.0/pow(stdevvals[m],2);

      for(kk=NPOP/EFAC;kk<NPOP;kk++)  // Now create rest of new generation
	{
	  randc = wsum*rand()/randmax;
	  for(m=0;m<NPOP;m++)
	    if(randc > wvals[m] && randc < wvals[m+1])
	      for(j=1;(j<nvar+1);j++)
		{
	          offspring[kk][j] = GRAND(fit_params[m][j],spread[j]);
		  while(offspring[kk][1] > (jdmin+2.0*paramslocal_array[kk][31])) offspring[kk][1] = GRAND(fit_params[kk][1],spread[1]);
		  while(offspring[kk][1] < jdmin) offspring[kk][1] = GRAND(fit_params[kk][1],spread[1]);
                  while(offspring[kk][2] < 1.5) offspring[kk][2] = GRAND(fit_params[kk][2],spread[2]);
                  while(offspring[kk][2] > 22) offspring[kk][2] = GRAND(fit_params[kk][2],spread[2]);
		  while(offspring[kk][3] > (jdmin+2.0*paramslocal_array[kk][31])) offspring[kk][3] = GRAND(fit_params[kk][3],spread[3]);
		  while(offspring[kk][3] < jdmin) offspring[kk][3] = GRAND(fit_params[kk][3],spread[3]);
                  while(offspring[kk][4] < 1.5) offspring[kk][4] = GRAND(fit_params[kk][4],spread[4]);
                  while(offspring[kk][4] > 22) offspring[kk][4] = GRAND(fit_params[kk][4],spread[4]);
                  while(offspring[kk][5] > (jdmin+2.0*paramslocal_array[kk][31])) offspring[kk][5] = GRAND(fit_params[kk][5],spread[5]);
                  while(offspring[kk][5] < jdmin) offspring[kk][5] = GRAND(fit_params[kk][5],spread[5]);
                  while(offspring[kk][6] < paramslocal_array[kk][31] || sqrt((1-paramslocal_array[kk][31]/offspring[kk][6])/paramslocal_array[kk][30]) > .99 || offspring[kk][6] > 1.2*paramslocal_array[kk][31]) offspring[kk][6] = GRAND(fit_params[kk][6],spread[6]);
                  while(offspring[kk][7] < 1.5) offspring[kk][7] = GRAND(fit_params[kk][7],spread[7]);
                  while(offspring[kk][7] > 20) offspring[kk][7] = GRAND(fit_params[kk][7],spread[7]);

		}
	}

      for(m=0;m<NPOP;m++) { // Now make current generaton elitist parents plus new offspring
	for(j=1;(j<nvar+1);j++) {
          fit_params[m][j] = offspring[m][j];
        }
      }
      //modulate generation range of individual parameters by chi2 of the individual parameters	

      for(j=1;(j<nvar+1);j++) { // Make spread equal to current stdev
        for(m=0;m<NPOP;m++) {
	  Jparam[m] = fit_params[m][j];
	}
	spread[j]=pow(0.1,1.0/help[0])*STDEV(Jparam,NPOP);
        if(j==2) spread[2]=pow(.1,1.0/help[0])*STDEV(Jparam,NPOP)*.1;
        if(j==4) spread[5]=pow(.1,1.0/help[0])*STDEV(Jparam,NPOP)*.1;
        if(j==7) spread[8]=pow(.1,1.0/help[0])*STDEV(Jparam,NPOP)*.1;

	//cout << spread[j] << endl;
	}
      //cout << "       " << help[0] << endl;;
    }  // Do next generation (switch to the correct fit_params)
  cout << nn << " " << help[0] << endl;
  }
  time2 = time(NULL);
  cout << time2-time1 << "secs" << endl;
  outfile.close();
//  outfile2.close();
  return 0;
}
///////////////////////////////////////////////////////////////////////////////
//The calling sequence and actual procedure is here
void runmodel_special(int printinfo_int, double paramslocal_array[MAXsmall], int Ntop, double jd_array[MAX], double flux_array[MAX], double merr_array[MAX], double model_array[MAX], int *Ntop_model_return, double jd_model_array[MAX], double model_plot_array[MAX], double resid_array[MAX], double flat_array[MAX], double Timeinc, int Degrees_of_Freedom, int numspots, double TimeMin, double TimeMax, double *ChiCare_return)
{
  CheckParams(paramslocal_array);
  int spots_overlap = 0;
  double degrees_overlap = 0, degrees_overlap_care = 0;
  //double U, Epoch, p, inc, l, b, g;
  
  if(paramslocal_array[13] > -9999){PhaseMin = paramslocal_array[13];}
  else{PhaseMin = TimeMin;}
  //we probably want to calculate more model points than data points in the case that the period is very fast compared to our orbital period
  for(int ii=0; ii<Ntop; ii++){
    double t = jd_array[ii];
    for(int bb=1; bb<=numspots; bb++){
      m_Kw_sigmac[bb] = sigmac_call(t,paramslocal_array[bb*100+1],paramslocal_array[bb*100+2],paramslocal_array[bb*100+3],paramslocal_array[3],paramslocal_array[bb*100+4],paramslocal_array[bb*100+5],paramslocal_array[bb*100+6]);
    }
    double lcq = lcq_submit(numspots,m_Kw_sigmac,paramslocal_array);
    model_array[ii] = lcq;
    if(EnsureSpotsDontOverlap(t,printinfo_int,numspots,paramslocal_array,&degrees_overlap)==1){
      spots_overlap = 1;
      degrees_overlap_care = degrees_overlap;
    }
    if(printinfo_int==1){
      //fprintf(Printz,"ii: %d, time: %lf, model: %lf",ii,t,model_array[ii]);
    }
  }
  //This segment produces the values needed for a pretty plot
  /*
  //NOW DO THE SAME THING WITH MORE POINTS
  //Now we determine how many extra points to plot
  double time_step = Timeinc/10.;
  int Ntop_model = int((TimeMax-TimeMin)/time_step);
  if(Ntop_model > MAXmed){
    Ntop_model = MAXmed;
    time_step = (TimeMax-TimeMin)/double(Ntop_model);
  }
  for(int jj=0; jj<Ntop_model; jj++){
    double t = TimeMin +(jj)*time_step;
    jd_model_array[jj] = t;
    
    //double sigmac1 = sigmac_b(t,U,Epoch1,p1,inc,l1,b1,g1);
    for(int bb=1; bb<=numspots; bb++){
      m_Kw_sigmac[bb] = sigmac_call(t,paramslocal_array[bb*100+1],paramslocal_array[bb*100+2],paramslocal_array[bb*100+3],paramslocal_array[3],paramslocal_array[bb*100+4],paramslocal_array[bb*100+5],paramslocal_array[bb*100+6]);
    }
    
    double lcq = lcq_submit(numspots,m_Kw_sigmac,paramslocal_array);
    model_plot_array[jj] = lcq;
    //ok we are going to care whether they overlap in this region as well, be careful however
    //as this could lead to problems
    if(EnsureSpotsDontOverlap(t,printinfo_int,numspots,paramslocal_array,&degrees_overlap)==1){
      spots_overlap = 1;
      degrees_overlap_care = degrees_overlap;
    }
  }
  */
  chicare = 9e9;
  if(spots_overlap==1){
    m_SpotsOverlap = 1;
    //printf("ATTENTION YOUR SPOTS OVERLAP AT SOME POINT ");
    m_RepeatOverlap++;

    //CALCULATE CHI2 but add a penalty
    chicare = ChiSquared(Ntop,jd_array,flux_array,merr_array,model_array,printinfo_int,Degrees_of_Freedom);
    //here is the penalty, this should prevent them from overlapping
    chicare = 1000*(1+degrees_overlap_care) + chicare;
    //obviously overlap gives a bad chi2
  }
  //more common occurence
  else{
    m_SpotsOverlap = 0;
    if(printinfo_int==1){
      //      printf("SPOT Configuration does not overlap");
    }
    chicare = ChiSquared(Ntop,jd_array,flux_array,merr_array,model_array,printinfo_int,Degrees_of_Freedom);
  }
  if(printinfo_int==1){
    //    fprintf(Printz,"Reduced Chi-Squared: %0.11lf, Chi^2: %lf, Ntop: %d, free: %d,",chicare,chicare*Degrees_of_Freedom,Ntop,Degrees_of_Freedom);
  }
  
  //make_residuals(Ntop,Ntop_model,jd_array,flux_array,model_array,resid_array,flat_array);
  //*Ntop_model_return = Ntop_model;
  *ChiCare_return = chicare;
  //	UpdateThings();
}
///////////////////////////////////////////////////////////////////////////////
void NormalizeParams(double params_local[MAXsmall])
//This function finds if k and Peq are declared, then it determines
//the appropriate periods
{
  if(params_local[31] >= 0.0 && params_local[30] >= -9999){
    
    //then we are using the k and Peq fitting method
    //find the number of spots
    int numspots = int(params_local[11]);
    double k = params_local[30];   //differential rotation parameter 0:1?
    double peq = params_local[31]; //equatorital period (days)

    params_local[30] = k;
    params_local[31] = peq;
    for(int jj=1; jj<=numspots; jj++){
      double pjj = params_local[jj*100+3];  //spot period
      double bjj = asin(sqrt((1-peq/pjj)/k))*180./m_ONE_PI;
      params_local[jj*100+5] = bjj;
      //double bjj = params_local[jj*100+5]; //spot latitude
      //double pjj = peq/(1-k*pow(sin(bjj*m_ONE_PI/180.),2.));  //spot period
      //params_local[jj*100+3] = pjj;
    }
  }
}
///////////////////////////////////////////////////////////////////////////////
void CheckParams(double params_local[MAXsmall])
{
  //this information should be checked in the GA, this routine is nearly pointless
  
  int numspots = int(params_local[11]);  //number of spots
  //checkparams
  
  //HERE TAKE PARAMS AND MAKE THEM NORMAL
  
  //make sure the inclination is below 90 degrees
  if(params_local[3] < 0.0){params_local[3] = .1;}
  if(params_local[3] > 90.0){params_local[3] = 89.9;}
  //A STARS OFF
  double b,l,g;
  double p;

  //check b,l,g
  for(int ii=1; ii<=numspots; ii++){
    //
    b = params_local[ii*100+5];  //latitude -90:90
    l = params_local[ii*100+4];  //longitude 0:360
    g = params_local[ii*100+6];  //angular radius 0:90
    
    while(b > 270){b = b - 360;}
    while(b < -90){b = b + 360;}
    
    if( b > 90){
      b = 180-b;
      l = 180+l;
    }
    //what happens
    
    while(l > 360){l = l - 360;}
    while(l < 0){l = l + 360;}
    l = 0;
    if(l > 0){
      params_local[ii*100+2] = params_local[ii*100+2] - (l/360.0)*params_local[ii*100+3];
      l = 0;
    }
    else if( l < 0){
      params_local[ii*100+2] = params_local[ii*100+2] + (l/360.0)*params_local[ii*100+3];
      l = 0;
    }
    
    g = fabs(g);
    if(g > 90){g = 90;}
    //check the period is greater than zero
    //params_local[ii*100+3] = fabs(params_local[ii*100+3]);
    //HELP HELP
    //add a size constraint
    //if(g > 50){g = 50;}
    
    p = params_local[ii*100+3];  //spot period (days)
    p = fabs(p);
    
    //redefine corrected parameters
    params_local[ii*100+3] = p;
    params_local[ii*100+4] = l;
    params_local[ii*100+5] = b;
    params_local[ii*100+6] = g;
  }
  //NormalizeParams(params_local);
}
///////////////////////////////////////////////////////////////////////////////
//Finds the minimum and maximum values of time it the data set
void returnminmax(int ntop,double value[MAX], double *min, double *max)
{
  double minner = 9e99;
  double maxxer = -9e99;
  for (int j=0; j<ntop; j++){
    if(value[j] < minner){
      minner = value[j];
    }
    if(value[j] > maxxer){
      maxxer = value[j];
    }
  }
  *min = minner;
  *max = maxxer;
  //ntop is the total number of data points
}
///////////////////////////////////////////////////////////////////////////////
double sigmac_call(double t, double u, double E, double P, double i, double l, double b, double g)
{
  double sigmac;
  if(m_BuddingBool==0){
    //use budding's model
    sigmac = sigmac_budding(t,u,E,P,i,l,b,g);
  }
  else if(m_BuddingBool==1){
    //use dorren's model
    sigmac = sigmac_dorren(t,u,E,P,i,l,b,g);
  }
  return sigmac;
}
///////////////////////////////////////////////////////////////////////////////
//Budding model for spot position/limb darkening?
double sigmac_budding(double t, double u, double E, double P, double i, double l, double b, double g)
//THIS CODE OBTAINED NEARLY VERBATIM FROM SPOTMODEL
{
	//GNUPLOT 
	double D2R = 0.017453292519943295769236907684886127;
	double TWOPI = 6.28318530717958647692528676655900577;
	double PI_PER_2 =  1.57079632679489661923132169163975144;
	double PI = 3.14159265358979323846264338327950288;
	double SIGMA_2_PER_3 = 0.66666666666666666666666666666666667;

	double Phi, cosi, sini, cosb, sinb, sing, sings;
	double sings1, rsings1, singz0, z0, z0s, z0s1, ds, d, s, rss1, n, acosn, acoss, acosns;
	double rns1, sigma00, sigma10, bd;
	if (b>90) {
		b=180-b;
		l=180+l;
	}
	//l*=D2R;
	//b*=D2R;
	g = fabs(g);
	if (g>90){g=180-g;}
	Phi=TWOPI*(t-E)/P;

	//NOW WE SHOULD CONVERT DEGREES TO RADIANS
	i = i*m_ONE_PI/180.;
	l = l*m_ONE_PI/180.;
	b = b*m_ONE_PI/180.;
	g = g*m_ONE_PI/180.;

	cosi=cos(i);
	sini=sin(i);
	cosb=cos(b);
	sinb=sin(b);
	sing=sin(g);
	sings=sing*sing;
	sings1=1.0-sings;
	rsings1=sqrt(sings1);
	z0=cos(l-Phi)*cosb*sini+sinb*cosi;
	z0s=z0*z0;
	z0s1=1.0-z0s;
	singz0=sing*z0;
	ds=sings1*z0s1;
	d=sqrt(ds);
	bd=acos(z0); /* distance of the center of the visible hemisphere and the center of the spot */
	
	if (bd>=PI_PER_2+g) {
		/* out of view */
		sigma00=0.0;
		sigma10=0.0;
	}
	else {
		if (d<=sings1) {
			/* completely in view */
			sigma00=sings*z0;
			sigma10=SIGMA_2_PER_3*(1.0-rsings1*(sings1+(3.0*ds*sings)/(2.0*sings1)));
		}
		else{
			/* partially in view */
			s=sings1/d;
			rss1=sqrt(1.0-s*s);
			n=(d-s)/singz0;
			rns1=sqrt(1.0-n*n);
			acosn=acos(n);
			acoss=acos(s);
			acosns=acos(n/s);
			sigma00=(acoss-s*rss1+sings*z0*(acosn-n*rns1))/PI;
			sigma10=(2.0/(3.0*PI))*(acosns+(rsings1/(2.0*s))*(singz0*(3.0*sings-1.0)*rns1-(2.0*s*sings1+3.0*d*sings)*acosn));
		}
	}
	return((3.0/(3.0-u))*((1.0-u)*sigma00+u*sigma10));
}
///////////////////////////////////////////////////////////////////////////////
//Dorren model for spot position and better limb darkening
double sigmac_dorren(double t, double u, double E, double P, double i, double l, double b, double g)
//THIS CODE OBTAINED NEARLY VERBATIM FROM SPOTMODEL
{
	double Phi,cosi,sini,cosb,sinb,sigma00,sigma10,ad,cosbd,bd,sinbd,tanbd,sinad,cosad,tanad,cosdd,dd,coszd,sinzd,zd,T,PIdd,sinad2,cosad2,cosad3,sinbd2,sin2bd,A,B;
	double D2R = 0.017453292519943295769236907684886127;
	double TWOPI = 6.28318530717958647692528676655900577;
	double PI_PER_2 =  1.57079632679489661923132169163975144;
	double PI = 3.14159265358979323846264338327950288;
	double SIGMA_2_PER_3 = 0.66666666666666666666666666666666667;
	double SIGMA_1_PER_3 = 0.33333333333333333333333333333333333;
	double SIGMA_1_PER_6 = 0.16666666666666666666666666666666667;
	double SIGMA_1_PER_PI = 0.318309886183790671537767526745028724;
	double PI_PER_3 = PI/3.;
	
	if (b>90){
		b=180-b;
		l=180+l;
	}
	g = fabs(g);
	if (g>90){g=180-g;}

	Phi=TWOPI*(t-E)/P;

	//NOW WE SHOULD CONVERT DEGREES TO RADIANS
	i = i*m_ONE_PI/180.;
	l = l*m_ONE_PI/180.;
	b = b*m_ONE_PI/180.;
	g = g*m_ONE_PI/180.;

	cosi=cos(i);
	sini=sin(i);
	cosb=cos(b);
	sinb=sin(b);
	ad=g; /* radius of the spot */
	cosbd=cos(l-Phi)*cosb*sini+sinb*cosi; /* z0 */
	bd=acos(cosbd); /* distance of the center of the visible hemisphere and	the center of the spot */

	if (bd>=PI_PER_2+ad) {
		/* out of view */
		A=0.0;
		B=0.0;
	}
	else if (bd<=PI_PER_2-ad) {
		/* completely in view */
		sinbd=sqrt(1.0-cosbd*cosbd);
		sinad=sin(ad);
		cosad=cos(ad);
		sinad2=sinad*sinad;
		cosad3=cosad*cosad*cosad;
		sinbd2=sinbd*sinbd;
		A=PI*cosbd*sinad2;
		B=PI_PER_3*(2.0-2.0*cosad3-3.0*sinbd2*cosad*sinad2);
	}
	else {
		/* partially in view */
		sinbd=sqrt(1.0-cosbd*cosbd);
		tanbd=sinbd/cosbd;
		sinad=sin(ad);
		cosad=cos(ad);
		tanad=sinad/cosad;
		cosdd=1.0/(tanad*tanbd);
		dd=acos(cosdd);
		coszd=cosad/sinbd;
		sinzd=sqrt(1.0-coszd*coszd);
		zd=acos(coszd);
		if (bd<=PI_PER_2){
			T=atan(sinzd*tanbd);
		}
		else{
		    T=PI-atan(-sinzd*tanbd);
		}
		PIdd=PI-dd;
		sinad2=sinad*sinad;
		cosad2=cosad*cosad;
		cosad3=cosad2*cosad;
		sinbd2=sinbd*sinbd;
		sin2bd=2.0*sinbd*cosbd;
		A=zd+PIdd*cosbd*sinad2-sinzd*sinbd*cosad;
		B=SIGMA_1_PER_3*PIdd*(-2.0*cosad3-3.0*sinbd2*cosad*sinad2)+SIGMA_2_PER_3*(PI-T)+SIGMA_1_PER_6*sinzd*sin2bd*(2.0-3.0*cosad2);
	}
	sigma00=SIGMA_1_PER_PI*A;
	sigma10=SIGMA_1_PER_PI*B;
	return((3.0/(3.0-u))*((1.0-u)*sigma00+u*sigma10));
}
///////////////////////////////////////////////////////////////////////////////
//This seems to calculate how much light is attenuated by a spot
double lcq_submit(int nspots, double kw_sigmac[MAXsmall], double paramcare[MAXsmall])
{
  double U = paramcare[1];//gParam_set, unspotted intensity star/system
  double L1 = paramcare[2];//gParam_set, unspotted intensity of star
  double USI = paramcare[0];//unspotted magnitude, don't care
  //kw is the flux ratio of a spot
  double minus = 0;
  int kwfind;
  for(int bb=1; bb<=nspots; bb++){
    kwfind = bb*100;//kw listing is #00
    minus += (1-paramcare[kwfind])*kw_sigmac[bb];//spot flux * spot area
  }
  double lcq_care;
  
  //don't need this part
  if(m_gMagnitude==1){
    lcq_care = USI - 2.5*log10(U-L1*( minus )) ;
  }
  else{
    lcq_care = (U-L1*( minus )) ;//flux after total attenuation of all spots
  }
  return lcq_care;
}
///////////////////////////////////////////////////////////////////////////////
//obvious function here
int EnsureSpotsDontOverlap(double time, int print_info, int numspots,double paramslocal_array[MAX], double * degrees_overlap)
{
  //FILE *Printz;
  //char garbage[100] = "garbage";
  //Printz = fopen(garbage,"w");
  //if the spots don't overlap return 0
  int retval = 0;
  LBG spotA; LBG spotB;
  //lbg_return();
  for(int bb=1; bb<numspots; bb++){
    lbg_return(time, paramslocal_array[bb*100+3], paramslocal_array[bb*100+2], paramslocal_array[bb*100+4],paramslocal_array[bb*100+5],paramslocal_array[bb*100+6], &spotA);
    for(int jj=bb+1; jj<=numspots; jj++){
      lbg_return(time, paramslocal_array[jj*100+3], paramslocal_array[jj*100+2], paramslocal_array[jj*100+4],paramslocal_array[jj*100+5],paramslocal_array[jj*100+6], &spotB);
      if(CheckIfSpotsOverlap(spotA,spotB,&(*degrees_overlap))==1){
	if(print_info==1){
	  //fprintf(Printz,"Spots %d and %d overlap at time: %lf",bb,jj,time);
	}
	//the spots overlap
	retval = 1;
      }
    }
  }
  //fclose(Printz);
  return retval;
}
///////////////////////////////////////////////////////////////////////////////
//function for EnsureSpotsDontOverlap
void lbg_return(double t, double P, double E, double l, double b, double g,  LBG *spot)
{
  double ll = l;
  double bb = b;
  double gg = g;
  if (bb>90) {
    bb=180-bb;
    ll=180+ll;
  }
  gg = fabs(gg);
  if (gg>90){gg=180-gg;}
  double Phi=2*m_ONE_PI*(t-E)/P;
  //need to convert to radians
  ll = Phi + ll*m_ONE_PI/180.;
  bb = bb*m_ONE_PI/180.;
  gg = gg*m_ONE_PI/180.;
  spot->l = ll;
  spot->b = bb;
  spot->g = gg;
}
///////////////////////////////////////////////////////////////////////////////
//function for EnsureSpotsDontOverlap
int CheckIfSpotsOverlap(LBG spotlbg1, LBG spotlbg2, double * degrees_overlap)
{
  int retval = 0;
  *degrees_overlap = 0;
  //need to translate LBG into xyz
  XYZ spotxyz1;
  spotxyz1.x = cos(spotlbg1.b)*cos(spotlbg1.l);
  spotxyz1.y = sin(spotlbg1.b);
  spotxyz1.z = -cos(spotlbg1.b)*sin(spotlbg1.l);
  
  XYZ spotxyz2;
  spotxyz2.x = cos(spotlbg2.b)*cos(spotlbg2.l);
  spotxyz2.y = sin(spotlbg2.b);
  spotxyz2.z = -cos(spotlbg2.b)*sin(spotlbg2.l);
  
  double ang1 = acos(dotproduct(spotxyz1,spotxyz2));
  double ang2 = (spotlbg1.g + spotlbg2.g);
  if( ang1 < ang2){
    *degrees_overlap = (ang2-ang1)*180/m_ONE_PI;
    retval = 1;
  }
  return retval;
}

//could be useful
/*
double vsini_return(double P1_day,double Rstar_rsun,double inc_deg)
//this is for a solid-rotator
{
	double Peq_seconds = P1_day*24*3600;
	double Rsun = 6.9599e8;//m
	double Rstar_m = Rstar_rsun*Rsun;
	double eq_speed_kms = 1e-3*2*g_pDlg->m_ONE_PI*Rstar_m/Peq_seconds;//m/s
	double inc_rad = inc_deg*g_pDlg->m_ONE_PI/180.; 
	double vsini = eq_speed_kms*sin(inc_rad);
	return vsini;
}
*/
///////////////////////////////////////////////////////////////////////////////
void make_residuals(int ntop, int ntop2, double time[MAX], double fluxa[MAX], double modelflux[MAX], double mag_residuals[MAX], double flatline[MAX])
{
	for(int jj=0; jj<ntop; jj++){
		mag_residuals[jj] = fluxa[jj]-modelflux[jj];
	}
	for(int jj=0; jj<ntop2; jj++){
		flatline[jj] = 0.0;
	}
}
/*
//might need
void UpdateThings()
{
  updateperiod(m_gParams[103]);
  updatetime(m_gTime);
  if(m_StripWant){
    UpdateStrip();
  }
  if(m_StarWant){
    OnPaintSpot();
    SaveDIB( );
  }
}

//useful?
double return_mean(int ntop,double mean[MAX])
{

	double meanmean = -9e9;
	double ongoing = 0;

	if(ntop!=0){
		for(int i=0;i<ntop;i++){
			//printf("%d, value: %f\n",ntop,mean[i]);
			ongoing = ongoing + mean[i];
		}
		meanmean = ongoing/(ntop);
	}
	return meanmean;
}
//useful?
double return_std_dev(int ntop,double mean,double value[MAX])
{
	
	double std_deviation = -9999;
	double ongoing = 0;
	if(ntop!=0){
		for(int i=0;i<ntop;i++){
			ongoing = ongoing + pow(mean-value[i],2);
		}
		std_deviation = sqrt(ongoing/(ntop) );
	}
	return std_deviation;
}

//need this?
void return_k_and_peq(double params_local_array[MAX], double *k_return, double *peq1_return, double *peq2_return)
{
  FILE *Printz;
  char garbage[100] = "garbage";
  Printz = fopen(garbage,"w");
	double Peq1_day,Peq2_day,k;
	int numspots = int(params_local_array[11]);

	if(numspots==3){
		double P3 = params_local_array[303];
		double P2 = params_local_array[203];
		double P1 = params_local_array[103];
		double theta1 = params_local_array[105]*m_ONE_PI/180.;
		double theta2 = params_local_array[205]*m_ONE_PI/180.;
		double theta3 = params_local_array[305]*m_ONE_PI/180.;
		
		//double k12 = (P2-P1)/(P2*pow( sin(theta2),2.) - P1*pow( sin(theta1),2.) );
		//double k23 = (P3-P2)/(P3*pow( sin(theta3),2.) - P2*pow( sin(theta2),2.) );
		//double k13 = (P3-P1)/(P3*pow( sin(theta3),2.) - P1*pow( sin(theta1),2.) );
		//k = (k12+k23+k13)/double(3.);
		k = (P2-P1)/(P2*pow( sin(theta2),2.) - P1*pow( sin(theta1),2.) );
		Peq1_day = P1*(1-k*pow(sin(theta1),2.));
		Peq2_day = P2*(1-k*pow(sin(theta2),2.));
		//Peq3_day = P3*(1-k*pow(sin(theta3),2.));
		
	}
	else if(numspots==2){
		double P1 = params_local_array[103];
		double P2 = params_local_array[203];
		double theta1 = params_local_array[105]*m_ONE_PI/180.;
		double theta2 = params_local_array[205]*m_ONE_PI/180.;
		
		fprintf(g_pDlg->Printz,"P1: %lf, P2: %lf, theta1: %lf, theta2: %lf,",P1,P2,theta1,theta2);
		
		k = (P2-P1)/(P2*pow( sin(theta2),2.) - P1*pow( sin(theta1),2.) );
		Peq1_day = P1*(1-k*pow(sin(theta1),2.));
		Peq2_day = P2*(1-k*pow(sin(theta2),2.));
		
	}
	else if(numspots==1){
		Peq1_day = params_local_array[103];
		Peq2_day = 0.00;
		k = 0.00;
		
	}
	else{
		Peq1_day = params_local_array[103];
		Peq2_day = 0.00;
		k = 0.00;
		
	}
	params_local_array[30] = k;
	params_local_array[31] = Peq1_day;
	*k_return = k;
	*peq1_return = Peq1_day;
	*peq2_return = Peq2_day;
	fprintf(Printz,"kdif: %lf, Peq: %lf, numspots: %d,",k,Peq1_day,numspots);
	
}
*/
///////////////////////////////////////////////////////////////////////////////
//chi square analysis routines
double ChiSquared(int ntop, double jd_array[MAX], double mag_array[MAX], double merr_array[MAX], double model_array[MAX], int printwant_int, int Degrees_of_Freedom)
{
	double chicare = ChiSquaredUnReduced(ntop,jd_array,mag_array,merr_array,model_array,printwant_int);
	chicare = chicare/double(Degrees_of_Freedom);
	return chicare;
}
///////////////////////////////////////////////////////////////////////////////
double ChiSquaredUnReduced(int ntop, double jd_array[MAX], double mag_array[MAX], double merr_array[MAX], double model_array[MAX], int printwant_int)
{
	double chicare = 0.0;
	for(int jj=0; jj<ntop; jj++){
		chicare +=  pow( (mag_array[jj]- model_array[jj])/merr_array[jj],2.);
	}
	return chicare;
}
///////////////////////////////////////////////////////////////////////////////
double dotproduct(XYZ vec1, XYZ vec2)
{
	double dotter = vec1.x*vec2.x + vec1.y*vec2.y + vec1.z*vec2.z;
	return dotter;
}
///////////////////////////////////////////////////////////////////////////////
double ERAND(double top, double bottom)  //create random number between top and bottom values
{
  double out;
  out = ((double)rand() / randmax)*(top-bottom) + bottom;
  return out;
}
///////////////////////////////////////////////////////////////////////////////
double STDEV(double y[],int a) {
/*-------------------------------------------------------------
  Calculate the standard deviation of an array y, given a terms
  -------------------------------------------------------------*/

double ymean,sum;
int z;
  
ymean=0;
for(z=0;z<a;z++)
  ymean+=y[z];
ymean/=a;

sum=0; 
for(z=0;z<a;z++)
  sum+=pow((y[z]-ymean),2.0); 

return sqrt(sum/a);
}
///////////////////////////////////////////////////////////////////////////////
double GRAND(double center, double sigma)  // Create random gaussian number with center and sigma values
  {
  int z;
  double u1,u2,out,randmax = RAND_MAX;

  if(sigma < 0)
  sigma = fabs(sigma);

  z=0;
  while(z==0) // Create Gaussian Noise Array of sigma S
    {
    u1 = 2.0*random()/randmax;
    if(u1>1.0)
      u1 = u1 - 2.0;
    u2 = 2.0*random()/randmax;
    if(u2>1.0)
      u2 = u2 - 2.0;
    if(u1*u1+u2*u2>0 && u1*u1+u2*u2<1)
      {
      out = center + sigma*u1*sqrt(-2*log(u1*u1+u2*u2)/(u1*u1+u2*u2));
      z=1;
      }
    }
  return out;
  }
