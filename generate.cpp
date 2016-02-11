#include <cmath>
#include <stdarg.h>	
#include <stdio.h>
#include <cstdlib>
#include <time.h>
#include "gen.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <omp.h>

using namespace std;

//these values need to be global parameters
//Spot model constants
const double m_Timeinc = .16;  //model will have five times as many data points
const double m_TimeMin = 0.0;
const double m_TimeMax = 40.0;
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
const int m_Number_Fitted_Params = 9;
const int m_Number_Of_Steps_Run_Fit = 0;
const double m_ONE_PI = 3.1415926535897932384626433832795;
const double randmax = RAND_MAX;

/*
paramslocal[#]
0 = USI unspotted magnitude of the star.
1 = U unspotted intensity of the star/system.
2 = L1 unspotted intensity of the spotted star.
3 = i inclination of the stellar rotation axis to the line of sight.
11 = numspots number of spots.
12 = RStar the radius of the star in solar radii.
13 = PHASEmin the eopch of minimum phase (HJD).
14 = NUMFit = the number of fitted parameters
30 = k differential rotation parameter
31 = Peq equatorital period (days).
#00 = kw# the flux ratio between spot # and the unspotted photosphere (0-1)
#01 = u# linear limb darkening coefficient (0-1).
#02 = E# spot epoch in HJD.
#03 = p# spot period in days.
#04 = l# spot longitude (0-360), lambda.
#05 = b# spot latitude (-90-90), beta.
#06 = g# spot angular radius (0-90), gamma.
*/

main() {
  FILE *fp;  //input file
  ofstream outfile1 ("model.params.dat",ofstream::binary);
  ofstream outfile2 ("model.data.dat",ofstream::binary);
  srand ( time(NULL) );
  printinfo_int = 0;

  char infile[100] = "temp.dat";  //model params from lightcurve.pro
  char new_line[100];
  i = 0;
  if ((fp = fopen(infile,"r")) == NULL)  {
    cout << "Could not open file for reading: " << infile;
    abort();
  }
  while (fgets(new_line, sizeof(new_line), fp))  {
    if (sscanf(new_line,"%lf",&param) == 1)
      {
        m_param[i] = param;
        i++;
      }
  }
  
  paramslocal_array[0] = m_param[0];
  paramslocal_array[1] = m_param[0];
  paramslocal_array[2] = m_param[0];
  paramslocal_array[3] = m_param[1];//inclination
  paramslocal_array[11] = 3;//numspots
  paramslocal_array[12] = m_param[2];//Rstar_solar
  paramslocal_array[13] = 0.0;//PhaseMin
  paramslocal_array[14] = 4*paramslocal_array[11];//NumFit
  paramslocal_array[30] = m_param[3];//k
  paramslocal_array[31] = m_param[4];//Peq
  numspots = paramslocal_array[11];

  paramslocal_array[100] = .3;
  paramslocal_array[100+1] = .684;
  paramslocal_array[100+2] = m_param[5];//epoch
  paramslocal_array[100+3] = m_param[6];//per
  paramslocal_array[100+4] = 0.0;
  paramslocal_array[100+5] = m_param[7];//lat
  paramslocal_array[100+6] = m_param[8];//rad
  paramslocal_array[2*100] = .3;
  paramslocal_array[2*100+1] = .684;
  paramslocal_array[2*100+2] = m_param[9];
  paramslocal_array[2*100+3] = m_param[10];
  paramslocal_array[2*100+4] = 0.0;
  paramslocal_array[2*100+5] = m_param[11];
  paramslocal_array[2*100+6] = m_param[12];
  paramslocal_array[3*100] = .3;
  paramslocal_array[3*100+1] = .684;
  paramslocal_array[3*100+2] = m_param[13];
  paramslocal_array[3*100+3] = m_param[14];
  paramslocal_array[3*100+4] = 0.0;
  paramslocal_array[3*100+5] = m_param[15];
  paramslocal_array[3*100+6] = m_param[16];

  /*
  for(int bb=1; bb<=numspots; bb++){
    paramslocal_array[bb*100] = .3;//kw
    paramslocal_array[bb*100+1] = .684;//u
    paramslocal_array[bb*100+2] = ERAND(10,0);//E1 spot epoch //bb*(paramslocal_array[31]/2.);
    //paramslocal_array[bb*100+3] = paramslocal_array[31];//spot period
    while(paramslocal_array[bb*100+3] < paramslocal_array[31]) paramslocal_array[bb*100+3] = GRAND(paramslocal_array[31],.5);//spot period
    paramslocal_array[bb*100+4] = 0.0;//spot longitude
    paramslocal_array[bb*100+5] = ERAND(90,-90);//spot latitude
    paramslocal_array[bb*100+6] = ERAND(30,5);//spot radius (angular)
    }
  */
  TimeMax=m_TimeMax;
  TimeMin=m_TimeMin;
  Timeinc=m_Timeinc;
  outfile1 << "inclination:" << " " << paramslocal_array[3] << endl;
  outfile1 << "period:"  << " " << paramslocal_array[31] << endl;
  for(int bb=1; bb<=numspots; bb++) {
    outfile1 << "epoch" << bb << ":" << " " << paramslocal_array[bb*100+2] << endl;
    outfile1 << "period" << bb << ":" << " " << paramslocal_array[bb*100+3] << endl;
    outfile1 << "latitude" << bb << ":" << " " << paramslocal_array[bb*100+5] << endl;
    outfile1 << "spot size" << bb << ":" << " " << paramslocal_array[bb*100+6] << endl;
  }
  outfile1.close();

  double time_step = Timeinc;
  int Ntop_model = int((TimeMax-TimeMin)/time_step);
  if(Ntop_model > MAXmed){
    Ntop_model = MAXmed;
    time_step = (TimeMax-TimeMin)/double(Ntop_model);
  }
  for(int jj=0; jj<Ntop_model; jj++){
    double t = TimeMin +(jj)*time_step;
    jd_model_array[jj] = t;
    
    //double sigmac1 = sigmac_b(t,U,Epoch1,p1,inc,l1,b1,g1);
    for(int bb=0; bb<=numspots; bb++){
      m_Kw_sigmac[bb] = sigmac_call(t,paramslocal_array[bb*100+1],paramslocal_array[bb*100+2],paramslocal_array[bb*100+3],paramslocal_array[3],paramslocal_array[bb*100+4],paramslocal_array[bb*100+5],paramslocal_array[bb*100+6]);
    }
    
    double lcq = lcq_submit(numspots,m_Kw_sigmac,paramslocal_array);
    model_plot_array[jj] = lcq;

    if(EnsureSpotsDontOverlap(t,printinfo_int,numspots,paramslocal_array,&degrees_overlap)==1){
      spots_overlap = 1;
      cout << "SPOTS OVERLAP, DATA INVALID!" << endl;
      degrees_overlap_care = degrees_overlap;
    }
  }
  for(int jj=0;jj<Ntop_model; jj++) {
    //    bump = GRAND(0.0,.0005);  //random noise generator
    error = .002;
    if(error<0) error = -1*error;
    outfile2 << jd_model_array[jj] << " " << bump+model_plot_array[jj] << " " <<error << endl;
  }
  outfile2.close();
  return 0;
}

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

//obvious function here
int EnsureSpotsDontOverlap(double time, int print_info, int numspots,double paramslocal_array[MAX], double * degrees_overlap)
{
  FILE *Printz;
  //  char garbage[100] = "garbage.dat";
  //  Printz = fopen(garbage,"w");
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
	  fprintf(Printz,"Spots %d and %d overlap at time: %lf",bb,jj,time);
	}
	//the spots overlap
	retval = 1;
      }
    }
  }
  //fclose(Printz);
  return retval;
}

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

double ERAND(double top, double bottom)  //create random number between top and bottom values
{
  double out;
  out = ((double)rand() / randmax)*(top-bottom) + bottom;
  return out;
}

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

double dotproduct(XYZ vec1, XYZ vec2)
{
	double dotter = vec1.x*vec2.x + vec1.y*vec2.y + vec1.z*vec2.z;
	return dotter;
}
