
double MEAN(int a, double x[MAX], double xerr[MAX]), MEANERR(int a, double x[MAX], double xerr[MAX]), STDEV(double y[MAX],int a), SUM(int a, double x[MAX]), SUMERR(int a, double xerr[MAX]), GRAND(double center, double sigma), WSTDEV(double y[MAX],double yerr[MAX], int a);

////////////////////////////////////////////////////////////////////////////////////////

double MEAN(int a, double x[], double xerr[]) {  // Return error-weighted mean of first a terms

int z;
double invvarsum,xmean;

invvarsum=0.0;
for(z=0;z<a;z++)
  invvarsum += 1.0/(xerr[z]*xerr[z]);
xmean=0;
for(z=0;z<a;z++)
  xmean += x[z]/(xerr[z]*xerr[z]);

return(xmean/invvarsum);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

double MEANERR(int a, double x[], double xerr[]) {  // Return the error of the error-weighted mean of first a terms

int z;
double invvarsum;

nvvarsum=0.0;
for(z=0;z<a;z++)
  invvarsum += 1.0/(xerr[z]*xerr[z]);

return(sqrt(1.0/invvarsum));
}

/////////////////////////////////////////////////////////////////////////////

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


////////////////////////////////////////////////////////////////////////////////////////////////////

double SUM(int a, double x[]) {  // Return sum of first a terms

int z;
double xsum;

xsum=0;
for(z=0;z<a;z++)
  xsum += x[z];

return(xsum);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

double SUMERR(int a, double xerr[]) {  // Return the error of the first a terms added in quadrature

int z;
double errsum;

errsum=0.0;
for(z=0;z<a;z++)
  errsum += (xerr[z]*xerr[z]);

return(sqrt(errsum));
}

//////////////////////////////////////////////////////////////////////////

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

///////////////////////////////////////////////////////////////////////////////

double WSTDEV(double y[],double yerr[], int a) {  // Return error-weighted standard deviation for an array y, with errors yerr, and number of terms a

double ymean,invvarsum,wstdev;
double wsum;
int z;

invvarsum=0;
for(z=0;z<a;z++)
  invvarsum = invvarsum + 1.0/(yerr[z]*yerr[z]);  // sum of weights

ymean=0;
for(z=0;z<a;z++)
  ymean=ymean + y[z]/(yerr[z]*yerr[z]);
ymean=ymean/invvarsum;

wsum=0;
for(z=0;z<a;z++)
  wsum = wsum + (1.0/(yerr[z]*yerr[z]))*pow((y[z]-ymean),2.0);

wstdev = sqrt(wsum/((a-1)*invvarsum/a));

return wstdev;
}

///////////////////////////////////////////////////////////////////////////////


