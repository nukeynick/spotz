
//from StarSpotzDlg.h
#define MAXmed 16000
#define MAXsmall 1600
#define MAX 1600000

using namespace std;

//structs from StarSpotzDlg.h
//XYZ is used to check for spot overlap, comes from LBG
struct XYZ
{
    double x;
    double y;
    double z;
};
//RTP was used to create the 3D image of the star
struct RTP
{
    double r;
    double t;
    double p;
};
//LBG defines the location and size of the spots
struct LBG
{
    double l;
    double b;
    double g;
};

void runmodel_special(int printinfo_int, double paramslocal_array[MAXmed], int Ntop, double jd_array[MAX], double flux_array[MAX], double merr_array[MAX], double model_array[MAX], int *Ntop_model_return, double jd_model_array[MAX], double model_plot_array[MAX], double resid_array[MAX], double flat_array[MAX], double Timeinc, int Degrees_of_Freedom, int numspots, double TimeMin, double TimeMax, double *ChiCare_return, double param, double m_param[MAXsmall]), NormalizeParams(double paramslocal_array[MAXmed]), CheckParams(double params_local[MAXmed]), returnminmax(int ntop,double value[MAX], double *min, double *max), lbg_return(double t, double P, double E, double l, double b, double g,  LBG *spot), make_residuals(int ntop, int ntop2, double time[MAX], double fluxa[MAX], double modelflux[MAX], double mag_residuals[MAX], double flatline[MAX]), UpdateThings(), return_k_and_peq(double params_local_array[MAXmed], double *k_return, double *peq1_return, double *peq2_return);

double sigmac_call(double t, double u, double E, double P, double i, double l, double b, double g), sigmac_budding(double t, double u, double E, double P, double i, double l, double b, double g), sigmac_dorren(double t, double u, double E, double P, double i, double l, double b, double g), lcq_submit(int nspots, double kw_sigmac[MAXsmall], double paramcare[MAXsmall]), vsini_return(double P1_day,double Rstar_rsun,double inc_deg), return_mean(int ntop,double mean[MAX]), return_std_dev(int ntop,double mean,double value[MAX]), ChiSquared(int ntop, double jd_array[MAX], double mag_array[MAX], double merr_array[MAX], double model_array[MAX], int printwant_int, int Degrees_of_Freedom), ChiSquaredUnReduced(int ntop, double jd_array[MAX], double mag_array[MAX], double merr_array[MAX], double model_array[MAX], int printwant_int), dotproduct(XYZ vec1, XYZ vec2);

int EnsureSpotsDontOverlap(double time, int print_info, int numspots,double paramslocal_array[MAXmed], double * degrees_overlap), CheckIfSpotsOverlap(LBG spotlbg1, LBG spotlbg2, double * degrees_overlap);

double m_ChiCare,m_ChiStart,jd,mag,err,jdmin,jdmax,PhaseMin,TimeMin,TimeMax,Timeinc,degrees_overlap,spots_overlap,degrees_overlap_care,bump,error,param;
char Printz[512];
char PrintChar[256];
int printinfo_int,m_Ntop,m_Ntop_model,m_Fitparams_int,m_RepeatOverlap,m_SpotsOverlap,i,m_Degrees_of_Freedom,numspots;
double chicare;
double m_gParams[MAXsmall];
double m_gParams_Errors[MAXsmall];
double m_gParams_new[MAXsmall];
double m_gParams_newest[MAXsmall];
int m_Fitparams[MAXsmall];
double m_Kw_sigmac[MAXsmall];
double m_modelTISMO[MAXmed];//after the fitting
double m_Phase[MAXmed];
double m_resid[MAXmed];
double m_flat[MAXmed];
double m_flux[MAXmed];
double m_jd[MAXmed];
double m_jd_model[MAXmed];
double m_model[MAXmed];
double m_model_plot[MAXmed];
double jd_model_array[MAXmed];
double model_plot_array[MAXmed];
double m_merr[MAXmed];
double paramslocal_array[MAXmed];
double m_param[MAXsmall];

//AGA initialization
double ERAND(double top, double bottom), GRAND(double center, double sigma), STDEV(double y[MAX],int a);

double spreadinit,centinit,wsum,randc;
int nvar,j,m,kk,n,time1,time2,nn;
double spread[MAXmed];
double stdevvals[MAXmed];
double help[MAXmed];
double wvals[MAXmed];
double Jparam[MAXmed];
