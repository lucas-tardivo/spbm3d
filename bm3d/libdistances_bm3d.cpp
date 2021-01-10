#include "mex.h"
#include <math.h>
#include <string.h>
#include "libdistances_bm3d.h"
#include <boost/math/distributions/chi_squared.hpp>

#include "lib_bm3d_mle.h"


#define min_val 0.0001f
#define min_sig 0.000000000001f

extern double *fExpPositiveLut;
#include <gsl/gsl_sf_hyperg.h> //Gauss Hypergeometric function

using namespace std;
#define Hypergeometric2F1(a,b,c,d) gsl_sf_hyperg_2F1(a,b,c,d)


//Reference https://www.physicsforums.com/threads/calculating-hypergeometric-function-2f1-for-z-1.501956/
double hyperg_z_GT1 (double a, double b, double c, double z) { 
//calculates 2F1 for z < -1 
    double coef1=tgamma(c)*tgamma(b-a)*pow(1-z,-a)/(tgamma(b)*tgamma(c-a)); 
    double coef2=tgamma(c)*tgamma(a-b)*pow(1-z,-b)/(tgamma(a)*tgamma(c-b)); 
    return coef1*gsl_sf_hyperg_2F1(a,c-b,a-b+1,1/(1-z))+coef2*gsl_sf_hyperg_2F1(b,c-a,b-a+1,1/(1-z)); 
}

double euclidean_distance(double a1, double a2, double s, double b) {
    if (a1 == a2) return 0.0f;
    
    double d = sqrt(pow((a1 - a2), 2));
    
    if (!isfinite(d)) {
        d = 3.0f; //If not finite, return any very high value
    } else if (d<0.0f && fabs(d)<min_val) {
        d = fabs(d);
    }
    return d;   
}

/************** POISSON DISTANCES ********************************/
double kullback_leibler_distance(double a1, double a2, double s, double b) {
    if (a1 == a2) return 0.0f;
    
    double d = (0.5f * (a1 - a2) * log(a1 / a2));
    
    if (!isfinite(d)) {
        d = 3.0f; //If not finite, return any very high value
    } else if (d<0.0f && fabs(d)<min_val) {
        d = fabs(d);
    }
    return d;   
}

double hellinger_distance(double a1, double a2, double s, double b) {
    if (a1 == a2) return 0.0f;
    
    double d = 1 - exp((-0.5 * (a1 + a2)) + sqrt(a1 * a2));
    
    if (!isfinite(d)) {
        d = 3.0f; //If not finite, return any very high value
    } else if (d<0.0f && fabs(d)<min_val) {
        d = fabs(d);
    }
    return d;   
}

double bhattacharryya_distance(double a1, double a2, double s, double b) {
    if (a1 == a2) return 0.0f;
    
    double d = 0.5 * (a1 + a2) - sqrt(a1 * a2);
    
    if (!isfinite(d)) {
        d = 3.0f; //If not finite, return any very high value
    } else if (d<0.0f && fabs(d)<min_val) {
        d = fabs(d);
    }
    return d;   
}

double renyi_distance(double a1, double a2, double s, double b) {
    if (a1 == a2) return 0.0f;
    
    double d = (1 / (b - 1)) * log(0.5 * ((exp((pow(a2,(1 - b))) * (pow(a1,b)) - (((1 - b) * a2) + (b * a1)))) + 
            (exp((pow(a1,(1 - b))) * (pow(a2,b)) - (((1 - b) * a1) + (b * a2))))));
    
    if (!isfinite(d)) {
        d = 3.0f; //If not finite, return any very high value
    } else if (d<0.0f && fabs(d)<min_val) {
        d = fabs(d);
    }
    return d;   
}

sd_distance get_distance_function(char * pStr) {
    if (strcmp(pStr,"kullback_leibler_distance")==0) {
         return &kullback_leibler_distance;                   
    } else  if (strcmp(pStr,"hellinger_distance")==0) {
         return &hellinger_distance;                   
    } else  if (strcmp(pStr,"bhattacharryya_distance")==0) {
         return &bhattacharryya_distance;                   
    }  else  if (strcmp(pStr,"renyi_distance")==0) {
         return &renyi_distance;                   
    }   else  if (strcmp(pStr,"euclidean_distance")==0) {
         return &euclidean_distance;                   
    } else {
         mexPrintf("<%s>\n", pStr);
         mexErrMsgTxt("Invalid distance function.");
    }    
}

/******** STATISTICS ***********************************************************/
double bhattacharyya_statistic(double dist, int M, int N,double s) {
    return ((8.0f*M*N)/(M+N))*dist;
}

double hellinger_statistic(double dist, int M, int N,double s) {
    return ((8.0f*M*N)/(M+N))*dist;
}

double kl_statistic(double dist, int M, int N,double s) {
    return ((2.0f*M*N)/(M+N))*dist;
}

double renyi_statistic(double dist, int M, int N,double s) {
    return ((1.0f/s)*(2.0f*M*N)/(M+N))*dist;
}

/*******************************************************************************/

double get_pvalue(double stat) {
    boost::math::chi_squared chi2dist(1);
    return (1.0f-boost::math::cdf(chi2dist,stat));
}

stat_ftype get_stat_function(char * pStr) {
         /************** FISHER-TIPPETT *************************/
    if (strcmp(pStr,"bhattacharyya_statistic")==0) {
        return &bhattacharyya_statistic; 
    } else  if (strcmp(pStr,"hellinger_statistic")==0) {
        return &hellinger_statistic;          
    } else  if (strcmp(pStr,"kl_statistic")==0) {
        return &kl_statistic;     
    } else  if (strcmp(pStr,"renyi_statistic")==0) {
        return &renyi_statistic;                       
    } else {
         mexErrMsgTxt("Invalid statistic function.");
    }    
}

