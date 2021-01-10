
#ifndef _LIBDISTANCES_BM3D_H_
#define _LIBDISTANCES_BM3D_H_

typedef double (*sd_distance)(double,double,double,double);
typedef double (*stat_ftype)(double,int, int,double);
typedef double (*speckle_index)(double,double);

/******* DISTANCES *******************************************************************************/
double euclidean_distance(double a1, double a2, double s, double b);

/************** POISSON DISTANCES *********************************************************/
double kullback_leibler_distance(double a1, double a2, double s, double b);
double hellinger_distance(double a1, double a2, double s, double b);
double bhattacharryya_distance(double a1, double a2, double s, double b);
double renyi_distance(double a1, double a2, double s, double b);

sd_distance get_distance_function(char * pStr);
stat_ftype get_stat_function(char * pStr);
speckle_index get_speckle_index_function(char * pStr);

double bhattacharyya_statistic(double dist, int M, int N,double beta);
double hellinger_statistic(double dist, int M, int N,double beta);
double kl_statistic(double dist, int M, int N,double beta);
double renyi_statistic(double dist, int M, int N,double beta);

double get_pvalue(double stat);

#endif


