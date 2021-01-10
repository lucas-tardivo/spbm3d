#ifndef LIB_BLOCK_MATCHING_H_INCLUDED
#define LIB_BLOCK_MATCHING_H_INCLUDED

#include "fftw3.h"
#include "lib_bm3d_mle.h"
#include "libdistances_bm3d.h"

 using namespace std;

typedef struct {
    sd_distance             distance_function;
	mle_estimator           mle_function;                
    char*                   distance_function_name;
	char*                   mle_function_name;                
    char*                   block_matching_function_name;  
    float                   s ; 
    float                   sd_thr ; 
    int                     use_std;
    int                     transform;
    float                   lambdaHard3D;
    float                   tauMatch;
    vector<float>    filtered_img;
    float                   b;
}sdbm3dPtype;

typedef void (*block_matching_ftype)(    
    std::vector<std::vector<unsigned> > &patch_table
,   const std::vector<float> &img      
,   const unsigned width
,   const unsigned height
,   const unsigned kHW
,   const unsigned NHW
,   const unsigned nHW
,   const unsigned pHW
,   const float    tauMatch_ec
,   const float    sigma
,   const          sdbm3dPtype parms    
);

block_matching_ftype get_block_matching_function(char * pStr);

void precompute_BM_sd(
    std::vector<std::vector<unsigned> > &patch_table
,   const std::vector<float> &img      
,   const unsigned width
,   const unsigned height
,   const unsigned kHW
,   const unsigned NHW
,   const unsigned nHW
,   const unsigned pHW
,   const float    tauMatch_ec
,   const float    sigma
,   const          sdbm3dPtype parms          
);

#endif // LIB_BLOCK_MATCHING_H_INCLUDED
