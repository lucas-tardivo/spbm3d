
#include <iostream>
#include <algorithm>
#include <math.h>
#include <mex.h>

#include "utilities.h"
#include "lib_transforms.h"
#include "libdistances_bm3d.h"
#include "lib_block_matching.h"

 using namespace std;

void dump_table_distance(std::vector<pair<float, unsigned> > table_distance);
 
 bool ComparaisonFirst(pair<float,unsigned> pair1, pair<float,unsigned> pair2)
{
	return pair1.first < pair2.first;
} 
 
extern int print_distance;

void dump_table_distance(vector<pair<float, unsigned> > table_distance) {
    int ii = table_distance.size();
} 

block_matching_ftype get_block_matching_function(char * pStr) {
    if (strcmp(pStr,"sd")==0) {
        return &precompute_BM_sd;            
    } else {
        mexPrintf("<%s>\n", pStr);
        mexErrMsgTxt("Invalid block maching function.");
    }   
}

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
){
    //! Declarations
    const unsigned  Ns    = 2 * nHW + 1;
    int disp_msg          = 0;
   
    if (patch_table.size() != width * height)
        patch_table.resize(width * height);
    
    vector<unsigned> row_ind;
    ind_initialize(row_ind, height - kHW + 1, nHW, pHW);
    
    vector<unsigned> column_ind;
    ind_initialize(column_ind, width - kHW + 1, nHW, pHW);

    //A threshold based on the distance, if you want to use a threshold
    float threshold = parms.distance_function(0.1, parms.sd_thr, parms.s, parms.b);
    
    if (parms.sd_thr < 0) {
        threshold = parms.sd_thr * -1;
    }
           
    vector<pair<float, unsigned> > table_distance;
    
    //! To avoid reallocation
    table_distance.reserve(Ns * Ns);
    
   
    vector<float> table_mle_sig2(width*height,-1.0f);
    
    for (unsigned ind_i = 0; ind_i < row_ind.size(); ind_i++)
    {
        for (unsigned ind_j = 0; ind_j < column_ind.size(); ind_j++)
        {
            //! Initialization
            const unsigned k_r = row_ind[ind_i] * width + column_ind[ind_j];
            table_distance.clear();
            patch_table[k_r].clear();


            if (table_mle_sig2[k_r] < 0)
                table_mle_sig2[k_r] = parms.mle_function(img, k_r, kHW, width, ind_i * width + ind_j, parms.filtered_img);
            
            for (int dj = -(int)nHW; dj <= (int) nHW; dj++) {
                for (int di = -(int)nHW; di <= (int) nHW; di++) {
                    const unsigned int k_r0 = k_r + di * width + dj;
                    
                    if (table_mle_sig2[k_r0] < 0)
                        table_mle_sig2[k_r0] = parms.mle_function(img, k_r0, kHW, width, ind_i * width + ind_j, parms.filtered_img);
                    
                    float dist = 0.0f;
                            
                    if (k_r != k_r0)
                        dist = parms.distance_function(table_mle_sig2[k_r], table_mle_sig2[k_r0], parms.s, parms.b);                   
                    
                    if (dist < threshold) {
                        table_distance.push_back(make_pair(dist, k_r0));
                    }
                }
            }        
            
            //! We need a power of 2 for the number of similar patches,
            //! because of the Welsh-Hadamard transform on the third dimension.
            //! We assume that NHW is already a power of 2
            const unsigned nSx_r = (NHW > table_distance.size() ?
                                    closest_power_of_2(table_distance.size()) : NHW);
                                    
			if (nSx_r == 1 && table_distance.size() == 0)
			{
				table_distance.push_back(make_pair(0, k_r));
			}
            //! Sort patches according to their distance to the reference one
            partial_sort(   table_distance.begin(), 
                            table_distance.begin() + nSx_r,
                            table_distance.end(), 
                            ComparaisonFirst);        
            
            //! Keep a maximum of NHW similar patches
            for (unsigned n = 0; n < nSx_r; n++)
                patch_table[k_r].push_back(table_distance[n].second);

			//! To avoid problem
			if (nSx_r == 1)
				patch_table[k_r].push_back(table_distance[0].second);
        }
    }
    dump_table_distance(table_distance);
}

