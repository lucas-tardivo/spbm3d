//mex validate_parameters_input.cpp
#include <math.h>
#include <vector>
#include <mex.h>

#include "libdistances_bm3d.h"
#include "lib_block_matching.h"

using namespace std;

void get_parameters(block_matching_ftype &blk_method, sdbm3dPtype &parms, const mxArray *pm) {  
    
    mxArray *tmp = mxGetField(pm, 0, "use_std");
    if (tmp != NULL)
        parms.use_std = (int)mxGetScalar(tmp);   
    else {
        mexErrMsgTxt("Could not read use_std parameter!!!");        
    }
    
    tmp = mxGetField(pm, 0, "s");
    if (tmp != NULL) {
        parms.s = (float)mxGetScalar(tmp);   
    } else {
        mexErrMsgTxt("Could not read s parameter!!!");        
    }
    
    tmp = mxGetField(pm, 0, "b");
    if (tmp != NULL) {
        parms.b = (float)mxGetScalar(tmp);   
    } else {
        mexErrMsgTxt("Could not read b parameter!!!");        
    }
    
    tmp = mxGetField(pm, 0, "transform");
    if (tmp != NULL) {
        int trans = (int)mxGetScalar(tmp);
        if (trans == 0) {
            parms.transform = 4; //DCT   
        } else if (trans == 1) {
            parms.transform = 5; //BIOR      
        } else {
            mexErrMsgTxt("Invalid transform!!!");     
        }
    } else {
        mexErrMsgTxt("Could not read transform parameter!!!");         
    }
    
    tmp = mxGetField(pm, 0, "lambdaHard3D");
    if (tmp != NULL) {
        parms.lambdaHard3D = (float)mxGetScalar(tmp);  
    } else {
        mexErrMsgTxt("Could not read lambdaHard3D parameter!!!");        
    }        
    
    tmp = mxGetField(pm, 0, "tauMatch");
    if (tmp != NULL)
        parms.tauMatch = (float)mxGetScalar(tmp);   
    else {
        mexErrMsgTxt("Could not read tauMatch parameter!!!");        
    }
    
    tmp = mxGetField(pm, 0, "sd_thr");
    if (tmp != NULL)
        parms.sd_thr = (float)mxGetScalar(tmp);   
    else {
        mexErrMsgTxt("Could not read sd_thr parameter!!!");        
    }
    
    tmp = mxGetField(pm, 0, "mle_func");
    if (tmp != NULL) {
        int lenStr = mxGetN(tmp); 
        parms.mle_function_name = new char[lenStr+1];
        mxGetString(tmp, parms.mle_function_name, lenStr+1); 
        parms.mle_function = get_mle_function(parms.mle_function_name);        
    } else {
        mexErrMsgTxt("Could not read mle function parameter!!!");        
    } 
    
    tmp = mxGetField(pm, 0, "distance_func");
    if (tmp != NULL) {
        int lenStr = mxGetN(tmp); 
        parms.distance_function_name = new char[lenStr+1];
        mxGetString(tmp, parms.distance_function_name, lenStr+1); 
        parms.distance_function = get_distance_function(parms.distance_function_name);        
    } else {
        mexErrMsgTxt("Could not read distance function parameter!!!");        
    } 
    
    tmp = mxGetField(pm, 0, "blk_func");
    if (tmp != NULL) {
        int lenStr = mxGetN(tmp); 
        parms.block_matching_function_name = new char[lenStr+1];
        mxGetString(tmp, parms.block_matching_function_name, lenStr+1); 
        blk_method = get_block_matching_function(parms.block_matching_function_name);        
    } else {
        mexErrMsgTxt("Could not read block_match_function parameter!!!");        
    }
}

void dump_parameters(sdbm3dPtype parms, int step) {
    mexPrintf("---------------------------------\n");    
    mexPrintf("Parameters for step<%d>\n", step);    
    mexPrintf(" *use std<%d>\n", parms.use_std);
    if (parms.transform == 4) {
        mexPrintf(" *transform<DCT>\n");
    } else { 
        if (parms.transform == 5) {
            mexPrintf(" *transform<BIOR>\n");
        } else {
            mexErrMsgTxt("ERROR with transform parameter!!!");
        }
    }    
    mexPrintf(" *lambdaHard3D<%f>\n", parms.lambdaHard3D);  
    mexPrintf(" *tauMatch<%f>\n", parms.tauMatch);      
    mexPrintf(" *distance_function<%s>\n", parms.distance_function_name);  
    mexPrintf(" *mle_function<%s>\n", parms.mle_function_name);  
    mexPrintf(" *block_match_function<%s>\n", parms.block_matching_function_name); 
    mexPrintf(" *sd_thr<%f>\n", parms.sd_thr);    
    mexPrintf(" *s<%f>\n", parms.s);    
    mexPrintf("---------------------------------\n\n");     
}

