#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <string.h>

#include <mex.h>
#include "sd_bm3d.h"
#include "utilities.h"
#include "lib_block_matching.h"
#include "libparameters.h"

#define YUV       0
#define YCBCR     1
#define OPP       2
#define RGB       3
#define DCT       4
#define BIOR      5
#define HADAMARD  6
#define NONE      7

#define IN_IMG         prhs[0]
#define FILTERED_IMG   prhs[4]
#define OUT_IMG        plhs[0]
#define OUT_BASIC_IMG  plhs[1]

using namespace std;

/**
 * @file   main.cpp
 * @brief  Main executable file. Do not use lib_fftw to
 *         process DCT.
 *
 * @author MARC LEBRUN  <marc.lebrun@cmla.ens-cachan.fr>
 */


void mexFunction(  int nlhs, mxArray *plhs[],         /* Output variables */
                    int nrhs, const mxArray *prhs[])   /* Input variables  */
{
	//! Declarations
	vector<float> img_noisy, img_basic, img_denoised, img_filtered;
    sdbm3dPtype parmsStep1, parmsStep2;
    block_matching_ftype blk_method1, blk_method2;
    unsigned width, height, chnls;
    
   /* Check for proper number of arguments */
    if (nrhs != 5) {
        mexErrMsgIdAndTxt("MATLAB:mexcpp:nargin", 
                "SD_BM3D requires 5 input arguments.");
    }       
    if (!mxIsDouble(IN_IMG)) {
        mexErrMsgIdAndTxt("MATLAB:mexcpp:nargin", 
                "SD_BM3D requires an double type input image.");
    } 
    
    //Get noise image
    chnls = 1;
    height = mxGetM(IN_IMG);
    width  = mxGetN(IN_IMG);
    double *fpMatInput = (double*)mxGetPr(IN_IMG);
    
    //Get filtered image
    double *fpFilteredMatInput = (double*)mxGetPr(FILTERED_IMG);
    
    //Sigma
    float fSigma = (float)mxGetScalar(prhs[1]);
    
    //Step parameters
    get_parameters(blk_method1, parmsStep1, prhs[2]);
    get_parameters(blk_method2, parmsStep2, prhs[3]);
    
    //dump_parameters(parmsStep1, 1);
    //dump_parameters(parmsStep2, 2);
   
   
    //alocate Output images
    OUT_IMG = mxCreateDoubleMatrix(height, width, mxREAL);
    double *fpMatOutput = (double*)mxGetPr(OUT_IMG);  
    
    OUT_BASIC_IMG = mxCreateDoubleMatrix(height, width, mxREAL);
    double *fpMatOutput1 = (double*)mxGetPr(OUT_BASIC_IMG);     
    
	unsigned       wh           = (unsigned) width * height;
	unsigned       whc          = (unsigned) wh * chnls;
	img_noisy.resize(whc);
    img_filtered.resize(whc);

    //mexPrintf("SIZE    :<%d><%d>\n", height, width);

    //! Load image
    for(int n=0 ; n < width ; n++) {
      for(int m=0 ; m < height ; m++) {
        img_noisy[m*width+n]=(float)fpMatInput[m*width+n];
      }      
    }
    
    //! Load image
    for(int n=0 ; n < width ; n++) {
      for(int m=0 ; m < height ; m++) {
        img_filtered[m*width+n]=(float)fpFilteredMatInput[m*width+n];
      }      
    }
    
    parmsStep1.filtered_img = img_filtered;
    parmsStep2.filtered_img = img_filtered;

    //! Denoising
    run_sd_bm3d(fSigma, 
                img_noisy, 
                img_basic, 
                img_denoised, 
                width, 
                height, 
                chnls,
                blk_method1,
                parmsStep1,
                blk_method2,
                parmsStep2);
    
      //Copy result in matlab output var
    for(int n=0 ; n < width ; n++) {
      for(int m=0 ; m < height ; m++) {  
        fpMatOutput[m*width+n]=(double)img_denoised[m*width+n];
        fpMatOutput1[m*width+n]=(double)img_basic[m*width+n];
      }      
    }  

}


