
#ifndef _LIB_MLE_H_
#define _LIB_MLE_H_

#include <algorithm>    // std::sort
#include <vector>
#include "mex.h"
#include <math.h>
#include <string.h>
#include "boost/math/distributions/chi_squared.hpp"


typedef float (*mle_estimator)(const std::vector<float> &img, int k_r, int kHW, int width, int row, const std::vector<float> filtered_img);

float he_greenshields_sigma_mle(const std::vector<float> &img, int k_r, int kHW, int width, int row, const std::vector<float> filtered_img);


mle_estimator get_mle_function(char * pStr);
#endif


