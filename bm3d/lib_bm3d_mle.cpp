/*
 * Copyright (c) 2009-2011, A. Buades <toni.buades@uib.es>,
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "mex.h"
#include <math.h>
#include <string.h>
#include <vector>
#include "libdistances_bm3d.h"
#include "lib_bm3d_mle.h"

using namespace std;

float he_greenshields_sigma_mle(const std::vector<float> &img, int k_r, int kHW, int width, int row, const std::vector<float> filtered_img) {
    float f = filtered_img[k_r];
  
    return f;
}

mle_estimator get_mle_function(char * pStr) {
    if (strcmp(pStr,"he_greenshields_sigma_mle")==0) {
        return &he_greenshields_sigma_mle;       
    } else {
         mexPrintf("<%s>\n", pStr);
         mexErrMsgTxt("Invalid mle function.");
    }    
}