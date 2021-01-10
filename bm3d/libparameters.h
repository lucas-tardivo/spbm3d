#ifndef LIB_PARAMETERS_INCLUDED
#define LIB_PARAMETERS_INCLUDED


#include "lib_block_matching.h"

void get_parameters(block_matching_ftype &blk_method, sdbm3dPtype &parms, const mxArray *pm);

void dump_parameters(sdbm3dPtype parms, int step);

#endif // LIB_PARAMETERS_INCLUDED