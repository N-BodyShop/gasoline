#ifndef OUTTYPE_HINCLUDED
#define OUTTYPE_HINCLUDED

#include "pkd.h"

#define OUT_COLOR_ARRAY		1
#define OUT_DENSITY_ARRAY	2
#define OUT_POT_ARRAY		3
#define OUT_AMAG_ARRAY		4
#define OUT_IMASS_ARRAY		5
#define OUT_RUNG_ARRAY		6

#define OUT_POS_VECTOR		1
#define OUT_VEL_VECTOR		2
#define OUT_ACCEL_VECTOR	3

void pkdOutArray(PKD,char *,int);
void pkdOutVector(PKD,char *,int,int);

#endif
