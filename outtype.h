#ifndef OUTTYPE_HINCLUDED
#define OUTTYPE_HINCLUDED

#include "pkd.h"

#define OUT_COLOR_ARRAY		1
#define OUT_DENSITY_ARRAY	2
#define OUT_POT_ARRAY		3
#define OUT_AMAG_ARRAY		4
#define OUT_IMASS_ARRAY		5
#define OUT_RUNG_ARRAY		6

#ifdef GASOLINE

#define OUT_U_ARRAY			7
#define OUT_UDOT_ARRAY		8
#define OUT_HSMDIVV_ARRAY	9

#endif

#ifdef PLANETS

#define OUT_DT_ARRAY		10
#define OUT_CT_ARRAY		11

#endif /* PLANETS */

#define OUT_POS_VECTOR		1
#define OUT_VEL_VECTOR		2
#define OUT_ACCEL_VECTOR	3

#ifdef GASOLINE

#define OUT_VPRED_VECTOR	4

#endif

void pkdOutArray(PKD,char *,int);
void pkdOutVector(PKD,char *,int,int);

#endif







