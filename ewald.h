#ifndef EWALD_HINCLUDED
#define EWALD_HINCLUDED

#include "pkd.h"

void pkdEwaldInit(PKD,double,int);
int pkdBucketEwald(PKD,int,int,double,int);

#endif
