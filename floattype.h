#ifndef FLOATTYPE_HINCLUDED
#define FLOATTYPE_HINCLUDED

#include <limits.h>
#include <float.h>

#ifndef FLT_MAX
#define FLT_MAX 3.402823466E+38F
#endif

#ifndef DBL_MAX
#define DBL_MAX 1.7976931348623157E+308
#endif

/*
 * If you change these, see also xdr_FLOAT() in outtype.c
 */
#ifndef SINGLE
#define FLOAT double
#define FLOAT_MAXVAL DBL_MAX
#else
#define FLOAT float
#define FLOAT_MAXVAL FLT_MAX
#endif

#endif
