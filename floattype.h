#ifndef FLOATTYPE_INCLUDED
#define FLOATTYPE_INCLUDED

#include <limits.h>
#include <float.h>

#ifndef FLT_MAX
#define FLT_MAX 3.402823466E+38F
#endif

#ifndef DBL_MAX
#define DBL_MAX 1.7976931348623157E+308
#endif

#define FLOAT double
#define FLOAT_MAXVAL DBL_MAX

#endif
