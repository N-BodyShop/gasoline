#ifndef FLOATTYPE_INCLUDED
#define FLOATTYPE_INCLUDED

#include <limits.h>
#include <values.h>

#ifndef SINGLE

#define FLOAT			double
#define FLOAT_MAXVAL	DBL_MAX

#else

#define FLOAT			float
#define FLOAT_MAXVAL	FLT_MAX

#endif

#endif
