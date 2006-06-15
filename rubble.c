#ifdef RUBBLE_ZML
/* NOTE: some of the rpu routines contain errors remember to update rpu routines
  once they have been fixed */

/* ZML 02.19.03 resolves collisions in a planetesimal disk */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h> /*for getpid()*/
#include <math.h>
#include <assert.h>
#include "rubble.h"
#include "collision.h"
#include "random.h"

#define SMALL_MASS 0.05

#ifndef PI
#define PI M_PI
#endif

#ifndef TWO_PI
#define TWO_PI (2*PI)
#endif

/* parameters needed for rubble pile creation */
#define PART_SCALING 0.99
#define PART_MAX_SPEED 0.1
#define PEMAX 0.74 /* theoretical maximum packing efficiency */
#define PEN(n) ((n) == 1 ? 1.0 : 2*PEMAX/PI*atan(0.496*pow((n),0.264)))
#define PER(r) ((r) == 1 ? 1.0 : 2*PEMAX/PI*atan(0.368*pow((r),0.843)))

/* put comment here */

#define N_MU 4
#define N_B 5
#define N_VEL 5
#define N_EPS 5

typedef struct {
	double b;
	double v;
	double mu;
	double eps;
	} INTERP_PARAMS;

/*data compiled using find_data - from simsII and sims dirs*/
/*b=1.0 and v=0.0 data was "pinned" by hand*/
/*interpolation grid data from collision simulations done 09.02*/
/*the mass of the largest remnant is stored as a preinitialized*/
/*4x4x4x5 constant double array to avoid having individual processors*/ 
/*access a file*/

double Y_prim[N_MU][N_B][N_VEL][N_EPS] = {
	{{{1.000000, 1.000000, 1.000000, 1.000000, 1.000000},
		  {1.000000, 1.000000, 1.000000, 0.998953, 0.999476},
		  {0.991100, 0.988482, 0.981675, 0.966492, 0.900524},
		  {0.900000, 0.891100, 0.878011, 0.830367, 0.624084},
		  {0.685864, 0.704189, 0.669630, 0.564398, 0.005759}},
		 {{1.000000, 1.000000, 1.000000, 1.000000, 1.000000},
			  {1.000000, 1.000000, 1.000000, 1.000000, 1.000000},
			  {0.989005, 0.986911, 0.983246, 0.969634, 0.908901},
			  {0.506283, 0.884817, 0.876963, 0.808901, 0.314136},
			  {0.274869, 0.257068, 0.256021, 0.230890, 0.011518}},
		 {{1.000000, 1.000000, 1.000000, 1.000000, 1.000000},
			  {1.000000, 1.000000, 1.000000, 1.000000, 1.000000},
			  {0.495288, 0.523560, 0.501571, 0.482199, 0.475393},
			  {0.450262, 0.451833, 0.434032, 0.410995, 0.362827},
			  {0.380628, 0.375393, 0.376440, 0.343456, 0.271204}},
		 {{1.000000, 1.000000, 1.000000, 1.000000, 1.000000},
			  {1.000000, 1.000000, 1.000000, 1.000000, 1.000000},
			  {0.490576, 0.491623, 0.495288, 0.493717, 0.487434},
			  {0.481675, 0.484817, 0.481152, 0.468587, 0.460733},
			  {0.472251, 0.470157, 0.467016, 0.461780, 0.426178}},
		 {{1.000000, 1.000000, 1.000000, 1.000000, 1.000000},
			  {0.500000, 0.500000, 0.500000, 0.500000, 0.500000},
			  {0.500000, 0.500000, 0.500000, 0.500000, 0.500000},
			  {0.500000, 0.500000, 0.500000, 0.500000, 0.500000},
			  {0.500000, 0.500000, 0.500000, 0.500000, 0.500000}}},
	{{{1.000000, 1.000000, 1.000000, 1.000000, 1.000000},
		  {1.000000, 1.000000, 0.999215, 0.998429, 0.992930},
		  {0.977219, 0.967793, 0.966221, 0.940299, 0.881383},
		  {0.898665, 0.897093, 0.877455, 0.824823, 0.684211},
		  {0.756481, 0.741555, 0.706991, 0.619796, 0.427337}},
		 {{1.000000, 1.000000, 1.000000, 1.000000, 1.000000},
			  {1.000000, 1.000000, 0.999215, 0.998429, 0.992930},
			  {0.962294, 0.959152, 0.943441, 0.906520, 0.834250},
			  {0.781618, 0.796544, 0.771406, 0.728987, 0.637863},
			  {0.646504, 0.676355, 0.668500, 0.589945, 0.443833}},
		 {{1.000000, 1.000000, 1.000000, 1.000000, 1.000000},
			  {1.000000, 1.000000, 1.000000, 1.000000, 0.996072},
			  {0.790259, 0.794972, 0.785546, 0.768264, 0.748625},
			  {0.725059, 0.722702, 0.718775, 0.717989, 0.664572},
			  {0.674784, 0.668500, 0.666928, 0.663001, 0.593873}},
		 {{1.000000, 1.000000, 1.000000, 1.000000, 1.000000},
			  {1.000000, 1.000000, 1.000000, 0.999215, 1.000000},
			  {0.745483, 0.755695, 0.756481, 0.757266, 0.749411},
			  {0.736057, 0.737628, 0.734485, 0.732129, 0.722702},
			  {0.712490, 0.716418, 0.712490, 0.709348, 0.685782}},
		 {{1.000000, 1.000000, 1.000000, 1.000000, 1.000000},
			  {0.666667,0.666667,0.666667,0.666667,0.666667},
			  {0.666667,0.666667,0.666667,0.666667,0.666667},
			  {0.666667,0.666667,0.666667,0.666667,0.666667},
			  {0.666667,0.666667,0.666667,0.666667,0.666667}}},
	{{{1.000000, 1.000000, 1.000000, 1.000000, 1.000000},
		  {0.999102, 0.999102, 0.999102, 0.997307, 0.989228},
		  {0.973968, 0.978456, 0.964093, 0.942549, 0.878815},
		  {0.901257, 0.908438, 0.903950, 0.833034, 0.746858},
		  {0.780969, 0.774686, 0.734291, 0.681329, 0.575404}},
		 {{1.000000, 1.000000, 1.000000, 1.000000, 1.000000},
			  {1.000000, 0.997307, 0.997307, 0.996409, 0.983842},
			  {0.936266, 0.946140, 0.939856, 0.914722, 0.873429},
			  {0.843806, 0.841113, 0.822262, 0.806104, 0.738779},
			  {0.761221, 0.732495, 0.739677, 0.720826, 0.593357}},
		 {{1.000000, 1.000000, 1.000000, 1.000000, 1.000000},
			  {0.998205, 0.994614, 0.991023, 0.989228, 0.968582},
			  {0.880610, 0.884201, 0.881508, 0.874327, 0.862657},
			  {0.834830, 0.820467, 0.833034, 0.806104, 0.756732},
			  {0.768402, 0.777379, 0.776481, 0.741472, 0.673249}},
		 {{1.000000, 1.000000, 1.000000, 1.000000, 1.000000},
			  {0.880610, 0.881508, 0.883303, 0.886894, 0.885996},
			  {0.856373, 0.854578, 0.858169, 0.854578, 0.852783},
			  {0.850090, 0.847397, 0.851885, 0.849192, 0.836625},
			  {0.843806, 0.842011, 0.834830, 0.826751, 0.809695}},
		 {{1.000000, 1.000000, 1.000000, 1.000000, 1.000000},
			  {0.833333,0.833333,0.833333,0.833333,0.833333},
			  {0.833333,0.833333,0.833333,0.833333,0.833333},
			  {0.833333,0.833333,0.833333,0.833333,0.833333},
			  {0.833333,0.833333,0.833333,0.833333,0.833333}}},
	{{{1.000000, 1.000000, 1.000000, 1.000000, 1.000000},
		  {0.999058, 1.000000, 0.999058, 0.994345, 0.983035},
		  {0.967012, 0.967012, 0.961357, 0.935909, 0.885014},
		  {0.893497, 0.894439, 0.869934, 0.846371, 0.773798},
		  {0.810556, 0.824694, 0.802073, 0.739868, 0.638077}},
		 {{1.000000, 1.000000, 1.000000, 1.000000, 1.000000},
			  {0.992460, 0.995287, 0.995287, 0.990575, 0.980207},
			  {0.950990, 0.944392, 0.940622, 0.920829, 0.873704},
			  {0.897267, 0.879359, 0.865221, 0.839774, 0.771913},
			  {0.808671, 0.802073, 0.800188, 0.773798, 0.672950}},
		 {{1.000000, 1.000000, 1.000000, 1.000000, 1.000000},
			  {0.964185, 0.960415, 0.967955, 0.956645, 0.947220},
			  {0.907634, 0.906692, 0.901979, 0.899152, 0.881244},
			  {0.869934, 0.876532, 0.872761, 0.853911, 0.806786},
			  {0.840716, 0.836946, 0.835061, 0.807728, 0.769086}},
		 {{1.000000, 1.000000, 1.000000, 1.000000, 1.000000},
			  {0.911404, 0.912347, 0.913289, 0.914232, 0.914232},
			  {0.901037, 0.902922, 0.901979, 0.898209, 0.892554},
			  {0.898209, 0.898209, 0.893497, 0.896324, 0.881244},
			  {0.890669, 0.887842, 0.889727, 0.888784, 0.875589}},
		 {{1.000000, 1.000000, 1.000000, 1.000000, 1.000000},
			  {0.888889,0.888889,0.888889,0.888889,0.888889},
			  {0.888889,0.888889,0.888889,0.888889,0.888889},
			  {0.888889,0.888889,0.888889,0.888889,0.888889},
			  {0.888889,0.888889,0.888889,0.888889,0.888889}}} 
	};

double Y_sec[N_MU][N_B][N_VEL][N_EPS] = {
	{{{0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
		  {0.000000, 0.000000, 0.000000, 0.000524, 0.000524},
		  {0.000524, 0.000524, 0.001047, 0.001047, 0.001047},
		  {0.002094, 0.003665, 0.001571, 0.002094, 0.002618},
		  {0.010471, 0.006806, 0.010471, 0.014660, 0.001571}},
		 {{0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
			  {0.000000, 0.000000, 0.000000, 0.000000, 0.000524},
			  {0.001047, 0.001047, 0.000524, 0.000524, 0.001047},
			  {0.378534, 0.002618, 0.006283, 0.004712, 0.203141},
			  {0.169634, 0.224607, 0.172775, 0.104188, 0.001571}},
		 {{0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
			  {0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
			  {0.486387, 0.469634, 0.483246, 0.468063, 0.448168},
			  {0.420419, 0.402094, 0.406283, 0.395288, 0.343456},
			  {0.363874, 0.346073, 0.362827, 0.301047, 0.240314}},
		 {{0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
			  {0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
			  {0.486387, 0.487958, 0.485340, 0.483770, 0.470681},
			  {0.471728, 0.468063, 0.471204, 0.457592, 0.452356},
			  {0.463351, 0.457068, 0.453927, 0.441885, 0.418848}},
		 {{0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
			  {0.500000,0.500000,0.500000,0.500000,0.500000},
			  {0.500000,0.500000,0.500000,0.500000,0.500000},
			  {0.500000,0.500000,0.500000,0.500000,0.500000},
			  {0.500000,0.500000,0.500000,0.500000,0.500000}}},
	{{{0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
		  {0.000000, 0.000000, 0.000786, 0.000786, 0.000786},
		  {0.000786, 0.001571, 0.001571, 0.001571, 0.000786},
		  {0.001571, 0.001571, 0.002357, 0.003142, 0.003142},
		  {0.010998, 0.010212, 0.006284, 0.005499, 0.003928}},
		 {{0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
			  {0.000000, 0.000000, 0.000786, 0.000786, 0.000786},
			  {0.002357, 0.002357, 0.004713, 0.010998, 0.008641},
			  {0.038492, 0.029065, 0.019639, 0.011783, 0.002357},
			  {0.028280, 0.008641, 0.007855, 0.017282, 0.009427}},
		 {{0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
			  {0.000000, 0.000000, 0.000000, 0.000786, 0.000786},
			  {0.131186, 0.113904, 0.138256, 0.113119, 0.077769},
			  {0.086410, 0.085625, 0.072270, 0.034564, 0.003928},
			  {0.022781, 0.020424, 0.029851, 0.014140, 0.002357}},
		 {{0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
			  {0.000000, 0.042419, 0.000000, 0.039277, 0.123331},
			  {0.219167, 0.215240, 0.218382, 0.208955, 0.199529},
			  {0.193244, 0.202671, 0.174391, 0.164179, 0.139827},
			  {0.161037, 0.162608, 0.162608, 0.124116, 0.025137}},
		 {{0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
			  {0.333333,0.333333,0.333333,0.333333,0.333333},
			  {0.333333,0.333333,0.333333,0.333333,0.333333},
			  {0.333333,0.333333,0.333333,0.333333,0.333333},
			  {0.333333,0.333333,0.333333,0.333333,0.333333}}},
	{{{0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
		  {0.000898, 0.000898, 0.000898, 0.000898, 0.000898},
		  {0.000898, 0.001795, 0.000898, 0.000898, 0.001795},
		  {0.003591, 0.001795, 0.001795, 0.002693, 0.001795},
		  {0.005386, 0.006284, 0.005386, 0.002693, 0.002693}},
		 {{0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
			  {0.000898, 0.000898, 0.000898, 0.000898, 0.000898},
			  {0.008977, 0.003591, 0.002693, 0.000898, 0.000898},
			  {0.007181, 0.006284, 0.009874, 0.004488, 0.001795},
			  {0.017056, 0.005386, 0.005386, 0.004488, 0.002693}},
		 {{0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
			  {0.001795, 0.002693, 0.002693, 0.003591, 0.006284},
			  {0.017953, 0.013465, 0.008977, 0.008977, 0.001795},
			  {0.007181, 0.004488, 0.003591, 0.002693, 0.001795},
			  {0.004488, 0.002693, 0.002693, 0.001795, 0.001795}},
		 {{0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
			  {0.118492, 0.115799, 0.114901, 0.110413, 0.096050},
			  {0.111311, 0.102334, 0.109515, 0.074506, 0.007181},
			  {0.061939, 0.059246, 0.054758, 0.046679, 0.001795},
			  {0.049372, 0.054758, 0.031418, 0.017056, 0.001795}},
		 {{0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
			  {0.166667,0.166667,0.166667,0.166667,0.166667},
			  {0.166667,0.166667,0.166667,0.166667,0.166667},
			  {0.166667,0.166667,0.166667,0.166667,0.166667},
			  {0.166667,0.166667,0.166667,0.166667,0.166667}}},
	{{{0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
		  {0.000943, 0.000943, 0.000943, 0.000943, 0.000943},
		  {0.000943, 0.000943, 0.000943, 0.000943, 0.001885},
		  {0.003770, 0.002828, 0.001885, 0.001885, 0.001885},
		  {0.004713, 0.004713, 0.003770, 0.005655, 0.001885}},
		 {{0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
			  {0.000943, 0.000943, 0.000943, 0.000943, 0.000943},
			  {0.002828, 0.001885, 0.000943, 0.001885, 0.001885},
			  {0.004713, 0.002828, 0.003770, 0.001885, 0.001885},
			  {0.006598, 0.003770, 0.002828, 0.002828, 0.000943}},
		 {{0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
			  {0.006598, 0.014138, 0.012253, 0.013195, 0.000943},
			  {0.005655, 0.008483, 0.001885, 0.001885, 0.000943},
			  {0.002828, 0.001885, 0.001885, 0.001885, 0.000943},
			  {0.001885, 0.002828, 0.002828, 0.000943, 0.000943}},
		 {{0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
			  {0.082941, 0.080113, 0.079171, 0.065975, 0.065033},
			  {0.085768, 0.084826, 0.077286, 0.056550, 0.009425},
			  {0.053723, 0.037700, 0.032045, 0.003770, 0.001885},
			  {0.019793, 0.001885, 0.002828, 0.001885, 0.000943}},
		 {{0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
			  {0.111111,0.111111,0.111111,0.111111,0.111111},
			  {0.111111,0.111111,0.111111,0.111111,0.111111},
			  {0.111111,0.111111,0.111111,0.111111,0.111111},
			  {0.111111,0.111111,0.111111,0.111111,0.111111}}}
	};

/* numerical recipe routines */

#define N_SUBDIM 2

void 
polin2(double x1a[],double x2a[],double x3a[],double x4a[],
			double ya[N_SUBDIM][N_SUBDIM][N_SUBDIM][N_SUBDIM],int l,int m,
			int n,int p,double x1,double x2,double x3,double x4,double *y,
			double *dy);
void 
polint(double xa[],double ya[],int n,double x,double *y,double *dy);
static double 
ran(void);

/* pkd routines */
void
rubLocateBin(PARTICLE *p,DUST_BINS_PARAMS *DB,double M,double *a,int *piBin)
{
	/*calculates semi-major axis and finds the bin that the particle belongs to*/

	double r,r_xy,v2;

	r = sqrt(SQ(p->r[0]) + SQ(p->r[1]) + SQ(p->r[2]));
	r_xy = sqrt(SQ(p->r[0]) + SQ(p->r[1]));
	v2 = SQ(p->v[0]) + SQ(p->v[1]) + SQ(p->v[2]);
	assert(r > 0.0);
	assert(r_xy > 0.0);
	assert(M > 0.0);
	/* semimajor axis */
	*a = 1.0/(2.0/r - v2/(M + p->fMass));
/*	printf("particle %d has semimajor axis = %f\n",p->iOrder, *a);*/
	if (*a <= 0.0 && p->iColor != PLANETESIMAL) {/*DEBUG*/
		*piBin = -1;
		printf("rubble particle %d has neg a = %f and will be removed\n",
			   p->iOrder, *a);
		printf("v2 = %e, r = %f, mass = %e\n", v2,r,p->fMass);
		printf("iColor = %i, position = %e %e %e\n", 
			   p->iColor, p->r[X], p->r[Y], p->r[Z]);
		printf("velocity = %e %e %e\n", p->v[X], p->v[Y], p->v[Z]);
		return;
	} else 	/*allow of rubble particles to be kicked out*/
		assert(*a > 0.0); /* orbits of non-rubble must be bound, with r > 0.0 */

	if (DB->nDustBins == 0) {
		*piBin = -1;
		return;
		}
	assert(DB->nDustBins > 0);
	/*
	 ** Now figure out what dust bin the planetesimal is in.
	 ** NOTE: for simplicity, we base this solely on semimajor axis.
	 ** To do this properly, we would want to take the planetesimal's
	 ** eccentricity into account, in case the planetesimal's orbit
	 ** takes it through more than one bin. Implicitly, we're assuming
	 ** the eccentricity is small, and the accretion interval is small
	 ** compared to the orbit period.
	 */
	/* 01.28.04 - I think that iBin should be determined by instantaneous
	 ** xy position not the semi-major axis of the planetesimal 
	 */

/*	*piBin = (*a - DB->dDustBinsInner)/DB->dDustBinsWidth;*/
	*piBin = (r_xy - DB->dDustBinsInner)/DB->dDustBinsWidth;
/*	printf("*a = %f, Inner = %f, Outer = %f\n",
		   *a,DB->dDustBinsInner, DB->dDustBinsWidth);
	printf("particle = %d rubLocateBin = %d\n",p->iOrder,*piBin);*/
/*	if (*a == DB->dDustBinsOuter) *piBin = DB->nDustBins - 1;*/ /* very unlikely! */
	if (r_xy == DB->dDustBinsOuter) *piBin = DB->nDustBins - 1; /* very unlikely! */
	}

void 
pkdDustBinsGetMass(PKD pkd,DUST_BINS_PARAMS *DB,DustBins pDustBins[],
				   double dt,double M,double pDustBinsMassLoss[])
{
	PARTICLE *p;
	double a;
	int i,nLocal = pkdLocal(pkd);

	for (i=0;i<nLocal;i++) {

		p = &pkd->pStore[i];

		/* first compute orbital parameters */

		rubLocateBin(p,DB,M,&a,&p->iBin);

		if (p->iBin < 0 || p->iBin >= DB->nDustBins)
			p->dDustMass = 0.0;
		else {
			double omega,P,h[3],h2,e,rho;
			if (a < 0.0) /*DEBUG*/
				printf("bin = %i, particle = %i, color = %i\n", 
					   p->iBin, p->iOrder, p->iColor); 
			assert(a > 0.0);
			/* angular speed */
			omega = sqrt((M + p->fMass)/CUBE(a));
			/* orbital period */
			P = TWO_PI/omega;
			assert(P > 0.0);
			/* eccentricity */
			h[0] = p->r[1]*p->v[2] - p->r[2]*p->v[1];
			h[1] = p->r[2]*p->v[0] - p->r[0]*p->v[2];
			h[2] = p->r[0]*p->v[1] - p->r[1]*p->v[0];
			h2 = SQ(h[0]) + SQ(h[1]) + SQ(h[2]);
			e = sqrt(1.0 - h2/(a*(M + p->fMass)));
			assert(e < 1.0);

			/* compute mass increase, ignoring particles outside dust region */

			assert(pDustBins[p->iBin].dVolume > 0.0);
			rho = pDustBins[p->iBin].dMass/pDustBins[p->iBin].dVolume;
			p->dDustMass = e*PI*SQ(RADIUS(p))*TWO_PI*a*rho*dt/P;
			
			/* accumulate mass loss */

			pDustBinsMassLoss[p->iBin] += p->dDustMass;
/*			if (p->dDustMass < 1e-19)*/ /*DEBUG*/
/*				printf("iBin = %i, dDustMass = %e, e=%f, h2=%e, M=%f, a=%f, diff=%e\n", 
					   p->iBin, p->dDustMass, e, h2, M, a, (1.0 - h2/(a*(M + p->fMass))));*/
			}
		}
/*	for (i=0;i<10;i++) 
		printf("iBin = %i, Mass loss = %e\n", i, pDustBinsMassLoss[p->iBin]);*/
	}

void 
pkdDustBinsApply(PKD pkd,double M,double pMassIncrFrac[],int nBins)
{
	/*
	 ** Applies dust accretion to particles, increasing mass
	 ** (by a maximum amount determined by the mass left in the
	 ** dust bin) and radius (assuming constant density), and
	 ** adjusting the velocity (in the manner of a momentum-
	 ** conserving drag force).
	 */

	PARTICLE *p;
	double dm,R,omega,vkx,vky;
	int i;

/*	return;*/ /*DEBUG*/
/*	printf("GOt to far/n");*/

	for (i=0;i<pkd->nLocal;i++) {
		p = &pkd->pStore[i];

		if (p->iBin >= 0 && p->iBin < nBins) {

			dm = pMassIncrFrac[p->iBin]*p->dDustMass;
/*			printf("dm = %e\n", dm);*/
			if (dm <= 0.0) {
				printf("dm = %e, p->dDustMass = %e, Incr = %e, iBin = %i\n", dm, p->dDustMass, 
					   pMassIncrFrac[p->iBin], p->iBin);
				}
			assert(dm >= 0.0);
			assert(p->fMass > dm); /* if it isn't, something is probably wrong */

		/* compute change to orbit, based on momentum conservation */

			omega = sqrt(M + p->fMass)/pow(SQ(p->r[0]) + SQ(p->r[1]),0.75); /* ang speed in plane */
			vkx = - omega*p->r[1]; /* Kepler components */
			vky =   omega*p->r[0];

		/*
		 ** If planetesimal's instantaneous velocity is larger than
		 ** Keplerian at that location the planetesimal should suffer
		 ** a drag. If the planetesimal's instantaneous velocity
		 ** is smaller than Keplerian the planetesimal should 
		 ** receive a kick. The dust should circularize the planetesimal's
		 ** orbit. 12.4.03
		 */
		
			p->v[0] = vkx + p->fMass*(p->v[0] - vkx)/(dm + p->fMass);
			p->v[1] = vky + p->fMass*(p->v[1] - vky)/(dm + p->fMass);
			p->v[2] = p->fMass*p->v[2]/(dm + p->fMass);
			/* grow radius and mass of particle */

			R = RADIUS(p)*pow(1.0 + dm/p->fMass,1.0/3.0);
			p->fSoft = 0.5*R; /* radius is twice particle softening */
			p->fMass += dm;

		} 
		else {
/*			double temp_pos, temp_speed;*/
/*			printf("WARNING: planetesimal %d is outside of bin range %d\n",
				   p->iOrder, p->iBin);*/ /*DEBUG*/
/*			printf("iBin = %d\n", p->iBin);
			printf("x = %f y = %f z = %f, mag_xy = %f\n", 
				   p->r[0], p->r[1], p->r[2], 
				   sqrt(SQ(p->r[0]) + SQ(p->r[1])));
			printf("iColor = %i, position = %e %e %e\n", p->iColor, 
				   p->r[X], p->r[Y], p->r[Z]);
			printf("velocity = %e %e %e, mass = %e\n", p->v[X],
				   p->v[Y], p->v[Z], p->fMass);
			temp_pos = sqrt(SQ(p->r[0])+SQ(p->r[1])+SQ(p->r[2]));
			temp_speed = SQ(p->v[0])+SQ(p->v[1])+SQ(p->v[2]);
			printf("a = %f\n", 
				   1.0/(2.0/temp_pos - temp_speed/(M+p->fMass)));*/
			} /* DEBUG 01.28.04 */
		}
	}

void 
pkdRubbleResetColFlag(PKD pkd)
{
	int i,nLocal = pkdLocal(pkd);

	for (i=0;i<nLocal;i++)
		pkd->pStore[i].bMayCollide = 0;
	}

int 
pkdRubbleCheckForKDKRestart(PKD pkd)
{
	int i,nLocal = pkdLocal(pkd);

	for (i=0;i<nLocal;i++)
		if (pkd->pStore[i].bMayCollide)
			return 1;

	return 0;
	}

void 
pkdRubbleStep(PKD pkd,double dMaxStep,double dMinStep)
{
	/*
	 ** We want rubble particles to be on the bottom rung and "normal"
	 ** (planetesimal) particles to be on the top rung. Note that a
	 ** particle is considered a "planetesimal" if its color is equal
	 ** to PLANETESIMAL (defined in ssdefs.h; usually 3, i.e. green).
	 ** Later we may want to allow a range of rungs; for now we stick
	 ** with a bimodal distribution.
	 **
	 ** ADDED: to keep rungs in sync (since collisions are determined
	 ** *during* the drift), planetesimals that *may* collide during
	 ** this interval should also be forced to the lowest rung.
	 */

	PARTICLE *p;
	int i,nLocal = pkdLocal(pkd);
	/*DEBUG*/
/*	printf("Got to pkdRubbleStep \n");*/

	for (i=0;i<nLocal;i++) {
		p = &pkd->pStore[i];
/*printf("%d bMaycollide = %d\n", p->iColor,p->bMayCollide); */
		if (p->iColor == PLANETESIMAL && !(p->bMayCollide)) 
			pkd->pStore[i].dt = dMaxStep;
		else
			pkd->pStore[i].dt = dMinStep;
/*		if (pkd->pStore[i].iOrder == 2536 || pkd->pStore[i].iOrder == 1834)
			printf("particle %d is on step %e Minstep = %e Maxstep = %e\n",
				   pkd->pStore[i].iOrder, pkd->pStore[i].dt,dMinStep,
				   dMaxStep);*/ /* DEBUG 01.30.04*/
		}
	}

void 
pkdRubCleanup(PKD pkd,int iColor,DUST_BINS_PARAMS *DB,double M,
				   double pDustBinsMassGain[],double *dDustBinsRubTrash)
{
	PARTICLE *p;
	double a;
	int i,iBin,nLocal = pkdLocal(pkd);

/*	printf("Got to RubCleanup - returning\n");
	return; *//*DEBUG*/
/*	pritnf("Got too far\n");*/

	for (i=0;i<nLocal;i++) {
		p = &pkd->pStore[i];
/*		assert(p->iOrder >= 0);*/ /*DEBUG 06.01.04*/
		if (p->iColor == iColor) {
/*			if (p->fMass >= RUB_MIN_MASS)*/
			if (p->fMass >= DB->dRubMinMass) {
				printf("original color %d\n", p->iColor);
				p->iColor = PLANETESIMAL;
/*				printf("particle %d,mass %e, color %d, rung %d is kept\n", 
					   p->iOrder, p->fMass, p->iColor, p->iRung);*/ /*DEBUG*/
				{ /*DEBUG -- CHECK FOR ABNORMAL DENSITY*/
				double rho;
				rho = p->fMass/(4.0/3.0*M_PI*CUBE(RADIUS(p)));
				/*if (rho > 2e4) for f=6*/
				if (rho > 4e6) /*for f=1*/
					(void) printf("DUST CHECK BAD DENSITY: iOrder=%i iOrgIdx=%i rho=%g\n",p->iOrder,p->iOrgIdx,rho);
				}
				}
			else {
/*				printf("particle %d, mass %e, color %d, rung %d went to dust\n", 
					   p->iOrder, p->fMass, p->iColor, p->iRung);*/ /*DEBUG*/
/*				printf("Rubble Cleanup\n");*/
				rubLocateBin(p,DB,M,&a,&iBin); /* "a" ignored */
				/* dust outside bin range goes into overflow variable */
				if (iBin >= 0 && iBin < DB->nDustBins && a > 0.0)
					pDustBinsMassGain[iBin] += p->fMass;
				else {
					*dDustBinsRubTrash += p->fMass;
/*					printf("particle %d, is outside dust bin range %f, mass %e has gone to trash\n",
						   p->iOrder, a, p->fMass);*/ /*DEBUG*/
					}
				pkdDeleteParticle(pkd,p);	
/*				printf("particle %d, iActive %d, iColor %d, rung %d\n",
					   p->iOrder, p->iActive, p->iColor, p->iRung); DEBUG */
				}	
			}
		}
	}

void 
pkdRubInterpCleanup(PKD pkd,DUST_BINS_PARAMS *DB,double M,int iOrder,
					int *iBin,double *dDustBinsInterpMass)
{
	PARTICLE *p;
	double a;
	int i,nLocal = pkdLocal(pkd);

	*iBin = INT_MAX; /* initialization */
	*dDustBinsInterpMass = 0.0;
/*	printf("Got to RubInterpCleanup - returning\n");
	return; *//*DEBUG*/
/*	printf("Got too far\n");*/

/*	printf("iOrder = %i, dDustBinsInterpMass = %e\n",iOrder,*dDustBinsInterpMass);
	printf("Looking through %i particles\n", pkdLocal(pkd));*/ /*DEBUG*/
	for (i=0;i<nLocal;i++) {
		p = &pkd->pStore[i];
		if (p->iOrder == iOrder) {
/*			printf("Interpolated Collision Cleanup\n");*/ /*DEBUG*/
			rubLocateBin(p,DB,M,&a,iBin); /* "a" and iBin ignored */
			assert(a > 0.0); /*paranoia check*/
			/* accept dust outside bin range -- let master handle it */
		/*	printf("p->dDustMass = %e\n",p->dDustMass);*/
			*dDustBinsInterpMass = p->dDustMass;
/*			printf("p->dDustMass = %e, dDustBinsInterpMass = %e\n",
				   p->dDustMass,*dDustBinsInterpMass);*/
/*			printf("particle %d, position: %f %f %f\n",
				   iOrder, p->r[0], p->r[1], p->r[2]);
			printf("velocity: %e %e %e\n", p->v[0],p->v[1],p->v[2]);*/
			return; /* no need to check other particles */
			}
		}
/*	printf("Didn't find particle with iOrder = %d\n",iOrder);*/
	}


/* functions below from vector.c dcr 94-08-24 */	
/* zml added to rubble.c 03-05-02 */

static void 
swap(double *x, double *y)
{
	/* Swaps double values pointed to by "x" and "y" */
  
	double t = *x;
  
	*x = *y;
	*y = t;
	}

void 
Transpose(MATRIX a)
{
	/* Applies transpose operator on matrix "a" */
  
	swap(&a[X][Y], &a[Y][X]);
	swap(&a[X][Z], &a[Z][X]);
	swap(&a[Y][Z], &a[Z][Y]);
	}

void 
Transform(MATRIX a, VECTOR u, VECTOR v)
{
	/* Applies matrix "a" to vector "u", returning vector "v" */

	v[X] = DOT(a[X], u);
	v[Y] = DOT(a[Y], u);
	v[Z] = DOT(a[Z], u);
	}
/* end of vector.c functions */

/* functions below from rpu.c dcr 98-09-24 */
/* zml added to rubble.c 03-05-02 */

BOOLEAN
rpuInEllipsoid(const VECTOR r0,double R,const VECTOR a)
{
	/*
	 ** Returns TRUE if a ball with radius R located at r0 lies *entirely*
	 ** within the ellipsoid defined by semi-axes "a" (measured along the
	 ** Cartesian axes and centred at the origin). To get this right, need
	 ** to compute direction from r0 to nearest point on ellipsoid surface.
	 ** This is too hard, so settle for more conservative boundary.
	 */
  
	VECTOR r;
	double d;
	int k;
  
	assert(a[X] > 0.0 && a[Y] > 0.0 && a[Z] > 0.0);
  
	/* Check easy cases first */
  
	d = 0;
  
	for (k=0;k<N_DIM;k++) {
		if (r0[k] == 0 && R > a[k]) return FALSE;
		d += SQ(r0[k]/a[k]);
		}
  
	if (d == 0) return TRUE;
	if (d > 1) return FALSE;
  
	/* Make a (conservative) guess for the general case */
  
	d = MAG(r0);
	COPY_VEC(r0,r);
	SCALE_VEC(r,1 + R/d);
	d = SQ(r[X]/a[X]) + SQ(r[Y]/a[Y]) + SQ(r[Z]/a[Z]);
	
	return (d <= 1);
	}

/*
 ** Following functions take rubble-pile arguments.
 */

void
rpuMalloc(RUBBLE_PILE *rp)
{
	assert(rp);
	assert(rp->n_particles > 0);
	rp->data = (SSDATA *) malloc(rp->n_particles*sizeof(SSDATA));
	assert(rp->data);
/*DEBUG (void) fprintf(stderr,"rpuMalloc(): %i bytes.\n",rp->n_particles*sizeof(SSDATA));*/
	}

void
rpuRealloc(RUBBLE_PILE *rp)
{
	assert(rp);
	assert(rp->n_particles > 0);
	rp->data =
		(SSDATA *) realloc((void *) rp->data,rp->n_particles*sizeof(SSDATA));
	assert(rp->data);
/*DEBUG (void) fprintf(stderr,"rpuRealloc(): %i bytes.\n",rp->n_particles*sizeof(SSDATA));*/
	}

void
rpuFree(RUBBLE_PILE *rp)
{
	assert(rp);
	assert(rp->data);
	free((void *) rp->data);
	rp->data = NULL;
/*DEBUG (void) fprintf(stderr,"rpuFree().\n");*/
	}

void rpuCalcPos(RUBBLE_PILE *rp)
{
	/*
	 ** Sets rp->pos to center-of-mass position (assuming particles are
	 ** equal mass).
	 */

	int i;

	assert(rp != NULL);
	assert(rp->n_particles > 0);
	assert(rp->data != NULL);

	ZERO_VEC(rp->pos);
	for (i=0;i<rp->n_particles;i++) {
		ADD_VEC(rp->pos,rp->data[i].pos,rp->pos);
		}
	NORM_VEC(rp->pos,rp->n_particles);
	}

void rpuCalcVel(RUBBLE_PILE *rp)
{
	/*
	 ** Sets rp->vel to center-of-mass velocity (assuming particles are
	 ** equal mass).
	 */

	int i;

	assert(rp != NULL);
	assert(rp->n_particles > 0);
	assert(rp->data != NULL);

	ZERO_VEC(rp->vel);
	for (i=0;i<rp->n_particles;i++) {
		ADD_VEC(rp->vel,rp->data[i].vel,rp->vel);
		}
	NORM_VEC(rp->vel,rp->n_particles);
	}

void
rpuTransPos(RUBBLE_PILE *rp)
{
	/* Sets new particle positions relative to rp->pos */
  
	int i;
  
	assert(rp);
	assert(rp->n_particles > 0);
	assert(rp->data);
  
	for (i=0;i<rp->n_particles;i++)
		ADD_VEC(rp->data[i].pos,rp->pos,rp->data[i].pos);
	}

void
rpuTransVel(RUBBLE_PILE *rp)
{
	/* Sets new particle velocities relative to rp->vel */
  
	int i;
  
	assert(rp);
   	assert(rp->n_particles > 0);
	assert(rp->data);
  
	for (i=0;i<rp->n_particles;i++)
		ADD_VEC(rp->data[i].vel,rp->vel,rp->data[i].vel);
	}

void rpuApplySpin(RUBBLE_PILE *rp)
{
	/*
	 ** Applies spin to particle velocities, ASSUMING RUBBLE
	 ** PILE POSITION AND VELOCITY ARE ZERO. 
	 */
	
	SSDATA *d;
	VECTOR v;
	int i;
	
	assert(rp != NULL);
	assert(rp->n_particles > 0);
	assert(rp->data != NULL);
	
	for (i=0;i<rp->n_particles;i++) {
		d = &rp->data[i];
		CROSS(rp->spin,d->pos,v);
		ADD_VEC(d->vel,v,d->vel);
		COPY_VEC(rp->spin,d->spin);
		}
}

/* end of rpu.c functions */

void
rubGetInterpParams(const COLLIDER *col_a,const COLLIDER *col_b,
				   INTERP_PARAMS *ip)
{
	/* based on a and b parameters finds b,v,mu, and eps to pass to interp */
	/* creates a new collider particle based on interp results */

	VECTOR vel,loc,cross;
	double loc_mag,v,R,M_tot,M_red,vcrit;

	/* ip->v must be in vcrit units */
	SUB_VEC(col_a->v,col_b->v,vel);
	v = MAG(vel);
	R = pow(col_a->fRadius,3.0)+pow(col_b->fRadius,3.0);
	R = pow(R,(1/3.0));
	M_tot = col_a->fMass+col_b->fMass;
	M_red = col_a->fMass*col_b->fMass/M_tot;
	vcrit = M_tot*sqrt(6.0/(5*M_red*R)); /* assumes G == 1 */
	ip->v = v/vcrit;
	if (col_a->fMass > col_b->fMass)
		ip->mu = col_b->fMass/col_a->fMass;
	else
		ip->mu = col_a->fMass/col_b->fMass;

/*	ip->eps = Z_EPS;*/
	/*use the velocity vector direction for impact parameter*/

	SUB_VEC(col_a->r,col_b->r,loc);
	loc_mag = MAG(loc);

	/*impact parameter in units of the seperation distance*/
	CROSS(vel,loc,cross);
	ip->b = MAG(cross);
	ip->b /= (loc_mag*v); 
	}

double
rubInterp(INTERP_PARAMS *ip, int n)
{
	int dim=N_SUBDIM;
	int i,j,k,l,ii,jj,kk,ll;
	int aa=0,bb=0,cc=0,dd=0;
	/* 4d subgrid of 2x2x2x2 - lin intp, extrap */

	double ysub[N_SUBDIM][N_SUBDIM][N_SUBDIM][N_SUBDIM]; 
	double mu_sub[N_SUBDIM],b_sub[N_SUBDIM],v_sub[N_SUBDIM],eps_sub[N_SUBDIM];
	double mu_int,b_int,v_int,eps_int,y_int,dy_int;

	const double mu[N_MU] = {1.0,1.0/3.0,1.0/6.0,1.0/9.0};
	const double b[N_B] = {0.0,0.25,0.50,0.75,1.00};
	const double v[N_VEL] = {0.50,1.00,1.50,2.00};
	const double eps[N_EPS] = {0.1,0.3,0.5,0.7,0.9};

	v_int = ip->v;
	b_int = ip->b;
	mu_int = ip->mu;
	eps_int = ip->eps;
	printf("v = %f, b = %f, mu = %f, eps = %f\n", v_int, b_int,mu_int,eps_int);
	/* interp (lin) mlrem if mu_int inside data grid extrap (lin) if outside */
	if (mu_int>mu[0]) {
		printf("Mu is greater than 1 - error\n");
		return 1;
	} else {
		for (i=0;i<(N_MU-1) && mu_int<=mu[i];i++)
			;
		}
	mu_sub[0]=mu[i-1];
	mu_sub[1]=mu[i];

	/* the impact parameter is fully contained - only interpolated */
	if (b_int>b[N_B-1]) {
		printf("The b value is too large for collision to take place\n");
		return 1;
	} else if (b_int<b[0]) {
		printf("The b value is less than zero -- error\n");
		return 1;
	} else {
		for (j=0;j<(N_B-1) && b_int>=b[j];j++)
			;
		}
	b_sub[0]=b[j-1];
	b_sub[1]=b[j];

	/* interp if v_int is inside data grid - extrap in both directions */
/*	if (v_int<=0.0) {*/
	if (v_int<0.0) {
		printf("The v value is less than zero - error.\n");
		return 1;
	} else {
		for (k=0;k<(N_VEL-1) && v_int>=v[k];k++)
			;
		}
	if (k==0) /* problem if extrapolating */
		k=1;
	v_sub[0]=v[k-1];
	v_sub[1]=v[k];
	/* did not add data for eps=0 or 1 - not straightforward and not nec. */
	if (eps_int>eps[N_EPS-1] || eps_int<eps[0]) {
		printf("The eps value is outside of the data grid\n");
		return 1;
	} else {
		for (l=0;l<(N_EPS-1) && eps_int>=eps[l];l++)
			;
		}
	eps_sub[0]=eps[l-1];
	eps_sub[1]=eps[l];

	for (ii=(i-1);ii<=i;ii++) {
		bb=0;
		for (jj=(j-1);jj<=j;jj++) {
			cc=0;
			for (kk=(k-1);kk<=k;kk++) {
				dd=0;
				for (ll=(l-1);ll<=l;ll++) {
					switch (n) {
					case 1:
						ysub[aa][bb][cc][dd] = Y_prim[ii][jj][kk][ll];
						break;
					case 2:
						ysub[aa][bb][cc][dd] = Y_sec[ii][jj][kk][ll];
						break;
					default:
						assert(0);
						}
					dd++;
					}
				cc++;
				}
			bb++;
			}
		aa++;
		}

	y_int = 0.0;

	polin2(mu_sub,b_sub,v_sub,eps_sub,ysub,dim,dim,dim,dim,
		   mu_int,b_int,v_int,eps_int,&y_int,&dy_int);

	if (y_int <= 0.0) /*to avoid returning neg mass*/
		y_int = 0.0;

	if (y_int > 1.0){ /*to avoid returning a mass greater than the sys mass*/
		y_int = 1.0;
		printf("y_int > 1 set to 1\n");
		}
	printf("y_int = %f\n", y_int);
	return y_int;
	}

void
rubCreatePart(const COLLIDER *col_a,const COLLIDER *col_b,
	      double mass_frac,COLLIDER **col_c)
{
	/*creates one particle that is used as the result of collision*/
	/*this function used when rubble_collide determines full resolution*/
	/*of collision is not needed*/

	VECTOR r_temp_a,r_temp_b,v_temp_a,v_temp_b,rcom,vcom;
	VECTOR rcom_a,rcom_b,vcom_a,vcom_b,L_a,L_b,Iw_a,Iw_b;
	double I_a,I_b,I_c,orig_den;

	*col_c = (COLLIDER *) malloc(sizeof(COLLIDER));
	assert(*col_c != NULL);

	(*col_c)->fMass = mass_frac*(col_a->fMass + col_b->fMass);
/*	(*col_c)->fRadius = pow(((*col_c)->fMass/(4.0/3.0*PI*DEN)),(1.0/3.0)); */
	orig_den = col_a->fMass/(4.0/3.0*PI*pow(col_a->fRadius,3.0));
	(*col_c)->fRadius = pow(((*col_c)->fMass/(4.0/3.0*PI*orig_den)),(1.0/3.0));
	COPY_VEC(col_a->r,r_temp_a);
	SCALE_VEC(r_temp_a,col_a->fMass);
	COPY_VEC(col_b->r,r_temp_b);
	SCALE_VEC(r_temp_b,col_b->fMass);
	ADD_VEC(r_temp_a,r_temp_b,rcom);
	COPY_VEC(rcom,(*col_c)->r);
	NORM_VEC((*col_c)->r,col_a->fMass + col_b->fMass);

	printf("mass_a = %e, mass_b = %e, mass_frac = %f\n",col_a->fMass,col_b->fMass,mass_frac);
	printf("col_a->v = %12e %12e %12e\n", col_a->v[0],col_a->v[1],col_a->v[2]);
	printf("col_b->v = %12e %12e %12e\n", col_b->v[0],col_b->v[1],col_b->v[2]);
	printf("col_a->r = %12e %12e %12e\n", col_a->r[0],col_a->r[1],col_a->r[2]);
	printf("col_b->r = %12e %12e %12e\n", col_b->r[0],col_b->r[1],col_b->r[2]);
	COPY_VEC(col_a->v,v_temp_a);
	SCALE_VEC(v_temp_a,col_a->fMass);
	COPY_VEC(col_b->v,v_temp_b);
	SCALE_VEC(v_temp_b,col_b->fMass);
	ADD_VEC(v_temp_a,v_temp_b,vcom);
	COPY_VEC(vcom,(*col_c)->v);
	NORM_VEC((*col_c)->v,col_a->fMass + col_b->fMass);
/*	NORM_VEC((*col_c)->v,(*col_c)->fMass);*/ /*conservation of linear momentum*/
/*	printf("(*col_c)->v = %12e %12e %12e\n", (*col_c)->v[0],(*col_c)->v[1],(*col_c)->v[2]);*/

	/* angular momentum conserved (spin+orbital) about com */
	SUB_VEC(col_a->r,(*col_c)->r,rcom_a);
	SUB_VEC(col_a->v,vcom,vcom_a);
	SUB_VEC(col_b->r,(*col_c)->r,rcom_b);
	SUB_VEC(col_b->v,vcom,vcom_b);

	CROSS(rcom_a,vcom_a,L_a);

	SCALE_VEC(L_a,col_a->fMass);
	I_a = 2.0/5.0*col_a->fMass*col_a->fRadius*col_a->fRadius;
	COPY_VEC(col_a->w,Iw_a);
	SCALE_VEC(Iw_a,I_a);
	ADD_VEC(L_a,Iw_a,L_a);

	CROSS(rcom_b,vcom_b,L_b);
	SCALE_VEC(L_b,col_b->fMass);
	I_b = 2.0/5.0*col_b->fMass*col_b->fRadius*col_b->fRadius;
	COPY_VEC(col_b->w,Iw_b);
	SCALE_VEC(Iw_b,I_b);
	ADD_VEC(L_b,Iw_b,L_b);

	ADD_VEC(L_a,L_b,(*col_c)->w);
	I_c = 2.0/5.0*(*col_c)->fMass*(*col_c)->fRadius*(*col_c)->fRadius;
	NORM_VEC((*col_c)->w,I_c);

	if (col_a->iRung >= col_b->iRung)
		(*col_c)->iRung = col_a->iRung;	
	else
		(*col_c)->iRung = col_b->iRung;
	}

static void
rubRanVel(double speed_max,VECTOR v)
{
	double v_max = speed_max/sqrt(N_DIM);

	do {
		SET_VEC(v,ran(),ran(),ran());
		SCALE_VEC(v,v_max);
	} while (MAG_SQ(v) > SQ(speed_max));
	}

void
rubCreateColliders(int nOut,RUBBLE_PILE *rp_target,int target_idx,
				   RUBBLE_PILE *rp_proj,int proj_idx,COLLIDER *col) 
{
	/*
	 ** Combines target & projectile rubble piles into single COLLIDER
	 ** structure for return to pkdDoCollisions(). Particles belonging
	 ** to each pile are identified by target_idx & proj_idx.
	 */

	int i,j;

	for (i=0;i<rp_target->n_particles;i++) {
		col[i].id.iIndex = rp_target->data[i].org_idx; /*should be iOrgIdx?*/
		col[i].id.iOrgIdx = target_idx;
		col[i].fMass = rp_target->data[i].mass;
		col[i].fRadius = rp_target->data[i].radius;
		COPY_VEC(rp_target->data[i].pos,col[i].r);
		COPY_VEC(rp_target->data[i].vel,col[i].v);
		COPY_VEC(rp_target->data[i].spin,col[i].w);
		}

	for (j=0;i<nOut;i++,j++) {
/*		col[i].id.iIndex = rp_proj->data[i].org_idx;*/
		col[i].id.iIndex = rp_proj->data[j].org_idx;
		col[i].id.iOrgIdx = proj_idx;
		col[i].fMass = rp_proj->data[j].mass;
		col[i].fRadius = rp_proj->data[j].radius;
		COPY_VEC(rp_proj->data[j].pos,col[i].r);
		COPY_VEC(rp_proj->data[j].vel,col[i].v);
		COPY_VEC(rp_proj->data[j].spin,col[i].w);
		}
	}

void
rubHcp(const COLLIDER *col,int max_npart,double *part_rad,RUBBLE_PILE *rp)
{
	VECTOR pos;
	double p_radius,dx,dy,dz,rmax,speed_esc;
	int np,nx,ny,nz,ix,iy,iz,i;

	rp->axis_len[0] = col->fRadius;
	rp->axis_len[1] = col->fRadius;
	rp->axis_len[2] = col->fRadius;
/*	printf("axis length = %e\n", col->fRadius);*/ /*DEBUG*/
/*	rp->density = DEN; *//*fix at some point*/
	rp->density = col->fMass/(4.0/3.0*PI*pow(col->fRadius,3.0));
	rp->radius = rp->axis_len[0];
	rp->mass = col->fMass;

	if (*part_rad == 0) { /* if rp = target start particle rad = zero*/
		rp->packing = PEN(rp->n_particles);
		*part_rad = rp->radius*pow(rp->packing/rp->n_particles,1.0/3);
	} else { /* rp = proj, particle radius set */
		rp->packing = PER(rp->radius/(*part_rad)); /* *p_radius = target_p_rad */
		rp->n_particles = rp->packing/CUBE(*part_rad/rp->radius);
/*		if (1.1*rp->n_particles > MAX_NPART) *//* a buffer to proj part num */
/*			rp->n_particles = MAX_NPART;	*//* 04.19.04 */
		if (1.1*rp->n_particles > max_npart)
			rp->n_particles = max_npart;
		else
			rp->n_particles = 1.1*rp->n_particles;

		/*printf("packing = %f, den = %f\n",
			   rp->packing, CUBE(*part_rad/rp->radius)); *//*DEBUG*/
		assert(rp->n_particles >= 0);
/*		printf("n_particles = %d\n", rp->n_particles); *//* DEBUG */
		if (rp->n_particles < 13) { /* min number of particles in rp is 13 */
			printf("WARNING: Projectile is small\n"); /* increased by 2 for */
			*part_rad = col->fRadius;	/* buffering n_particles */
			rp->n_particles = 1;
			}
		}

	p_radius = *part_rad; /* need a dummy var so *part_rad is not changed*/
/*	printf("rubble density = %f, particle density = %f\n",
		   rp->density, p_density); *//*DEBUG*/
/*	printf("rubble mass = %e, rubble radius = %f, particle radius = %f\n", 
		   rp->mass, rp->radius, p_radius);*//*DEBUG*/

	/* hexagonal closest packing - target */
	
	dx = 2*p_radius;
	dy = sqrt(3.0)*p_radius;
	dz = (2.0/3)*sqrt(6.0)*p_radius;
	
	rmax = rp->radius - p_radius;
	nx = rmax/dx;
	ny = rmax/dy;
	nz = rmax/dz;
	
	np = 0;
  
	for (ix=-nx;ix<=nx;ix++) 
		for (iy=-ny;iy<=ny;iy++) 
			for (iz=-nz;iz<=nz;iz++) {
				pos[X] = ix*dx + ((iy + iz)%2)*0.5*dx;
				pos[Y] = iy*dy - ((iz%4 - 2*SGN(iz))%2)*(dy/3);
				pos[Z] = iz*dz;
				if (rpuInEllipsoid(pos,p_radius,rp->axis_len)) {
					if (np == rp->n_particles) {
						(void) fprintf(stderr,"hcp(): Too many particles.\n");
						(void) fprintf(stderr, "np = %d, n = %d\n",
									   np, rp->n_particles);
						exit(1);
						}
					COPY_VEC(pos,rp->data[np].pos);
					++np;
					}
				}
  
	rp->n_particles = np;
/*	printf("n_particles = %d\n", rp->n_particles); *//*DEBUG*/

	/* Fill in remaining particle data fields */
  
	p_radius *= PART_SCALING; /*need to use a temp var here*/
  
	for (i=0;i<rp->n_particles;i++) {
		rp->data[i].mass = rp->mass/rp->n_particles;
		rp->data[i].radius = p_radius;
		ZERO_VEC(rp->data[i].vel); /* recalculated below */
		ZERO_VEC(rp->data[i].spin); /* taken care of with new rpuApplySpin */
		rp->data[i].org_idx = i;
		}

	/*
	 ** NOTE: many of the bulk property parameters in the rubble pile
	 ** are left uninitialized since we never need to use them.
	 */

 	/* Set center-of-mass position to zero */
 
	rpuCalcPos(rp);
	SCALE_VEC(rp->pos,-1);
	rpuTransPos(rp);
  
	/* Escape speed from particle surface */
	/*DEBUG: assumes equal sized particles*/
	speed_esc = sqrt(2*rp->data[0].mass/rp->data[0].radius);
  
	/* Set new particle velocities (random to prevent simultaneous collisions) */

	for (i=0;i<rp->n_particles;i++)
		rubRanVel(PART_MAX_SPEED*speed_esc,rp->data[i].vel);

	/* Set center-of-mass velocity to zero */

	rpuCalcVel(rp);
	SCALE_VEC(rp->vel,-1);
	rpuTransVel(rp);

	/*
	 ** Set spin ASSUMING RANDOM VELOCITIES CONTRIBUTE NEGLIGIBLY
	 ** (otherwise would need to compute angular momentum,
	 ** multiply by inverse of inertia tensor, and remove this
	 ** spin component first -- yuck).
	 */

	COPY_VEC(col->w,rp->spin);
	rpuApplySpin(rp);

	/* Now translate to collider frame */

	COPY_VEC(col->r,rp->pos);
	rpuTransPos(rp);
	COPY_VEC(col->v,rp->vel);
	rpuTransVel(rp);
	}

int
rubCreateRubblePiles(const COLLIDER *col_a,const COLLIDER *col_b,
					 const COLLISION_PARAMS *cp,COLLIDER **col_c,int *nOut)
{
	/*creates rubble piles in place of single particles col_a and col_b*/	

	const COLLIDER *target, *proj;
	RUBBLE_PILE rp_target, rp_proj;
	double p_radius=0;
	int max_npart,redo;

	if (col_a->fMass > col_b->fMass) { 
		target = col_a;
		proj = col_b;
	} else {
		target = col_b;
		proj = col_a;
		}

	/*
	 ** We require that dDensity be set to a sensible value
	 ** so that we don't get "puffy" projectiles...
	 */
	assert(target->fRadius >= 0.9/*DEBUG some slack*/*proj->fRadius);

	redo = 0;
	max_npart = MAX_NPART;

	rp_target.data = rp_proj.data = NULL; /* force malloc() on first rpRealloc() call */
	do {
		redo = 0;
		rp_target.n_particles = max_npart;
		rpuRealloc(&rp_target); 
/*		printf("Target Radius = %e\n", target->fRadius);*/ /*DEBUG*/
		printf("max_npart=%i,rp_target.n_particles=%i, p_radius=%e\n",
			   max_npart,rp_target.n_particles, p_radius);
		rubHcp(target,max_npart,&p_radius,&rp_target); /*initially p_radius = 0*/
		printf("max_npart=%i,rp_target.n_particles=%i, p_radius=%e\n",
			   max_npart,rp_target.n_particles, p_radius);
/*		printf("Returned from first call to rubHcp()\n");*/ /*DEBUG*/
/*		printf("Particle radius = %e, target radius = %e\n",
			   p_radius,target->fRadius);*/ /*DEBUG*/
		assert(rp_target.n_particles <= max_npart);
		rp_proj.n_particles = rp_target.n_particles; /*to allow malloc*/
		printf("rp_proj=%i,rp_target.n_particles=%i\n",
			   rp_proj.n_particles,rp_target.n_particles);
		rpuRealloc(&rp_proj); 
/*		printf("Particle radius = %e, target radius = %e\n",
			   p_radius,proj->fRadius);*//*DEBUG*/
		rubHcp(proj,max_npart,&p_radius,&rp_proj); /*to ensure proj and target = proj part*/
/*		printf("Returned from second call to rubHcp() \n");*/
		assert(rp_proj.n_particles <= max_npart); /*DEBUG 07.17.04*/
		/*
		 ** If the rubble pieces are larger than the simulation resolution,
		 ** then any survivors could be turned into planetesimals.  But these
		 ** bodies tend to have higher densities because of the packing
		 ** inefficiency of a rubble pile.  This in turn means the timesteps
		 ** should decrease.  Because of these issues, we DO NOT ALLOW rubble
		 ** pieces to be smaller than the resolution.  The only way around this
		 ** then is to use more particles in the rubble piles when the
		 ** planetesimal masses start getting large...
		 */
		if (rp_target.mass/rp_target.n_particles > cp->DB.dRubMinMass ||
			rp_proj.mass/rp_proj.n_particles > cp->DB.dRubMinMass) { /* assumes equal masses! */
			printf("max mass = %e, target part %e, proj part %e,targ = %i, proj = %i\n", 
				   cp->DB.dRubMinMass, rp_target.mass/rp_target.n_particles,
				   rp_proj.mass/rp_proj.n_particles,rp_target.n_particles,
				   rp_proj.n_particles); /*DUBUG*/
			printf("rubCreateRubblePiles(): growing max_npart\n");
			assert(max_npart < SUPER_NPART); /* for now */
			max_npart = SUPER_NPART;
			printf("super_npart = %i\n", max_npart);
			p_radius = 0.0; /*DEBUG*/
			redo = 1;
			}
		} while (redo);
	rpuRealloc(&rp_target);
	rpuRealloc(&rp_proj);	

	/* need to modify the usage of this so that when two part collide */
	/* pkdgrav does not call rubble.c again */
  
	*nOut = rp_target.n_particles+rp_proj.n_particles;
	printf("nOut = %d\n", *nOut); 
	*col_c = (COLLIDER *) malloc((*nOut)*sizeof(COLLIDER));
	assert(*col_c != NULL);

	rubCreateColliders((*nOut),&rp_target,target->id.iOrgIdx,
					   &rp_proj,proj->id.iOrgIdx,*col_c);
	{
	/*DEBUG is there a better way/place to do this?*/
	/*(needed for PutColliderInfo())*/
	/*
	 ** Basically, we would like colliding particles to be
	 ** on the same (smallest?) rung, so that we can be sure
	 ** the encounter was treated accurately.  The rung
	 ** inheritance below is overridden by pkdRubbleStep(),
	 ** but it's important to initialize the rungs of the
	 ** collision remnants to *something* to avoid FPE issues.
	 */
	int i;
	if (col_a->iRung != col_b->iRung)
		(void) printf("WARNING: %i on rung %i but %i on rung %i\n",
					  col_a->id.iOrder,col_a->iRung,
					  col_b->id.iOrder,col_b->iRung);
	/*assert(col_a->iRung == col_b->iRung);*/
	for (i=0;i<*nOut;i++)
		(*col_c)[i].iRung = col_a->iRung; /* arbitrary -- reset by pkdRubbleStep() */
	}

	rpuFree(&rp_proj);
	rpuFree(&rp_target);

	return 0;
	}     
#define TOLERANCE 1.0e-16

void
rubRubbleCollide(const COLLIDER *col_a,const COLLIDER *col_b,
				 const COLLISION_PARAMS *cp,COLLIDER **col_c,
				 double *dMassInDust,int *nOut)
{
	/*determines whether full resolution of planetesimal collision is needed*/

	/*
	 ** NOTE: particle colors and other identifying info handled by caller.
	 ** This means, for example, the rubble pile structure has NO color field!
	 */
	
	static int first_call = 1;

	INTERP_PARAMS ip;
/* if f > 1 no need to define w2 and w2max */
	double m,m_prim; /*w2,w2max;*/
/*	double m,m_prim;*/

	*dMassInDust = 0.0;

	if (first_call) { /* initialize random number generator on first call */ 
		(void) ran(); 
		first_call = 0; 
		} 

	printf("m %e r %e dens1 = %e, m %e r %e dens2 = %e\n", 
		   col_a->fMass, col_a->fRadius,
		   col_a->fMass/(4.0/3.0*M_PI*col_a->fRadius*col_a->fRadius*col_a->fRadius),
		   col_b->fMass, col_b->fRadius,
		   col_b->fMass/(4.0/3.0*M_PI*col_b->fRadius*col_b->fRadius*col_b->fRadius));
	ip.eps = cp->dEpsN;
	rubGetInterpParams(col_a,col_b,&ip); 
	/*printf("%f %f %f %f\n", ip.b, ip.v, ip.mu, ip.eps);	   */ /*DEBUG*/
	m = rubInterp(&ip,1); /* get mass of largest remnant */ 
	m_prim = m;

/*	printf("primary mass = %f, min_mass = %f\n",m,cp->dRubbleMinFracMass);*/ 

	if (1.0 - m >= cp->dRubbleMinFracMass) {
		if (m < SMALL_MASS) { /* 8.28.03 if interp puts m to zero go to rp */
			printf("WARNING: Primary mass very small\n");
			rubCreateRubblePiles(col_a,col_b,cp,col_c,nOut);
			assert(*nOut > 1);
			return;
		} else {	
		m = rubInterp(&ip,2); /* get mass of second-largest remnant */
		printf("Returned from secondary mass test, m = %f\n",m);
		}
		if (m > cp->dRubbleMinFracMass) {
			/*did prim and sec change mass much?*/
			/*how should this be done?*/
			printf("Creating rubble piles\n");
			rubCreateRubblePiles(col_a,col_b,cp,col_c,nOut);
			assert(*nOut > 1);
			return;
			}
		}
	/*	printf("Got to the else condition\n");*/
	/*primary represents a significant mass of the entire system*/
	rubCreatePart(col_a,col_b,m_prim,col_c); /*rest goes to dust*/	
	printf("col_c->v[0] = %e, col_c->v[1] = %e, col_c->v[2] = %e\n",
		   (*col_c)->v[0], (*col_c)->v[1], (*col_c)->v[2]); 
	printf("col_c->r[0] = %e, col_c->r[1] = %e, col_c->r[2] = %e\n",
		   (*col_c)->r[0], (*col_c)->r[1], (*col_c)->r[2]); 
	*nOut = 1;

	/*test largest rem is not spinning too fast*/
	/*if is force direct simulation of collision*/
/* the following is commented out when f > 1 */
/*	w2max = (*col_c)->fMass/((*col_c)->fRadius*(*col_c)->fRadius*(*col_c)->fRadius);
	w2 = 0.0;
	for (i=0;i<3;i++)
		w2 += (*col_c)->w[i]*(*col_c)->w[i];
	if (w2 > 0.9*w2max) {*/
/*	if (w2 > w2max) {*/
/*		free((void *) *col_c);*/ /*free space used when part created*/
/*		printf("WARNING: Primary spinning too fast forcing resolved collision\n");
		rubCreateRubblePiles(col_a,col_b,cp,col_c,nOut);
		assert(*nOut > 1);
	} else { */
/* end of spin test that is commented out when f > 1 */
/*DEBUG 05.26.04*/
	/*printf("No spin test\n");*/
/*		printf("dMassInDust = %e, before perfect merge test\n",*dMassInDust);*/
		if ((col_c[0]->fMass >= (col_a->fMass+col_b->fMass)) &&
			(col_c[0]->fMass < 1.000001*(col_a->fMass+col_b->fMass))){
			printf("col_c[0]->fMass = %1.20e\n", col_c[0]->fMass);
			col_c[0]->fMass = col_a->fMass+col_b->fMass;
			*dMassInDust = 0.0;
			}
		else
			*dMassInDust = col_a->fMass + col_b->fMass - col_c[0]->fMass;

/*		*dMassInDust = col_c[0]->fMass - (col_a->fMass + col_b->fMass);
		if (fabs(*dMassInDust)/col_c[0]->fMass < TOLERANCE) *dMassInDust = 0.0;*/
		printf("dMassInDust = %e\n",*dMassInDust);
/*		col_c[0]->id.iOrder = col_a->id.iOrder;*/
/*		printf("After rubCreatePart\n");*/
/*		printf("col_a.id.iOrder=%d,col_b.id.iOrder=%d\n",
			   col_a->id.iOrder,col_b->id.iOrder);*/
/*		for (i=0;i<*nOut;i++) 
			printf("col_c[%d].id.iOrder=%d,iIndex=%d,iOrgIdx=%d\n",
				   i,col_c[i]->id.iOrder,col_c[i]->id.iIndex,
				   col_c[i]->id.iOrgIdx);*/
		assert(*dMassInDust >= 0.0);
/*		}*/ /*DEBUG 05.26.04*/
	}
#undef TOLERANCE
/*
 ** Following routines adapted from Numerical Recipes in C (2nd ed).
 ** from rpu.c dcr 98-09-24
 ** zml added to rubble.c 03-05-02 
 */

/* zml modified from NumRec in C */
void polin2(double x1a[], double x2a[], double x3a[], double x4a[], 
			double ya[N_SUBDIM][N_SUBDIM][N_SUBDIM][N_SUBDIM], int l,
			int m, int n, int p, double x1, 
			double x2,double x3, double x4, double *y, double *dy)
{
	void polint(double xa[],double ya[],int p,double x,double *y,double *dy);
	int h,i,j;
	double ***ymtmp,**yymtmp,*yyymtmp;

	/*converting this num rec routine back to normal C 11.06.02*/
	ymtmp = (double ***) malloc(l*sizeof(double **));
	assert(ymtmp != NULL);
	for (h=0;h<l;h++) {
		ymtmp[h] = (double **) malloc(m*sizeof(double *));
		assert(ymtmp[h] != NULL);
		for (i=0;i<m;i++) {
			ymtmp[h][i] = (double *) malloc(n*sizeof(double));
			assert(ymtmp[h][i] != NULL);
			} 
		}

	yymtmp = (double **) malloc(l*sizeof(double *));
	assert(yymtmp);
	for (i=0;i<l;i++) {
		yymtmp[i] = (double *) malloc(m*sizeof(double));
		assert(yymtmp[i] != NULL);
		}

	yyymtmp = (double *) malloc(l*sizeof(double));
	assert(yyymtmp != NULL);
  
	for (h=0;h<l;h++)
		for (i=0;i<m;i++)
			for (j=0;j<n;j++) 
				polint(x4a,ya[h][i][j],p,x4,&ymtmp[h][i][j],dy);

	for (i=0;i<l;i++)
		for (j=0;j<m;j++) 
			polint(x3a,ymtmp[i][j],n,x3,&yymtmp[i][j],dy);

	for (j=0;j<l;j++) 
		polint(x2a,yymtmp[j],m,x2,&yyymtmp[j],dy);

	polint(x1a,yyymtmp,l,x1,y,dy);
	for (i=0;i<l;i++) {
		for (j=0;j<m;j++)
			free((void *) ymtmp[i][j]);
		free((void *) ymtmp[i]);
		}
	free((void *) ymtmp);

	for (i=0;i<l;i++)
		free((void *) yymtmp[i]);
	free((void *) yymtmp);

	free((void *) yyymtmp);
	}

void polint(double xa[], double ya[], int n, double x, double *y, double *dy)
{
	/*converted this num rec routine to normal C 11.06.02*/
	int i,m,ns=0;
	double den,dif,dift,ho,hp,w;
	double *c,*d;
  
	dif=fabs(x-xa[0]);
	c = (double *) malloc(n*sizeof(double));
	assert(c != NULL);
	d = (double *) malloc(n*sizeof(double));
	assert(d != NULL);

	for (i=0;i<n;i++) {
		if ( (dift=fabs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
			}
		c[i]=ya[i];
		d[i]=ya[i];
		}
	*y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=0;i<(n-m);i++) {
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ( (den=ho-hp) == 0.0) {
				printf("Error in routine polint\n");
				exit (1);
				}
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
			}
		*y += (*dy=(2*ns < (n-m-1) ? c[ns+1] : d[ns--]));
		}

	free((void *) d);
	free((void *) c);
	}
/* NRiC 2nd ed */

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

static double
ran(void)
{
	static unsigned long n_times=0;
	static long seed;

	if (!n_times) {
		seed = getpid();
		printf("ran seed = %ld\n", seed);
		srangen(seed);
		}

	n_times++;
	return rangen();
	}

#ifdef RANOLD /*DEBUG*/
static double
ran_old(void)
{
	static long idum;
	
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (!iy) { /* throw away returned value in this case! */
		idum = getpid();
		idum = 421; /* DEBUG: want to run same simulation (11.25.03)*/
		printf("idum = %ld\n", idum);
		for (j=NTAB+7;j>=0;j--) {
			k=idum/IQ;
			idum=IA*(idum-k*IQ)-IR*k;
			if (idum < 0) idum += IM;
			if (j < NTAB) iv[j] = idum;
			}
		iy=iv[0];
		}
	k=idum/IQ;
	idum=IA*(idum-k*IQ)-IR*k;
	if (idum < 0) idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return (double) temp;
	}
#endif

#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

/*#undef NRANSI*/
/* (C) Copr. 1986-92 Numerical Recipes Software ?421.1-9. */

#endif /* RUBBLE_ZML */
