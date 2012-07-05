/*
 ** Function to solve the equation M = e*sinh(F) - F for the value of F.
 ** This is used in solving the analogue of Kepler's equation for 
 ** for Hyperbolic orbits. It is mainly called by conversion from the
 ** Delaunay elements to cartesian coordinates.
 **
 ** Joachim Stadel, Jan. 11, 1995
 */
#include <math.h>
#include <assert.h>

#ifdef CRAY_T3D
#include "hyperlib.h"
#endif

#define MAX_HYPERBOLA_ITTR		32
/*
 ** This upper bound should be plenty because it corresponds to
 ** a mean anomally of greater than 1e+42! Also assures that the 
 ** hyperbolic funtions don't blow up!
 */
#define MAX_F					100.0

double dHypAnom(double M,double e)
{
	double w,w1,w2,w3;
	double Y,crit,sgm;
#ifdef __i386__
        long double dF,F=0,lb,ub=0;
#else
        double dF,F=0,lb,ub=0;
#endif
	int i;

	/*
	 ** Find a good first approximation to F!
	 */
	if (M < 0.0) {
		sgm = -1.0;
		M = -M;
		}
	else {
		sgm = 1.0;
		}
	crit = 5.0*e - 2.5;
	if (e > 1.000001) {
		if (M < crit) {
			Y = sqrt(8.0*(e - 1.0)/e);
			w = asinh(3.0*M/(Y*(e - 1.0)))/3.0;
			F = Y*sinh(w);
			if (F > MAX_F) F = MAX_F;
			ub = F;
			}
		else {
			ub = MAX_F;
			F = log(2.0*M/e);
			}
		}
	else if (e >= 1.0) {
		F = pow(6.0*M/e,1.0/3.0);
		if (F > MAX_F) F = MAX_F;
		ub = F;
		}
	else {
		/*
		 ** Bad value of e, not appropriate for Hyperbolic orbit!
		 */
		assert(0);
		}
	lb = 0.0;
	for (i=0;i<MAX_HYPERBOLA_ITTR;++i) {
		w2 = e*sinh(F);
		w3 = e*cosh(F);
		w = M + F - w2;
		w1 = w3 - 1.0;
		dF = w/w1;		
		dF = w/(w1 - dF*w2);
		dF = w/(w1 - dF*(0.5*w2 + dF*w3/6.0));
		if (dF < 0.0) ub = F;
		else lb = F;
		dF += F;
		if (F == dF) return(sgm*F);
		if (dF > lb && dF < ub) F = dF;
		else F = 0.5*(lb + ub);
		if (F == lb || F == ub) return(sgm*F);
		}
	/*
	 ** Hyperbola solution has failed. Please note M and e which caused
	 ** the failure and report to the author.
	 */
	assert(0);
	return(0.0);
	}
