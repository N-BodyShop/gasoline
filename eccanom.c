/*
 ** This function solves Kepler's equation, M = E - e*sin(E), for E.
 ** It is mainly called by conversion from the Delaunay elements to 
 ** cartesian coordinates.
 ** Original code from the Saha and Tremaine long term planetary integrator. 
 **
 ** Joachim Stadel, Jan. 11, 1995
 */
#include <math.h>
#include <assert.h>

#define MAX_KEPLER_ITTR		32

double dEccAnom(double M,double e)
{
	double w,w1,w2,w3,twopi;
#ifdef __i386__
	long double dE,E,lb,ub;
#else
	double dE,E,lb,ub;
#endif
	int i;

	twopi = 2.0*M_PI;
	M -= ((int)(M/twopi))*twopi;
	if (sin(M) > 0.0) E = M - 0.85*e;  /* sin(M) is -ve */
	else E = M + 0.85*e;  /* sin(M) is +ve */
	lb = -twopi;
	ub = twopi;
	for (i=0;i<MAX_KEPLER_ITTR;++i) {
		w2 = e*sin(E);
		w3 = e*cos(E);
		w = M + w2 - E;
		w1 = 1.0 - w3;
		dE = w/w1;
		dE = w/(w1 + 0.5*dE*w2);
		dE = w/(w1 + dE*(0.5*w2 + dE*w3/6.0));
		if (dE < 0.0) ub = E;
		else lb = E;
		dE += E;
		if (E == dE) return(E);
		if (dE > lb && dE < ub) E = dE;
		else E = 0.5*(lb + ub);
		if (E == lb || E == ub) return(E);
		}
	/*
	 ** Kepler solution has failed.
	 */
	assert(0);
	return(0.0);
	}
