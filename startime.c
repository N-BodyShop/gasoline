#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <assert.h>
#include <stdlib.h>
#include "startime.h"

#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))

void PadovaInitialize(PDVAPARAM *ppdva) 
{
    PDVAPARAM pdva;
/* parameters from Raiteri, Villata, & Navarro, A&A, 315, 105, 1996 */
    struct PadovaContext initpdva = 
	{ 	10.13,	0.07547, -0.008084, 
		-4.424,	-0.7939, -0.1187,
		1.262,	0.3385, 0.05417,
                0.0, 0.0, 0.0,
                7e-5, 3e-2, 0.02, 9e-3};

    pdva = (PDVAPARAM) malloc(sizeof(struct PadovaContext));
    assert(pdva != NULL);
    *pdva = initpdva;
    *ppdva = pdva;
    }


void PadovaCoefInit (PDVAPARAM ppdva, double dMetals)
{
    double logZ;
    double Z = dMetals;
    
    Z = min(Z, ppdva->zmax);
    Z = max(Z, ppdva->zmin);
    logZ = log10(Z);
    ppdva->a0 = ppdva->a00 + ppdva->a01*logZ + ppdva->a02*logZ*logZ;
    ppdva->a1 = ppdva->a10 + ppdva->a11*logZ + ppdva->a12*logZ*logZ;
    ppdva->a2 = ppdva->a20 + ppdva->a21*logZ + ppdva->a22*logZ*logZ;
}

double dSTLtimeMStar (PDVAPARAM ppdva, double dStarMass,
		      double dMetals) 
{
    /* finds stellar lifetime in yr corresponding to stellar 
       mass in solar masses */
    double logStarMass, logLtime, Ltime;
    
    logStarMass = log10(dStarMass);
    PadovaCoefInit (ppdva, dMetals);
    logLtime = ppdva->a0 + ppdva->a1 *logStarMass + ppdva->a2*logStarMass*logStarMass;
    Ltime = pow(10.0, logLtime);
    
    return Ltime;
}

double dSTMStarLtime (PDVAPARAM ppdva, double dStarLtime, double dMetals) 
{
    /* finds stellar mass in solar masses corresponding to stellar 
       lifetime dStarTime in yr */
    double logStarMass, StarMass;
    double a, b, c;
    
    if(dStarLtime <= 0.0)	/* Time can be zero */
	return DBL_MAX;
    
    PadovaCoefInit (ppdva, dMetals);
    c = ppdva->a0;
    c -= log10(dStarLtime);
    b = ppdva->a1;
    a = ppdva->a2;
    if(b*b - 4*a*c < 0.0) 	/* time is too small for fitting
				   formula */
	return DBL_MAX;

    logStarMass = (-b - sqrt(b*b - 4*a*c))/(2*a);
    StarMass = pow(10., logStarMass);
    return StarMass;
}

    
#if 0
int
main(int argc, char **argv)
{
    int nsamp;
    int i;
    double dlgm;
    double lgm;
    double Ntot;
    double Mtot;

    PDVAPARAM ppdva;
    SN sn;
    double dStarMass;
    double mass[30], time[30];
    PARTICLE *p;
    
    assert(argc == 2);
    
    nsamp = atoi(argv[1]);
    dlgm = (2.0 + 1.0)/nsamp;
    snInitialize (&sn);
//    Ntot = dMSCumNumber(&MSparam, 0.0);
//    Mtot = dMSCumMass(&MSparam, 0.0);
    
//    dSTLtimeMStar (ppdva, sn, dStarMass, p);
//    dSTMStarLtime (ppdva, sn, dStarLtime, p);
    
    PadovaInitialize(&ppdva);
//    mass = 8;
//    printf ("first %g %g\n", mass, dSTLtimeMStar (ppdva, mass, p));

    for(i = 0; i < nsamp; i++) {
	
//	lgm = -1 + i*dlgm;
	lgm = i*dlgm;

	mass[i] = pow(10.0, lgm);
        time[i] = dSTLtimeMStar (ppdva, sn, mass[i], p);
	printf("%g %g\n", mass[i], time[i]);        

	}

    for(i = 0; i < nsamp; i++) {

        mass[i] = dSTMStarLtime (ppdva, sn, time[i], p);
	printf("%g %g\n", mass[i], time[i]);        
    }
    
    return 0;
    }
#endif
