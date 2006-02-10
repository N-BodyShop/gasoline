#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include "millerscalo.h"


/* Miller-Scalo IMF (Miller & Scalo, Ap.J. Supp., 41, 513, 1979) in
   stars per unit logarithmic mass.  Divide by M (mass) for IMF in
   stars per unit mass.  Also IMF is defined per yer per pc^2,
   integrated over a cylinder that extends "several hundred parsecs on
   either side of the plane of the galaxy" */

void MSInitialize(MSPARAM *pms) 
{
    MSPARAM ms;
    
#ifdef KROUPA
    struct MillerScaloContext initms = 
         {       0.3029*1.86606, -0.3, .08, /* parameters from A&A, 315,1996 */
                 0.3029, -1.2, 0.5,
                 0.3029, -1.7, 1.0,
                 100.0};
#else
    struct MillerScaloContext initms = 
	{ 	0.5652,	-0.3, .08, /* parameters from A+A 315 105 Raiteri et al */
		0.3029,	-1.2, 0.5,
		0.3029,	-1.7, 1.0, 
		100.0};
/*	{ 	42.0,	-0.4, .1, /* parameters from Ap.J. Supp., 41,1979 /
		42.0,	-1.5, 1.0,
		240.0,	-2.3, 10.0, /*This is discontinuous, but is what /
		100.0};		    /*they report in paper, s we leave it.*/
#endif

    ms = (MSPARAM) malloc(sizeof(struct MillerScaloContext));
    
    assert(ms != NULL);
    *ms = initms;
    *pms = ms;
    }

double dMSIMF(MSPARAM p, double mass)
{
    double dIMF;
    
    if(mass > p->mmax)
	return 0.0;
    if(mass > p->m3)
	dIMF = p->a3*pow(mass, p->b3);
    else if(mass > p->m2)
	dIMF = p->a2*pow(mass, p->b2);
    else if(mass > p->m1)
	dIMF = p->a1*pow(mass, p->b1);
    else
	dIMF = 0.0;
    return dIMF;
    }


/*
 * Cumulative number of stars with mass greater than "mass".
 */
double dMSCumNumber(MSPARAM p, double mass)
{
    double dCumN;
    
    if(mass > p->mmax)
	return 0;
    if(mass > p->m3)
	return p->a3/p->b3*(pow(p->mmax, p->b3) - pow(mass, p->b3));
    else
	dCumN = p->a3/p->b3*(pow(p->mmax, p->b3) - pow(p->m3, p->b3)); 
    if(mass > p->m2) {
	dCumN += p->a2/p->b2*(pow(p->m3, p->b2) - pow(mass, p->b2));
	return dCumN;
    }
    else {
	dCumN += p->a2/p->b2*(pow(p->m3, p->b2) - pow(p->m2, p->b2));
	}
    if(mass > p->m1)
	dCumN += p->a1/p->b1*(pow(p->m2, p->b1) - pow(mass, p->b1));
    else
	dCumN += p->a1/p->b1*(pow(p->m2, p->b1) - pow(p->m1, p->b1));
    
    return dCumN;
    }

/*
 * Cumulative mass of stars with mass greater than "mass".
 */
double dMSCumMass(MSPARAM p, double mass)
{
    double dCumM;
    
    if(mass > p->mmax)
	return 0;
    if(mass > p->m3)
	return p->a3/(p->b3 + 1)*(pow(p->mmax, p->b3 + 1)
				  - pow(mass, p->b3 + 1));
    else
	dCumM = p->a3/(p->b3 + 1)*(pow(p->mmax, p->b3 + 1)
				   - pow(p->m3, p->b3 + 1)); 
    if(mass > p->m2) {
	dCumM += p->a2/(p->b2 + 1)*(pow(p->m3, p->b2 + 1)
				    - pow(mass, p->b2 + 1));
	return dCumM;
	}
    else {
	dCumM += p->a2/(p->b2 + 1)*(pow(p->m3, p->b2 + 1)
				    - pow(p->m2, p->b2 + 1));
	}
    if(mass > p->m1)
	dCumM += p->a1/(p->b1 + 1)*(pow(p->m2, p->b1 + 1)
				    - pow(mass, p->b1 + 1));
    else
	dCumM += p->a1/(p->b1 + 1)*(pow(p->m2, p->b1 + 1)
					- pow(p->m1, p->b1 + 1));
    
    return dCumM;
    }

double imf1to8Exp(MSPARAM p)
{
  if(p->m3 < 3.0) return p->b3;
  else if(p->m2 < 3.0) return p->b2;
  else return p->b1;
}

double imf1to8PreFactor(MSPARAM p)
{
  if(p->m3 < 3.0) return p->a3;
  else if(p->m2 < 3.0) return p->a2;
  else return p->a1;
}

#if 0
int
main(int argc, char **argv)
{
    int nsamp;
    int i, imax;
    double dlgm;
    double lgm;
    double Ntot;
    double Mtot;
    double dMassT1, dMassT2, dMass, part, sum;
    
    MSPARAM MSparam;
    MSInitialize (&MSparam);

    sum = 0;
    dMassT1 = 3;
    dMassT2 = 8;
    imax = 101;
    for (i = 0; i < imax; i++) {
        dMass = dMassT1 + i*(dMassT2-dMassT1)/(float) imax;
        part = dMSIMFSec (MSparam, dMass);
        printf ("%g %g\n", dMass, part);
        sum += part*0.05;
    }
    printf ("sum = %g\n", sum);

    printf ("number SN Ia = %g\n", dNSNIa (MSparam, dMassT1, dMassT2));
    printf ("mass in SN Ia %g\n", dMSNIa (MSparam, dMassT1, dMassT2));

    assert(argc == 2);
    
    nsamp = atoi(argv[1]);
    dlgm = (2.0 + 1.0)/nsamp;
    
    Ntot = dMSCumNumber(MSparam, 0.0);
    Mtot = dMSCumMass(MSparam, 0.0);
    
    printf("Ntot, Mtot = %g %g\n", Ntot, Mtot);
    
    for(i = 0; i < nsamp; i++) {
	double mass;
	
	lgm = -1 + i*dlgm;

	mass = pow(10.0, lgm);


	printf("%g %g %g %g\n", mass, dMSIMF(MSparam, mass),
	       dMSCumNumber(MSparam, mass),
	       dMSCumMass(MSparam, mass));
	}
    return 0;
    }
#endif
