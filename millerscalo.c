#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include "millerscalo.h"

void MSInitialize(MSPARAM *pms) 
{
    MSPARAM ms;
    struct MillerScaloContext initms = 
	{ 	42.0,	-0.4, .1, /* parameters from Ap.J. Supp., 41,1979 */
		42.0,	-1.5, 1.0,
		240.0,	-2.3, 10.0,
		100.0};
    
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
    
    assert(argc == 2);
    
    nsamp = atoi(argv[1]);
    dlgm = (2.0 + 1.0)/nsamp;
    
    Ntot = dMSCumNumber(&MSparam, 0.0);
    Mtot = dMSCumMass(&MSparam, 0.0);
    
    printf("%g %g\n", Ntot, Mtot);
    
    for(i = 0; i < nsamp; i++) {
	double mass;
	
	lgm = -1 + i*dlgm;

	mass = pow(10.0, lgm);
	printf("%g %g %g %g\n", mass, dMSIMF(&MSparam, mass),
	       dMSCumNumber(&MSparam, mass)/Ntot,
	       dMSCumMass(&MSparam, mass)/Mtot);
	}
    return 0;
    }
#endif
