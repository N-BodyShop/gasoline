#ifdef STARFORM
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "supernova.h"
#include "millerscalo.h"
#include "supernovaia.h"

/* Calculation of number and mass of stars that go SN Type Ia in a
   given stellar mass range.  Must be in separate file from
   supernova.c and millerscalo.c because functions require MSSN
   structure that packages together MillerScaloContext structure and
   snContext structure.  SN Ia functions need both the Miller-Scalo
   IMF parameteres and the SN Ia parameters from supernova.h.  
   
   XXX - Might want to move constants relating to SN Ia to new
   structure in supernovaia.  This would require separate
   initialization, etc. same changes to msrInitialize as before. */

#define EPSSNIA 1e-7

#include "romberg.h"

double dNSNIa (MSSN mssn, double dMassT1, double dMassT2)
{
    assert (dMassT1 < dMassT2 && dMassT2 <= mssn->sn.dMBmax/2.);
    
	/* calculate number of SN Type Ia a la Raiteri, Villata, Navarro, A&A
	   315, 105, 1996) Returns number of SN Type Ia that occur during
	   timestep in which dMassT1 and dMassT2 are masses of stars that end
	   their lives at the end and beginning of timestep, respectively */
    return dRombergO(mssn, (double (*)(void *, double)) dMSIMFSec, dMassT1, dMassT2, EPSSNIA);
	}

#if 0
/* XXX The following is wrong and not used */
double dMSNIa (MSSN mssn, double dMassT1, double dMassT2)
{
    assert (dMassT1 < dMassT2 && dMassT1 >= mssn->sn.dMBmin && dMassT2 <= mssn->sn.dMBmax/2.);
    
	/* calculate mass of stars that go SN Type Ia a la Raiteri, Villata,
	   Navarro, A&A 315, 105, 1996) Returns total mass in stars that go SN
	   Type Ia during timestep in which dMassT1 and dMassT2 are masses of
	   stars that end their lives at the end and beginning of timestep,
	   respectively */
    return dRombergO(mssn, (double (*)(void *, double)) dMSIMFSecM, dMassT1, dMassT2, EPSSNIA);
	}

/* XXX - move dMBmin, dMBmax in dNSNIa and dMSNIa
   XXX - move dFracBinSNIa in dMSIMFSec and dMSIMFSecM
   */

/* Not really used any more because mass SNIa ejecta = NSNIa*1.4 */
double dMSIMFSecM(MSSN mssn, double dMass2)
{
	/* Mass times IMF of secondary in binary system that goes SN Ia */
    double dIMFSecExp, dIMFSecIntExp, dIMFSec;
    double dMass2_2, Msup, Minf;
    
    dIMFSecIntExp = imf1to8Exp(&mssn->ms) - 2.; /* subtract 1 from exponent in integrand
                                      because MS IMF is per unit log mass,
                                      subtract 2 because of square of ratio of mass
                                      of secondary to mass of binary, add 1 to multiply 
                                      IMF by mass */
    /*dIMFSecIntExp = dIMFSecExp + 1.;  exponent of anti-derivative is +1 */
    Msup = dMass2 + 8;
    dMass2_2 = 2.*dMass2;    
    Minf = (dMass2_2 > 3)?dMass2_2:3;
    dIMFSec = pow (Msup, dIMFSecIntExp) - pow(Minf, dIMFSecIntExp);
    dIMFSec *= mssn->sn.dFracBinSNIa * imf1to8PreFactor(&mssn->ms) * dMass2*dMass2*dMass2 / dIMFSecIntExp;
    return dIMFSec;
    }
#endif

/** IMF of secondary in binary system that goes SN Ia
 * The distribution of secondary mass ratios is assumed to be a power
 * law, mu^dGamma as in Matteucci & Greggio (1986)
 */

double dMSIMFSec(MSSN mssn, double dMass2)
{
    double dIMFSecExp, dIMFSecIntExp, dIMFSec;
    double dMass2_2, Msup, Minf;
    const double dGamma_2nd = 2.0;
    /* Normalization of 2ndary */
    const double dNorm_2nd = pow(2.0, 1 + dGamma_2nd)*(1 + dGamma_2nd);

    
    /* subtract 1 from exponent in integrand because MS IMF is per
       unit log mass.  Also multiply in ratio of masses */
    dIMFSecExp = imf1to8Exp(&mssn->ms) - (1. + dGamma_2nd);
    dIMFSecIntExp = dIMFSecExp + 1.; /* exponent of anti-derivative is +1 */
    Msup = dMass2 + 8; /* Mass where primary would have gone
			  supernova. */
    dMass2_2 = 2.*dMass2;  /* Minimum mass of binary */    
    Minf = (dMass2_2 > 3.0)?dMass2_2:3.0;
    dIMFSec = pow (Msup, dIMFSecIntExp) - pow(Minf, dIMFSecIntExp);
    /* Factor of log(10) because mass functions are per log10 mass
       interval. */
    dIMFSec *= mssn->sn.dFracBinSNIa * dNorm_2nd
	* imf1to8PreFactor(&mssn->ms) * dMass2*dMass2 / dIMFSecIntExp
	/pow(log(10.0), 2.0);
    return dIMFSec;
    }


#if 0
#define STARFORM_HINCLUDED
#define STARFORM
#define GASOLINE
int main ()
{
    SN sn;
    MSPARAM MSparam;
    struct inMSIMFSec mssn;
	//    MSSN mssn;
    
    MSInitialize (&MSparam);
    // crap();
    
    snInitialize (&sn);
	//    snInitConstants (sn);
    mssn.ms = *MSparam;
    mssn.sn = *sn;
    printf ("%g %g\n", mssn.ms.a1, mssn.sn.dFracBinSNIa);
    
	}
#endif

#endif
