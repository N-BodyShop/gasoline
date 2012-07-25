#ifndef MILLERSCALO_HINCLUDED
#define MILLERSCALO_HINCLUDED
#ifdef CHABRIER
#include "romberg.h"
/*
  Use the log normal + power law fit of Chabrier, 2003, Galactic
  Stellar and Substellar Initial Mass Function", PASP 115, 763.
*/

typedef struct MillerScaloContext
{
    /*
      Chabrier low mass formula:
       \xi(log m) = A exp [ - (log m - log m_c)^2/2 \sigma^2]
     */
    double a1, sigma, mc;
    /* For high mass: normalization, index, minimum mass */
    double a2, b2, m2;
    double mmax;
    } * MSPARAM;
#else
/*
    Uses the 3 segment power law fit for the Miller-Scalo IMF
    (Ap.J. Supp., 41,1979).

                                a1*(M**(b1))          0.1<M<1.
              IMF(Log10(M))=    a2*(M**(b2))         1.<M<10.
                                a3*(M**(b3))         10.<M
*/

typedef struct MillerScaloContext
{
    /* normalization, index, minimum mass */
    double a1, b1, m1;
    double a2, b2, m2;
    double a3, b3, m3;
    double mmax;
    } * MSPARAM;
#endif


void MSInitialize(MSPARAM *pms);

double dMSIMF(MSPARAM p, double mass);
double dMSCumNumber(MSPARAM p, double mass);
double dMSCumMass(MSPARAM p, double mass);
double imfCumLuminosity(MSPARAM p, double mass);
double imf1to8Exp(MSPARAM p);
double imf1to8PreFactor(MSPARAM p);

#endif
