#ifndef MILLERSCALO_HINCLUDED
#define MILLERSCALO_HINCLUDED
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


void MSInitialize(MSPARAM *pms);

double dMSIMF(MSPARAM p, double mass);
double dMSCumNumber(MSPARAM p, double mass);
double dMSCumMass(MSPARAM p, double mass);
double imf1to8Exp(MSPARAM p);
double imf1to8PreFactor(MSPARAM p);

#endif
