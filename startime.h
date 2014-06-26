
#ifndef STARTIME_HINCLUDED
#define STARTIME_HINCLUDED
/*
  Uses Raiteri, Villata, and Navarro (A&A, 315, 105, 1996) fit to 
  Padova group stellar models to find stellar lifetimes as a 
  function of mass, M,and metallicity, Z.

  log t = a0(Z) + a1(Z)*log(M) + a2(Z)*(log(M))**2

  where t is stellar lifetime in years, M is stellar mass 
  in solar masses, and Z, defined as the total mass fracion in
  elements heavier than He, is in the range 0.0004-0.05.  
  Coefficients are given as follows:

  a0(Z) = 10.13 + 0.07547*log(Z) - 0.008084*(log(Z))**2

  a1(Z) = -4.424 + 0.7939*log(Z) - 0.1187*(log(Z))**2

  a2(Z) = 1.262 + 0.3385*log(Z) - 0.1187*(log(Z))**2

  zmin, zmax = minimum and maximum metallicity in stellar lifetimes equation 
  zsol = solar metal abundance 
  xoxsol = solar oxygen abundance */

typedef struct PadovaContext
{
    double a00, a01, a02;
    double a10, a11, a12;
    double a20, a21, a22;
    double a0, a1, a2;    
    double zmin, zmax, zsol, xoxsol;    
    } * PDVAPARAM;


void PadovaInitialize(PDVAPARAM *ppdva);
void PadovaCoefInit (PDVAPARAM ppdva, double dMetals);
double dSTLtimeMStar(PDVAPARAM ppdva, double dStarMass, double dMetals);
double dSTMStarLtime(PDVAPARAM ppdva, double dStarLtime, double dMetals);
#endif
