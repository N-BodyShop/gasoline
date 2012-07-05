#ifndef COOLING_COSMO_HINCLUDED
#define COOLING_COSMO_HINCLUDED

/* Global consts */
#if defined(COOLDEBUG) || defined(STARFORM)
#include "mdl.h"
#else
#endif
#include "floattype.h"
#include "param.h"

/* Constants */
#define CL_B_gm         (6.022e23*(938.7830/931.494))
#define CL_k_Boltzmann  1.38066e-16
#define CL_eV_erg       1.60219e-12
#define CL_eV_per_K     (CL_k_Boltzmann/CL_eV_erg)
/* 
 * Work around for Dec ev6 flawed
 * treatment of sub-normal numbers 
 */
#define CL_MAX_NEG_EXP_ARG  -500.

#define CL_NMAXBYTETABLE   56000


typedef struct CoolingParametersStruct {
	double    BaseT;
	double    dParam2;
	double    dParam3;
	double    dParam4;
	double    Y_Total;
	double dCoolingTmin;
	double dCoolingTmax;
	} COOLPARAM;

/* How to do a dummy? --> can't I guess */
typedef struct CoolingParticleStruct {
	FLOAT Y_Total;
	} COOLPARTICLE;

/* Heating Cooling Context */

typedef struct CoolingPKDStruct { 
   double     z; /* Redshift */
   double     dTime;

   double     dGmPerCcUnit;
   double     dComovingGmPerCcUnit;
   double     dErgPerGmUnit;
   double     dSecUnit;
   double     dErgPerGmPerSecUnit;
   double     diErgPerGmUnit;
   double     dKpcUnit;

   double    BaseT;
   double    dParam2;
   double    dParam3;
   double    dParam4;
   double    Y_Total;
   double    Tmin;
   double    Tmax;
   void       *DerivsData;

   int nTableRead; /* Internal Tables read from Files */

#if defined(COOLDEBUG) || defined(STARFORM)
   MDL        mdl; /* For diag/debug outputs */
   struct particle *p; /* particle pointer needed for SN feedback */
#endif
   
} COOL;

typedef struct clDerivsDataStruct {
  void *IntegratorContext;
  COOL *cl;
  double rho,PdV,E,T,Y_Total,rFactor;
  double     dlnE;
  int        its;  /* Debug */
} clDerivsData;


COOL *CoolInit( );
void CoolFinalize( COOL *cl );

void clInitConstants( COOL *cl, double dGMPerCcunit, double dComovingGmPerCcUnit,
					 double dErgPerGmUnit, double dSecUnit, double dKpcUnit, COOLPARAM CoolParam);
void clInitUV(COOL *cl, int nTableColumns, int nTableRows, double *dTableData );
void clInitRatesTable( COOL *cl, double TMin, double TMax, int nTable );
void CoolInitRatesTable( COOL *cl, COOLPARAM CoolParam);

double clThermalEnergy( double Y_Total, double T );
double clTemperature( double Y_Total, double E );

double clEdotInstant( COOL *cl, double E, double T, double rho, double r );
void clIntegrateEnergy(COOL *cl, double *E, 
		       double PdV, double rho, double Y_Total, double radius, double tStep );


void clDerivs(void *Data, double x, double *y, double *dydx) ;

void clJacobn(void *Data, double x, double y[], double dfdx[], double **dfdy) ;
  
void CoolAddParams( COOLPARAM *CoolParam, PRM );
void CoolLogParams( COOLPARAM *CoolParam, FILE *fp );
void CoolOutputArray( COOLPARAM *CoolParam, int, int *, char * );

#define COOL_ARRAY0_EXT "Y"
FLOAT COOL_ARRAY0(COOL *cl,COOLPARTICLE *cp, float *Z);
#define COOL_ARRAY0( cl_, cp, Z_ ) ((cp)->Y_Total)

#define COOL_ARRAY1_EXT "junk"
FLOAT COOL_ARRAY1(COOL *cl,COOLPARTICLE *cp, float *Z);
#define COOL_ARRAY1( cl_, cp, Z_ ) (0)

#define COOL_ARRAY2_EXT "junk2"
FLOAT COOL_ARRAY2(COOL *cl,COOLPARTICLE *cp, float *Z);
#define COOL_ARRAY2( cl_, cp, Z_ ) (0)

FLOAT COOL_EDOT( COOL *cl_, COOLPARTICLE *cp_, double ECode_, double rhoCode_, double ZMetal_, double *posCode_ );
#define COOL_EDOT( cl_, cp_, ECode_, rhoCode_, ZMetal_, posCode_) (CoolCodeWorkToErgPerGmPerSec( cl_, CoolEdotInstantCode( cl_, cp_, ECode_, rhoCode_, ZMetal_, posCode_ )))

FLOAT COOL_COOLING( COOL *cl_, COOLPARTICLE *cp_, double ECode_, double rhoCode_, double ZMetal_, double *posCode_ );
#define COOL_COOLING( cl_, cp_, ECode_, rhoCode_, ZMetal_, posCode_) (CoolCodeWorkToErgPerGmPerSec( cl_, CoolEdotInstantCode( cl_, cp_, ECode_, rhoCode_, ZMetal_, posCode_ )))

FLOAT COOL_HEATING( COOL *cl_, COOLPARTICLE *cp_, double ECode_, double rhoCode_, double ZMetal_, double *posCode_ );
#define COOL_HEATING( cl_, cp_, ECode_, rhoCode_, ZMetal_, posCode_) (0)

double CoolCodeEnergyToTemperature( COOL *Cool, COOLPARTICLE *cp, double E, double ZMetal );

/* Note: nod to cosmology (z parameter) unavoidable unless we want to access cosmo.[ch] from here */
void CoolSetTime( COOL *Cool, double dTime, double z );

double CoolCodeTimeToSeconds( COOL *Cool, double dCodeTime );

#define CoolCodeTimeToSeconds( Cool, dCodeTime ) ((Cool)->dSecUnit*(dCodeTime))

double CoolSecondsToCodeTime( COOL *Cool, double dTime );

#define CoolSecondsToCodeTime( Cool, dTime ) ((dTime)/(Cool)->dSecUnit)

double CoolCodeEnergyToErgPerGm( COOL *Cool, double dCodeEnergy );

#define CoolCodeEnergyToErgPerGm( Cool, dCodeEnergy ) ((Cool)->dErgPerGmUnit*(dCodeEnergy))

double CoolErgPerGmToCodeEnergy( COOL *Cool, double dEnergy );

#define CoolErgPerGmToCodeEnergy( Cool, dEnergy ) ((Cool)->diErgPerGmUnit*(dEnergy))

double CoolCodeWorkToErgPerGmPerSec( COOL *Cool, double dCodeWork );

#define CoolCodeWorkToErgPerGmPerSec( Cool, dCodeWork ) ((Cool)->dErgPerGmPerSecUnit*(dCodeWork))

double CoolErgPerGmPerSecToCodeWork( COOL *Cool, double dWork );

#define CoolErgPerGmPerSecToCodeWork( Cool, dWork ) ((dWork)/(Cool)->dErgPerGmPerSecUnit)

double CodeDensityToComovingGmPerCc( COOL *Cool, double dCodeDensity );

#define CodeDensityToComovingGmPerCc( Cool, dCodeDensity )  ((Cool)->dComovingGmPerCcUnit*(dCodeDensity))

void CoolIntegrateEnergy(COOL *cl, COOLPARTICLE *cp, double *E, 
		       double PdV, double rho, double ZMetal, double tStep );

void CoolIntegrateEnergyCode(COOL *cl, COOLPARTICLE *cp, double *E, 
		       double PdV, double rho, double ZMetal, double *r, double tStep );

void CoolDefaultParticleData( COOLPARTICLE *cp );

void CoolInitEnergyAndParticleData( COOL *cl, COOLPARTICLE *cp, double *E, double dDensity, double dTemp, double ZMetal );

/* Deprecated */
double CoolHeatingRate( COOL *cl, COOLPARTICLE *cp, double E, double dDensity );

double CoolEdotInstantCode(COOL *cl, COOLPARTICLE *cp, double ECode, 
			   double rhoCode, double ZMetal, double *posCode );

void CoolCodePressureOnDensitySoundSpeed( COOL *cl, COOLPARTICLE *cp, double uPred, double fDensity, double gamma, double gammam1, double *PoverRho, double *c );

/* Taken from Bate, ApJ, 508, L95 (1998) 
   rho = powerlaw with 4 gamma values 
   PoverRho has units of ergs per gm (same as speed squared)
*/

#ifdef MODBATEPOLY
/* close to "opacity limit" of Bate et al. -- corrected up to actual T when used */
//#define RHOMIN 1e-16
//#define RHOMIN 3.5e-19
#define RHOMIN 1.5e-15
  /* mu=2.33, T=1 K  P/rho = kT/(mu*mH) */
#define PONRHOMIN (1.38066e-16/(2.33*1.67e-24)*1)
#define GetGammaEff(rho__,gammaEff__, PonRho__) { \
  if ( rho__ <= RHOMIN ) { \
     gammaEff__ = 1.0; \
     PonRho__  = PONRHOMIN; \
  } else if ( rho__ <= RHOMIN*10 ) { \
     double x_ = log10(rho__/RHOMIN); \
     if (x_ < 0.5) { \
       gammaEff__ = 1+2*x_*x_; \
       PonRho__  = PONRHOMIN*pow(10.,2/3.*x_*x_*x_); \
       } \
     else { \
       gammaEff__ = 4*x_-2*x_*x_; \
       PonRho__  = PONRHOMIN*pow(10.,2*x_*x_-2./3.*x_*x_*x_+4./24.-x_); \
       } \
  } else { \
     gammaEff__ = 2; \
     PonRho__  = rho__*PONRHOMIN/RHOMIN*0.316228; \
  } \
}

#define CoolCodePressureOnDensitySoundSpeed( cl__, cp__, uPred__, fDensity__, gamma__, gammam1__, PoverRho__, c__ ) { \
  double rho__ = (fDensity__)*(cl__)->dGmPerCcUnit; \
  double gammaEff__, PonRhoCGS__; \
  GetGammaEff(rho__,gammaEff__, PonRhoCGS__); \
  PonRhoCGS__ *= cl->BaseT; \
  *(PoverRho__) = PonRhoCGS__*(cl__)->diErgPerGmUnit; \
  *(c__) = sqrt((gammaEff__)*(*(PoverRho__))); } 

#else 
#define GetGammaEff(rho__,gammaEff__, KEff__) { \
  double coeff__ = 0.1*18321.359*18321.359; /* mu=2.46, T=1 K, gamma=7/5   2.0e4*2.0e4; */ \
  if ( rho__ <= 1e-13 ) { \
     gammaEff__ = 1.0; \
     KEff__  = coeff__; \
  } else if ( rho__ <= 5.7e-8 ) { \
     gammaEff__ = 1.4; \
     KEff__  = coeff__*158489.; \
  } else if ( rho__ <= 1e-3 ) { \
     gammaEff__ = 1.15; \
     KEff__  = coeff__*158489.*0.0154514; \
  } else { \
     gammaEff__ = 5./3.; \
     KEff__  = coeff__*158489.*0.0154514*35.4813; \
  } \
}

#define CoolCodePressureOnDensitySoundSpeed( cl__, cp__, uPred__, fDensity__, gamma__, gammam1__, PoverRho__, c__ ) { \
  double rho__ = (fDensity__)*(cl__)->dGmPerCcUnit; \
  double gammaEff__, KEff__; \
  GetGammaEff(rho__,gammaEff__, KEff__); \
  KEff__ *= cl->BaseT; \
  *(PoverRho__) = KEff__*pow(rho__,gammaEff__-1.)*(cl__)->diErgPerGmUnit; \
  *(c__) = sqrt((gammaEff__)*(*(PoverRho__))); } 

#endif
/*
double CoolCodePressureOnDensity( COOL *cl, COOLPARTICLE *cp, double uPred, double fDensity, double gammam1 );

#define CoolCodePressureOnDensity( cl, cp, uPred, fDensity, gammam1 ) ((gammam1)*(uPred))
*/

struct inInitCooling {
    double dGmPerCcUnit;
    double dComovingGmPerCcUnit;
    double dErgPerGmUnit;
    double dSecUnit;
	double dKpcUnit;
    double z;
	double dTime;
	COOLPARAM CoolParam;
    };

struct inInitEnergy {
	double dTuFac;
	double z;
	double dTime;
	};

void CoolTableReadInfo( COOLPARAM *CoolParam, int cntTable, int *nTableColumns, char *suffix );

void CoolTableRead( COOL *Cool, int nData, void *vData);

#endif


