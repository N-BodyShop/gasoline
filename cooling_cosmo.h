#ifndef COOLING_COSMO_HINCLUDED
#define COOLING_COSMO_HINCLUDED

/* Global consts */
#if defined(COOLDEBUG) || defined(STARFORM)
#include "mdl.h"
#endif
#include "floattype.h"
#include "param.h"

/* Constants */
#define CL_B_gm         (6.022e23*(938.7830/931.494))
#define CL_k_Boltzmann  1.38066e-16
#define CL_eV_erg       1.60219e-12
#define CL_eV_per_K     (CL_k_Boltzmann/CL_eV_erg)
/*
#define CL_RT_FLOAT      float
#define CL_RT_MIN        1e-38
#define CL_RT_MIN        FLT_MIN
*/

#define CL_RT_FLOAT      double
#define CL_RT_MIN        1e-100

/*
#define CL_RT_MIN        DBL_MIN
*/
/* 
 * Work around for Dec ev6 flawed
 * treatment of sub-normal numbers 
 */
#define CL_MAX_NEG_EXP_ARG  -500.

#define CL_NMAXBYTETABLE   56000


typedef struct CoolingParametersStruct {
	int    bIonNonEqm;
	int    nCoolingTable;
	int    bUV;
	int    bUVTableUsesTime;
	int    bDoIonOutput;
	double dMassFracHelium;
	double dCoolingTmin;
	double dCoolingTmax;
	} COOLPARAM;

typedef struct CoolingParticleStruct {
	FLOAT Y_HI,Y_HeI,Y_HeII;	/* Abundance of ions */
	} COOLPARTICLE;

typedef struct { 
  double e,Total;
  double HI,HII,HeI,HeII,HeIII;
} PERBARYON;

typedef struct { 
  double   zTime;

  double   Rate_Phot_HI;
  double   Rate_Phot_HeI;
  double   Rate_Phot_HeII;

  double   Heat_Phot_HI;
  double   Heat_Phot_HeI;
  double   Heat_Phot_HeII;
} UVSPECTRUM;

typedef struct { 
  double   Rate_Phot_HI;
  double   Rate_Phot_HeI;
  double   Rate_Phot_HeII;

  double   Heat_Phot_HI;
  double   Heat_Phot_HeI;
  double   Heat_Phot_HeII;
 
  double   Cool_Coll_HI;
  double   Cool_Coll_HeI;
  double   Cool_Coll_HeII;
  double   Cool_Diel_HeII;
  double   Cool_Comp;

} RATES_NO_T;

typedef struct { 
  CL_RT_FLOAT   Rate_Coll_HI;
  CL_RT_FLOAT   Rate_Coll_HeI;
  CL_RT_FLOAT   Rate_Coll_HeII;
  CL_RT_FLOAT   Rate_Radr_HII;
  CL_RT_FLOAT   Rate_Radr_HeII;
  CL_RT_FLOAT   Rate_Radr_HeIII;
  CL_RT_FLOAT   Rate_Diel_HeII;

  CL_RT_FLOAT   Cool_Brem_1;
  CL_RT_FLOAT   Cool_Brem_2;
  CL_RT_FLOAT   Cool_Radr_HII;
  CL_RT_FLOAT   Cool_Radr_HeII;
  CL_RT_FLOAT   Cool_Radr_HeIII;
  CL_RT_FLOAT   Cool_Line_HI;
  CL_RT_FLOAT   Cool_Line_HeI;
  CL_RT_FLOAT   Cool_Line_HeII;

} RATES_T;

/* Heating Cooling Context */

typedef struct CoolingPKDStruct { 
   double     z; /* Redshift */
   double     dTime;
 /* Rates independent of Temperature */ 
   RATES_NO_T  R;
 /* Table for Temperature dependent rates */ 
   int        nTable;
   double     TMin;
   double     TMax;
   double     TlnMin;
   double     TlnMax;
   double     rDeltaTln;
   RATES_T     *RT;

   int        nTableRead; /* number of Tables read from files */

   int        bUV;
   int        nUV;
   UVSPECTRUM *UV;
   int        bUVTableUsesTime;
   int        bUVTableLinear;

   double     dGmPerCcUnit;
   double     dComovingGmPerCcUnit;
   double     dErgPerGmUnit;
   double     dSecUnit;
   double     dErgPerGmPerSecUnit;
   double     diErgPerGmUnit;
   double     dKpcUnit;
   double     Y_H;
   double     Y_He;
   double     Y_eMAX;
   void       *DerivsData;

#if defined(COOLDEBUG) || defined(STARFORM)
   MDL        mdl; /* For diag/debug outputs */
   struct particle *p; /* particle pointer needed for SN feedback */
#endif
   
} COOL;

typedef struct {
  double   T, Tln;
  double   Coll_HI;
  double   Coll_HeI;
  double   Coll_HeII;
  double   Radr_HII;
  double   Radr_HeII;
  double   Diel_HeII;
  double   Totr_HeII;
  double   Radr_HeIII;
} RATE;

typedef struct {
  double compton;
  double bremHII;
  double bremHeII;
  double bremHeIII;
  double radrecHII;
  double radrecHeII;
  double radrecHeIII;
  double collionHI; 
  double collionHeI;
  double collionHeII;
  double dielrecHeII;
  double lineHI;
  double lineHeI;
  double lineHeII;
} COOL_ERGPERSPERGM;


typedef struct clDerivsDataStruct {
  void *IntegratorContext;
  COOL *cl;
  double rho,PdV,E;
  RATE Rate;
  PERBARYON Y;
  double     Y_Total0, Y_Total1;
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

void clRatesTableError( COOL *cl );
void clRatesRedshift( COOL *cl, double z, double dTime );
double clHeatTotal ( COOL *cl, PERBARYON *Y );
void clRates ( COOL *cl, RATE *Rate, double T );
double clCoolTotal ( COOL *cl, PERBARYON *Y, RATE *Rate, double rho );
COOL_ERGPERSPERGM  clTestCool ( COOL *cl, PERBARYON *Y, RATE *Rate, double rho );
void clPrintCool ( COOL *cl, PERBARYON *Y, RATE *Rate, double rho );

void clAbunds( COOL *cl, PERBARYON *Y, RATE *Rate, double rho);
double clThermalEnergy( double Y_Total, double T );
double clTemperature( double Y_Total, double E );
double clRateCollHI( double T );
double clRateCollHeI( double T );
double clRateCollHeII( double T );
double clRateRadrHII( double T );
double clRateRadrHeII( double T );
double clRateDielHeII( double T );
double clRateRadrHeIII( double T );
double clCoolBrem1( double T );
double clCoolBrem2( double T );
double clCoolRadrHII( double T );
double clCoolRadrHeII( double T );
double clCoolRadrHeIII( double T );
double clCoolLineHI( double T );
double clCoolLineHeI( double T );
double clCoolLineHeII( double T );
double clEdotInstant ( COOL *cl, PERBARYON *Y, RATE *Rate, double rho );
void clIntegrateEnergy(COOL *cl, PERBARYON *Y, double *E, 
		       double PdV, double rho, double dt );
void clIntegrateEnergyDEBUG(COOL *cl, PERBARYON *Y, double *E, 
		       double PdV, double rho, double dt );


void clDerivs(void *Data, double x, double *y, double *dydx) ;

void clJacobn(void *Data, double x, double y[], double dfdx[], double **dfdy) ;
  
void CoolAddParams( COOLPARAM *CoolParam, PRM );
void CoolLogParams( COOLPARAM *CoolParam, FILE *fp );
void CoolOutputArray( COOLPARAM *CoolParam, int, int *, char * );

FLOAT COOL_ARRAY0(COOLPARTICLE *cp);
#define COOL_ARRAY0( cp ) ((cp)->Y_HI)

FLOAT COOL_ARRAY1(COOLPARTICLE *cp);
#define COOL_ARRAY1( cp ) ((cp)->Y_HeI)

FLOAT COOL_ARRAY2(COOLPARTICLE *cp);
#define COOL_ARRAY2( cp ) ((cp)->Y_HeII)

void CoolPARTICLEtoPERBARYON(PERBARYON *Y, COOLPARTICLE *cp, double HTotal, double HeTotal);

#define CoolPARTICLEtoPERBARYON(Y, cp, HTotal, HeTotal) { \
    (Y)->HI = (cp)->Y_HI; \
    (Y)->HII = HTotal - (Y)->HI; \
    (Y)->HeI = (cp)->Y_HeI; \
    (Y)->HeII = (cp)->Y_HeII; \
    (Y)->HeIII = HeTotal - (Y)->HeI - (Y)->HeII; \
    (Y)->e = (Y)->HII + (Y)->HeII + 2*(Y)->HeIII; \
    (Y)->Total = (Y)->e + HTotal + HeTotal; }

void CoolPERBARYONtoPARTICLE(PERBARYON *Y, COOLPARTICLE *cp);

#define CoolPERBARYONtoPARTICLE(Y, cp) { \
    (cp)->Y_HI = (Y)->HI; \
    (cp)->Y_HeI = (Y)->HeI; \
    (cp)->Y_HeII = (Y)->HeII; }


double CoolCodeEnergyToTemperature( COOL *Cool, COOLPARTICLE *cp, double E );

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

double CodeDensityToComovingGmPerCc( COOL *Cool, double dCodeDensity );

#define CodeDensityToComovingGmPerCc( Cool, dCodeDensity )  ((Cool)->dComovingGmPerCcUnit*(dCodeDensity))

void CoolIntegrateEnergy(COOL *cl, COOLPARTICLE *cp, double *E, 
		       double PdV, double rho, double tStep );

void CoolIntegrateEnergyEPDRCode(COOL *cl, COOLPARTICLE *cp, double *E, 
		       double PdV, double rho, double *r, double tStep );

void CoolDefaultParticleData( COOLPARTICLE *cp );

void CoolInitEnergyAndParticleData( COOL *cl, COOLPARTICLE *cp, double *E, double dDensity, double dTemp );

double CoolHeatingRate( COOL *cl, COOLPARTICLE *cp, double E, double dDensity );

double CoolCodePressureOnDensity( COOL *cl, COOLPARTICLE *cp, double uPred, double fDensity, double gammam1 );

#define CoolCodePressureOnDensity( cl, cp, uPred, fDensity, gammam1 ) ((gammam1)*(uPred))

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


