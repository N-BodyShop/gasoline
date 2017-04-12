#ifndef COOLING_GRACKLE_HINCLUDED
#define COOLING_GRACKLE_HINCLUDED

#ifndef LOG_HINCLUDED
#include "log.h"
#endif

/* Global consts */
#include "floattype.h"
#include "param.h"
#include <sys/param.h> /* for MAXPATHLEN */
#ifndef MAXPATHLEN
#define MAXPATHLEN 256
#endif

// needed for PST
#define CL_NMAXBYTETABLE   56000

// double for variables -- must be consistent with GRACKLE compile
#define CONFIG_BFLOAT_8

#include "grackle.h"

// Default to tabular only version unless compiled in.  Max sensible value for this is 3
#ifndef GRACKLE_PRIMORDIAL_CHEMISTRY_MAX
#define GRACKLE_PRIMORDIAL_CHEMISTRY_MAX 1
#endif

typedef struct CoolingParametersStruct {
    int bDoIonOutput;

// Note many more possible parameters:  see chemistry_data.h
// Note some are probably reset by internal code
    int grackle_verbose; // verbose flag
    int use_grackle; // = 1;            // chemistry on
    int with_radiative_cooling; // = 1; // cooling on
    int primordial_chemistry; // = 3;   // molecular network with H, He, D
    int metal_cooling; // = 1;          // metal cooling on
    int UVbackground; // = 1;           // UV background on

    int bComoving; // part of units
    char grackle_data_file[MAXPATHLEN]; // "../../input/CloudyData_UVB=HM2012.h5"; // data file
} COOLPARAM;


typedef struct CoolingParticleStruct {
#if (GRACKLE_PRIMORDIAL_CHEMISTRY_MAX<1)
    float dummy;
#endif
#if (GRACKLE_PRIMORDIAL_CHEMISTRY_MAX>=1)
    gr_float HI, HII, HeI, HeII, HeIII, e;
#if (GRACKLE_PRIMORDIAL_CHEMISTRY_MAX>=2)
    gr_float HM, H2I, H2II;
#if (GRACKLE_PRIMORDIAL_CHEMISTRY_MAX>=3)
    gr_float DI, DII, HDI
#endif
#endif
#endif
} COOLPARTICLE;

/* Heating Cooling Context */

typedef struct CoolingPKDStruct { 
// not official grackle stuff    
    double     z; 
    double     a;
    double     dTime;
    double     dSecUnit; 
    double     dComovingGmPerCcUnit; 
    double     dErgPerGmUnit; 
    double     diErgPerGmUnit;
    double     dErgPerGmPerSecUnit;
// Grackle data
    char grackle_data_file[MAXPATHLEN]; // "../../input/CloudyData_UVB=HM2012.h5"; // data file
    chemistry_data *pgrackle_data;  // defined in chemistry_data.h, points at global grackle_data
    code_units my_units;     // defined in code_units.h
} COOL;

COOL *CoolInit( );
void CoolFinalize( COOL *cl );

void clInitConstants( COOL *cl, double dGMPerCcunit, double dComovingGmPerCcUnit,
		      double dErgPerGmUnit, double dSecUnit, double dKpcUnit, COOLPARAM CoolParam);
void CoolInitRatesTable( COOL *cl, COOLPARAM CoolParam);

void CoolAddParams( COOLPARAM *CoolParam, PRM );
void CoolLogParams( COOLPARAM *CoolParam, LOGGER *lgr );
void CoolOutputArray( COOLPARAM *CoolParam, int, int *, char * );

// GRACKLE_PRIMORDIAL_CHEMISTRY_MAX >=1
#define COOL_ARRAY0_EXT  "HI"
FLOAT COOL_ARRAY0(COOL *cl, COOLPARTICLE *cp, double ZMetal);
void COOL_IN_ARRAY0(COOL *cl, COOLPARTICLE *cp, double ZMetal, double Data);
#define COOL_IN_ARRAY0(w,x,y,z)  

#define COOL_ARRAY1_EXT  "HII"
FLOAT COOL_ARRAY1(COOL *cl, COOLPARTICLE *cp, double ZMetal);
void COOL_IN_ARRAY1(COOL *cl, COOLPARTICLE *cp, double ZMetal, double Data);
#define COOL_IN_ARRAY1(w,x,y,z)  

#define COOL_ARRAY2_EXT  "HeI"
FLOAT COOL_ARRAY2(COOL *cl, COOLPARTICLE *cp, double ZMetal);
void COOL_IN_ARRAY2(COOL *cl, COOLPARTICLE *cp, double ZMetal, double Data);
#define COOL_IN_ARRAY2(w,x,y,z)  

#define COOL_ARRAY3_EXT  "HeII"
FLOAT COOL_ARRAY3(COOL *cl, COOLPARTICLE *cp, double ZMetal);
void COOL_IN_ARRAY3(COOL *cl, COOLPARTICLE *cp, double ZMetal, double Data);
#define COOL_IN_ARRAY3(w,x,y,z)  

#define COOL_ARRAY4_EXT  "HeIII"
FLOAT COOL_ARRAY4(COOL *cl, COOLPARTICLE *cp, double ZMetal);
void COOL_IN_ARRAY4(COOL *cl, COOLPARTICLE *cp, double ZMetal, double Data);
#define COOL_IN_ARRAY4(w,x,y,z)  

#define COOL_ARRAY5_EXT  "e"
FLOAT COOL_ARRAY5(COOL *cl, COOLPARTICLE *cp, double ZMetal);
void COOL_IN_ARRAY5(COOL *cl, COOLPARTICLE *cp, double ZMetal, double Data);
#define COOL_IN_ARRAY5(w,x,y,z)  

// GRACKLE_PRIMORDIAL_CHEMISTRY_MAX >=2
#define COOL_ARRAY6_EXT  "HM"
FLOAT COOL_ARRAY6(COOL *cl, COOLPARTICLE *cp, double ZMetal);
void COOL_IN_ARRAY6(COOL *cl, COOLPARTICLE *cp, double ZMetal, double Data);
#define COOL_IN_ARRAY6(w,x,y,z)  

#define COOL_ARRAY7_EXT  "H2I"
FLOAT COOL_ARRAY7(COOL *cl, COOLPARTICLE *cp, double ZMetal);
void COOL_IN_ARRAY7(COOL *cl, COOLPARTICLE *cp, double ZMetal, double Data);
#define COOL_IN_ARRAY7(w,x,y,z)  

#define COOL_ARRAY8_EXT  "H2II"
FLOAT COOL_ARRAY8(COOL *cl, COOLPARTICLE *cp, double ZMetal);
void COOL_IN_ARRAY8(COOL *cl, COOLPARTICLE *cp, double ZMetal, double Data);
#define COOL_IN_ARRAY8(w,x,y,z)  

// GRACKLE_PRIMORDIAL_CHEMISTRY_MAX >=3
#define COOL_ARRAY9_EXT  "DI"
FLOAT COOL_ARRAY9(COOL *cl, COOLPARTICLE *cp, double ZMetal);
void COOL_IN_ARRAY9(COOL *cl, COOLPARTICLE *cp, double ZMetal, double Data);
#define COOL_IN_ARRAY9(w,x,y,z)  

#define COOL_ARRAY10_EXT  "DII"
FLOAT COOL_ARRAY10(COOL *cl, COOLPARTICLE *cp, double ZMetal);
void COOL_IN_ARRAY10(COOL *cl, COOLPARTICLE *cp, double ZMetal, double Data);
#define COOL_IN_ARRAY10(w,x,y,z)  

#define COOL_ARRAY11_EXT  "HDI"
FLOAT COOL_ARRAY11(COOL *cl, COOLPARTICLE *cp, double ZMetal);
void COOL_IN_ARRAY11(COOL *cl, COOLPARTICLE *cp, double ZMetal, double Data);
#define COOL_IN_ARRAY11(w,x,y,z)  

#define COOL_ARRAY12_EXT  "dummy"
#define COOL_ARRAY12(x,y,z)  0
#define COOL_IN_ARRAY12(w,x,y,z)  

#define COOL_ARRAY13_EXT  "dummy"
#define COOL_ARRAY13(x,y,z)  0
#define COOL_IN_ARRAY13(w,x,y,z)  

#define COOL_ARRAY14_EXT  "dummy"
#define COOL_ARRAY14(x,y,z)  0
#define COOL_IN_ARRAY14(w,x,y,z)  

#define COOL_ARRAY15_EXT  "dummy"
#define COOL_ARRAY15(x,y,z)  0
#define COOL_IN_ARRAY15(w,x,y,z)  

FLOAT COOL_EDOT( COOL *cl_, COOLPARTICLE *cp_, double ECode_, double rhoCode_, double ZMetal_, double *posCode_ );
#define COOL_EDOT( cl_, cp_, ECode_, rhoCode_, ZMetal_, posCode_) (CoolCodeWorkToErgPerGmPerSec( cl_, CoolEdotInstantCode( cl_, cp_, ECode_, rhoCode_, ZMetal_, posCode_ )))

FLOAT COOL_COOLING( COOL *cl_, COOLPARTICLE *cp_, double ECode_, double rhoCode_, double ZMetal_, double *posCode_ );
#define COOL_COOLING( cl_, cp_, ECode_, rhoCode_, ZMetal_, posCode_) (CoolCodeWorkToErgPerGmPerSec( cl_, CoolCoolingCode( cl_, cp_, ECode_, rhoCode_, ZMetal_, posCode_ )))

FLOAT COOL_HEATING( COOL *cl_, COOLPARTICLE *cp_, double ECode_, double rhoCode_, double ZMetal_, double *posCode_ );
#define COOL_HEATING( cl_, cp_, ECode_, rhoCode_, ZMetal_, posCode_) (CoolCodeWorkToErgPerGmPerSec( cl_, CoolHeatingCode( cl_, cp_, ECode_, rhoCode_, ZMetal_, posCode_ ))) 

//void CoolPARTICLEtoPERBARYON(COOL *cl_, PERBARYON *Y, COOLPARTICLE *cp, double ZMetal);
//void CoolPERBARYONtoPARTICLE(COOL *cl_, PERBARYON *Y, COOLPARTICLE *cp, double ZMetal);

double CoolEnergyToTemperature( COOL *Cool, COOLPARTICLE *cp, double E, double, double ZMetal);
double CoolCodeEnergyToTemperature( COOL *Cool, COOLPARTICLE *cp, double E, double, double ZMetal);

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
			 double ExternalHeating, double rho, double ZMetal, double *rp,  double tStep );

void CoolIntegrateEnergyCode(COOL *cl, COOLPARTICLE *cp, double *E, 
			     double ExternalHeating, double rho, double ZMetal, double *r, double tStep );

void CoolDefaultParticleData( COOLPARTICLE *cp );

void CoolInitEnergyAndParticleData( COOL *cl, COOLPARTICLE *cp, double *E, double dDensity, double dTemp, double ZMetal);

/* Deprecated */
double CoolHeatingRate( COOL *cl, COOLPARTICLE *cp, double E, double dDensity, double ZMetal, double rkpc);

double CoolEdotInstantCode(COOL *cl, COOLPARTICLE *cp, double ECode, 
			   double rhoCode, double ZMetal, double *posCode );
double CoolCoolingCode(COOL *cl, COOLPARTICLE *cp, double ECode, 
		       double rhoCode, double ZMetal, double *posCode );
double CoolHeatingCode(COOL *cl, COOLPARTICLE *cp, double ECode, 
		       double rhoCode, double ZMetal, double *posCode );

void CoolCodePressureOnDensitySoundSpeed( COOL *cl, COOLPARTICLE *cp, double uPred, double fDensity, double gamma, double gammam1, double *PoverRho, double *c );

/* Note: gamma should be 5/3 for this to be consistent! */
#define CoolCodePressureOnDensitySoundSpeed( cl__, cp__, uPred__, fDensity__, gamma__, gammam1__, PoverRho__, c__ ) { \
  *(PoverRho__) = ((5./3.-1)*(uPred__)); \
  *(c__) = sqrt((5./3.)*(*(PoverRho__))); }

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
