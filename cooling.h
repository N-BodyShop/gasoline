/* Global consts */
#include "mdl.h"

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


typedef struct { 
  double e,Total;
  double HI,HII,HeI,HeII,HeIII;
} PERBARYON;

typedef struct { 
  double   z;

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

typedef struct { 
   double     z; /* Redshift */
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

   int        nUV;
   UVSPECTRUM *UV;

   double     dGmPerCcUnit;
   double     dComovingGmPerCcUnit;
   double     dErgPerGmUnit;
   double     dSecUnit;
   double     dErgPerGmPerSecUnit;
   double     diErgPerGmUnit;
   double     Y_H;
   double     Y_He;
   double     Y_eMAX;

   MDL        mdl; /* For diag/debug outputs */
   struct particle *p; /* particle pointer needed for SN feedback */
   
} CL;

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


void clInitConstants(CL *cl, double dGMPerCcunit, double dComovingGmPerCcUnit,
		     double dErgPerGmUnit, double dSecUnit, double dMassFracHelium);
void clInitUV( CL *cl, UVSPECTRUM *UV, int nUV);
void clInitRatesTable( CL *cl, double TMin, double TMax, int nTable );
void clRatesTableError( CL *cl );
void clRatesRedshift( CL *cl, double z );
double clHeatTotal ( CL *cl, PERBARYON *Y );
void clRates ( CL *cl, RATE *Rate, double T );
double clCoolTotal ( CL *cl, PERBARYON *Y, RATE *Rate, double rho );
COOL_ERGPERSPERGM  clTestCool ( CL *cl, PERBARYON *Y, RATE *Rate, double rho );
void clPrintCool ( CL *cl, PERBARYON *Y, RATE *Rate, double rho );

void clAbunds( CL *cl, PERBARYON *Y, RATE *Rate, double rho);
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
double clEdot( CL *cl, PERBARYON *Y, RATE *Rate, double rho, 
	       PERBARYON *Y1, PERBARYON *Y2, double dt);
double clEdot_new( CL *cl, RATE *R1, RATE *R2, double rho, 
		           PERBARYON *Y1, PERBARYON *Y2, double dt);
double clEdotHarmonic( CL *cl, RATE *R1, RATE *R2, PERBARYON *Y1, PERBARYON *Y2, 
		       double rho, double PdV, double dt);
void clIntegrateEnergy(CL *cl, PERBARYON *Y, double *E, 
		       double PdV, double rho, double dt );
void clIntegrateEnergyDEBUG(CL *cl, PERBARYON *Y, double *E, 
		       double PdV, double rho, double dt );


  



