#ifdef GASOLINE
#ifndef NOCOOLING

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

/*
#include "cooling.h"
*/

/* debug */
#include "pkd.h"

#define EPS 1e-5
#define MAXABUNDITERATIONS 20
#define NDOT_MIN  1e-60
  
#define CL_eHI     (13.60*CL_eV_erg)
#define CL_eHeI    (24.59*CL_eV_erg)
#define CL_eHeII   (54.42*CL_eV_erg)
#define CL_E2HeII  (3.0*13.6*CL_eV_erg)

void clInitConstants(CL *cl, double dGmPerCcUnit, double dComovingGmPerCcUnit, 
		double dErgPerGmUnit, double dSecUnit, double dMassFracHelium)
{
    cl->dGmPerCcUnit = dGmPerCcUnit;
    cl->dComovingGmPerCcUnit = dComovingGmPerCcUnit;
    cl->dErgPerGmUnit = dErgPerGmUnit;
    cl->dSecUnit = dSecUnit;
    cl->dErgPerGmPerSecUnit = cl->dErgPerGmUnit / cl->dSecUnit;
    cl->diErgPerGmUnit = 1./dErgPerGmUnit;

    cl->nUV = 0;
    cl->UV = NULL;

    cl->Y_H = 1.0-dMassFracHelium;
    cl->Y_He = dMassFracHelium/4;
    cl->Y_eMAX = cl->Y_H + cl->Y_He*2;
    }

void clInitUV(CL *cl, UVSPECTRUM *UV, int nUV)
{
    int i;

    cl->nUV = nUV;
    cl->UV = (UVSPECTRUM *) malloc(nUV*sizeof(UVSPECTRUM));
    
    for (i=0;i<nUV;i++) (cl->UV)[i] = UV[i];
    }
   
void clInitRatesTable( CL *cl, double TMin, double TMax, int nTable ) {
/*-----------------------------------------------------------------
 *      Function of ln(T) tables:
 * cl assumed to be predefined (allocated)
 * A table spacing of 1.23e-3 in log(T) eg. nTable=15001 10->1e9 K
 * Maximum 1% errors except at discontinuities in the functions.
 * Storage for such a table is (15001*15*4|8b) => 0.9|1.7 Mb
 *-----------------------------------------------------------------*/
  int i;
  double DeltaTln, Tln, T;

  cl->R.Cool_Coll_HI = CL_eHI*CL_B_gm;
  cl->R.Cool_Coll_HeI = CL_eHeI*CL_B_gm;
  cl->R.Cool_Coll_HeII = CL_eHeII*CL_B_gm;
  cl->R.Cool_Diel_HeII = (CL_E2HeII+CL_eHeI)*CL_B_gm;

  cl->nTable = nTable;
  cl->TMin = TMin;
  cl->TMax = TMax;
  cl->TlnMin = log( TMin );
  cl->TlnMax = log( TMax );
  DeltaTln = ( cl->TlnMax-cl->TlnMin )/( nTable - 1 );
  cl->rDeltaTln = 1./DeltaTln;
  cl->RT = (RATES_T *) malloc( nTable * sizeof(RATES_T) );

  for ( i=0; i<nTable; i++ ) {
    Tln = cl->TlnMin + DeltaTln*(i-1);
    T=exp(Tln);

    (cl->RT+i)->Rate_Coll_HI = clRateCollHI( T );
    if ( (cl->RT+i)->Rate_Coll_HI < CL_RT_MIN ) (cl->RT+i)->Rate_Coll_HI = CL_RT_MIN;
    (cl->RT+i)->Rate_Coll_HeI = clRateCollHeI( T );
    if ( (cl->RT+i)->Rate_Coll_HeI < CL_RT_MIN ) (cl->RT+i)->Rate_Coll_HeI = CL_RT_MIN;
    (cl->RT+i)->Rate_Coll_HeII = clRateCollHeII( T );
    if ( (cl->RT+i)->Rate_Coll_HeII < CL_RT_MIN ) (cl->RT+i)->Rate_Coll_HeII = CL_RT_MIN;

    (cl->RT+i)->Rate_Radr_HII = clRateRadrHII( T );
    (cl->RT+i)->Rate_Radr_HeII = clRateRadrHeII( T );
    (cl->RT+i)->Rate_Radr_HeIII = clRateRadrHeIII( T );
    (cl->RT+i)->Rate_Diel_HeII = clRateDielHeII( T );

    (cl->RT+i)->Cool_Brem_1 = clCoolBrem1( T );
    (cl->RT+i)->Cool_Brem_2 = clCoolBrem2( T );
    (cl->RT+i)->Cool_Radr_HII = clCoolRadrHII( T );
    (cl->RT+i)->Cool_Radr_HeII = clCoolRadrHeII( T );
    (cl->RT+i)->Cool_Radr_HeIII = clCoolRadrHeIII( T );
    (cl->RT+i)->Cool_Line_HI = clCoolLineHI( T );
    (cl->RT+i)->Cool_Line_HeI = clCoolLineHeII( T );
    (cl->RT+i)->Cool_Line_HeII = clCoolLineHeII( T );
  }    
}

void clRatesTableError( CL *cl ) {
  /* This estimates the error for a table of half the size of
   * the current one.
   * The collisional, dielectric and line cooling rates can all go to
   * zero (minimum) and will even for double precision (exp(-1e5/T)) 
   * the critical values that must never get down to zero
   * are the radiative recombination rates: Check this before using
   * CL_RT_FLOAT float for the table.
   */
  int i,j,ierr[15];
  double min[15],max[15];
  double err,maxerr[15];
  CL_RT_FLOAT *p,*p1,*p2;

  for (j=0;j<15;j++) {
    maxerr[j]=0.0;
    min[j]=1e300;
    max[j]=-1e300;
  }
  for ( i=1; i<cl->nTable-1; i+=2 ) {
    p = (CL_RT_FLOAT *) &((cl->RT+i)->Rate_Coll_HI);
    p1 = (CL_RT_FLOAT *) &((cl->RT+i-1)->Rate_Coll_HI);
    p2 = (CL_RT_FLOAT *) &((cl->RT+i+1)->Rate_Coll_HI);
    if (i==10001) {
      printf(" Comp %i %i %e %e\n",i,(int) sizeof(CL_RT_FLOAT),p[1],(cl->RT+i)->Rate_Coll_HeI);
    }
    for (j=0;j<15;j++) {
      if (p[j] < min[j]) min[j]=p[j];
      if (p[j] > max[j]) max[j]=p[j];
      err= fabs((0.5*(p1[j]+p2[j])-p[j])/(p[j]+1e-30));
      if (err > maxerr[j] ) {
	maxerr[j]=err;
	ierr[j]=i;
      }
    }
  }  
  for (j=0;j<15;j++) {
    printf("Col %i  Max Error: %e at T=%e dlnT=%e (min %e max %e)\n",j,
	 maxerr[j],exp(cl->TlnMin+ierr[j]/cl->rDeltaTln),2/cl->rDeltaTln,min[j],max[j]);
  }
}

#define CL_Ccomp0 0.565e-9
#define CL_Tcmb0  2.735
#define CL_Ccomp  (CL_Ccomp0*CL_Tcmb0)

void clRatesRedshift( CL *cl, double z ) {
  int i;
  double xx;
  UVSPECTRUM *UV,*UV0;

  cl->z = z;
  cl->dComovingGmPerCcUnit = cl->dGmPerCcUnit*pow(1.+z,3.);
     
  cl->R.Cool_Comp = pow((1+z)*CL_Ccomp,4.0)*CL_B_gm;

  /* Photo-Ionization rates 
  ** Table in order of high to low redshift 
   */

  UV = cl->UV;
  for ( i=0; i < cl->nUV && z <= UV->z ; i++,UV++ );

  if (i == 0) {
       cl->R.Rate_Phot_HI = CL_RT_MIN;
       cl->R.Rate_Phot_HeI = CL_RT_MIN;
       cl->R.Rate_Phot_HeII = CL_RT_MIN;

       cl->R.Heat_Phot_HI = 0.0;
       cl->R.Heat_Phot_HeI = 0.0;
       cl->R.Heat_Phot_HeII = 0.0;
       return;
  }

  UV0=UV-1;
  if (i == cl->nUV ) {
       cl->R.Rate_Phot_HI = UV0->Rate_Phot_HI;
       cl->R.Rate_Phot_HeI = UV0->Rate_Phot_HeI;
       cl->R.Rate_Phot_HeII = UV0->Rate_Phot_HeII;

       cl->R.Heat_Phot_HI = UV0->Heat_Phot_HI*CL_B_gm;
       cl->R.Heat_Phot_HeI = UV0->Heat_Phot_HeI*CL_B_gm;
       cl->R.Heat_Phot_HeII = UV0->Heat_Phot_HeII*CL_B_gm;
  }
  else {
       xx = log((1+z)/(1+UV0->z))/log((1+UV->z)/(1+UV0->z));
       cl->R.Rate_Phot_HI = pow(UV0->Rate_Phot_HI,1-xx)*pow(UV->Rate_Phot_HI,xx);
       cl->R.Rate_Phot_HeI = pow(UV0->Rate_Phot_HeI,1-xx)*pow(UV->Rate_Phot_HeI,xx);
       cl->R.Rate_Phot_HeII = pow(UV0->Rate_Phot_HeII,1-xx)*pow(UV->Rate_Phot_HeII,xx);

       cl->R.Heat_Phot_HI = pow(UV0->Heat_Phot_HI,1-xx)*pow(UV->Heat_Phot_HI,xx)*CL_B_gm;
       cl->R.Heat_Phot_HeI = pow(UV0->Heat_Phot_HeI,1-xx)*pow(UV->Heat_Phot_HeI,xx)*CL_B_gm;
       cl->R.Heat_Phot_HeII = pow(UV0->Heat_Phot_HeII,1-xx)*pow(UV->Heat_Phot_HeII,xx)*CL_B_gm;
  }
  if (cl->R.Rate_Phot_HI < CL_RT_MIN) cl->R.Rate_Phot_HI = CL_RT_MIN;
  if (cl->R.Rate_Phot_HeI < CL_RT_MIN) cl->R.Rate_Phot_HeI = CL_RT_MIN;
  if (cl->R.Rate_Phot_HeII < CL_RT_MIN) cl->R.Rate_Phot_HeII = CL_RT_MIN;
  return;
}

void clRates ( CL *cl, RATE *Rate, double T ) {
  double Tln;
  double xTln,wTln0,wTln1;
  RATES_T *RT0,*RT1;
  int iTln;

  if (T >= cl->TMax) T=cl->TMax*(1.0 - EPS);   
  if (T < cl->TMin) T=cl->TMin;
  Tln = log(T);

  Rate->T = T;
  Rate->Tln = Tln; 

  xTln = (Tln-cl->TlnMin)*cl->rDeltaTln;
  iTln = xTln;
  RT0 = (cl->RT+iTln);
  RT1 = RT0+1; 
  wTln1 = xTln-iTln;
  wTln0 = 1-wTln1;
  
  Rate->Coll_HI = (wTln0*RT0->Rate_Coll_HI+wTln1*RT1->Rate_Coll_HI);
  Rate->Coll_HeI = (wTln0*RT0->Rate_Coll_HeI+wTln1*RT1->Rate_Coll_HeI);
  Rate->Coll_HeII = (wTln0*RT0->Rate_Coll_HeII+wTln1*RT1->Rate_Coll_HeII);

  Rate->Radr_HII = (wTln0*RT0->Rate_Radr_HII+wTln1*RT1->Rate_Radr_HII);
  Rate->Radr_HeII = (wTln0*RT0->Rate_Radr_HeII+wTln1*RT1->Rate_Radr_HeII);
  Rate->Diel_HeII = (wTln0*RT0->Rate_Diel_HeII+wTln1*RT1->Rate_Diel_HeII);
  Rate->Totr_HeII = Rate->Radr_HeII + Rate->Diel_HeII;
  Rate->Radr_HeIII = (wTln0*RT0->Rate_Radr_HeIII+wTln1*RT1->Rate_Radr_HeIII);
}

double clHeatTotal ( CL *cl, PERBARYON *Y ) {
  /* erg /gram /sec
     Note: QQ_* premultiplied by (CL_B_gm*erg_ev) */

  return Y->HI   * cl->R.Heat_Phot_HI * cl->R.Rate_Phot_HI +
         Y->HeI  * cl->R.Heat_Phot_HeI * cl->R.Rate_Phot_HeI +
         Y->HeII * cl->R.Heat_Phot_HeII * cl->R.Rate_Phot_HeII;
}

double clCoolTotal ( CL *cl, PERBARYON *Y, RATE *Rate, double rho ) {
  /* Assumes clRates called previously */
  /* erg /gram /sec */

  double en_B=rho*CL_B_gm;
  double xTln,wTln0,wTln1;
  RATES_T *RT0,*RT1;
  int iTln;

  xTln = (Rate->Tln-cl->TlnMin)*cl->rDeltaTln;
  iTln = xTln;
  RT0 = (cl->RT+iTln);
  RT1 = RT0+1; 
  wTln1 = xTln-iTln;
  wTln0 = 1-wTln1;

  /* PUT INTO erg/gm/sec */
  return Y->e * ( 
    cl->R.Cool_Comp * ( Rate->T-CL_Tcmb0 ) * ( 1 + cl->z ) + 
    en_B * (
    (wTln0*RT0->Cool_Brem_1+wTln1*RT1->Cool_Brem_1) * ( Y->HII + Y->HeII ) +
    (wTln0*RT0->Cool_Brem_2+wTln1*RT1->Cool_Brem_2) * Y->HeIII +

    (wTln0*RT0->Cool_Radr_HII+wTln1*RT1->Cool_Radr_HII) * Y->HII * Rate->Radr_HII +
    (wTln0*RT0->Cool_Radr_HeII+wTln1*RT1->Cool_Radr_HeII) * Y->HeII * Rate->Radr_HeII +
    (wTln0*RT0->Cool_Radr_HeIII+wTln1*RT1->Cool_Radr_HeIII) * Y->HeIII * Rate->Radr_HeIII +
 
    cl->R.Cool_Coll_HI * Y->HI * Rate->Coll_HI +
    cl->R.Cool_Coll_HeI * Y->HeI * Rate->Coll_HeI +
    cl->R.Cool_Coll_HeII * Y->HeII * Rate->Coll_HeII +

    cl->R.Cool_Diel_HeII * Y->HeII * Rate->Diel_HeII +

    (wTln0*RT0->Cool_Line_HI+wTln1*RT1->Cool_Line_HI) * Y->HI +
    (wTln0*RT0->Cool_Line_HeI+wTln1*RT1->Cool_Line_HeI) * Y->HeI +
    (wTln0*RT0->Cool_Line_HeII+wTln1*RT1->Cool_Line_HeII) * Y->HeII ) );
}

COOL_ERGPERSPERGM  clTestCool ( CL *cl, PERBARYON *Y, RATE *Rate, double rho ) {
  /* Assumes clRates called previously */
  /* erg /gram /sec */

  double en_B=rho*CL_B_gm;
  double xTln,wTln0,wTln1;
  RATES_T *RT0,*RT1;
  int iTln;
  COOL_ERGPERSPERGM ret;

  xTln = (Rate->Tln-cl->TlnMin)*cl->rDeltaTln;
  iTln = xTln;
  RT0 = (cl->RT+iTln);
  RT1 = RT0+1; 
  wTln1 = xTln-iTln;
  wTln0 = 1-wTln1;

  /* PUT INTO erg/gm/sec */
  ret.compton = Y->e * ( 
    cl->R.Cool_Comp * ( Rate->T-CL_Tcmb0 ) * ( 1 + cl->z ));
  ret.bremHII = Y->e * en_B * (
    (wTln0*RT0->Cool_Brem_1+wTln1*RT1->Cool_Brem_1) * ( Y->HII ));
  ret.bremHeII = Y->e * en_B * (
    (wTln0*RT0->Cool_Brem_1+wTln1*RT1->Cool_Brem_1) * ( Y->HeII ));
  ret.bremHeIII = Y->e * en_B * (
    (wTln0*RT0->Cool_Brem_2+wTln1*RT1->Cool_Brem_2) * Y->HeIII );
  ret.radrecHII = Y->e * en_B * 
    (wTln0*RT0->Cool_Radr_HII+wTln1*RT1->Cool_Radr_HII) * Y->HII * Rate->Radr_HII;
  ret.radrecHeII = Y->e * en_B * 
    (wTln0*RT0->Cool_Radr_HeII+wTln1*RT1->Cool_Radr_HeII) * Y->HeII * Rate->Radr_HeII;
  ret.radrecHeIII = Y->e * en_B * 
    (wTln0*RT0->Cool_Radr_HeIII+wTln1*RT1->Cool_Radr_HeIII) * Y->HeIII * Rate->Radr_HeIII;
  ret.collionHI = Y->e * en_B * 
    cl->R.Cool_Coll_HI * Y->HI * Rate->Coll_HI;
  ret.collionHeI = Y->e * en_B * 
    cl->R.Cool_Coll_HeI * Y->HeI * Rate->Coll_HeI;
  ret.collionHeII = Y->e * en_B * 
    cl->R.Cool_Coll_HeII * Y->HeII * Rate->Coll_HeII;
  ret.dielrecHeII = Y->e * en_B * 
    cl->R.Cool_Diel_HeII * Y->HeII * Rate->Diel_HeII;
  ret.lineHI = Y->e * en_B * 
    (wTln0*RT0->Cool_Line_HI+wTln1*RT1->Cool_Line_HI) * Y->HI;
  ret.lineHeI = Y->e * en_B * 
    (wTln0*RT0->Cool_Line_HeI+wTln1*RT1->Cool_Line_HeI) * Y->HeI;
  ret.lineHeII = Y->e * en_B * 
    (wTln0*RT0->Cool_Line_HeII+wTln1*RT1->Cool_Line_HeII) * Y->HeII;

  return ret;
}

void clPrintCool ( CL *cl, PERBARYON *Y, RATE *Rate, double rho ) {
  /* Assumes clRates called previously */
  /* erg /gram /sec */

  double en_B=rho*CL_B_gm;
  double xTln,wTln0,wTln1;
  RATES_T *RT0,*RT1;
  int iTln;

  xTln = (Rate->Tln-cl->TlnMin)*cl->rDeltaTln;
  iTln = xTln;
  RT0 = (cl->RT+iTln);
  RT1 = RT0+1; 
  wTln1 = xTln-iTln;
  wTln0 = 1-wTln1;

  /* PUT INTO erg/gm/sec */
  printf("Compton:  %e\n",
    Y->e * ( 
    cl->R.Cool_Comp * ( Rate->T-CL_Tcmb0 ) * ( 1 + cl->z )));
  printf("Cool Brem HII    %e\n",
    Y->e * en_B * (
    (wTln0*RT0->Cool_Brem_1+wTln1*RT1->Cool_Brem_1) * ( Y->HII )) );
  printf("Cool Brem HeII   %e\n",
    Y->e * en_B * (
    (wTln0*RT0->Cool_Brem_1+wTln1*RT1->Cool_Brem_1) * ( Y->HeII )) );
  printf("Cool Brem HeIII  %e\n",
    Y->e * en_B * (
    (wTln0*RT0->Cool_Brem_2+wTln1*RT1->Cool_Brem_2) * Y->HeIII ) );
  printf("Radiative Recombination  HII    %e\n",
    Y->e * en_B * 
    (wTln0*RT0->Cool_Radr_HII+wTln1*RT1->Cool_Radr_HII) * Y->HII * Rate->Radr_HII );
  printf("Radiative Recombination  HeII   %e\n",
    Y->e * en_B * 
    (wTln0*RT0->Cool_Radr_HeII+wTln1*RT1->Cool_Radr_HeII) * Y->HeII * Rate->Radr_HeII);
  printf("Radiative Recombination  HeIII  %e\n",
    Y->e * en_B * 
    (wTln0*RT0->Cool_Radr_HeIII+wTln1*RT1->Cool_Radr_HeIII) * Y->HeIII * Rate->Radr_HeIII);
  printf("Collisional Ionization  HI    %e\n",
    Y->e * en_B * 
    cl->R.Cool_Coll_HI * Y->HI * Rate->Coll_HI);
  printf("Collisional Ionization  HeI   %e\n",
    Y->e * en_B * 
    cl->R.Cool_Coll_HeI * Y->HeI * Rate->Coll_HeI);
  printf("Collisional Ionization  HeII  %e\n",
    Y->e * en_B * 
    cl->R.Cool_Coll_HeII * Y->HeII * Rate->Coll_HeII);
  printf("Dielectric Recombination HeII %e\n",
    Y->e * en_B * cl->R.Cool_Diel_HeII * Y->HeII * Rate->Diel_HeII);

  printf("Line cooling HI   %e\n",
    Y->e * en_B * 
    (wTln0*RT0->Cool_Line_HI+wTln1*RT1->Cool_Line_HI) * Y->HI);
  printf("Line cooling HeI  %e\n",
    Y->e * en_B * 
    (wTln0*RT0->Cool_Line_HeI+wTln1*RT1->Cool_Line_HeI) * Y->HeI);
  printf("Line cooling HeII  %e\n",
    Y->e * en_B * 
    (wTln0*RT0->Cool_Line_HeII+wTln1*RT1->Cool_Line_HeII) * Y->HeII );
}

void clAbunds( CL *cl, PERBARYON *Y, RATE *Rate, double rho ) {
  double en_B=rho*CL_B_gm;
  double rcirrHI   = (Rate->Coll_HI)/(Rate->Radr_HII);
  double rcirrHeI  = (Rate->Coll_HeI)/(Rate->Totr_HeII);
  double rcirrHeII = (Rate->Coll_HeII)/(Rate->Radr_HeIII);
  double rpirrHI   = (cl->R.Rate_Phot_HI)/(Rate->Radr_HII * en_B);
  double rpirrHeI  = (cl->R.Rate_Phot_HeI)/(Rate->Totr_HeII * en_B);
  double rpirrHeII = (cl->R.Rate_Phot_HeII)/(Rate->Radr_HeIII * en_B);
  double yHI = 0.;
  double yHeI = 0.;
  double yHeII = 0.;
  double yH = cl->Y_H;
  double yHe = cl->Y_He;
  double yeMAX = cl->Y_eMAX;
  double rye,ye;
  double fHI,fHeI,fHeII,rfHe,yHI_old,yHeII_old;
  int i;  

  for ( i=0 ; i<MAXABUNDITERATIONS ; i++ ) {
    yHI_old = yHI;
    yHeII_old = yHeII;

    ye = (yeMAX-(yHI + 2 * yHeI + yHeII));
    if (ye <= 0) {
      Y->e = 0;
      Y->HI = yH;
      Y->HII = 0;
      Y->HeI = yHe;
      Y->HeII = 0;
      Y->HeIII = 0;
      Y->Total = yH + yHe;
      return;
    }
    rye = 1/ye;
    fHI = rcirrHI + rpirrHI * rye;
    fHeI = rcirrHeI + rpirrHeI * rye;
    rfHe = 1 / ( 1 + fHeI * (1+rcirrHeII+rpirrHeII*rye) );

    yHI = yH / (1.0+fHI);
    yHeI = yHe * rfHe;
    yHeII = yHe * fHeI * rfHe;
   
    if ( fabs(yHeII_old-yHeII) < EPS * yHeII && fabs(yHI_old-yHI) < EPS * yHI ) break;
  }

  Y->e = ye;
  Y->HI = yHI;
  Y->HII = yH / (1.0/fHI+1.0);
  Y->HeI = yHeI;
  Y->HeII = yHeII;
  fHeII = rcirrHeII + rpirrHeII*rye;
  Y->HeIII = yHe / ((1.0/fHeI+1.0)/fHeII+1.0);
  Y->Total = Y->e + yH + yHe;
}

#define CL_Rgascode         8.2494e7
#define CL_Eerg_gm_degK     CL_Rgascode
#define CL_ev_degK          1.0/1.1604e4
#define CL_Eerg_gm_ev       CL_Eerg_gm_degK/CL_ev_degK
#define CL_Eerg_gm_degK3_2  1.5*CL_Eerg_gm_degK

/* 
 * Though 13.6eV is lost to the Gas as radiation during H recombination, calculating the
 * Energy using u = E per unit mass = 3/2 n/rho k T requires we don't subtract it there.
 * This formulation is useful because then pressure = 2/3 u rho.
 * Instead we subtract the 13.6eV for H collisional ionization which actually
 * removes no energy from the Gas ( similarly for Helium ) 
 * It also means photoionization doesn't add the 13.6eV, only the excess.
 */
double clThermalEnergy( double Y_Total, double T ) {
  return Y_Total*CL_Eerg_gm_degK3_2*T;

}

double clTemperature( double Y_Total, double E ) {
  return E/(Y_Total*CL_Eerg_gm_degK3_2);
}

/*-----------------------------------------------------------------
 *     Collisional Ionization rates
 *-----------------------------------------------------------------*/
/*     H + e- -> H+ + 2e-  Janev et al. 1987 (Abel 1996) */
double clRateCollHI( double T ) {
  double TL,arg;
    
  TL = log(T*CL_eV_per_K);
  arg = -32.713967867 + TL*(13.536556     + TL*(-5.73932875 +
      TL*(1.56315498 +
      TL*(-0.2877056     + TL*(3.48255977e-2 + TL*(-2.63197617e-3 +
      TL*(1.11954395e-4  + TL*(-2.03914985e-6))))))));
  if (arg < CL_MAX_NEG_EXP_ARG) return 0;
  return exp ( arg );
}
  
/*     He + e- -> He+ + 2e-  Janev et al. 1987 (Abel 1996) */
double clRateCollHeI( double T ) {
  double TL,arg;
    
  TL = log(T*CL_eV_per_K);
  arg = -44.09864886  + TL*(23.91596563   + TL*(-10.7532302 + 
      TL*(3.05803875 + 
      TL*(-0.56851189    + TL*(6.79539123e-2 + TL*(-5.00905610e-3 +
      TL*(2.06723616e-4  + TL*(-3.64916141e-6))))))));
  if (arg < CL_MAX_NEG_EXP_ARG) return 0;
  return exp( arg );
}

/*     He+ + e- -> He++ + 2e- Aladdin Database 1989 (Abel 1996) */
double clRateCollHeII( double T ) {
  double TL,arg;

  TL = log(T*CL_eV_per_K);
  arg = -68.71040990  + TL*(43.93347633   + TL*(-18.4806699 + 
      TL*(4.70162649 +
      TL*(-0.76924663    + TL*(8.113042e-2   + TL*(-5.32402063e-3 + 
      TL*(1.97570531e-4  + TL*(-3.16558106e-6))))))));
  if (arg < CL_MAX_NEG_EXP_ARG) return 0;
  return exp( arg );
}

/*-----------------------------------------------------------------
 *     Radiative Recombination rates
 *-----------------------------------------------------------------*/
/*     H+ + e- -> H + gam  Verner & Ferland 1996 */
double clRateRadrHII( double T ) {
  double Tsq = sqrt(T);

  return 7.982e-11/( Tsq*0.563615 *
		pow(1+Tsq*0.563615,0.252) * pow(1+Tsq*1.192167e-3,1.748));
}

/*     He+ + e- -> He + gam  radiative  Verner & Ferland 1996 */
double clRateRadrHeII( double T ) {
  /* 
   * Note that these functions do not meet perfectly at 1e6 -- 2% difference 
   * The derivatives are different there also: So the apparent error is large 
   */
  double Tsq = sqrt(T);
  if (T < 1e6)
    return  3.294e-11/( Tsq*0.253673 *
             pow(1+Tsq*0.253673,0.309) * pow(1+Tsq*1.649348e-4,1.691));
  else
    return  9.356e-10/( Tsq*4.841607 *
             pow(1+Tsq*4.841607,0.2108) * pow(1+Tsq*4.628935e-4,1.7892));
}

/*     He+ + e- -> He + gam  dielectronic  Aldovandi&Pequignot 1973 (Black 1981) */
double clRateDielHeII( double T ) {
  double T_inv = 1.0/T,arg;

  arg = -4.7e5*T_inv;
  if (arg < CL_MAX_NEG_EXP_ARG) return 0;
  return 1.9e-3*pow(T,-1.5)*exp(arg)*(1+0.3*exp(-9.4e4*T_inv));
}

/*     He++ + e- -> He+ + gam  Verner & Ferland 1996 */
double clRateRadrHeIII( double T ) {
  double Tsq = sqrt(T);

  return 1.891e-10/( Tsq*0.326686 *
          pow(1+Tsq*0.326686,0.2476) * pow(1+Tsq*6.004084e-4,1.7524));
}

/*-----------------------------------------------------------------
 *     Bremsstrahlung   
 *-----------------------------------------------------------------*/
#define CL_Cbremss1 1.426e-27
#define CL_al       0.79464
#define CL_bl       0.1243
#define CL_ar       2.13164
#define CL_br       (-0.1240)

double clCoolBrem1( double T ) {
  double Tlog10, Tsq;

  Tlog10 = log10(T);
  Tsq = sqrt(T);
  if (T < 3.2e5) 
     return Tsq*CL_Cbremss1*(CL_al+CL_bl*Tlog10)*CL_B_gm;
  else   
     return Tsq*CL_Cbremss1*(CL_ar+CL_br*Tlog10)*CL_B_gm;
}

#define CL_alog4   0.602059991
#define CL_alII    (4.0*(CL_al-CL_bl*CL_alog4))
#define CL_blII    (4.0*CL_bl)
#define CL_arII    (4.0*(CL_ar-CL_br*CL_alog4))
#define CL_brII    (4.0*CL_br)

double clCoolBrem2( double T ) {
  double Tlog10, Tsq;

  Tlog10 = log10(T);
  Tsq = sqrt(T);

  if (T<12.8e5) 
    return Tsq*CL_Cbremss1*(CL_alII+CL_blII*Tlog10)*CL_B_gm;
  else
    return Tsq*CL_Cbremss1*(CL_arII+CL_brII*Tlog10)*CL_B_gm;
}

/*-----------------------------------------------------------------
 *     Cooling multiplier for radiative recombination
 *-----------------------------------------------------------------*/
#define CL_aHII  0.0215964
#define CL_b     0.270251

double clCoolRadrHII( double T ) {

  double Tpow;
    
  Tpow=pow(T,CL_b);
  /* return CL_B_gm*(CL_eHI+exp(-CL_aHII*Tpow)*CL_k_Boltzmann*T); */
  /* Though 13.6eV is lost to the Gas as radiation, calculating the
   * Energy using u = 3/2 k T requires we don't subtract it here.
   */
  return CL_B_gm*(exp(-CL_aHII*Tpow)*CL_k_Boltzmann*T);
}
 
double clCoolRadrHeII( double T ) {

  double Tpow;
    
  Tpow=pow(T,CL_b);
  /* return CL_B_gm*(CL_eHeI+exp(-(CL_aHII*pow(13.6/24.59,CL_b))*Tpow)*CL_k_Boltzmann*T); */
  return CL_B_gm*(exp(-(CL_aHII*pow(13.6/24.59,CL_b))*Tpow)*CL_k_Boltzmann*T);
}

double clCoolRadrHeIII( double T ) {
  double Tpow;
    
  Tpow=pow(T,CL_b);
  /* return CL_B_gm*(CL_eHeII+exp(-(CL_aHII*pow(13.6/54.42,CL_b))*Tpow)*CL_k_Boltzmann*T); */
  return CL_B_gm*(exp(-(CL_aHII*pow(13.6/54.42,CL_b))*Tpow)*CL_k_Boltzmann*T);
}

/*-----------------------------------------------------------------
 *     Line Cooling
 *-----------------------------------------------------------------*/
/*      CEN (1992, Ap.J.Suppl 78,341) ADVOCATES MULTIPLYING EACH OF 
 *      THESE RATES BY Cen_correctn - HE CLAIMS THIS GIVES THE RIGHT
 *      HIGH T LIMIT FOR PROCESSES INVOLVING A FREE EL INTERACTING 
 *      WITH AN ORBITAL ELECTRON ?? */
#define CL_aHI  7.5e-19
#define CL_bHI  1.18348e05

double clCoolLineHI( double T ) {
  double T_inv, arg;
  double Cen_correctn = 1.0/(1.0+sqrt(T/1.0e5));

  T_inv=1.0/T;
  arg = -CL_bHI*T_inv;
  if (arg < CL_MAX_NEG_EXP_ARG) return 0;
  return CL_B_gm*CL_aHI*exp( arg )*Cen_correctn;
}

#define CL_aHeI   9.10e-27
#define CL_bHeI   1.3179e04
#define CL_p_HeI  0.1687

double clCoolLineHeI( double T ) {
  double T_inv,arg;
  double Cen_correctn = 1.0/(1.0+sqrt(T/1.0e5));

  T_inv=1.0/T;
  arg = -CL_bHeI*T_inv;
  if (arg < CL_MAX_NEG_EXP_ARG) return 0;
  return CL_B_gm*CL_aHeI*exp(-CL_bHeI*T_inv)*pow(T_inv,CL_p_HeI)*Cen_correctn;
}

#define CL_aHeII   5.54e-17
#define CL_bHeII   4.73638e05
#define CL_p_HeII  0.397

double clCoolLineHeII( double T ) {

  double T_inv,arg;
  double Cen_correctn = 1.0/(1.0+sqrt(T/1.0e5));

  T_inv=1.0/T;
  arg = -CL_bHeII*T_inv;
  if (arg < CL_MAX_NEG_EXP_ARG) return 0;
  return CL_B_gm*CL_aHeII*exp(-CL_bHeII*T_inv)*pow(T_inv,CL_p_HeII)*Cen_correctn;
}

/* 
 * This ugly function is an huge improvement over the standard clCoolTotal
 * because it adds in the heating and cooling associated with
 * sudden changes in ionization (unlike most Equlibrium codes).
 */

#define SIGFRAC 1e-7

double clEdot( CL *cl, PERBARYON *Y, RATE *Rate, double rho, 
		      PERBARYON *Y1, PERBARYON *Y2, double dt) {
  double en_B = rho*CL_B_gm;
  double xTln,wTln0,wTln1;
  RATES_T *RT0,*RT1;
  int iTln;

  double Edot,ne,ndot,ndot_Phot,ndot_Coll,ndot_Radr_in,ndot_Radr,ndot_Diel,add;

  ne = Y->e*en_B;

  xTln = (Rate->Tln-cl->TlnMin)*cl->rDeltaTln;
  iTln = xTln;
  RT0 = (cl->RT+iTln);
  RT1 = RT0+1; 
  wTln1 = (xTln-iTln);
  wTln0 = (1-wTln1);

  Edot = - Y->e * cl->R.Cool_Comp * ( Rate->T-CL_Tcmb0 ) * ( 1 + cl->z ) -
    ne * ((wTln0*RT0->Cool_Brem_1+wTln1*RT1->Cool_Brem_1) * ( Y->HII + Y->HeII ) +
	   (wTln0*RT0->Cool_Brem_2+wTln1*RT1->Cool_Brem_2) * Y->HeIII +
	   cl->R.Cool_Diel_HeII * Y->HeII * Rate->Diel_HeII +
	   (wTln0*RT0->Cool_Line_HI+wTln1*RT1->Cool_Line_HI) * Y->HI +
	   (wTln0*RT0->Cool_Line_HeI+wTln1*RT1->Cool_Line_HeI) * Y->HeI +
	   (wTln0*RT0->Cool_Line_HeII+wTln1*RT1->Cool_Line_HeII) * Y->HeII );

  /* 
   * In equilibrium ndot_Phot + ndot_Coll - ndot_Radr = 0
   * However, there has been a net change ndot so we adjust the rates
   * so that ndot_Phot + ndot_Coll - ndot_Radr = ndot
   * and the change in Sum_i (log(ndot_i)-log(ndot_i_in))^2 is minimized 
   * We do this for HI,HeI and HeII 
   */
  ndot = (Y1->HI - Y2->HI)/dt;
  ndot_Phot = cl->R.Rate_Phot_HI * Y->HI + NDOT_MIN;
  ndot_Coll = Rate->Coll_HI * ne * Y->HI;
  ndot_Radr_in = Rate->Radr_HII * ne * Y->HII + NDOT_MIN;
  add = ndot_Phot*ndot_Radr_in + ndot_Coll*ndot_Radr_in;
  if ( add < ndot*fabs(ndot)*SIGFRAC )
    ndot_Radr = add/ndot;
  else
    ndot_Radr = (-.5*ndot + sqrt( 0.25*ndot*ndot + add ));
  /* 
  printf("HI %e ndot %e %e %e -> %e %e %e\n",ndot,ndot_Phot,ndot_Coll,ndot_Radr_in, 
	 ndot_Phot*ndot_Radr_in/ndot_Radr,ndot_Coll*ndot_Radr_in/ndot_Radr, ndot_Radr );
  */
  ndot_Phot = ndot_Phot*ndot_Radr_in/ndot_Radr;
  ndot_Coll = ndot_Coll*ndot_Radr_in/ndot_Radr;

  Edot += ndot_Phot * cl->R.Heat_Phot_HI
    - ndot_Coll * cl->R.Cool_Coll_HI
    - ndot_Radr * (wTln0*RT0->Cool_Radr_HII+wTln1*RT1->Cool_Radr_HII);

  ndot = (Y1->HeI - Y2->HeI)/dt;
  ndot_Phot = cl->R.Rate_Phot_HeI * Y->HeI + NDOT_MIN;
  ndot_Coll = Rate->Coll_HeI * ne * Y->HeI;
  ndot_Radr_in = Rate->Radr_HeII * ne * Y->HeII + NDOT_MIN;
  ndot_Diel = Rate->Diel_HeII * ne * Y->HeII;
  add = ndot_Phot*ndot_Radr_in + ndot_Coll*ndot_Radr_in;
  if ( add < ndot*fabs(ndot)*SIGFRAC )
    ndot_Radr = add/( ndot* ( 1 + ndot_Diel/ndot_Radr_in ) );
  else
    ndot_Radr = (-.5*ndot + sqrt( 0.25*ndot*ndot + add ))/( 1 + ndot_Diel/ndot_Radr_in );
  /*
  printf("HeI %e ndot %e %e %e %e -> %e %e %e %e\n",ndot,ndot_Phot,ndot_Coll,ndot_Radr_in,ndot_Diel, 
	 ndot_Phot*ndot_Radr_in/ndot_Radr,ndot_Coll*ndot_Radr_in/ndot_Radr, ndot_Radr, ndot_Radr*ndot_Diel/ndot_Radr_in);
  */
  ndot_Phot = ndot_Phot*ndot_Radr_in/ndot_Radr;
  ndot_Coll = ndot_Coll*ndot_Radr_in/ndot_Radr;
  ndot_Diel = ndot_Diel*ndot_Radr/ndot_Radr_in;

  Edot += ndot_Phot * cl->R.Heat_Phot_HeI
    - ndot_Coll * cl->R.Cool_Coll_HeI
    - ndot_Radr * (wTln0*RT0->Cool_Radr_HeII+wTln1*RT1->Cool_Radr_HeII)
    - ndot_Diel * cl->R.Cool_Diel_HeII;

  ndot = (Y2->HeIII - Y1->HeIII)/dt; /* note: Use HeIII to get HeII rates */
  ndot_Phot = cl->R.Rate_Phot_HeII * Y->HeII + NDOT_MIN;
  ndot_Coll = Rate->Coll_HeII * ne * Y->HeII;
  ndot_Radr_in = Rate->Radr_HeIII * ne * Y->HeIII + NDOT_MIN;
  add = ndot_Phot*ndot_Radr_in + ndot_Coll*ndot_Radr_in;
  if ( add < ndot*fabs(ndot)*SIGFRAC )
    ndot_Radr = add/ndot;
  else
    ndot_Radr = (-.5*ndot + sqrt( 0.25*ndot*ndot + add ));
  /*
  printf("HeII %e ndot %e %e %e -> %e %e %e\n",ndot,ndot_Phot,ndot_Coll,ndot_Radr_in, 
	 ndot_Phot*ndot_Radr_in/ndot_Radr,ndot_Coll*ndot_Radr_in/ndot_Radr, ndot_Radr );
  */
  ndot_Phot = ndot_Phot*ndot_Radr_in/ndot_Radr;
  ndot_Coll = ndot_Coll*ndot_Radr_in/ndot_Radr;

  Edot += ndot_Phot * cl->R.Heat_Phot_HeII
    - ndot_Coll * cl->R.Cool_Coll_HeII
    - ndot_Radr * (wTln0*RT0->Cool_Radr_HeIII+wTln1*RT1->Cool_Radr_HeIII);

  /* 
  printf("Edot %e %e\n",Edot,Edot*dt); 
  */
  return Edot;
}

/*
 *  We solve the implicit equation:  
 *  Eout = Ein + dt * Cooling ( Yin, Yout, Ein, Eout )
 *
 *  E erg/gm, PdV erg/gm/s, rho gm/cc, dt sec
 * 
 *  Brent's method is fine but you need a bracketing
 *  I bracket using Ein and Ein+Edot(Ein)*dt*A 
 *  where A is some power of 2
 */

#define MAXBRACKET 10

#define ITMAX 100
#define BRENTEPS 1e-5
#define BRENTTOL 1e6

typedef struct {
  PERBARYON Y;
  double E,T,F;
} BRENT;
  
double clEdotHarmonic( CL *cl, RATE *R1, RATE *R2, PERBARYON *Y1, PERBARYON *Y2, 
		       double rho, double PdV, double dt) {
  double Edot1,Edot2;
  double Prod;

  Edot1 =  clEdot( cl, Y1, R1, rho, Y1, Y2, dt ) + PdV;
  Edot2 =  clEdot( cl, Y2, R2, rho, Y1, Y2, dt ) + PdV;

  Prod = Edot1*Edot2;
  if (Prod < 0) return Edot2; /* Revert to First order always stable endpoint derivative */

  return Prod/(0.5*(Edot1+Edot2)); /* Second order Harmonic mean derivative */
}

void clIntegrateEnergy(CL *cl, PERBARYON *Y, double *E, 
		       double PdV, double rho, double dt ) {
  BRENT a,b,c;
  PERBARYON Yin = *Y, YY;
  RATE Rate,Ratein;
  double Ein = *E, dE, Tin, Ttmp;
  int i;
  /* brent */
  int iter;
  double d=0,e=0,min1,min2;
  double p,q,r,s,tol1,xm;

  /* debug */  /*
  PARTICLE *part = cl->p;

  char ach[256];
  */

  if (dt <= 0) return;

  /* 
   *  We recalculate Y based on the current table (rather than use Yin)
   *  to take into account changes in photoionization for dE 
   */
  Tin = clTemperature( Yin.Total, Ein ); 
  clRates( cl, &Ratein, Tin );
  clAbunds( cl, &YY, &Ratein, rho );
  Ttmp = clTemperature( YY.Total, Ein ); 
  clRates( cl, &Rate, Ttmp );
  clAbunds( cl, &YY, &Rate, rho );

  dE  = dt * clEdotHarmonic( cl, &Ratein, &Rate, &Yin, &YY, rho, PdV, dt );
  if (dE < -0.2*Ein) dE = -0.2*Ein;
  /*
  if (part->iOrder==880556) {
     sprintf(ach,"%d: T: %e E %e dE: %e dt %e   %e %e %e\n",part->iOrder,Tin,Ein,dE,dt,
	 -clCoolTotal( cl, &YY, &Ratein, rho ),clHeatTotal( cl, &YY ),PdV);
     mdlDiag(cl->mdl, ach);
  }
  */
  a.Y = Yin;
  a.E = Ein;
  a.T = Tin;
  a.F = - dE;
  
  /* Bracket the root */
  b.Y = YY;

  for ( i = 1; i <= MAXBRACKET; i++) {
    b.E = Ein + dE;
    b.T = clTemperature( b.Y.Total, b.E );
    if (b.T > cl->TMin) {
      clRates( cl, &Rate, b.T );
      clAbunds( cl, &b.Y, &Rate, rho );
      b.T = clTemperature( b.Y.Total, b.E ); 
    }
    if (b.T <= cl->TMin) {
      b.T = cl->TMin;
      clRates( cl, &Rate, b.T );
      clAbunds( cl, &b.Y, &Rate, rho );
      b.E = clThermalEnergy( b.Y.Total, b.T );
      b.F = b.E - Ein - dt * clEdotHarmonic( cl, &Ratein, &Rate, &Yin, &b.Y, rho, PdV, dt );
      if ( b.F * a.F < 0 ) break;
      /* Stick to the minimum Temperature */
      /*
      if (part->iOrder==880556) {
	sprintf(ach,"Floored it after %i brackets %g %g %g %g\n   (%g %g %g %g %g)\n",i,a.T,a.F,b.T,b.F,b.E,Ein,
	     dt*clEdotHarmonic( cl, &Ratein, &Rate, &Yin, &b.Y, rho, PdV, dt ),
	     dt*( clEdot( cl, &Yin, &Ratein, rho, &Yin, &b.Y, dt ) + PdV),
	     dt*( clEdot( cl, &b.Y, &Rate, rho, &Yin, &b.Y, dt ) + PdV) );
	mdlDiag(cl->mdl, ach);
      }
      */
      *Y = b.Y;
      *E = b.E;
      return;
    }
    clRates( cl, &Rate, b.T ); 
    clAbunds( cl, &b.Y, &Rate, rho );

    b.F = b.E - Ein - dt * clEdotHarmonic( cl, &Ratein, &Rate, &Yin, &b.Y, rho, PdV, dt );
    /*
    if (part->iOrder==880556) {
      sprintf(ach,"%d B %d T: %e E %e dE: %e dt %g Edot %g %e %e PdV %e F %g\n",part->iOrder,i,b.T,b.E,dE,dt,
	      clEdotHarmonic( cl, &Ratein, &Rate, &Yin, &b.Y, rho, PdV, dt ),
	      -clCoolTotal( cl, &b.Y, &Rate, rho ),clHeatTotal( cl, &b.Y ),PdV,b.F);
      mdlDiag(cl->mdl, ach);
    }
    */

    if ( b.F * a.F < 0 ) break;
    if (i<5 && Ein+2*dE > 0) dE *= 2;
    else {
      if (dE < 0) dE -= Ein*0.2;
      else dE = dE*2 + Ein*0.2;
    }
  }

#if (0) 
  if (i>MAXBRACKET) {
    /* Try going the other direction */
    dE = a.F;
    b.Y = Yin;
    
    for ( i = 1; i <= MAXBRACKET; i++) {
      b.E = Ein + dE;
      b.T = clTemperature( b.Y.Total, b.E );
      if (b.T > cl->TMin) {
	clRates( cl, &Rate, b.T );
	clAbunds( cl, &b.Y, &Rate, rho );
	b.T = clTemperature( b.Y.Total, b.E ); 
      }
      if (b.T <= cl->TMin) {
	b.T = cl->TMin;
	clRates( cl, &Rate, b.T );
	clAbunds( cl, &b.Y, &Rate, rho );
	b.E = clThermalEnergy( b.Y.Total, b.T );
	b.F = b.E - Ein - dt * clEdotHarmonic( cl, &Ratein, &Rate, &Yin, &b.Y, rho, PdV, dt );
	if ( b.F * a.F < 0 ) break;
	/* Stick to the minimum Temperature */
	/*
	printf("Floored it after %i brackets\n",i);
	*/
	*Y = b.Y;
	*E = b.E;
	return;
      }
      clRates( cl, &Rate, b.T ); 
      clAbunds( cl, &b.Y, &Rate, rho );
      
      b.F = b.E - Ein - dt * clEdotHarmonic( cl, &Ratein, &Rate, &Yin, &b.Y, rho, PdV, dt );
      if ( b.F * a.F < 0 ) break;
      if (i<5) dE *= 2;
      else {
	if (dE < 0) dE -= Ein*0.2;
	else dE = dE*2 + Ein*0.2;
      }
    }
  }
#endif
  /*
  if (i>MAXBRACKET) {
    sprintf(ach,"COOL: %d %d MaxBracket Y %g %g %g E %g PdV %g rho %g dt %g\n", part->iOrder,part->iRung,Y->HI,Y->HeI,Y->HeII, *E, PdV, rho, dt);
    mdlDiag(cl->mdl, ach);
  }  
  */
  assert(i<=MAXBRACKET);
  /*
    printf("Bracket Iterations %i\n",i);
  */

  /* Brent it - Algorithm taken from Numerical Recipes */  
  c = b;

  for (iter=1;iter<=ITMAX;iter++) {
    if ((b.F > 0.0 && c.F > 0.0) || (b.F < 0.0 && c.F < 0.0)) {
      c = a;
      e = d = b.E - a.E;
    }
    if (fabs(c.F) < fabs(b.F)) {
      a = b;
      b = c;
      c = a;
    }
    tol1 = 2*BRENTEPS*fabs(b.E) + 0.5*BRENTTOL;
    xm = 0.5 * (c.E - b.E);
    
    if (fabs(xm) <= tol1 || b.F == 0.0) {
      /*
      printf("Brent iterations: %i\n",iter);
      */
      /* finished -> b.E */
      *E = b.E;
      *Y = b.Y;
      return;
    }
    if (fabs(e) >= tol1 && fabs(a.F) > fabs(b.F)) {
      s=b.F/a.F;
      if (a.E == c.E) {
        p = 2*xm*s;
        q = 1-s;
      } else {
        q = a.F/c.F;
        r = b.F/c.F;
        p = s*(2.0*xm*q*(q-r)-(b.E-a.E)*(r-1.0));
        q = (q-1.0)*(r-1.0)*(s-1.0);
      }
      if (p > 0.0) q = -q;
      p=fabs(p);
      min1 = 3.0*xm*q - fabs(tol1*q);
      min2 = fabs(e*q);
      if (2.0*p < (min1 < min2 ? min1 : min2)) {
        e=d;
        d=p/q;
      } else {
        d=xm;
        e=d;
      }
    } else {
      d=xm;
      e=d;
    }
    a=b;
    if (fabs(d) > tol1)
      b.E += d;
    else
      b.E += (xm > 0.0 ? fabs(tol1) : -fabs(tol1));

    b.T = clTemperature( b.Y.Total, b.E ); 
    clRates( cl, &Rate, b.T );
    clAbunds( cl, &b.Y, &Rate, rho );
    b.T = clTemperature( b.Y.Total, b.E ); 
    clRates( cl, &Rate, b.T ); 
    clAbunds( cl, &b.Y, &Rate, rho );

    b.F = b.E - Ein - dt * clEdotHarmonic( cl, &Ratein, &Rate, &Yin, &b.Y, rho, PdV, dt );
  }
  assert( 0 ); /* Brent failed to converge */
}

#ifdef COOLDEBUG
void clIntegrateEnergyDEBUG(CL *cl, PERBARYON *Y, double *E, 
		       double PdV, double rho, double dt ) {
  BRENT a,b,c;
  PERBARYON Yin = *Y, YY;
  RATE Rate,Ratein;
  double Ein = *E, dE, Tin;
  int i;
  /* brent */
  int iter;
  double d,e,min1,min2;
  double p,q,r,s,tol1,xm;
  
  printf("Ein: %g dt: %g enB: %g Pdt: %g ",
	 Ein,dt,rho*CL_B_gm,dt*PdV);

  if (dt <= 0) return;

  /* 
   *  We recalculate Y based on the current table (rather than use Yin)
   *  to take into account changes in photoionization for dE 
   */
  Tin = clTemperature( Yin.Total, Ein ); 
  clRates( cl, &Ratein, Tin );
  clAbunds( cl, &YY, &Ratein, rho );
  Tin = clTemperature( Yin.Total, Ein ); 
  clRates( cl, &Ratein, Tin );
  clAbunds( cl, &YY, &Ratein, rho );

  printf("Cdt: %g Hdt: %g Tin: %g ",dt*clCoolTotal(cl,&YY,&Ratein,rho),dt*clHeatTotal(cl,&YY),Tin);

  dE  = dt * clEdotHarmonic( cl, &Ratein, &Ratein, &Yin, &YY, rho, PdV, dt );
  if (dE < -0.2*Ein) dE = -0.2*Ein;
  /*
    printf("T: %e E %e dE: %e dt %e   %e %e %e\n",Tin,Ein,dE,dt,
	 -clCoolTotal( cl, &YY, &Ratein, rho ),clHeatTotal( cl, &YY ),PdV);
  */
  a.Y = Yin;
  a.E = Ein;
  a.T = Tin;
  a.F = - dE;
  
  /* Bracket the root */
  b.Y = Yin;

  for ( i = 1; i <= MAXBRACKET; i++) {
    b.E = Ein + dE;
    b.T = clTemperature( b.Y.Total, b.E );
    if (b.T > cl->TMin) {
      clRates( cl, &Rate, b.T );
      clAbunds( cl, &b.Y, &Rate, rho );
      b.T = clTemperature( b.Y.Total, b.E ); 
    }
    if (b.T <= cl->TMin) {
      b.T = cl->TMin;
      clRates( cl, &Rate, b.T );
      clAbunds( cl, &b.Y, &Rate, rho );
      b.E = clThermalEnergy( b.Y.Total, b.T );
      b.F = b.E - Ein - dt * clEdotHarmonic( cl, &Ratein, &Rate, &Yin, &b.Y, rho, PdV, dt );
      if ( b.F * a.F < 0 ) break;
      /* Stick to the minimum Temperature */
      /*
      printf("Floored it after %i brackets %g %g %g %g\n   (%g %g %g %g %g)\n",i,a.T,a.F,b.T,b.F,b.E,Ein,
	     dt*clEdotHarmonic( cl, &Ratein, &Rate, &Yin, &b.Y, rho, PdV, dt ),
	     dt*( clEdot( cl, &Yin, &Ratein, rho, &Yin, &b.Y, dt ) + PdV),
	     dt*( clEdot( cl, &b.Y, &Rate, rho, &Yin, &b.Y, dt ) + PdV) );
      */
      *Y = b.Y;
      *E = b.E;
      printf("Floored Eout: %g, Tout: %g ",b.E,b.T);
      return;
    }
    clRates( cl, &Rate, b.T ); 
    clAbunds( cl, &b.Y, &Rate, rho );

    b.F = b.E - Ein - dt * clEdotHarmonic( cl, &Ratein, &Rate, &Yin, &b.Y, rho, PdV, dt );
    if ( b.F * a.F < 0 ) break;
    if (i<5) dE *= 2;
    else {
      if (dE < 0) dE -= Ein*0.2;
      else dE = dE*2 + Ein*0.2;
    }
  }
  
  if (i>MAXBRACKET) {
    /* Try going the other direction */
    dE = a.F; /* -dE original */
    b.Y = Yin;
    
    for ( i = 1; i <= MAXBRACKET; i++) {
      b.E = Ein + dE;
      b.T = clTemperature( b.Y.Total, b.E );
      if (b.T > cl->TMin) {
	clRates( cl, &Rate, b.T );
	clAbunds( cl, &b.Y, &Rate, rho );
	b.T = clTemperature( b.Y.Total, b.E ); 
      }
      if (b.T <= cl->TMin) {
	b.T = cl->TMin;
	clRates( cl, &Rate, b.T );
	clAbunds( cl, &b.Y, &Rate, rho );
	b.E = clThermalEnergy( b.Y.Total, b.T );
	b.F = b.E - Ein - dt * clEdotHarmonic( cl, &Ratein, &Rate, &Yin, &b.Y, rho, PdV, dt );
	if ( b.F * a.F < 0 ) break;
	/* Stick to the minimum Temperature */
	/*
	printf("Floored it after %i brackets\n",i);
	*/
	*Y = b.Y;
	*E = b.E;
	printf("Floored(2) Eout: %g, Tout: %g ",b.E,b.T);
	return;
      }
      clRates( cl, &Rate, b.T ); 
      clAbunds( cl, &b.Y, &Rate, rho );
      
      b.F = b.E - Ein - dt * clEdotHarmonic( cl, &Ratein, &Rate, &Yin, &b.Y, rho, PdV, dt );
      if ( b.F * a.F < 0 ) break;
      if (i<5) dE *= 2;
      else {
	if (dE < 0) dE -= Ein*0.2;
	else dE = dE*2 + Ein*0.2;
      }
    }
    assert(i<=MAXBRACKET);
  }
  /*
    printf("Bracket Iterations %i\n",i);
  */

  /* Brent it - Algorithm taken from Numerical Recipes */  
  c = b;

  for (iter=1;iter<=ITMAX;iter++) {
    if ((b.F > 0.0 && c.F > 0.0) || (b.F < 0.0 && c.F < 0.0)) {
      c = a;
      e = d = b.E - a.E;
    }
    if (fabs(c.F) < fabs(b.F)) {
      a = b;
      b = c;
      c = a;
    }
    tol1 = 2*BRENTEPS*fabs(b.E) + 0.5*BRENTTOL;
    xm = 0.5 * (c.E - b.E);
    
    if (fabs(xm) <= tol1 || b.F == 0.0) {
      /*
      printf("Brent iterations: %i\n",iter);
      */
      /* finished -> b.E */
      *E = b.E;
      *Y = b.Y;
      printf("Eout: %g, Tout: %g ",b.E,b.T);
      return;
    }
    if (fabs(e) >= tol1 && fabs(a.F) > fabs(b.F)) {
      s=b.F/a.F;
      if (a.E == c.E) {
        p = 2*xm*s;
        q = 1-s;
      } else {
        q = a.F/c.F;
        r = b.F/c.F;
        p = s*(2.0*xm*q*(q-r)-(b.E-a.E)*(r-1.0));
        q = (q-1.0)*(r-1.0)*(s-1.0);
      }
      if (p > 0.0) q = -q;
      p=fabs(p);
      min1 = 3.0*xm*q - fabs(tol1*q);
      min2 = fabs(e*q);
      if (2.0*p < (min1 < min2 ? min1 : min2)) {
        e=d;
        d=p/q;
      } else {
        d=xm;
        e=d;
      }
    } else {
      d=xm;
      e=d;
    }
    a=b;
    if (fabs(d) > tol1)
      b.E += d;
    else
      b.E += (xm > 0.0 ? fabs(tol1) : -fabs(tol1));

    b.T = clTemperature( b.Y.Total, b.E ); 
    clRates( cl, &Rate, b.T );
    clAbunds( cl, &b.Y, &Rate, rho );
    b.T = clTemperature( b.Y.Total, b.E ); 
    clRates( cl, &Rate, b.T ); 
    clAbunds( cl, &b.Y, &Rate, rho );

    b.F = b.E - Ein - dt * clEdotHarmonic( cl, &Ratein, &Rate, &Yin, &b.Y, rho, PdV, dt );
  }
  assert( 0 ); /* Brent failed to converge */
}
#endif

#endif /* NOCOOLING */
#endif /* GASOLINE */
