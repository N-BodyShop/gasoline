
#ifdef GASOLINE
#ifndef NOCOOLING

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

/* General accuracy target */
#define EPS 1e-5
#define MAXABUNDITERATIONS 20
/* Accuracy for Temperature estimation for E,rho assuming eqm abundances */
#define EPSTEMP 1e-5
#define ETAMINTIMESTEP 1e-4

#include "stiff.h"
#if defined(COOLDEBUG) || defined(STARFORM) || defined(SIMPLESF)
#include "pkd.h"
#endif
#include "cooling.h"
#include "outtype.h"

/* Integrator constants */

/* When to use to just do a first order step and quit */
/* Good value to use I think */
#define ECHANGEFRACSMALL 1e-4  

/* Accuracy target for intergrators */
#define EPSINTEG  1e-5
#define MAXINTEGITS 20000

#define EMUL (1.01)

COOL *CoolInit( )
{
  COOL *cl;
  cl = (COOL *) malloc(sizeof(COOL));
  assert(cl!=NULL);

  cl->nTableRead = 0; /* Internal Tables read from Files */

  cl->DerivsData = malloc(sizeof(clDerivsData));
  assert(cl->DerivsData != NULL);
  ((clDerivsData *) (cl->DerivsData))->IntegratorContext = 
    StiffInit( EPSINTEG, 1, cl->DerivsData, clDerivs, clJacobn );
  
  return cl;
}

void CoolFinalize(COOL *cl ) 
{
  StiffFinalize( ((clDerivsData *) (cl->DerivsData))->IntegratorContext );
  free(cl->DerivsData);
  free(cl);
}

void clInitConstants( COOL *cl, double dGmPerCcUnit, double dComovingGmPerCcUnit, 
		double dErgPerGmUnit, double dSecUnit, double dKpcUnit, COOLPARAM CoolParam) 
{
  assert(cl!=NULL);
  cl->dGmPerCcUnit = dGmPerCcUnit;
  cl->dComovingGmPerCcUnit = dComovingGmPerCcUnit;
  cl->dErgPerGmUnit = dErgPerGmUnit;
  cl->dSecUnit = dSecUnit;
  cl->dErgPerGmPerSecUnit = cl->dErgPerGmUnit / cl->dSecUnit;
  cl->diErgPerGmUnit = 1./dErgPerGmUnit;
  cl->dKpcUnit = dKpcUnit;
  
  cl->BaseT = CoolParam.BaseT;
  cl->dParam2 = CoolParam.dParam2;
  cl->dParam3 = CoolParam.dParam3;
  cl->dParam4 = CoolParam.dParam4;
  cl->Y_Total = CoolParam.Y_Total;
  cl->Tmin = CoolParam.dCoolingTmin;
  cl->Tmax = CoolParam.dCoolingTmax;

  /* Derivs Data Struct */
  {
    clDerivsData *Data = cl->DerivsData;

    assert(Data != NULL);

    Data->cl = cl;
    Data->dlnE = (log(EMUL)-log(1/EMUL));
  }
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
    assert(0);
    return Y_Total*CL_Eerg_gm_degK3_2*T;

}

double clTemperature( double Y_Total, double E ) {
    assert(0);
    return E/(Y_Total*CL_Eerg_gm_degK3_2);
}

double clEdotInstant( COOL *cl, double E, double T, double rho, double rFactor )
{
	double Edot;

	Edot = 0;

	return Edot;
	}

/*
 *  We solve the implicit equation:  
 *  Eout = Ein + dt * Cooling ( Yin, Yout, Ein, Eout )
 *
 *  E erg/gm, PdV erg/gm/s, rho gm/cc, dt sec
 * 
 */

void clDerivs(void *Data, double x, double *y, double *dydx) {
  clDerivsData *d = Data;

  d->E = y[1];
  d->T = clTemperature( d->Y_Total, d->E );
  /*
  dydx[1] = clEdotInstant( d->cl, d->E, d->T, d->rho, d->rFactor ) + d->PdV;
  */
  dydx[1] = 0.0;
}

void clJacobn(void *Data, double x, double y[], double dfdx[], double **dfdy) {
  clDerivsData *d = Data;
  double E = y[1],dE;

  dfdx[1] = 0;

  /* Approximate dEdt/dE */
  d->E = E*(EMUL);
  d->T = clTemperature( d->Y_Total, d->E );
  dE = clEdotInstant( d->cl, d->E, d->T, d->rho, d->rFactor );

  d->E = E*(1/EMUL);
  d->T = clTemperature( d->Y_Total, d->E );
  dE -= clEdotInstant( d->cl, d->E, d->T, d->rho, d->rFactor );

  dfdy[1][1] = dE/(E*d->dlnE);
}

void clIntegrateEnergy(COOL *cl, double *E, 
		       double PdV, double rho, double Y_Total, double radius, double tStep ) {

  double dEdt,dE,Ein = *E,EMin;
  double t=0,dtused,dtnext,tstop = tStep*(1-1e-8),dtEst;
  clDerivsData *d = cl->DerivsData;
  STIFF *sbs = d->IntegratorContext;
  
  if (tStep <= 0) return;

  d->rho = rho;
  d->PdV = PdV;
  d->Y_Total = Y_Total;
/*  d->rFactor = cl->dParam1*pow(radius,-3./2.);*/
  
  EMin = clThermalEnergy( d->Y_Total, cl->Tmin );

  dtnext = tStep;
  {
    int its = 0;
    while (t<tstop) {
      double Eold;
      if (its++ > MAXINTEGITS) assert(0);
      d->E = *E;
	  d->T = clTemperature( d->Y_Total, d->E );
      clDerivs( d, t, E-1, (&dEdt)-1 );
      if (fabs(dEdt) > 0) {
		  dtEst = fabs(*E/dEdt);

      /* 
	 Since there is no time dependence and the function is smooth
	 if the changes become very small it must have reached a saddle or
	 an equilibrium.  I'll put my money on Equilibrium and abort 
      */
      /*
      if (tStep-t < ECHANGEFRACSMALL*dtEst) {
	fprintf(stderr,"Aborting -- changes too small\n");
	*E += (tStep-t)*dEdt;
	break;
      }
      if (dtEst < tStep-t) sbs->hMin = dtEst*ETAMINTIMESTEP;
      else sbs->hMin = (tStep-t)*ETAMINTIMESTEP;
      */

	if (dtnext > 0.5*dtEst) dtnext = 0.5*dtEst;
      }
      if (dtnext >= tStep-t) dtnext = tStep-t;
      StiffStep( sbs, (E-1), (&dEdt)-1,  &t, dtnext, (&Ein)-1, &dtused, &dtnext );
      Eold = *E;
#ifdef ASSERTENEG      
      assert(*E > 0);
#else
      if (*E < EMin) {
	*E = EMin;
	break;
      }
#endif    
    }
  }
  /* 
     Note Stored abundances are not necessary with
     this eqm integrator therefore the following operations
     could be moved to output routines 
  */

  d->E = *E;
  d->T = clTemperature( d->Y_Total, d->E );
  if (d->T < cl->Tmin ) {
	  d->T = cl->Tmin;
	  *E = clThermalEnergy( d->Y_Total, d->T );
	  }
}

/* Module Interface routines */

void CoolAddParams( COOLPARAM *CoolParam, PRM prm ) {
	CoolParam->BaseT = 10;
	prmAddParam(prm,"CoolBaseT",2,&CoolParam->BaseT,
				sizeof(double),"BaseT",
				"<Param1> = 10.0");
	CoolParam->dParam2 = 0;
	prmAddParam(prm,"CooldParam2",2,&CoolParam->dParam2,
				sizeof(double),"dParam2",
				"<Param2> = 0.0");
	CoolParam->dParam3 = 0;
	prmAddParam(prm,"CooldParam3",2,&CoolParam->dParam3,
				sizeof(double),"dParam3",
				"<Param3> = 0.0");
	CoolParam->dParam4 = 0;
	prmAddParam(prm,"CooldParam4",2,&CoolParam->dParam4,
				sizeof(double),"dParam4",
				"<Param4> = 0.0");
	CoolParam->Y_Total = 0.5;
	prmAddParam(prm,"dY_Total",0.5,&CoolParam->Y_Total,
				sizeof(double),"Y_Total",
                "<Y_Total> = 0.4166667");  /* Solar metallicity, H2 Y_T=1/2.4 */
	CoolParam->dCoolingTmin = 10;
	prmAddParam(prm,"dCoolingTmin",2,&CoolParam->dCoolingTmin,
				sizeof(double),"ctmin",
				"<Minimum Temperature for Cooling> = 10K");
	CoolParam->dCoolingTmax = 1e9;
	prmAddParam(prm,"dCoolingTmax",2,&CoolParam->dCoolingTmax,
				sizeof(double),"ctmax",
				"<Maximum Temperature for Cooling> = 1e9K");
	}
	
void CoolLogParams( COOLPARAM *CoolParam, LOGGER *lgr) {
    LogParams(lgr, "COOLING", "CoolBaseT: %g",CoolParam->BaseT); 
    LogParams(lgr, "COOLING", "CooldParam2: %g",CoolParam->dParam2); 
    LogParams(lgr, "COOLING", "CooldParam3: %g",CoolParam->dParam3); 
    LogParams(lgr, "COOLING", "ColldParam4: %g",CoolParam->dParam4); 
    LogParams(lgr, "COOLING", "Y_Total: %g",CoolParam->Y_Total); 
    LogParams(lgr, "COOLING", "dCoolingTmin: %g",CoolParam->dCoolingTmin); 
    LogParams(lgr, "COOLING", "dCoolingTmax: %g",CoolParam->dCoolingTmax); 
#ifdef MODBATEPOLY
    LogParams(lgr, "COOLING", " Polytrope RHOMIN %g",RHOMIN); 
#endif
	}

void CoolOutputArray( COOLPARAM *CoolParam, int cnt, int *type, char *suffix ) {
	*type = OUT_NULL;

	}

/* Initialization Routines */

void CoolTableReadInfo( COOLPARAM *CoolParam, int cntTable, int *nTableColumns, char *suffix )
{
   int localcntTable = 0;

   *nTableColumns = 0;
   }

void CoolTableRead( COOL *Cool, int nData, void *vData)
{
   fprintf(stderr," Attempt to initialize non-exitent table in cooling\n");
   assert(0);

   }

void CoolDefaultParticleData( COOLPARTICLE *cp )
{
	cp->Y_Total = 0.5;
	}

void CoolInitEnergyAndParticleData( COOL *cl, COOLPARTICLE *cp, double *E, double dDensity, double dTemp, double ZMetal )
{
	cp->Y_Total = cl->Y_Total;
#ifdef WOLFIRE
        {
        double u,PonRho,T;
        WolfirePressureEnergySoundSpeed(dDensity*cl->dGmPerCcUnit,&u,&PonRho,&T);
        *E = CoolErgPerGmToCodeEnergy(cl, u);
        }
#else
	*E = clThermalEnergy(cp->Y_Total,dTemp)*cl->diErgPerGmUnit;
#endif
	}

void CoolInitRatesTable( COOL *cl, COOLPARAM CoolParam ) {
	}

void CoolSetTime( COOL *cl, double dTime, double z ) {
	cl->z = z;
	cl->dTime = dTime;
	}

/* Output Conversion Routines */

double CoolEnergyToTemperature( COOL *Cool, COOLPARTICLE *cp, double E, double rho ) {
#ifdef WOLFIRE
    double u,PonRho,T;
    WolfirePressureEnergySoundSpeed(rho,&u,&PonRho,&T);
    return T;
#else
    assert(0);
#endif
    }

double CoolCodeEnergyToTemperature( COOL *Cool, COOLPARTICLE *cp, double E, double rho, double ZMetal ) {
	return CoolEnergyToTemperature( Cool, cp, E*Cool->dErgPerGmUnit, rho*Cool->dGmPerCcUnit );
	}

#ifdef WOLFIRE
/* Solar metallicity Y_H = 0.676 Y_He = 0.31 Y_Z = 0.014 => mmw = 2.4 for all H2
   nH is hydrogen nuclei = rho*0.676/mH
   T(nH) Based on Wolfire et al. 2003 R = 8.5 kpc T(nH)   
   fH2(nH) Based on Gnedin et al 2009 (solar case)
   n = nH*(1-0.5*fH2(nH)+(0.31*0.24+0.014/12.)/0.676)
   P = n k_B T(nH)
   double dTdn = -1e4*2*nH/(0.8*0.8)/(iT1*iT1) - 50*(1/30.)/(iT2*iT2);  
   double dfH2dn = 2*30.*30./(nH*nH*nH)*fH2*fH2;                        
   approximate cs max as sqrt(dPdrho*max(gamma_eff)*1.4)  max(gamma_eff) = 1
   add in factor of 1.4 to help code stability since c_s is just for timesteps
   Note: dP/drho is negative (min ~ -0.75) near nH ~ 1-10 H/cc anyway -- useless
*/
void WolfirePressureEnergySoundSpeed(double rho,double *u,double *PonRho,double *T) 
    {                                                                   
    double nB=rho/1.672e-24;                                           
    double nH=nB*0.676;                                                 
    double iT1 = (nH*nH/(0.8*0.8)+1), iT2 = (nH/30.+1);                 
    double Temp = 1e4/iT1 + 50/iT2 + 10;                                   
    double fH2 = 1/(1+(30.*30.)/(nH*nH));                               
    double nmono = nB*((1-fH2)*0.676+0.31*0.25 + 0.014/12.);            
    double ndi = nB*fH2*0.676*0.5;                                      
    double kB = 1.38066e-16,Epergm;                                           
#ifdef WOLFIRE_NOTWOPHASE
    if (nH > .8093477708 && (nmono+ndi)*Temp < 4516.73) Temp=4516.73/(nmono+ndi);
#endif
    Epergm = kB*Temp/rho;                                         
    *T = Temp;                                                             
    *u = (1.5*nmono+2.5*ndi)*Epergm;                                    
    *PonRho = (nmono+ndi)*Epergm;                                       
    }

void CoolCodePressureOnDensitySoundSpeed( COOL *cl, COOLPARTICLE *cp, double uPred, double fDensity, double gamma, double gammam1, double *PoverRho, double *c ) {
    double rho = (fDensity)*(cl)->dGmPerCcUnit;                 
    double uCGS, PonRhoCGS, T, gammaEff = 1.4; 
                 
    WolfirePressureEnergySoundSpeed(rho,&uCGS,&PonRhoCGS,&T);            
    *PoverRho = PonRhoCGS*(cl)->diErgPerGmUnit;                 
    *c = sqrt((gammaEff)*(*(PoverRho))); 
    } 

#endif
#ifdef MODBATEPOLY
/* Taken from Bate, ApJ, 508, L95 (1998) 
   rho = powerlaw with 4 gamma values 
   PoverRho has units of ergs per gm (same as speed squared)
*/
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

#endif
/*default*/
#ifdef BATEPOLY
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

/* Integration Routines */

#define CONVERT_CMPERKPC (3.0857e21)

void CoolIntegrateEnergyCode(COOL *cl, COOLPARTICLE *cp, double *ECode, 
		       double PdVCode, double rhoCode, double ZMetal, double *posCode, double tStep ) {
	double radius;

	double rho = rhoCode*cl->dGmPerCcUnit, Ephys; 
	double gammaEffco; 
#ifdef WOLFIRE
	double PonRho,T; 
	WolfirePressureEnergySoundSpeed(rho,&Ephys,&PonRho,&T);
	*ECode = CoolErgPerGmToCodeEnergy(cl, Ephys);
#endif
#ifdef MODBATEPOLY
	double PonRho; 
	GetGammaEff(rho,gammaEff,PonRho);
	PonRho *= cl->BaseT;
	
	Ephys = PonRho*(7./5.-1); /* cgs Erg per g */

	*ECode = CoolErgPerGmToCodeEnergy(cl, Ephys);
#endif
#ifdef BATEPOLY
	double KEff,gammaEff;
	GetGammaEff(rho,gammaEff,KEff);
	KEff *= cl->BaseT;
	
/* Bate specifies P(rho),  e = P/rho /(gamma-1) */
	if (gammaEff < 1.66) 
	    Ephys = KEff*pow(rho,gammaEff-1.)/(7./5.-1);
	else
	    Ephys = KEff*pow(rho,gammaEff-1.)/(5./3.-1);

	*ECode = CoolErgPerGmToCodeEnergy(cl, Ephys);
#endif
	/*
	radius= sqrt(posCode[0]*posCode[0]+posCode[1]*posCode[1]+posCode[2]*posCode[2])
		*cl->dKpcUnit*CONVERT_CMPERKPC;

	*ECode = CoolCodeEnergyToErgPerGm( cl, *ECode );
	clIntegrateEnergy(cl,  ECode, CoolCodeWorkToErgPerGmPerSec( cl, PdVCode ), 
					  CodeDensityToComovingGmPerCc(cl, rhoCode), cp->Y_Total, radius, tStep);
	*ECode = CoolErgPerGmToCodeEnergy(cl, *ECode);
	*/

	}

/* Star form function -- do not use */
double CoolHeatingRate( COOL *cl, COOLPARTICLE *cp, double T, double dDensity ) {
	assert(0);
	return 1.0;
	}

/* Not implemented */
double CoolEdotInstantCode(COOL *cl, COOLPARTICLE *cp, double ECode, 
			  double rhoCode, double ZMetal, double *posCode ) {
    double T,E,rho,Edot;

    return 0; /*CoolErgPerGmPerSecToCodeWork( cl, 0. );*/
    }

#endif /* NOCOOLING */
#endif /* GASOLINE */

