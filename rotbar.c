#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "rotbar.h"
#include "romberg.h"
#include "pkd.h"

#define EPS 1e-3

void rotbarAddParams(ROTBAR rotbar, PRM prm)
{
	rotbar->bFixedBar = 0;
	prmAddParam(prm,"bFixedBar",0,&rotbar->bFixedBar, sizeof(int),
		    "fixedbar", "bar with fixed rotation speed = -fixedbar");
	rotbar->bMonopole = 1;
	prmAddParam(prm,"bBarMonopole",0,&rotbar->bMonopole,
		    sizeof(int),"barmonopole",
		    "include monopole potential of bar = +barmonopole");
	rotbar->dMass = 0;
	prmAddParam(prm,"dRotBarMass",2,&rotbar->dMass,
		    sizeof(double), "rbmass","mass of bar = 0");
	rotbar->dAmplitude = 0.0;
	prmAddParam(prm,"dRotBarAmplitude",2,&rotbar->dAmplitude,
		    sizeof(double), "rbamp","quadrapole amplitude of bar = 0");
	rotbar->dLength = 0;
	prmAddParam(prm,"dRotBarLength",2,&rotbar->dLength,
		    sizeof(double), "rblen","length of bar = 0");
	rotbar->dCorotFac = 0;
	prmAddParam(prm,"dRotBarCorotFac",2,&rotbar->dCorotFac,
		    sizeof(double), "rbcorot",
		    "fraction of corotation of bar = 0");
	rotbar->dTurnOn = 0;
	prmAddParam(prm,"dRotBarTurnOn",2,&rotbar->dTurnOn,
		    sizeof(double), "rbturnon", "turn on time of bar = 0");
	rotbar->dTurnOff = 0;
	prmAddParam(prm,"dRotBarTurnOff",2,&rotbar->dTurnOff,
		    sizeof(double), "rbturnoff", "turn off time of bar = 0");
	rotbar->dDuration = 0;
	prmAddParam(prm,"dRotBarDuration",2,&rotbar->dDuration, sizeof(double),
		    "rbduration", "off/on time duration of bar = 0");
	rotbar->bratio = 0.5;
	prmAddParam(prm,"dRotBarBratio",2,&rotbar->bratio,
		    sizeof(double), "rbbratio", "b/a ratio of bar = 0.5");
	rotbar->cratio = 0.1;
	prmAddParam(prm,"dRotBarCratio",2,&rotbar->cratio,
		    sizeof(double), "rbcratio", "c/b ratio of bar = 0.1");
    }

void rotbarLogParams(ROTBAR rotbar, FILE *fp )
{
	fprintf(fp," dRotBarMass: %g",rotbar->dMass );
	fprintf(fp," dRotBarLength: %g",rotbar->dLength );
	fprintf(fp," dRotBarCorotFac: %g",rotbar->dCorotFac );
	fprintf(fp," dRotBarTurnOff: %g",rotbar->dTurnOff );
	fprintf(fp," dRotBarTurnOn: %g",rotbar->dTurnOn );
	fprintf(fp," dRotBarAmpFac: %g",rotbar->dAmpFac );
	fprintf(fp," dRotBarDuration: %g",rotbar->dDuration );
	fprintf(fp,"\n# dRotBarPosAng: %g",rotbar->dPosAng );
	fprintf(fp," dRotBarTime0: %g",rotbar->dTime0 );
	fprintf(fp," dRotBarOmega: %g",rotbar->dOmega );
	fprintf(fp," dRotBarB5: %g",rotbar->dB5 );
	fprintf(fp," dRotBarIz: %g",rotbar->dIz );
	fprintf(fp," dRotBarLz: %g",rotbar->dLz );
	fprintf(fp," dRotBarLz0: %g",rotbar->dLz0 );
	fprintf(fp," bFixedBar: %d",rotbar->bFixedBar);
	fprintf(fp," bBarMonopole: %d",rotbar->bMonopole);
    }

void rotbarCheckWrite(ROTBAR rotbar, FDL_CTX *fdl)
{
	FDL_write(fdl,"dRotBarMass", &rotbar->dMass);
	FDL_write(fdl,"dRotBarLength", &rotbar->dLength);
	FDL_write(fdl,"dRotBarCorotFac", &rotbar->dCorotFac);
	FDL_write(fdl,"dRotBarTurnOff", &rotbar->dTurnOff);
	FDL_write(fdl,"dRotBarTurnOn", &rotbar->dTurnOn);
	FDL_write(fdl,"dRotBarAmpFac", &rotbar->dAmpFac);
	FDL_write(fdl,"dRotBarDuration", &rotbar->dDuration);
	FDL_write(fdl,"dRotBarPosAng", &rotbar->dPosAng);
	FDL_write(fdl,"dRotBarPosX", &rotbar->dPos[0]);
	FDL_write(fdl,"dRotBarPosY", &rotbar->dPos[1]);
	FDL_write(fdl,"dRotBarPosZ", &rotbar->dPos[2]);
	FDL_write(fdl,"dRotBarVelX", &rotbar->dVel[0]);
	FDL_write(fdl,"dRotBarVelY", &rotbar->dVel[1]);
	FDL_write(fdl,"dRotBarVelZ", &rotbar->dVel[2]);
	FDL_write(fdl,"dRotBarTime0", &rotbar->dTime0);
	FDL_write(fdl,"dRotBarOmega", &rotbar->dOmega);
	FDL_write(fdl,"dRotBarB5", &rotbar->dB5);
	FDL_write(fdl,"dRotBarIz", &rotbar->dIz);
	FDL_write(fdl,"dRotBarLz", &rotbar->dLz);
	FDL_write(fdl,"dRotBarLz0", &rotbar->dLz0);
    }

void rotbarCheckRead(ROTBAR rotbar, FDL_CTX *fdl)
{
    FDL_read(fdl,"dRotBarMass", &rotbar->dMass);
    FDL_read(fdl,"dRotBarLength", &rotbar->dLength);	
    FDL_read(fdl,"dRotBarCorotFac", &rotbar->dCorotFac);
    FDL_read(fdl,"dRotBarTurnOff", &rotbar->dTurnOff);
    FDL_read(fdl,"dRotBarTurnOn", &rotbar->dTurnOn);
    FDL_read(fdl,"dRotBarAmpFac", &rotbar->dAmpFac);
    FDL_read(fdl,"dRotBarDuration", &rotbar->dDuration);
    FDL_read(fdl,"dRotBarPosAng", &rotbar->dPosAng);
    FDL_read(fdl,"dRotBarPosX", &rotbar->dPos[0]);
    FDL_read(fdl,"dRotBarPosY", &rotbar->dPos[1]);
    FDL_read(fdl,"dRotBarPosZ", &rotbar->dPos[2]);
    FDL_read(fdl,"dRotBarVelX", &rotbar->dVel[0]);
    FDL_read(fdl,"dRotBarVelY", &rotbar->dVel[1]);
    FDL_read(fdl,"dRotBarVelZ", &rotbar->dVel[2]);
    FDL_read(fdl,"dRotBarTime0", &rotbar->dTime0);
    FDL_read(fdl,"dRotBarOmega", &rotbar->dOmega);
    FDL_read(fdl,"dRotBarB5", &rotbar->dB5);
    FDL_read(fdl,"dRotBarIz", &rotbar->dIz);
    FDL_read(fdl,"dRotBarLz", &rotbar->dLz);
    FDL_read(fdl,"dRotBarLz0", &rotbar->dLz0);
    }

void rotbarInitialize(ROTBAR *protbar)
{
    ROTBAR rotbar;
    
    rotbar = (ROTBAR) malloc(sizeof(struct rotbarContext));
    assert(rotbar != NULL);
    /* maybe initialize stuff here */
    *protbar = rotbar;
    }

typedef struct intContext 
{
    double a, b, c;
    double x, y, z;
    double xfac, yfac, zfac;
    } INTCTX;

double ellipden(void *vctx, double z)
{
    INTCTX *ctx = (INTCTX *)vctx;
    double x = ctx->x;
    double y = ctx->y;
    double a = ctx->a;
    double b = ctx->b;
    double c = ctx->c;

    double ans = 0.0;
    /* Is inside ellipsoid? */
    if (x*x/(a*a) + y*y/(b*b) + z*z/(c*c) < 1.0) 
	ans = 1.0;

    return ans;
    }

double ellipden2d(void *vctx, double y)
{
    INTCTX *ctx = (INTCTX *)vctx;
    double x = ctx->x;
    double a = ctx->a;
    double b = ctx->b;
    double c = ctx->c;
    double spheremax;
    double ellipsemax;
    double zmax;
    
    if( (ctx->xfac*ctx->xfac - x*x - y*y) > 0.0)
    	spheremax = sqrt(ctx->xfac*ctx->xfac - x*x - y*y);
    else 
	spheremax = 0.0;
    if((1.0 - x*x/(a*a) - y*y/(b*b)) > 0.0)
	ellipsemax = c*sqrt(1.0 - x*x/(a*a) - y*y/(b*b));
    else 
	return 0.0;
    
    zmax = (spheremax < ellipsemax ) ? spheremax : ellipsemax;
    return zmax;
    }

double ellipden3d(void *vctx, double x)
{
    INTCTX *ctx = (INTCTX *)vctx;
    
    ctx->x = x;
    if( (ctx->xfac*ctx->xfac - x*x) > 0.0)
    	ctx->yfac = sqrt(ctx->xfac*ctx->xfac - x*x);
    else
	return 0.0;
    return dRombergC(ctx, ellipden2d, 0.0, ctx->yfac, EPS);
    }

/* Linear interpolation derivative */
double linderive(double xinterp, double *x, double *y, int nMax)
{
    int i;
    assert(xinterp >= x[0]);
    for(i = 0; i < nMax-1 && xinterp > x[i];  i++);
    if(i == 0) i++;
    
    return((y[i] - y[i-1])/(x[i] - x[i-1]));
    }

/* Linear interpolation */
double lininterp(double xinterp, double *x, double *y, int nMax)
{
    int i;
    assert(xinterp >= x[0]);
    for(i = 0; (i < nMax-1) && (xinterp > x[i]);  i++);
    if(i == 0) return y[0];
    
    return(y[i-1] + (xinterp - x[i-1])*(y[i] - y[i-1])/(x[i] - x[i-1]));
    }

double rotbarGetMass(ROTBAR rotbar, double dR)
{
    if(dR >= rotbar->A) return rotbar->pdMmono[rotbar->nBinMono - 1];
    return lininterp(dR, rotbar->pdRmono, rotbar->pdMmono, rotbar->nBinMono);
    }

double rotbarGetPot(ROTBAR rotbar, double dR)
{
    if(dR >= rotbar->A) return rotbar->pdPmono[rotbar->nBinMono - 1];
    return lininterp(dR, rotbar->pdRmono, rotbar->pdPmono, rotbar->nBinMono);
    }


void rotbarInitValues(ROTBAR rotbar)
{
    const int nInterp = 100;
  /* Work vectors */
    double w1[nInterp], w2[nInterp];
    double mfac;
    double mass;
    double a, b, c;
    int numr;
    double ans;
    double dr;
    int v;
    INTCTX intctx;

    rotbar->nBinMono = nInterp;
    rotbar->pdRmono = malloc(rotbar->nBinMono*sizeof(*(rotbar->pdRmono)));
    rotbar->pdMmono = malloc(rotbar->nBinMono*sizeof(*(rotbar->pdMmono)));
    rotbar->pdPmono = malloc(rotbar->nBinMono*sizeof(*(rotbar->pdPmono)));
    
    a = rotbar->A;
    b = rotbar->B;
    c = rotbar->C;
    numr = nInterp;
    mass = rotbar->dMass;

    dr = a/numr;

    mfac = mass/(4.0*M_PI/3.0*a*b*c);

    for (v=0; v<numr; v++) {

	rotbar->pdRmono[v] = dr*v;

	intctx.a = rotbar->A;
	intctx.b = rotbar->B;
	intctx.c = rotbar->C;
	intctx.xfac = (rotbar->pdRmono[v] < a) ? rotbar->pdRmono[v] : a;
	ans = dRombergC(&intctx, ellipden3d, 0.0, intctx.xfac, EPS);
	
    rotbar->pdMmono[v] = 8.0*ans*mfac;	/* Tabulate total mass */
  }

				/* External potential integrand: (dM/dr)/r */
  for (v=0; v<numr; v++) {
    if (rotbar->pdRmono[v] <= 0.0) w1[v] = 0.0;
    else w1[v] = linderive(rotbar->pdRmono[v], rotbar->pdRmono,
			   rotbar->pdMmono, rotbar->nBinMono)
	    	/rotbar->pdRmono[v];
  }

				/* Integrate external potential */
  w2[0] = 0.0;			/* using trapezoidal rule */
  for (v=1; v<numr; v++) 
    w2[v] = 0.5*(rotbar->pdRmono[v] - rotbar->pdRmono[v-1])*(w1[v] + w1[v-1])
	+ w2[v-1];

				/* Compute the total gravitational potential */
  for (v=0; v<numr; v++) {
      if (rotbar->pdRmono[v] <= 0.0) rotbar->pdPmono[v] = -w2[numr-1];
      else rotbar->pdPmono[v] = -rotbar->pdMmono[v]/rotbar->pdRmono[v] - (w2[numr-1] - w2[v]);
      }
    
  rotbar->getMass = &rotbarGetMass;
  rotbar->getPot = &rotbarGetPot;
    }

void rotbarDrift(ROTBAR rotbar, double dTime, double dDelta)
{
    double tnow = dTime;
    int j;

    printf("rotbar: %g %g %g %g %g %g %g\n", dTime, rotbar->dPos[0],
	   rotbar->dPos[1], rotbar->dPos[2], rotbar->dPosAng, rotbar->dOmega,
	   rotbar->amplitude);
    
    if(rotbar->bMonopole) {
	for (j=0;j<3;++j) {
	    rotbar->dPos[j] += dDelta*rotbar->dVel[j];
	    }
	}
    
    rotbar->amplitude = rotbar->dAmplitude*rotbar->dAmpFac
	*0.5*(1.0 + erf((tnow - rotbar->dTurnOn)/rotbar->dDuration))
	*0.5*(1.0 - erf((tnow - rotbar->dTurnOff)/rotbar->dDuration)) ;
    rotbar->dPosAng += dDelta*rotbar->dOmega;
    }

void rotbarKick(ROTBAR rotbar, double dvFacOne, double dvFacTwo)
{
    int j;
    
    if(rotbar->bMonopole) {
	for (j=0;j<3;++j) {
	    rotbar->dVel[j] = rotbar->dVel[j]*dvFacOne
		+ rotbar->dAcc[j]*dvFacTwo;
	    }
	}

    if(!rotbar->bFixedBar) {
	/*
	rotbar->dOmega = (rotbar->dLz + rotbar->dLz0 -
			  rotbar->dLzPart)/rotbar->dIz;
	*/
	fprintf(stderr, "Torque: %g %g %g %g\n", rotbar->dTorque[0],
		rotbar->dTorque[1], rotbar->dTorque[2], rotbar->dIz);
	rotbar->dOmega += dvFacTwo*rotbar->dTorque[2]/rotbar->dIz;
	}
    }

void pkdInitRotBar(PKD pkd, ROTBAR rotbar)
{
    ROTBAR rbNode;
    
    rotbarInitialize(&rbNode);
    
    *rbNode = *rotbar;
    
    rotbarInitValues(rbNode);
    
    pkd->rotbar = rbNode;
    }

/*  U1 (0.5)	Bar length		*/
/*  U2 (0.3)	Bar amplitude		*/
/*  U3 (-20.0)	Turn on start time	*/
/*  U4 (1.0)	Turn on duration	*/
/*  U5 (0.0)	Corotation factor	*/
/*  U6 (1000.0)	Turn off start time	*/

void pkdRotatingBar(PKD pkd, double amp, /* relative amplitude of bar */
		    double posang, /* position angle of bar */
		    double b5,	/* radial scale length (^5) */
		    FLOAT *aCom, /* Center of mass */
		    double *accCom, /* acceleration (returned) */
		    double *dTorque) /* Torque (returned) */
{
  double fac, ffac;
  const double numfac = 3.86274202023190e-01; /* sqrt(15/(32*Pi)) */

  double xx, yy, zz, rr, nn, pp; 
  double cos2p = cos(2.0*posang);
  double sin2p = sin(2.0*posang);
  double acc[3];
  double pos[3];
  double dTorqueTmp[3];
  double M0;
  int soft = pkd->rotbar->soft;
  
  int i, k;

  PARTICLE *p = pkd->pStore;
  int nLocal = pkd->nLocal;
  for(k = 0; k < 3; k++) {
      accCom[k] = 0.0;
      dTorqueTmp[k] = 0.0;
      }

  for (i=0; i < nLocal; i++)
    {
	if (TYPEQueryACTIVE(&(p[i]))) {
	    for(k = 0; k < 3; k++)
		pos[k] = p[i].r[k] - aCom[k];	
	    xx = pos[0];	
	    yy = pos[1];
	    zz = pos[2];
	    rr = sqrt( xx*xx + yy*yy + zz*zz );

	    if (soft) {
	      fac = 1.0 + rr/b5;

	      ffac = -amp*numfac/pow(fac, 6.0);

	      pp = (xx*xx - yy*yy)*cos2p + 2.0*xx*yy*sin2p;
	      nn = pp /( b5*rr ) ;
	    } else {
		fac = 1.0 + pow(rr/b5, 5.0);

		ffac = -amp*numfac/(fac*fac); 

		pp = (xx*xx - yy*yy)*cos2p + 2.0*xx*yy*sin2p;
		nn = pp * pow(rr/b5, 3.0)/(b5*b5);
		}

	    acc[0] = ffac* ( 2.0*( xx*cos2p + yy*sin2p)*fac - 5.0*nn*xx );
	    acc[1] = ffac* ( 2.0*(-yy*cos2p + xx*sin2p)*fac - 5.0*nn*yy );
	    acc[2] = ffac* ( -5.0*nn*zz );
	    
	    dTorqueTmp[0] -= p[i].fMass*(pos[1] * acc[2] - pos[2]*acc[1]);
	    dTorqueTmp[1] -= p[i].fMass*(pos[2] * acc[0] - pos[0]*acc[2]);
	    dTorqueTmp[2] -= p[i].fMass*(pos[0] * acc[1] - pos[1]*acc[0]);

	    if(pkd->rotbar->bMonopole)
		M0 = pkd->rotbar->getMass(pkd->rotbar, rr);
	    else
		M0 = 0.0;

	    for (k=0; k<3; k++) {
		/* Add monopole acceleration */
		acc[k] += -M0*pos[k]/(rr*rr*rr);

	      /* Add bar acceleration to particle */
		p[i].a[k] += acc[k];
		p[i].fPot += -ffac*pp*fac
		    + pkd->rotbar->getPot(pkd->rotbar, rr);

		/* Force on bar (via Newton's 3rd law) */
		accCom[k] -= p[i].fMass*acc[k];
		}
	    }
	}
  for (k=0; k<3; k++) {
      /* divide by total mass of bar here */
      accCom[k] /= pkd->rotbar->dMass;
      dTorque[k] = dTorqueTmp[k];
      }
  }

#include "master.h"

void msrInitRotatingBar(MSR msr, double dTime)
{
    struct inMassInR in;
    struct outMassInR out;
    struct outCalcEandL outL;
    struct inRotBar inRotBar;
    ROTBAR rotbar = msr->param.rotbar;
    int iDum;
    double a1, a2, a3, geom, A12, A22, A32;
    double u, d, t, denom, ans1=0.0, ans2=0.0;
    double mass = rotbar->dMass;
    const int N = 4000;
    double dt = 1.0/N;
    double rho, b1, b25;
    int i;
    
    if(!msr->param.bRestart) {
	in.R = rotbar->dLength*rotbar->dCorotFac;
	pstMassInR(msr->pst, &in, sizeof(in), &out, &iDum);
	rotbar->dOmega = sqrt(out.dMass/pow(in.R, 3.0));
	}

    rotbar->A = a1 = rotbar->dLength;
    rotbar->B = a2 = rotbar->bratio*a1;		
    rotbar->C = a3 = rotbar->cratio*a2;
    rotbar->soft = 0;

    geom = pow(a1*a2*a3, 1.0/3.0);

    A12 = a1*a1/geom/geom;
    A22 = a2*a2/geom/geom;
    A32 = a3*a3/geom/geom;

    for (i=1; i<=N; i++) { /* simple centered rectangle
				  integration; can be replaced by
				  Romberg.
			       */
	t = 0.5*M_PI*dt*((double)i-0.5);
	u = tan(t);
	d = cos(t);
	d = 1.0/(d*d);

	denom = sqrt( (A12+u)*(A22+u)*(A32+u) );
	ans1 += d/( (A12+u)*denom );
	ans2 += d/( (A22+u)*denom );
	}
    ans1 *= 0.5*M_PI*dt;
    ans2 *= 0.5*M_PI*dt;

    printf("Computed quadrupole fit to homogenous ellipsoid\n");
    printf("with Mass=%g A_1=%g A_2=%g A_3=%g\n", mass, a1, a2, a3);
    printf("with an exact fit to asymptotic quadrupole.\n");
    printf("V_1=%g\n", ans1);
    printf("V_2=%g\n", ans2);
    printf("I_3=%g\n", 0.2*mass*(a1*a1 + a2*a2));

    rho = mass/(4.0*M_PI/3.0*a1*a2*a3);
    b1 = M_PI*rho*sqrt(2.0*M_PI/15.0)*(ans1 - ans2);
    b25 = 0.4*a1*a2*a3*(a2*a2 - a1*a1)/(ans1 - ans2);

    rotbar->dB5 = pow(b25, 0.2);
    rotbar->dAmpFac = b1;
    printf("b1=%g\n", b1);
    printf("b5=%g\n", rotbar->dB5);
    printf("afac=%g\n", rotbar->dAmpFac);

    if(!msr->param.bRestart) {
	pstCalcEandL(msr->pst, NULL, 0, &outL, &iDum);

	rotbar->dIz = 0.2*mass*(a1*a1 + a2*a2);
	rotbar->dLz = 0.2*mass*(a1*a1 + a2*a2)*rotbar->dOmega;
	rotbar->dLz0 = outL.L[2];
	rotbar->dLzPart = outL.L[2];
	rotbar->dPosAng = 0.0;
	rotbar->dTime0 = dTime;
	}
    
    rotbar->amplitude = rotbar->dAmplitude*rotbar->dAmpFac
	*0.5*(1.0 + erf((dTime - rotbar->dTurnOn)/rotbar->dDuration))
	*0.5*(1.0 - erf((dTime - rotbar->dTurnOff)/rotbar->dDuration)) ;
    inRotBar.rotbar = *rotbar;
    pstInitRotBar(msr->pst, &inRotBar, sizeof(struct inRotBar), NULL, NULL);
    }

