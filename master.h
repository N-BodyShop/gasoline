#ifndef MASTER_HINCLUDED
#define MASTER_HINCLUDED

#include "param.h"
#include "pst.h"
#include "mdl.h"
#include "parameters.h"

#define MSR_INIT_ECOSMO		1
#define MSR_STEP_ECOSMO		0

#define MAX_REDSHIFTS	1000

typedef struct msrContext {
	PRM prm;
	PST pst;
	MDL mdl;
	LCL lcl;
	float fCenter[3];
	struct parameters param;
	int nThreads;
	int N;
	int iOpenType;
	double dCrit;
	/*
	 ** Comoving coordinate variables.
	 */
	double dRedshift;
	double dHubble;
	double dCosmoFac;
	double dEcosmo;
	double dHubbleOld,dUOld,dTimeOld;
	/*
	 ** Redshift output points.
	 */
	int nRed;
	double dRedOut[MAX_REDSHIFTS];
	} * MSR;



void msrInitialize(MSR *,MDL,int,char **,char *);
void msrLogParams(MSR msr, FILE *fp);
void msrFinish(MSR);
double msrReadTipsy(MSR);
void msrWriteTipsy(MSR,char *,double);
void msrSetSoft(MSR msr);
void msrBuildTree(MSR);
void msrDomainColor(MSR);
void msrReorder(MSR);
void msrOutArray(MSR,char *,int);
void msrOutVector(MSR,char *,int);
void msrDensity(MSR);
void msrGravity(MSR,int *,double *,double *,double *);
void msrCalcE(MSR,int,double,double *,double *,double *);
void msrDrift(MSR,double);
void msrKick(MSR,double);
double msrReadCheck(MSR,int *);
void msrWriteCheck(MSR,double,int);
void msrStepCosmo(MSR,double);
double msrRedOut(MSR,int);
void msrReadRed(MSR);
/*
 ** Interface functions.
 */
int msrSteps(MSR);
char *msrOutName(MSR);
double msrDelta(MSR);
double msrRedshift(MSR);
int msrLogInterval(MSR);
int msrCheckInterval(MSR);
int msrOutInterval(MSR);
int msrRestart(MSR);
int msrComove(MSR);

#endif




