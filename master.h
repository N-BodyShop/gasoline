#ifndef MASTER_HINCLUDED
#define MASTER_HINCLUDED

#include "param.h"
#include "pst.h"
#include "mdl.h"
#include "parameters.h"

#define MSR_INIT_ECOSMO		1
#define MSR_STEP_ECOSMO		0

/*
 ** An integer marking the type of tree currently in use.
 ** MSR_TREE_NONE: undefined tree type
 ** MSR_TREE_SPATIAL: spatial binary tree
 ** MSR_TREE_DENSITY: density binary tree (the old style KD-tree!)
 */
#define MSR_TREE_NONE		0
#define MSR_TREE_SPATIAL	1
#define MSR_TREE_DENSITY	2

typedef struct msrContext {
	PRM prm;
	PST pst;
	MDL mdl;
	LCL lcl;
	float fCenter[3];
	/*
	 ** Parameters.
	 */
	struct parameters param;	   
	/*
	 ** Other stuff...
	 */
	int nThreads;
	int N;
	int iCurrMaxRung;
	int bOpenSpec;	/* was an opening parameter specified (used by +restart) */
	int iOpenType;
	double dCrit;
	/*
	 ** Comoving coordinate variables.
	 */
	double dEcosmo;
	double dUOld;
	double dTimeOld;
	/*
	 ** Redshift output points.
	 */
	int nMaxOuts;
	int nOuts;
	double *pdOutTime;
	int iOut;
	/*
	 ** Processor mapping for one-node-output functions.
	 */
	int *pMap;
	/*
	 ** An integer marking the type of tree currently in use.
	 ** MSR_TREE_NONE: undefined tree type
	 ** MSR_TREE_SPATIAL: spatial binary tree
	 ** MSR_TREE_DENSITY: density binary tree (the old style KD-tree!)
	 */
	int iTreeType;
	} * MSR;



void msrInitialize(MSR *,MDL,int,char **,char *);
void msrLogParams(MSR msr, FILE *fp);
void msrFinish(MSR);
double msrTime2Exp(MSR,double);
double msrExp2Time(MSR,double);
double msrTime2Hub(MSR,double);
double msrComoveDriftFac(MSR,double,double);
double msrComoveKickFac(MSR,double,double);
double msrReadTipsy(MSR);
void msrWriteTipsy(MSR,char *,double);
void msrSetSoft(MSR msr,double);
void msrBuildTree(MSR,int,double,int);
void msrDomainColor(MSR);
void msrReorder(MSR);
void msrOutArray(MSR,char *,int);
void msrOutVector(MSR,char *,int);
void msrDensity(MSR);
void msrGravity(MSR,double,int *,double *,double *,double *,int *);
void msrCalcE(MSR,int,double,double *,double *,double *);
void msrDrift(MSR,double,double);
void msrKick(MSR,double,double);
double msrReadCheck(MSR,int *);
void msrWriteCheck(MSR,double,int);
int msrOutTime(MSR,double);
void msrReadOuts(MSR,double);
double msrMassCheck(MSR,double,char *);
void msrTopStep(MSR msr, double dStep, double dTime, double dDelta, 
				double *pdMultiEff);
void msrRungStats(MSR);
void msrActiveRung(MSR msr, int iRung, int bGreater);
void msrVelocityRung(MSR msr, int iRung, double dDelta, double dTime,
		     int bAll);
/*
 ** Interface functions.
 */
int msrSteps(MSR);
char *msrOutName(MSR);
double msrDelta(MSR);
int msrLogInterval(MSR);
int msrCheckInterval(MSR);
int msrOutInterval(MSR);
int msrRestart(MSR);
int msrComove(MSR);
int msrKDK(MSR);
double msrSoft(MSR);
int msrDoDensity(MSR);

#endif







