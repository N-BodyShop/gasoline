#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <malloc.h>
#include <assert.h>
#include <time.h>
#include <math.h>

#include <rpc/types.h>
#include <rpc/xdr.h>

#ifdef CRAY_T3D
#include "hyperlib.h"
#endif

#include "master.h"
#include "tipsydefs.h"
#include "opentype.h"
#include "fdl.h"
#include "outtype.h"
#include "smoothfcn.h"

void _msrLeader(void)
{
#ifdef GASOLINE
    puts("USAGE: gasoline [SETTINGS | FLAGS] [SIM_FILE]");
#else
    puts("USAGE: pkdgrav [SETTINGS | FLAGS] [SIM_FILE]");
#endif
    puts("SIM_FILE: Configuration file of a particular simulation, which");
    puts("          includes desired settings and relevant input and");
    puts("          output files. Settings specified in this file override");
    puts("          the default settings.");
    puts("SETTINGS");
    puts("or FLAGS: Command line settings or flags for a simulation which");
    puts("          will override any defaults and any settings or flags");
    puts("          specified in the SIM_FILE.");
	}


void _msrTrailer(void)
{
	puts("(see man page for more information)");
	}


void _msrExit(MSR msr)
{
	MDL mdl=msr->mdl;
	msrFinish(msr);
	mdlFinish(mdl);
	exit(1);
	}


void msrInitialize(MSR *pmsr,MDL mdl,int argc,char **argv)
{
	MSR msr;
	int j,ret;
	int id;
	struct inSetAdd inAdd;
	struct inLevelize inLvl;
	struct inGetMap inGM;

	msr = (MSR)malloc(sizeof(struct msrContext));
	assert(msr != NULL);
	msr->mdl = mdl;
	msr->pst = NULL;
	*pmsr = msr;
	/*
	 ** Now setup for the input parameters.
	 */
	prmInitialize(&msr->prm,_msrLeader,_msrTrailer);
	msr->param.nThreads = 1;
	prmAddParam(msr->prm,"nThreads",1,&msr->param.nThreads,sizeof(int),"sz",
				"<nThreads>");
	msr->param.bDiag = 0;
	prmAddParam(msr->prm,"bDiag",0,&msr->param.bDiag,sizeof(int),"d",
				"enable/disable per thread diagnostic output");
	msr->param.bVerbose = 0;
	prmAddParam(msr->prm,"bVerbose",0,&msr->param.bVerbose,sizeof(int),"v",
				"enable/disable verbose output");
	msr->param.bPeriodic = 0;
	prmAddParam(msr->prm,"bPeriodic",0,&msr->param.bPeriodic,sizeof(int),"p",
				"periodic/non-periodic = -p");
	msr->param.bRestart = 0;
	prmAddParam(msr->prm,"bRestart",0,&msr->param.bRestart,sizeof(int),"restart",
				"restart from checkpoint");
	msr->param.bParaRead = 1;
	prmAddParam(msr->prm,"bParaRead",0,&msr->param.bParaRead,sizeof(int),"par",
				"enable/disable parallel reading of files = +par");
	msr->param.bParaWrite = 1;
	prmAddParam(msr->prm,"bParaWrite",0,&msr->param.bParaWrite,sizeof(int),"paw",
				"enable/disable parallel writing of files = +paw");
	msr->param.bCannonical = 1;
	prmAddParam(msr->prm,"bCannonical",0,&msr->param.bCannonical,sizeof(int),"can",
				"enable/disable use of cannonical momentum = +can");
	msr->param.bKDK = 1;
	prmAddParam(msr->prm,"bKDK",0,&msr->param.bKDK,sizeof(int),"kdk",
				"enable/disable use of kick-drift-kick integration = +kdk");
	msr->param.bBinary = 1;
	prmAddParam(msr->prm,"bBinary",0,&msr->param.bBinary,sizeof(int),"bb",
				"spatial/density binary trees = +bb");
	msr->param.bDoDensity = 1;
	prmAddParam(msr->prm,"bDoDensity",0,&msr->param.bDoDensity,sizeof(int),
				"den","enable/disable density outputs = +den");
	msr->param.nBucket = 8;
	prmAddParam(msr->prm,"nBucket",1,&msr->param.nBucket,sizeof(int),"b",
				"<max number of particles in a bucket> = 8");
	msr->param.nSteps = 0;
	prmAddParam(msr->prm,"nSteps",1,&msr->param.nSteps,sizeof(int),"n",
				"<number of timesteps> = 0");
	msr->param.iOutInterval = 0;
	prmAddParam(msr->prm,"iOutInterval",1,&msr->param.iOutInterval,sizeof(int),
				"oi","<number of timesteps between snapshots> = 0");
	msr->param.iLogInterval = 10;
	prmAddParam(msr->prm,"iLogInterval",1,&msr->param.iLogInterval,sizeof(int),
				"ol","<number of timesteps between logfile outputs> = 10");
	msr->param.iCheckInterval = 10;
	prmAddParam(msr->prm,"iCheckInterval",1,&msr->param.iCheckInterval,sizeof(int),
				"oc","<number of timesteps between checkpoints> = 10");
	msr->param.iOrder = 4;
	prmAddParam(msr->prm,"iOrder",1,&msr->param.iOrder,sizeof(int),"or",
				"<multipole expansion order: 1, 2, 3 or 4> = 4");
	msr->param.iEwOrder = 4;
	prmAddParam(msr->prm,"iEwOrder",1,&msr->param.iEwOrder,sizeof(int),"ewo",
				"<Ewald multipole expansion order: 1, 2, 3 or 4> = 4");
	msr->param.nReplicas = 0;
	prmAddParam(msr->prm,"nReplicas",1,&msr->param.nReplicas,sizeof(int),"nrep",
				"<nReplicas> = 0 for -p, or 1 for +p");
	msr->param.dSoft = 0.0;
	prmAddParam(msr->prm,"dSoft",2,&msr->param.dSoft,sizeof(double),"e",
				"<gravitational softening length> = 0.0");
	msr->param.dDelta = 0.0;
	prmAddParam(msr->prm,"dDelta",2,&msr->param.dDelta,sizeof(double),"dt",
				"<time step>");
	msr->param.dEta = 0.1;
	prmAddParam(msr->prm,"dEta",2,&msr->param.dEta,sizeof(double),"eta",
				"<time step criterion>");
	msr->param.bEpsVel = 1;
	prmAddParam(msr->prm,"bEpsVel",0,&msr->param.bEpsVel,sizeof(int),
		    "ev", "<Epsilon on V timestepping>");
	msr->param.bNonSymp = 0;
	prmAddParam(msr->prm,"bNonSymp",0,&msr->param.bNonSymp,sizeof(int),
		    "ns", "<Non-symplectic density stepping>");
	msr->param.iMaxRung = 1;
	prmAddParam(msr->prm,"iMaxRung",1,&msr->param.iMaxRung,sizeof(int),
		    "mrung", "<maximum timestep rung>");
	msr->param.dEwCut = 2.6;
	prmAddParam(msr->prm,"dEwCut",2,&msr->param.dEwCut,sizeof(double),"ew",
				"<dEwCut> = 2.6");
	msr->param.dEwhCut = 2.8;
	prmAddParam(msr->prm,"dEwhCut",2,&msr->param.dEwhCut,sizeof(double),"ewh",
				"<dEwhCut> = 2.8");
	msr->param.dTheta = 0.8;
	prmAddParam(msr->prm,"dTheta",2,&msr->param.dTheta,sizeof(double),"theta",
				"<Barnes opening criterion> = 0.8");
	msr->param.dAbsPartial = 0.0;
	prmAddParam(msr->prm,"dAbsPartial",2,&msr->param.dAbsPartial,sizeof(double),"ap",
				"<absolute partial error opening criterion>");
	msr->param.dRelPartial = 0.0;
	prmAddParam(msr->prm,"dRelPartial",2,&msr->param.dRelPartial,sizeof(double),"rp",
				"<relative partial error opening criterion>");
	msr->param.dAbsTotal = 0.0;
	prmAddParam(msr->prm,"dAbsTotal",2,&msr->param.dAbsTotal,sizeof(double),"at",
				"<absolute total error opening criterion>");
	msr->param.dRelTotal = 0.0;
	prmAddParam(msr->prm,"dRelTotal",2,&msr->param.dRelTotal,sizeof(double),"rt",
				"<relative total error opening criterion>");
	msr->param.dPeriod = 1.0;
	prmAddParam(msr->prm,"dPeriod",2,&msr->param.dPeriod,sizeof(double),"L",
				"<periodic box length> = 1.0");
	msr->param.achInFile[0] = 0;
	prmAddParam(msr->prm,"achInFile",3,msr->param.achInFile,256,"I",
				"<input file name> (file in TIPSY binary format)");
#ifdef GASOLINE
	strcpy(msr->param.achOutName,"gasoline");
	prmAddParam(msr->prm,"achOutName",3,&msr->param.achOutName,256,"o",
				"<output name for snapshots and logfile> = \"gasoline\"");
#else
	strcpy(msr->param.achOutName,"pkdgrav");
	prmAddParam(msr->prm,"achOutName",3,&msr->param.achOutName,256,"o",
				"<output name for snapshots and logfile> = \"pkdgrav\"");
#endif
	msr->param.bComove = 0;
	prmAddParam(msr->prm,"bComove",0,&msr->param.bComove,sizeof(int),"cm",
				"enable/disable comoving coordinates = -cm");
	msr->param.dHubble0 = 0.0;
	prmAddParam(msr->prm,"dHubble0",2,&msr->param.dHubble0,sizeof(double),"Hub",
				"<dHubble0> = 0.0");
	msr->param.dOmega0 = 1.0;
	prmAddParam(msr->prm,"dOmega0",2,&msr->param.dOmega0,sizeof(double),"Om",
				"<dOmega0> = 1.0");
	strcpy(msr->param.achDataSubPath,".");
	prmAddParam(msr->prm,"achDataSubPath",3,&msr->param.achDataSubPath,256,
				NULL,NULL);
	msr->param.dExtraStore = 0.1;
	prmAddParam(msr->prm,"dExtraStore",2,&msr->param.dExtraStore,
				sizeof(double),NULL,NULL);
	msr->param.nSmooth = 64;
	prmAddParam(msr->prm,"nSmooth",1,&msr->param.nSmooth,sizeof(int),"s",
				"<number of particles to smooth over> = 64");
	msr->param.bStandard = 0;
	prmAddParam(msr->prm,"bStandard",0,&msr->param.bStandard,sizeof(int),"std",
				"output in standard TIPSY binary format = -std");
	msr->param.dRedTo = 0.0;	
	prmAddParam(msr->prm,"dRedTo",2,&msr->param.dRedTo,sizeof(double),"zto",
				"specifies final redshift for the simulation");
	msr->param.nSuperCool = 0;
	prmAddParam(msr->prm,"nSuperCool",1,&msr->param.nSuperCool,sizeof(int),
				"scn","<number of supercooled particles> = 0");
	msr->param.dCoolFac = 0.95;
	prmAddParam(msr->prm,"dCoolFac",2,&msr->param.dCoolFac,sizeof(double),
				"scf","<Velocity Cooling factor> = 0.95 (no cool = 1.0)");
	msr->param.dCoolDens = 50.0;
	prmAddParam(msr->prm,"dCoolDens",2,&msr->param.dCoolDens,sizeof(double),
				"scd","<Velocity Cooling Critical Density> = 50");
	/*
	 ** Set the box center to (0,0,0) for now!
	 */
	for (j=0;j<3;++j) msr->fCenter[j] = 0.0;
	/*
	 ** Define any "LOCAL" parameters (LCL)
	 */
	msr->lcl.pszDataPath = getenv("PTOOLS_DATA_PATH");
	/*
	 ** Process command line arguments.
	 */
	ret = prmArgProc(msr->prm,argc,argv);
	if (!ret) {
		_msrExit(msr);
		}
	if (msr->param.bPeriodic && !prmSpecified(msr->prm,"nReplicas")) {
		msr->param.nReplicas = 1;
		}
	if (!msr->param.achInFile[0] && !msr->param.bRestart) {
		printf("ERROR: no input file specified\n");
		_msrExit(msr);
		}
	/*
	 ** Should we have restarted, maybe?
	 */
	if (!msr->param.bRestart) {
		}
	msr->nThreads = mdlThreads(mdl);
	/*
	 ** Determine opening type.
	 */
	msr->iOpenType = 0;
	msr->bOpenSpec = 1;
	if (prmFileSpecified(msr->prm,"dAbsPartial")) msr->iOpenType = OPEN_ABSPAR;
	if (prmFileSpecified(msr->prm,"dRelPartial")) msr->iOpenType = OPEN_RELPAR;
	if (prmFileSpecified(msr->prm,"dAbsTotal")) msr->iOpenType = OPEN_ABSTOT;
	if (prmFileSpecified(msr->prm,"dRelTotal")) msr->iOpenType = OPEN_RELTOT;
	if (prmArgSpecified(msr->prm,"dTheta")) msr->iOpenType = OPEN_JOSH;
	if (prmArgSpecified(msr->prm,"dAbsPartial")) msr->iOpenType = OPEN_ABSPAR;
	if (prmArgSpecified(msr->prm,"dRelPartial")) msr->iOpenType = OPEN_RELPAR;
	if (prmArgSpecified(msr->prm,"dAbsTotal")) msr->iOpenType = OPEN_ABSTOT;
	if (prmArgSpecified(msr->prm,"dRelTotal")) msr->iOpenType = OPEN_RELTOT;
	if (!msr->iOpenType) {
		msr->iOpenType = OPEN_JOSH;
		msr->bOpenSpec = 0;
		}
	switch (msr->iOpenType) {
	case OPEN_JOSH:
		msr->dCrit = msr->param.dTheta;
		break;
	case OPEN_ABSPAR:
		msr->dCrit = msr->param.dAbsPartial;
		break;
	case OPEN_RELPAR:
		msr->dCrit = msr->param.dRelPartial;
		break;
	case OPEN_ABSTOT:
		msr->dCrit = msr->param.dAbsTotal;
		break;
	case OPEN_RELTOT:
		msr->dCrit = msr->param.dRelTotal;
		break;
	default:
		msr->dCrit = msr->param.dTheta;
		}
	/*
	 ** Initialize comove variables.
	 */
	msr->nMaxOuts = 100;
	msr->pdOutTime = malloc(msr->nMaxOuts*sizeof(double));
	assert(msr->pdOutTime != NULL);
	msr->nOuts = 0;
        
	pstInitialize(&msr->pst,msr->mdl,&msr->lcl);

	pstAddServices(msr->pst,msr->mdl);
	/*
	 ** Create the processor subset tree.
	 */
	for (id=1;id<msr->nThreads;++id) {
		if (msr->param.bVerbose) printf("Adding %d to the pst\n",id);
		inAdd.id = id;
		pstSetAdd(msr->pst,&inAdd,sizeof(inAdd),NULL,NULL);
		}
	if (msr->param.bVerbose) printf("\n");
	/*
	 ** Levelize the PST.
	 */
	inLvl.iLvl = 0;
	pstLevelize(msr->pst,&inLvl,sizeof(inLvl),NULL,NULL);
	/*
	 ** Create the processor mapping array for the one-node output
	 ** routines.
	 */
	msr->pMap = malloc(msr->nThreads*sizeof(int));
	assert(msr->pMap != NULL);
	inGM.nStart = 0;
	pstGetMap(msr->pst,&inGM,sizeof(inGM),msr->pMap,NULL);
	/*
	 ** Initialize tree type to none.
	 */
	msr->iTreeType = MSR_TREE_NONE;
	msr->iCurrMaxRung = 0;
	}


void msrLogParams(MSR msr,FILE *fp)
{
	double z;
	int i;
  	
	fprintf(fp,"# N: %d",msr->N);
	fprintf(fp," nThreads: %d",msr->param.nThreads);
	fprintf(fp," bDiag: %d",msr->param.bDiag);
	fprintf(fp," bVerbose: %d",msr->param.bVerbose);
	fprintf(fp,"\n# bPeriodic: %d",msr->param.bPeriodic);
	fprintf(fp," bRestart: %d",msr->param.bRestart);
	fprintf(fp," bComove: %d",msr->param.bComove);
	fprintf(fp,"\n# bParaRead: %d",msr->param.bParaRead);
	fprintf(fp," bParaWrite: %d",msr->param.bParaWrite);
	fprintf(fp," bCannonical: %d",msr->param.bCannonical);
	fprintf(fp," bStandard: %d",msr->param.bStandard);
	fprintf(fp,"\n# bKDK: %d",msr->param.bKDK);
	fprintf(fp," nBucket: %d",msr->param.nBucket);
	fprintf(fp," iOutInterval: %d",msr->param.iOutInterval);
	fprintf(fp," iLogInterval: %d",msr->param.iLogInterval);
	fprintf(fp,"\n# iCheckInterval: %d",msr->param.iCheckInterval);
	fprintf(fp," iOrder: %d",msr->param.iOrder);
	fprintf(fp," iEwOrder: %d",msr->param.iEwOrder);
	fprintf(fp," nReplicas: %d",msr->param.nReplicas);
	fprintf(fp,"\n# dEwCut: %f",msr->param.dEwCut);
	fprintf(fp," dEwhCut: %f",msr->param.dEwhCut);
	fprintf(fp,"\n# nSteps: %d",msr->param.nSteps);
	fprintf(fp," nSmooth: %d",msr->param.nSmooth);
	fprintf(fp," dExtraStore: %f",msr->param.dExtraStore);
	if (prmSpecified(msr->prm,"dSoft"))
		fprintf(fp," dSoft: %g",msr->param.dSoft);
	else
		fprintf(fp," dSoft: input");
	fprintf(fp,"\n# dDelta: %g",msr->param.dDelta);
	fprintf(fp," dEta: %g",msr->param.dEta);
	fprintf(fp," iMaxRung: %d",msr->param.iMaxRung);
	fprintf(fp," bEpsVel: %d",msr->param.bEpsVel);
	fprintf(fp," bNonSymp: %d",msr->param.bNonSymp);
	switch (msr->iOpenType) {
	case OPEN_JOSH:
		fprintf(fp,"\n# iOpenType: JOSH");
		break;
	case OPEN_ABSPAR:
		fprintf(fp,"\n# iOpenType: ABSPAR");
		break;
	case OPEN_RELPAR:
		fprintf(fp,"\n# iOpenType: RELPAR");
		break;
	case OPEN_ABSTOT:
		fprintf(fp,"\n# iOpenType: ABSTOT");
		break;
	case OPEN_RELTOT:
		fprintf(fp,"\n# iOpenType: RELTOT");
		break;
	default:
		fprintf(fp,"\n# iOpenType: NONE?");
		}
	fprintf(fp," dTheta: %f",msr->param.dTheta);
	fprintf(fp,"\n# dAbsPartial: %g",msr->param.dAbsPartial);
	fprintf(fp," dRealPartial: %g",msr->param.dRelPartial);
	fprintf(fp," dAbsTotal: %g",msr->param.dAbsTotal);
	fprintf(fp," dRelTotal: %g",msr->param.dRelTotal);
	fprintf(fp,"\n# dPeriod: %g",msr->param.dPeriod);
	fprintf(fp," dHubble0: %g",msr->param.dHubble0);
	fprintf(fp," dOmega0: %g",msr->param.dOmega0);
	fprintf(fp,"\n# achInFile: %s",msr->param.achInFile);
	fprintf(fp,"\n# achOutName: %s",msr->param.achOutName); 
	fprintf(fp,"\n# achDataSubPath: %s",msr->param.achDataSubPath);
	if (msr->param.bComove) {
		fprintf(fp,"\n# RedOut:");
		if (msr->nOuts == 0) fprintf(fp," none");
		for (i=0;i<msr->nOuts;i++) {
			if (i%5 == 0) fprintf(fp,"\n#   ");
			z = 1.0/msrTime2Exp(msr,msr->pdOutTime[i])-1.0;
			fprintf(fp," %f",z);
			}
		fprintf(fp,"\n");
		fflush(fp);
		}
	else {
		fprintf(fp,"\n# TimeOut:");
		if (msr->nOuts == 0) fprintf(fp," none");
		for (i=0;i<msr->nOuts;i++) {
			if (i%5 == 0) fprintf(fp,"\n#   ");
			fprintf(fp," %f",msr->pdOutTime[i]);
			}
		fprintf(fp,"\n");
		fflush(fp);
		}
	}


void msrFinish(MSR msr)
{
	int id;

	for (id=1;id<msr->nThreads;++id) {
		if (msr->param.bVerbose) printf("Stopping thread %d\n",id);		
		mdlReqService(msr->mdl,id,SRV_STOP,NULL,0);
		mdlGetReply(msr->mdl,id,NULL,NULL);
		}
	pstFinish(msr->pst);
	/*
	 ** finish with parameter stuff, dealocate and exit.
	 */
	prmFinish(msr->prm);
	free(msr->pMap);
	free(msr);
	}


/*
 ** Link with code from file eccanom.c and hypanom.c.
 */
double dEccAnom(double,double);
double dHypAnom(double,double);

double msrTime2Exp(MSR msr,double dTime)
{
	double dOmega0 = msr->param.dOmega0;
	double dHubble0 = msr->param.dHubble0;
	double a0,A,B,eta;

	if (!msr->param.bComove) return(1.0);
	if (dOmega0 == 1.0) {
		assert(dHubble0 > 0.0);
		if (dTime == 0.0) return(0.0);
		return(pow(3.0*dHubble0*dTime/2.0,2.0/3.0));
		}
	else if (dOmega0 > 1.0) {
		assert(dHubble0 >= 0.0);
		if (dHubble0 == 0.0) {
			B = 1.0/sqrt(dOmega0);
			eta = dEccAnom(dTime/B,1.0);
			return(1.0-cos(eta));
			}
		if (dTime == 0.0) return(0.0);
		a0 = 1.0/dHubble0/sqrt(dOmega0-1.0);
		A = 0.5*dOmega0/(dOmega0-1.0);
		B = A*a0;
		eta = dEccAnom(dTime/B,1.0);
		return(A*(1.0-cos(eta)));
		}
	else if (dOmega0 > 0.0) {
		assert(dHubble0 > 0.0);
		if (dTime == 0.0) return(0.0);
		a0 = 1.0/dHubble0/sqrt(1.0-dOmega0);
		A = 0.5*dOmega0/(1.0-dOmega0);
		B = A*a0;
		eta = dHypAnom(dTime/B,1.0);
		return(A*(cosh(eta)-1.0));
		}
	else if (dOmega0 == 0.0) {
		assert(dHubble0 > 0.0);
		if (dTime == 0.0) return(0.0);
		return(dTime*dHubble0);
		}
	else {
		/*
		 ** Bad value.
		 */
		assert(0);
		return(0.0);
		}
	}


double msrExp2Time(MSR msr,double dExp)
{
	double dOmega0 = msr->param.dOmega0;
	double dHubble0 = msr->param.dHubble0;
	double a0,A,B,eta;

	if (!msr->param.bComove) {
		/*
		 ** Invalid call!
		 */
		assert(0);
		}
	if (dOmega0 == 1.0) {
		assert(dHubble0 > 0.0);
		if (dExp == 0.0) return(0.0);
		return(2.0/(3.0*dHubble0)*pow(dExp,1.5));
		}
	else if (dOmega0 > 1.0) {
		assert(dHubble0 >= 0.0);
		if (dHubble0 == 0.0) {
			B = 1.0/sqrt(dOmega0);
			eta = acos(1.0-dExp); 
			return(B*(eta-sin(eta)));
			}
		if (dExp == 0.0) return(0.0);
		a0 = 1.0/dHubble0/sqrt(dOmega0-1.0);
		A = 0.5*dOmega0/(dOmega0-1.0);
		B = A*a0;
		eta = acos(1.0-dExp/A);
		return(B*(eta-sin(eta)));
		}
	else if (dOmega0 > 0.0) {
		assert(dHubble0 > 0.0);
		if (dExp == 0.0) return(0.0);
		a0 = 1.0/dHubble0/sqrt(1.0-dOmega0);
		A = 0.5*dOmega0/(1.0-dOmega0);
		B = A*a0;
		eta = acosh(dExp/A+1.0);
		return(B*(sinh(eta)-eta));
		}
	else if (dOmega0 == 0.0) {
		assert(dHubble0 > 0.0);
		if (dExp == 0.0) return(0.0);
		return(dExp/dHubble0);
		}
	else {
		/*
		 ** Bad value.
		 */
		assert(0);
		return(0.0);
		}	
	}


double msrTime2Hub(MSR msr,double dTime)
{
	double a = msrTime2Exp(msr,dTime);
	double z;

	assert(a > 0.0);
	z = 1.0/a-1.0;
	return(msr->param.dHubble0*(1.0+z)*sqrt(1.0+msr->param.dOmega0*z));
	}


/*
 ** This function integrates the time dependence of the "drift"-Hamiltonian.
 */
double msrComoveDriftFac(MSR msr,double dTime,double dDelta)
{
	double dOmega0 = msr->param.dOmega0;
	double dHubble0 = msr->param.dHubble0;
	double a0,A,B,a1,a2,eta1,eta2;

	if (!msr->param.bComove) return(dDelta);
	else {
		a1 = msrTime2Exp(msr,dTime);
		a2 = msrTime2Exp(msr,dTime+dDelta);
		if (dOmega0 == 1.0) {
			return((2.0/dHubble0)*(1.0/sqrt(a1) - 1.0/sqrt(a2)));
			}
		else if (dOmega0 > 1.0) {
			assert(dHubble0 >= 0.0);
			if (dHubble0 == 0.0) {
				A = 1.0;
				B = 1.0/sqrt(dOmega0);
				}
			else {
				a0 = 1.0/dHubble0/sqrt(dOmega0-1.0);
				A = 0.5*dOmega0/(dOmega0-1.0);
				B = A*a0;
				}
			eta1 = acos(1.0-a1/A);
			eta2 = acos(1.0-a2/A);
			return(B/A/A*(1.0/tan(0.5*eta1) - 1.0/tan(0.5*eta2)));
			}
		else if (dOmega0 > 0.0) {
			assert(dHubble0 > 0.0);
			a0 = 1.0/dHubble0/sqrt(1.0-dOmega0);
			A = 0.5*dOmega0/(1.0-dOmega0);
			B = A*a0;
			eta1 = acosh(a1/A+1.0);
			eta2 = acosh(a2/A+1.0);
			return(B/A/A*(1.0/tanh(0.5*eta1) - 1.0/tanh(0.5*eta2)));
			}
		else if (dOmega0 == 0.0) {
			/*
			 ** YOU figure this one out!
			 */
			assert(0);
			return(0.0);
			}
		else {
			/*
			 ** Bad value?
			 */
			assert(0);
			return(0.0);
			}
		}
	}


/*
 ** This function integrates the time dependence of the "kick"-Hamiltonian.
 */
double msrComoveKickFac(MSR msr,double dTime,double dDelta)
{
	double dOmega0 = msr->param.dOmega0;
	double dHubble0 = msr->param.dHubble0;
	double a0,A,B,a1,a2,eta1,eta2;

	if (!msr->param.bComove) return(dDelta);
	else {
		a1 = msrTime2Exp(msr,dTime);
		a2 = msrTime2Exp(msr,dTime+dDelta);
		if (dOmega0 == 1.0) {
			return((2.0/dHubble0)*(sqrt(a2) - sqrt(a1)));
			}
		else if (dOmega0 > 1.0) {
			assert(dHubble0 >= 0.0);
			if (dHubble0 == 0.0) {
				A = 1.0;
				B = 1.0/sqrt(dOmega0);
				}
			else {
				a0 = 1.0/dHubble0/sqrt(dOmega0-1.0);
				A = 0.5*dOmega0/(dOmega0-1.0);
				B = A*a0;
				}
			eta1 = acos(1.0-a1/A);
			eta2 = acos(1.0-a2/A);
			return(B/A*(eta2 - eta1));
			}
		else if (dOmega0 > 0.0) {
			assert(dHubble0 > 0.0);
			a0 = 1.0/dHubble0/sqrt(1.0-dOmega0);
			A = 0.5*dOmega0/(1.0-dOmega0);
			B = A*a0;
			eta1 = acosh(a1/A+1.0);
			eta2 = acosh(a2/A+1.0);
			return(B/A*(eta2 - eta1));
			}
		else if (dOmega0 == 0.0) {
			/*
			 ** YOU figure this one out!
			 */
			assert(0);
			return(0.0);
			}
		else {
			/*
			 ** Bad value?
			 */
			assert(0);
			return(0.0);
			}
		}
	}

void msrOneNodeReadTipsy(MSR msr, struct inReadTipsy *in)
{
    int i,id;
    int *nParts;		/* number of particles for each processor */
    int nStart;
    PST pst0;
    LCL *plcl;
    char achInFile[PST_FILENAME_SIZE];
    int nid;
    int inswap;
	
    nParts = malloc(msr->nThreads*sizeof(*nParts));
    for (id=0;id<msr->nThreads;++id) {
		nParts[id] = -1;
		}
	
    pstOneNodeReadInit(msr->pst, in, sizeof(*in), nParts, &nid);
    assert(nid/sizeof(*nParts) == msr->nThreads);
    for (id=0;id<msr->nThreads;++id) {
		assert(nParts[id] > 0);
		}
	
    pst0 = msr->pst;
    while(pst0->nLeaves > 1)
		pst0 = pst0->pstLower;
    plcl = pst0->plcl;
    /*
     ** Add the local Data Path to the provided filename.
     */
    achInFile[0] = 0;
    if (plcl->pszDataPath) {
	    strcat(achInFile,plcl->pszDataPath);
	    strcat(achInFile,"/");
	    }
    strcat(achInFile,in->achInFile);
	
    nStart = nParts[0];
	assert(msr->pMap[0] == 0);
    for (i=1;i<msr->nThreads;++i) {
		id = msr->pMap[i];
		/* 
		 * Read particles into the local storage.
		 */
		assert(plcl->pkd->nStore >= nParts[id]);
		pkdReadTipsy(plcl->pkd,achInFile,nStart,nParts[id],in->dvFac);
		nStart += nParts[id];
		/* 
		 * Now shove them over to the remote processor.
		 */
		inswap = 0;
		mdlReqService(pst0->mdl,id,PST_SWAPALL,&inswap,sizeof(inswap));
		pkdSwapAll(plcl->pkd, id);
		mdlGetReply(pst0->mdl,id,NULL,NULL);
    	}
    assert(nStart == msr->N);
    /* 
     * Now read our own particles.
     */
    pkdReadTipsy(plcl->pkd,achInFile,0,nParts[0],in->dvFac);
    }

double msrReadTipsy(MSR msr)
{
	FILE *fp;
	struct dump h;
	struct inReadTipsy in;
	char achInFile[PST_FILENAME_SIZE];
	int j;
	LCL *plcl = msr->pst->plcl;
	double dTime,aTo,tTo,z;
	
	if (msr->param.achInFile[0]) {
		/*
		 ** Add Data Subpath for local and non-local names.
		 */
		achInFile[0] = 0;
		strcat(achInFile,msr->param.achDataSubPath);
		strcat(achInFile,"/");
		strcat(achInFile,msr->param.achInFile);
		strcpy(in.achInFile,achInFile);
		/*
		 ** Add local Data Path.
		 */
		achInFile[0] = 0;
		if (plcl->pszDataPath) {
			strcat(achInFile,plcl->pszDataPath);
			strcat(achInFile,"/");
			}
		strcat(achInFile,in.achInFile);

		fp = fopen(achInFile,"r");
		if (!fp) {
			printf("Could not open InFile:%s\n",achInFile);
			_msrExit(msr);
			}
		}
	else {
		printf("No input file specified\n");
		_msrExit(msr);
		}
	/*
	 ** Assume tipsy format for now, and dark matter only.
	 */
	fread(&h,sizeof(struct dump),1,fp);
	fclose(fp);
	msr->N = h.nbodies;
	assert(msr->N == h.ndark);
	if (msr->param.bComove) {
		if(msr->param.dHubble0 == 0.0) {
			printf("No hubble constant specified\n");
			_msrExit(msr);
			}
		dTime = msrExp2Time(msr,h.time);
		z = 1.0/h.time - 1.0;
		printf("Input file, Time:%g Redshift:%g Expansion factor:%g\n",
			   dTime,z,h.time);
		if (prmSpecified(msr->prm,"dRedTo")) {
			if (!prmArgSpecified(msr->prm,"nSteps") &&
				prmArgSpecified(msr->prm,"dDelta")) {
				aTo = 1.0/(msr->param.dRedTo + 1.0);
				tTo = msrExp2Time(msr,aTo);
				printf("Simulation to Time:%g Redshift:%g Expansion factor:%g\n",
					   tTo,1.0/aTo-1.0,aTo);
				if (tTo < dTime) {
					printf("Badly specified final redshift, check -zto parameter.\n");
					_msrExit(msr);
					}
				msr->param.nSteps = (int)ceil((tTo-dTime)/msr->param.dDelta);
				}
			else if (!prmArgSpecified(msr->prm,"dDelta") &&
					 prmArgSpecified(msr->prm,"nSteps")) {
				aTo = 1.0/(msr->param.dRedTo + 1.0);
				tTo = msrExp2Time(msr,aTo);
				printf("Simulation to Time:%g Redshift:%g Expansion factor:%g\n",
					   tTo,1.0/aTo-1.0,aTo);
				if (tTo < dTime) {
					printf("Badly specified final redshift, check -zto parameter.\n");	
					_msrExit(msr);
					}
				if(msr->param.nSteps != 0)
				    msr->param.dDelta = (tTo-dTime)/msr->param.nSteps;
				else
				    msr->param.dDelta = 0.0;
				}
			else if (!prmSpecified(msr->prm,"nSteps") &&
				prmFileSpecified(msr->prm,"dDelta")) {
				aTo = 1.0/(msr->param.dRedTo + 1.0);
				tTo = msrExp2Time(msr,aTo);
				printf("Simulation to Time:%g Redshift:%g Expansion factor:%g\n",
					   tTo,1.0/aTo-1.0,aTo);
				if (tTo < dTime) {
					printf("Badly specified final redshift, check -zto parameter.\n");
					_msrExit(msr);
					}
				msr->param.nSteps = (int)ceil((tTo-dTime)/msr->param.dDelta);
				}
			else if (!prmSpecified(msr->prm,"dDelta") &&
					 prmFileSpecified(msr->prm,"nSteps")) {
				aTo = 1.0/(msr->param.dRedTo + 1.0);
				tTo = msrExp2Time(msr,aTo);
				printf("Simulation to Time:%g Redshift:%g Expansion factor:%g\n",
					   tTo,1.0/aTo-1.0,aTo);
				if (tTo < dTime) {
					printf("Badly specified final redshift, check -zto parameter.\n");	
					_msrExit(msr);
					}
				if(msr->param.nSteps != 0)
				    msr->param.dDelta = (tTo-dTime)/msr->param.nSteps;
				else
				    msr->param.dDelta = 0.0;
				}
			}
		else {
			tTo = dTime + msr->param.nSteps*msr->param.dDelta;
			aTo = msrTime2Exp(msr,tTo);
			printf("Simulation to Time:%g Redshift:%g Expansion factor:%g\n",
				   tTo,1.0/aTo-1.0,aTo);
			}
		printf("Reading file...\nN:%d\n",msr->N);
		if (msr->param.bCannonical) {
			in.dvFac = h.time*h.time;
			}
		else {
			in.dvFac = 1.0;
			}
		}
	else {
		dTime = h.time;
		printf("Input file, Time:%g\n",dTime);
		tTo = dTime + msr->param.nSteps*msr->param.dDelta;
		printf("Simulation to Time:%g\n",tTo);
		printf("Reading file...\nN:%d Time:%g\n",msr->N,dTime);
		in.dvFac = 1.0;
		}
	in.nStart = 0;
	in.nEnd = msr->N - 1;
	in.iOrder = msr->param.iOrder;
	/*
	 ** Since pstReadTipsy causes the allocation of the local particle
	 ** store, we need to tell it the percentage of extra storage it
	 ** should allocate for load balancing differences in the number of
	 ** particles.
	 */
	in.fExtraStore = msr->param.dExtraStore;
	/*
	 ** Provide the period.
	 */
	for (j=0;j<3;++j) {
		in.fPeriod[j] = msr->param.dPeriod;
		}
	if(msr->param.bParaRead)
	    pstReadTipsy(msr->pst,&in,sizeof(in),NULL,NULL);
	else
	    msrOneNodeReadTipsy(msr, &in);
	if (msr->param.bVerbose) puts("Input file has been successfully read.");
	/*
	 ** Now read in the output points, passing the initial time.
	 ** We do this only if nSteps is not equal to zero.
	 */
	if (msrSteps(msr) > 0) msrReadOuts(msr,dTime);
	/*
	 ** Set up the output counter.
	 */
	for (msr->iOut=0;msr->iOut<msr->nOuts;++msr->iOut) {
		if (dTime < msr->pdOutTime[msr->iOut]) break;
		}
	return(dTime);
	}


int xdrWriteHeader(XDR *pxdrs,struct dump *ph)
{
	int pad = 0;
	
	if (!xdr_double(pxdrs,&ph->time)) return 0;
	if (!xdr_int(pxdrs,&ph->nbodies)) return 0;
	if (!xdr_int(pxdrs,&ph->ndim)) return 0;
	if (!xdr_int(pxdrs,&ph->nsph)) return 0;
	if (!xdr_int(pxdrs,&ph->ndark)) return 0;
	if (!xdr_int(pxdrs,&ph->nstar)) return 0;
	if (!xdr_int(pxdrs,&pad)) return 0;
	return 1;
	}

void msrOneNodeWriteTipsy(MSR msr, struct inWriteTipsy *in)
{
    int i,id;
    int nStart;
    PST pst0;
    LCL *plcl;
    char achOutFile[PST_FILENAME_SIZE];
    int inswap;

    pst0 = msr->pst;
    while(pst0->nLeaves > 1)
		pst0 = pst0->pstLower;
    plcl = pst0->plcl;
    /*
     ** Add the local Data Path to the provided filename.
     */
    achOutFile[0] = 0;
    if (plcl->pszDataPath) {
	    strcat(achOutFile,plcl->pszDataPath);
	    strcat(achOutFile,"/");
	    }
    strcat(achOutFile,in->achOutFile);

    /* 
     * First write our own particles.
     */
    pkdWriteTipsy(plcl->pkd,achOutFile,0,plcl->pkd->nLocal - 1,in->bStandard,
				  in->dvFac); 
    nStart = plcl->pkd->nLocal;
	assert(msr->pMap[0] == 0);
    for (i=1;i<msr->nThreads;++i) {
		id = msr->pMap[i];
		/* 
		 * Swap particles with the remote processor.
		 */
		inswap = 0;
		mdlReqService(pst0->mdl,id,PST_SWAPALL,&inswap,sizeof(inswap));
		pkdSwapAll(plcl->pkd, id);
		mdlGetReply(pst0->mdl,id,NULL,NULL);
		/* 
		 * Write the swapped particles.
		 */
		pkdWriteTipsy(plcl->pkd,achOutFile,nStart,
					  nStart + plcl->pkd->nLocal - 1,
					  in->bStandard, in->dvFac); 
		nStart += plcl->pkd->nLocal;
		/* 
		 * Swap them back again.
		 */
		inswap = 0;
		mdlReqService(pst0->mdl,id,PST_SWAPALL,&inswap,sizeof(inswap));
		pkdSwapAll(plcl->pkd, id);
		mdlGetReply(pst0->mdl,id,NULL,NULL);
    	}
    assert(nStart == msr->N);
    }

void msrWriteTipsy(MSR msr,char *pszFileName,double dTime)
{
	FILE *fp;
	struct dump h;
	struct inWriteTipsy in;
	char achOutFile[PST_FILENAME_SIZE];
	LCL *plcl = msr->pst->plcl;
	
	/*
	 ** Add Data Subpath for local and non-local names.
	 */
	achOutFile[0] = 0;
	strcat(achOutFile,msr->param.achDataSubPath);
	strcat(achOutFile,"/");
	strcat(achOutFile,pszFileName);
	strcpy(in.achOutFile,achOutFile);
	/*
	 ** Add local Data Path.
	 */
	achOutFile[0] = 0;
	if (plcl->pszDataPath) {
		strcat(achOutFile,plcl->pszDataPath);
		strcat(achOutFile,"/");
		}
	strcat(achOutFile,in.achOutFile);
	
	fp = fopen(achOutFile,"w");
	if (!fp) {
		printf("Could not open OutFile:%s\n",achOutFile);
		_msrExit(msr);
		}
	in.bStandard = msr->param.bStandard;
	/*
	 ** Assume tipsy format for now, and dark matter only.
	 */
	h.nbodies = msr->N;
	h.ndark = msr->N;
	h.nsph = 0;
	h.nstar = 0;
	if (msr->param.bComove) {
		h.time = msrTime2Exp(msr,dTime);
		if (msr->param.bCannonical) {
			in.dvFac = 1.0/(h.time*h.time);
			}
		else {
			in.dvFac = 1.0;
			}
		}
	else {
		h.time = dTime;
		in.dvFac = 1.0;
		}
	h.ndim = 3;
	if (msr->param.bVerbose) {
		if (msr->param.bComove) {
			printf("Writing file...\nTime:%g Redshift:%g\n",
				   dTime,(1.0/h.time - 1.0));
			}
		else {
			printf("Writing file...\nTime:%g\n",dTime);
			}
		}
	if (in.bStandard) {
		XDR xdrs;

		xdrstdio_create(&xdrs,fp,XDR_ENCODE);
		xdrWriteHeader(&xdrs,&h);
		xdr_destroy(&xdrs);
		}
	else {
		fwrite(&h,sizeof(struct dump),1,fp);
		}
	fclose(fp);
	if(msr->param.bParaWrite)
	    pstWriteTipsy(msr->pst,&in,sizeof(in),NULL,NULL);
	else
	    msrOneNodeWriteTipsy(msr, &in);
	if (msr->param.bVerbose) {
		puts("Output file has been successfully written.");
		}
	}


void msrSetSoft(MSR msr,double dSoft)
{
	struct inSetSoft in;
  
	if (msr->param.bVerbose) printf("Set Softening...\n");
	in.dSoft = dSoft;
	pstSetSoft(msr->pst,&in,sizeof(in),NULL,NULL);
	}


void msrBuildTree(MSR msr,int bActiveOnly,double dMass,int bSmooth)
{
	struct inBuildTree in;
	struct outBuildTree out;
	struct inColCells inc;
	struct ioCalcRoot root;
	KDN *pkdn;
	int sec,dsec;
	int iDum,nCell;
	int i;

	if (msr->param.bVerbose) printf("Domain Decomposition...\n");
	sec = time(0);
	pstDomainDecomp(msr->pst,NULL,0,NULL,NULL);
	msrMassCheck(msr,dMass,"After pstDomainDecomp in msrBuildTree");
	dsec = time(0) - sec;
	if (msr->param.bVerbose) {
		printf("Domain Decomposition complete, Wallclock: %d secs\n\n",dsec);
		}
	if (msr->param.bVerbose) printf("Building local trees...\n");
	/*
	 ** First make sure the particles are in Active/Inactive order.
	 */
	pstActiveOrder(msr->pst,NULL,0,NULL,NULL);
	in.nBucket = msr->param.nBucket;
	in.iOpenType = msr->iOpenType;
	in.iOrder = (msr->param.iOrder >= msr->param.iEwOrder)?
		msr->param.iOrder:msr->param.iEwOrder;
	in.dCrit = msr->dCrit;
	if (bSmooth) {
		in.bBinary = 0;
		msr->iTreeType = MSR_TREE_DENSITY;
		}
	else {
		in.bBinary = msr->param.bBinary;
		if (msr->param.bBinary) {
			msr->iTreeType = MSR_TREE_SPATIAL;
			}
		else {
			msr->iTreeType = MSR_TREE_DENSITY;
			}
		}
	in.bActiveOnly = bActiveOnly;
	sec = time(0);
	pstBuildTree(msr->pst,&in,sizeof(in),&out,&iDum);
	msrMassCheck(msr,dMass,"After pstBuildTree in msrBuildTree");
	dsec = time(0) - sec;
	if (msr->param.bVerbose) {
		printf("Tree built, Wallclock: %d secs\n\n",dsec);
		}
	nCell = 1<<(1+(int)ceil(log((double)msr->nThreads)/log(2.0)));
	pkdn = malloc(nCell*sizeof(KDN));
	assert(pkdn != NULL);
	inc.iCell = ROOT;
	inc.nCell = nCell;
	pstColCells(msr->pst,&inc,sizeof(inc),pkdn,NULL);
	msrMassCheck(msr,dMass,"After pstColCells in msrBuildTree");
#if (0)
	for (i=1;i<nCell;++i) {
		struct pkdCalcCellStruct *m;

		printf("\nLTTO:%d\n",i);
		printf("    iDim:%1d fSplit:%g pLower:%d pUpper:%d\n",
			   pkdn[i].iDim,pkdn[i].fSplit,pkdn[i].pLower,pkdn[i].pUpper);
		printf("    bnd:(%g,%g) (%g,%g) (%g,%g)\n",
			   pkdn[i].bnd.fMin[0],pkdn[i].bnd.fMax[0],
			   pkdn[i].bnd.fMin[1],pkdn[i].bnd.fMax[1],
			   pkdn[i].bnd.fMin[2],pkdn[i].bnd.fMax[2]);
		printf("    fMass:%g fSoft:%g\n",pkdn[i].fMass,pkdn[i].fSoft);
		printf("    rcm:%g %g %g fOpen2:%g\n",pkdn[i].r[0],pkdn[i].r[1],
			   pkdn[i].r[2],pkdn[i].fOpen2);
		m = &pkdn[i].mom;
		printf("    xx:%g yy:%g zz:%g xy:%g xz:%g yz:%g\n",
			   m->Qxx,m->Qyy,m->Qzz,m->Qxy,m->Qxz,m->Qyz);
		printf("    xxx:%g xyy:%g xxy:%g yyy:%g xxz:%g yyz:%g xyz:%g\n",
			   m->Oxxx,m->Oxyy,m->Oxxy,m->Oyyy,m->Oxxz,m->Oyyz,m->Oxyz);
		printf("    xxxx:%g xyyy:%g xxxy:%g yyyy:%g xxxz:%g\n",
			   m->Hxxxx,m->Hxyyy,m->Hxxxy,m->Hyyyy,m->Hxxxz);
		printf("    yyyz:%g xxyy:%g xxyz:%g xyyz:%g\n",
			   m->Hyyyz,m->Hxxyy,m->Hxxyz,m->Hxyyz);
		}
#endif
	pstDistribCells(msr->pst,pkdn,nCell*sizeof(KDN),NULL,NULL);
	msrMassCheck(msr,dMass,"After pstDistribCells in msrBuildTree");
	free(pkdn);
	pstCalcRoot(msr->pst,NULL,0,&root,&iDum);
	msrMassCheck(msr,dMass,"After pstCalcRoot in msrBuildTree");
	pstDistribRoot(msr->pst,&root,sizeof(struct ioCalcRoot),NULL,NULL);
	msrMassCheck(msr,dMass,"After pstDistribRoot in msrBuildTree");
	}


void msrDomainColor(MSR msr)
{
	pstDomainColor(msr->pst,NULL,0,NULL,NULL);
	}


void msrReorder(MSR msr)
{
	int sec,dsec;

	if (msr->param.bVerbose) printf("Ordering...\n");
	sec = time(0);
	pstDomainOrder(msr->pst,NULL,0,NULL,NULL);
	pstLocalOrder(msr->pst,NULL,0,NULL,NULL);
	dsec = time(0) - sec;
	if (msr->param.bVerbose) {
		printf("Order established, Wallclock: %d secs\n\n",dsec);
		}
	/*
	 ** Mark tree as void.
	 */
	msr->iTreeType = MSR_TREE_NONE;
	}


void msrOutArray(MSR msr,char *pszFile,int iType)
{
	struct inOutArray in;
	char achOutFile[PST_FILENAME_SIZE];
	LCL *plcl;
	PST pst0;
	FILE *fp;
	int id,i;
	int inswap;

	pst0 = msr->pst;
	while(pst0->nLeaves > 1)
	    pst0 = pst0->pstLower;
	plcl = pst0->plcl;
	if (pszFile) {
		/*
		 ** Add Data Subpath for local and non-local names.
		 */
		achOutFile[0] = 0;
		strcat(achOutFile,msr->param.achDataSubPath);
		strcat(achOutFile,"/");
		strcat(achOutFile,pszFile);
		strcpy(in.achOutFile,achOutFile);
		/*
		 ** Add local Data Path.
		 */
		achOutFile[0] = 0;
		if (plcl->pszDataPath) {
			strcat(achOutFile,plcl->pszDataPath);
			strcat(achOutFile,"/");
			}
		strcat(achOutFile,in.achOutFile);

		fp = fopen(achOutFile,"w");
		if (!fp) {
			printf("Could not open Array Output File:%s\n",achOutFile);
			_msrExit(msr);
			}
		}
	else {
		printf("No Array Output File specified\n");
		_msrExit(msr);
		}
	/*
	 ** Write the Header information and close the file again.
	 */
	fprintf(fp,"%d\n",msr->N);
	fclose(fp);
	/* 
	 * First write our own particles.
	 */
	assert(msr->pMap[0] == 0);
	pkdOutArray(plcl->pkd,achOutFile,iType); 
	for (i=1;i<msr->nThreads;++i) {
		id = msr->pMap[i];
	    /* 
	     * Swap particles with the remote processor.
	     */
	    inswap = 0;
	    mdlReqService(pst0->mdl,id,PST_SWAPALL,&inswap,sizeof(inswap));
	    pkdSwapAll(plcl->pkd, id);
	    mdlGetReply(pst0->mdl,id,NULL,NULL);
	    /* 
	     * Write the swapped particles.
	     */
	    pkdOutArray(plcl->pkd,achOutFile,iType); 
	    /* 
	     * Swap them back again.
	     */
	    inswap = 0;
	    mdlReqService(pst0->mdl,id,PST_SWAPALL,&inswap,sizeof(inswap));
	    pkdSwapAll(plcl->pkd, id);
	    mdlGetReply(pst0->mdl,id,NULL,NULL);
	    }
	}


void msrOutVector(MSR msr,char *pszFile,int iType)
{
	struct inOutVector in;
	char achOutFile[PST_FILENAME_SIZE];
	LCL *plcl;
	PST pst0;
	FILE *fp;
	int id,i;
	int inswap;
	int iDim;

	pst0 = msr->pst;
	while(pst0->nLeaves > 1)
	    pst0 = pst0->pstLower;
	plcl = pst0->plcl;
	if (pszFile) {
		/*
		 ** Add Data Subpath for local and non-local names.
		 */
		achOutFile[0] = 0;
		strcat(achOutFile,msr->param.achDataSubPath);
		strcat(achOutFile,"/");
		strcat(achOutFile,pszFile);
		strcpy(in.achOutFile,achOutFile);
		/*
		 ** Add local Data Path.
		 */
		achOutFile[0] = 0;
		if (plcl->pszDataPath) {
			strcat(achOutFile,plcl->pszDataPath);
			strcat(achOutFile,"/");
			}
		strcat(achOutFile,in.achOutFile);

		fp = fopen(achOutFile,"w");
		if (!fp) {
			printf("Could not open Vector Output File:%s\n",achOutFile);
			_msrExit(msr);
			}
		}
	else {
		printf("No Vector Output File specified\n");
		_msrExit(msr);
		}
	/*
	 ** Write the Header information and close the file again.
	 */
	fprintf(fp,"%d\n",msr->N);
	fclose(fp);
	/* 
	 * First write our own particles.
	 */
	assert(msr->pMap[0] == 0);
	for (iDim=0;iDim<3;++iDim) {
		pkdOutVector(plcl->pkd,achOutFile,iDim,iType); 
		for (i=1;i<msr->nThreads;++i) {
			id = msr->pMap[i];
			/* 
			 * Swap particles with the remote processor.
			 */
			inswap = 0;
			mdlReqService(pst0->mdl,id,PST_SWAPALL,&inswap,sizeof(inswap));
			pkdSwapAll(plcl->pkd, id);
			mdlGetReply(pst0->mdl,id,NULL,NULL);
			/* 
			 * Write the swapped particles.
			 */
			pkdOutVector(plcl->pkd,achOutFile,iDim,iType); 
			/* 
			 * Swap them back again.
			 */
			inswap = 0;
			mdlReqService(pst0->mdl,id,PST_SWAPALL,&inswap,sizeof(inswap));
			pkdSwapAll(plcl->pkd,id);
			mdlGetReply(pst0->mdl,id,NULL,NULL);
			}
		}
	}


void msrSmooth(MSR msr,int iSmoothType,int bSymmetric)
{
	struct inSmooth in;
	int sec,dsec;

	/*
	 ** Make sure that the type of tree is a density binary tree!
	 */
	assert(msr->iTreeType == MSR_TREE_DENSITY);
	in.nSmooth = msr->param.nSmooth;
	in.bPeriodic = msr->param.bPeriodic;
	in.bSymmetric = bSymmetric;
	in.iSmoothType = iSmoothType;
	if (msr->param.bVerbose) printf("Smoothing...\n");
	sec = time(0);
	pstSmooth(msr->pst,&in,sizeof(in),NULL,NULL);
	dsec = time(0) - sec;
	printf("Smooth Calculated, Wallclock:%d secs\n\n",dsec);
	}


void msrReSmooth(MSR msr,int iSmoothType,int bSymmetric)
{
	struct inReSmooth in;
	int sec,dsec;

	/*
	 ** Make sure that the type of tree is a density binary tree!
	 */
	assert(msr->iTreeType == MSR_TREE_DENSITY);
	in.nSmooth = msr->param.nSmooth;
	in.bPeriodic = msr->param.bPeriodic;
	in.bSymmetric = bSymmetric;
	in.iSmoothType = iSmoothType;
	if (msr->param.bVerbose) printf("ReSmoothing...\n");
	sec = time(0);
	pstReSmooth(msr->pst,&in,sizeof(in),NULL,NULL);
	dsec = time(0) - sec;
	printf("ReSmooth Calculated, Wallclock:%d secs\n\n",dsec);
	}


void msrGravity(MSR msr,double dStep,
				int *piSec,double *pdWMax,double *pdIMax,double *pdEMax,
				int *pnActive)
{
	struct inGravity in;
	struct outGravity out;
	int iDum;
	int sec,dsec;
	double dPartAvg,dCellAvg;
	double dMFlops;
	double dWAvg,dWMax,dWMin;
	double dIAvg,dIMax,dIMin;
	double dEAvg,dEMax,dEMin;
	double iP;

	assert(msr->iTreeType == MSR_TREE_SPATIAL || 
		   msr->iTreeType == MSR_TREE_DENSITY);
	printf("Calculating Gravity, Step:%.1f\n",dStep);
    in.nReps = msr->param.nReplicas;
    in.bPeriodic = msr->param.bPeriodic;
	in.iOrder = msr->param.iOrder;
	in.iEwOrder = msr->param.iEwOrder;
    in.dEwCut = msr->param.dEwCut;
    in.dEwhCut = msr->param.dEwhCut;
	sec = time(0);
	pstGravity(msr->pst,&in,sizeof(in),&out,&iDum);
	dsec = time(0) - sec;
	dMFlops = out.dFlop/dsec*1e-6;
	printf("Gravity Calculated, Wallclock:%d secs, MFlops:%.1f, Flop:%.3g\n",
		   dsec,dMFlops,out.dFlop);
	*piSec = dsec;
	dPartAvg = out.dPartSum/out.nActive;
	dCellAvg = out.dCellSum/out.nActive;
	*pnActive = out.nActive;
	iP = 1.0/msr->nThreads;
	dWAvg = out.dWSum*iP;
	dIAvg = out.dISum*iP;
	dEAvg = out.dESum*iP;
	dWMax = out.dWMax;
	*pdWMax = dWMax;
	dIMax = out.dIMax;
	*pdIMax = dIMax;
	dEMax = out.dEMax;
	*pdEMax = dEMax;
	dWMin = out.dWMin;
	dIMin = out.dIMin;
	dEMin = out.dEMin;
	printf("dPartAvg:%f dCellAvg:%f\n",dPartAvg,dCellAvg);
	printf("Walk CPU     Avg:%10f Max:%10f Min:%10f\n",dWAvg,dWMax,dWMin);
	printf("Interact CPU Avg:%10f Max:%10f Min:%10f\n",dIAvg,dIMax,dIMin);
	printf("Ewald CPU    Avg:%10f Max:%10f Min:%10f\n",dEAvg,dEMax,dEMin);	
	if (msr->nThreads > 1) {
		printf("Particle Cache Statistics (average per processor):\n");
		printf("    Accesses:    %10g\n",out.dpASum*iP);
		printf("    Miss Ratio:  %10g\n",out.dpMSum*iP);
		printf("    Min Ratio:   %10g\n",out.dpTSum*iP);
		printf("    Coll Ratio:  %10g\n",out.dpCSum*iP);
		printf("Cell Cache Statistics (average per processor):\n");
		printf("    Accesses:    %10g\n",out.dcASum*iP);
		printf("    Miss Ratio:  %10g\n",out.dcMSum*iP);
		printf("    Min Ratio:   %10g\n",out.dcTSum*iP);
		printf("    Coll Ratio:  %10g\n",out.dcCSum*iP);
		printf("\n");
		}
	}


void msrCalcE(MSR msr,int bFirst,double dTime,double *E,double *T,double *U)
{
	struct outCalcE out;
	double a;

	pstCalcE(msr->pst,NULL,0,&out,NULL);
	*T = out.T;
	*U = out.U;
	/*
	 ** Do the comoving coordinates stuff.
	 */
	a = msrTime2Exp(msr,dTime);
	*U *= a;
	if (!msr->param.bCannonical) *T *= pow(a,4.0);
	if (msr->param.bComove && !bFirst) {
		msr->dEcosmo += 0.5*(dTime - msr->dTimeOld)*
			(msrTime2Hub(msr,dTime)*(*U) + 
			 msrTime2Hub(msr,msr->dTimeOld)*msr->dUOld);
		}
	else {
		msr->dEcosmo = 0.0;
		}
	msr->dTimeOld = dTime;
	msr->dUOld = *U;
	*E = (*T) + (*U) - msr->dEcosmo;
	}


void msrDrift(MSR msr,double dTime,double dDelta)
{
	struct inDrift in;
	int j;

	if (msr->param.bCannonical) {
		in.dDelta = msrComoveDriftFac(msr,dTime,dDelta);
		}
	else {
		in.dDelta = dDelta;
		}
	for (j=0;j<3;++j) {
		in.fCenter[j] = msr->fCenter[j];
		}
	in.bPeriodic = msr->param.bPeriodic;
	pstDrift(msr->pst,&in,sizeof(in),NULL,NULL);
	/*
	 ** Once we move the particles the tree should no longer be considered 
	 ** valid.
	 */
	msr->iTreeType = MSR_TREE_NONE;
	}


void msrCoolVelocity(MSR msr,double dMass)
{
	struct inActiveCool ina;
	struct inCoolVelocity in;
	
	if (msr->param.nSuperCool > 0) {
		/*
		 ** First calculate densities for the supercool particles.
		 */
		ina.nSuperCool = msr->param.nSuperCool;
		pstActiveCool(msr->pst,&ina,sizeof(ina),NULL,NULL);
		/*
		 ** Build the tree, but with active only particles, just like 
		 ** sph see!
		 */
		msrBuildTree(msr,1,dMass,1);
		msrSmooth(msr,SMX_DENSITY,1);
#ifdef SUPERCOOL
		msrReSmooth(msr,SMX_MEANVEL,1);
#endif
		/*
		 ** Now cool them.
		 */
		in.nSuperCool = msr->param.nSuperCool;
		in.dCoolFac = msr->param.dCoolFac;
		in.dCoolDens = msr->param.dCoolDens;
		pstCoolVelocity(msr->pst,&in,sizeof(in),NULL,NULL);
		}
	}


void msrKick(MSR msr,double dTime,double dDelta)
{
	double H,a;
	struct inKick in;
	
	if (msr->param.bCannonical) {
		in.dvFacOne = 1.0; /* no hubble drag, man! */
		in.dvFacTwo = msrComoveKickFac(msr,dTime,dDelta);
		}
	else {
		/*
		 ** Careful! For non-cannonical we want H and a at the 
		 ** HALF-STEP! This is a bit messy but has to be special
		 ** cased in some way.
		 */
		dTime += dDelta/2.0;
		a = msrTime2Exp(msr,dTime);
		H = msrTime2Hub(msr,dTime);
		in.dvFacOne = (1.0 - H*dDelta)/(1.0 + H*dDelta);
		in.dvFacTwo = dDelta/pow(a,3.0)/(1.0 + H*dDelta);
		}
	pstKick(msr->pst,&in,sizeof(in),NULL,NULL);
	}


void msrOneNodeReadCheck(MSR msr, struct inReadCheck *in)
{
    int i,id;
    int *nParts;				/* number of particles for each processor */
    int nStart;
    PST pst0;
    LCL *plcl;
    char achInFile[PST_FILENAME_SIZE];
    int nid;
    int inswap;
    struct inReadTipsy tin;
    int j;

    nParts = malloc(msr->nThreads*sizeof(*nParts));
    for (id=0;id<msr->nThreads;++id) {
		nParts[id] = -1;
		}

    tin.nStart = in->nStart;
    tin.nEnd = in->nEnd;
    tin.iOrder = in->iOrder;
    tin.fExtraStore = in->fExtraStore;
    for(j = 0; j < 3; j++)
		tin.fPeriod[j] = in->fPeriod[j];

    pstOneNodeReadInit(msr->pst, &tin, sizeof(tin), nParts, &nid);
    assert(nid == msr->nThreads);
    for (id=0;id<msr->nThreads;++id) {
		assert(nParts[id] > 0);
		}

    pst0 = msr->pst;
    while(pst0->nLeaves > 1)
		pst0 = pst0->pstLower;
    plcl = pst0->plcl;
    /*
     ** Add the local Data Path to the provided filename.
     */
    achInFile[0] = 0;
    if (plcl->pszDataPath) {
	    strcat(achInFile,plcl->pszDataPath);
	    strcat(achInFile,"/");
	    }
    strcat(achInFile,in->achInFile);

    nStart = nParts[0];
	assert(msr->pMap[0] == 0);
    for (i=1;i<msr->nThreads;++i) {
		id = msr->pMap[i];
		/* 
		 * Read particles into the local storage.
		 */
		assert(plcl->pkd->nStore >= nParts[id]);
		pkdReadCheck(plcl->pkd,achInFile,in->iVersion, in->iOffset,
					 in->nStart, nParts[id]);
		nStart += nParts[id];
		/* 
		 * Now shove them over to the remote processor.
		 */
		inswap = 0;
		mdlReqService(pst0->mdl,id,PST_SWAPALL,&inswap,sizeof(inswap));
		pkdSwapAll(plcl->pkd, id);
		mdlGetReply(pst0->mdl,id,NULL,NULL);
    	}
    assert(nStart == msr->N);
    /* 
     * Now read our own particles.
     */
    pkdReadCheck(plcl->pkd,achInFile,in->iVersion, in->iOffset, 0, nParts[0]);
    }

double msrReadCheck(MSR msr,int *piStep)
{
	struct inReadCheck in;
	char achInFile[PST_FILENAME_SIZE];
	int i,j;
	LCL *plcl = msr->pst->plcl;
	double dTime;
	int iVersion,iNotCorrupt;
	FDL_CTX *fdl;

	/*
	 ** Add Data Subpath for local and non-local names.
	 */
	sprintf(achInFile,"%s/%s.chk",msr->param.achDataSubPath,
			msr->param.achOutName);
	strcpy(in.achInFile,achInFile);
	/*
	 ** Add local Data Path.
	 */
	if (plcl->pszDataPath) {
		sprintf(achInFile,"%s/%s",plcl->pszDataPath,in.achInFile);
		}
	fdl = FDL_open(achInFile);
	FDL_read(fdl,"version",&iVersion);
	if (msr->param.bVerbose) 
		printf("Reading Version-%d Checkpoint file.\n",iVersion);
	FDL_read(fdl,"not_corrupt_flag",&iNotCorrupt);
	if (iNotCorrupt != 1) {
	printf("Sorry the checkpoint file is corrupted.\n");
		_msrExit(msr);
		}
	FDL_read(fdl,"number_of_particles",&msr->N);
	FDL_read(fdl,"current_timestep",piStep);
	FDL_read(fdl,"current_time",&dTime);
	FDL_read(fdl,"current_ecosmo",&msr->dEcosmo);
	FDL_read(fdl,"old_time",&msr->dTimeOld);
	FDL_read(fdl,"old_potentiale",&msr->dUOld);
	if (!msr->bOpenSpec) {
		FDL_read(fdl,"opening_type",&msr->iOpenType);
		FDL_read(fdl,"opening_criterion",&msr->dCrit);
		}
	FDL_read(fdl,"number_of_outs",&msr->nOuts);
	if (msr->nOuts > msr->nMaxOuts) {
		msr->nMaxOuts = msr->nOuts;
		msr->pdOutTime = realloc(msr->pdOutTime,msr->nMaxOuts*sizeof(double));
		assert(msr->pdOutTime != NULL);
		}
	for (i=0;i<msr->nOuts;++i) {
		FDL_index(fdl,"out_time_index",i);
		FDL_read(fdl,"out_time",&msr->pdOutTime[i]);
		}
	/*
	 ** Read the old parameters.
	 */
	FDL_read(fdl,"bPeriodic",&msr->param.bPeriodic);
	FDL_read(fdl,"bComove",&msr->param.bComove);
	if (!prmSpecified(msr->prm,"bParaRead"))
		FDL_read(fdl,"bParaRead",&msr->param.bParaRead);
	if (!prmSpecified(msr->prm,"bParaWrite"))
		FDL_read(fdl,"bParaWrite",&msr->param.bParaWrite);
	/*
	 ** Checkpoints can NOT be switched to a different coordinate system!
	 */
	FDL_read(fdl,"bCannonical",&msr->param.bCannonical);
	if (!prmSpecified(msr->prm,"bStandard"))
		FDL_read(fdl,"bStandard",&msr->param.bStandard);
	FDL_read(fdl,"bKDK",&msr->param.bKDK);
	if (!prmSpecified(msr->prm,"nBucket"))
		FDL_read(fdl,"nBucket",&msr->param.nBucket);
	if (!prmSpecified(msr->prm,"iOutInterval"))
		FDL_read(fdl,"iOutInterval",&msr->param.iOutInterval);
	if (!prmSpecified(msr->prm,"iLogInterval"))
		FDL_read(fdl,"iLogInterval",&msr->param.iLogInterval);
	if (!prmSpecified(msr->prm,"iCheckInterval"))
		FDL_read(fdl,"iCheckInterval",&msr->param.iCheckInterval);
	if (!prmSpecified(msr->prm,"iOrder"))
		FDL_read(fdl,"iExpOrder",&msr->param.iOrder);
	if (!prmSpecified(msr->prm,"iEwOrder"))
		FDL_read(fdl,"iEwOrder",&msr->param.iEwOrder);
	if (!prmSpecified(msr->prm,"nReplicas"))
		FDL_read(fdl,"nReplicas",&msr->param.nReplicas);
	if (!prmSpecified(msr->prm,"nSteps"))
		FDL_read(fdl,"nSteps",&msr->param.nSteps);
	if (!prmSpecified(msr->prm,"dExtraStore"))
		FDL_read(fdl,"dExtraStore",&msr->param.dExtraStore);
	if (!prmSpecified(msr->prm,"dDelta"))
		FDL_read(fdl,"dDelta",&msr->param.dDelta);
	if (!prmSpecified(msr->prm,"dEta"))
		FDL_read(fdl,"dEta",&msr->param.dEta);
	if (!prmSpecified(msr->prm,"bEpsVel"))
		FDL_read(fdl,"bEpsVel",&msr->param.bEpsVel);
	if (!prmSpecified(msr->prm,"bNonSymp"))
		FDL_read(fdl,"bNonSymp",&msr->param.bNonSymp);
	if (!prmSpecified(msr->prm,"iMaxRung"))
		FDL_read(fdl,"iMaxRung",&msr->param.iMaxRung);
	if (!prmSpecified(msr->prm,"dEwCut"))
		FDL_read(fdl,"dEwCut",&msr->param.dEwCut);
	if (!prmSpecified(msr->prm,"dEwhCut"))
		FDL_read(fdl,"dEwhCut",&msr->param.dEwhCut);
	FDL_read(fdl,"dPeriod",&msr->param.dPeriod);
	FDL_read(fdl,"dHubble0",&msr->param.dHubble0);
	FDL_read(fdl,"dOmega0",&msr->param.dOmega0);
	if (msr->param.bVerbose) {
		printf("Reading checkpoint file...\nN:%d Time:%g\n",msr->N,dTime);
		}
	in.nStart = 0;
	in.nEnd = msr->N - 1;
	in.iOrder = msr->param.iOrder;
	/*
	 ** Since pstReadCheck causes the allocation of the local particle
	 ** store, we need to tell it the percentage of extra storage it
	 ** should allocate for load balancing differences in the number of
	 ** particles.
	 */
	in.fExtraStore = msr->param.dExtraStore;
	/*
	 ** Provide the period.
	 */
	for (j=0;j<3;++j) {
		in.fPeriod[j] = msr->param.dPeriod;
		}
	in.iVersion = iVersion;
	in.iOffset = FDL_offset(fdl,"particle_array");
	FDL_finish(fdl);
	if(msr->param.bParaRead)
	    pstReadCheck(msr->pst,&in,sizeof(in),NULL,NULL);
	else
	    msrOneNodeReadCheck(msr, &in);
	if (msr->param.bVerbose) puts("Checkpoint file has been successfully read.");
	/*
	 ** Set up the output counter.
	 */
	for (msr->iOut=0;msr->iOut<msr->nOuts;++msr->iOut) {
		if (dTime < msr->pdOutTime[msr->iOut]) break;
		}
	return(dTime);
	}

void msrOneNodeWriteCheck(MSR msr, struct inWriteCheck *in)
{
    int i,id;
    int nStart;
    PST pst0;
    LCL *plcl;
    char achOutFile[PST_FILENAME_SIZE];
    int inswap;

    pst0 = msr->pst;
    while(pst0->nLeaves > 1)
		pst0 = pst0->pstLower;
    plcl = pst0->plcl;
    /*
     ** Add the local Data Path to the provided filename.
     */
    achOutFile[0] = 0;
    if (plcl->pszDataPath) {
	    strcat(achOutFile,plcl->pszDataPath);
	    strcat(achOutFile,"/");
	    }
    strcat(achOutFile,in->achOutFile);

    /* 
     * First write our own particles.
     */
    pkdWriteCheck(plcl->pkd,achOutFile,in->iOffset, 0); 
    nStart = plcl->pkd->nLocal;
	assert(msr->pMap[0] == 0);
    for (i=1;i<msr->nThreads;++i) {
		id = msr->pMap[i];
		/* 
		 * Swap particles with the remote processor.
		 */
		inswap = 0;
		mdlReqService(pst0->mdl,id,PST_SWAPALL,&inswap,sizeof(inswap));
		pkdSwapAll(plcl->pkd, id);
		mdlGetReply(pst0->mdl,id,NULL,NULL);
		/* 
		 * Write the swapped particles.
		 */
		pkdWriteCheck(plcl->pkd,achOutFile,in->iOffset, nStart); 
		nStart += plcl->pkd->nLocal;
		/* 
		 * Swap them back again.
		 */
		inswap = 0;
		mdlReqService(pst0->mdl,id,PST_SWAPALL,&inswap,sizeof(inswap));
		pkdSwapAll(plcl->pkd, id);
		mdlGetReply(pst0->mdl,id,NULL,NULL);
    	}
    assert(nStart == msr->N);
    }

void msrWriteCheck(MSR msr,double dTime,int iStep)
{
	struct inWriteCheck in;
	struct outSetTotal oute;
	char achOutFile[PST_FILENAME_SIZE];
	int iDum,i;
	LCL *plcl = msr->pst->plcl;
	FDL_CTX *fdl;
	char *pszFdl;
	int iVersion,iNotCorrupt;
	static int first = 1;
	
	/*
	 ** Add Data Subpath for local and non-local names.
	 */
	achOutFile[0] = 0;
	strcat(achOutFile,msr->param.achDataSubPath);
	strcat(achOutFile,"/");
	if(first) {
	    sprintf(achOutFile,"%s.chk0",msr->param.achOutName);
	    first = 0;
	} else {
	    sprintf(achOutFile,"%s.chk1",msr->param.achOutName);
	    first = 1;
	}
	strcpy(in.achOutFile,achOutFile);
	/*
	 ** Add local Data Path.
	 */
	achOutFile[0] = 0;
	if (plcl->pszDataPath) {
		strcat(achOutFile,plcl->pszDataPath);
		strcat(achOutFile,"/");
		}
	strcat(achOutFile,in.achOutFile);
	pszFdl = getenv("PKDGRAV_CHECKPOINT_FDL");
	assert(pszFdl != NULL);
	fdl = FDL_create(achOutFile,pszFdl);
	iVersion = CHECKPOINT_VERSION;
	FDL_write(fdl,"version",&iVersion);
	iNotCorrupt = 0;
	FDL_write(fdl,"not_corrupt_flag",&iNotCorrupt);
	/*
	 ** Checkpoint header.
	 */
	FDL_write(fdl,"number_of_particles",&msr->N);
	FDL_write(fdl,"current_timestep",&iStep);
	FDL_write(fdl,"current_time",&dTime);
	FDL_write(fdl,"current_ecosmo",&msr->dEcosmo);
	FDL_write(fdl,"old_time",&msr->dTimeOld);
	FDL_write(fdl,"old_potentiale",&msr->dUOld);
	FDL_write(fdl,"opening_type",&msr->iOpenType);
	FDL_write(fdl,"opening_criterion",&msr->dCrit);
	FDL_write(fdl,"number_of_outs",&msr->nOuts);
	for (i=0;i<msr->nOuts;++i) {
		FDL_index(fdl,"out_time_index",i);
		FDL_write(fdl,"out_time",&msr->pdOutTime[i]);
		}
	/*
	 ** Write the old parameters.
	 */
	FDL_write(fdl,"bPeriodic",&msr->param.bPeriodic);
	FDL_write(fdl,"bComove",&msr->param.bComove);
	FDL_write(fdl,"bParaRead",&msr->param.bParaRead);
	FDL_write(fdl,"bParaWrite",&msr->param.bParaWrite);
	FDL_write(fdl,"bCannonical",&msr->param.bCannonical);
	FDL_write(fdl,"bStandard",&msr->param.bStandard);
	FDL_write(fdl,"bKDK",&msr->param.bKDK);
	FDL_write(fdl,"nBucket",&msr->param.nBucket);
	FDL_write(fdl,"iOutInterval",&msr->param.iOutInterval);
	FDL_write(fdl,"iLogInterval",&msr->param.iLogInterval);
	FDL_write(fdl,"iCheckInterval",&msr->param.iCheckInterval);
	FDL_write(fdl,"iExpOrder",&msr->param.iOrder);
	FDL_write(fdl,"iEwOrder",&msr->param.iEwOrder);
	FDL_write(fdl,"nReplicas",&msr->param.nReplicas);
	FDL_write(fdl,"nSteps",&msr->param.nSteps);
	FDL_write(fdl,"dExtraStore",&msr->param.dExtraStore);
	FDL_write(fdl,"dDelta",&msr->param.dDelta);
	FDL_write(fdl,"dEta",&msr->param.dEta);
	FDL_write(fdl,"bEpsVel",&msr->param.bEpsVel);
	FDL_write(fdl,"bNonSymp",&msr->param.bNonSymp);
	FDL_write(fdl,"iMaxRung",&msr->param.iMaxRung);
	FDL_write(fdl,"dEwCut",&msr->param.dEwCut);
	FDL_write(fdl,"dEwhCut",&msr->param.dEwhCut);
	FDL_write(fdl,"dPeriod",&msr->param.dPeriod);
	FDL_write(fdl,"dHubble0",&msr->param.dHubble0);
	FDL_write(fdl,"dOmega0",&msr->param.dOmega0);
	if (msr->param.bVerbose) {
		printf("Writing checkpoint file...\nTime:%g\n",dTime);
		}
	/*
	 ** Do a parallel or serial write to the output file.
	 */
	pstSetTotal(msr->pst,NULL,0,&oute,&iDum);
	in.iOffset = FDL_offset(fdl,"particle_array");
	in.nStart = 0;
	if(msr->param.bParaWrite)
	    pstWriteCheck(msr->pst,&in,sizeof(in),NULL,NULL);
	else
	    msrOneNodeWriteCheck(msr, &in);
	if (msr->param.bVerbose) {
		puts("Checkpoint file has been successfully written.");
		}
	iNotCorrupt = 1;
	FDL_write(fdl,"not_corrupt_flag",&iNotCorrupt);
	FDL_finish(fdl);
	}


int msrOutTime(MSR msr,double dTime)
{	
	if (msr->iOut < msr->nOuts) {
		if (dTime >= msr->pdOutTime[msr->iOut]) {
			++msr->iOut;
			return(1);
			}
		else return(0);
		}
	else return(0);
	}


int cmpTime(const void *v1,const void *v2) 
{
	double *d1 = (double *)v1;
	double *d2 = (double *)v2;

	if (*d1 < *d2) return(-1);
	else if (*d1 == *d2) return(0);
	else return(1);
	}

void msrReadOuts(MSR msr,double dTime)
{
	char achFile[PST_FILENAME_SIZE];
	char ach[PST_FILENAME_SIZE];
	LCL *plcl = &msr->lcl;
	FILE *fp;
	int i,ret;
	double z,a,n;
	char achIn[80];
	
	/*
	 ** Add Data Subpath for local and non-local names.
	 */
	achFile[0] = 0;
	sprintf(achFile,"%s/%s.red",msr->param.achDataSubPath,
			msr->param.achOutName);
	/*
	 ** Add local Data Path.
	 */
	if (plcl->pszDataPath) {
		strcpy(ach,achFile);
		sprintf(achFile,"%s/%s",plcl->pszDataPath,ach);
		}
	fp = fopen(achFile,"r");
	if (!fp) {
		printf("WARNING: Could not open redshift input file:%s\n",achFile);
		msr->nOuts = 0;
		return;
		}
	i = 0;
	while (1) {
		if (!fgets(achIn,80,fp)) goto NoMoreOuts;
		switch (achIn[0]) {
		case 'z':
			ret = sscanf(&achIn[1],"%lf",&z);
			if (ret != 1) goto NoMoreOuts;
			a = 1.0/(z+1.0);
			msr->pdOutTime[i] = msrExp2Time(msr,a);
			break;
		case 'a':
			ret = sscanf(&achIn[1],"%lf",&a);
			if (ret != 1) goto NoMoreOuts;
			msr->pdOutTime[i] = msrExp2Time(msr,a);
			break;
		case 't':
			ret = sscanf(&achIn[1],"%lf",&msr->pdOutTime[i]);
			if (ret != 1) goto NoMoreOuts;
			break;
		case 'n':
			ret = sscanf(&achIn[1],"%lf",&n);
			if (ret != 1) goto NoMoreOuts;
			msr->pdOutTime[i] = dTime + (n-0.5)*msrDelta(msr);
			break;
		default:
			ret = sscanf(achIn,"%lf",&z);
			if (ret != 1) goto NoMoreOuts;
			a = 1.0/(z+1.0);
			msr->pdOutTime[i] = msrExp2Time(msr,a);
			}
		++i;
		if(i > msr->nMaxOuts) {
			msr->nMaxOuts *= 2;
			msr->pdOutTime = realloc(msr->pdOutTime,
									 msr->nMaxOuts*sizeof(double));
			assert(msr->pdOutTime != NULL);
		    }
		}
 NoMoreOuts:
	msr->nOuts = i;
	/*
	 ** Now sort the array of output times into ascending order.
	 */
	qsort(msr->pdOutTime,msr->nOuts,sizeof(double),cmpTime);
	fclose(fp);
	}


int msrSteps(MSR msr)
{
	return(msr->param.nSteps);
	}


char *msrOutName(MSR msr)
{
	return(msr->param.achOutName);
	}


double msrDelta(MSR msr)
{
	return(msr->param.dDelta);
	}


int msrCheckInterval(MSR msr)
{
	return(msr->param.iCheckInterval);
	}


int msrLogInterval(MSR msr)
{
	return(msr->param.iLogInterval);
	}


int msrOutInterval(MSR msr)
{
	return(msr->param.iOutInterval);
	}


int msrRestart(MSR msr)
{
	return(msr->param.bRestart);
	}


int msrComove(MSR msr)
{
	return(msr->param.bComove);
	}


int msrKDK(MSR msr)
{
	return(msr->param.bCannonical && msr->param.bKDK);
	}


double msrSoft(MSR msr)
{
	return(msr->param.dSoft);
	}


double msrMassCheck(MSR msr,double dMass,char *pszWhere)
{
	struct outMassCheck out;

	if (msr->param.bVerbose) puts("doing mass check...");
	pstMassCheck(msr->pst,NULL,0,&out,NULL);
	if (dMass < 0.0) return(out.dMass);
	else if (dMass != out.dMass) {
		printf("ERROR: Mass not conserved (%s): %.15f != %.15f!\n",
			   pszWhere,dMass,out.dMass);
		_msrExit(msr);
		return(0.0);
		}
	else return(out.dMass);
	}

void
msrInitStep(MSR msr)
{
    struct inSetRung in;

    in.iRung = msr->param.iMaxRung - 1;
    pstSetRung(msr->pst, &in, sizeof(in), NULL, NULL);
    msr->iCurrMaxRung = in.iRung;
    }


int msrMaxRung(MSR msr)
{
    return msr->param.iMaxRung;
    }


int msrCurrMaxRung(MSR msr)
{
    return msr->iCurrMaxRung;
    }


double msrEta(MSR msr)
{
    return msr->param.dEta;
    }


void msrActiveRung(MSR msr, int iRung, int bGreater)
{
    struct inActiveRung in;

    in.iRung = iRung;
    in.bGreater = bGreater;
    pstActiveRung(msr->pst, &in, sizeof(in), NULL, NULL);
    }


int msrCurrRung(MSR msr, int iRung)
{
    struct inCurrRung in;
    struct outCurrRung out;

    in.iRung = iRung;
    pstCurrRung(msr->pst, &in, sizeof(in), &out, NULL);
    return out.iCurrent;
    }


void msrDensityRung(MSR msr, int iRung, double dDelta, double dTime,
		    int bAll)
{
    struct inDensityRung in;
    struct outDensityRung out;
    double expand;
    struct inSmooth inden;
    int sec,dsec;

    inden.nSmooth = msr->param.nSmooth;
    inden.bPeriodic = msr->param.bPeriodic;
	inden.bSymmetric = 0;
	inden.iSmoothType = SMX_DENSITY;
    if (msr->param.bVerbose) printf("Calculating Rung Densities...\n");
    sec = time(0);
    pstSmooth(msr->pst,&inden,sizeof(inden),NULL,NULL);
    dsec = time(0) - sec;
    printf("Rung Densities Calculated, Wallclock:%d secs\n\n",dsec);
    in.iRung = iRung;
    in.dDelta = dDelta;
    in.dEta = msrEta(msr);
    expand = msrTime2Exp(msr, dTime);
    in.dRhoFac = 1.0/(expand*expand*expand);
    in.bAll = bAll;
    pstDensityRung(msr->pst, &in, sizeof(in), &out, NULL);
    msr->iCurrMaxRung = out.iMaxRung;
    }


void msrVelocityRung(MSR msr, int iRung, double dDelta, double dTime, int bAll)
{
    struct inVelocityRung in;
    struct outVelocityRung out;
    double a;

    in.iRung = iRung;
    in.dDelta = dDelta;
    in.dEta = msrEta(msr);
    in.iMaxRung = msrMaxRung(msr);
    a = msrTime2Exp(msr,dTime);
    if(msr->param.bCannonical)
	in.dVelFac = 1.0/(a*a);
    else
	in.dVelFac = 1.0;
    in.dAccFac = 1.0/(a*a*a);
    in.bAll = bAll;
    pstVelocityRung(msr->pst, &in, sizeof(in), &out, NULL);
    msr->iCurrMaxRung = out.iMaxRung;
    }


void msrTopStepDen(MSR msr, double dStep, double dTime, double dDelta, 
				   int iRung, double *pdActiveSum)
{
    double dMass = -1.0;
    int iSec;
    double dWMax, dIMax, dEMax;
	int nActive;

	if(msrCurrMaxRung(msr) >= iRung) { /* particles to be kicked? */
	    if(iRung < msrMaxRung(msr)-1) {
			if (msr->param.bVerbose) {
				printf("Adjust, iRung: %d\n", iRung);
				}
			msrDrift(msr,dTime,0.5*dDelta);
			msrActiveRung(msr, iRung, 1);
			msrBuildTree(msr,0,dMass,1);
			dTime += 0.5*dDelta;
			msrDensityRung(msr,iRung, dDelta, dTime, 0);
			msrDrift(msr,dTime,-0.5*dDelta);
			dTime -= 0.5*dDelta;
			}
		/*
		 ** Actual Stepping.
		 */
		msrTopStepDen(msr, dStep, dTime, 0.5*dDelta,iRung+1,pdActiveSum);
		dStep += 1.0/(2 << iRung);
		if(msrCurrRung(msr, iRung)) {
			if (msr->param.bVerbose) {
			    printf("Kick, iRung: %d\n", iRung);
			    }
			msrActiveRung(msr, iRung, 0);
			msrBuildTree(msr,0,dMass,0);
			msrGravity(msr,dStep,&iSec,&dWMax,&dIMax,&dEMax,&nActive);
			*pdActiveSum += (double)nActive/msr->N;
			msrKick(msr, dTime, dDelta);
			msrRungStats(msr);
			}
		msrTopStepDen(msr,dStep,dTime+0.5*dDelta,0.5*dDelta,iRung+1,
					  pdActiveSum);
		}
	else {    
		if (msr->param.bVerbose) {
		    printf("Drift, iRung: %d\n",iRung-1);
		    }
		msrDrift(msr,dTime,dDelta);
		}
	}


void msrRungStats(MSR msr)
{
	struct inRungStats in;
	struct outRungStats out;
	int i;

	printf("Rung distribution: (");
	for (i=0;i<msr->param.iMaxRung;++i) {
		if (i != 0) printf(",");
		in.iRung = i;
		pstRungStats(msr->pst,&in,sizeof(in),&out,NULL);
		printf("%d",out.nParticles);
		}
	printf(")\n");
	}


void msrTopStepNS(MSR msr, double dStep, double dTime, double dDelta, int
				  iRung, int iAdjust, double *pdActiveSum)
{
    double dMass = -1.0;
    int iSec;
    double dWMax, dIMax, dEMax;
	int nActive;

	if(msrCurrMaxRung(msr) >= iRung) { /* particles to be kicked? */
		if(iAdjust && (iRung < msrMaxRung(msr)-1)) {
			if (msr->param.bVerbose) {
				printf("Adjust, iRung: %d\n", iRung);
				}

			msrActiveRung(msr, iRung, 1);
			msrBuildTree(msr,1,dMass,1);
			msrDensityRung(msr,iRung, dDelta, dTime, 1);
			}
		/*
		 ** Actual Stepping.
		 */
		msrTopStepNS(msr, dStep, dTime, 0.5*dDelta,iRung+1, 0, pdActiveSum);
		dStep += 1.0/(2 << iRung);
		if(msrCurrRung(msr, iRung)) {
			if (msr->param.bVerbose) {
				printf("Kick, iRung: %d\n", iRung);
				}
			msrActiveRung(msr, iRung, 0);
			msrBuildTree(msr,0,dMass,0);
			msrGravity(msr,dStep,&iSec,&dWMax,&dIMax,&dEMax,&nActive);
			*pdActiveSum += (double)nActive/msr->N;
			msrKick(msr, dTime, dDelta);
			}
		msrTopStepNS(msr, dStep, dTime+0.5*dDelta,
					 0.5*dDelta,iRung+1, 1, pdActiveSum);
		}
	else {    
		if (msr->param.bVerbose) {
			printf("Drift, iRung: %d\n", iRung-1);
			}
		msrDrift(msr,dTime,dDelta);
		}
	}


void msrTopStepVel(MSR msr, double dStep, double dTime, double dDelta, int
				   iRung, int iAdjust, double *pdActiveSum)
{
    double dMass = -1.0;
    int iSec;
    double dWMax, dIMax, dEMax;
	int nActive;

	if(msrCurrMaxRung(msr) >= iRung) { /* particles to be kicked? */
		if(iAdjust && (iRung < msrMaxRung(msr)-1)) {
			if (msr->param.bVerbose) {
				printf("Adjust, iRung: %d\n", iRung);
				}

			msrActiveRung(msr, iRung, 1);
			msrVelocityRung(msr,iRung, dDelta, dTime, 1);
			}
		/*
		 ** Actual Stepping.
		 */
		msrTopStepVel(msr, dStep, dTime, 0.5*dDelta,iRung+1, 0, pdActiveSum);
		dStep += 1.0/(2 << iRung);
		if(msrCurrRung(msr, iRung)) {
			if (msr->param.bVerbose) {
				printf("Kick, iRung: %d\n", iRung);
				}
			msrActiveRung(msr, iRung, 0);
			msrBuildTree(msr,0,dMass,0);
			msrGravity(msr,dStep,&iSec,&dWMax,&dIMax,&dEMax,&nActive);
			*pdActiveSum += (double)nActive/msr->N;
			msrKick(msr, dTime, dDelta);
			}
		msrTopStepVel(msr, dStep, dTime+0.5*dDelta, 0.5*dDelta,
					  iRung+1, 1, pdActiveSum);
		}
	else {    
		if (msr->param.bVerbose) {
			printf("Drift, iRung: %d\n", iRung-1);
			}
		msrDrift(msr,dTime,dDelta);
		}
	}


void msrTopStep(MSR msr, double dStep, double dTime, double dDelta, 
				double *pdMultiEff)
{
	int iRung = 0;

	*pdMultiEff = 0.0;
    if(msr->param.bEpsVel)
		msrTopStepVel(msr,dStep,dTime,dDelta,iRung,1,pdMultiEff);
    else {
		if(msr->param.bNonSymp)
			msrTopStepNS(msr,dStep,dTime,dDelta,iRung,1,pdMultiEff);
		else
			msrTopStepDen(msr,dStep,dTime,dDelta,iRung,pdMultiEff);
		}
	printf("Multistep Efficiency (average number of microsteps per step):%f\n",
		   *pdMultiEff);
	}


int msrDoDensity(MSR msr)
{
	return(msr->param.bDoDensity);
	}

