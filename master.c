#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <malloc.h>
#include <assert.h>
#include <time.h>
#include <math.h>

#ifndef CRAY_T3D
#include <rpc/types.h>
#include <rpc/xdr.h>
#else
#include "hyperlib.h"
#endif

#include "master.h"
#include "tipsydefs.h"
#include "opentype.h"
#include "checkdefs.h"


void _msrLeader(void)
{
    puts("USAGE: main [SETTINGS | FLAGS] [SIM_FILE]");
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


void msrInitialize(MSR *pmsr,MDL mdl,int argc,char **argv,char *pszDefaultName)
{
	MSR msr;
	int j,ret;
	int id;
	struct inSetAdd inAdd;
	struct inLevelize inLvl;
	
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
	prmAddParam(msr->prm,"nThreads",1,&msr->param.nThreads,"sz",
				"<nThreads>");
	msr->param.bDiag = 0;
	prmAddParam(msr->prm,"bDiag",0,&msr->param.bDiag,"d",
				"enable/disable per thread diagnostic output");
	msr->param.bVerbose = 0;
	prmAddParam(msr->prm,"bVerbose",0,&msr->param.bVerbose,"v",
				"enable/disable verbose output");
	msr->param.bPeriodic = 0;
	prmAddParam(msr->prm,"bPeriodic",0,&msr->param.bPeriodic,"p",
				"periodic/non-periodic = -p");
	msr->param.bRestart = 0;
	prmAddParam(msr->prm,"bRestart",0,&msr->param.bRestart,"restart",
				"restart from checkpoint");
	msr->param.bParaRead = 1;
	prmAddParam(msr->prm,"bParaRead",0,&msr->param.bParaRead,"par",
				"enable/disable parallel reading of files = +par");
	msr->param.bParaWrite = 1;
	prmAddParam(msr->prm,"bParaWrite",0,&msr->param.bParaWrite,"paw",
				"enable/disable parallel writing of files = +paw");
	msr->param.nBucket = 8;
	prmAddParam(msr->prm,"nBucket",1,&msr->param.nBucket,"b",
				"<max number of particles in a bucket> = 8");
	msr->param.nSteps = 0;
	prmAddParam(msr->prm,"nSteps",1,&msr->param.nSteps,"n",
				"<number of timesteps> = 0");
	msr->param.iOutInterval = 20;
	prmAddParam(msr->prm,"iOutInterval",1,&msr->param.iOutInterval,
				"oi","<number of timesteps between snapshots> = 20");
	msr->param.iLogInterval = 10;
	prmAddParam(msr->prm,"iLogInterval",1,&msr->param.iLogInterval,
				"ol","<number of timesteps between logfile outputs> = 10");
	msr->param.iCheckInterval = 10;
	prmAddParam(msr->prm,"iCheckInterval",1,&msr->param.iCheckInterval,
				"oc","<number of timesteps between checkpoints> = 10");
	msr->param.iOrder = 4;
	prmAddParam(msr->prm,"iOrder",1,&msr->param.iOrder,"or",
				"<multipole expansion order: 1, 2, 3 or 4> = 4");
	msr->param.iEwOrder = 4;
	prmAddParam(msr->prm,"iEwOrder",1,&msr->param.iEwOrder,"ewo",
				"<Ewald multipole expansion order: 1, 2, 3 or 4> = 4");
	msr->param.nReplicas = 0;
	prmAddParam(msr->prm,"nReplicas",1,&msr->param.nReplicas,"nrep",
				"<nReplicas> = 0 for -p, or 1 for +p");
	msr->param.dSoft = 0.0;
	prmAddParam(msr->prm,"dSoft",2,&msr->param.dSoft,"e",
				"<gravitational softening length> = 0.0");
	msr->param.dDelta = 0.0;
	prmAddParam(msr->prm,"dDelta",2,&msr->param.dDelta,"dt",
				"<time step>");
	msr->param.dEwCut = 2.6;
	prmAddParam(msr->prm,"dEwCut",2,&msr->param.dEwCut,"ew",
				"<dEwCut> = 2.6");
	msr->param.dEwhCut = 2.8;
	prmAddParam(msr->prm,"dEwhCut",2,&msr->param.dEwhCut,"ewh",
				"<dEwhCut> = 2.8");
	msr->param.dTheta = 0.8;
	prmAddParam(msr->prm,"dTheta",2,&msr->param.dTheta,"theta",
				"<Barnes opening criterion> = 0.8");
	msr->param.dAbsPartial = 0.0;
	prmAddParam(msr->prm,"dAbsPartial",2,&msr->param.dAbsPartial,"ap",
				"<absolute partial error opening criterion>");
	msr->param.dRelPartial = 0.0;
	prmAddParam(msr->prm,"dRelPartial",2,&msr->param.dRelPartial,"rp",
				"<relative partial error opening criterion>");
	msr->param.dAbsTotal = 0.0;
	prmAddParam(msr->prm,"dAbsTotal",2,&msr->param.dAbsTotal,"at",
				"<absolute total error opening criterion>");
	msr->param.dRelTotal = 0.0;
	prmAddParam(msr->prm,"dRelTotal",2,&msr->param.dRelTotal,"rt",
				"<relative total error opening criterion>");
	msr->param.dPeriod = 1.0;
	prmAddParam(msr->prm,"dPeriod",2,&msr->param.dPeriod,"L",
				"<periodic box length> = 1.0");
	msr->param.achInFile[0] = 0;
	prmAddParam(msr->prm,"achInFile",3,msr->param.achInFile,"I",
				"<input file name> (file in TIPSY binary format)");
	strcpy(msr->param.achOutName,pszDefaultName); 
	prmAddParam(msr->prm,"achOutName",3,&msr->param.achOutName,"o",
				"<output name for snapshots and logfile> = \"pkdgrav\"");
	msr->param.bComove = 0;
	prmAddParam(msr->prm,"bComove",0,&msr->param.bComove,"cm",
				"enable/disable comoving coordinates = -cm");
	msr->param.dHubble0 = 0.0;
	prmAddParam(msr->prm,"dHubble0",2,&msr->param.dHubble0,"Hub",
				"<dHubble0> = 0.0");
	msr->param.dOmega0 = 1.0;
	prmAddParam(msr->prm,"dOmega0",2,&msr->param.dOmega0,"Om",
				"<dOmega0> = 1.0");
	strcpy(msr->param.achDataSubPath,".");
	prmAddParam(msr->prm,"achDataSubPath",3,&msr->param.achDataSubPath,
				NULL,NULL);
	msr->param.dExtraStore = 1.0;
	prmAddParam(msr->prm,"dExtraStore",2,&msr->param.dExtraStore,NULL,NULL);
	msr->param.nSmooth = 64;
	prmAddParam(msr->prm,"nSmooth",1,&msr->param.nSmooth,"s",
				"<number of particles to smooth over> = 64");
	msr->param.bGatherScatter = 1;
	prmAddParam(msr->prm,"bGatherScatter",0,&msr->param.bGatherScatter,"gs",
				"gather-scatter/gather-only smoothing kernel = +gs");
	/*
	 ** Added bStandard as a new parameter, but can't add it to 
	 ** msr->params because this would change the checkpoint 
	 ** format. Need to generalize the checkpoint format to handle
	 ** these types of changes!
	 ** Now also dRedTo...
	 */
	msr->bStandard = 0;
	prmAddParam(msr->prm,"bStandard",0,&msr->bStandard,"std",
				"output in standard TIPSY binary format");
	msr->dRedTo = 0.0;	
	prmAddParam(msr->prm,"dRedTo",2,&msr->dRedTo,"zto",
				"specifies final redshift for the simulation");
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
		msrFinish(msr);
		mdlFinish(mdl);
		exit(1);
		}
	if (msr->param.bPeriodic && !prmSpecified(msr->prm,"nReplicas")) {
		msr->param.nReplicas = 1;
		}
	if (!msr->param.achInFile[0] && !msr->param.bRestart) {
		printf("pkdgrav ERROR: no input file specified\n");
		msrFinish(msr);
		mdlFinish(mdl);
		exit(1);
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
	msr->iOpenType = OPEN_JOSH;
	if (prmFileSpecified(msr->prm,"dAbsPartial")) msr->iOpenType = OPEN_ABSPAR;
	if (prmFileSpecified(msr->prm,"dRelPartial")) msr->iOpenType = OPEN_RELPAR;
	if (prmFileSpecified(msr->prm,"dAbsTotal")) msr->iOpenType = OPEN_ABSTOT;
	if (prmFileSpecified(msr->prm,"dRelTotal")) msr->iOpenType = OPEN_RELTOT;
	if (prmArgSpecified(msr->prm,"dTheta")) msr->iOpenType = OPEN_JOSH;
	if (prmArgSpecified(msr->prm,"dAbsPartial")) msr->iOpenType = OPEN_ABSPAR;
	if (prmArgSpecified(msr->prm,"dRelPartial")) msr->iOpenType = OPEN_RELPAR;
	if (prmArgSpecified(msr->prm,"dAbsTotal")) msr->iOpenType = OPEN_ABSTOT;
	if (prmArgSpecified(msr->prm,"dRelTotal")) msr->iOpenType = OPEN_RELTOT;
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
	msr->dRedshift = 0.0;
	msr->dCosmoFac = 1.0;
	msr->dHubble = 0.0;
	msr->nRed = 0;
	if (msrComove(msr)) msrReadRed(msr);
        
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
	}


void msrLogParams(MSR msr,FILE *fp)
{
	int i;
  
	fprintf(fp,"# nThreads: %d",msr->param.nThreads);
	fprintf(fp," bDiag: %d",msr->param.bDiag);
	fprintf(fp," bVerbose: %d",msr->param.bVerbose);
	fprintf(fp," bPeriodic: %d",msr->param.bPeriodic);
	fprintf(fp," bRestart: %d",msr->param.bRestart);
	fprintf(fp," nBucket: %d",msr->param.nBucket);
	fprintf(fp,"\n# nSteps: %d",msr->param.nSteps);
	fprintf(fp," iOutInterval: %d",msr->param.iOutInterval);
	fprintf(fp," iLogInterval: %d",msr->param.iLogInterval);
	fprintf(fp," iCheckInterval: %d",msr->param.iCheckInterval);
	fprintf(fp," iOrder: %d",msr->param.iOrder);
	fprintf(fp,"\n# nReplicas: %d",msr->param.nReplicas);
	if (prmArgSpecified(msr->prm,"dSoft"))
		fprintf(fp," dSoft: %g",msr->param.dSoft);
	else
		fprintf(fp," dSoft: input");
	fprintf(fp," dDelta: %g",msr->param.dDelta);
	fprintf(fp," dEwCut: %f",msr->param.dEwCut);
	fprintf(fp," dEwhCut: %f",msr->param.dEwhCut);
	fprintf(fp,"\n#");
	switch (msr->iOpenType) {
	case OPEN_JOSH:
		fprintf(fp," iOpenType: JOSH");
		break;
	case OPEN_ABSPAR:
		fprintf(fp," iOpenType: ABSPAR");
		break;
	case OPEN_RELPAR:
		fprintf(fp," iOpenType: RELPAR");
		break;
	case OPEN_ABSTOT:
		fprintf(fp," iOpenType: ABSTOT");
		break;
	case OPEN_RELTOT:
		fprintf(fp," iOpenType: RELTOT");
		break;
	default:
		fprintf(fp," iOpenType: NONE?");
		}
	fprintf(fp," dTheta: %f",msr->param.dTheta);
	fprintf(fp," dAbsPartial: %g",msr->param.dAbsPartial);
	fprintf(fp," dRealPartial: %g",msr->param.dRelPartial);
	fprintf(fp," dAbsTotal: %g",msr->param.dAbsTotal);
	fprintf(fp," dRelTotal: %g",msr->param.dRelTotal);
	fprintf(fp,"\n# dPeriod: %g",msr->param.dPeriod);
	fprintf(fp," achInFile: %s",msr->param.achInFile);
	fprintf(fp," achOutName: %s",msr->param.achOutName); 
	fprintf(fp," bComove: %d",msr->param.bComove);
	fprintf(fp," dHubble0: %g",msr->param.dHubble0);
	fprintf(fp," dOmega0: %g",msr->param.dOmega0);
	fprintf(fp,"\n# achDataSubPath: %s",msr->param.achDataSubPath);
	fprintf(fp," dExtraStore: %f",msr->param.dExtraStore);
	fprintf(fp," nSmooth: %d",msr->param.nSmooth);
	fprintf(fp," bGatherScatter: %d\n",msr->param.bGatherScatter);
	fprintf(fp,"# RedOut:");
	for (i=0;i<msr->nRed;i++)
		fprintf(fp," %f",msrRedOut(msr, i));
	fprintf(fp,"\n");
	fflush(fp);
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


double msrReadTipsy(MSR msr)
{
	FILE *fp;
	struct dump h;
	struct inReadTipsy in;
	char achInFile[PST_FILENAME_SIZE];
	int j;
	LCL *plcl = msr->pst->plcl;
	double dTime,aTo,tTo;
	
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
			msrFinish(msr);
			mdlFinish(msr->mdl);
			exit(1);
			}
		}
	else {
		printf("No input file specified\n");
		msrFinish(msr);
		mdlFinish(msr->mdl);
		exit(1);
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
			msrFinish(msr);
			mdlFinish(msr->mdl);
			exit(1);
			}
		dTime = msrExp2Time(msr,h.time);
		msr->dRedshift = 1.0/h.time - 1.0;
		printf("Input file, Time:%g Redshift:%g Expansion factor:%g\n",
			   dTime,msr->dRedshift,h.time);
		if (prmSpecified(msr->prm,"dRedTo")) {
			if (!prmArgSpecified(msr->prm,"nSteps") &&
				prmArgSpecified(msr->prm,"dDelta")) {
				aTo = 1.0/(msr->dRedTo + 1.0);
				tTo = msrExp2Time(msr,aTo);
				printf("Simulation to Time:%g Redshift:%g Expansion factor:%g\n",
					   tTo,1.0/aTo-1.0,aTo);
				if (tTo < dTime) {
					printf("Badly specified final redshift, check -zto parameter.\n");
					msrFinish(msr);
					mdlFinish(msr->mdl);
					exit(1);
					}
				msr->param.nSteps = (int)ceil((tTo-dTime)/msr->param.dDelta);
				}
			else if (!prmArgSpecified(msr->prm,"dDelta") &&
					 prmArgSpecified(msr->prm,"nSteps")) {
				aTo = 1.0/(msr->dRedTo + 1.0);
				tTo = msrExp2Time(msr,aTo);
				printf("Simulation to Time:%g Redshift:%g Expansion factor:%g\n",
					   tTo,1.0/aTo-1.0,aTo);
				if (tTo < dTime) {
					printf("Badly specified final redshift, check -zto parameter.\n");	
					msrFinish(msr);
					mdlFinish(msr->mdl);
					exit(1);
					}
				msr->param.dDelta = (tTo-dTime)/msr->param.nSteps;
				}
			else if (!prmSpecified(msr->prm,"nSteps") &&
				prmFileSpecified(msr->prm,"dDelta")) {
				aTo = 1.0/(msr->dRedTo + 1.0);
				tTo = msrExp2Time(msr,aTo);
				printf("Simulation to Time:%g Redshift:%g Expansion factor:%g\n",
					   tTo,1.0/aTo-1.0,aTo);
				if (tTo < dTime) {
					printf("Badly specified final redshift, check -zto parameter.\n");
					msrFinish(msr);
					mdlFinish(msr->mdl);
					exit(1);
					}
				msr->param.nSteps = (int)ceil((tTo-dTime)/msr->param.dDelta);
				}
			else if (!prmSpecified(msr->prm,"dDelta") &&
					 prmFileSpecified(msr->prm,"nSteps")) {
				aTo = 1.0/(msr->dRedTo + 1.0);
				tTo = msrExp2Time(msr,aTo);
				printf("Simulation to Time:%g Redshift:%g Expansion factor:%g\n",
					   tTo,1.0/aTo-1.0,aTo);
				if (tTo < dTime) {
					printf("Badly specified final redshift, check -zto parameter.\n");	
					msrFinish(msr);
					mdlFinish(msr->mdl);
					exit(1);
					}
				msr->param.dDelta = (tTo-dTime)/msr->param.nSteps;
				}
			}
		else {
			tTo = dTime + msr->param.nSteps*msr->param.dDelta;
			aTo = msrTime2Exp(msr,tTo);
			printf("Simulation to Time:%g Redshift:%g Expansion factor:%g\n",
				   tTo,1.0/aTo-1.0,aTo);
			}
		printf("Reading file...\nN:%d\n",msr->N);
		}
	else {
		dTime = h.time;
		msr->dRedshift = 0.0;
		printf("Input file, Time:%g\n",dTime);
		tTo = dTime + msr->param.nSteps*msr->param.dDelta;
		printf("Simulation to Time:%g\n",tTo);
		printf("Reading file...\nN:%d Time:%g\n",msr->N,dTime);
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
	in.bParaRead = msr->param.bParaRead;
	pstReadTipsy(msr->pst,&in,sizeof(in),NULL,NULL);
	if (msr->param.bVerbose) puts("Input file has been successfully read.");
	return(dTime);
	}


#ifndef CRAY_T3D
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
#endif


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
		msrFinish(msr);
		mdlFinish(msr->mdl);
		exit(1);
		}
	in.bStandard = msr->bStandard;
	/*
	 ** Assume tipsy format for now, and dark matter only.
	 */
	h.nbodies = msr->N;
	h.ndark = msr->N;
	h.nsph = 0;
	h.nstar = 0;
	if (msr->param.bComove)
	  h.time = msrTime2Exp(msr,dTime);
	else
	  h.time = dTime;
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
#ifndef CRAY_T3D
		XDR xdrs;

		assert(sizeof(struct dump) == sizeof(double)+6*sizeof(int));
		xdrstdio_create(&xdrs,fp,XDR_ENCODE);
		xdrWriteHeader(&xdrs,&h);
		xdr_destroy(&xdrs);
#endif
		}
	else {
		fwrite(&h,sizeof(struct dump),1,fp);
		}
	fclose(fp);
	in.bParaWrite = msr->param.bParaWrite;
	pstWriteTipsy(msr->pst,&in,sizeof(in),NULL,NULL);
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


void msrBuildTree(MSR msr)
{
	struct inBuildTree in;
	struct outBuildTree out;
	struct inColCells inc;
	struct ioCalcRoot root;
	KDN *pkdn;
	int sec,dsec;
	int iDum,nCell;

	if (msr->param.bVerbose) printf("Domain Decomposition...\n");
	sec = time(0);
	pstDomainDecomp(msr->pst,NULL,0,NULL,NULL);
	dsec = time(0) - sec;
	if (msr->param.bVerbose) {
		printf("Domain Decomposition complete, Wallclock: %d secs\n\n",dsec);
		}
	if (msr->param.bVerbose) printf("Building local trees...\n");
	in.nBucket = msr->param.nBucket;
	in.iOpenType = msr->iOpenType;
	in.iOrder = (msr->param.iOrder >= msr->param.iEwOrder)?
		msr->param.iOrder:msr->param.iEwOrder;
	in.dCrit = msr->dCrit;
	sec = time(0);
	pstBuildTree(msr->pst,&in,sizeof(in),&out,&iDum);
	dsec = time(0) - sec;
	if (msr->param.bVerbose) {
		printf("Tree built, Wallclock: %d secs\n\n",dsec);
		}
	nCell = 1<<(1+(int)ceil(log(msr->nThreads)/log(2.0)));
	pkdn = malloc(nCell*sizeof(KDN));
	assert(pkdn != NULL);
	inc.iCell = ROOT;
	inc.nCell = nCell;
	pstColCells(msr->pst,&inc,sizeof(inc),pkdn,NULL);
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
	free(pkdn);
	pstCalcRoot(msr->pst,NULL,0,&root,&iDum);
	pstDistribRoot(msr->pst,&root,sizeof(struct ioCalcRoot),NULL,NULL);
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
	}


void msrOutArray(MSR msr,char *pszFile,int iType)
{
	struct inOutArray in;
	char achOutFile[PST_FILENAME_SIZE];
	LCL *plcl = msr->pst->plcl;
	FILE *fp;

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
			msrFinish(msr);
			mdlFinish(msr->mdl);
			exit(1);
			}
		}
	else {
		printf("No Array Output File specified\n");
		msrFinish(msr);
		mdlFinish(msr->mdl);
		exit(1);
		}
	/*
	 ** Write the Header information and close the file again.
	 */
	fprintf(fp,"%d\n",msr->N);
	fclose(fp);
	in.iType = iType;
	pstOutArray(msr->pst,&in,sizeof(in),NULL,NULL);
	}


void msrOutVector(MSR msr,char *pszFile,int iType)
{
	struct inOutVector in;
	char achOutFile[PST_FILENAME_SIZE];
	LCL *plcl = msr->pst->plcl;
	FILE *fp;

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
			msrFinish(msr);
			mdlFinish(msr->mdl);
			exit(1);
			}
		}
	else {
		printf("No Vector Output File specified\n");
		msrFinish(msr);
		mdlFinish(msr->mdl);
		exit(1);
		}
	/*
	 ** Write the Header information and close the file again.
	 */
	fprintf(fp,"%d\n",msr->N);
	fclose(fp);
	in.iType = iType;
	in.iDim = 0;
	pstOutVector(msr->pst,&in,sizeof(in),NULL,NULL);
	in.iDim = 1;
	pstOutVector(msr->pst,&in,sizeof(in),NULL,NULL);
	in.iDim = 2;
	pstOutVector(msr->pst,&in,sizeof(in),NULL,NULL);
	}


void msrDensity(MSR msr)
{
	struct inDensity in;
	int sec,dsec;

	in.nSmooth = msr->param.nSmooth;
	in.bGatherScatter = msr->param.bGatherScatter;
	if (msr->param.bVerbose) printf("Calculating Densities...\n");
	sec = time(0);
	pstDensity(msr->pst,&in,sizeof(in),NULL,NULL);
	dsec = time(0) - sec;
	printf("Densities Calculated, Wallclock:%d secs\n\n",dsec);
	}


void msrGravity(MSR msr,int *piSec,double *pdWMax,double *pdIMax,
				double *pdEMax)
{
	struct inGravity in;
	struct outGravity out;
	int iDum;
	int sec,dsec;
	double dPartAvg,dCellAvg;
	double dWAvg,dWMax,dWMin;
	double dIAvg,dIMax,dIMin;
	double dEAvg,dEMax,dEMin;
	double iP;

	if (msr->param.bVerbose) printf("Calculating Gravity...\n");
	sec = time(0);
    in.nReps = msr->param.nReplicas;
    in.bPeriodic = msr->param.bPeriodic;
	in.iOrder = msr->param.iOrder;
	in.iEwOrder = msr->param.iEwOrder;
    in.dEwCut = msr->param.dEwCut;
    in.dEwhCut = msr->param.dEwhCut;
	pstGravity(msr->pst,&in,sizeof(in),&out,&iDum);
	dsec = time(0) - sec;
	printf("Gravity Calculated, Wallclock:%d secs\n",dsec);
	*piSec = dsec;
	dPartAvg = out.dPartSum/msr->N;
	dCellAvg = out.dCellSum/msr->N;
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


void msrStepCosmo(MSR msr,double dTime)
{
	double a;
	
	a = msrTime2Exp(msr,dTime);
	msr->dRedshift = 1.0/a - 1.0;
	msr->dHubble = msr->param.dHubble0*(1.0+msr->dRedshift)*
		sqrt(1.0+msr->param.dOmega0*msr->dRedshift);
	if (msr->param.bComove) msr->dCosmoFac = a;	
	}


void msrCalcE(MSR msr,int bFirst,double dTime,double *E,double *T,double *U)
{
	struct outCalcE out;

	pstCalcE(msr->pst,NULL,0,&out,NULL);
	*T = out.T;
	*U = out.U;
	/*
	 ** Do the comoving coordinates stuff.
	 */
	*U *= msr->dCosmoFac;
	*T *= pow(msr->dCosmoFac,4.0);
	if (msr->param.bComove && !bFirst) {
		msr->dEcosmo += 0.5*(dTime - msr->dTimeOld)*
			(msr->dHubbleOld*msr->dUOld + msr->dHubble*(*U));
		}
	else {
		msr->dEcosmo = 0.0;
		}
	msr->dTimeOld = dTime;
	msr->dUOld = *U;
	msr->dHubbleOld = msr->dHubble;
	*E = (*T) + (*U) - msr->dEcosmo;
	}


void msrDrift(MSR msr,double dDelta)
{
	struct inDrift in;
	int j;

	in.dDelta = dDelta;
	for (j=0;j<3;++j) {
		in.fCenter[j] = msr->fCenter[j];
		}
	in.bPeriodic = msr->param.bPeriodic;
	pstDrift(msr->pst,&in,sizeof(in),NULL,NULL);
	}


void msrKick(MSR msr,double dDelta)
{
	struct inKick in;
	
	in.dvFacOne = (1.0 - msr->dHubble*dDelta)/(1.0 + msr->dHubble*dDelta);
	in.dvFacTwo = dDelta/pow(msr->dCosmoFac,3.0)/(1.0 + msr->dHubble*dDelta);
	pstKick(msr->pst,&in,sizeof(in),NULL,NULL);
	}


double msrReadCheck(MSR msr,int *piStep)
{
	FILE *fp;
	struct msrCheckPointHeader h;
	struct inReadCheck in;
	char achInFile[PST_FILENAME_SIZE];
	int i,j,bNewCheck;
	LCL *plcl = msr->pst->plcl;
	double dTime;
	
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
	fp = fopen(achInFile,"r");
	if (!fp) {
		printf("Could not open checkpoint file:%s\n",achInFile);
		/*
		 ** Try opening a .ochk
		 ** Add Data Subpath for local and non-local names.
		 */
		sprintf(achInFile,"%s/%s.ochk",msr->param.achDataSubPath,
				msr->param.achOutName);
		strcpy(in.achInFile,achInFile);
		/*
		 ** Add local Data Path.
		 */
		if (plcl->pszDataPath) {
			sprintf(achInFile,"%s/%s",plcl->pszDataPath,in.achInFile);
			}
		fp = fopen(achInFile,"r");
		if (!fp) {
			printf("Could not open checkpoint file:%s\n",achInFile);
			msrFinish(msr);
			mdlFinish(msr->mdl);
			exit(1);
			}
		bNewCheck = 0;
		}
	else {
		bNewCheck = 1;
		}
	fread(&h,sizeof(struct msrCheckPointHeader),1,fp);
	fclose(fp);
	dTime = h.dTime;
	*piStep = h.iStep;
	msr->N = h.N;
	msr->param = h.param;
	msr->iOpenType = h.iOpenType;
	msr->dCrit = h.dCrit;
	msr->dRedshift = h.dRedshift;
	msr->dHubble = h.dHubble;
	msr->dCosmoFac = h.dCosmoFac;
	msr->dEcosmo = h.dEcosmo;
	msr->dHubbleOld = h.dHubbleOld;
	msr->dUOld = h.dUOld;
	msr->dTimeOld = h.dTimeOld;
	msr->nRed = h.nRed;
	for (i=0;i<msr->nRed;++i) {
		msr->dRedOut[i] = h.dRedOut[i];
		}
	if (msr->param.bVerbose) {
		if (msr->param.bComove) {
			printf("Reading checkpoint file...\nN:%d Time:%g Redshift:%g\n",
				   msr->N,dTime,msr->dRedshift);
			}
		else {
			printf("Reading checkpoint file...\nN:%d Time:%g\n",msr->N,dTime);
			}
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
	in.bNewCheck = bNewCheck;
	in.bParaRead = msr->param.bParaRead;
	pstReadCheck(msr->pst,&in,sizeof(in),NULL,NULL);
	if (msr->param.bVerbose) puts("Checkpoint file has been successfully read.");
	return(dTime);
	}


void msrWriteCheck(MSR msr,double dTime,int iStep)
{
	FILE *fp;
	struct msrCheckPointHeader h;
	struct inWriteCheck in;
	struct outSetTotal oute;
	char achOutFile[PST_FILENAME_SIZE];
	char achOutTmp[PST_FILENAME_SIZE];
	char achOutBak[PST_FILENAME_SIZE];
	char ach[4*PST_FILENAME_SIZE];
	int iDum,i;
	LCL *plcl = msr->pst->plcl;
	
	/*
	 ** Add Data Subpath for local and non-local names.
	 */
	achOutFile[0] = 0;
	strcat(achOutFile,msr->param.achDataSubPath);
	strcat(achOutFile,"/");
	sprintf(achOutFile,"%s.chk",msr->param.achOutName);
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
	sprintf(achOutTmp,"%s%s",achOutFile,".tmp");
	sprintf(achOutBak,"%s%s",achOutFile,".bak");
#ifdef SAFE_CHECK
	strcat(in.achOutFile,".tmp");	
	fp = fopen(achOutTmp,"w");
#else
	fp = fopen(achOutFile,"w");
#endif
	if (!fp) {
		printf("Could not open OutFile:%s\n",achOutTmp);
		msrFinish(msr);
		mdlFinish(msr->mdl);
		exit(1);
		}
	h.dTime = dTime;
	h.iStep = iStep;
	h.N = msr->N;
	h.param = msr->param;
	h.iOpenType = msr->iOpenType;
	h.dCrit = msr->dCrit;
	h.dRedshift = msr->dRedshift;
	h.dHubble = msr->dHubble;
	h.dCosmoFac = msr->dCosmoFac;
	h.dEcosmo = msr->dEcosmo;
	h.dHubbleOld = msr->dHubbleOld;
	h.dUOld = msr->dUOld;
	h.dTimeOld = msr->dTimeOld;
	h.nRed = msr->nRed;
	for (i=0;i<msr->nRed;++i) {
		h.dRedOut[i] = msr->dRedOut[i];
		}
	if (msr->param.bVerbose) {
		if (msr->param.bComove) {
			printf("Writing checkpoint file...\nTime:%g Redshift:%g\n",
				   dTime,(1.0/dTime - 1.0));
			}
		else {
			printf("Writing checkpoint file...\nTime:%g\n",dTime);
			}
		}
	fwrite(&h,sizeof(struct msrCheckPointHeader),1,fp);
	fclose(fp);
	/*
	 ** Do a parallel write to the output file.
	 */
	pstSetTotal(msr->pst,NULL,0,&oute,&iDum);
	in.nStart = 0;
	in.bParaWrite = msr->param.bParaWrite;
	pstWriteCheck(msr->pst,&in,sizeof(in),NULL,NULL);
	if (msr->param.bVerbose) {
		puts("Checkpoint file has been successfully written.");
		}
	sprintf(ach,"mv -f %s %s; mv %s %s",achOutFile,achOutBak,achOutTmp,
			achOutFile);
#ifdef SAFE_CHECK
	system(ach);
#endif
	}


double msrRedOut(MSR msr,int iRed)
{
	if (iRed < msr->nRed) {
		return(msr->dRedOut[iRed]);
		}
	else {
		return(-1.0);
		}
	}


void msrReadRed(MSR msr)
{
	char achFile[PST_FILENAME_SIZE];
	char ach[PST_FILENAME_SIZE];
	LCL *plcl = &msr->lcl;
	FILE *fp;
	int iRed,ret;
	
	/*
	 ** Add Data Subpath for local and non-local names.
	 */
	achFile[0] = 0;
	strcat(achFile,msr->param.achDataSubPath);
	strcat(achFile,"/");
	sprintf(achFile,"%s.red",msr->param.achOutName);
	/*
	 ** Add local Data Path.
	 */
	if (plcl->pszDataPath) {
		strcpy(ach,achFile);
		sprintf(achFile,"%s/%s",plcl->pszDataPath,ach);
		}
	fp = fopen(achFile,"r");
	if (!fp) {
		printf("Could not open redshift input file:%s\n",achFile);
		msrFinish(msr);
		mdlFinish(msr->mdl);
		exit(1);
		}
	iRed = 0;
	while (1) {
		ret = fscanf(fp,"%lf",&msr->dRedOut[iRed]);
		if (ret != 1) break;
		++iRed;
		if(iRed > MAX_REDSHIFTS){
		    printf("Too many output redshifts, recompile with greater MAX_REDSHIFTS\n");
		    msrFinish(msr);
		    mdlFinish(msr->mdl);
		    exit(1);
		    }
		}
	msr->nRed = iRed;
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


double msrRedshift(MSR msr)
{
	return(msr->dRedshift);
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


double msrSoft(MSR msr)
{
	return(msr->param.dSoft);
	}













