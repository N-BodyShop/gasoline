#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "mdl.h"
#include "master.h"
#include "outtype.h"


void main_ch(MDL mdl)
{
	PST pst;
	LCL lcl;
	int nThreads;

	lcl.pszDataPath = (char *)getenv("PTOOLS_DATA_PATH");
	pstInitialize(&pst,mdl,&lcl);

	nThreads = mdlThreads(mdl);
	mdlAddService(mdl,PST_SETADD,pst,pstSetAdd,
				  sizeof(struct inSetAdd),0);
	mdlAddService(mdl,PST_LEVELIZE,pst,pstLevelize,
				  sizeof(struct inLevelize),0);
	mdlAddService(mdl,PST_READTIPSY,pst,pstReadTipsy,
				  sizeof(struct inReadTipsy),0);
	mdlAddService(mdl,PST_DOMAINDECOMP,pst,pstDomainDecomp,
				  0,0);
	mdlAddService(mdl,PST_CALCBOUND,pst,pstCalcBound,
				  0,sizeof(struct outCalcBound));
	mdlAddService(mdl,PST_WEIGHT,pst,pstWeight,
				  sizeof(struct inWeight),sizeof(struct outWeight));
	mdlAddService(mdl,PST_FREESTORE,pst,pstFreeStore,
				  0,sizeof(struct outFreeStore));
	mdlAddService(mdl,PST_COLREJECTS,pst,pstColRejects,
				  sizeof(struct inColRejects),nThreads*sizeof(OREJ));
	mdlAddService(mdl,PST_SWAPREJECTS,pst,pstSwapRejects,
				  nThreads*sizeof(int),nThreads*sizeof(OREJ));
	mdlAddService(mdl,PST_DOMAINCOLOR,pst,pstDomainColor,
				  0,0);
	mdlAddService(mdl,PST_COLORDREJECTS,pst,pstColOrdRejects,
				  sizeof(struct inColOrdRejects),nThreads*sizeof(OREJ));
	mdlAddService(mdl,PST_DOMAINORDER,pst,pstDomainOrder,
				  0,0);
	mdlAddService(mdl,PST_LOCALORDER,pst,pstLocalOrder,
				  0,0);
	mdlAddService(mdl,PST_OUTARRAY,pst,pstOutArray,
				  sizeof(struct inOutArray),0);
	mdlAddService(mdl,PST_OUTVECTOR,pst,pstOutVector,
				  sizeof(struct inOutVector),0);
	mdlAddService(mdl,PST_WRITETIPSY,pst,pstWriteTipsy,
				  sizeof(struct inWriteTipsy),0);
	mdlAddService(mdl,PST_BUILDTREE,pst,pstBuildTree,
				  sizeof(struct inBuildTree),0);
	mdlAddService(mdl,PST_DENSITY,pst,pstDensity,
				  sizeof(struct inDensity),0);
	mdlAddService(mdl,PST_GRAVITY,pst,pstGravity,
				  sizeof(struct inGravity),sizeof(struct outGravity));
	mdlAddService(mdl,PST_CALCE,pst,pstCalcE,
				  0,sizeof(struct outCalcE));
	mdlAddService(mdl,PST_DRIFT,pst,pstDrift,
				  sizeof(struct inDrift),0);
	mdlAddService(mdl,PST_KICK,pst,pstKick,
				  sizeof(struct inKick),0);
	mdlAddService(mdl,PST_READCHECK,pst,pstReadCheck,
				  sizeof(struct inReadCheck),0);
	mdlAddService(mdl,PST_WRITECHECK,pst,pstWriteCheck,
				  sizeof(struct inWriteCheck),0);
	mdlAddService(mdl,PST_SETSOFT,pst,pstSetSoft,
				  sizeof(struct inSetSoft),0);
	mdlAddService(mdl,PST_SETTOTAL,pst,pstSetTotal,
				  0,sizeof(struct outSetTotal));

	mdlHandler(mdl);

	pstFinish(pst);
	}


void main(int argc,char **argv)
{
	MDL mdl;
	MSR msr;
	int iStep,iSec,iRed, i;
	double dTime,E,T,U,dWMax,dIMax,dEMax;
	char achFile[256];
	FILE *fpLog;

	mdlInitialize(&mdl,argv,main_ch);
	msrInitialize(&msr,mdl,argc,argv,"pkdgrav");
	iRed = 0;	
	/*
	 ** Check if a restart has been requested.
	 ** Or if it might be required.
	 */
	if (msrRestart(msr)) {
		dTime = msrReadCheck(msr,&iStep);
		printf("Restart Time:%g Redshift:%g Step:%d\n",dTime,msrRedshift(msr),
			   iStep);
		sprintf(achFile,"%s.log",msrOutName(msr));
		fpLog = fopen(achFile,"a");
		assert(fpLog != NULL);
		while (msrRedshift(msr) <= msrRedOut(msr,iRed)) ++iRed;
		goto Restart;
		}
	sprintf(achFile,"%s.log",msrOutName(msr));
	fpLog = fopen(achFile,"w");
	assert(fpLog != NULL);
	/*
	 ** Include a comment at the start of the log file showing the
	 ** command line options.
	 */
	fprintf(fpLog,"#");
	for (i=0;i<argc;++i) fprintf(fpLog,"%s ",argv[i]);
	fprintf(fpLog,"\n");
	fflush(fpLog);

	msrLogParams(msr, fpLog);

	dTime = msrReadTipsy(msr);
	if (msrComove(msr)) msrStepCosmo(msr,dTime);
	if(prmArgSpecified(msr->prm,"dSoft")) msrSetSoft(msr,msrSoft(msr));

	msrBuildTree(msr);
	msrGravity(msr,&iSec,&dWMax,&dIMax,&dEMax);
	msrCalcE(msr,MSR_INIT_ECOSMO,dTime,&E,&T,&U);

	fprintf(fpLog,"%10g %10g %10g %10g %10g %5d %7.1f %7.1f %7.1f\n",
			dTime,msrRedshift(msr),E,T,U,iSec,dWMax,dIMax,dEMax);
	fflush(fpLog);

	if (msrSteps(msr) > 0) {
		for (iStep=1;iStep<=msrSteps(msr);++iStep) {
			msrDrift(msr,msrDelta(msr)/2.0);
			dTime += msrDelta(msr)/2.0;
			msrStepCosmo(msr,dTime);
			msrBuildTree(msr);
			msrGravity(msr,&iSec,&dWMax,&dIMax,&dEMax);
			msrKick(msr,msrDelta(msr));
			msrDrift(msr,msrDelta(msr)/2.0);
			dTime += msrDelta(msr)/2.0;
			msrStepCosmo(msr,dTime);
			if (iStep%msrLogInterval(msr) == 0) {
				/*
				 ** Output a log file line.
				 */
				msrBuildTree(msr);
				msrGravity(msr,&iSec,&dWMax,&dIMax,&dEMax);
				msrCalcE(msr,MSR_STEP_ECOSMO,dTime,&E,&T,&U);

				fprintf(fpLog,"%10g %10g %10g %10g %10g %5d %7.1f %7.1f %7.1f\n",
						dTime,msrRedshift(msr),E,T,U,iSec,dWMax,dIMax,dEMax);
				fflush(fpLog);
				
				}
			if (msrRedshift(msr) <= msrRedOut(msr,iRed)) {
				msrReorder(msr);
				sprintf(achFile,"%s.%05d",msrOutName(msr),iStep);
				msrWriteTipsy(msr,achFile,dTime);
				++iRed;
				}
			if (iStep%msrCheckInterval(msr) == 0) {
				/*
				 ** Write a checkpoint.
				 */
				msrWriteCheck(msr,dTime,iStep);
			Restart:
				;
				}
			}
		/*
		 ** Final output.
		 */
		msrReorder(msr);
		sprintf(achFile,"%s.%05d",msrOutName(msr),iStep);
		msrWriteTipsy(msr,achFile,dTime);
		} 
	else {
		msrReorder(msr);
		sprintf(achFile,"%s.acc",msrOutName(msr));
		msrOutVector(msr,achFile,OUT_ACCEL_VECTOR);
		sprintf(achFile,"%s.pot",msrOutName(msr));
		msrOutArray(msr,achFile,OUT_POT_ARRAY);
		}
	msrFinish(msr);
	mdlFinish(mdl);
	}

