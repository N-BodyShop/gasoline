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

	lcl.pszDataPath = (char *)getenv("PTOOLS_DATA_PATH");
	pstInitialize(&pst,mdl,&lcl);

	mdlAddService(mdl,PST_SETADD,(void *)pst,pstSetAdd);
	mdlAddService(mdl,PST_LEVELIZE,(void *)pst,pstLevelize);
	mdlAddService(mdl,PST_SHOWPST,(void *)pst,pstShowPst);
	mdlAddService(mdl,PST_READTIPSY,(void *)pst,pstReadTipsy);
	mdlAddService(mdl,PST_DOMAINDECOMP,(void *)pst,pstDomainDecomp);
	mdlAddService(mdl,PST_CALCBOUND,(void *)pst,pstCalcBound);
	mdlAddService(mdl,PST_WEIGHT,(void *)pst,pstWeight);
	mdlAddService(mdl,PST_FREESTORE,(void *)pst,pstFreeStore);
	mdlAddService(mdl,PST_COLREJECTS,(void *)pst,pstColRejects);
	mdlAddService(mdl,PST_SWAPREJECTS,(void *)pst,pstSwapRejects);
	mdlAddService(mdl,PST_DOMAINCOLOR,(void *)pst,pstDomainColor);
	mdlAddService(mdl,PST_COLORDREJECTS,(void *)pst,pstColOrdRejects);
	mdlAddService(mdl,PST_DOMAINORDER,(void *)pst,pstDomainOrder);
	mdlAddService(mdl,PST_LOCALORDER,(void *)pst,pstLocalOrder);
	mdlAddService(mdl,PST_OUTARRAY,(void *)pst,pstOutArray);
	mdlAddService(mdl,PST_OUTVECTOR,(void *)pst,pstOutVector);
	mdlAddService(mdl,PST_WRITETIPSY,(void *)pst,pstWriteTipsy);
	mdlAddService(mdl,PST_BUILDTREE,(void *)pst,pstBuildTree);
	mdlAddService(mdl,PST_DENSITY,(void *)pst,pstDensity);
	mdlAddService(mdl,PST_GRAVITY,(void *)pst,pstGravity);
	mdlAddService(mdl,PST_CALCE,(void *)pst,pstCalcE);
	mdlAddService(mdl,PST_DRIFT,(void *)pst,pstDrift);
	mdlAddService(mdl,PST_KICK,(void *)pst,pstKick);
	mdlAddService(mdl,PST_READCHECK,(void *)pst,pstReadCheck);
	mdlAddService(mdl,PST_WRITECHECK,(void *)pst,pstWriteCheck);
	mdlAddService(mdl,PST_SETSOFT,(void *)pst,pstSetSoft);
	mdlAddService(mdl,PST_SETTOTAL,(void *)pst,pstSetTotal);

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

        if(prmArgSpecified(msr->prm, "dSoft")) msrSetSoft(msr);

	msrBuildTree(msr);
	msrGravity(msr,&iSec,&dWMax,&dIMax,&dEMax);
	msrCalcE(msr,MSR_INIT_ECOSMO,dTime,&E,&T,&U);

	fprintf(fpLog,"%10g %10g %10g %10g %10g %5d %7.1f %7.1f %7.1f\n",
			dTime,msrRedshift(msr),E,T,U,iSec,dWMax,dIMax,dEMax);
	fflush(fpLog);

	if (msrSteps(msr) > 0) {
		for (iStep=1;iStep<=msrSteps(msr);++iStep) {
			msrDrift(msr,msrDelta(msr)/2);
			dTime += msrDelta(msr)/2;
			msrStepCosmo(msr,dTime);
			msrBuildTree(msr);
			msrGravity(msr,&iSec,&dWMax,&dIMax,&dEMax);
			msrKick(msr,msrDelta(msr));
			msrDrift(msr,msrDelta(msr)/2);
			dTime += msrDelta(msr)/2;
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
#if 0
		sprintf(achFile,"%s.imass",msrOutName(msr));
		msrOutArray(msr,achFile,OUT_IMASS_ARRAY);
#endif
		}
	msrFinish(msr);
	mdlFinish(mdl);
	}

