#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "mdl.h"
#include "pst.h"
#include "master.h"
#include "outtype.h"


void main_ch(MDL mdl)
{
	PST pst;
	LCL lcl;

	lcl.pszDataPath = (char *)getenv("PTOOLS_DATA_PATH");
	pstInitialize(&pst,mdl,&lcl);

	pstAddServices(pst,mdl);

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
	/*
	 ** Read in the binary file, this may set the number of timesteps or
	 ** the size of the timestep when the zto parameter is used.
	 */
	dTime = msrReadTipsy(msr);
	if(prmArgSpecified(msr->prm,"dSoft")) msrSetSoft(msr,msrSoft(msr));
	/*
	 ** Now we have all the parameters for the simulation we can make a 
	 ** log file entry.
	 */
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

	if (msrComove(msr)) msrStepCosmo(msr,dTime);

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

