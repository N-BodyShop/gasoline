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
	int iStep,iSec, i;
	double dTime,E,T,U,dWMax,dIMax,dEMax;
	char achFile[256];
	FILE *fpLog;

	mdlInitialize(&mdl,argv,main_ch);
	msrInitialize(&msr,mdl,argc,argv,"pkdgrav");
	/*
	 ** Check if a restart has been requested.
	 ** Or if it might be required.
	 */
	if (msrRestart(msr)) {
		dTime = msrReadCheck(msr,&iStep);
		printf("Restart Step:%d\n",iStep);
		sprintf(achFile,"%s.log",msrOutName(msr));
		fpLog = fopen(achFile,"a");
		assert(fpLog != NULL);
		goto Restart;
		}
	/*
	 ** Read in the binary file, this may set the number of timesteps or
	 ** the size of the timestep when the zto parameter is used.
	 */
	dTime = msrReadTipsy(msr);
	if(prmArgSpecified(msr->prm,"dSoft")) msrSetSoft(msr,msrSoft(msr));
	/*
	 ** If the simulation is periodic make sure to wrap all particles into
	 ** the "unit" cell. Doing a drift of 0.0 will always take care of this.
	 */
	msrDrift(msr,dTime,0.0);
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

	msrBuildTree(msr);
	msrGravity(msr,0.0,&iSec,&dWMax,&dIMax,&dEMax);
	msrCalcE(msr,MSR_INIT_ECOSMO,dTime,&E,&T,&U);

	fprintf(fpLog,"%10g %10g %10g %10g %10g %5d %7.1f %7.1f %7.1f\n",
			dTime,1.0/msrTime2Exp(msr,dTime)-1.0,E,T,U,iSec,dWMax,dIMax,dEMax);
	fflush(fpLog);

	if (msrSteps(msr) > 0) {
		for (iStep=1;iStep<=msrSteps(msr);++iStep) {
			if (msrKDK(msr)) {
				msrKick(msr,dTime,0.5*msrDelta(msr));
				msrDrift(msr,dTime,msrDelta(msr));
				msrBuildTree(msr);
				msrGravity(msr,iStep,&iSec,&dWMax,&dIMax,&dEMax);
				msrKick(msr,dTime+0.5*msrDelta(msr),0.5*msrDelta(msr));
				dTime += msrDelta(msr);
				/*
				 ** Output a log file line at each step.
				 ** Note: no extra gravity calculation required.
				 */
				msrCalcE(msr,MSR_STEP_ECOSMO,dTime,&E,&T,&U);
				fprintf(fpLog,"%10g %10g %10g %10g %10g %5d %7.1f %7.1f %7.1f\n",
						dTime,1.0/msrTime2Exp(msr,dTime)-1.0,E,T,U,iSec,dWMax,dIMax,dEMax);
				fflush(fpLog);			
				}
			else {
				msrDrift(msr,dTime,0.5*msrDelta(msr));
				msrBuildTree(msr);
				msrGravity(msr,iStep-0.5,&iSec,&dWMax,&dIMax,&dEMax);
				msrKick(msr,dTime,msrDelta(msr));
				msrDrift(msr,dTime+0.5*msrDelta(msr),0.5*msrDelta(msr));
				dTime += msrDelta(msr);
				if (iStep%msrLogInterval(msr) == 0) {
					/*
					 ** Output a log file line.
					 */
					msrBuildTree(msr);
					msrGravity(msr,iStep,&iSec,&dWMax,&dIMax,&dEMax);
					msrCalcE(msr,MSR_STEP_ECOSMO,dTime,&E,&T,&U);
					fprintf(fpLog,"%10g %10g %10g %10g %10g %5d %7.1f %7.1f %7.1f\n",
							dTime,1.0/msrTime2Exp(msr,dTime)-1.0,E,T,U,iSec,dWMax,dIMax,dEMax);
					fflush(fpLog);		
					}
				}
			if (msrOutTime(msr,dTime)) {
				msrReorder(msr);
				sprintf(achFile,"%s.%05d",msrOutName(msr),iStep);
				msrWriteTipsy(msr,achFile,dTime);
				/*
				 ** Don't allow duplicate outputs.
				 */
				while (msrOutTime(msr,dTime));
				}
			else if (iStep == msrSteps(msr)) {
				/*
				 ** Final output always produced.
				 */
				msrReorder(msr);
				sprintf(achFile,"%s.%05d",msrOutName(msr),iStep);
				msrWriteTipsy(msr,achFile,dTime);
				}
			else if (msrOutInterval(msr) > 0) {
				if (iStep%msrOutInterval(msr) == 0) {
					msrReorder(msr);
					sprintf(achFile,"%s.%05d",msrOutName(msr),iStep);
					msrWriteTipsy(msr,achFile,dTime);
					}
				}
			if (iStep%msrCheckInterval(msr) == 0 && iStep != msrSteps(msr)) {
				/*
				 ** Write a checkpoint.
				 */
				msrWriteCheck(msr,dTime,iStep);
			Restart:
				;
				}
			}
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

