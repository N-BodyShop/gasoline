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
	double dMass;

	mdlInitialize(&mdl,argv,main_ch);
	msrInitialize(&msr,mdl,argc,argv,"pkdgrav");
	/*
	 ** Check if a restart has been requested.
	 ** Or if it might be required.
	 */
	if (msrRestart(msr)) {
		dTime = msrReadCheck(msr,&iStep);
		dMass = msrMassCheck(msr,-1.0,"Initial");
		printf("Restart Step:%d\n",iStep);
		sprintf(achFile,"%s.log",msrOutName(msr));
		fpLog = fopen(achFile,"a");
		assert(fpLog != NULL);
		if(msrKDK(msr)) {
			msrBuildTree(msr,0,dMass);
			msrMassCheck(msr,dMass,"After msrBuildTree");
			msrGravity(msr,iStep,&iSec,&dWMax,&dIMax,&dEMax);
			}
		goto Restart;
		}
	/*
	 ** Read in the binary file, this may set the number of timesteps or
	 ** the size of the timestep when the zto parameter is used.
	 */
	dTime = msrReadTipsy(msr);
	dMass = msrMassCheck(msr,-1.0,"Initial");
	if(prmArgSpecified(msr->prm,"dSoft")) msrSetSoft(msr,msrSoft(msr));
	msrMassCheck(msr,dMass,"After msrSetSoft");
	/*
	 ** If the simulation is periodic make sure to wrap all particles into
	 ** the "unit" cell. Doing a drift of 0.0 will always take care of this.
	 */
	msrDrift(msr,dTime,0.0);
	msrMassCheck(msr,dMass,"After initial msrDrift");
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
	msrMassCheck(msr,dMass,"After msrLogParams");

	msrBuildTree(msr,0,dMass);
	msrMassCheck(msr,dMass,"After msrBuildTree");
	msrGravity(msr,0.0,&iSec,&dWMax,&dIMax,&dEMax);
	msrMassCheck(msr,dMass,"After msrGravity");
	msrCalcE(msr,MSR_INIT_ECOSMO,dTime,&E,&T,&U);
	msrMassCheck(msr,dMass,"After msrCalcE");

	fprintf(fpLog,"%10g %10g %10g %10g %10g %5d %7.1f %7.1f %7.1f\n",
			dTime,1.0/msrTime2Exp(msr,dTime)-1.0,E,T,U,iSec,dWMax,dIMax,dEMax);
	fflush(fpLog);

	if (msrSteps(msr) > 0) {
		for (iStep=1;iStep<=msrSteps(msr);++iStep) {
			if (msrKDK(msr)) {
				msrKick(msr,dTime,0.5*msrDelta(msr));
				msrMassCheck(msr,dMass,"After msrKick-1 in KDK");
				msrDrift(msr,dTime,msrDelta(msr));
				msrMassCheck(msr,dMass,"After msrDrift in KDK");
				msrBuildTree(msr,0,dMass);
				msrMassCheck(msr,dMass,"After msrBuildTree in KDK");
				msrGravity(msr,iStep,&iSec,&dWMax,&dIMax,&dEMax);
				msrMassCheck(msr,dMass,"After msrGravity in KDK");
				msrKick(msr,dTime+0.5*msrDelta(msr),0.5*msrDelta(msr));
				msrMassCheck(msr,dMass,"After msrKick-2 in KDK");
				dTime += msrDelta(msr);
				/*
				 ** Output a log file line at each step.
				 ** Note: no extra gravity calculation required.
				 */
				msrCalcE(msr,MSR_STEP_ECOSMO,dTime,&E,&T,&U);
				msrMassCheck(msr,dMass,"After msrCalcE in KDK");
				fprintf(fpLog,"%10g %10g %10g %10g %10g %5d %7.1f %7.1f %7.1f\n",
						dTime,1.0/msrTime2Exp(msr,dTime)-1.0,E,T,U,iSec,dWMax,dIMax,dEMax);
				fflush(fpLog);			
				}
			else {
				msrDrift(msr,dTime,0.5*msrDelta(msr));
				msrMassCheck(msr,dMass,"After msrDrift-1 in DKD");
				msrBuildTree(msr,0,dMass);
				msrMassCheck(msr,dMass,"After msrBuildTree in DKD");
				msrGravity(msr,iStep-0.5,&iSec,&dWMax,&dIMax,&dEMax);
				msrMassCheck(msr,dMass,"After msrGravity in DKD");
				msrKick(msr,dTime,msrDelta(msr));
				msrMassCheck(msr,dMass,"After msrKick in DKD");
				msrDrift(msr,dTime+0.5*msrDelta(msr),0.5*msrDelta(msr));
				msrMassCheck(msr,dMass,"After msrDrift-2 in DKD");
				dTime += msrDelta(msr);
				if (iStep%msrLogInterval(msr) == 0) {
					/*
					 ** Output a log file line.
					 */
					msrBuildTree(msr,0,dMass);
					msrMassCheck(msr,dMass,"After msrBuildTree in DKD-log");
					msrGravity(msr,iStep,&iSec,&dWMax,&dIMax,&dEMax);
					msrMassCheck(msr,dMass,"After msrGravity in DKD-log");
					msrCalcE(msr,MSR_STEP_ECOSMO,dTime,&E,&T,&U);
					msrMassCheck(msr,dMass,"After msrCalcE in DKD-log");
					fprintf(fpLog,"%10g %10g %10g %10g %10g %5d %7.1f %7.1f %7.1f\n",
							dTime,1.0/msrTime2Exp(msr,dTime)-1.0,E,T,U,iSec,dWMax,dIMax,dEMax);
					fflush(fpLog);		
					}
				}
			if (msrOutTime(msr,dTime)) {
				msrReorder(msr);
				msrMassCheck(msr,dMass,"After msrReorder in OutTime");
				sprintf(achFile,"%s.%05d",msrOutName(msr),iStep);
				msrWriteTipsy(msr,achFile,dTime);
				msrMassCheck(msr,dMass,"After msrWriteTipsy in OutTime");
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
				msrMassCheck(msr,dMass,"After msrReorder in OutFinal");
				sprintf(achFile,"%s.%05d",msrOutName(msr),iStep);
				msrWriteTipsy(msr,achFile,dTime);
				msrMassCheck(msr,dMass,"After msrWriteTipsy in OutFinal");
				}
			else if (msrOutInterval(msr) > 0) {
				if (iStep%msrOutInterval(msr) == 0) {
					msrReorder(msr);
					msrMassCheck(msr,dMass,"After msrReorder in OutInt");
					sprintf(achFile,"%s.%05d",msrOutName(msr),iStep);
					msrWriteTipsy(msr,achFile,dTime);
					msrMassCheck(msr,dMass,"After msrWriteTipsy in OutInt");
					}
				}
			if (iStep%msrCheckInterval(msr) == 0 && iStep != msrSteps(msr)) {
				/*
				 ** Write a checkpoint.
				 */
				msrWriteCheck(msr,dTime,iStep);
				msrMassCheck(msr,dMass,"After msrWriteCheck");
			Restart:
				;
				}
			}
		} 
	else {
		msrReorder(msr);
		msrMassCheck(msr,dMass,"After msrReorder in OutArray");
		sprintf(achFile,"%s.acc",msrOutName(msr));
		msrOutVector(msr,achFile,OUT_ACCEL_VECTOR);
		msrMassCheck(msr,dMass,"After msrOutVector");
		sprintf(achFile,"%s.pot",msrOutName(msr));
		msrOutArray(msr,achFile,OUT_POT_ARRAY);
		msrMassCheck(msr,dMass,"After msrOutArray");

		sprintf(achFile,"%s.buc",msrOutName(msr));
		msrOutArray(msr,achFile,OUT_COLOR_ARRAY);
		}
	msrFinish(msr);
	mdlFinish(mdl);
	}

