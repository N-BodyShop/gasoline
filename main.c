#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "mdl.h"
#include "pst.h"
#include "master.h"
#include "outtype.h"
#include "smoothfcn.h"

#ifdef FORCE_CORE_DUMPS
#include <signal.h>
#include <sys/resource.h>
#endif

void main_ch(MDL mdl)
{
	PST pst;
	LCL lcl;

	lcl.pszDataPath = (char *)getenv("PTOOLS_DATA_PATH");
	lcl.pkd = NULL;
	pstInitialize(&pst,mdl,&lcl);

	pstAddServices(pst,mdl);

	mdlHandler(mdl);

	pstFinish(pst);
	}


int main(int argc,char **argv)
{
	MDL mdl;
	MSR msr;
	int iStep,iSec, i;
	double dTime,E,T,U,dWMax,dIMax,dEMax;
	char achFile[256];
	FILE *fpLog = NULL;
	double dMass;
	int nActive;
	double dMultiEff;

#ifdef TINY_PTHREAD_STACK
	static int first = 1;
	static char **save_argv;

	/*
	 * Hackery to get around SGI's tiny pthread stack.
	 * Main will be called twice.  The second time, argc and argv
	 * will be garbage, so we have to save them from the first.
	 * Another way to do this would involve changing the interface
	 * to mdlInitialize(), so that this hackery could be hidden
	 * down there.
	 */
	if(first) {
	    save_argv = malloc((argc+1)*sizeof(*save_argv));
	    for(i = 0; i < argc; i++)
	        save_argv[i] = strdup(argv[i]);
	    save_argv[argc] = NULL;
	    }
	else {
	    argv = save_argv;
	    }
	first = 0;
#endif /* TINY_PTHREAD_STACK */

#ifdef FORCE_CORE_DUMPS
	{
	/* Force core dumps on segv */

	struct rlimit rlim;

	rlim.rlim_cur = rlim.rlim_max = RLIM_INFINITY;
	setrlimit(RLIMIT_CORE,&rlim);
	signal(SIGSEGV,SIG_DFL);
	}
#endif

#ifdef NO_STDOUT_BUF
	setbuf(stdout,(char *)NULL);
#endif

	mdlInitialize(&mdl,argv,main_ch);
	for(argc = 0; argv[argc]; argc++); /* some MDLs can trash argv */
	msrInitialize(&msr,mdl,argc,argv);
	/*
	 ** Check if a restart has been requested.
	 ** Or if it might be required.
	 */
	if (msrRestart(msr)) {
		dTime = msrReadCheck(msr,&iStep);
		msrInitStep(msr);
		dMass = msrMassCheck(msr,-1.0,"Initial");
		printf("Restart Step:%d\n",iStep);
		if (msrLogInterval(msr)) {
			sprintf(achFile,"%s.log",msrOutName(msr));
			fpLog = fopen(achFile,"a");
			assert(fpLog != NULL);
			}
		if(msrKDK(msr) || msr->param.bEpsVel) {
			msrBuildTree(msr,0,dMass,0);
			msrMassCheck(msr,dMass,"After msrBuildTree");
			msrInitAccel(msr);
			if (msrDoGravity(msr)) {
				msrGravity(msr,iStep,msrDoSun(msr),&iSec,&dWMax,&dIMax,&dEMax,&nActive);
				}
			}
		goto Restart;
		}
	/*
	 ** Read in the binary file, this may set the number of timesteps or
	 ** the size of the timestep when the zto parameter is used.
	 */
#ifdef PLANETS
	dTime = msrReadSS(msr);
	if (msr->param.nSmooth > msr->N)
		printf("WARNING: nSmooth > N\n");
#else /* PLANETS */
	dTime = msrReadTipsy(msr);
#endif /* !PLANETS */
	msrInitStep(msr);
#ifdef GASOLINE
	msrInitSph(msr,dTime);
#endif
	dMass = msrMassCheck(msr,-1.0,"Initial");
	if (prmSpecified(msr->prm,"dSoft")) msrSetSoft(msr,msrSoft(msr));
	msrMassCheck(msr,dMass,"After msrSetSoft");
	/*
	 ** If the simulation is periodic make sure to wrap all particles into
	 ** the "unit" cell. Doing a drift of 0.0 will always take care of this.
	 */
#ifdef PLANETS
	/* This also checks for initial particle overlaps! */
#endif /* PLANETS */
	msrDrift(msr,dTime,0.0);
	msrMassCheck(msr,dMass,"After initial msrDrift");

	if (msrSteps(msr) > 0) {
		/*
		 ** Now we have all the parameters for the simulation we can make a 
		 ** log file entry.
		 */
		if (msrComove(msr)) {
			msrSwitchTheta(msr,dTime);
			}
		if (msrLogInterval(msr)) {
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
			}
		/*
		 ** Build tree, activating all particles first (just in case).
		 */
		msrActiveRung(msr,0,1);
		msrBuildTree(msr,0,dMass,0);
		msrMassCheck(msr,dMass,"After msrBuildTree");

		msrInitAccel(msr);
		if (msrDoGravity(msr)) {
			msrGravity(msr,0.0,msrDoSun(msr),&iSec,&dWMax,&dIMax,&dEMax,&nActive);
			msrMassCheck(msr,dMass,"After msrGravity");
			msrCalcE(msr,MSR_INIT_ECOSMO,dTime,&E,&T,&U);
			msrMassCheck(msr,dMass,"After msrCalcE");
			dMultiEff = 1.0;
			if (msrLogInterval(msr)) {
				(void) fprintf(fpLog,"%e %e %e %e %e %i %e %e %e %e\n",dTime,
						1.0/msrTime2Exp(msr,dTime)-1.0,E,T,U,iSec,dWMax,dIMax,
						dEMax,dMultiEff);
				(void) fflush(fpLog);
				}
			}

		for (iStep=1;iStep<=msrSteps(msr);++iStep) {
			if (msrComove(msr)) {
				msrSwitchTheta(msr,dTime);
				}
			if (msrKDK(msr)) {
				dMultiEff = 0.0;
				msrTopStepKDK(msr, iStep-1, dTime,
					      msrDelta(msr), 0, 0, 1,
					      &dMultiEff, &dWMax,
					      &dIMax, &dEMax, &iSec);
				msrRungStats(msr);
				msrCoolVelocity(msr,dTime,dMass);	/* Supercooling if specified */
				msrMassCheck(msr,dMass,"After CoolVelocity in KDK");
				dTime += msrDelta(msr);
				/*
				 ** Output a log file line at each step.
				 ** Note: no extra gravity calculation required.
				 */
				msrCalcE(msr,MSR_STEP_ECOSMO,dTime,&E,&T,&U);
				msrMassCheck(msr,dMass,"After msrCalcE in KDK");
				if (msrLogInterval(msr) && iStep%msrLogInterval(msr) == 0) {
					(void) fprintf(fpLog,"%e %e %e %e %e %i %e %e %e %e\n",
								   dTime,1.0/msrTime2Exp(msr,dTime)-1.0,E,T,U,
								   iSec,dWMax,dIMax,dEMax,dMultiEff);
					(void) fflush(fpLog);
					}
				}
			else {
				msrTopStepDKD(msr, iStep-1, dTime,
					      msrDelta(msr), &dMultiEff);
				msrRungStats(msr);
				msrCoolVelocity(msr,dTime,dMass);	/* Supercooling if specified */
				msrMassCheck(msr,dMass,"After CoolVelocity in DKD");
				dTime += msrDelta(msr);
				if (iStep%msrLogInterval(msr) == 0) {
					/*
					 ** Output a log file line.
					 ** Reactivate all particles.
					 */
					msrActiveRung(msr,0,1);
					msrBuildTree(msr,0,dMass,0);
					msrMassCheck(msr,dMass,"After msrBuildTree in DKD-log");
					msrInitAccel(msr);
					if (msrDoGravity(msr)) {
						msrGravity(msr,iStep,msrDoSun(msr),&iSec,&dWMax,&dIMax,&dEMax,&nActive);
						msrMassCheck(msr,dMass,"After msrGravity in DKD-log");
						msrCalcE(msr,MSR_STEP_ECOSMO,dTime,&E,&T,&U);
						msrMassCheck(msr,dMass,"After msrCalcE in DKD-log");
						if (msrLogInterval(msr) && msrLogInterval(msr)%iStep == 0) {
							(void) fprintf(fpLog,"%e %e %e %e %e %i %e %e %e %e\n",dTime,
										   1.0/msrTime2Exp(msr,dTime)-1.0,E,T,U,iSec,
										   dWMax,dIMax,dEMax,dMultiEff);
							(void) fflush(fpLog);
							}
						}
					}
				}
			if (msrOutTime(msr,dTime)) {
				if (msrDoDensity(msr)) {
					msrActiveRung(msr,0,1);
					msrBuildTree(msr,0,dMass,1);
					msrMassCheck(msr,dMass,"After msrBuildTree in OutTime");
					msrSmooth(msr,dTime,SMX_DENSITY,1);
					msrMassCheck(msr,dMass,"After msrSmooth in OutTime");
					}
				msrReorder(msr);
				msrMassCheck(msr,dMass,"After msrReorder in OutTime");
				sprintf(achFile,"%s.%05d",msrOutName(msr),iStep);
#ifdef PLANETS
				msrWriteSS(msr,achFile,dTime);
				msrMassCheck(msr,dMass,"After msrWriteSS in OutTime");
#ifdef OUTPUT_DT
				{
				char achFileDt[256/*MAXPATHLEN*/];
				(void) sprintf(achFileDt,"%sdt.%05d",msrOutName(msr),iStep);
				msrOutArray(msr,achFileDt,OUT_DT_ARRAY);
				}
#endif
#else /* PLANETS */
				msrWriteTipsy(msr,achFile,dTime);
				msrMassCheck(msr,dMass,"After msrWriteTipsy in OutTime");
#endif /* !PLANETS */
				if (msrDoDensity(msr)) {
					sprintf(achFile,"%s.%05d.den",msrOutName(msr),iStep);
					msrOutArray(msr,achFile,OUT_DENSITY_ARRAY);
					msrMassCheck(msr,dMass,"After msrOutArray in OutTime");
					}
				/*
				 ** Don't allow duplicate outputs.
				 */
				while (msrOutTime(msr,dTime));
				}
			else if (iStep == msrSteps(msr)) {
				/*
				 ** Final output always produced.
				 */
				if (msrDoDensity(msr)) {
					msrActiveRung(msr,0,1);
					msrBuildTree(msr,0,dMass,1);
					msrMassCheck(msr,dMass,"After msrBuildTree in OutFinal");
					msrSmooth(msr,dTime,SMX_DENSITY,1);
					msrMassCheck(msr,dMass,"After msrSmooth in OutFinal");
					}
				msrReorder(msr);
				msrMassCheck(msr,dMass,"After msrReorder in OutFinal");
				sprintf(achFile,"%s.%05d",msrOutName(msr),iStep);
#ifdef PLANETS
				msrWriteSS(msr,achFile,dTime);
				msrMassCheck(msr,dMass,"After msrWriteSS in OutFinal");
#ifdef OUTPUT_DT
				{
				char achFileDt[256/*MAXPATHLEN*/];
				(void) sprintf(achFileDt,"%sdt.%05d",msrOutName(msr),iStep);
				msrOutArray(msr,achFileDt,OUT_DT_ARRAY);
				}
#endif
#else /* PLANETS */
				msrWriteTipsy(msr,achFile,dTime);
				msrMassCheck(msr,dMass,"After msrWriteTipsy in OutFinal");
#endif /* !PLANETS */
				if (msrDoDensity(msr)) {
					sprintf(achFile,"%s.%05d.den",msrOutName(msr),iStep);
					msrOutArray(msr,achFile,OUT_DENSITY_ARRAY);
					msrMassCheck(msr,dMass,"After msrOutArray in OutFinal");
					}
				}
			else if (msrOutInterval(msr) > 0) {
				if (iStep%msrOutInterval(msr) == 0) {
					if (msrDoDensity(msr)) {
						msrActiveRung(msr,0,1);
						msrBuildTree(msr,0,dMass,1);
						msrMassCheck(msr,dMass,"After msrBuildTree in OutInt");
						msrSmooth(msr,dTime,SMX_DENSITY,1);
						msrMassCheck(msr,dMass,"After msrSmooth in OutInt");
						}
					msrReorder(msr);
					msrMassCheck(msr,dMass,"After msrReorder in OutInt");
					sprintf(achFile,"%s.%05d",msrOutName(msr),iStep);
#ifdef PLANETS
					msrWriteSS(msr,achFile,dTime);
					msrMassCheck(msr,dMass,"After msrWriteSS in OutInt");
#ifdef OUTPUT_DT
					{
					char achFileDt[256/*MAXPATHLEN*/];
					(void) sprintf(achFileDt,"%sdt.%05d",msrOutName(msr),iStep);
					msrOutArray(msr,achFileDt,OUT_DT_ARRAY);
					}
#endif
#else /* PLANETS */
					msrWriteTipsy(msr,achFile,dTime);
					msrMassCheck(msr,dMass,"After msrWriteTipsy in OutInt");
#endif /* !PLANETS */
					if (msrDoDensity(msr)) {
						sprintf(achFile,"%s.%05d.den",msrOutName(msr),iStep);
						msrOutArray(msr,achFile,OUT_DENSITY_ARRAY);
						msrMassCheck(msr,dMass,"After msrOutArray in OutInt");
						}
					}
				}
			if (msrCheckInterval(msr) && iStep%msrCheckInterval(msr) == 0 &&
				iStep != msrSteps(msr)) {
				/*
				 ** Write a checkpoint.
				 */
				msrWriteCheck(msr,dTime,iStep);
				msrMassCheck(msr,dMass,"After msrWriteCheck");
			Restart:
				;
				}
			}
		if (msrLogInterval(msr)) (void) fclose(fpLog);
		} 
	else {
		/*
		 ** Build tree, activating all particles first (just in case).
		 */
		if (msrDoDensity(msr)) {
			msrActiveRung(msr,0,1);
			msrBuildTree(msr,0,dMass,1);
			msrMassCheck(msr,dMass,"After msrBuildTree in OutSingle Density");
			msrSmooth(msr,dTime,SMX_DENSITY,1);
			msrMassCheck(msr,dMass,"After msrSmooth in OutSingle Density");
			msrReorder(msr);
			msrMassCheck(msr,dMass,"After msrReorder in OutSingle Density");
			sprintf(achFile,"%s.den",msrOutName(msr));
			msrOutArray(msr,achFile,OUT_DENSITY_ARRAY);
			msrMassCheck(msr,dMass,"After msrOutArray in OutSingle Density");
#ifdef PLANETS
			{
			    struct inSmooth smooth;
			    
			    msrBuildQQTree(msr, 0, dMass);
			    smooth.nSmooth = 32;
			    smooth.bPeriodic = 0;
			    smooth.bSymmetric = 0;
			    smooth.iSmoothType = SMX_ENCOUNTER;
			    smooth.smf.pkd = NULL;
			    pstQQSmooth(msr->pst, &smooth,
					sizeof(smooth), NULL, NULL);
			    msrReorder(msr);
			    sprintf(achFile,"%s.enc",msrOutName(msr));
			    msrOutArray(msr,achFile,OUT_DENSITY_ARRAY);
			    }
#endif
			}
		else if (msrDoGravity(msr)) {
			msrActiveRung(msr,0,1);
			msrBuildTree(msr,0,dMass,0);
			msrMassCheck(msr,dMass,"After msrBuildTree in OutSingle Gravity");
			msrInitAccel(msr);
			msrGravity(msr,0.0,msrDoSun(msr),&iSec,&dWMax,&dIMax,&dEMax,&nActive);
			msrMassCheck(msr,dMass,"After msrGravity in OutSingle Gravity");
			msrReorder(msr);
			msrMassCheck(msr,dMass,"After msrReorder in OutSingle Gravity");
			sprintf(achFile,"%s.acc",msrOutName(msr));
			msrOutVector(msr,achFile,OUT_ACCEL_VECTOR);
			msrMassCheck(msr,dMass,"After msrOutVector in OutSingle Gravity");
			sprintf(achFile,"%s.pot",msrOutName(msr));
			msrOutArray(msr,achFile,OUT_POT_ARRAY);
			msrMassCheck(msr,dMass,"After msrOutArray in OutSingle Gravity");
			msrInitDt(msr);
			msrAccelStep(msr, dTime);
			msrDtToRung(msr, 0, msrDelta(msr), 1);
			sprintf(achFile,"%s.dt",msrOutName(msr));
			msrOutArray(msr,achFile,OUT_DT_ARRAY);
			}
		}
	msrFinish(msr);
	mdlFinish(mdl);
	return 0;
	}
