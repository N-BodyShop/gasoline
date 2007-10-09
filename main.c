/*#include <fenv.h>*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include "mdl.h"
#include "master.h"
#include "outtype.h"
#include "smoothfcn.h"

#ifdef COLLISIONS
#include <rpc/xdr.h> /* needed for time stamping terse collision log on restart */
#endif

/*DEBUG for FPE trapping... (not defined on most systems)
#include <fpu_control.h>
*/

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

/*DEBUG Should opaque nature of "msr" be enforced in main()? -- DCR*/ 

#ifdef AMPI
#include <converse.h>
#define printf CmiPrintf
/* Charm MPI requires this name as "main" */
int AMPI_Main(int argc,char **argv)
#else
int main(int argc,char **argv)
#endif
{
	MDL mdl;
	MSR msr;
	FILE *fpLog = NULL;
	FILE *fpLogTiming = NULL;
	char achFile[256]; /*DEBUG use MAXPATHLEN here (& elsewhere)? -- DCR*/
	double dTime;
	double E=0,T=0,U=0,Eth=0,L[3]={0,0,0};
	double dWMax=0,dIMax=0,dEMax=0,dMass=0,dMultiEff=0;
	long lSec=0,lStart;
	int i,iStep,iSec=0,iStop=0,nActive,iNumOutputs, OutputList[NUMOUTPUTS];
	char achBaseMask[256];

#ifdef COLLISIONS
	double sec,dsec;
#endif

	/* code to make gasoline core dump if there is a floating point exception 
	feenableexcept(FE_OVERFLOW | FE_DIVBYZERO | FE_INVALID);*/ 

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

#ifndef CCC
	/* no stdout buffering */
	setbuf(stdout,(char *) NULL);
#endif

#ifdef __FAST_MATH__
	assert(0); /* too dangerous! (gcc compile option) */
#endif

/*DEBUG following may work to explicitly trap FPEs -- also uncomment #include above
	{
	fpu_control_t cw = 0x1372;
	_FPU_SETCW(cw);
	}
*/

	lStart=time(0);
	mdlInitialize(&mdl,argv,main_ch);
	for(argc = 0; argv[argc]; argc++); /* some MDLs can trash argv */
	msrInitialize(&msr,mdl,argc,argv);

	(void) strncpy(achBaseMask,msr->param.achDigitMask,256);

	/*
	 Look for checkpoint files.  If not found, we start as normal.
	 If found, msrFindCheck() will move most recent to .chk, and 
	 we restart. bOverwrite means start from beginning, even if 
	 checkpoints exist.
	 */
	if(!msr->param.bOverwrite && msrFindCheck(msr)) {
                msr->param.bRestart = 1;
		dTime = msrReadCheck(msr,&iStep);
		msr->param.bRestart = 1;
#ifdef COLLISIONS
		if (msr->param.nSmooth > msr->N) {
			msr->param.nSmooth = msr->N;
			if (msr->param.bVWarnings)
				printf("WARNING: nSmooth reduced to %i\n",msr->N);
			}
#endif /* COLLISIONS */

#ifdef AGGS
		/*
		 ** Aggregate info not currently stored in checkpoints, so
		 ** reconstruct now.
		 */
		msrAggsFind(msr);
#endif
#ifdef GASOLINE
#ifndef NOCOOLING
		if (msr->param.iGasModel == GASMODEL_COOLING
			|| msr->param.bStarForm) 
		    msrInitCooling(msr);
		if(msr->param.bStarForm)
		    msrInitStarLog(msr);
#endif
#endif
		if(msr->param.bRotatingBar) {
		    msrInitRotatingBar(msr, dTime);
		    }
		msrInitStep(msr);
		dMass = msrMassCheck(msr,-1.0,"Initial");
		msrSetSink(msr);
		if (msr->param.bVStart) printf("Restart Step:%d\n",iStep);
		if (msrLogInterval(msr)) {
			sprintf(achFile,"%s.log",msrOutName(msr));
			fpLog = fopen(achFile,"a");
			assert(fpLog != NULL);
			setbuf(fpLog,(char *) NULL); /* no buffering */
			fprintf(fpLog,"# RESTART (dTime = %g)\n# ",dTime);
			for (i=0;i<argc;++i) fprintf(fpLog,"%s ",argv[i]);
			fprintf(fpLog,"\n");
			msrLogParams(msr,fpLog);
			/* Timing data, if requested */
			if ((fpLogTiming = LogTimingInit( msr, "a" ))) {
			    fprintf(fpLogTiming,"# RESTART (dTime = %g)\n# ",dTime);
			    }
		        }
#ifdef COLLISIONS
		if (msr->param.iCollLogOption != COLL_LOG_NONE) {
			FILE *fp = fopen(msr->param.achCollLog,"r");
			if (fp) { /* add RESTART tag only if file already exists */
				fclose(fp);
				fp = fopen(msr->param.achCollLog,"a");
				assert(fp != NULL);
				switch (msr->param.iCollLogOption) {
				case COLL_LOG_VERBOSE:
					fprintf(fp,"RESTART:T=%e\n",dTime);
					break;
				case COLL_LOG_TERSE:
					{
					XDR xdrs;
					int dum = -1;
					xdrstdio_create(&xdrs,fp,XDR_ENCODE);
					(void) xdr_double(&xdrs,&dTime);
					(void) xdr_int(&xdrs,&dum);
					(void) xdr_int(&xdrs,&dum);
					(void) xdr_int(&xdrs,&dum);
					xdr_destroy(&xdrs);
					break;
					}
				default:
					assert(0); /* should never happen */
					}
				fclose(fp);
				}
			}
#endif
		if(msrKDK(msr) || msr->param.bGravStep || msr->param.bAccelStep) {
			msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE);
			msrDomainDecomp(msr,0,1);
			msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE);
			msrInitAccel(msr);

			msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE);
			msrUpdateSoft(msr,dTime);
			msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE);
			msrBuildTree(msr,0,dMass,0);
			msrMassCheck(msr,dMass,"After msrBuildTree");
			if (msrDoGravity(msr)) {
				msrGravity(msr,iStep,msrDoSun(msr),&iSec,&dWMax,&dIMax,&dEMax,&nActive);
				}
			}
#ifdef GASOLINE
		msrInitSph(msr,dTime);
#endif
		if (msr->param.bDoSinksAtStart) msrDoSinks(msr, dTime, 0.0, 0);
		/* 
		 ** Dump Frame Initialization
		 */
		/* Bring frame count up to correct place for restart. */
		if( msrDumpFrameInit( msr, dTime, 1.0*msr->param.iStartStep, 1 )
                    && msr->df->dDumpFrameStep > 0) {
			while(msr->df->dStep + msr->df->dDumpFrameStep < iStep) {
				msr->df->dStep += msr->df->dDumpFrameStep;
				msr->df->nFrame++;
			}
		}

		if (msrSteps(msr) == 0) goto CheckForDiagnosticOutput;
		goto Restart;
		}
	if(msr->param.bRestart) {
	    printf("Error: restart requested and no checkpoint file found\n");
	    msrFinish(msr);
	    mdlFinish(mdl);
	    return 1;
	    }
	    
	/*
	 ** Read in the binary file, this may set the number of timesteps or
	 ** the size of the timestep when the zto parameter is used.
	 */
#ifndef COLLISIONS
	dTime = msrReadTipsy(msr);
#else
	dTime = msrReadSS(msr); /* must use "Solar System" (SS) I/O format... */
	if (msr->param.nSmooth > msr->N) {
		msr->param.nSmooth = msr->N;
		if (msr->param.bVWarnings)
			printf("WARNING: nSmooth reduced to %i\n",msr->N);
		}
	if (msr->param.iCollLogOption != COLL_LOG_NONE) {
		FILE *fp;
		if (msr->param.iStartStep > 0) { /* append if non-zero start step */
			fp = fopen(msr->param.achCollLog,"a");
			assert(fp != NULL);
			fprintf(fp,"START:T=%e\n",dTime);
			}
		else { /* otherwise erase any old log */
			fp = fopen(msr->param.achCollLog,"w");
			assert(fp != NULL);
			}
		fclose(fp);
		}
#ifdef SLIDING_PATCH
	if (msr->param.iRandStep) {
	    FILE *rfp = fopen("random.log","w");
	    assert(rfp);
	    fclose(rfp);
	    msr->param.iNextRandomization=msrGetNextRandomTime(msr->param.iRandStep,msr->param.iStartStep+1);
	    }
	
#endif /* SLIDING_PATCH */

#endif
#ifdef GASOLINE
#ifndef NOCOOLING
	if (msr->param.iGasModel == GASMODEL_COOLING ||
		msr->param.bStarForm)
		msrInitCooling(msr);
	if(msr->param.bStarForm)
	    msrInitStarLog(msr);
#endif
#endif
	if(msr->param.bRotatingBar)
	    msrInitRotatingBar(msr, dTime);
	msrInitStep(msr);
#ifdef GLASS
	msrInitGlass(msr);
#endif
	dMass = msrMassCheck(msr,-1.0,"Initial");
	if (prmSpecified(msr->prm,"dSoft")) msrSetSoft(msr,msrSoft(msr));
	msrMassCheck(msr,dMass,"After msrSetSoft");

	msrSetSink(msr);
#ifdef COLLISIONS
	if (msr->param.bFindRejects) msrFindRejects(msr);
#endif
#ifdef AGGS
	/* find and initialize any aggregates */
	msrAggsFind(msr);
	msrMassCheck(msr,dMass,"After msrAggsFind");
#endif
	/*
	 ** If the simulation is periodic make sure to wrap all particles into
	 ** the "unit" cell. Doing a drift of 0.0 will always take care of this.
	 */
	msrDrift(msr,dTime,0.0); /* also finds initial overlaps for COLLISIONS */
	msrMassCheck(msr,dMass,"After initial msrDrift");

 CheckForDiagnosticOutput:
	if (msrSteps(msr) > 0) {
		if (msrComove(msr)) {
			msrSwitchTheta(msr,dTime);
			}
		/*
		 ** Now we have all the parameters for the simulation we can make a 
		 ** log file entry.
		 */
		if (msrLogInterval(msr)) {
			sprintf(achFile,"%s.log",msrOutName(msr));
			fpLog = fopen(achFile,"w");
			assert(fpLog != NULL);
			setbuf(fpLog,(char *) NULL); /* no buffering */
			/*
			 ** Include a comment at the start of the log file showing the
			 ** command line options.
			 */
			fprintf(fpLog,"# ");
			for (i=0;i<argc;++i) fprintf(fpLog,"%s ",argv[i]);
			fprintf(fpLog,"\n");
			msrLogParams(msr,fpLog);
			/* Timing data, if requested */
			fpLogTiming = LogTimingInit( msr, "w" );
			}
		/*
		 ** Build tree, activating all particles first (just in case).
		 */
		msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE);
		msrDomainDecomp(msr,0,1);
		msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE);
		msrInitAccel(msr);
		msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE);
		msrUpdateSoft(msr,dTime);
		msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE);
		msrBuildTree(msr,0,dMass,0);
		msrMassCheck(msr,dMass,"After msrBuildTree");
		if (msrDoGravity(msr)) {
			msrGravity(msr,0.0,msrDoSun(msr),&iSec,&dWMax,&dIMax,&dEMax,&nActive);
			msrMassCheck(msr,dMass,"After msrGravity");
			msrCalcEandL(msr,MSR_INIT_E,dTime,&E,&T,&U,&Eth,L);
			msrMassCheck(msr,dMass,"After msrCalcEandL");
			dMultiEff = 1.0;
			if (msrLogInterval(msr)) {
				(void) fprintf(fpLog,"%e %e %.16e %e %e %e %.16e %.16e %.16e "
							   "%i %e %e %e %e\n",dTime,
							   1.0/csmTime2Exp(msr->param.csm,dTime)-1.0,
							   E,T,U,Eth,L[0],L[1],L[2],iSec,dWMax,dIMax,dEMax,
							   dMultiEff);
				}
			/* LogTimingOutput( msr, fpLogTiming, dTime, 0 ); */
			}
#ifdef GASOLINE
		msrInitSph(msr,dTime);
#endif
		if (msr->param.bDoSinksAtStart) msrDoSinks(msr, dTime, 0.0, 0);
		/* 
		 ** Dump Frame Initialization
		 */
		msrDumpFrameInit( msr, dTime, 1.0*msr->param.iStartStep, 0);

		LogTimingZeroCounters( msr );
		for (iStep=msr->param.iStartStep+1;iStep<=msr->param.iStopStep;++iStep) {
			if (msrComove(msr)) {
				msrSwitchTheta(msr,dTime);
				}
			if (msrKDK(msr)) {
				dMultiEff = 0.0;
				lSec = time(0);
#ifdef OLD_KEPLER
				if (msr->param.bFandG) {
					msrPlanetsKDK(msr,iStep - 1,dTime,msrDelta(msr),
								  &dWMax,&dIMax,&dEMax,&iSec);
					continue;
				}
#endif
#ifdef COLLISIONS
				if (msr->param.iMinBinaryRung > 0 && 
					msr->iCurrMaxRung >= msr->param.iMinBinaryRung) {
					if (msr->param.bVDetails) {
					  sec = msrTime();
					  printf("\nSearching for binaries...\n");
					  msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
					  msrDomainDecomp(msr,0,1);
					  msrBuildTree(msr,0,dMass,1);
					  msrCheckForBinary(msr,dTime);
					  dsec=msrTime() - sec;
					  printf("Binary search complete, Wallclock: %f sec\n\n",dsec);
					} else {
					  msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
					  msrDomainDecomp(msr,0,1);
					  msrBuildTree(msr,0,dMass,1);

					  msrCheckForBinary(msr,dTime);
					}
				}

#ifdef SLIDING_PATCH
				if (msr->param.iRandStep) {
					if (iStep >= msr->param.iNextRandomization) {
					    msrRandomizeLargeMasses(msr,iStep,dTime);
					}
				}
#endif /* SLIDING_PATCH */

#endif /* COLLISIONS */

#ifdef RUBBLE_ZML
				{
				  int j;
				  msrMassCheck(msr,dMass,"Before msrRubCleanup"); /*DEBUG*/

				  /*
				  ** Are there any dust particles that need to be added 
				  ** to the dust bins?
				  ** Skip step 0 so initial conditions are preserved.
				  */

				  if (iStep > 0)
					msrRubCleanup(msr,dTime);
				  msrMassCheck(msr,dMass,"After msrRubCleanup"); /*DEBUG*/
				  /*
				  ** Is it time to add dust to planetesimals?
				  ** Skip step 0 so initial conditions are preserved.
				  */

				  if (iStep > 0 && msr->param.CP.DB.nDustBins > 0 && 
					  iStep%msr->param.CP.DB.iDustBinsApplyInt == 0) {
					msrDustBinsApply(msr);
					if (iStep%msr->param.iOutInterval == 0) {
					  printf("iStep = %i\n", iStep);
					  for (j=0;j<msr->param.CP.DB.nDustBins;j++)
						printf("DustBin[%i] = %e\n",j,
							   msr->aDustBins[j].dMass);
					}
				  }
				  msrMassCheck(msr,dMass,"After dust applied to planetesimals"); /*DEBUG*/
				  /*
				  ** The rubble routines need to know if two
				  ** planetesimals will collide during the drift
				  ** interval so that they can be forced to the
				  ** smallest rung. But this may actually result in
				  ** the two planetesimals *not* colliding (since
				  ** their orbits will be better integrated), so
				  ** it's necessary before each top step to reset
				  ** the flags warning of imminent collision.
				  */
				  msrRubbleResetColFlag(msr);
				  msrMassCheck(msr,dMass,"Before msrTopStepKDK"); /*DEBUG*/
				}
#endif
				{
				  msrTopStepKDK(msr,iStep-1,dTime,
								msrDelta(msr),0,0,1,
								&dMultiEff,&dWMax,&dIMax,
								&dEMax,&iSec);
				}
				
				msrRungStats(msr);
				msrCoolVelocity(msr,dTime,dMass);	/* Supercooling if specified */
				msrMassCheck(msr,dMass,"After CoolVelocity in KDK");
				dTime += msrDelta(msr);
				lSec = time(0) - lSec;
				/*
				** Output a log file line if requested.
				** Note: no extra gravity calculation required.
				*/
				if (msrLogInterval(msr) && iStep%msrLogInterval(msr) == 0) {
				  msrCalcEandL(msr,MSR_STEP_E,dTime,&E,&T,&U,&Eth,L);
				  msrMassCheck(msr,dMass,"After msrCalcEandL in KDK");
				  (void) fprintf(fpLog,"%e %e %.16e %e %e %e %.16e %.16e "
								 "%.16e %li %e %e %e %e\n",dTime,
								 1.0/csmTime2Exp(msr->param.csm,dTime)-1.0,
								 E,T,U,Eth,L[0],L[1],L[2],lSec,dWMax,dIMax,
								 dEMax,dMultiEff);
				}
				LogTimingOutput( msr, fpLogTiming, dTime, 0 );
			}
			else {
				lSec = time(0);
				msr->bDoneDomainDecomp = 0;
				msrTopStepDKD(msr,iStep-1,dTime,msrDelta(msr),&dMultiEff);
				msrRungStats(msr);
				msrCoolVelocity(msr,dTime,dMass); /* Supercooling if specified */
				msrMassCheck(msr,dMass,"After CoolVelocity in DKD");
				msrGrowMass(msr,dTime,msrDelta(msr)); /* Grow Masses if specified */
				dTime += msrDelta(msr);
				if (msrLogInterval(msr) && iStep%msrLogInterval(msr) == 0) {
					/*
					 ** Output a log file line.
					 ** Reactivate all particles.
					 */
					if (msrDoGravity(msr)) {
						msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE);
						msrDomainDecomp(msr,0,1);
						msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE);
						msrUpdateSoft(msr,dTime);
						msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE);
						msrBuildTree(msr,0,dMass,0);
						msrMassCheck(msr,dMass,"After msrBuildTree in DKD-log");
						msrInitAccel(msr);
						msrGravity(msr,iStep,msrDoSun(msr),&iSec,&dWMax,&dIMax,&dEMax,&nActive);
						msrMassCheck(msr,dMass,"After msrGravity in DKD-log");
						}
					msrCalcEandL(msr,MSR_STEP_E,dTime,&E,&T,&U,&Eth,L);
					msrMassCheck(msr,dMass,"After msrCalcEandL in DKD-log");
					(void) fprintf(fpLog,"%e %e %.16e %e %e %e %.16e %.16e "
								   "%.16e %li %e %e %e %e\n",dTime,
								   1.0/csmTime2Exp(msr->param.csm,dTime)-1.0,
								   E,T,U,Eth,L[0],L[1],L[2],time(0)-lSec,dWMax,
								   dIMax,dEMax,dMultiEff);

					}
				LogTimingOutput( msr, fpLogTiming, dTime, 0 );
				lSec = time(0) - lSec;
				}
			/*
			 ** Check for user interrupt.
			 */
			iStop = msrCheckForStop(msr);
			/*
			** Output if 1) we've hit an output time
			**           2) We are stopping
			**           3) we're at an output interval
			*/
#ifndef BENCHMARK
			if (msrOutTime(msr,dTime) || iStep == msr->param.iStopStep || iStop ||
				(msrOutInterval(msr) > 0 && iStep%msrOutInterval(msr) == 0)) {
				if (msr->nGas && !msr->param.bKDK) {
					msrActiveType(msr,TYPE_GAS,TYPE_ACTIVE|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
					msrBuildTree(msr,1,-1.0,1);
					msrSmooth(msr,dTime,SMX_DENSITY,1);
					}
#ifdef DENSITYU
				msrGetDensityU(msr);
#endif
				msrReorder(msr);
				msrMassCheck(msr,dMass,"After msrReorder in OutTime");
				sprintf(achFile,msr->param.achDigitMask,msrOutName(msr),iStep);
                                msrCreateGasOutputList(msr, &iNumOutputs, OutputList);
                                msrWriteOutputs(msr, achFile, OutputList, iNumOutputs, dTime);

				if (msrDoDensity(msr) || msr->param.bDohOutput) {
					msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
					msrDomainDecomp(msr,0,1);
					msrBuildTree(msr,0,dMass,1);
					msrSmooth(msr,dTime,SMX_DENSITY,1);
				        msrReorder(msr);
				        }
                                msrCreateAllOutputList(msr, &iNumOutputs, OutputList);
                                msrWriteOutputs(msr, achFile, OutputList, iNumOutputs, dTime);
				msrFlushStarLog(msr);
				/*
				 ** Don't allow duplicate outputs.
				 */
				while (msrOutTime(msr,dTime));
				}
#endif
			if (!iStop && msr->param.iWallRunTime > 0) {
			    if (msr->param.iWallRunTime*60 - (time(0)-lStart) < ((int) (lSec*1.5)) ) {
					printf("RunTime limit exceeded.  Writing checkpoint and exiting.\n");
					printf("    iWallRunTime(sec): %d   Time running: %ld   Last step: %ld\n",
						   msr->param.iWallRunTime*60,time(0)-lStart,lSec);
					iStop = 1;
					}
				}
			if (iStop || iStep == msr->param.iStopStep ||
				(msrCheckInterval(msr) && iStep%msrCheckInterval(msr) == 0)) {
				/*
				 ** Write a checkpoint.
				 */
#ifndef BENCHMARK
				msrFlushStarLog(msr);
				msrWriteCheck(msr,dTime,iStep);
				msrMassCheck(msr,dMass,"After msrWriteCheck");
#endif
			Restart:
				;
				}
			if (iStop) break;
			}
		if (msrLogInterval(msr)) {
		    (void) fclose(fpLog);
		    LogTimingFinish( msr, fpLogTiming, dTime );
		    }
		if (msr->param.bVStart) printf("Integration complete\n");
		}
	else {
            /* Do DiagnosticOutput */
            struct inInitDt in;
            msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);

            in.dDelta = 1e37; /* large number */
            pstInitDt(msr->pst,&in,sizeof(in),NULL,NULL);
            msrInitAccel(msr);
            puts("Initialized Accel and dt\n");
            sprintf(achFile,"%s",msrOutName(msr));

            if (msrRestart(msr)) {
                msrReorder(msr);
                sprintf(achFile,"%s",msrOutName(msr));
#ifndef COLLISIONS
                msrWriteTipsy(msr,achFile,dTime);
#else
                msrWriteSS(msr,achFile,dTime);
#endif
                }

            if (msrDoGravity(msr)) {
                msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE );
                msrDomainDecomp(msr,0,1);
                msrUpdateSoft(msr,dTime);
                msrBuildTree(msr,0,dMass,0);
                msrMassCheck(msr,dMass,"After msrBuildTree in OutSingle Gravity");
                msrGravity(msr,0.0,msrDoSun(msr),&iSec,&dWMax,&dIMax,&dEMax,&nActive);
                msrMassCheck(msr,dMass,"After msrGravity in OutSingle Gravity");
                msrReorder(msr);
                msrMassCheck(msr,dMass,"After msrReorder in OutSingle Gravity");
                iNumOutputs = 0;
                OutputList[iNumOutputs++]=OUT_ACCELG_VECTOR;
                OutputList[iNumOutputs++]=OUT_POT_ARRAY;
                msrWriteOutputs(msr, achFile, OutputList, iNumOutputs, dTime);
                msrMassCheck(msr,dMass,"After msrOutArray in OutSingle Gravity");
                }
#ifdef GASOLINE
            if (msr->nGas > 0) {
		msrInitSph(msr,dTime);
                msrCreateGasStepZeroOutputList(msr, &iNumOutputs,OutputList);
#ifdef DENSITYU
		OutputList[(iNumOutputs)++]=OUT_DENSITYU_ARRAY;
#endif
#ifdef STARFORM
#ifdef CHECKSF
		if(msr->param.bStarForm){
		    msrFormStars(msr, dTime, 1e-6); /* dDelta = 1e-6 not 1e37 */
		    OutputList[(iNumOutputs)++]=OUT_TOFF_YR_ARRAY;
		    OutputList[(iNumOutputs)++]=OUT_TCOOL_YR_ARRAY;
		    OutputList[(iNumOutputs)++]=OUT_TDYN_YR_ARRAY;
		    OutputList[(iNumOutputs)++]=OUT_RATIOSOUNDDYN_ARRAY;
		    OutputList[(iNumOutputs)++]=OUT_L_JEANS_ARRAY;
		    OutputList[(iNumOutputs)++]=OUT_ISMALL_JEANS_ARRAY;
		    }
#endif		
#endif
#ifndef NOCOOLING
		if (msr->param.bGasCooling) {
		    OutputList[iNumOutputs++]=OUT_COOL_EDOT_ARRAY;
		    OutputList[iNumOutputs++]=OUT_COOL_COOLING_ARRAY;
		    OutputList[iNumOutputs++]=OUT_COOL_HEATING_ARRAY;
		    }
#endif
                msrReorder(msr);
                if (msr->param.bSphStep) {
                    fprintf(stdout,"Adding SphStep dt\n");
                    msrSphStep(msr,dTime);
                    }
                msrWriteOutputs(msr, achFile, OutputList, iNumOutputs, dTime);
		msrFlushStarLog(msr);
                }
#endif
            /*
             ** Build tree, activating all particles first (just in case).
             */
            if (msrDoDensity(msr) || msr->param.bDensityStep) {
                    msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
                    msrDomainDecomp(msr,0,1);
                    msrBuildTree(msr,0,-1.0,1);
                    msrMassCheck(msr,dMass,"After msrBuildTree in OutSingle Density");
                    msrSmooth(msr,dTime,SMX_DENSITY,1);
                    msrMassCheck(msr,dMass,"After msrSmooth in OutSingle Density");
                    } 
            if (msrDoDensity(msr)) {
                    msrReorder(msr);
                    msrMassCheck(msr,dMass,"After msrReorder in OutSingle Density");
                    iNumOutputs = 0;
                    OutputList[iNumOutputs++]=OUT_DENSITY_ARRAY;
                    msrWriteOutputs(msr, achFile, OutputList, iNumOutputs, dTime);
            /*	sprintf(achFile,"%s.den",msrOutName(msr));
                    msrReorder(msr);
                    msrOutArray(msr,achFile,OUT_DENSITY_ARRAY);*/
                    msrMassCheck(msr,dMass,"After msrOutArray in OutSingle Density");

#if OLD_KEPLER/*DEBUG*/
                    {
                    struct inSmooth smooth;
                    int sec,dsec;

                    sec = time(0);
                    (void) printf("Building QQ tree...\n");
                    msrBuildQQTree(msr, 0, dMass);
                    dsec = time(0) - sec;
                    (void) printf("QQ tree built, time = %i sec\n",dsec);
                    smooth.nSmooth = 32;
                    smooth.bPeriodic = 0;
                    smooth.bSymmetric = 0;
                    smooth.iSmoothType = SMX_ENCOUNTER;
                    smooth.smf.pkd = NULL;
                    smooth.smf.dTime = dTime;
                    smooth.smf.dDelta = msrDelta(msr);
                    smooth.smf.dCentMass = msr->param.dCentMass;
                    sec = time(0);
                    (void) printf("Beginning encounter search...\n");
                    pstQQSmooth(msr->pst,&smooth,sizeof(smooth),NULL,NULL);
                    dsec = time(0) - sec;
                    (void) printf("Encounter search completed, time = %i sec\n",dsec);
                    msrReorder(msr);
                    iNumOutputs = 0;
                    OutputList[iNumOutputs++]=OUT_DENSITY_ARRAY;
                    msrWriteOutputs(msr, achFile, OutputList, iNumOutputs, dTime);
    /*		sprintf(achFile,"%s.enc",msrOutName(msr));

                    msrOutArray(msr,achFile,OUT_DENSITY_ARRAY);*/
                    }
#endif
                    }

            if (msrDoGravity(msr)) {
                    if (msr->param.bGravStep) {
                            fprintf(stdout,"Adding GravStep dt\n");
                            msrGravStep(msr,dTime);
                            }
                    if (msr->param.bAccelStep) {
                            fprintf(stdout,"Adding AccelStep dt\n");
                            msrAccelStep(msr,dTime);
                            }
                    }

            if (msr->param.bDensityStep) {
                fprintf(stdout,"Adding DensStep dt\n");
                msrDensityStep(msr,dTime);
                    }

            if (msr->param.bDeltaAccelStep) {
                fprintf(stdout,"Adding DeltaAccelStep dt\n");
                    if (!msr->param.bDeltaAccelStepGasTree) {
                        msrActiveType(msr,TYPE_ALL,TYPE_TREEACTIVE);
                        msrBuildTree(msr,0,-1.0,1);
                        }
                    else {
                        msrActiveType(msr,TYPE_GAS,TYPE_TREEACTIVE);
                        msrBuildTree(msr,0,-1.0,1);
                    }
                    msrActiveType(msr,TYPE_ALL,TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
                    msrBuildTree(msr,0,-1.0,1);
                    msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
                    /* This smooth sets dt directly -- hardwired coefficient */
                    msrSmooth(msr,dTime,SMX_DELTAACCEL,0);
                }
            msrReorder(msr);
            iNumOutputs = 0;
            OutputList[iNumOutputs++]=OUT_DT_ARRAY;
            if(msr->param.iMaxRung > 1
	       && (msr->param.bDensityStep || msrDoGravity(msr))) {
		msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
                msrDtToRung(msr,0,msrDelta(msr),1);
                msrRungStats(msr);
		OutputList[iNumOutputs++]=OUT_RUNG_ARRAY;
                }
            msrWriteOutputs(msr, achFile, OutputList, iNumOutputs, dTime);
            }
	
	dfFinalize( msr->df );
	msrFinish(msr);
	mdlFinish(mdl);
	return 0;
	}
