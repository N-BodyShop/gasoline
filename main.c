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

/*DEBUG! for FPE trapping... (not defined on most systems)
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
	int i,iStep,iSec=0,iStop=0,nActive;

	char achBaseMask[256];

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
#endif
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
#endif
#endif
		msrInitStep(msr);
		dMass = msrMassCheck(msr,-1.0,"Initial");
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
		if (msr->param.bDoSinksAtStart) msrDoSinks(msr,0.0);
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
#endif
#ifdef GASOLINE
#ifndef NOCOOLING
	if (msr->param.iGasModel == GASMODEL_COOLING ||
		msr->param.bStarForm)
		msrInitCooling(msr);
#endif
#endif
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
		if (msr->param.bDoSinksAtStart) msrDoSinks(msr,0.0);
		/* 
		 ** Dump Frame Initialization
		 */
		msrDumpFrameInit( msr, dTime, 1.0*msr->param.iStartStep, 0);

		LogTimingZeroCounters( msr );
		for (iStep=msr->param.iStartStep+1;iStep<=msrSteps(msr);++iStep) {
			if (msrComove(msr)) {
				msrSwitchTheta(msr,dTime);
				}
			if (msrKDK(msr)) {
				dMultiEff = 0.0;
				lSec = time(0);
#ifdef OLD_KEPLER
				if (msr->param.bFandG)
					msrPlanetsKDK(msr,iStep - 1,dTime,msrDelta(msr),
								  &dWMax,&dIMax,&dEMax,&iSec);
				else
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
				/*
				 ** Output a log file line if requested.
				 ** Note: no extra gravity calculation required.
				 */
				if (msrLogInterval(msr) && iStep%msrLogInterval(msr) == 0) {
					msrCalcEandL(msr,MSR_STEP_E,dTime,&E,&T,&U,&Eth,L);
					msrMassCheck(msr,dMass,"After msrCalcEandL in KDK");
					lSec = time(0) - lSec;
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
			if (msrOutTime(msr,dTime) || iStep == msrSteps(msr) || iStop ||
				(msrOutInterval(msr) > 0 && iStep%msrOutInterval(msr) == 0)) {
				if (msr->nGas && !msr->param.bKDK) {
					msrActiveType(msr,TYPE_GAS,TYPE_ACTIVE|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
					msrBuildTree(msr,1,-1.0,1);
					msrSmooth(msr,dTime,SMX_DENSITY,1);
					}
				msrReorder(msr);
				msrMassCheck(msr,dMass,"After msrReorder in OutTime");
				sprintf(achFile,msr->param.achDigitMask,msrOutName(msr),iStep);
#ifndef COLLISIONS
				msrWriteTipsy(msr,achFile,dTime);
#else
				msrWriteSS(msr,achFile,dTime);
#endif
				if(msr->param.bDoIOrderOutput) {
				    sprintf(achFile,achBaseMask,
					    msrOutName(msr),iStep);
				    strncat(achFile,".iord",256);
				    msrOutArray(msr,achFile,OUT_IORDER_ARRAY);
				    }
				if (msrDoDensity(msr) || msr->param.bDohOutput) {
					msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
					msrDomainDecomp(msr,0,1);
					msrBuildTree(msr,0,dMass,1);
					msrSmooth(msr,dTime,SMX_DENSITY,1);
				        }
				if (msrDoDensity(msr)) {
				        msrReorder(msr);
					sprintf(achFile,achBaseMask,msrOutName(msr),iStep);
					strncat(achFile,".den",256);
					msrOutArray(msr,achFile,OUT_DENSITY_ARRAY);
				        }
				if (msr->param.bDoSoftOutput) {
				        msrReorder(msr);
					sprintf(achFile,achBaseMask,msrOutName(msr),iStep);
					strncat(achFile,".soft",256);
					msrOutArray(msr,achFile,OUT_SOFT_ARRAY);
				        }
				if (msr->param.bDohOutput) {
					msrReorder(msr);
					sprintf(achFile,achBaseMask,msrOutName(msr),iStep);
					strncat(achFile,".H",256);
					msrOutArray(msr,achFile,OUT_H_ARRAY);
					}
#ifdef GASOLINE				
 				if (msr->param.bDoSphhOutput) {
					msrReorder(msr);
 					sprintf(achFile,achBaseMask,msrOutName(msr),iStep);
					strncat(achFile,".SPHH",256);
					msrOutArray(msr,achFile,OUT_H_ARRAY);
					}
				if (msrDoDensity(msr) || msr->param.bDohOutput || msr->param.bDoSphhOutput) {
					msrActiveType(msr,TYPE_GAS,TYPE_ACTIVE|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
					msrDomainDecomp(msr,0,1);
					msrBuildTree(msr,1,-1.0,1);
					msrSmooth(msr,dTime,SMX_DENSITY,1);
					}
#ifdef PDVDEBUG
				msrReorder(msr);
				sprintf(achFile,achBaseMask,msrOutName(msr),iStep);
				strncat(achFile,".PdVpres",256);
 				msrOutArray(msr,achFile,OUT_PDVPRES_ARRAY);
				sprintf(achFile,achBaseMask,msrOutName(msr),iStep);
				strncat(achFile,".PdVvisc",256);
 				msrOutArray(msr,achFile,OUT_PDVVISC_ARRAY);
#endif
				if (msr->param.bShockTracker) {
 					msrReorder(msr);
					sprintf(achFile,achBaseMask,msrOutName(msr),iStep);
					strncat(achFile,".ST",256);
					msrOutArray(msr,achFile,OUT_SHOCKTRACKER_ARRAY);

					sprintf(achFile,achBaseMask,msrOutName(msr),iStep);
					strncat(achFile,".BSw",256);
					msrOutArray(msr,achFile,OUT_BALSARASWITCH_ARRAY);

					sprintf(achFile,achBaseMask,msrOutName(msr),iStep);
					strncat(achFile,".SPHH",256);
					msrOutArray(msr,achFile,OUT_H_ARRAY);

					sprintf(achFile,achBaseMask,msrOutName(msr),iStep);
					strncat(achFile,".divv",256);
					msrOutArray(msr,achFile,OUT_DIVV_ARRAY);
					
					sprintf(achFile,achBaseMask,msrOutName(msr),iStep);
					strncat(achFile,".divrhov",256);
					msrOutArray(msr,achFile,OUT_DIVRHOV_ARRAY);

					sprintf(achFile,achBaseMask,msrOutName(msr),iStep);
					strncat(achFile,".gradrho",256);
					msrOutVector(msr,achFile,OUT_GRADRHO_VECTOR);

				}
				if (msr->param.bSN) {
					msrReorder(msr);
					sprintf(achFile,achBaseMask,msrOutName(msr),iStep);
					strncat(achFile,".PdVSN",256);
					msrOutArray(msr,achFile,OUT_PDVSN_ARRAY);
					sprintf(achFile,achBaseMask,msrOutName(msr),iStep);
					strncat(achFile,".uSN",256);
					msrOutArray(msr,achFile,OUT_USN_ARRAY);
				        }
#ifndef NOCOOLING				
				{
				int ArrayCnt = 0;
				char OutSuffix[20];
				int OutType;
				
				msrReorder(msr);
				for (;;) {	
					CoolOutputArray( &msr->param.CoolParam, ArrayCnt, &OutType, OutSuffix );
					if (OutType == OUT_NULL) break;
					sprintf(achFile,achBaseMask,msrOutName(msr),iStep);
					strncat(achFile,OutSuffix,256);
					msrOutArray(msr,achFile, OutType);
                    ArrayCnt++;
					}
				}
#endif

				if(msr->param.bStarForm || msr->param.bFeedBack) {
					msrReorder(msr);
					sprintf(achFile,achBaseMask,msrOutName(msr),iStep);
					strncat(achFile,".igasorder",256);
					msrOutArray(msr,achFile,OUT_IGASORDER_ARRAY);
					sprintf(achFile,achBaseMask,msrOutName(msr),iStep);
					strncat(achFile,".massform",256);
					msrOutArray(msr,achFile,OUT_MASSFORM_ARRAY);
					sprintf(achFile,achBaseMask,msrOutName(msr),iStep);
					strncat(achFile,".coolontime",256);
					msrOutArray(msr,achFile,OUT_COOLTURNONTIME_ARRAY);
					sprintf(achFile,achBaseMask,msrOutName(msr),iStep);
					strncat(achFile,".OxMassFrac",256);
					msrOutArray(msr,achFile,OUT_OXYGENMASSFRAC_ARRAY);
					sprintf(achFile,achBaseMask,msrOutName(msr),iStep);
					strncat(achFile,".FeMassFrac",256);
					msrOutArray(msr,achFile,OUT_IRONMASSFRAC_ARRAY);
					sprintf(achFile,achBaseMask,msrOutName(msr),iStep);
					strncat(achFile,".rhoform",256);
					msrOutArray(msr,achFile,OUT_DENSITY_ARRAY);
					sprintf(achFile,achBaseMask,msrOutName(msr),iStep);
					strncat(achFile,".Tform",256);
					msrOutArray(msr,achFile,OUT_U_ARRAY);
					sprintf(achFile,achBaseMask,msrOutName(msr),iStep);
					strncat(achFile,".rform",256);
					msrOutVector(msr,achFile,OUT_RFORM_VECTOR);
					sprintf(achFile,achBaseMask,msrOutName(msr),iStep);
					strncat(achFile,".vform",256);
					msrOutVector(msr,achFile,OUT_VFORM_VECTOR);
#ifdef SIMPLESF
					sprintf(achFile,achBaseMask,msrOutName(msr),iStep);
					strncat(achFile,".divv",256);
					msrOutArray(msr,achFile,OUT_DIVV_ARRAY);
					sprintf(achFile,achBaseMask,msrOutName(msr),iStep);
					strncat(achFile,".tCoolAgain",256);
					msrOutArray(msr,achFile,OUT_TCOOLAGAIN_ARRAY);
					sprintf(achFile,achBaseMask,msrOutName(msr),iStep);
					strncat(achFile,".mStar",256);
					msrOutArray(msr,achFile,OUT_MSTAR_ARRAY);
#endif
				    }
#endif
				if (msr->param.bDodtOutput) {
					msrReorder(msr);
					sprintf(achFile,achBaseMask,msrOutName(msr),iStep);
					strncat(achFile,".dt",256);
					msrOutArray(msr,achFile,OUT_DT_ARRAY);
					}
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
			if (iStop || iStep == msrSteps(msr) ||
				(msrCheckInterval(msr) && iStep%msrCheckInterval(msr) == 0)) {
				/*
				 ** Write a checkpoint.
				 */
#ifndef BENCHMARK
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
    
		if (msrRestart(msr)) {
			msrReorder(msr);
			sprintf(achFile,"%s.diag",msrOutName(msr));
#ifndef COLLISIONS
			msrWriteTipsy(msr,achFile,dTime);
#else
			msrWriteSS(msr,achFile,dTime);
#endif
			}

		if (msrDoGravity(msr)) {
			msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE );
			msrDomainDecomp(msr,0,1);
			msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE );
			msrUpdateSoft(msr,dTime);
			msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE );
			msrBuildTree(msr,0,dMass,0);
			msrMassCheck(msr,dMass,"After msrBuildTree in OutSingle Gravity");
			msrGravity(msr,0.0,msrDoSun(msr),&iSec,&dWMax,&dIMax,&dEMax,&nActive);
			msrMassCheck(msr,dMass,"After msrGravity in OutSingle Gravity");
			msrReorder(msr);
			msrMassCheck(msr,dMass,"After msrReorder in OutSingle Gravity");
			sprintf(achFile,"%s.accg",msrOutName(msr));
			msrOutVector(msr,achFile,OUT_ACCEL_VECTOR);
			msrMassCheck(msr,dMass,"After msrOutVector in OutSingle Gravity");
			sprintf(achFile,"%s.pot",msrOutName(msr));
			msrReorder(msr);
			msrOutArray(msr,achFile,OUT_POT_ARRAY);
			msrMassCheck(msr,dMass,"After msrOutArray in OutSingle Gravity");
		        }
#ifdef GASOLINE				
		if (msr->nGas > 0) {
		        msrActiveType(msr,TYPE_GAS,TYPE_ACTIVE|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
			msrDomainDecomp(msr,0,1);
			msrBuildTree(msr,1,-1.0,1);
			msrSmooth(msr,dTime,SMX_DENSITY,1);
			if (msr->param.bBulkViscosity) {
			  msrReSmooth(msr,dTime,SMX_DIVVORT,1);
			  msrGetGasPressure(msr);
			  msrReSmooth(msr,dTime,SMX_HKPRESSURETERMS,1);
			  } 
			else {
			  if (msr->param.bViscosityLimiter || msr->param.bShockTracker)
			    msrReSmooth(msr,dTime,SMX_DIVVORT,1);

			  msrSphViscosityLimiter(msr, dTime);
			  msrGetGasPressure(msr);
			  /*
			  msrReSmooth(msr,dTime,SMX_SPHPRESSURETERMS,1);
			  */
			  msrReSmooth(msr,dTime,SMX_SPHPRESSURE,1);
			  msrUpdateShockTracker(msr, 0.0);
			  msrReSmooth(msr,dTime,SMX_SPHVISCOSITY,1);
			  /*
 			  if (msr->param.bShockTracker)
			    msrReSmooth(msr,dTime,SMX_SHOCKTRACK,1);
			  */
			  /*
			  msrReSmooth(msr,dTime,SMX_SPHPRESSURETERMS,1);
			  */
			  msrReorder(msr);
			  sprintf(achFile,"%s.BSw",msrOutName(msr));
			  msrOutArray(msr,achFile,OUT_BALSARASWITCH_ARRAY);
			  sprintf(achFile,"%s.divv",msrOutName(msr));
			  msrOutArray(msr,achFile,OUT_DIVV_ARRAY);
			  sprintf(achFile,"%s.mumax",msrOutName(msr));
			  msrOutArray(msr,achFile,OUT_MUMAX_ARRAY);
			  if (msr->param.bShockTracker) {
			       sprintf(achFile,"%s.ST",msrOutName(msr));
			       msrOutArray(msr,achFile,OUT_SHOCKTRACKER_ARRAY);

			       sprintf(achFile,"%s.dch",msrOutName(msr));
			       msrOutArray(msr,achFile,OUT_DIVONCONH_ARRAY);
			       sprintf(achFile,"%s.dcx",msrOutName(msr));
			       msrOutArray(msr,achFile,OUT_DIVONCONX_ARRAY);

			       sprintf(achFile,"%s.divrhov",msrOutName(msr));
			       msrOutArray(msr,achFile,OUT_DIVRHOV_ARRAY);
			       sprintf(achFile,"%s.gradrho",msrOutName(msr));
			       msrOutVector(msr,achFile,OUT_GRADRHO_VECTOR);

			       sprintf(achFile,"%s.accp",msrOutName(msr));
			       msrOutVector(msr,achFile,OUT_ACCELPRES_VECTOR);
			  }

			  sprintf(achFile,"%s.acc",msrOutName(msr));
			  msrOutVector(msr,achFile,OUT_ACCEL_VECTOR);
			  sprintf(achFile,"%s.PdV",msrOutName(msr));
			  msrOutArray(msr,achFile,OUT_PDV_ARRAY);
 			  sprintf(achFile,"%s.PdVpres",msrOutName(msr));
			  msrOutArray(msr,achFile,OUT_PDVPRES_ARRAY);
			  sprintf(achFile,"%s.PdVvisc",msrOutName(msr));
			  msrOutArray(msr,achFile,OUT_PDVVISC_ARRAY);
			  }

			if (msr->param.bSphStep) {
		          fprintf(stderr,"Adding SphStep dt\n");
			  msrSphStep(msr,dTime);
			  msrReorder(msr);
			  sprintf(achFile,"%s.SPHdt",msrOutName(msr));
			  msrOutArray(msr,achFile,OUT_DT_ARRAY);
			  }
		        }

		if (msr->param.bDoSphhOutput) {
		        msrReorder(msr);
			sprintf(achFile,"%s.SPHH",msrOutName(msr));
			msrOutArray(msr,achFile,OUT_H_ARRAY);
		        }
#ifndef NOCOOLING				
				{
				int ArrayCnt = 0;
				char OutSuffix[20];
				int OutType;
				
				for (;;) {	
					CoolOutputArray( &msr->param.CoolParam, ArrayCnt, &OutType, OutSuffix );
					if (OutType == OUT_NULL) break;
					sprintf(achFile,"%s.%s",msrOutName(msr),OutSuffix);
					msrOutArray(msr,achFile, OutType);
                    ArrayCnt++;
					}
				}
#endif
#endif
		/*
		 ** Build tree, activating all particles first (just in case).
		 */
		if (msrDoDensity(msr) || msr->param.bDensityStep) {
			msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
			msrDomainDecomp(msr,0,1);
			msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
			msrBuildTree(msr,0,-1.0,1);
			msrMassCheck(msr,dMass,"After msrBuildTree in OutSingle Density");
			msrSmooth(msr,dTime,SMX_DENSITY,1);
			msrMassCheck(msr,dMass,"After msrSmooth in OutSingle Density");
		        } 
		if (msrDoDensity(msr)) {
			msrReorder(msr);
			msrMassCheck(msr,dMass,"After msrReorder in OutSingle Density");
			sprintf(achFile,"%s.den",msrOutName(msr));
			msrReorder(msr);
			msrOutArray(msr,achFile,OUT_DENSITY_ARRAY);
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
			sprintf(achFile,"%s.enc",msrOutName(msr));

			msrOutArray(msr,achFile,OUT_DENSITY_ARRAY);
			}
#endif
		        }

		if (msrDoGravity(msr)) {
			if (msr->param.bGravStep) {
			        fprintf(stderr,"Adding GravStep dt\n");
				msrGravStep(msr,dTime);
				}
			if (msr->param.bAccelStep) {
			        fprintf(stderr,"Adding AccelStep dt\n");
				msrAccelStep(msr,dTime);
				}
			}

		if (msr->param.bDensityStep) {
		    fprintf(stderr,"Adding DensStep dt\n");
		    msrDensityStep(msr,dTime);
		        }

		if (msr->param.bDeltaAccelStep) {
 		    fprintf(stderr,"Adding DeltaAccelStep dt\n");
			if (!msr->param.bDeltaAccelStepGasTree) {
			    msrActiveType(msr,TYPE_ALL,TYPE_TREEACTIVE);
			    msrBuildTree(msr,0,-1.0,1);
			    }
			else {
			    msrActiveType(msr,TYPE_GAS,TYPE_TREEACTIVE);
			    msrBuildTree(msr,1,-1.0,1);
			}
			msrActiveType(msr,TYPE_ALL,TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
			msrBuildTree(msr,0,-1.0,1);
		     msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
			/* This smooth sets dt directly -- hardwired coefficient */
			msrSmooth(msr,dTime,SMX_DELTAACCEL,0);
		    }
		/*
		msrDtToRung(msr,0,msrDelta(msr),1);
		msrRungStats(msr);
		*/
		msrReorder(msr);
		sprintf(achFile,"%s.dt",msrOutName(msr));
		msrOutArray(msr,achFile,OUT_DT_ARRAY);
		if(msr->param.iMaxRung > 1 && (msr->param.bDensityStep || msrDoGravity(msr))) {
			msrDtToRung(msr,0,msrDelta(msr),1);
			msrRungStats(msr);
			}
		}
	
	dfFinalize( msr->df );
	msrFinish(msr);
	mdlFinish(mdl);
	return 0;
	}
