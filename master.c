#include <stdio.h>
#include <stdlib.h>
#include <unistd.h> /* for unlink() */
#include <stddef.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <sys/stat.h>

#define max(A,B) ((A) > (B) ? (A) : (B))

#include <sys/param.h> /* for MAXPATHLEN and MAXHOSTNAMELEN, if available */
#ifndef MAXPATHLEN
#define MAXPATHLEN 256
#endif

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

#ifdef AMPI
#include "charm.h"
#define printf CmiPrintf
#endif

#ifdef COLLISIONS
#include "ssdefs.h" /* in turn includes ssio.h */
#include "collision.h"
#endif

#define LOCKFILE ".lockfile"	/* for safety lock */
#define STOPFILE "STOP"			/* for user interrupt */

#define NEWTIME
#ifdef NEWTIME 
double msrTime() {
	struct timeval tv;

	gettimeofday(&tv,NULL);
	return tv.tv_sec + tv.tv_usec*1e-6;
	}
#else
double msrTime() {
	return time(NULL);
	}
#endif

#define DEN_CGS_SYS 1.6831e6 /* multiply density in cgs by this to get
								density in system units (AU, M_Sun) */

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


void _msrExit(MSR msr,int status)
{
	MDL mdl=msr->mdl;
	msrFinish(msr);
	mdlFinish(mdl);
	exit(status);
	}


void
_msrMakePath(const char *dir,const char *base,char *path)
{
	/*
	 ** Prepends "dir" to "base" and returns the result in "path". It is the
	 ** caller's responsibility to ensure enough memory has been allocated
	 ** for "path".
	 */

	if (!path) return;
	path[0] = 0;
	if (dir) {
		strcat(path,dir);
		strcat(path,"/");
		}
	if (!base) return;
	strcat(path,base);
	}

void
_msrStripQuotes(const char achIn[],char achOut[])
{
	/* Removes leading and trailing double quotes from string */

	assert(achIn && achOut);
	if (strlen(achIn) > 0 && achIn[0] == '"')
		(void) strcpy(achOut,achIn + 1);
	else
		(void) strcpy(achOut,achIn);
	if (strlen(achOut) > 0 && achOut[strlen(achOut) - 1] == '"')
		achOut[strlen(achOut) - 1] = '\0';
	}

#ifdef SAND_PILE

void
_msrGetWallData(MSR msr,const char achFilenameWithQuotes[])
{
	FILE *fp;
	WALLS *w = &msr->param.CP.walls;
	char achFilename[MAXPATHLEN],achTmp[MAXPATHLEN];
	double dd;
	int i,di;

	assert(msr && achFilenameWithQuotes);
	_msrStripQuotes(achFilenameWithQuotes,achFilename);
	if (!strlen(achFilename)) {
		w->nWalls = 0;
		return;
		}
	_msrMakePath(msr->param.achDataSubPath,achFilename,achTmp);
	_msrMakePath(msr->lcl.pszDataPath,achTmp,achFilename);
	if (!(fp = fopen(achFilename,"r"))) {
		(void) fprintf(stderr,"Unable to open \"%s\"\n",achFilename);
		goto abort;
		}
	if (fscanf(fp,"%i",&w->nWalls) != 1) {
		(void) fprintf(stderr,"Expected no. walls in \"%s\"\n",achFilename);
		goto abort;
		}
	if (w->nWalls <= 0) {
		(void) fprintf(stderr,"Invalid no. walls in \"%s\"\n",achFilename);
		goto abort;
		}
	if (w->nWalls > MAX_NUM_WALLS) {
		(void) fprintf(stderr,"Number of walls (%i) exceeds maximum (%i)\n",
					   w->nWalls,MAX_NUM_WALLS);
		goto abort;
		}
	for (i=0;i<w->nWalls;i++) {
#ifdef TUMBLER
		if (fscanf(fp,"%lf%lf%lf%lf%lf%lf%lf%lf%lf%i",
			&w->wall[i].n[0],&w->wall[i].n[1],&w->wall[i].n[2],
			&w->wall[i].ndotp,&w->wall[i].radius,&w->wall[i].omega,
			&w->wall[i].dEpsN,&w->wall[i].dEpsT,
			&w->wall[i].hotParam,&w->wall[i].type) != 10) {
			(void) fprintf(stderr,"Invalid/missing data in \"%s\" (wall %i)\n",
						   achFilename,i);
			goto abort;
			}
		if (w->wall[i].radius < 0) {
			(void) fprintf(stderr,"Invalid radius (%g) in \"%s\", wall %i\n",
						   w->wall[i].radius,achFilename,i);
			goto abort;
			}
		if (w->wall[i].hotParam < 0) {
			(void) fprintf(stderr,"Invalid hotParam (%g) in \"%s\", wall %i\n",
						   w->wall[i].hotParam,achFilename,i);
			goto abort;
			}
		if (w->wall[i].type < 0 || w->wall[i].type > 1) {
			(void) fprintf(stderr,"Invalid wall type (%i) in \"%s\", wall %i\n",
						   w->wall[i].type,achFilename,i);
			goto abort;
			}
#else
		if (fscanf(fp,"%lf%lf%lf%lf%lf%lf%lf%lf%i",
				   &w->wall[i].x1,&dd,&w->wall[i].z1,
				   &w->wall[i].x2,&dd,&w->wall[i].z2,
				   &w->wall[i].dEpsN,&w->wall[i].dEpsT,&di) != 9) {
			(void) fprintf(stderr,"Invalid/missing data in \"%s\" (wall %i)\n",
						   achFilename,i);
			goto abort;
			}
#endif
		if (w->wall[i].dEpsN < 0 || w->wall[i].dEpsN > 1) {
			(void) fprintf(stderr,"Invalid epsn (%g) in \"%s\", wall %i\n",
						   w->wall[i].dEpsN,achFilename,i);
			goto abort;
			}
		if (w->wall[i].dEpsT < -1 || w->wall[i].dEpsT > 1) {
			(void) fprintf(stderr,"Invalid epst (%g) in \"%s\", wall %i\n",
						   w->wall[i].dEpsT,achFilename,i);
			goto abort;
			}
		}
	(void) fclose(fp);
	return;
 abort:
	if (fp) (void) fclose(fp);
	_msrExit(msr,1);
	}

#endif /* SAND_PILE */

#ifdef SPECIAL_PARTICLES

void
_msrGetSpecialData(MSR msr,const char achFilenameWithQuotes[])
{
	FILE *fp;
	struct parameters *p;
	SPECIAL_PARTICLE_DATA *s;
	char achFilename[MAXPATHLEN],achTmp[MAXPATHLEN];
	int i,n=0;

	assert(msr && achFilenameWithQuotes);
	p = &msr->param;
	_msrStripQuotes(achFilenameWithQuotes,achFilename);
	if (!strlen(achFilename)) {
		p->nSpecial = 0;
		return;
		}
	_msrMakePath(p->achDataSubPath,achFilename,achTmp);
	_msrMakePath(msr->lcl.pszDataPath,achTmp,achFilename);
	if (!(fp = fopen(achFilename,"r"))) {
		(void) fprintf(stderr,"Unable to open \"%s\"\n",achFilename);
		goto abort;
		}
	if (fscanf(fp,"%i",&p->nSpecial) != 1) {
		(void) fprintf(stderr,"Expected no. special particles in \"%s\"\n",achFilename);
		goto abort;
		}
	if (p->nSpecial <= 0) {
		(void) fprintf(stderr,"Invalid no. special particles in \"%s\"\n",achFilename);
		goto abort;
		}
	if (p->nSpecial > MAX_NUM_SPECIAL_PARTICLES) {
		(void) fprintf(stderr,"Number of special particles (%i) exceeds maximum (%i)\n",
					   p->nSpecial,MAX_NUM_SPECIAL_PARTICLES);
		goto abort;
		}
	s = p->sSpecialData;
	for (i=0;i<p->nSpecial;i++) {
		if (fscanf(fp,"%i%i",&p->iSpecialId[i],&s[i].iType) != 2) {
			(void) fprintf(stderr,"Invalid/missing data in \"%s\" (entry %i)\n",
						   achFilename,i);
			goto abort;
			}
		if (p->iSpecialId[i] < -1) {
			(void) fprintf(stderr,"Invalid original index (%i) in \"%s\", entry %i\n",
						   p->iSpecialId[i],achFilename,i);
			goto abort;
			}
		if (p->iSpecialId[i] == -1) {
			if (!msr->param.bHeliocentric) {
				puts("ERROR: Negative special particle index requires non-inertial frame");
				goto abort;
				}
			if (++n > 1) {
				puts("ERROR: Only one special particle can have index -1");
				goto abort;
				}
			}
		if (!((s[i].iType & SPECIAL_OBLATE) || (s[i].iType & SPECIAL_GR))) {
			(void) fprintf(stderr,"Invalid data type (%i) in \"%s\", entry %i\n",
						   s[i].iType,achFilename,i);
			goto abort;
			}
		if (s[i].iType & SPECIAL_OBLATE) {
			if (fscanf(fp,"%lf%lf%lf",&s[i].oblate.dRadEq,&s[i].oblate.J2,
					   &s[i].oblate.J4) != 3) {
				(void) fprintf(stderr,"Invalid/missing data in \"%s\" (entry %i)\n",
							   achFilename,i);
				goto abort;
				}
			}
		if (s[i].iType & SPECIAL_GR) {
			(void) fprintf(stderr,"GR not supported yet.\n");
			goto abort;
			}
		}
	(void) fclose(fp);
	return;
 abort:
	if (fp) (void) fclose(fp);
	_msrExit(msr,1);
	}

#endif /* SPECIAL_PARTICLES */

void msrInitialize(MSR *pmsr,MDL mdl,int argc,char **argv)
{
	MSR msr;
	int j,ret;
	int id,nDigits;
	struct inSetAdd inAdd;
	struct inLevelize inLvl;
	struct inGetMap inGM;

#ifdef COLLISIONS
	char achWallFile[MAXPATHLEN]; /* needed for SAND_PILE only */
	char achSpecialFile[MAXPATHLEN]; /* needed for SPECIAL_PARTICLES only */
#endif

	msr = (MSR)malloc(sizeof(struct msrContext));
	assert(msr != NULL);

	msr->bDumpFrame = 0;
	msr->df = NULL;

	msr->mdl = mdl;
	msr->pst = NULL;
	msr->lcl.pkd = NULL;
	*pmsr = msr;
	/*
	 ** default dTimeOld value
	 */
	msr->dTimeOld = 1e-20;
	csmInitialize(&msr->param.csm);
	/*
	 ** Now setup for the input parameters.
	 **
	 ** NOTE: nThreads & bDiag are parsed here, but the actual values are
	 ** read from the command line via mdlInitialize(). This means the
	 ** values of nThreads & bDiag read by prmAddParam() are ignored!
	 */
	prmInitialize(&msr->prm,_msrLeader,_msrTrailer);
	msr->param.nThreads = 1;
	prmAddParam(msr->prm,"nThreads",1,&msr->param.nThreads,sizeof(int),"sz",
				"<nThreads>");
	msr->param.bDiag = 0;
	prmAddParam(msr->prm,"bDiag",0,&msr->param.bDiag,sizeof(int),"d",
				"enable/disable per thread diagnostic output");
	msr->param.bOverwrite = 0;
	prmAddParam(msr->prm,"bOverwrite",0,&msr->param.bOverwrite,sizeof(int),
				"overwrite","enable/disable overwrite safety lock = +overwrite");
	msr->param.bVWarnings = 1;
	prmAddParam(msr->prm,"bVWarnings",0,&msr->param.bVWarnings,sizeof(int),
				"vwarnings","enable/disable warnings = +vwarnings");
	msr->param.bVStart = 1;
	prmAddParam(msr->prm,"bVStart",0,&msr->param.bVStart,sizeof(int),
				"vstart","enable/disable verbose start = +vstart");
	msr->param.bVStep = 1;
	prmAddParam(msr->prm,"bVStep",0,&msr->param.bVStep,sizeof(int),
				"vstep","enable/disable verbose step = +vstep");
	msr->param.bVRungStat = 1;
	prmAddParam(msr->prm,"bVRungStat",0,&msr->param.bVRungStat,sizeof(int),
				"vrungstat","enable/disable rung statistics = +vrungstat");
	msr->param.bVDetails = 0;
	prmAddParam(msr->prm,"bVDetails",0,&msr->param.bVDetails,sizeof(int),
				"vdetails","enable/disable verbose details = +vdetails");
	nDigits = 5;
	prmAddParam(msr->prm,"nDigits",1,&nDigits,sizeof(int),"nd",
				"<number of digits to use in output filenames> = 5");
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
	msr->param.iBinaryOutput = 0;
	prmAddParam(msr->prm,"iBinaryOutput",1,&msr->param.iBinaryOutput,sizeof(int),
				"binout","<array outputs 0 ascii, 1 float, 2 double, 3 FLOAT(internal)> = 0");
	msr->param.bPackedVector = 0;
	prmAddParam(msr->prm,"bPackedVector",0,&msr->param.bPackedVector,sizeof(int),
				"pvec","enable/disable packed vector outputs = +pvec");
	msr->param.bDoDensity = 1;
	prmAddParam(msr->prm,"bDoDensity",0,&msr->param.bDoDensity,sizeof(int),
				"den","enable/disable density outputs = +den");
	msr->param.iReadIOrder = 0;
	prmAddParam(msr->prm,"iReadIOrder",1,&msr->param.iReadIOrder,sizeof(int),
				"iordin","<array outputs 0 NO, 1 int, 2 long, 3 int (internal)> = 0");
	msr->param.bDoIOrderOutput = 0;
	prmAddParam(msr->prm,"bDoIOrderOutput",0,&msr->param.bDoIOrderOutput,
		    sizeof(int), "iordout","enable/disable iOrder outputs = -iordout");
	msr->param.bDohOutput = 0;
	prmAddParam(msr->prm,"bDohOutput",0,&msr->param.bDohOutput,sizeof(int),
				"hout","enable/disable h outputs = -hout");
	msr->param.bDoSphhOutput = 0;
	prmAddParam(msr->prm,"bDoSphhOutput",0,&msr->param.bDoSphhOutput,sizeof(int),
				"sphhout","enable/disable Sph h outputs = -sphhout");
	msr->param.bDodtOutput = 0;
	prmAddParam(msr->prm,"bDodtOutput",0,&msr->param.bDodtOutput,sizeof(int),
				"dtout","enable/disable dt outputs = -dtout");
	msr->param.nBucket = 8;
	prmAddParam(msr->prm,"nBucket",1,&msr->param.nBucket,sizeof(int),"b",
				"<max number of particles in a bucket> = 8");
	msr->param.iStartStep = 0;
	prmAddParam(msr->prm,"iStartStep",1,&msr->param.iStartStep,
				sizeof(int),"nstart","<initial step numbering> = 0");
	msr->param.nSteps = 0;
	prmAddParam(msr->prm,"nSteps",1,&msr->param.nSteps,sizeof(int),"n",
				"<number of timesteps> = 0");
	msr->param.iOutInterval = 0;
	prmAddParam(msr->prm,"iOutInterval",1,&msr->param.iOutInterval,sizeof(int),
				"oi","<number of timesteps between snapshots> = 0");
	msr->param.dDumpFrameStep = 0;
	prmAddParam(msr->prm,"dDumpFrameStep",2,&msr->param.dDumpFrameStep,sizeof(double),
				"dfi","<number of steps between dumped frames> = 0");
	msr->param.dDumpFrameTime = 0;
	prmAddParam(msr->prm,"dDumpFrameTime",2,&msr->param.dDumpFrameTime,sizeof(double),
				"dft","<number of timesteps between dumped frames> = 0");
	msr->param.iLogInterval = 10;
	prmAddParam(msr->prm,"iLogInterval",1,&msr->param.iLogInterval,sizeof(int),
				"ol","<number of timesteps between logfile outputs> = 10");
	msr->param.iCheckInterval = 10;
	prmAddParam(msr->prm,"iCheckInterval",1,&msr->param.iCheckInterval,sizeof(int),
				"oc","<number of timesteps between checkpoints> = 10");
	msr->param.iOrder = 4;
	prmAddParam(msr->prm,"iOrder",1,&msr->param.iOrder,sizeof(int),"or",
				"<multipole expansion order: 1, 2, 3 or 4> = 4");
	msr->param.bEwald = 1;
	prmAddParam(msr->prm,"bEwald",0,&msr->param.bEwald,sizeof(int),"ewald",
				"enable/disable Ewald correction = +ewald");
	msr->param.iEwOrder = 4;
	prmAddParam(msr->prm,"iEwOrder",1,&msr->param.iEwOrder,sizeof(int),"ewo",
				"<Ewald multipole expansion order: 1, 2, 3 or 4> = 4");
	msr->param.nReplicas = 0;
	prmAddParam(msr->prm,"nReplicas",1,&msr->param.nReplicas,sizeof(int),"nrep",
				"<nReplicas> = 0 for -p, or 1 for +p");
	msr->param.dSoft = 0.0;
	prmAddParam(msr->prm,"dSoft",2,&msr->param.dSoft,sizeof(double),"e",
				"<gravitational softening length> = 0.0");
	msr->param.dSoftMax = 0.0;
	prmAddParam(msr->prm,"dSoftMax",2,&msr->param.dSoftMax,sizeof(double),"eMax",
				"<maximum comoving gravitational softening length (abs or multiplier)> = 0.0");
	msr->param.bPhysicalSoft = 0;
	prmAddParam(msr->prm,"bPhysicalSoft",0,&msr->param.bPhysicalSoft,sizeof(int),"PhysSoft",
				"<Physical gravitational softening length> -PhysSoft");
	msr->param.bSoftMaxMul = 1;
	prmAddParam(msr->prm,"bSoftMaxMul",0,&msr->param.bSoftMaxMul,sizeof(int),"SMM",
				"<Use maximum comoving gravitational softening length as a multiplier> +SMM");
	msr->param.bVariableSoft = 0;
	prmAddParam(msr->prm,"bVariableSoft",0,&msr->param.bVariableSoft,sizeof(int),"VarSoft",
				"<Variable gravitational softening length> -VarSoft");
	msr->param.nSoftNbr = 32;
	prmAddParam(msr->prm,"nSoftNbr",1,&msr->param.nSoftNbr,sizeof(int),"VarSoft",
				"<Neighbours for Variable gravitational softening length> 32");
	msr->param.bSoftByType = 1;
	prmAddParam(msr->prm,"bSoftByType",0,&msr->param.bSoftByType,sizeof(int),"SBT",
				"<Variable gravitational softening length by Type> +SBT");
	msr->param.bDoSoftOutput = 0;
	prmAddParam(msr->prm,"bDoSoftOutput",0,&msr->param.bDoSoftOutput,sizeof(int),
				"softout","enable/disable soft outputs = -softout");

	msr->param.dDelta = 0.0;
	prmAddParam(msr->prm,"dDelta",2,&msr->param.dDelta,sizeof(double),"dt",
				"<time step>");
	msr->param.dEta = 0.1;
	prmAddParam(msr->prm,"dEta",2,&msr->param.dEta,sizeof(double),"eta",
				"<time step criterion> = 0.1");
	msr->param.dEtaDeltaAccel = 0.2;
	prmAddParam(msr->prm,"dEtaDeltaAccel",2,&msr->param.dEtaDeltaAccel,sizeof(double),"etadrda",
				"<drda time step criterion> = 0.2");
	msr->param.dEtaCourant = 0.4;
	prmAddParam(msr->prm,"dEtaCourant",2,&msr->param.dEtaCourant,sizeof(double),"etaC",
				"<Courant criterion> = 0.4");
	msr->param.dEtauDot = 0.25;
	prmAddParam(msr->prm,"dEtauDot",2,&msr->param.dEtauDot,sizeof(double),"etau",
				"<uDot criterion> = 0.25");
	msr->param.duDotLimit = -0.2;
	prmAddParam(msr->prm,"duDotLimit",2,&msr->param.duDotLimit,sizeof(double),"uDL",
				"<uDotLimit:  Treat udot/u < duDotLimit specially> = -0.2 < 0");
	msr->param.bGravStep = 0;
	prmAddParam(msr->prm,"bGravStep",0,&msr->param.bGravStep,sizeof(int),
				"gs","<Symmetric gravity timestepping (sqrt(r^3/(mi + mj)))>");
	msr->param.bEpsAccStep = 1;
	prmAddParam(msr->prm,"bEpsAccStep",0,&msr->param.bEpsAccStep,sizeof(int),
				"ea", "<Sqrt(Epsilon on a) timestepping>");
	msr->param.bSqrtPhiStep = 0;
	prmAddParam(msr->prm,"bSqrtPhiStep",0,&msr->param.bSqrtPhiStep,sizeof(int),
				"sphi", "<Sqrt(Phi) on a timestepping>");
	msr->param.bDensityStep = 0;
	prmAddParam(msr->prm,"bDensityStep",0,&msr->param.bDensityStep,sizeof(int),
				"isrho", "<Sqrt(1/Rho) timestepping>");
	msr->param.bDeltaAccelStep = 0;
	prmAddParam(msr->prm,"bDeltaAccelStep",0,&msr->param.bDeltaAccelStep,sizeof(int),
				"isdrda", "<Sqrt(dr/da) timestepping>");
	msr->param.bDeltaAccelStepGasTree = 0;
	prmAddParam(msr->prm,"bDeltaAccelStepGasTree",0,&msr->param.bDeltaAccelStepGasTree,sizeof(int),
				"isdrdagt", "<Sqrt(dr/da) timestepping via gas tree>");
	msr->param.nTruncateRung = 0;
	prmAddParam(msr->prm,"nTruncateRung",1,&msr->param.nTruncateRung,sizeof(int),"nTR",
				"<number of MaxRung particles to delete MaxRung> = 0");
	msr->param.bNonSymp = 1;
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
	msr->param.dTheta2 = msr->param.dTheta;
	prmAddParam(msr->prm,"dTheta",2,&msr->param.dTheta,sizeof(double),"theta",
				"<Barnes opening criterion> = 0.8");
	prmAddParam(msr->prm,"dTheta2",2,&msr->param.dTheta2,sizeof(double),
				"theta2","<Barnes opening criterion after a < daSwitchTheta> = 0.8");
	msr->param.daSwitchTheta = 1./3.;
	prmAddParam(msr->prm,"daSwitchTheta",2,&msr->param.daSwitchTheta,sizeof(double),"aSwitchTheta",
				"<a to switch theta at> = 1./3.");
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

	msr->param.bDoSinks = 0;
	prmAddParam(msr->prm,"bDoSinks",0,&msr->param.bDoSinks,sizeof(int),
				"sinks","enable/disable sinks = -sinks");
	msr->param.bSinkThermal = 0;
	prmAddParam(msr->prm,"bSinkThermal",0,&msr->param.bSinkThermal,sizeof(int),
				"tsinks","enable/disable thermal energy in sink calcs = -tsinks");
	msr->param.dSinkRadius = 0.0;
	prmAddParam(msr->prm,"dSinkRadius",2,&msr->param.dSinkRadius,sizeof(double),"sinkr",
				"<Sink Radius>");
	msr->param.dSinkBoundOrbitRadius = 0.0;
	prmAddParam(msr->prm,"dSinkBoundOrbitRadius",2,&msr->param.dSinkBoundOrbitRadius,sizeof(double),"sinkbor",
				"<Sink Bound Orbit Radius>");
	msr->param.dDeltaSink = msr->param.dDelta;
	prmAddParam(msr->prm,"dDeltaSink", 2, &msr->param.dDeltaSink,
		    sizeof(double), "dDeltaSink",
		    "<Minimum SF timestep in years> = dDelta");
	msr->param.iSinkRung = 0;
	prmAddParam(msr->prm,"iSinkRung", 2, &msr->param.iSinkRung,
		    sizeof(double), "iSinkRung",
		    "<Star Formation Rung> = 0");

	msr->param.dPeriod = 1.0;
	prmAddParam(msr->prm,"dPeriod",2,&msr->param.dPeriod,sizeof(double),"L",
				"<periodic box length> = 1.0");
	msr->param.dxPeriod = 1.0;
	prmAddParam(msr->prm,"dxPeriod",2,&msr->param.dxPeriod,sizeof(double),"Lx",
				"<periodic box length in x-dimension> = 1.0");
	msr->param.dyPeriod = 1.0;
	prmAddParam(msr->prm,"dyPeriod",2,&msr->param.dyPeriod,sizeof(double),"Ly",
				"<periodic box length in y-dimension> = 1.0");
	msr->param.dzPeriod = 1.0;
	prmAddParam(msr->prm,"dzPeriod",2,&msr->param.dzPeriod,sizeof(double),"Lz",
				"<periodic box length in z-dimension> = 1.0");
	msr->param.achInFile[0] = 0;
	prmAddParam(msr->prm,"achInFile",3,msr->param.achInFile,256,"I",
				"<input file name> (file in TIPSY binary format)");
#ifdef GASOLINE
	strcpy(msr->param.achOutName,"gasoline");
	prmAddParam(msr->prm,"achOutName",3,msr->param.achOutName,256,"o",
				"<output name for snapshots and logfile> = \"gasoline\"");
#else
	strcpy(msr->param.achOutName,"pkdgrav");
	prmAddParam(msr->prm,"achOutName",3,msr->param.achOutName,256,"o",
				"<output name for snapshots and logfile> = \"pkdgrav\"");
#endif
	msr->param.csm->bComove = 0;
	prmAddParam(msr->prm,"bComove",0,&msr->param.csm->bComove,sizeof(int),
				"cm", "enable/disable comoving coordinates = -cm");
	msr->param.csm->dHubble0 = 0.0;
	prmAddParam(msr->prm,"dHubble0",2,&msr->param.csm->dHubble0, 
				sizeof(double),"Hub", "<dHubble0> = 0.0");
	msr->param.csm->dOmega0 = 1.0;
	prmAddParam(msr->prm,"dOmega0",2,&msr->param.csm->dOmega0,
				sizeof(double),"Om", "<dOmega0> = 1.0");
	msr->param.csm->dLambda = 0.0;
	prmAddParam(msr->prm,"dLambda",2,&msr->param.csm->dLambda,
				sizeof(double),"Lambda", "<dLambda> = 0.0");
	msr->param.csm->dOmegaRad = 0.0;
	prmAddParam(msr->prm,"dOmegaRad",2,&msr->param.csm->dOmegaRad,
				sizeof(double),"Omrad", "<dOmegaRad> = 0.0");
	msr->param.csm->dOmegab = 0.0;
	prmAddParam(msr->prm,"dOmegab",2,&msr->param.csm->dOmegab,
				sizeof(double),"Omb", "<dOmegab> = 0.0");
	msr->param.csm->dQuintess = 0.0;
	prmAddParam(msr->prm,"dQuintess",2,&msr->param.csm->dQuintess,
				sizeof(double),"Quint",
		    "<dQuintessence (constant w = -1/2) > = 0.0");
	strcpy(msr->param.achDataSubPath,".");
	prmAddParam(msr->prm,"achDataSubPath",3,msr->param.achDataSubPath,256,
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
	msr->param.dCoolMaxDens = 1e8;
	prmAddParam(msr->prm,"dCoolMaxDens",2,&msr->param.dCoolMaxDens,
				sizeof(double),"scmd",
				"<Velocity Cooling Maximum Density> = 1e8");
	msr->param.bSymCool = 0;
	prmAddParam(msr->prm,"bSymCool",0,&msr->param.bSymCool,sizeof(int),NULL,
				NULL);
	msr->param.nGrowMass = 0;
	prmAddParam(msr->prm,"nGrowMass",1,&msr->param.nGrowMass,sizeof(int),
				"gmn","<number of particles to increase mass> = 0");
	msr->param.dGrowDeltaM = 0.0;
	prmAddParam(msr->prm,"dGrowDeltaM",2,&msr->param.dGrowDeltaM,
				sizeof(double),"gmdm","<Total growth in mass/particle> = 0.0");
	msr->param.dGrowStartT = 0.0;
	prmAddParam(msr->prm,"dGrowStartT",2,&msr->param.dGrowStartT,
				sizeof(double),"gmst","<Start time for growing mass> = 0.0");
	msr->param.dGrowEndT = 1.0;
	prmAddParam(msr->prm,"dGrowEndT",2,&msr->param.dGrowEndT,
				sizeof(double),"gmet","<End time for growing mass> = 1.0");
	msr->param.dFracNoDomainDecomp = 0.002;
	prmAddParam(msr->prm,"dFracNoDomainDecomp",2,&msr->param.dFracNoDomainDecomp,
				sizeof(double),"fndd",
				"<Fraction of Active Particles for no new DD> = 0.002");
	msr->param.dFracNoDomainDimChoice = 0.1;
	prmAddParam(msr->prm,"dFracNoDomainDimChoice",2,&msr->param.dFracNoDomainDimChoice,
				sizeof(double),"fnddc",
				"<Fraction of Active Particles for no new DD dimension choice> = 0.1");
	msr->param.dFracFastGas = 0.2;
	prmAddParam(msr->prm,"dFracFastGas",2,&msr->param.dFracFastGas,
				sizeof(double),"fndd",
				"<Fraction of Active Particles for Fast Gas> = 0.01");
	msr->param.dhMinOverSoft = 0.0;
	prmAddParam(msr->prm,"dhMinOverSoft",2,&msr->param.dhMinOverSoft,
				sizeof(double),"hmin",
				"<Minimum h as a fraction of Softening> = 0.0");
	msr->param.bDoGravity = 1;
	prmAddParam(msr->prm,"bDoGravity",0,&msr->param.bDoGravity,sizeof(int),"g",
				"enable/disable gravity (interparticle and external potentials) = +g");
	msr->param.bDoSelfGravity = 1;
	prmAddParam(msr->prm,"bDoSelfGravity",0,&msr->param.bDoSelfGravity,sizeof(int),"sg",
				"enable/disable interparticle self gravity = +sg");
	msr->param.bRungDD = 0;
	prmAddParam(msr->prm,"bRungDomainDecomp",0,&msr->param.bRungDD,sizeof(int),
				"RungDD","<Rung Domain Decomp> = 0");
	msr->param.dRungDDWeight = 1.0;
	prmAddParam(msr->prm,"dRungDDWeight",2,&msr->param.dRungDDWeight,sizeof(int),
				"RungDDWeight","<Rung Domain Decomp Weight> = 1.0");
	msr->param.bFandG = 0;
	prmAddParam(msr->prm,"bFandG",0,&msr->param.bFandG,sizeof(int),"fg",
				"use/don't use Kepler orbit drifts = -fg");
	msr->param.bHeliocentric = 0;
	prmAddParam(msr->prm,"bHeliocentric",0,&msr->param.bHeliocentric,
				sizeof(int),"hc","use/don't use Heliocentric coordinates = -hc");
	msr->param.dCentMass = 1.0;
	prmAddParam(msr->prm,"dCentMass",2,&msr->param.dCentMass,sizeof(double),
				"fgm","specifies the central mass for Keplerian orbits");
	msr->param.bLogHalo = 0;
	prmAddParam(msr->prm,"bLogHalo",0,&msr->param.bLogHalo,
				sizeof(int),"halo","use/don't use galaxy halo = -halo");
	msr->param.bHernquistSpheroid = 0;
	prmAddParam(msr->prm,"bHernquistSpheroid",0,&msr->param.bHernquistSpheroid,
				sizeof(int),"hspher","use/don't use galaxy Hernquist Spheroid = -hspher");
	msr->param.bNFWSpheroid = 0;
	prmAddParam(msr->prm,"bNFWSpheroid",0,&msr->param.bNFWSpheroid,
				sizeof(int),"NFWspher","use/don't use galaxy NFW Spheroid = -NFWspher");
	msr->param.bElliptical=0;
	prmAddParam(msr->prm,"bElliptical",0,&msr->param.bElliptical,
				sizeof(int),"elliptical","use/don't");
	msr->param.bEllipticalDarkNFW=0;
	prmAddParam(msr->prm,"bEllipticalDarkNFW",0,&msr->param.bEllipticalDarkNFW,
		    sizeof(int),"ellipticaldarknfw","use/dont");
	msr->param.bHomogSpheroid = 0;
	prmAddParam(msr->prm,"bHomogSpheroid",0,&msr->param.bHomogSpheroid,
				sizeof(int),"hspher","use/don't use galaxy Homog Spheroid = -homogspher");
	msr->param.bBodyForce = 0;
	prmAddParam(msr->prm,"bBodyForce",0,&msr->param.bBodyForce,
				sizeof(int),"bodyforce","use/don't use body force = -bf");
	msr->param.bMiyamotoDisk = 0;
	prmAddParam(msr->prm,"bMiyamotoDisk",0,&msr->param.bMiyamotoDisk,
				sizeof(int),"mdisk","use/don't use galaxy Miyamoto Disk = -mdisk");
	msr->param.bTimeVarying = 0;
	prmAddParam(msr->prm,"bTimeVarying",0,&msr->param.bTimeVarying,
				sizeof(int),"tvar","use/don't use the time varying externalpotential = -tvar");
	msr->param.bRotFrame = 0;
	prmAddParam(msr->prm,"bRotFrame",0,&msr->param.bRotFrame,
				sizeof(int),"rframe","use/don't use rotating frame = -rframe");
	msr->param.dOmega = 0;
	prmAddParam(msr->prm,"dOmega",2,&msr->param.dOmega,sizeof(double),
				"omega","<dOmega> = 0");
	msr->param.dOmegaDot = 0;
	prmAddParam(msr->prm,"dOmegaDot",2,&msr->param.dOmegaDot,sizeof(double),
				"omegadot","<dOmegaDot> = 0");
	msr->param.iWallRunTime = 0;
	prmAddParam(msr->prm,"iWallRunTime",1,&msr->param.iWallRunTime,
				sizeof(int),"wall",
				"<Maximum Wallclock time (in minutes) to run> = 0 = infinite");
	msr->param.dSunSoft = 0.0;
	prmAddParam(msr->prm,"dSunSoft", 2, &msr->param.dSunSoft,
		    sizeof(double), "sunSoft",
		    "<Softening length of the Sun in heliocentric coordinates> = 0.0");
#ifdef GASOLINE
	msr->param.bSphStep = 1;
	prmAddParam(msr->prm,"bSphStep",0,&msr->param.bSphStep,sizeof(int),
				"ss","<SPH timestepping>");
	msr->param.bDoGas = 1;
	prmAddParam(msr->prm,"bDoGas",0,&msr->param.bDoGas,sizeof(int),"gas",
				"calculate gas/don't calculate gas = +gas");
#ifndef GASOLINE
	msr->param.bDoGas = 0;
#endif
	msr->param.bGeometric = 0;
	prmAddParam(msr->prm,"bGeometric",0,&msr->param.bGeometric,sizeof(int),
				"geo","geometric/arithmetic mean to calc Grad(P/rho) = +geo");
	msr->param.bGasAdiabatic = 0;
	prmAddParam(msr->prm,"bGasAdiabatic",0,&msr->param.bGasAdiabatic,
				sizeof(int),"GasAdiabatic",
				"<Gas is Adiabatic> = +GasAdiabatic");
	msr->param.bGasIsothermal = 0;
	prmAddParam(msr->prm,"bGasIsothermal",0,&msr->param.bGasIsothermal,
				sizeof(int),"GasIsothermal",
				"<Gas is Isothermal> = +GasIsothermal");
	msr->param.bGasCooling = 0;
	prmAddParam(msr->prm,"bGasCooling",0,&msr->param.bGasCooling,
				sizeof(int),"GasCooling",
				"<Gas is Cooling> = +GasCooling");
	msr->param.iGasModel = GASMODEL_UNSET; /* Deprecated in for backwards compatibility */
	prmAddParam(msr->prm,"iGasModel",0,&msr->param.iGasModel,
				sizeof(int),"GasModel",
				"<Gas model employed> = 0 (Adiabatic)");
#ifndef NOCOOLING
	CoolAddParams( &msr->param.CoolParam, msr->prm );
#endif
	msr->param.dShockTrackerA = 0.16; 
	prmAddParam(msr->prm,"dShockTrackerA",2,&msr->param.dShockTrackerA,
				sizeof(double),"STA",
				"<Shock Tracker A constant> = 0.16");
	msr->param.dShockTrackerB = 0.4; 
	prmAddParam(msr->prm,"dShockTrackerB",2,&msr->param.dShockTrackerB,
				sizeof(double),"STB",
				"<Shock Tracker B constant> = 0.4");
	msr->param.dConstAlpha = 1.0; 	/* Default changed to 0.5 later if bBulkViscosity */
	prmAddParam(msr->prm,"dConstAlpha",2,&msr->param.dConstAlpha,
				sizeof(double),"alpha",
				"<Alpha constant in viscosity> = 1.0 or 0.5 (bBulkViscosity)");
	msr->param.dConstBeta = 2.0; 	/* Default changed to 0.5 later if bBulkViscosity */
	prmAddParam(msr->prm,"dConstBeta",2,&msr->param.dConstBeta,
				sizeof(double),"beta",
				"<Beta constant in viscosity> = 2.0 or 0.5 (bBulkViscosity)");
	msr->param.dConstGamma = 5.0/3.0;
	prmAddParam(msr->prm,"dConstGamma",2,&msr->param.dConstGamma,
				sizeof(double),"gamma",
				"<Ratio of specific heats> = 5/3");
	msr->param.dMeanMolWeight = 1.0;
	prmAddParam(msr->prm,"dMeanMolWeight",2,&msr->param.dMeanMolWeight,
				sizeof(double),"mmw",
				"<Mean molecular weight in amu> = 1.0");
	msr->param.dGasConst = 1.0;
	prmAddParam(msr->prm,"dGasConst",2,&msr->param.dGasConst,
				sizeof(double),"gcnst",
				"<Gas Constant>");
	msr->param.dKBoltzUnit = 1.0;
	prmAddParam(msr->prm,"dKBoltzUnit",2,&msr->param.dKBoltzUnit,
				sizeof(double),"gcnst",
				"<Boltzmann Constant in System Units>");
	msr->param.dMsolUnit = 1.0;
	prmAddParam(msr->prm,"dMsolUnit",2,&msr->param.dMsolUnit,
				sizeof(double),"msu",
				"<Solar mass/system mass unit>");
	msr->param.dKpcUnit = 1000.0;
	prmAddParam(msr->prm,"dKpcUnit",2,&msr->param.dKpcUnit,
				sizeof(double),"kpcu",
				"<Kiloparsec/system length unit>");
	msr->param.ddHonHLimit = 0.1;
	prmAddParam(msr->prm,"ddHonHLimit",2,&msr->param.ddHonHLimit,
				sizeof(double),"dhonh",
				"<|dH|/H Limiter> = 0.1");
	msr->param.bViscosityLimiter = 0;
	prmAddParam(msr->prm,"bViscosityLimiter",0,&msr->param.bViscosityLimiter,sizeof(int),
				"vlim","<Balsara Viscosity Limiter> = 0");
	msr->param.bViscosityLimitdt = 0;
	prmAddParam(msr->prm,"bViscosityLimitdt",0,&msr->param.bViscosityLimitdt,sizeof(int),
				"vlim","<Balsara Viscosity Limit dt> = 0");
	msr->param.bShockTracker = 0;
	prmAddParam(msr->prm,"bShockTracker",0,&msr->param.bShockTracker,sizeof(int),
				"st","<Shock Tracker> = 0");
	msr->param.bBulkViscosity = 0;
	prmAddParam(msr->prm,"bBulkViscosity",0,&msr->param.bBulkViscosity,sizeof(int),
				"bulk","<Bulk Viscosity> = 0");
	msr->param.bGasDomainDecomp = 0;
	prmAddParam(msr->prm,"bGasDomainDecomp",0,&msr->param.bGasDomainDecomp,sizeof(int),
				"gasDD","<Gas Domain Decomp> = 0");
	msr->param.bLowerSoundSpeed = 0;
	prmAddParam(msr->prm,"bLowerSoundSpeed",0,&msr->param.bLowerSoundSpeed,sizeof(int),
				"bLowc","<Lower Sound Speed> = 0");
	msr->param.bFastGas = 1;
	prmAddParam(msr->prm,"bFastGas",0,&msr->param.bFastGas,sizeof(int),
				"Fgas","<Fast Gas Method> = 1");
	msr->param.bStarForm = 0;
	prmAddParam(msr->prm,"bStarForm",0,&msr->param.bStarForm,sizeof(int),
				"stfm","<Star Forming> = 0");
	msr->param.bFeedBack = 0;
	prmAddParam(msr->prm,"bFeedBack",0,&msr->param.bFeedBack,sizeof(int),
				"fdbk","<Stars provide feedback> = 0");
#ifdef SIMPLESF
	msr->param.SSF_dComovingDenMin = 200.0;
	prmAddParam(msr->prm,"SSF_dComovingDenMin", 2, &msr->param.SSF_dComovingDenMin,
		    sizeof(double), "stODmin",
		    "<Minimum overdensity for forming stars> = 2");
	msr->param.SSF_dPhysDenMin =  7e-26;
	prmAddParam(msr->prm,"SSF_dPhysDenMin", 2, &msr->param.SSF_dPhysDenMin,
		    sizeof(double), "stPDmin",
		    "<Minimum physical density for forming stars (gm/cc)> =  7e-26");
    msr->param.SSF_dInitStarMass = 0;
    prmAddParam(msr->prm,"SSF_dInitStarMass", 2, &msr->param.SSF_dInitStarMass,
				sizeof(double), "stm0",
				"<Initial star mass> = 0");
	msr->param.SSF_dESNPerStarMass = 1.25e16;
	prmAddParam(msr->prm,"SSF_dESNPerStarMass", 2, &msr->param.SSF_dESNPerStarMass,
		    sizeof(double), "ESNPerStarMass",
		    "<ESN per star mass, erg per g of stars> = 1.25e16");
	msr->param.SSF_dTMax = 3e4;
	prmAddParam(msr->prm,"SSF_dTMax", 2, &msr->param.SSF_dTMax,
		    sizeof(double), "SSF_dTMax",
		    "<Maximum temperature for forming stars, K> = 3e4");
	msr->param.SSF_dEfficiency = 0.1;
	prmAddParam(msr->prm,"SSF_dEfficiency", 2, &msr->param.SSF_dEfficiency,
		    sizeof(double), "SSF_dEfficiency",
		    "<SF Efficiency> = 0.1");
	msr->param.SSF_dtCoolingShutoff = 30e6;
	prmAddParam(msr->prm,"SSF_dtCoolingShutoff", 2, &msr->param.SSF_dtCoolingShutoff,
		    sizeof(double), "SSF_dtCoolingShutoff",
		    "<SF Cooling Shutoff duration> = 30e6");
	msr->param.SSF_bdivv = 1;
	prmAddParam(msr->prm,"SSF_bdivv", 0, &msr->param.SSF_bdivv,
		    sizeof(int), "SSF_bdivv",
		    "<SF Use div v for star formation> = 1");
#endif /* SIMPLESF */

#ifdef STARFORM
	stfmInitialize(&msr->param.stfm);
	msr->param.stfm->dOverDenMin = 2.0;
	prmAddParam(msr->prm,"dOverDenMin", 2, &msr->param.stfm->dOverDenMin,
		    sizeof(double), "stODmin",
		    "<Minimum overdensity for forming stars> = 2");
	msr->param.stfm->dPhysDenMin = 0.1;
	prmAddParam(msr->prm,"dPhysDenMin", 2, &msr->param.stfm->dPhysDenMin,
		    sizeof(double), "stPDmin",
		    "<Minimum physical density for forming stars (atoms/cc)> = .1");
	msr->param.stfm->dStarEff = .3333;
	prmAddParam(msr->prm,"dStarEff", 2, &msr->param.stfm->dStarEff,
		    sizeof(double), "stEff",
		    "<Fraction of gas converted into stars per timestep> = .3333");
    msr->param.stfm->dInitStarMass = 0;
    prmAddParam(msr->prm,"dInitStarMass", 2, &msr->param.stfm->dInitStarMass,
				sizeof(double), "stm0",
				"<Initial star mass> = 0");
	msr->param.stfm->dMinMassFrac = .1;
	prmAddParam(msr->prm,"dMinMassFrac", 2, &msr->param.stfm->dMinMassFrac,
		    sizeof(double), "stMinFrac",
		    "<Minimum fraction of average mass of neighbour particles required for gas particles to avoid deletion> = .1");
	msr->param.stfm->dMinGasMass = 0.0;
	prmAddParam(msr->prm,"dMinGasMass", 2, &msr->param.stfm->dMinGasMass,
		    sizeof(double), "stMinGas",
		    "<Minimum mass of a gas particle> = 0.0");
	msr->param.stfm->dMaxStarMass = 0.0;
	prmAddParam(msr->prm,"dMaxStarMass", 2, &msr->param.stfm->dMaxStarMass,
		    sizeof(double), "stMaxStarMass",
		    "<Maximum amount of star mass a hybrid particle can contain = 0.0");
	msr->param.stfm->dCStar = 0.05;
	prmAddParam(msr->prm,"dCStar", 2, &msr->param.stfm->dCStar,
		    sizeof(double), "stCStar",
		    "<Star formation coefficient> = 0.1");
	msr->param.stfm->dTempMax = 3.0e4;
	prmAddParam(msr->prm,"dTempMax", 2, &msr->param.stfm->dTempMax,
		    sizeof(double), "stTempMax",
		    "<Maximum temperature at which star formation occurs> = 0.0");
	msr->param.stfm->dSoftMin = 1.0;
	prmAddParam(msr->prm,"dSoftMin", 2, &msr->param.stfm->dSoftMin,
		    sizeof(double), "stSoftMin",
		    "<Minimum softening for star formation> = 0.0");
	msr->param.dDeltaStarForm = msr->param.dDelta;
	prmAddParam(msr->prm,"dDeltaStarForm", 2, &msr->param.dDeltaStarForm,
		    sizeof(double), "dDeltaStarForm",
		    "<Minimum SF timestep in years> = dDelta");
	msr->param.bShortCoolShutoff = 0;
	prmAddParam(msr->prm,"bShortCoolShutoff", 0, &msr->param.bShortCoolShutoff,
		    sizeof(int), "bShortCoolShutoff",
		    "<Which cooling shutoff time to use> = long one");
	msr->param.bSmallSNSmooth = 1;
	prmAddParam(msr->prm,"bSmallSNSmooth", 0, &msr->param.bSmallSNSmooth,
		    sizeof(int), "bSmallSNSmooth",
		    "<smooth SN ejecta over blast or smoothing radius> = blast radius");
	msr->param.iStarFormRung = 0;
	prmAddParam(msr->prm,"iStarFormRung", 2, &msr->param.iStarFormRung,
		    sizeof(double), "iStarFormRung",
		    "<Star Formation Rung> = 0");

/* supernova constants */
	fbInitialize(&msr->param.fb);
        snInitialize(&msr->param.sn);
        snInitConstants(msr->param.sn);
	msr->param.sn->dESN = 1e51;
	prmAddParam(msr->prm,"dESN", 2, &msr->param.sn->dESN,
		    sizeof(double), "snESN",
		    "<Energy of supernova in ergs> = 1e51");

#endif /* STARFORM */
#endif /* GASOLINE */
#ifdef GLASS
	msr->param.dGlassDamper = 0.0;
	prmAddParam(msr->prm,"dGlassDamper",2,&msr->param.dGlassDamper,
		    sizeof(double),"dGlassDamper",
		    "Lose 0.0 dt velocity per step (no damping)");
	msr->param.dGlassPoverRhoL = 1.0;
	prmAddParam(msr->prm,"dGlassPoverRhoL",2,&msr->param.dGlassPoverRhoL,
				sizeof(double),"dGlassPoverRhoL","Left P on Rho = 1.0");
	msr->param.dGlassPoverRhoL = 1.0;
	prmAddParam(msr->prm,"dGlassPoverRhoR",2,&msr->param.dGlassPoverRhoR,
				sizeof(double),"dGlassPoverRhoR","Right P on Rho = 1.0");
	msr->param.dGlassxL = 0.0;
	prmAddParam(msr->prm,"dGlassxL",2,&msr->param.dGlassxL,sizeof(double),
				"dGlassxL","Left smoothing range = 0.0");
	msr->param.dGlassxR = 0.0;
	prmAddParam(msr->prm,"dGlassxR",2,&msr->param.dGlassxR,sizeof(double),
				"dGlassxR","Right smoothing range = 0.0");
	msr->param.dGlassVL = 0.0;
	prmAddParam(msr->prm,"dGlassVL",2,&msr->param.dGlassVL,sizeof(double),
				"dGlassVL","Left Max Random Velocity = 0.0");
	msr->param.dGlassVR = 0.0;
	prmAddParam(msr->prm,"dGlassVR",2,&msr->param.dGlassVR,sizeof(double),
				"dGlassVR","Right Max Random Velocity = 0.0");
#endif /* GLASS */
	msr->param.bPatch = 0;
	prmAddParam(msr->prm,"bPatch",0,&msr->param.bPatch,sizeof(int),
				"patch","enable/disable patch reference frame = -patch");
	msr->param.dOrbFreq = 0.0;
	prmAddParam(msr->prm,"dOrbFreq",2,&msr->param.dOrbFreq,
				sizeof(double),"orbfreq","<Patch orbit frequency>");
	msr->param.bSimpleGasDrag = 0;
	prmAddParam(msr->prm,"bSimpleGasDrag",0,&msr->param.bSimpleGasDrag,
				sizeof(int),"sgd","enable/disable simple gas drag = -sgd");
	msr->param.bEpstein = 1;
	prmAddParam(msr->prm,"bEpstein",0,&msr->param.bEpstein,sizeof(int),
				"epstein","enable/disable Epstein regime for gas drag");
	msr->param.dGamma = 1.0e-11;
	prmAddParam(msr->prm,"dGamma",2,&msr->param.dGamma,sizeof(double),
				"gamma","coefficient for inverse stopping time = 1.0e-11");
#ifdef COLLISIONS
	msr->param.bFindRejects = 0;
	prmAddParam(msr->prm,"bFindRejects",0,&msr->param.bFindRejects,
				sizeof(int),"rejects","enable/disable check for rejected ICs = -rejects");
	msr->param.iCollLogOption = 0;
	prmAddParam(msr->prm,"iCollLogOption",1,&msr->param.iCollLogOption,
				sizeof(int),"clog","<Collision log option> = 0");
	msr->param.dSmallStep = 0.0;
	prmAddParam(msr->prm,"dSmallStep",2,&msr->param.dSmallStep,
				sizeof(double),"sstep","<Rectilinear time-step>");
	msr->param.dxUnifGrav = msr->param.dyUnifGrav = msr->param.dzUnifGrav = 0;
	prmAddParam(msr->prm,"dxUnifGrav",2,&msr->param.dxUnifGrav,sizeof(double),
				"gx","<x component of uniform gravity field> = 0");
	prmAddParam(msr->prm,"dyUnifGrav",2,&msr->param.dyUnifGrav,sizeof(double),
				"gy","<y component of uniform gravity field> = 0");
	prmAddParam(msr->prm,"dzUnifGrav",2,&msr->param.dzUnifGrav,sizeof(double),
				"gz","<z component of uniform gravity field> = 0");
	msr->param.CP.iOutcomes = BOUNCE;
	prmAddParam(msr->prm,"iOutcomes",1,&msr->param.CP.iOutcomes,
				sizeof(int),"outcomes","<Allowed collision outcomes> = 0");
	msr->param.CP.dDensity = 0.0;
	prmAddParam(msr->prm,"dDensity",2,&msr->param.CP.dDensity,
				sizeof(double),"density","<Merged particle density> = 0");
	msr->param.CP.iBounceOption = ConstEps;
	prmAddParam(msr->prm,"iBounceOption",1,&msr->param.CP.iBounceOption,
				sizeof(int),"bopt","<Bounce option> = 0");
	msr->param.CP.dEpsN = 1.0;
	prmAddParam(msr->prm,"dEpsN",2,&msr->param.CP.dEpsN,
				sizeof(double),"epsn","<Coefficient of restitution> = 1");
	msr->param.CP.dEpsT = 1.0;
	prmAddParam(msr->prm,"dEpsT",2,&msr->param.CP.dEpsT,
				sizeof(double),"epst","<Coefficient of surface friction> = 1");
	msr->param.CP.iSlideOption = EscVel;
	prmAddParam(msr->prm,"iSlideOption",1,&msr->param.CP.iSlideOption,
				sizeof(int),"sopt","<Slide option> = 0");
	msr->param.CP.dSlideLimit = 0;
	prmAddParam(msr->prm,"dSlideLimit",2,&msr->param.CP.dSlideLimit,
				sizeof(double),"slide","<Sliding motion parameter> = 0");
	msr->param.CP.dSlideEpsN = 1.0;
	prmAddParam(msr->prm,"dSlideEpsN",2,&msr->param.CP.dSlideEpsN,
				sizeof(double),"sepsn","<epsn if speed less than minimum> = 1");
	msr->param.CP.dSlideEpsT = 1.0;
	prmAddParam(msr->prm,"dSlideEpsT",2,&msr->param.CP.dSlideEpsT,
				sizeof(double),"sepst","<epst if speed less than minimum> = 1");
	msr->param.CP.dCollapseLimit = 0;
	prmAddParam(msr->prm,"dCollapseLimit",2,&msr->param.CP.dCollapseLimit,
				sizeof(double),"collapse","<Inelastic collapse parameter> = 0");
	msr->param.CP.dCollapseEpsN = 1.0;
	prmAddParam(msr->prm,"dCollapseEpsN",2,&msr->param.CP.dCollapseEpsN,
				sizeof(double),"cepsn","<epsn for inelastic collapse> = 1");
	msr->param.CP.dCollapseEpsT = 1.0;
	prmAddParam(msr->prm,"dCollapseEpsT",2,&msr->param.CP.dCollapseEpsT,
				sizeof(double),"cepst","<epst for inelastic collapse> = 1");
	msr->param.CP.dCrushLimit = 0.0; /* i.e. no limit, set to DBL_MAX below */
	prmAddParam(msr->prm,"dCrushLimit",2,&msr->param.CP.dCrushLimit,
				sizeof(double),"crush","<Maximum impact speed squared> = 0");
	msr->param.CP.dCrushEpsN = 1.0;
	prmAddParam(msr->prm,"dCrushEpsN",2,&msr->param.CP.dCrushEpsN,
				sizeof(double),"cepsn","<epsn if speed greater than maximum> = 1");
	msr->param.CP.dCrushEpsT = 1.0;
	prmAddParam(msr->prm,"dCrushEpsT",2,&msr->param.CP.dCrushEpsT,
				sizeof(double),"cepst","<epst if speed greater than maximum> = 1");
	/*
	 ** The following parameters are only relevant to SAND_PILE and SPECIAL_PARTICLES,
	 ** but the parser is picky about unrecognized commands in parameter files, so
	 ** we'll include them in the general COLLISIONS model for now.
	 */
	strcpy(achWallFile,"");
	prmAddParam(msr->prm,"achWallFile",3,achWallFile,MAXPATHLEN,"walls",
				"<name of wall data file> = \"\"");
	strcpy(achSpecialFile,"");
	prmAddParam(msr->prm,"achSpecialFile",3,achSpecialFile,MAXPATHLEN,"special",
				"<name of special particle file> = \"\"");
#endif /* COLLISIONS */
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
		_msrExit(msr,1);
		}

#if defined(GASOLINE) || defined(COLLISIONS)
	assert(msr->param.bKDK); /*DEBUG DKD broken at the moment...*/
#endif

	if (nDigits < 1 || nDigits > 9) {
		(void) fprintf(stderr,"Unreasonable number of filename digits.\n");
		_msrExit(msr,1);
		}

	(void) sprintf(msr->param.achDigitMask,"%%s.%%0%ii",nDigits);

	if ( getenv("PKDGRAV_CHECKPOINT_FDL") == NULL ) { 
          fprintf(stderr,"PKDGRAV_CHECKPOINT_FDL environment variable not set\n");
	}
	/*
	 ** Don't allow periodic BC's for Kepler orbital problems.
	 ** It just doesn't make sense, does it?
	 */
	if (msr->param.bFandG) msr->param.bPeriodic = 0;
	/*
	 ** Make sure that we have some setting for nReplicas if bPeriodic is set.
	 */
	if (msr->param.bPeriodic && !prmSpecified(msr->prm,"nReplicas")) {
		msr->param.nReplicas = 1;
		}
	/*
	 ** Warn that we have a setting for nReplicas if bPeriodic NOT set.
	 */
	if (!msr->param.bPeriodic && msr->param.nReplicas != 0) {
		printf("WARNING: nReplicas set to non-zero value for non-periodic!\n");
		}

	if (!msr->param.achInFile[0] && !msr->param.bRestart) {
		puts("ERROR: no input file specified");
		_msrExit(msr,1);
		}

	if (msr->param.dTheta <= 0) {
		if (msr->param.dTheta == 0 && msr->param.bVWarnings)
			fprintf(stderr,"WARNING: Zero opening angle may cause numerical problems\n");
		else if (msr->param.dTheta < 0) {
			fprintf(stderr,"ERROR: Opening angle must be non-negative\n");
			_msrExit(msr,1);
			}
		}

	msr->nThreads = mdlThreads(mdl);

	if (msr->param.bDoSelfGravity && !msr->param.bDoGravity)
		fprintf(stderr,"WARNING: May need bDoGravity on for bDoSelfGravity to work\n");

	/*
	 ** Always set bCannonical = 1 if bComove == 0
	 */
	if (!msrComove(msr)) {
		if (!msr->param.bCannonical)
			printf("WARNING: bCannonical reset to 1 for non-comoving (bComove == 0)\n");
		msr->param.bCannonical = 1;
		} 
	/* 
	 * Softening 
	 */
	
	if (msr->param.bPhysicalSoft || msr->param.bVariableSoft) {
	  if (msr->param.bPhysicalSoft && !msrComove(msr)) {
	    printf("WARNING: bPhysicalSoft reset to 0 for non-comoving (bComove == 0)\n");
	    msr->param.bPhysicalSoft = 0;
	  }
#ifndef CHANGESOFT
	  fprintf(stderr,"ERROR: You must compile with -DCHANGESOFT to use changing softening options\n");
	  _msrExit(msr,1);
#endif
	  if (msr->param.bVariableSoft && !prmSpecified(msr->prm,"bDoSoftOutput")) msr->param.bDoSoftOutput=1;
  
	  if (msr->param.bPhysicalSoft && msr->param.bVariableSoft) {
	    fprintf(stderr,"ERROR: You may only choose one of Physical or Variable softening\n");
	    _msrExit(msr,1);
	  }
	}

	/*
	 ** Determine the period of the box that we are using.
	 ** Set the new d[xyz]Period parameters which are now used instead
	 ** of a single dPeriod, but we still want to have compatibility
	 ** with the old method of setting dPeriod.
	 */
	if (prmSpecified(msr->prm,"dPeriod") && 
		!prmSpecified(msr->prm,"dxPeriod")) {
		msr->param.dxPeriod = msr->param.dPeriod;
		}
	if (prmSpecified(msr->prm,"dPeriod") && 
		!prmSpecified(msr->prm,"dyPeriod")) {
		msr->param.dyPeriod = msr->param.dPeriod;
		}
	if (prmSpecified(msr->prm,"dPeriod") && 
		!prmSpecified(msr->prm,"dzPeriod")) {
		msr->param.dzPeriod = msr->param.dPeriod;
		}
	/*
	 ** Periodic boundary conditions can be disabled along any of the
	 ** x,y,z axes by specifying a period of zero for the given axis.
	 ** Internally, the period is set to infinity (Cf. pkdBucketWalk()
	 ** and pkdDrift(); also the INTERSECT() macro in smooth.h).
	 */
	if (msr->param.dPeriod  == 0) msr->param.dPeriod  = FLOAT_MAXVAL;
	if (msr->param.dxPeriod == 0) msr->param.dxPeriod = FLOAT_MAXVAL;
	if (msr->param.dyPeriod == 0) msr->param.dyPeriod = FLOAT_MAXVAL;
	if (msr->param.dzPeriod == 0) msr->param.dzPeriod = FLOAT_MAXVAL;
#ifdef GASOLINE
	assert(msr->param.duDotLimit <= 0);
	if (msr->param.bBulkViscosity) {
		if (!prmSpecified(msr->prm,"dConstAlpha"))
			msr->param.dConstAlpha=0.5;
		if (!prmSpecified(msr->prm,"dConstBeta"))
			msr->param.dConstBeta=0.5;
		}
#ifndef SHOCKTRACK
	if (msr->param.bShockTracker != 0) {
	        fprintf(stderr,"Compile with -DSHOCKTRACK for Shock Tracking.\n");
			assert(0);
	        }
#endif
#endif
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
		if (!prmSpecified(msr->prm,"dTheta2")) 
			msr->param.dTheta2 = msr->param.dTheta;
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
		if (!prmSpecified(msr->prm,"dTheta2")) 
			msr->param.dTheta2 = msr->param.dTheta;
		}
	/*
	 ** Initialize comove variables.
	 */
	msr->nMaxOuts = 100;
	msr->pdOutTime = malloc(msr->nMaxOuts*sizeof(double));
	assert(msr->pdOutTime != NULL);
	msr->nOuts = 0;

	/*
	 ** Check timestepping.
	 */

	if (msr->param.iMaxRung < 1) {
		msr->param.iMaxRung = 1;
		if (msr->param.bVWarnings)
			(void) fprintf(stderr,"WARNING: iMaxRung set to 1\n");
		}
		if(msr->param.iMaxRung > 29) {
			fprintf(stderr, "WARNING: Cannot have %d rungs, reducing to 29 (the maximum)\n", msr->param.iMaxRung);
			msr->param.iMaxRung = 29;
		}

	if (msr->param.bGravStep && !msr->param.bDoGravity) {
		puts("ERROR: need gravity to use gravity stepping...");
		_msrExit(msr,1);
		}
	if (msr->param.bEpsAccStep || msr->param.bSqrtPhiStep) {
		msr->param.bAccelStep = 1;
		}
#ifdef GASOLINE
	assert(msr->param.bSphStep);
#endif

#define SECONDSPERYEAR   31557600.
#ifdef GASOLINE
#define msrSetGasModel( iModel ) { \
  if (msr->param.iGasModel == GASMODEL_UNSET) msr->param.iGasModel = iModel; \
  else { \
    fprintf( stderr, "More than one Gas Model specified [%i and %i]", msr->param.iGasModel, iModel ); \
    assert(0); \
  } \
}
	if (msr->param.bGasAdiabatic    ) msrSetGasModel( GASMODEL_ADIABATIC );
	if (msr->param.bGasIsothermal   ) msrSetGasModel( GASMODEL_ISOTHERMAL );
	if (msr->param.bGasCooling      ) msrSetGasModel( GASMODEL_COOLING );

	if (msr->param.iGasModel == GASMODEL_UNSET) {
		msr->param.iGasModel = GASMODEL_ADIABATIC;
		fprintf( stderr, "Defaulting to Adiabatic Gas Model."); 
		}	

	switch(msr->param.iGasModel) {
	case GASMODEL_ADIABATIC:
		msr->param.bGasAdiabatic = 1;	
		break;
	case GASMODEL_ISOTHERMAL:
		msr->param.bGasIsothermal = 1;
		break;
	case GASMODEL_COOLING:
		msr->param.bGasCooling = 1;
		break;
		}

	if(msr->param.iGasModel == GASMODEL_COOLING) {
		assert (prmSpecified(msr->prm, "dMsolUnit") &&
				prmSpecified(msr->prm, "dKpcUnit"));
		}
	/* bolzman constant in cgs */
#define KBOLTZ	1.38e-16
	/* mass of hydrogen atom in grams */
#define MHYDR 1.67e-24
	/* solar mass in grams */
#define MSOLG 1.99e33
	/* G in cgs */
#define GCGS 6.67e-8
	/* kiloparsec in centimeters */
#define KPCCM 3.085678e21
	/*
	 ** Convert kboltz/mhydrogen to system units, assuming that
	 ** G == 1.
	 */
	if(prmSpecified(msr->prm, "dMsolUnit") &&
	   prmSpecified(msr->prm, "dKpcUnit")) {
		msr->param.dGasConst = msr->param.dKpcUnit*KPCCM*KBOLTZ
			/MHYDR/GCGS/msr->param.dMsolUnit/MSOLG;
		/* code energy per unit mass --> erg per g */
		msr->param.dErgPerGmUnit = GCGS*msr->param.dMsolUnit*MSOLG/(msr->param.dKpcUnit*KPCCM);
		/* code density --> g per cc */
		msr->param.dGmPerCcUnit = (msr->param.dMsolUnit*MSOLG)/pow(msr->param.dKpcUnit*KPCCM,3.0);
		/* code time --> seconds */
		msr->param.dSecUnit = sqrt(1/(msr->param.dGmPerCcUnit*GCGS));
		/* code comoving density --> g per cc = msr->param.dGmPerCcUnit (1+z)^3 */
		msr->param.dComovingGmPerCcUnit = msr->param.dGmPerCcUnit;
		}

	if(msr->param.bStarForm) {
	    assert (prmSpecified(msr->prm, "dMsolUnit") &&
		    prmSpecified(msr->prm, "dKpcUnit"));

		if (!(msr->param.iGasModel == GASMODEL_COOLING)) {
			fprintf(stderr,"Warning: You are not running a cooling"
					"EOS with starformation\n");
			}

#ifdef STARFORM
	    assert((msr->param.stfm->dStarEff > 0 && 
                msr->param.stfm->dStarEff < 1) ||
			   msr->param.stfm->dInitStarMass > 0);
	    assert(msr->param.stfm->dMinMassFrac > 0 && 
                msr->param.stfm->dMinMassFrac < 1);
		if (msr->param.stfm->dInitStarMass > 0) {
/*			if (msr->param.stfm->dMinGasMass <= 0) */
 			  /* Only allow 10% underweight star particles */
				msr->param.stfm->dMinGasMass = 0.9*msr->param.stfm->dInitStarMass;
			}
		else
			assert(msr->param.stfm->dMinGasMass > 0);

	    msr->param.stfm->dSecUnit = msr->param.dSecUnit;
	    msr->param.stfm->dGmPerCcUnit = msr->param.dGmPerCcUnit;
	    msr->param.stfm->dGmUnit = msr->param.dMsolUnit*MSOLG;
	    msr->param.stfm->dErgUnit =
		GCGS*pow(msr->param.dMsolUnit*MSOLG, 2.0)
		/(msr->param.dKpcUnit*KPCCM);
            msr->param.dKBoltzUnit = msr->param.dKpcUnit*KPCCM*KBOLTZ
			/GCGS/msr->param.dMsolUnit/MSOLG/msr->param.dMsolUnit/MSOLG;
	    /* convert to system units */
	    msr->param.stfm->dPhysDenMin *= MHYDR/msr->param.stfm->dGmPerCcUnit;
            msr->param.dtCoolingShutoff *= SECONDSPERYEAR/msr->param.dSecUnit;
            if( prmSpecified(msr->prm, "dDeltaStarForm") ){
                msr->param.dDeltaStarForm *= SECONDSPERYEAR/msr->param.dSecUnit;
                msr->param.stfm->dDeltaT = msr->param.dDeltaStarForm;
                }

	    msr->param.fb->dSecUnit = msr->param.dSecUnit;
	    msr->param.fb->dGmUnit = msr->param.dMsolUnit*MSOLG;
	    msr->param.fb->dErgPerGmUnit = msr->param.dErgPerGmUnit;
	    msr->param.fb->dInitStarMass = msr->param.stfm->dInitStarMass;
#endif /* STARFORM */
#ifdef SIMPLESF		
		assert(msr->param.SSF_dInitStarMass > 0);
#endif
	    }
	
	    
#endif /* GASOLINE */

#ifdef ROT_FRAME
	if (msr->param.bVWarnings && !msr->param.bRotFrame) {
		puts("WARNING: ROT_FRAME set without bRotFrame");
		}
#else
	if (msr->param.bRotFrame) {
		puts("ERROR: bRotFrame set without ROT_FRAME");
		_msrExit(msr,1);
		}
#endif

#ifdef COLLISIONS
	/*
	 ** Parameter checks and initialization.
	 */
#ifdef GASOLINE
	puts("ERROR: can't mix COLLISIONS and GASOLINE!");
	_msrExit(msr,1);
#endif
	if (!msr->param.bKDK) {
		puts("ERROR: must use KDK scheme for collisions");
		_msrExit(msr,1);
		}
#ifdef SAND_PILE
	if (msr->param.bVWarnings && msr->param.nSmooth < 1)
		puts("WARNING: collision detection disabled (nSmooth < 1)");
#else
	if (msr->param.nSmooth < 1) {
		puts("ERROR: nSmooth must be positive");
		_msrExit(msr,1);
		}
	if (msr->param.bVWarnings && msr->param.nSmooth < 2)
		puts("WARNING: collision detection disabled (nSmooth < 2)");
#endif
	if (!msr->param.bHeliocentric) msr->param.dCentMass = 0; /* just in case */
	assert(msr->param.dCentMass >= 0);
	if (msr->param.bFandG) {
#ifdef SLIDING_PATCH
		assert(0);
#endif
#ifdef SAND_PILE
		assert(0);
#endif
		assert(msr->param.bHeliocentric);
		if (!msr->param.bCannonical) {
			puts("ERROR: must use cannonical momentum in FandG collision model");
			_msrExit(msr,1);
			}
#ifdef OLD_KEPLER /* override for now */
		if (msr->param.iMaxRung > 1) {
			puts("ERROR: multistepping not currently supported in FandG collision model");
			_msrExit(msr,1);
			}
#endif
		if (msr->param.dSmallStep < 0) {
			puts("ERROR: dSmallStep cannot be negative");
			_msrExit(msr,1);
			}
		if (msr->param.bVWarnings && msr->param.dSmallStep == 0)
			puts("WARNING: encounter detection disabled (dSmallStep = 0)");
		if (msr->param.dSmallStep > msr->param.dDelta) {
			puts("ERROR: inner step must be less than dDelta");
			_msrExit(msr,1);
			}
		}
	if (msr->param.bPatch) {
		if (msr->param.dOrbFreq <= 0) {
			puts("ERROR: must specify positive patch orbit frequency");
			_msrExit(msr,1);
			}
		if (!msr->param.bPeriodic) {
			puts("ERROR: must use periodic BCs for patch model");
			_msrExit(msr,1);
			}
		assert(msr->param.dxPeriod > 0);
		assert(msr->param.dyPeriod > 0);
		assert(msr->param.dzPeriod > 0);
		if (msr->param.dyPeriod == FLOAT_MAXVAL) {
			puts("ERROR: must specify positive y period");
			_msrExit(msr,1);
			}
		assert(msr->param.nReplicas > 0);
		if (msr->param.bEwald) {
			puts("ERROR: cannot use Ewald correction in patch model");
			_msrExit(msr,1);
			}
		}
#ifdef SAND_PILE
	if (msr->param.iCheckInterval) {
		puts("ERROR: Unable to checkpoint with SAND_PILE"); /* bStuck... */
		_msrExit(msr,1);
		}
#endif
#ifdef SLIDING_PATCH
	if (!msr->param.bPatch) {
		puts("WARNING: SLIDING_PATCH set without bPatch");
		}
	if (msr->param.bPatch && msr->param.dxPeriod == FLOAT_MAXVAL) {
		puts("ERROR: must specify positive x period");
		_msrExit(msr,1);
		}
	if (msr->param.bPatch && msr->param.bHeliocentric) {
		puts("ERROR: Patch and Helocentric are incompatible");
		_msrExit(msr,1);
		}
#else
	if (msr->param.bPatch) {
		puts("ERROR: bPatch set without SLIDING_PATCH");
		_msrExit(msr,1);
		}
#endif
	switch (msr->param.iCollLogOption) {
	case COLL_LOG_NONE:
		break;
	case COLL_LOG_VERBOSE:
		(void) strcpy(msr->param.achCollLog,COLL_LOG_TXT);
		break;
	case COLL_LOG_TERSE:
		(void) strcpy(msr->param.achCollLog,COLL_LOG_BIN);
		break;
	default:
		puts("ERROR: Invalid collision log option");
		_msrExit(msr,1);
		}
	if (msr->param.iCollLogOption && msr->param.bVStart)
		printf("Collision log: \"%s\"\n",msr->param.achCollLog);
#ifdef SAND_PILE
	if (msr->param.nSmooth >= 1) {
#else
	if (msr->param.nSmooth >= 2) {
#endif
		COLLISION_PARAMS *CP = &msr->param.CP;
		if (!(CP->iOutcomes & (MERGE | BOUNCE | FRAG))) {
			puts("ERROR: must specify one of MERGE/BOUNCE/FRAG");
			_msrExit(msr,1);
			}
		if (!(CP->iOutcomes & MERGE) && CP->dDensity)
			puts("WARNING: merged particle density ignored (no merging)");
		if ((CP->iOutcomes & MERGE) && CP->dDensity < 0) {
			puts("ERROR: invalid merge density");
			_msrExit(msr,1);
			}
		if (!(CP->iOutcomes & BOUNCE) && CP->dCollapseLimit) {
			puts("WARNING: collapse detection disabled (no bouncing)");
			CP->dCollapseLimit = 0;
			}
		if (CP->iOutcomes & BOUNCE) {
			if (CP->iBounceOption < ConstEps ||
				CP->iBounceOption > Glancing) {
				puts("ERROR: invalid bounce option");
				_msrExit(msr,1);
				}
			if (CP->dSlideLimit < 0 || CP->dCollapseLimit < 0 ||
				CP->dCrushLimit < 0) {
				puts("ERROR: invalid slide, collapse, and/or crush limit");
				_msrExit(msr,1);
				}
			if (CP->dCrushLimit == 0)
				CP->dCrushLimit = DBL_MAX;
			if (CP->dEpsN <= 0 || CP->dEpsN > 1 ||
				(CP->dSlideLimit > 0 &&
				 (CP->dSlideEpsN <= 0 ||
				  CP->dSlideEpsN > 1)) ||
				(CP->dCrushLimit < DBL_MAX &&
				 (CP->dCrushEpsN <= 0 ||
				  CP->dCrushEpsN > 1))) {
				puts("ERROR: coef of rest must be > 0 and <= 1");
				_msrExit(msr,1);
				}
			if (CP->dEpsT < -1 || CP->dEpsT > 1 ||
				(CP->dSlideLimit > 0 &&
				 (CP->dSlideEpsT < -1 ||
				  CP->dSlideEpsT > 1)) ||
				(CP->dCrushLimit < DBL_MAX &&
				 (CP->dCrushEpsT < -1 ||
				  CP->dCrushEpsT > 1))) {
				puts("ERROR: coef of surf frict must be >= -1 and <= 1");
				_msrExit(msr,1);
				}
			}
		if (CP->iOutcomes != MERGE && (CP->iOutcomes & MERGE) &&
			!msr->param.bDoSelfGravity) {
			puts("ERROR: Need interparticle gravity for conditional merging");
			_msrExit(msr,1);
			}
		}
	msr->param.CP.dDensity *= DEN_CGS_SYS; /* convert: cgs to system units */
	if (msr->param.CP.iSlideOption == MaxTrv) {
		double g = sqrt(msr->param.dxUnifGrav*msr->param.dxUnifGrav +
						msr->param.dyUnifGrav*msr->param.dyUnifGrav +
						msr->param.dzUnifGrav*msr->param.dzUnifGrav);
		msr->param.CP.dSlideLimit = sqrt(2*g*msr->param.CP.dSlideLimit);
		}
	msr->param.CP.dSlideLimit2 =
		msr->param.CP.dSlideLimit*msr->param.CP.dSlideLimit;
#endif /* COLLISIONS */

#ifdef AGGS
#ifndef COLLISIONS
	puts("ERROR: For aggregates, collisions must be defined.");
	_msrExit(msr, 1);
#endif
#endif

#ifdef SAND_PILE
	_msrGetWallData(msr,achWallFile);
#endif
#ifdef SPECIAL_PARTICLES
	_msrGetSpecialData(msr,achSpecialFile);
#endif

	pstInitialize(&msr->pst,msr->mdl,&msr->lcl);

	pstAddServices(msr->pst,msr->mdl);
	/*
	 ** Create the processor subset tree.
	 */
	for (id=1;id<msr->nThreads;++id) {
		if (msr->param.bVDetails) printf("Adding %d to the pst\n",id);
		inAdd.id = id;
		pstSetAdd(msr->pst,&inAdd,sizeof(inAdd),NULL,NULL);
		}
	if (msr->param.bVDetails) printf("\n");
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
	/*
	 ** Mark the Domain Decompositon as not done
	 */
	msr->bDoneDomainDecomp = 0;
	msr->iLastRungDomainDecomp = 0;
	msr->nRung = (int *) malloc( (msr->param.iMaxRung+1)*sizeof(int) );
	}


void msrLogParams(MSR msr,FILE *fp)
{
	double z, testDelta;
	int i;

#ifdef __DATE__
#ifdef __TIME__
	fprintf(fp,"# Compiled: %s %s\n",__DATE__,__TIME__);
#endif
#endif
	fprintf(fp,"# Preprocessor macros:");
#ifdef GASOLINE
	fprintf(fp," GASOLINE");
#endif
#ifdef STARFORM
	fprintf(fp," STARFORM");
#endif
#ifdef SIMPLESF
	fprintf(fp," SIMPLESF");
#endif
#ifdef LARGEFBALL
	fprintf(fp," LARGEFBALL");
#endif
#ifdef SHOCKTRACK
	fprintf(fp," SHOCKTRACK");
#endif
#ifdef PEAKEDKERNEL
	fprintf(fp," PEAKEDKERNEL");
#endif
#ifdef CHANGESOFT
 	fprintf(fp," CHANGESOFT");
#endif
#ifdef NOCOOLING
 	fprintf(fp," NOCOOLING");
#endif
#ifdef COOLING_COSMO
 	fprintf(fp," COOLING_COSMO");
#endif
#ifdef COOLING_PLANET
 	fprintf(fp," COOLING_PLANET");
#endif
#ifdef GLASS
	fprintf(fp," GLASS");
#endif
#ifdef HSHRINK
	fprintf(fp," HSHRINK");
#endif
#ifdef DEBUG
	fprintf(fp," DEBUG");
#endif
#ifdef ALTSPH
	fprintf(fp," ALTSPH");
#endif
#ifdef SUPERCOOL
	fprintf(fp," SUPERCOOL");
#endif
#ifdef ROT_FRAME
	fprintf(fp," ROT_FRAME");
#endif
#ifdef COLLISIONS
	fprintf(fp," COLLISIONS");
#endif
#ifdef SPECIAL_PARTICLES
	fprintf(fp," SPECIAL_PARTICLES");
#endif
#ifdef FIX_COLLAPSE
	fprintf(fp," FIX_COLLAPSE");
#endif
#ifdef SLIDING_PATCH
	fprintf(fp," SLIDING_PATCH");
#endif
#ifdef SAND_PILE
	fprintf(fp," SAND_PILE");
#endif
#ifdef TUMBLER
	fprintf(fp," TUMBLER");
#endif
#ifdef SIMPLE_GAS_DRAG
	fprintf(fp," SIMPLE_GAS_DRAG");
#endif
#ifdef _REENTRANT
	fprintf(fp," _REENTRANT");
#endif
#ifdef CRAY_T3D
	fprintf(fp," CRAY_T3D");
#endif
#ifdef PRES_MONAGHAN
	fprintf(fp," PRES_MONAGHAN");
#endif 
#ifdef PRES_HK
	fprintf(fp," PRES_HK");
#endif 
#ifdef MAXHOSTNAMELEN
	{
	char hostname[MAXHOSTNAMELEN];
	fprintf(fp,"\n# Master host: ");
	if (gethostname(hostname,MAXHOSTNAMELEN))
		fprintf(fp,"unknown");
	else
		fprintf(fp,"%s",hostname);
	}
#endif
	fprintf(fp,"\n# N: %d",msr->N);
	fprintf(fp," nThreads: %d",msr->param.nThreads);
	fprintf(fp," bDiag: %d",msr->param.bDiag);
	fprintf(fp," Verbosity flags: (%d,%d,%d,%d,%d)",msr->param.bVWarnings,
			msr->param.bVStart,msr->param.bVStep,msr->param.bVRungStat,
			msr->param.bVDetails);
	fprintf(fp,"\n# bPeriodic: %d",msr->param.bPeriodic);
	fprintf(fp," bRestart: %d",msr->param.bRestart);
	fprintf(fp," bComove: %d",msr->param.csm->bComove);
	fprintf(fp,"\n# bParaRead: %d",msr->param.bParaRead);
	fprintf(fp," bParaWrite: %d",msr->param.bParaWrite);
	fprintf(fp," bCannonical: %d",msr->param.bCannonical);
	fprintf(fp," bStandard: %d",msr->param.bStandard);
	fprintf(fp,"\n# bKDK: %d",msr->param.bKDK);
	fprintf(fp," nBucket: %d",msr->param.nBucket);
	fprintf(fp," iOutInterval(%d,%d): %d",msr->param.iBinaryOutput,msr->param.bPackedVector,msr->param.iOutInterval);
	fprintf(fp," dDumpFrameStep: %g",msr->param.dDumpFrameStep);
	fprintf(fp," dDumpFrameTime: %g",msr->param.dDumpFrameTime);
	fprintf(fp," iLogInterval: %d",msr->param.iLogInterval);
	fprintf(fp,"\n# iCheckInterval: %d",msr->param.iCheckInterval);
	fprintf(fp," iOrder: %d",msr->param.iOrder);
	fprintf(fp," iEwOrder: %d",msr->param.iEwOrder);
	fprintf(fp," nReplicas: %d",msr->param.nReplicas);
	fprintf(fp,"\n# dEwCut: %f",msr->param.dEwCut);
	fprintf(fp," dEwhCut: %f",msr->param.dEwhCut);
	fprintf(fp,"\n# iStartStep: %d",msr->param.iStartStep);
	fprintf(fp," nSteps: %d",msr->param.nSteps);
	fprintf(fp," nSmooth: %d",msr->param.nSmooth);
	fprintf(fp," dExtraStore: %f",msr->param.dExtraStore);
	if (prmSpecified(msr->prm,"dSoft"))
		fprintf(fp," dSoft: %g",msr->param.dSoft);
	else
		fprintf(fp," dSoft: input");
	fprintf(fp,"\n# bPhysicalSoft: %d",msr->param.bPhysicalSoft);
	fprintf(fp," bVariableSoft: %d",msr->param.bVariableSoft);
	fprintf(fp," nSoftNbr: %d",msr->param.nSoftNbr);
	fprintf(fp," bSoftByType: %d",msr->param.bSoftByType);
	fprintf(fp," bSoftMaxMul: %d",msr->param.bSoftMaxMul);
	fprintf(fp," dSoftMax: %g",msr->param.dSoftMax);
	fprintf(fp," bDoSoftOutput: %d",msr->param.bDoSoftOutput);
	fprintf(fp,"\n# dDelta: %g",msr->param.dDelta);
	fprintf(fp," dEta: %g",msr->param.dEta);
	fprintf(fp," dEtaDeltaAccel: %g",msr->param.dEtaDeltaAccel);
	fprintf(fp," dEtaCourant: %g",msr->param.dEtaCourant);
	fprintf(fp," iMaxRung: %d",msr->param.iMaxRung);
	fprintf(fp,"\n# bGravStep: %d",msr->param.bGravStep);
	fprintf(fp," bEpsAccStep: %d",msr->param.bEpsAccStep);
	fprintf(fp," bSqrtPhiStep: %d",msr->param.bSqrtPhiStep);
	fprintf(fp," bDensityStep: %d",msr->param.bDensityStep);
	fprintf(fp," bDeltaAccelStep: %d",msr->param.bDeltaAccelStep);
	fprintf(fp," (gt): %d",msr->param.bDeltaAccelStepGasTree);
	fprintf(fp," nTruncateRung: %d",msr->param.nTruncateRung);
	fprintf(fp," bNonSymp: %d",msr->param.bNonSymp);
	fprintf(fp,"\n# bDoGravity: %d",msr->param.bDoGravity);
	fprintf(fp," bDoSelfGravity: %d",msr->param.bDoSelfGravity);
	fprintf(fp," bFandG: %d",msr->param.bFandG);
	fprintf(fp," bHeliocentric: %d",msr->param.bHeliocentric);
	fprintf(fp," dCentMass: %g",msr->param.dCentMass);
	fprintf(fp," dSunSoft: %g",msr->param.dSunSoft);
	fprintf(fp,"\n# bLogHalo: %d",msr->param.bLogHalo );
	fprintf(fp," bHernquistSpheroid: %d",msr->param.bHernquistSpheroid );
	fprintf(fp," bNFWSpheroid: %d",msr->param.bNFWSpheroid );
	fprintf(fp," bHomogSpheroid: %d",msr->param.bHomogSpheroid );
	fprintf(fp," bBodyForce: %d",msr->param.bBodyForce );
	fprintf(fp," bMiyamotoDisk: %d",msr->param.bMiyamotoDisk );
	fprintf(fp," bTimeVarying: %d",msr->param.bTimeVarying );
	fprintf(fp,"\n# bRotFrame: %d",msr->param.bRotFrame);
	fprintf(fp," dOmega: %g",msr->param.dOmega);
	fprintf(fp," dOmegaDot: %g",msr->param.dOmegaDot);
	fprintf(fp," bDoSinks: %d",msr->param.bDoSinks );
	fprintf(fp," bSinksThermal: %d",msr->param.bSinkThermal );
	fprintf(fp," dSinkRadius: %g",msr->param.dSinkRadius);
	fprintf(fp," dSinkBoundOrbitRadius: %g",msr->param.dSinkBoundOrbitRadius);
	fprintf(fp,"\n# dFracNoDomainDecomp: %g",msr->param.dFracNoDomainDecomp);
	fprintf(fp," dFracNoDomainDimChoice: %g",msr->param.dFracNoDomainDimChoice);
	fprintf(fp," bFastGas: %d",msr->param.bFastGas);
	fprintf(fp," dFracFastGas: %g",msr->param.dFracFastGas);
	fprintf(fp," dhMinOverSoft: %g",msr->param.dhMinOverSoft);
	fprintf(fp," bRungDD: %d",msr->param.bRungDD);
	fprintf(fp," dRungDDWeight: %g ",msr->param.dRungDDWeight);
	fprintf(fp,"\n# nTruncateRung: %d",msr->param.nTruncateRung);
	fprintf(fp," bLowerSoundSpeed: %d",msr->param.bLowerSoundSpeed);
	fprintf(fp," bShockTracker: %d",msr->param.bShockTracker);
	fprintf(fp," dShockTrackerA: %f",msr->param.dShockTrackerA);
	fprintf(fp," dShockTrackerB: %f",msr->param.dShockTrackerB);
	fprintf(fp,"\n# GROWMASS: nGrowMass: %d",msr->param.nGrowMass);
	fprintf(fp," dGrowDeltaM: %g",msr->param.dGrowDeltaM);
	fprintf(fp," dGrowStartT: %g",msr->param.dGrowStartT);
	fprintf(fp," dGrowEndT: %g",msr->param.dGrowEndT);
#ifdef GASOLINE
	fprintf(fp,"\n# SPH: bDoGas: %d",msr->param.bDoGas);	
	fprintf(fp," bGeometric: %d",msr->param.bGeometric);
	/* fprintf(fp," iGasModel: %d",msr->param.iGasModel); // Deprecated usage */
	fprintf(fp," bGasAdiabatic: %d",msr->param.bGasAdiabatic);	
	fprintf(fp," bGasIsothermal: %d",msr->param.bGasIsothermal);	
	fprintf(fp," bGasCooling: %d",msr->param.bGasCooling);	
	fprintf(fp," dConstAlpha: %g",msr->param.dConstAlpha);
	fprintf(fp," dConstBeta: %g",msr->param.dConstBeta);
	fprintf(fp,"\n# dConstGamma: %g",msr->param.dConstGamma);
	fprintf(fp," dMeanMolWeight: %g",msr->param.dMeanMolWeight);
	fprintf(fp," dGasConst: %g",msr->param.dGasConst);
	fprintf(fp," dKBoltzUnit: %g",msr->param.dKBoltzUnit);
	fprintf(fp," dMsolUnit: %g",msr->param.dMsolUnit);
	fprintf(fp," dKpcUnit: %g",msr->param.dKpcUnit);
	fprintf(fp," ddHonHLimit: %g",msr->param.ddHonHLimit);
	fprintf(fp,"\n# bViscosityLimiter: %d",msr->param.bViscosityLimiter);
	fprintf(fp," bBulkViscosity: %d",msr->param.bBulkViscosity);
	fprintf(fp," bGasDomainDecomp: %d",msr->param.bGasDomainDecomp);
	fprintf(fp," bSphStep: %d",msr->param.bSphStep);
	fprintf(fp,"\n#bSN: %d",msr->param.bSN);
	fprintf(fp," dSNRhoCut: %g",msr->param.dSNRhoCut);
 	fprintf(fp," dSNTMin: %g",msr->param.dSNTMin);
        fprintf(fp," dSNTMax: %g",msr->param.dSNTMax);
	fprintf(fp," dSNMetalCut: %g",msr->param.dSNMetalCut);
	fprintf(fp," dSNHeatFraction: %g",msr->param.dSNHeatFraction);
#ifndef NOCOOLING
	CoolLogParams( &msr->param.CoolParam, fp );
#endif
#endif
#ifdef STARFORM
	fprintf(fp,"\n# Star Formation: bStarForm: %d",msr->param.bStarForm);
	fprintf(fp," bFeedBack: %d",msr->param.bFeedBack);	
	fprintf(fp," dOverDenMin: %g",msr->param.stfm->dOverDenMin);
	fprintf(fp," dPhysDenMin: %g",msr->param.stfm->dPhysDenMin);
	fprintf(fp," dStarEff: %g",msr->param.stfm->dStarEff);
	fprintf(fp," dCStar: %g",msr->param.stfm->dCStar);
	fprintf(fp," dTempMax: %g",msr->param.stfm->dTempMax);
	fprintf(fp," dSoftMin: %g",msr->param.stfm->dSoftMin);
	fprintf(fp," dMinMassFrac: %g",msr->param.stfm->dMinMassFrac);
	fprintf(fp," dMinGasMass: %g",msr->param.stfm->dMinGasMass);
	fprintf(fp," dMaxStarMass: %g",msr->param.stfm->dMaxStarMass);
	fprintf(fp," dESN: %g",msr->param.sn->dESN);
	fprintf(fp," bShortCoolShutoff: %i",msr->param.bShortCoolShutoff);
	fprintf(fp," bSmallSNSmooth: %i",msr->param.bSmallSNSmooth);
        for ( testDelta = msr->param.dDelta; 
            testDelta >= msr->param.dDeltaStarForm && 
            msr->param.dDeltaStarForm != 0.0; testDelta *= 0.5 ){
                    if ( !(prmSpecified(msr->prm,"iStarFormRung")) )
                        msr->param.iStarFormRung++;
                    }
        if ( testDelta <= msr->param.dDelta ){ 
            fprintf(fp," dDeltaStarForm (set): %g, effectively: %g = %g yrs, iStarFormRung: %i",
                    msr->param.dDeltaStarForm, testDelta,
                    testDelta*msr->param.dSecUnit/SECONDSPERYEAR,
                    msr->param.iStarFormRung );
            msr->param.stfm->dDeltaT = msr->param.dDeltaStarForm = testDelta;
            }
        else if ( msr->param.dDeltaStarForm == 0.0 ) {
            fprintf(fp," dDeltaStarForm (set): %g, effectively: 0.0 = 0.0 yrs, iStarFormRung: maxRung",
                    msr->param.dDeltaStarForm );
            msr->param.iStarFormRung = msr->param.iMaxRung;
            }
        else {
            fprintf(fp," dDeltaStarForm (set): %g, effectively:  NO STARS WILL FORM", msr->param.dDeltaStarForm);
            }

#endif
  	  if (msr->param.bDoSinks) {
        for ( testDelta = msr->param.dDelta; 
            testDelta >= msr->param.dDeltaSink && 
            msr->param.dDeltaSink != 0.0; testDelta *= 0.5 ){
                    if ( !(prmSpecified(msr->prm,"iSinkRung")) )
                        msr->param.iSinkRung++;
                    }
        if ( msr->param.dDeltaSink == 0.0 ) {
            fprintf(fp," dDeltaSink (set): %g, effectively: 0.0 = 0.0 yrs, iSinkRung: maxRung",
                    msr->param.dDeltaSink );
            msr->param.iSinkRung = msr->param.iMaxRung;
            }
        else if ( testDelta <= msr->param.dDelta ){ 
            fprintf(fp," dDeltaSink (set): %g, effectively: %g = %g yrs, iSinkRung: %i",
                    msr->param.dDeltaSink, testDelta,
                    testDelta*msr->param.dSecUnit/SECONDSPERYEAR,
                    msr->param.iSinkRung );
            msr->param.dDeltaSink = testDelta;
            }
        else {
            fprintf(fp," dDeltaSink (set): %g, effectively:  NO SINK ACCRETION", msr->param.dDeltaSink);
            }
	  }

#ifdef SIMPLESF
	fprintf(fp,"\n# SSF: bStarForm: %d",msr->param.bStarForm);
	fprintf(fp," bFeedBack: %d",msr->param.bFeedBack);	
	fprintf(fp," SSF_dEfficiency: %g",msr->param.SSF_dEfficiency);
	fprintf(fp," SSF_dTMax: %g",msr->param.SSF_dTMax);
	fprintf(fp," SSF_dPhysDenMin: %g",msr->param.SSF_dPhysDenMin);
	fprintf(fp," SSF_dComovingDenMin: %g",msr->param.SSF_dComovingDenMin);
	fprintf(fp," SSF_dESNPerStarMass: %g",msr->param.SSF_dESNPerStarMass);
	fprintf(fp," SSF_dInitStarMass: %g",msr->param.SSF_dInitStarMass);
	fprintf(fp," SSF_dtCoolingShutoff: %g",msr->param.SSF_dtCoolingShutoff);
	fprintf(fp," SSF_bdivv: %d",msr->param.SSF_bdivv);
#endif
	fprintf(fp,"\n# bPatch: %d",msr->param.bPatch);
	fprintf(fp," dOrbFreq: %g",msr->param.dOrbFreq);
	fprintf(fp,"\n# bSimpleGasDrag: %d",msr->param.bSimpleGasDrag);
	fprintf(fp," bEpstein: %d",msr->param.bEpstein);
	fprintf(fp," dGamma: %g",msr->param.dGamma);
#ifdef COLLISIONS
	fprintf(fp,"\n# Collisions...");
	fprintf(fp," bFindRejects: %d",msr->param.bFindRejects);
	fprintf(fp," iCollLogOption: %d",msr->param.iCollLogOption);
	fprintf(fp," dSmallStep: %g",msr->param.dSmallStep);
	fprintf(fp,"\n# dxUnifGrav: %g",msr->param.dxUnifGrav);
	fprintf(fp," dyUnifGrav: %g",msr->param.dyUnifGrav);
	fprintf(fp," dzUnifGrav: %g",msr->param.dzUnifGrav);
    fprintf(fp,"\n# iOutcomes: %d",msr->param.CP.iOutcomes);
	fprintf(fp," dDensity: %g",msr->param.CP.dDensity/DEN_CGS_SYS);
	fprintf(fp," iBounceOption: %d",msr->param.CP.iBounceOption);
    fprintf(fp," dEpsN: %g",msr->param.CP.dEpsN);
    fprintf(fp," dEpsT: %g",msr->param.CP.dEpsT);
    fprintf(fp,"\n# iSlideOption: %d",msr->param.CP.iSlideOption);
	fprintf(fp," dSlideLimit: %g",msr->param.CP.dSlideLimit);
    fprintf(fp," dSlideEpsN: %g",msr->param.CP.dSlideEpsN);
    fprintf(fp," dSlideEpsT: %g",msr->param.CP.dSlideEpsT);
    fprintf(fp,"\n# dCollapseLimit: %g",msr->param.CP.dCollapseLimit);
    fprintf(fp," dCollapseEpsN: %g",msr->param.CP.dCollapseEpsN);
    fprintf(fp," dCollapseEpsT: %g",msr->param.CP.dCollapseEpsT);
    fprintf(fp,"\n# dCrushLimit: %g",msr->param.CP.dCrushLimit);
    fprintf(fp," dCrushEpsN: %g",msr->param.CP.dCrushEpsN);
    fprintf(fp," dCrushEpsT: %g",msr->param.CP.dCrushEpsT);
#endif
#ifdef SPECIAL_PARTICLES
	{
	SPECIAL_PARTICLE_DATA *s;
	int i;
	fprintf(fp,"\n# nSpecial: %i",msr->param.nSpecial);
	for (i=0;i<msr->param.nSpecial;i++) {
		s = &msr->param.sSpecialData[i];
		fprintf(fp,"\n# Special %i: org_idx: %i type: %i ",i,
				msr->param.iSpecialId[i],s->iType);
		if (s->iType & SPECIAL_OBLATE)
			fprintf(fp,"dRadEq: %g J2: %g J4: %g",s->oblate.dRadEq,
					s->oblate.J2,s->oblate.J4);
		if (s->iType & SPECIAL_GR)
			fprintf(fp,"UNSUPPORTED");
		}
	}
#endif
#ifdef SAND_PILE
	{
	WALLS *w = &msr->param.CP.walls;
	int i;
	fprintf(fp,"\n# nWalls: %i",w->nWalls);
	for (i=0;i<w->nWalls;i++) {
#ifdef TUMBLER
		fprintf(fp,"\n# Wall %i: nx: %g ny: %g nz: %g",
			i,w->wall[i].n[0],w->wall[i].n[1],w->wall[i].n[2]);
		fprintf(fp," ndotp: %g radius: %g omega:%g",
			w->wall[i].ndotp,w->wall[i].radius,w->wall[i].omega);
		fprintf(fp," dEpsN: %g dEpsT: %g hotParam: %g type: %i",
			w->wall[i].dEpsN,w->wall[i].dEpsT,w->wall[i].hotParam,w->wall[i].type);
#else
		fprintf(fp,"\n# Wall %i: x1: %g z1: %g x2: %g z2: %g",
				i,w->wall[i].x1,w->wall[i].z1,w->wall[i].x2,w->wall[i].z2);
		fprintf(fp," dEpsN: %g dEpsT: %g",w->wall[i].dEpsN,w->wall[i].dEpsT);
#endif
		}
	}
#endif
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
	fprintf(fp," dxPeriod: %g",
			msr->param.dxPeriod >= FLOAT_MAXVAL ? 0 : msr->param.dxPeriod);
	fprintf(fp," dyPeriod: %g",
			msr->param.dyPeriod >= FLOAT_MAXVAL ? 0 : msr->param.dyPeriod);
	fprintf(fp," dzPeriod: %g",
			msr->param.dzPeriod >= FLOAT_MAXVAL ? 0 : msr->param.dzPeriod);
	fprintf(fp,"\n# dHubble0: %g",msr->param.csm->dHubble0);
	fprintf(fp," dOmega0: %g",msr->param.csm->dOmega0);
	fprintf(fp," dLambda: %g",msr->param.csm->dLambda);
	fprintf(fp," dOmegaRad: %g",msr->param.csm->dOmegaRad);
	fprintf(fp," dOmegab: %g",msr->param.csm->dOmegab);
	fprintf(fp," dQuintess: %g",msr->param.csm->dQuintess);
	fprintf(fp,"\n# achInFile: %s",msr->param.achInFile);
	fprintf(fp,"\n# achOutName: %s",msr->param.achOutName); 
	fprintf(fp,"\n# achDataSubPath: %s",msr->param.achDataSubPath);
	if (msrComove(msr)) {
		fprintf(fp,"\n# RedOut:");
		if (msr->nOuts == 0) fprintf(fp," none");
		for (i=0;i<msr->nOuts;i++) {
			if (i%5 == 0) fprintf(fp,"\n#   ");
			z = 1.0/csmTime2Exp(msr->param.csm, msr->pdOutTime[i]) - 1.0;
			fprintf(fp," %f",z);
			}
		fprintf(fp,"\n");
		}
	else {
		fprintf(fp,"\n# TimeOut:");
		if (msr->nOuts == 0) fprintf(fp," none");
		for (i=0;i<msr->nOuts;i++) {
			if (i%5 == 0) fprintf(fp,"\n#   ");
			fprintf(fp," %f",msr->pdOutTime[i]);
			}
		fprintf(fp,"\n");
		}
	}

int
msrGetLock(MSR msr)
{
	/*
	 ** Attempts to lock run directory to prevent overwriting. If an old lock
	 ** is detected with the same achOutName, an abort is signaled. Otherwise
	 ** a new lock is created. The bOverwrite parameter flag can be used to
	 ** suppress lock checking.
	 */

	FILE *fp = NULL;
	char achTmp[256],achFile[256];

	_msrMakePath(msr->param.achDataSubPath,LOCKFILE,achTmp);
	_msrMakePath(msr->lcl.pszDataPath,achTmp,achFile);
	if (!msr->param.bOverwrite && (fp = fopen(achFile,"r"))) {
		(void) fscanf(fp,"%s",achTmp);
		(void) fclose(fp);
		if (!strcmp(msr->param.achOutName,achTmp)) {
			(void) printf("ABORT: %s detected.\nPlease ensure data is safe to "
						  "overwrite. Delete lockfile and try again.\n",achFile);
			return 0;
			}
		}
	if (!(fp = fopen(achFile,"w"))) {
		if (msr->param.bOverwrite && msr->param.bVWarnings) {
			(void) printf("WARNING: Unable to create %s...ignored.\n",achFile);
			return 1;
			}
		else {
			(void) printf("Unable to create %s\n",achFile);
			return 0;
			}
		}
	(void) fprintf(fp,msr->param.achOutName);
	(void) fclose(fp);
	return 1;
	}

int
msrCheckForStop(MSR msr)
{
	/*
	 ** Checks for existence of STOPFILE in run directory. If found, the file
	 ** is removed and the return status is set to 1, otherwise 0.
	 */

	static char achFile[256];
	static int first_call = 1;

	FILE *fp = NULL;

	if (first_call) {
		char achTmp[256];
		_msrMakePath(msr->param.achDataSubPath,STOPFILE,achTmp);
		_msrMakePath(msr->lcl.pszDataPath,achTmp,achFile);
		first_call = 0;
		}
	if ((fp = fopen(achFile,"r"))) {
		(void) printf("User interrupt detected.\n");
		(void) fclose(fp);
		(void) unlink(achFile);
		return 1;
		}
	return 0;
	}

void msrFinish(MSR msr)
{
	int id;

	for (id=1;id<msr->nThreads;++id) {
		if (msr->param.bVDetails) printf("Stopping thread %d\n",id);		
		mdlReqService(msr->mdl,id,SRV_STOP,NULL,0);
		mdlGetReply(msr->mdl,id,NULL,NULL);
		}
	pstFinish(msr->pst);
	/*
	 ** finish with parameter stuff, deallocate and exit.
	 */
	prmFinish(msr->prm);
	free(msr->pMap);
	free(msr);
	}

void msrOneNodeReadTipsy(MSR msr, struct inReadTipsy *in)
{
    int i,id;
    int *nParts;				/* number of particles for each processor */
    int nStart;
    PST pst0;
    LCL *plcl;
    char achInFile[PST_FILENAME_SIZE];
    int nid;
    int inswap;
    struct inSetParticleTypes intype;

    nParts = malloc(msr->nThreads*sizeof(*nParts));
    for (id=0;id<msr->nThreads;++id) {
		nParts[id] = -1;
		}

    pstOneNodeReadInit(msr->pst, in, sizeof(*in), nParts, &nid);
    assert(nid == msr->nThreads*sizeof(*nParts));
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
	_msrMakePath(plcl->pszDataPath,in->achInFile,achInFile);

    nStart = nParts[0];
	assert(msr->pMap[0] == 0);
    for (i=1;i<msr->nThreads;++i) {
		id = msr->pMap[i];
		/* 
		 * Read particles into the local storage.
		 */
		assert(plcl->pkd->nStore >= nParts[id]);
		pkdReadTipsy(plcl->pkd,achInFile,nStart,nParts[id],
					 in->bStandard,in->iReadIOrder,in->dvFac,in->dTuFac);
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
    pkdReadTipsy(plcl->pkd,achInFile,0,nParts[0],in->bStandard,in->iReadIOrder,in->dvFac,
				 in->dTuFac);

/* I think this code can be removed -- this is done later  JW Sept 2002 */
/*
	if (msr->param.iReadIOrder) {
		struct outGetNParts outget;
		struct inSetNParts inset;
		
		pstGetNParts(msr->pst,NULL,0,&outget,NULL);
		assert(outget.nGas == msr->nGas);
		assert(outget.nDark == msr->nDark);
		assert(outget.nStar == msr->nStar);
		inset.nGas = outget.nGas;
		inset.nDark = outget.nDark;
		inset.nStar = outget.nStar;
		msr->nMaxOrderGas = inset.nMaxOrderGas = outget.iMaxOrderGas;
		msr->nMaxOrderDark = inset.nMaxOrderDark = outget.iMaxOrderDark;
        msr->nMaxOrder = inset.nMaxOrder     = outget.iMaxOrderStar;
		pstSetNParts(msr->pst,&inset,sizeof(inset),NULL,NULL);
		}

    intype.nSuperCool = msr->param.nSuperCool;
    pstSetParticleTypes(msr->pst,&intype,sizeof(intype),NULL,NULL);
*/
    }

int xdrHeader(XDR *pxdrs,struct dump *ph)
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


double msrReadTipsy(MSR msr)
{
	FILE *fp;
	struct dump h;
	struct inReadTipsy in;
	char achInFile[PST_FILENAME_SIZE];
	LCL *plcl = msr->pst->plcl;
	double dTime,aTo,tTo,z;
	struct inSetParticleTypes intype;
	
	if (msr->param.achInFile[0]) {
		/*
		 ** Add Data Subpath for local and non-local names.
		 */
		_msrMakePath(msr->param.achDataSubPath,msr->param.achInFile,in.achInFile);
		/*
		 ** Add local Data Path.
		 */
		_msrMakePath(plcl->pszDataPath,in.achInFile,achInFile);

		fp = fopen(achInFile,"r");
		if (!fp) {
			printf("Could not open InFile:%s\n",achInFile);
			_msrExit(msr,1);
			}
		}
	else {
		printf("No input file specified\n");
		_msrExit(msr,1);
		return -1.0;
		}
	/*
	 ** Assume tipsy format for now, and dark matter only.
	 */
	if (msr->param.bStandard) {
		XDR xdrs;

		xdrstdio_create(&xdrs,fp,XDR_DECODE);
		xdrHeader(&xdrs,&h);
		xdr_destroy(&xdrs);
		}
	else {
		fread(&h,sizeof(struct dump),1,fp);
		}
	fclose(fp);

	msr->N = h.nbodies;
	msr->nDark = h.ndark;
	msr->nGas = h.nsph;
	msr->nStar = h.nstar;
	msr->nMaxOrder = msr->N - 1;
	msr->nMaxOrderGas = msr->nGas - 1;
	msr->nMaxOrderDark = msr->nGas + msr->nDark - 1;

	assert(msr->N == msr->nDark+msr->nGas+msr->nStar);
#ifndef GASOLINE
	if (msr->nGas != 0) fprintf(stderr,"GASOLINE compile flag not set:  Treating %d Gas particles as Dark\n",msr->nGas);
#endif
	if (msrComove(msr)) {
		if(msr->param.csm->dHubble0 == 0.0) {
			printf("No hubble constant specified\n");
			_msrExit(msr,1);
			}
		dTime = csmExp2Time(msr->param.csm,h.time);
		z = 1.0/h.time - 1.0;
		if (msr->param.bVStart)
			printf("Input file, Time:%g Redshift:%g Expansion factor:%g\n",
				   dTime,z,h.time);
		if (prmSpecified(msr->prm,"dRedTo")) {
			if (!prmArgSpecified(msr->prm,"nSteps") &&
				prmArgSpecified(msr->prm,"dDelta")) {
				aTo = 1.0/(msr->param.dRedTo + 1.0);
				tTo = csmExp2Time(msr->param.csm,aTo);
				if (msr->param.bVStart)
					printf("Simulation to Time:%g Redshift:%g Expansion factor:%g\n",
						   tTo,1.0/aTo-1.0,aTo);
				if (tTo < dTime) {
					printf("Badly specified final redshift, check -zto parameter.\n");
					_msrExit(msr,1);
					}
				msr->param.nSteps = (int)ceil((tTo-dTime)/msr->param.dDelta);
				}
			else if (!prmArgSpecified(msr->prm,"dDelta") &&
					 prmArgSpecified(msr->prm,"nSteps")) {
				aTo = 1.0/(msr->param.dRedTo + 1.0);
				tTo = csmExp2Time(msr->param.csm,aTo);
				if (msr->param.bVStart)
					printf("Simulation to Time:%g Redshift:%g Expansion factor:%g\n",
						   tTo,1.0/aTo-1.0,aTo);
				if (tTo < dTime) {
					printf("Badly specified final redshift, check -zto parameter.\n");	
					_msrExit(msr,1);
					}
				if(msr->param.nSteps != 0)
				    msr->param.dDelta =
					(tTo-dTime)/(msr->param.nSteps -
						     msr->param.iStartStep);
				
				else
				    msr->param.dDelta = 0.0;
				}
			else if (!prmSpecified(msr->prm,"nSteps") &&
					 prmFileSpecified(msr->prm,"dDelta")) {
				aTo = 1.0/(msr->param.dRedTo + 1.0);
				tTo = csmExp2Time(msr->param.csm,aTo);
				if (msr->param.bVStart)
					printf("Simulation to Time:%g Redshift:%g Expansion factor:%g\n",
						   tTo,1.0/aTo-1.0,aTo);
				if (tTo < dTime) {
					printf("Badly specified final redshift, check -zto parameter.\n");
					_msrExit(msr,1);
					}
				msr->param.nSteps = (int)ceil((tTo-dTime)/msr->param.dDelta);
				}
			else if (!prmSpecified(msr->prm,"dDelta") &&
					 prmFileSpecified(msr->prm,"nSteps")) {
				aTo = 1.0/(msr->param.dRedTo + 1.0);
				tTo = csmExp2Time(msr->param.csm,aTo);
				if (msr->param.bVStart)
					printf("Simulation to Time:%g Redshift:%g Expansion factor:%g\n",
						   tTo,1.0/aTo-1.0,aTo);
				if (tTo < dTime) {
					printf("Badly specified final redshift, check -zto parameter.\n");	
					_msrExit(msr,1);
					}
				if(msr->param.nSteps != 0)
				    msr->param.dDelta =	(tTo-dTime)/(msr->param.nSteps
													 - msr->param.iStartStep);
				else
				    msr->param.dDelta = 0.0;
				}
			}
		else {
			tTo = dTime + msr->param.nSteps*msr->param.dDelta;
			aTo = csmTime2Exp(msr->param.csm,tTo);
			if (msr->param.bVStart)
				printf("Simulation to Time:%g Redshift:%g Expansion factor:%g\n",
					   tTo,1.0/aTo-1.0,aTo);
			}
		if (msr->param.bVStart)
			printf("Reading file...\nN:%d nDark:%d nGas:%d nStar:%d\n",msr->N,
				   msr->nDark,msr->nGas,msr->nStar);
		if (msr->param.bCannonical) {
			in.dvFac = h.time*h.time;
			}
		else {
			in.dvFac = 1.0;
			}
		}
	else {
		dTime = h.time;
		if (msr->param.bVStart) printf("Input file, Time:%g\n",dTime);
		tTo = dTime + msr->param.nSteps*msr->param.dDelta;
		if (msr->param.bVStart) {
			printf("Simulation to Time:%g\n",tTo);
			printf("Reading file...\nN:%d nDark:%d nGas:%d nStar:%d Time:%g\n",
				   msr->N,msr->nDark,msr->nGas,msr->nStar,dTime);
			}
		in.dvFac = 1.0;
		}
	in.nFileStart = 0;
	in.nFileEnd = msr->N - 1;
	in.nDark = msr->nDark;
	in.nGas = msr->nGas;
	in.nStar = msr->nStar;
	in.iOrder = msr->param.iOrder;
	in.bStandard = msr->param.bStandard;
	in.iReadIOrder = msr->param.iReadIOrder;
#ifdef GASOLINE
	in.dTuFac = msr->param.dGasConst/(msr->param.dConstGamma - 1)/
		msr->param.dMeanMolWeight;
#else
	in.dTuFac = 1.0;
#endif
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
	in.fPeriod[0] = msr->param.dxPeriod;
	in.fPeriod[1] = msr->param.dyPeriod;
	in.fPeriod[2] = msr->param.dzPeriod;

	if(msr->param.bParaRead)
	    pstReadTipsy(msr->pst,&in,sizeof(in),NULL,NULL);
	else
	    msrOneNodeReadTipsy(msr, &in);

	if (msr->param.iReadIOrder) {
		struct outGetNParts outget;
		struct inSetNParts inset;
		
		pstGetNParts(msr->pst,NULL,0,&outget,NULL);
		assert(outget.nGas == msr->nGas);
		assert(outget.nDark == msr->nDark);
		assert(outget.nStar == msr->nStar);
		inset.nGas = outget.nGas;
		inset.nDark = outget.nDark;
		inset.nStar = outget.nStar;
		msr->nMaxOrderGas = inset.nMaxOrderGas = outget.iMaxOrderGas;
		msr->nMaxOrderDark = inset.nMaxOrderDark = outget.iMaxOrderDark;
        msr->nMaxOrder = inset.nMaxOrder     = outget.iMaxOrderStar;
		pstSetNParts(msr->pst,&inset,sizeof(inset),NULL,NULL);
		if (msr->param.bVDetails) puts("IOrder file has been successfully read.");
		}

	intype.nSuperCool = msr->param.nSuperCool;
	pstSetParticleTypes(msr->pst, &intype, sizeof(intype), NULL, NULL);
	if (msr->param.bVDetails) puts("Input file has been successfully read.");
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


/*
 ** This function makes some DANGEROUS assumptions!!!
 ** Main problem is that it calls pkd level routines, bypassing the
 ** pst level. It uses plcl pointer which is not desirable.
 */
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
	_msrMakePath(plcl->pszDataPath,in->achOutFile,achOutFile);

    /* 
     * First write our own particles.
     */
    pkdWriteTipsy(plcl->pkd,achOutFile,plcl->nWriteStart,in->bStandard,
				  in->dvFac,in->duTFac,in->iGasModel); 
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
					  in->bStandard, in->dvFac, in->duTFac,in->iGasModel); 
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


void msrCalcWriteStart(MSR msr) 
{
	struct outSetTotal out;
	struct inSetWriteStart in;

	pstSetTotal(msr->pst,NULL,0,&out,NULL);
	assert(out.nTotal == msr->N);
	in.nWriteStart = 0;
	pstSetWriteStart(msr->pst,&in,sizeof(in),NULL,NULL);
	}


void msrWriteTipsy(MSR msr,char *pszFileName,double dTime)
{
	FILE *fp;
	struct dump h;
	struct inWriteTipsy in;
	char achOutFile[PST_FILENAME_SIZE];
	LCL *plcl = msr->pst->plcl;

	/*
	 ** Calculate where to start writing.
	 ** This sets plcl->nWriteStart.
	 */
	msrCalcWriteStart(msr);
	/*
	 ** Add Data Subpath for local and non-local names.
	 */
	_msrMakePath(msr->param.achDataSubPath,pszFileName,in.achOutFile);
	/*
	 ** Add local Data Path.
	 */
	_msrMakePath(plcl->pszDataPath,in.achOutFile,achOutFile);
	
	fp = fopen(achOutFile,"w");
	if (!fp) {
		printf("Could not open OutFile:%s\n",achOutFile);
		_msrExit(msr,1);
		}
	in.bStandard = msr->param.bStandard;
#ifdef GASOLINE
	in.duTFac = (msr->param.dConstGamma - 1)*msr->param.dMeanMolWeight/
		msr->param.dGasConst;
#else
	in.duTFac = 1.0;
#endif
	in.iGasModel = msr->param.iGasModel;
	/*
	 ** Assume tipsy format for now.
	 */
	h.nbodies = msr->N;
	h.ndark = msr->nDark;
	h.nsph = msr->nGas;
	h.nstar = msr->nStar;
	if (msrComove(msr)) {
		h.time = csmTime2Exp(msr->param.csm,dTime);
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
	if (msr->param.bVDetails) {
		if (msrComove(msr)) {
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
		xdrHeader(&xdrs,&h);
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
	if (msr->param.bVDetails)
		puts("Output file has been successfully written.");
	}


void msrSetSoft(MSR msr,double dSoft)
{
	struct inSetSoft in;
  
	if (msr->param.bVDetails) printf("Set Softening...\n");
	in.dSoft = dSoft;
	pstSetSoft(msr->pst,&in,sizeof(in),NULL,NULL);
	}


void msrDomainDecomp(MSR msr, int iRung, int bGreater)
{
	struct inDomainDecomp in;
	int iRungDD,iRungSD,nActive;

#ifdef GASOLINE
	/* Sanity check on Gas particles being present */
	if (msr->nGas==0 && (msr->param.bDoGas==1 || msr->param.bGasDomainDecomp)) {
		if (msr->param.bGasDomainDecomp) {
			printf("MDD: switching bGasDomainDecomp off.\n");
			msr->param.bGasDomainDecomp=0;
			msr->bDoneDomainDecomp=0;
			}
		if (msr->param.bDoGas==1) {
			printf("MDD: switching bDoGas off.\n");
			msr->param.bDoGas=0;
			}
		}
#endif

	in.bDoRootFind = 1;
	in.bDoSplitDimFind = 1;
	
	nActive=0;
	if (bGreater) {
		iRungDD=msr->iCurrMaxRung+1; 
		while (iRungDD > iRung) {
			iRungDD--;
			nActive+=msr->nRung[iRungDD];
			}
		while(iRungDD > 0 && nActive < msr->N*msr->param.dFracNoDomainDecomp) {
			iRungDD--;
			nActive+=msr->nRung[iRungDD];
			}
		iRungSD = iRungDD;
		while(iRungSD > 0 && nActive < msr->N*msr->param.dFracNoDomainDimChoice) {
			iRungSD--;
			nActive+=msr->nRung[iRungSD];
			}
		}
	else {
		iRungDD = iRung;
		while(iRungDD > 0 && msr->nRung[iRungDD] < msr->N*msr->param.dFracNoDomainDecomp) {
			iRungDD--;
			}
		iRungSD = iRungDD;
		while(iRungSD > 0 && msr->nRung[iRungSD] < msr->N*msr->param.dFracNoDomainDimChoice) {
			iRungSD--;
			}
		}

	if (msr->nActive < msr->N*msr->param.dFracNoDomainDecomp) {
		if (msr->bDoneDomainDecomp && msr->iLastRungDomainDecomp >= iRungDD) {
			if (msr->param.bVRungStat) printf("Skipping Root Finder (nActive = %d/%d, iRung %d/%d/%d)\n",msr->nActive,msr->N,iRung,iRungDD,msr->iLastRungDomainDecomp);
			in.bDoRootFind = 0;
			in.bDoSplitDimFind = 0;
			}
		else if (iRungDD < iRung) {
			/* Set up the DD for the highest rung that still gets one */
			msrActiveRung(msr,iRungDD,bGreater);
			}
		}
	else iRungDD = iRung;

	if (in.bDoRootFind && msr->bDoneDomainDecomp && iRungDD > iRungSD && msr->iLastRungDomainDecomp >= iRungSD) {
		if (msr->param.bVRungStat) printf("Skipping Split Dim Finding (nDDActive = %d/%d, iRung %d/%d/%d/%d)\n",msr->nActive,msr->N,iRung,iRungDD,iRungSD,msr->iLastRungDomainDecomp);
		in.bDoSplitDimFind = 0;
		}

	if (msr->param.bGasDomainDecomp) {
		pstGasWeight(msr->pst,NULL,0,NULL,NULL);
		}

	if (msr->param.bRungDD) {
	        struct inRungDDWeight inRDD;
		inRDD.iMaxRung = msr->iCurrMaxRung;
                inRDD.dWeight = msr->param.dRungDDWeight;
		pstRungDDWeight(msr->pst,&inRDD,sizeof(struct inRungDDWeight),NULL,NULL);
		}

	if (msr->param.bVDetails) {
		double sec,dsec;
		printf("nActive %d nTreeActive %d nSmoothActive %d\n",msr->nActive,msr->nTreeActive,msr->nSmoothActive);
		printf("Domain Decomposition...\n");
		sec = msrTime();

		pstDomainDecomp(msr->pst,&in,sizeof(in),NULL,NULL);
		msr->bDoneDomainDecomp = 1; 
		dsec = msrTime() - sec;
		printf("Domain Decomposition complete, Wallclock: %f secs\n\n",dsec);
		}
	else {
		pstDomainDecomp(msr->pst,&in,sizeof(in),NULL,NULL);
		msr->bDoneDomainDecomp = 1; 
		}

	msr->iLastRungDomainDecomp = iRungDD;
	if (iRungDD < iRung) {
	        /* Restore Active data */
		msrActiveRung(msr,iRung,bGreater);
		}
	}

void msrBuildTree(MSR msr,int bActiveOnly, double dMass,int bSmooth)
{
	struct inBuildTree in;
	struct outBuildTree out;
	struct inColCells inc;
	struct ioCalcRoot root;
	KDN *pkdn;
	int iDum,nCell;

	if (msr->param.bVDetails) printf("Building local trees...\n");

	/*
	 ** First make sure the particles are in Active/Inactive order.
	 */
	msrActiveTypeOrder(msr, TYPE_ACTIVE|TYPE_TREEACTIVE );
	in.nBucket = msr->param.nBucket;
	in.iOpenType = msr->iOpenType;
	in.iOrder = (msr->param.iOrder >= msr->param.iEwOrder)?
		msr->param.iOrder:msr->param.iEwOrder;
	in.dCrit = msr->dCrit;
	if (bSmooth) {
		in.bBinary = 0;
		in.bGravity = 0;
		msr->iTreeType = MSR_TREE_DENSITY;
		msr->bGravityTree = 0;
		}
	else {
		in.bBinary = msr->param.bBinary;
		in.bGravity = 1;
		msr->bGravityTree = 1;
		if (msr->param.bBinary) {
			msr->iTreeType = MSR_TREE_SPATIAL;
			}
		else {
			msr->iTreeType = MSR_TREE_DENSITY;
			}
		}
	in.bActiveOnly = bActiveOnly;
	in.bTreeActiveOnly = bActiveOnly;
	if (msr->param.bVDetails) {
		double sec,dsec;
		sec = msrTime();
		pstBuildTree(msr->pst,&in,sizeof(in),&out,&iDum);
		msrMassCheck(msr,dMass,"After pstBuildTree in msrBuildTree");
		dsec = msrTime() - sec;
		printf("Tree built, Wallclock: %f secs\n\n",dsec);
		}
	else {
		pstBuildTree(msr->pst,&in,sizeof(in),&out,&iDum);
		msrMassCheck(msr,dMass,"After pstBuildTree in msrBuildTree");
		}
	nCell = 1<<(1+(int)ceil(log((double)msr->nThreads)/log(2.0)));
	pkdn = malloc(nCell*sizeof(KDN));
	assert(pkdn != NULL);
	inc.iCell = ROOT;
	inc.nCell = nCell;
	pstColCells(msr->pst,&inc,sizeof(inc),pkdn,NULL);
	msrMassCheck(msr,dMass,"After pstColCells in msrBuildTree");
#if (0)
	{ int i;
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
	}
#endif
	pstDistribCells(msr->pst,pkdn,nCell*sizeof(KDN),NULL,NULL);
	msrMassCheck(msr,dMass,"After pstDistribCells in msrBuildTree");
	free(pkdn);
	if (!bSmooth) {
		pstCalcRoot(msr->pst,NULL,0,&root,&iDum);
		msrMassCheck(msr,dMass,"After pstCalcRoot in msrBuildTree");
		pstDistribRoot(msr->pst,&root,sizeof(struct ioCalcRoot),NULL,NULL);
		msrMassCheck(msr,dMass,"After pstDistribRoot in msrBuildTree");
	    }
    }


void msrDomainColor(MSR msr)
{
	pstDomainColor(msr->pst,NULL,0,NULL,NULL);
	}


void msrReorder(MSR msr)
{
	struct inDomainOrder in;

	in.iMaxOrder = msrMaxOrder(msr);
	if (msr->param.bVDetails) {
		double sec,dsec;
		printf("Ordering...\n");
		sec = msrTime();
		pstDomainOrder(msr->pst,&in,sizeof(in),NULL,NULL);
		pstLocalOrder(msr->pst,NULL,0,NULL,NULL);
		dsec = msrTime() - sec;
		printf("Order established, Wallclock: %f secs\n\n",dsec);
		}
	else {
		pstDomainOrder(msr->pst,&in,sizeof(in),NULL,NULL);
		pstLocalOrder(msr->pst,NULL,0,NULL,NULL);
		}
 	/*
	 ** Mark tree as void.
	 */
	msr->iTreeType = MSR_TREE_NONE;
 	/*
	 ** Mark domain decomp as not done.
	 */
	msr->bDoneDomainDecomp = 0;
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
		_msrMakePath(msr->param.achDataSubPath,pszFile,in.achOutFile);
		/*
		 ** Add local Data Path.
		 */
		_msrMakePath(plcl->pszDataPath,in.achOutFile,achOutFile);

		fp = fopen(achOutFile,"w");
		if (!fp) {
			printf("Could not open Array Output File:%s\n",achOutFile);
			_msrExit(msr,1);
			}
		}
	else {
		printf("No Array Output File specified\n");
		_msrExit(msr,1);
		return;
		}
	/*
	 ** Write the Header information and close the file again.
	 */
	if (msr->param.iBinaryOutput) {
		fwrite(&msr->N,sizeof(int),1,fp);
		}
	else {
		fprintf(fp,"%d\n",msr->N);
		}
	fclose(fp);
	/* 
	 * First write our own particles.
	 */
	assert(msr->pMap[0] == 0);
	pkdOutArray(plcl->pkd,achOutFile,iType,msr->param.iBinaryOutput); 
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
#ifdef SIMPLESF
	    pkdOutArray(plcl->pkd,achOutFile,iType,((iType==OUT_TCOOLAGAIN_ARRAY || iType==OUT_MSTAR_ARRAY) ? 1 : msr->param.iBinaryOutput));
#else
	    pkdOutArray(plcl->pkd,achOutFile,iType,msr->param.iBinaryOutput); 
#endif
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
		_msrMakePath(msr->param.achDataSubPath,pszFile,in.achOutFile);
		/*
		 ** Add local Data Path.
		 */
		_msrMakePath(plcl->pszDataPath,in.achOutFile,achOutFile);

		fp = fopen(achOutFile,"w");
		if (!fp) {
			printf("Could not open Vector Output File:%s\n",achOutFile);
			_msrExit(msr,1);
			}
		}
	else {
		printf("No Vector Output File specified\n");
		_msrExit(msr,1);
		return;
		}
	/*
	 ** Write the Header information and close the file again.
	 */
	if (msr->param.iBinaryOutput) {
		fwrite(&msr->N,sizeof(int),1,fp);
		}
	else {
		fprintf(fp,"%d\n",msr->N);
		}
	fclose(fp);

	/* 
	 * First write our own particles.
	 */
	assert(msr->pMap[0] == 0);
	
	if (msr->param.bPackedVector) {
		pkdOutVector(plcl->pkd,achOutFile,-3,iType,msr->param.iBinaryOutput); 
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
			pkdOutVector(plcl->pkd,achOutFile,-3,iType,msr->param.iBinaryOutput); 
			/* 
			 * Swap them back again.
			 */
			inswap = 0;
			mdlReqService(pst0->mdl,id,PST_SWAPALL,&inswap,sizeof(inswap));
			pkdSwapAll(plcl->pkd,id);
			mdlGetReply(pst0->mdl,id,NULL,NULL);
			}
		}
	else {
		for (iDim=0;iDim<3;++iDim) {
			pkdOutVector(plcl->pkd,achOutFile,iDim,iType,msr->param.iBinaryOutput); 
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
				pkdOutVector(plcl->pkd,achOutFile,iDim,iType,msr->param.iBinaryOutput); 
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
	}


void msrSmooth(MSR msr,double dTime,int iSmoothType,int bSymmetric)
{
	struct inSmooth in;
	struct outSmooth out;

	/*
	 ** Make sure that the type of tree is a density binary tree!
	 */
	assert(msr->iTreeType == MSR_TREE_DENSITY);
	in.nSmooth = msr->param.nSmooth;
	in.bPeriodic = msr->param.bPeriodic;
	in.bSymmetric = bSymmetric;
	in.iSmoothType = iSmoothType;
	in.dfBall2OverSoft2 = (msr->param.bLowerSoundSpeed ? 0 :
			       4.0*msr->param.dhMinOverSoft*msr->param.dhMinOverSoft);
	if (msrComove(msr)) {
		in.smf.H = csmTime2Hub(msr->param.csm,dTime);
		in.smf.a = csmTime2Exp(msr->param.csm,dTime);
		}
	else {
		in.smf.H = 0.0;
		in.smf.a = 1.0;
		}
	{
	    double dAccFac = 1.0/(in.smf.a*in.smf.a*in.smf.a);
	    in.smf.dDeltaAccelFac = msr->param.dEtaDeltaAccel/sqrt(dAccFac);
	    }
	in.smf.bSinkThermal = msr->param.bSinkThermal;
	in.smf.dSinkRadius = msr->param.dSinkRadius;
	in.smf.dSinkBoundOrbitRadius = msr->param.dSinkBoundOrbitRadius;
#ifdef GASOLINE
	in.smf.alpha = msr->param.dConstAlpha;
	in.smf.beta = msr->param.dConstBeta;
	in.smf.gamma = msr->param.dConstGamma;
	in.smf.algam = in.smf.alpha*sqrt(in.smf.gamma*(in.smf.gamma - 1));
	in.smf.bGeometric = msr->param.bGeometric;
	in.smf.bCannonical = msr->param.bCannonical;
	in.smf.bGrowSmoothList = 0;
#endif
#ifdef STARFORM
        in.smf.dMinMassFrac = msr->param.stfm->dMinMassFrac;
	in.smf.dTime = dTime;
	in.smf.bSmallSNSmooth = msr->param.bSmallSNSmooth;
	in.smf.bShortCoolShutoff = msr->param.bShortCoolShutoff;
        /* from McKee and Ostriker (1977) ApJ 218 148 */
        in.smf.dRadPreFactor = pow(10,1.74)/(msr->param.dKpcUnit*1000.0)*
    pow(MSOLG*msr->param.dMsolUnit/(msr->param.dMeanMolWeight*MHYDR*pow(KPCCM*msr->param.dKpcUnit,3)),-0.16)*
    pow(0.0001*GCGS*pow(MSOLG*msr->param.dMsolUnit,2)/(pow(KPCCM*msr->param.dKpcUnit,4)*KBOLTZ),-0.2);
	if (msr->param.bShortCoolShutoff){        /* end of snowplow */
        	in.smf.dTimePreFactor = SECONDSPERYEAR*pow(10,5.92)/(msr->param.dSecUnit)*
    pow(MSOLG*msr->param.dMsolUnit/(msr->param.dMeanMolWeight*MHYDR*pow(KPCCM*msr->param.dKpcUnit,3)),0.27)*
    pow(0.0001*GCGS*pow(MSOLG*msr->param.dMsolUnit,2)/(pow(KPCCM*msr->param.dKpcUnit,4)*KBOLTZ),-0.64);
		} else {       /* t_{max}*/
        	in.smf.dTimePreFactor = SECONDSPERYEAR*pow(10,6.85)/(msr->param.dSecUnit)*
    pow(MSOLG*msr->param.dMsolUnit/(msr->param.dMeanMolWeight*MHYDR*pow(KPCCM*msr->param.dKpcUnit,3)),0.32)*
    pow(0.0001*GCGS*pow(MSOLG*msr->param.dMsolUnit,2)/(pow(KPCCM*msr->param.dKpcUnit,4)*KBOLTZ),-0.70);
		}
#endif /*STARFORM*/
#ifdef COLLISIONS
	in.smf.dCentMass = msr->param.dCentMass; /* for Hill sphere checks */
#endif
#ifdef SLIDING_PATCH /* called by msrFindRejects() only */
	in.smf.dOrbFreq = msr->param.dOrbFreq;
	in.smf.dTime = dTime;
	in.smf.dCentMass = 0; /* to disable Hill sphere checks */
#endif
	if (msr->param.bVStep) {
		double sec,dsec;
		printf("Smoothing...\n");
		sec = msrTime();
		pstSmooth(msr->pst,&in,sizeof(in),&out,NULL);
		dsec = msrTime() - sec;
		printf("Smooth Calculated, Wallclock: %f secs\n\n",dsec);
		if (msr->nThreads > 1) {
		    double iP = 1.0/msr->nThreads;
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
	else {
		pstSmooth(msr->pst,&in,sizeof(in),&out,NULL);
		}
	}


void msrReSmooth(MSR msr,double dTime,int iSmoothType,int bSymmetric)
{
	struct inReSmooth in;

	/*
	 ** Make sure that the type of tree is a density binary tree!
	 */
	assert(msr->iTreeType == MSR_TREE_DENSITY);
	in.nSmooth = msr->param.nSmooth;
	in.bPeriodic = msr->param.bPeriodic;
	in.bSymmetric = bSymmetric;
	in.iSmoothType = iSmoothType;
	in.dfBall2OverSoft2 = (msr->param.bLowerSoundSpeed ? 0 :
			       4.0*msr->param.dhMinOverSoft*msr->param.dhMinOverSoft);
	if (msrComove(msr)) {
		in.smf.H = csmTime2Hub(msr->param.csm,dTime);
		in.smf.a = csmTime2Exp(msr->param.csm,dTime);
		}
	else {
		in.smf.H = 0.0;
		in.smf.a = 1.0;
		}
	{
	    double dAccFac = 1.0/(in.smf.a*in.smf.a*in.smf.a);
	    in.smf.dDeltaAccelFac = msr->param.dEtaDeltaAccel/sqrt(dAccFac);
	    }
	in.smf.dSinkRadius = msr->param.dSinkRadius;
	in.smf.dSinkBoundOrbitRadius = msr->param.dSinkBoundOrbitRadius;
#ifdef GASOLINE
	in.smf.alpha = msr->param.dConstAlpha;
	in.smf.beta = msr->param.dConstBeta;
	in.smf.gamma = msr->param.dConstGamma;
	in.smf.algam = in.smf.alpha*sqrt(in.smf.gamma*(in.smf.gamma - 1));
	in.smf.bGeometric = msr->param.bGeometric;
	in.smf.bCannonical = msr->param.bCannonical;
	in.smf.bGrowSmoothList = 0;
#endif
	if (msr->param.bVStep) {
		double sec,dsec;
		struct outSmooth out;

		printf("ReSmoothing...\n");
		sec = msrTime();
		pstReSmooth(msr->pst,&in,sizeof(in),&out,NULL);
		dsec = msrTime() - sec;
		printf("ReSmooth Calculated, Wallclock: %f secs\n\n",dsec);
		if (msr->nThreads > 1) {
		    double iP = 1.0/msr->nThreads;
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
	else {
		pstReSmooth(msr->pst,&in,sizeof(in),NULL,NULL);
		}
	}

void msrMarkSmooth(MSR msr,double dTime,int bSymmetric,int iMarkType)
{
	struct inMarkSmooth in;

	/*
	 ** Make sure that the type of tree is a density binary tree!
	 */
	assert(msr->iTreeType == MSR_TREE_DENSITY);
	in.nSmooth = msr->param.nSmooth;
	in.bPeriodic = msr->param.bPeriodic;
	in.bSymmetric = bSymmetric;
	in.iSmoothType = SMX_MARK; 
	in.iMarkType = iMarkType;
	if (msrComove(msr)) {
		in.smf.H = csmTime2Hub(msr->param.csm,dTime);
		in.smf.a = csmTime2Exp(msr->param.csm,dTime);
		}
	else {
		in.smf.H = 0.0;
		in.smf.a = 1.0;
		}
#ifdef GASOLINE
	in.smf.alpha = msr->param.dConstAlpha;
	in.smf.beta = msr->param.dConstBeta;
	in.smf.gamma = msr->param.dConstGamma;
	in.smf.algam = in.smf.alpha*sqrt(in.smf.gamma*(in.smf.gamma - 1));
	in.smf.bGeometric = msr->param.bGeometric;
	in.smf.bCannonical = msr->param.bCannonical;
	in.smf.bGrowSmoothList = 0;
#endif
	if (msr->param.bVStep) {
		double sec,dsec;
		printf("MarkSmoothing...\n");
		sec = msrTime();
		pstMarkSmooth(msr->pst,&in,sizeof(in),NULL,NULL);
		dsec = msrTime() - sec;
		printf("MarkSmooth Calculated, Wallclock: %f secs\n\n",dsec);
		}
	else {
		pstMarkSmooth(msr->pst,&in,sizeof(in),NULL,NULL);
		}
	}

void msrUpdateSoft(MSR msr,double dTime) {
#ifdef CHANGESOFT
       if (!(msr->param.bPhysicalSoft || msr->param.bVariableSoft)) return;
       if (msr->param.bPhysicalSoft) {
	 struct inPhysicalSoft in;

	 in.dFac = 1./csmTime2Exp(msr->param.csm,dTime);
	 in.bSoftMaxMul = msr->param.bSoftMaxMul;
	 in.dSoftMax = msr->param.dSoftMax;

	 if (msr->param.bSoftMaxMul && in.dFac > in.dSoftMax) in.dFac = in.dSoftMax;

	 pstPhysicalSoft(msr->pst,&in,sizeof(in),NULL,NULL);
       }
       else {
	 struct inPostVariableSoft inPost;

	 pstPreVariableSoft(msr->pst,NULL,0,NULL,NULL);

	 if (msr->param.bSoftByType) {
	   if (msr->nDark) {
	     msrActiveType(msr,TYPE_DARK,TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
	     msrBuildTree(msr,1,-1.0,1);
	     msrSmooth(msr,dTime,SMX_NULL,0);
	   }
	   if (msr->nGas) {
	     msrActiveType(msr,TYPE_GAS,TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
	     msrBuildTree(msr,1,-1.0,1);
	     msrSmooth(msr,dTime,SMX_NULL,0);
	   }
	   if (msr->nStar) {
	     msrActiveType(msr,TYPE_STAR,TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
	     msrBuildTree(msr,1,-1.0,1);
	     msrSmooth(msr,dTime,SMX_NULL,0);
	   }
	 }
	 else {
	   msrActiveType(msr,TYPE_ALL,TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
	   msrBuildTree(msr,1,-1.0,1);
	   msrSmooth(msr,dTime,SMX_NULL,0);
	 }

	 inPost.dSoftMax = msr->param.dSoftMax;
	 inPost.bSoftMaxMul = msr->param.bSoftMaxMul;
	 pstPostVariableSoft(msr->pst,&inPost,sizeof(inPost),NULL,NULL);
       }
#endif
}


void msrGravity(MSR msr,double dStep,int bDoSun,
				int *piSec,double *pdWMax,double *pdIMax,double *pdEMax,
				int *pnActive)
{
	struct inGravity in;
	struct outGravity out;
	struct inGravExternal inExt;
	int iDum,j;
	double sec,dsec;

	if (msr->param.bDoSelfGravity) {
		assert(msr->bGravityTree == 1);
		assert(msr->iTreeType == MSR_TREE_SPATIAL || 
			   msr->iTreeType == MSR_TREE_DENSITY);
		if (msr->param.bVStep) printf("Calculating Gravity, Step:%f\n",dStep);
		in.nReps = msr->param.nReplicas;
		in.bPeriodic = msr->param.bPeriodic;
		in.iOrder = msr->param.iOrder;
		in.bEwald = msr->param.bEwald;
		in.iEwOrder = msr->param.iEwOrder;
#ifdef SLIDING_PATCH
		in.dOrbFreq = msr->param.dOrbFreq;
		in.dTime = dStep*msr->param.dDelta;
#endif
		/*
		 ** The meaning of 'bDoSun' here is that we want the accel on (0,0,0)
		 ** to use in creating the indirect acceleration on each particle. This
		 ** is why it looks inconsistent with call to pstGravExternal() below.
		 */
		in.bDoSun = msr->param.bHeliocentric;
		in.dEwCut = msr->param.dEwCut;
		in.dEwhCut = msr->param.dEwhCut;
		in.dSunSoft = msr->param.dSunSoft;
		sec = msrTime();
		pstGravity(msr->pst,&in,sizeof(in),&out,&iDum);
		dsec = msrTime() - sec;
#ifdef SPECIAL_PARTICLES
		if (msr->param.nSpecial) {
			/*
			 ** Handle contributions from "special" particles now,
			 ** e.g. oblateness effects, GR, etc.
			 */
			struct inGetSpecial inGetSpec;
			struct outGetSpecial outGetSpec;
			struct inDoSpecial inDoSpec;
			struct outDoSpecial outDoSpec;
			inGetSpec.nSpecial = msr->param.nSpecial;
			for (j=0;j<msr->param.nSpecial;j++)
				inGetSpec.iId[j] = msr->param.iSpecialId[j]; /* orig indices */
			inGetSpec.dCentMass = msr->param.dCentMass; /* if special frame */
			pstGetSpecialParticles(msr->pst,&inGetSpec,sizeof(inGetSpec),
								   &outGetSpec,NULL);
			inDoSpec.nSpecial = msr->param.nSpecial;
			for (j=0;j<msr->param.nSpecial;j++) {
				inDoSpec.sData[j] = msr->param.sSpecialData[j]; /* struct cp */
				inDoSpec.sInfo[j] = outGetSpec.sInfo[j]; /* ditto */
				}
			inDoSpec.bNonInertial = msr->param.bHeliocentric;
			pstDoSpecialParticles(msr->pst,&inDoSpec,sizeof(inDoSpec),
								  &outDoSpec,NULL);
			for (j=0;j<3;j++)
				out.aSun[j] += outDoSpec.aFrame[j]; /* add non-inertial term */
			}
#endif
		}
	else { /* pstGravity() not called, so zero all output parameters */
		dsec = out.nActive = out.nTreeActive = out.aSun[0] =
		out.aSun[1] = out.aSun[2] = out.dPartSum = out.dCellSum =
		out.dSoftSum = out.dFlop = out.dWSum = out.dWMax = out.dWMin =
		out.dISum = out.dIMax = out.dIMin = out.dESum = out.dEMax =
		out.dEMin = out.dpASum = out.dpMSum = out.dpCSum = out.dpTSum
		= out.dcASum = out.dcMSum = out.dcCSum = out.dcTSum = 0;
		}
	/* enforced initialization */
#ifdef ROT_FRAME
	inExt.bRotFrame = 0;
#endif
#ifdef SLIDING_PATCH
	inExt.bPatch = 0;
#endif
#ifdef SIMPLE_GAS_DRAG
	inExt.bSimpleGasDrag = 0;
#endif
	/*
	 ** Calculate any external potential stuff.
	 ** This may contain a huge list of flags in the future, so we may want
	 ** to replace this test with something like bAnyExternal.
	 */
	if (msr->param.bHeliocentric || msr->param.bLogHalo ||
		msr->param.bHernquistSpheroid || msr->param.bNFWSpheroid ||
		msr->param.bElliptical ||
		msr->param.bHomogSpheroid || msr->param.bBodyForce ||
        msr->param.bMiyamotoDisk || msr->param.bTimeVarying) {
		/*
		 ** Provide the time.
		 */
		inExt.dTime = dStep*msr->param.dDelta;
		/*
		 ** Only allow inclusion of solar terms if we are in Heliocentric 
		 ** coordinates.
		 */
		inExt.bIndirect = msr->param.bHeliocentric;
		inExt.bDoSun = bDoSun;  /* Treat the Sun explicitly. */
		inExt.dSunMass = msr->param.dCentMass;
		inExt.dSunSoft = msr->param.dSunSoft;
		if(inExt.bIndirect)
		    for (j=0;j<3;++j) inExt.aSun[j] = out.aSun[j];
		inExt.bLogHalo = msr->param.bLogHalo;
		inExt.bHernquistSpheroid = msr->param.bHernquistSpheroid;
		inExt.bNFWSpheroid = msr->param.bNFWSpheroid;
		inExt.bElliptical= msr->param.bElliptical;
		inExt.bEllipticalDarkNFW= msr->param.bEllipticalDarkNFW;
		inExt.bHomogSpheroid = msr->param.bHomogSpheroid;
		inExt.bBodyForce = msr->param.bBodyForce;
		inExt.bMiyamotoDisk = msr->param.bMiyamotoDisk;
		inExt.bTimeVarying = msr->param.bTimeVarying;
		pstGravExternal(msr->pst,&inExt,sizeof(inExt),NULL,NULL);
		}
#ifdef ROT_FRAME
	if (msr->param.bRotFrame) { /* general rotating frame */
		inExt.bIndirect = inExt.bDoSun = inExt.bLogHalo =
			inExt.bNFWSpheroid = inExt.bHernquistSpheroid =
				inExt.bMiyamotoDisk = 0;
		inExt.bRotFrame = msr->param.bRotFrame;
		inExt.dOmega = msr->param.dOmega +
			dStep*msr->param.dDelta*msr->param.dOmegaDot;
		inExt.dOmegaDot = msr->param.dOmegaDot;
		pstGravExternal(msr->pst,&inExt,sizeof(inExt),NULL,NULL);
		}
#endif
#if defined(SLIDING_PATCH) || defined(SIMPLE_GAS_DRAG)
	if (msr->param.bPatch || msr->param.bSimpleGasDrag) {
		inExt.bIndirect = inExt.bDoSun = inExt.bLogHalo =
			inExt.bNFWSpheroid = inExt.bHernquistSpheroid =
				inExt.bBodyForce =
				inExt.bMiyamotoDisk = 0;
#ifdef SLIDING_PATCH
		if (msr->param.bPatch) {
			static int bFirstCall = 1;
			static double dOrbFreqZ2 = 0;
			inExt.bPatch = msr->param.bPatch;
			inExt.dOrbFreq = msr->param.dOrbFreq;
			if (bFirstCall) {
				dOrbFreqZ2 = msr->param.dOrbFreq*msr->param.dOrbFreq;
				if (msr->param.bDoSelfGravity)
					dOrbFreqZ2 += 2*M_PI*msrMassCheck(msr,-1.0,"")/msr->param.nReplicas/
						pow(msr->param.dxPeriod*msr->param.dyPeriod,1.5);
				bFirstCall = 0;
				if (msr->param.bVStart)
					(void) printf("Vert. freq. enhancement factor = %g\n",
								  sqrt(dOrbFreqZ2)/msr->param.dOrbFreq);
				}
			inExt.dOrbFreqZ2 = dOrbFreqZ2;
			}
#endif
#ifdef SIMPLE_GAS_DRAG
		if (msr->param.bSimpleGasDrag) {
			inExt.bSimpleGasDrag = msr->param.bSimpleGasDrag;
			inExt.iFlowOpt	= 1; /* temporary */
			inExt.bEpstein	= msr->param.bEpstein;
			inExt.dGamma	= msr->param.dGamma;
			inExt.dOmegaZ	= msr->param.dOrbFreq;
			}
#endif
		pstGravExternal(msr->pst,&inExt,sizeof(inExt),NULL,NULL);
		}
#endif /* if (SLIDING_PATCH || SIMPLE_GAS_DRAG) */
#ifdef AGGS
	msrAggsGravity(msr);
#endif
	/*
	 ** Output.
	 */
	*piSec = dsec;
	*pnActive = out.nActive;
	*pdWMax = out.dWMax;
	*pdIMax = out.dIMax;
	*pdEMax = out.dEMax;
	if (msr->param.bVStep) {
		double dPartAvg,dCellAvg,dSoftAvg;
		double dWAvg,dWMax,dWMin;
		double dIAvg,dIMax,dIMin;
		double dEAvg,dEMax,dEMin;
		double iP;

		/*
		 ** Output some info...
		 */
		if (dsec > 0.0) {
			double dMFlops = out.dFlop/dsec*1e-6;
			printf("Gravity Calculated, Wallclock: %f secs, MFlops:%.1f, Flop:%.3g\n",
				   dsec,dMFlops,out.dFlop);
			}
		else {
			printf("Gravity Calculated, Wallclock: %f secs, MFlops:unknown, Flop:%.3g\n",
				   dsec,out.dFlop);
			}
		if (out.nActive > 0) {
			dPartAvg = out.dPartSum/out.nActive;
			dCellAvg = out.dCellSum/out.nActive;
			dSoftAvg = out.dSoftSum/out.nActive;
			}
		else {
			dPartAvg = dCellAvg = dSoftAvg = 0;
#ifdef COLLISIONS
			/*
			 ** This is allowed to happen in the collision model because a
			 ** time-step rung may be left vacant following a merger event.
			 */
#else
			if (msr->param.bVWarnings)
				printf("WARNING: no particles found!\n");
#endif
			}
		iP = 1.0/msr->nThreads;
		dWAvg = out.dWSum*iP;
		dIAvg = out.dISum*iP;
		dEAvg = out.dESum*iP;
		dWMax = out.dWMax;
		dIMax = out.dIMax;
		dEMax = out.dEMax;
		dWMin = out.dWMin;
		dIMin = out.dIMin;
		dEMin = out.dEMin;
		printf("dPartAvg:%f dCellAvg:%f dSoftAvg:%f\n",
			   dPartAvg,dCellAvg,dSoftAvg);
		printf("Walk CPU     Avg:%10f Max:%10f Min:%10f\n",dWAvg,dWMax,dWMin);
		printf("Interact CPU Avg:%10f Max:%10f Min:%10f\n",dIAvg,dIMax,dIMin);
		if (msr->param.bEwald) printf("Ewald CPU    Avg:%10f Max:%10f Min:%10f\n",dEAvg,dEMax,dEMin);
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
	}


void msrCalcEandL(MSR msr,int bFirst,double dTime,double *E,double *T,
				  double *U,double *Eth,double L[])
{
	struct outCalcEandL out;
	struct inCalcEandLExt inExt;
	struct outCalcEandLExt outExt;
	double a;
	int k;

	pstCalcEandL(msr->pst,NULL,0,&out,NULL);
	*T = out.T;
	*U = out.U;
	*Eth = out.Eth;
	for (k=0;k<3;k++) L[k] = out.L[k];
	/*
	 ** Calculate E & L contribution from external potential and/or
	 ** reference frame. Currently only heliocentric frame implemented.
	 ** NOTE: In some cases the correction terms may be comparable to
	 ** the original values, resulting in precision errors.
	 */
	if (msr->param.bHeliocentric) {
		double dSunPos[3],dSunVel[3],dRxV[3],dTotMass,dSunVel2;
		inExt.bHeliocentric = msr->param.bHeliocentric;
		pstCalcEandLExt(msr->pst,&inExt,sizeof(inExt),&outExt,NULL);
		dTotMass = outExt.dMass + msr->param.dCentMass;
		assert(dTotMass > 0);
		dSunVel2 = 0;
		for (k=0;k<3;k++) {
			dSunPos[k] = - outExt.dSumMR[k]/dTotMass;
			dSunVel[k] = - outExt.dSumMV[k]/dTotMass;
			dSunVel2 += dSunVel[k]*dSunVel[k];
			}
		dRxV[0] = dSunPos[1]*dSunVel[2] - dSunPos[2]*dSunVel[1];
		dRxV[1] = dSunPos[2]*dSunVel[0] - dSunPos[0]*dSunVel[2];
		dRxV[2] = dSunPos[0]*dSunVel[1] - dSunPos[1]*dSunVel[0];
		*T -= 0.5*dTotMass*dSunVel2;
		*U -= 0.5*msr->param.dCentMass*outExt.dPot;
		for (k=0;k<3;k++) L[k] -= dTotMass*dRxV[k];
		}
	/*
	 ** Do the comoving coordinates stuff.
	 ** Currently L is not adjusted for this. Should it be?
	 */
	a = csmTime2Exp(msr->param.csm,dTime);
	if (!msr->param.bCannonical) *T *= pow(a,4.0);
	/*
	 * Estimate integral (\dot a*U*dt) over the interval.
	 * Note that this is equal to integral (W*da) and the latter
	 * is more accurate when a is changing rapidly.
	 */
	if (msrComove(msr) && !bFirst) {
		msr->dEcosmo += 0.5*(a - csmTime2Exp(msr->param.csm, msr->dTimeOld))
			*((*U) + msr->dUOld);
		}
	else {
		msr->dEcosmo = 0.0;
		}
	msr->dTimeOld = dTime;
	msr->dUOld = *U;
	*U *= a;
#ifdef COLLISIONS
	if (bFirst)
		msr->dTcoll = 0;
	else
		*T -= msr->dTcoll;
#endif
	*E = (*T) + (*U) - msr->dEcosmo + a*a*(*Eth);
	}


void msrDrift(MSR msr,double dTime,double dDelta)
{
	struct inDrift in;
	int j;
#ifdef NEED_VPRED
	struct inKickVpred invpr;
	double a;
#endif

	/*
	 ** This only gets done if growmass parameters have been set!
	 */
	msrGrowMass(msr,dTime,dDelta);

#ifdef AGGS
	msrAggsAdvanceOpen(msr);
#endif

#ifdef COLLISIONS
	msrDoCollisions(msr,dTime,dDelta);
#endif

	if (msr->param.bCannonical) {
		in.dDelta = csmComoveDriftFac(msr->param.csm,dTime,dDelta);
		}
	else {
		in.dDelta = dDelta;
		}
	for (j=0;j<3;++j) {
		in.fCenter[j] = msr->fCenter[j];
		}
	in.bPeriodic = msr->param.bPeriodic;
	in.bFandG = msr->param.bFandG;
	in.fCentMass = msr->param.dCentMass;
#ifdef SLIDING_PATCH
	in.dOrbFreq = msr->param.dOrbFreq;
	in.dTime = dTime;
#endif
	pstDrift(msr->pst,&in,sizeof(in),NULL,NULL);
	/*
	 ** Once we move the particles the tree should no longer be considered 
	 ** valid.
	 */
	msr->iTreeType = MSR_TREE_NONE;

#ifdef NEED_VPRED

#ifdef PREDRHO
	if (msr->param.bPredRho == 2) {
		struct inKickRhopred inRhop;
		if (msrComove(msr)) 
			inRhop.dHubbFac = 3*csmTime2Hub(msr->param.csm, dTime + dDelta/2.0);
		else
			inRhop.dHubbFac = 0.0;
		inRhop.dDelta = dDelta;
		/* Non Active Gas particles need to have updated densities */
		pstKickRhopred(msr->pst,&inRhop,sizeof(inRhop),NULL,NULL);
		}
#endif

	if (msr->param.bCannonical) {
#ifdef GLASS	  
		invpr.dvFacOne = 1.0 - fabs(dDelta)*msr->param.dGlassDamper; /* Damp velocities */
#else
		invpr.dvFacOne = 1.0; /* no hubble drag, man! */
#endif
		invpr.dvFacTwo = csmComoveKickFac(msr->param.csm,dTime,dDelta);
		invpr.duDelta  = dDelta;
		invpr.iGasModel = msr->param.iGasModel;
		dTime += dDelta/2.0;
		a = csmTime2Exp(msr->param.csm,dTime);
		invpr.z = 1/a - 1;
		invpr.duDotLimit = msr->param.duDotLimit;
		}
	else {
		double H;
		/*
		 ** Careful! For non-cannonical we want H and a at the 
		 ** HALF-STEP! This is a bit messy but has to be special
		 ** cased in some way.
		 */
		dTime += dDelta/2.0;
		a = csmTime2Exp(msr->param.csm,dTime);
		H = csmTime2Hub(msr->param.csm,dTime);
		invpr.dvFacOne = (1.0 - H*dDelta)/(1.0 + H*dDelta);
		invpr.dvFacTwo = dDelta/pow(a,3.0)/(1.0 + H*dDelta);
		invpr.duDelta  = dDelta;
		invpr.iGasModel = msr->param.iGasModel;
		invpr.z = 1/a - 1;
		invpr.duDotLimit = msr->param.duDotLimit;
		}
	if (dDelta != 0.0) {
		struct outKick out;
		int nout;
	    pstKickVpred(msr->pst,&invpr,sizeof(invpr),&out,&nout);
		if (nout) printf("Drift (Vpred): Avg Wallclock %f, Max Wallclock %f\n",
						 out.SumTime/out.nSum,out.MaxTime);
	    }
#endif /* NEED_VPRED */

#ifdef AGGS
	msrAggsAdvanceClose(msr,dDelta);
#endif
	}

/* Untested */
void msrCoolVelocity(MSR msr,double dTime,double dMass)
{
#ifdef SUPERCOOL
	struct inCoolVelocity in;
	
	if (msr->param.nSuperCool > 0) {
		/*
		 ** Activate all for densities if bSymCool == 0
		 */
		if (msr->param.bSymCool) {
			msrActiveType(msr,TYPE_SUPERCOOL,TYPE_ACTIVE);
			msrActiveType(msr,TYPE_ALL,TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
			msrDomainDecomp(msr,0,1);
			msrActiveType(msr,TYPE_SUPERCOOL,TYPE_ACTIVE);
			msrBuildTree(msr,1,dMass,1);
			msrSmooth(msr,dTime,SMX_DENSITY,1);
			msrReSmooth(msr,dTime,SMX_MEANVEL,1);
			}
		else {
			/*
			 ** Note, here we need to calculate the densities of all
			 ** the particles so that the mean velocity can be 
			 ** calculated.
			 */
			/* activate all */
			msrActiveTypeRung(msr,TYPE_SUPERCOOL,TYPE_ACTIVE,0,1); 
			msrActiveType(msr,TYPE_ALL,TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
			msrDomainDecomp(msr,0,1);
			msrActiveTypeRung(msr,TYPE_SUPERCOOL,TYPE_ACTIVE,0,1); 
			msrBuildTree(msr,0,SMX_DENSITY,1);
			msrSmooth(msr,dTime,SMX_DENSITY,0);
			msrReSmooth(msr,dTime,SMX_MEANVEL,0);
			}
		/*
		 ** Now cool them.
		 */
		in.nSuperCool = msr->param.nSuperCool;
		in.dCoolFac = msr->param.dCoolFac;
		in.dCoolDens = msr->param.dCoolDens;
		in.dCoolMaxDens = msr->param.dCoolMaxDens;
		pstCoolVelocity(msr->pst,&in,sizeof(in),NULL,NULL);
		}
#endif
	}

void msrGrowMass(MSR msr, double dTime, double dDelta)
{
    struct inGrowMass in;
    
    if (msr->param.nGrowMass > 0 && dTime > msr->param.dGrowStartT &&
		dTime <= msr->param.dGrowEndT) {
		in.nGrowMass = msr->param.nGrowMass;
		in.dDeltaM = msr->param.dGrowDeltaM*dDelta/
			(msr->param.dGrowEndT - msr->param.dGrowStartT);
		pstGrowMass(msr->pst, &in, sizeof(in), NULL, NULL);
		}
    }

/*
 * For gasoline, updates predicted velocities to middle of timestep.
 */
void msrKickDKD(MSR msr,double dTime,double dDelta)
{
	double H,a;
	struct inKick in;
	struct outKick out;
	
	if (msr->param.bCannonical) {
#ifdef GLASS	  
		in.dvFacOne = 1.0 - fabs(dDelta)*msr->param.dGlassDamper; /* Damp velocities */
#else
		in.dvFacOne = 1.0; /* no hubble drag, man! */
#endif
#ifdef NEED_VPRED
#ifdef GLASS	  
		in.dvPredFacOne = 1.0 - fabs(dDelta)*msr->param.dGlassDamper; /* Damp velocities */
#else
		in.dvPredFacOne = 1.0;
#endif
		in.dvPredFacTwo = csmComoveKickFac(msr->param.csm,dTime,0.5*dDelta);
		in.duDelta      = dDelta;
		in.duPredDelta  = 0.5*dDelta;
		in.iGasModel = msr->param.iGasModel;
		dTime += dDelta/2.0;
		a = csmTime2Exp(msr->param.csm,dTime);
		in.z = 1/a - 1;
		in.duDotLimit = msr->param.duDotLimit;
#endif /* NEED_VPRED */
		}
	else {
		/*
		 ** Careful! For non-cannonical we want H and a at the 
		 ** HALF-STEP! This is a bit messy but has to be special
		 ** cased in some way.
		 */
		dTime += dDelta/2.0;
		a = csmTime2Exp(msr->param.csm,dTime);
		H = csmTime2Hub(msr->param.csm,dTime);
		in.dvFacOne = (1.0 - H*dDelta)/(1.0 + H*dDelta);
		in.dvFacTwo = dDelta/pow(a,3.0)/(1.0 + H*dDelta);
#ifdef NEED_VPRED
		dTime -= dDelta/4.0;
		a = csmTime2Exp(msr->param.csm,dTime);
		H = csmTime2Hub(msr->param.csm,dTime);
		in.dvPredFacOne = (1.0 - 0.5*H*dDelta)/(1.0 + 0.5*H*dDelta);
		in.dvPredFacTwo = 0.5*dDelta/pow(a,3.0)/(1.0 + 0.5*H*dDelta);
		in.duDelta      = dDelta;
		in.duPredDelta  = 0.5*dDelta;
		in.iGasModel = msr->param.iGasModel;
		in.z = 1/a - 1;
		in.duDotLimit = msr->param.duDotLimit;
#endif /* NEED_VPRED */
		}
	pstKick(msr->pst,&in,sizeof(in),&out,NULL);
	printf("Kick: Avg Wallclock %f, Max Wallclock %f\n",
	       out.SumTime/out.nSum,out.MaxTime);
	}

/*
 * For gasoline, updates predicted velocities to beginning of timestep.
 */
void msrKickKDKOpen(MSR msr,double dTime,double dDelta)
{
	/* NOTE: dDelta should be 1/2 the drift step here... */

	double H,a;
	struct inKick in;
	struct outKick out;
	
	if (msr->param.bCannonical) {
#ifdef GLASS	  
		in.dvFacOne = 1.0 - fabs(dDelta)*msr->param.dGlassDamper;  /* Damp velocities */
#else
		in.dvFacOne = 1.0;		/* no hubble drag, man! */
#endif
		in.dvFacTwo = csmComoveKickFac(msr->param.csm,dTime,dDelta);
#ifdef NEED_VPRED
#ifdef GLASS	  
		in.dvPredFacOne = 1.0 - fabs(dDelta)*msr->param.dGlassDamper;  /* Damp velocities */
#else
		in.dvPredFacOne = 1.0;
#endif
		in.dvPredFacTwo = 0.0;
		in.duDelta      = dDelta;
		in.duPredDelta  = 0.0;
		in.iGasModel = msr->param.iGasModel;
		dTime += dDelta/2.0;
		a = csmTime2Exp(msr->param.csm,dTime);
		in.z = 1/a - 1;
		in.duDotLimit = msr->param.duDotLimit;
#endif /* NEED_VPRED */
		}
	else {
		/*
		 ** Careful! For non-cannonical we want H and a at the 
		 ** HALF-STEP! This is a bit messy but has to be special
		 ** cased in some way.
		 */
		dTime += dDelta/2.0;
		a = csmTime2Exp(msr->param.csm,dTime);
		H = csmTime2Hub(msr->param.csm,dTime);
		in.dvFacOne = (1.0 - H*dDelta)/(1.0 + H*dDelta);
		in.dvFacTwo = dDelta/pow(a,3.0)/(1.0 + H*dDelta);
#ifdef NEED_VPRED
		in.dvPredFacOne = 1.0;
		in.dvPredFacTwo = 0.0;
		in.duDelta      = dDelta;
		in.duPredDelta  = 0.0;
		in.iGasModel = msr->param.iGasModel;
		in.z = 1/a - 1;
		in.duDotLimit = msr->param.duDotLimit;
#endif /* NEED_VPRED */
		}
#ifdef AGGS
	msrAggsDeactivate(msr);
#endif
	pstKick(msr->pst,&in,sizeof(in),&out,NULL);
#ifdef AGGS
	msrAggsActivate(msr);
#endif
	if (msr->param.bVDetails)
		printf("KickOpen: Avg Wallclock %f, Max Wallclock %f\n",
			   out.SumTime/out.nSum,out.MaxTime);
#ifdef COLLISIONS
	/*
	 ** This would be better as an external potential call, but usually
	 ** we want uniform gravity *without* interparticle gravity, in which
	 ** case pstGravExternal() is never called.
	 */
	{
	struct inKickUnifGrav inkug;
	inkug.dvx = msr->param.dxUnifGrav*dDelta;
	inkug.dvy = msr->param.dyUnifGrav*dDelta;
	inkug.dvz = msr->param.dzUnifGrav*dDelta;
	pstKickUnifGrav(msr->pst,&inkug,sizeof(inkug),NULL,NULL);
	}
#endif
#ifdef AGGS
	msrAggsKick(msr,dDelta);
#endif
	}

/*
 * For gasoline, updates predicted velocities to end of timestep.
 */
void msrKickKDKClose(MSR msr,double dTime,double dDelta)
{
	/* NOTE: dDelta should be 1/2 the drift step here... */

	double H,a;
	struct inKick in;
	struct outKick out;
	
	if (msr->param.bCannonical) {
#ifdef GLASS	  
		in.dvFacOne = 1.0 - fabs(dDelta)*msr->param.dGlassDamper; /* Damp velocities */
#else
		in.dvFacOne = 1.0; /* no hubble drag, man! */
#endif
		in.dvFacTwo = csmComoveKickFac(msr->param.csm,dTime,dDelta);
#ifdef NEED_VPRED
#ifdef GLASS	  
		in.dvPredFacOne = 1.0 - fabs(dDelta)*msr->param.dGlassDamper; /* Damp velocities */
#else
		in.dvPredFacOne = 1.0;
#endif
		in.dvPredFacTwo = in.dvFacTwo;
		in.duDelta      = dDelta;
		in.duPredDelta  = dDelta;
		in.iGasModel = msr->param.iGasModel;
		dTime += dDelta/2.0;
		a = csmTime2Exp(msr->param.csm,dTime);
		in.z = 1/a - 1;
		in.duDotLimit = msr->param.duDotLimit;
#endif /* NEED_VPRED */
		}
	else {
		/*
		 ** Careful! For non-cannonical we want H and a at the 
		 ** HALF-STEP! This is a bit messy but has to be special
		 ** cased in some way.
		 */
		dTime += dDelta/2.0;
		a = csmTime2Exp(msr->param.csm,dTime);
		H = csmTime2Hub(msr->param.csm,dTime);
		in.dvFacOne = (1.0 - H*dDelta)/(1.0 + H*dDelta);
		in.dvFacTwo = dDelta/pow(a,3.0)/(1.0 + H*dDelta);
#ifdef NEED_VPRED
		in.dvPredFacOne = in.dvFacOne;
		in.dvPredFacTwo = in.dvFacTwo;
		in.duDelta      = dDelta;
		in.duPredDelta  = dDelta;
		in.iGasModel = msr->param.iGasModel;
		in.z = 1/a - 1;
		in.duDotLimit = msr->param.duDotLimit;
#endif /* NEED_VPRED */
		}
#ifdef AGGS
	msrAggsDeactivate(msr);
#endif
	pstKick(msr->pst,&in,sizeof(in),&out,NULL);
#ifdef AGGS
	msrAggsActivate(msr);
#endif
	if (msr->param.bVDetails)
		printf("KickClose: Avg Wallclock %f, Max Wallclock %f\n",
			   out.SumTime/out.nSum,out.MaxTime);
#ifdef COLLISIONS
	{
	struct inKickUnifGrav inkug;
	inkug.dvx = msr->param.dxUnifGrav*dDelta;
	inkug.dvy = msr->param.dyUnifGrav*dDelta;
	inkug.dvz = msr->param.dzUnifGrav*dDelta;
	pstKickUnifGrav(msr->pst,&inkug,sizeof(inkug),NULL,NULL);
	}
#endif
#ifdef AGGS
	msrAggsKick(msr,dDelta);
#endif
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

    tin.nFileStart = in->nFileStart;
    tin.nFileEnd = in->nFileEnd;
    tin.iOrder = in->iOrder;
    tin.fExtraStore = in->fExtraStore;
    tin.nDark = in->nDark;
    tin.nGas = in->nGas;
    tin.nStar = in->nStar;
    tin.bStandard = 0;
	tin.iReadIOrder = 0;
    tin.dvFac = 1;
    tin.dTuFac = 1;

    for(j = 0; j < 3; j++)
		tin.fPeriod[j] = in->fPeriod[j];

    pstOneNodeReadInit(msr->pst, &tin, sizeof(tin), nParts, &nid);
    assert(nid == msr->nThreads*sizeof(*nParts));
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
	_msrMakePath(plcl->pszDataPath,in->achInFile,achInFile);

    nStart = nParts[0];
    assert(msr->pMap[0] == 0);
    for (i=1;i<msr->nThreads;++i) {
		id = msr->pMap[i];
		/* 
		 * Read particles into the local storage.
		 */
		assert(plcl->pkd->nStore >= nParts[id]);
		pkdReadCheck(plcl->pkd,achInFile,in->iVersion,
					 in->iOffset, nStart, nParts[id]);
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
    pkdReadCheck(plcl->pkd,achInFile,in->iVersion,in->iOffset,0,nParts[0]);
    }

int msrFindCheck(MSR msr) {
	char achCheckFile[PST_FILENAME_SIZE];
	char achLatestCheckFile[PST_FILENAME_SIZE];
	time_t latestTime;
	char extraString[PST_FILENAME_SIZE];
	int iNotCorrupt;
	FDL_CTX* fdl;
	struct stat statbuf;
	char* suffixes[3] = {"chk0", "chk1", "chk"};
	int i;
	
	latestTime = 0;
	
	for(i = 0; i < 3; ++i) {
		/* Get the filename of the checkpoint file formatted properly. */
		sprintf(achCheckFile, "%s/%s.%s", msr->param.achDataSubPath, msr->param.achOutName, suffixes[i]);
		if(msr->pst->plcl->pszDataPath) {
			strcpy(extraString, achCheckFile);
			sprintf(achCheckFile,"%s/%s", msr->pst->plcl->pszDataPath, extraString);
		}

		/* Check if the checkpoint file exists. */
		if(!stat(achCheckFile, &statbuf)) {
			/* Check if the checkpoint file is corrupt. */
			fdl = FDL_open(achCheckFile);
			/* Read and ignore the version number. */
			FDL_read(fdl, "version" ,&iNotCorrupt);
			FDL_read(fdl, "not_corrupt_flag", &iNotCorrupt);
			if(iNotCorrupt == 1) {
				/* Get the file creation time. */
				if(statbuf.st_mtime > latestTime) {
					/* The checkpoint file exists, isn't corrupt, and is more recent, so mark it. */
					latestTime = statbuf.st_mtime;
					strcpy(achLatestCheckFile, achCheckFile);
				}
			}
			FDL_finish(fdl);
		}
	}
	
	/* Did we find any checkpoint files? */
	if(latestTime == 0)
		return 0;
	
	/* Rename latest valid checkpoint file to .chk, return success. */
	if(!rename(achLatestCheckFile, achCheckFile)) /* The last filename checked is the one we want. */
		return 1;
	else
		return 0;
}

double msrReadCheck(MSR msr,int *piStep)
{
	struct inReadCheck in;
	struct inSetNParts inset;
	struct inSetParticleTypes intype;
	char achInFile[PST_FILENAME_SIZE];
	int i;
	LCL *plcl = msr->pst->plcl;
	double dTime;
	int iVersion,iNotCorrupt;
	FDL_CTX *fdl;
	int nMaxOutsTmp,nOutsTmp;
	double *pdOutTimeTmp;

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
	if (msr->param.bVDetails)
		printf("Reading Version-%d Checkpoint file.\n",iVersion);
	FDL_read(fdl,"not_corrupt_flag",&iNotCorrupt);
	if (iNotCorrupt != 1) {
		printf("Sorry the checkpoint file is corrupted.\n");
		_msrExit(msr,1);
		}
	FDL_read(fdl,"number_of_particles",&msr->N);
	/*
	 ** As of checkpoint version 3 we include numbers of dark gas and star 
	 ** particles to support GASOLINE.
	 */
	if (iVersion > 2) {
		FDL_read(fdl,"number_of_gas_particles",&msr->nGas);
		FDL_read(fdl,"number_of_dark_particles",&msr->nDark);
		FDL_read(fdl,"number_of_star_particles",&msr->nStar);
		assert(msr->N == msr->nDark+msr->nGas+msr->nStar);
		}
	else {
		msr->nDark = msr->N;
		msr->nGas = 0;
		msr->nStar = 0;
		}
	FDL_read(fdl,"current_timestep",piStep);
	FDL_read(fdl,"current_time",&dTime);
	FDL_read(fdl,"current_ecosmo",&msr->dEcosmo);
	FDL_read(fdl,"old_time",&msr->dTimeOld);
	FDL_read(fdl,"old_potentiale",&msr->dUOld);
#ifdef COLLISIONS
	msr->dTcoll = 0; /*DEBUG*/
#endif
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
	FDL_read(fdl,"bComove",&msr->param.csm->bComove);
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
	/*
	 ** bBinary somehow got left out of the previous versions of the
	 ** checkpoint file. We fix this as of checkpoint version 5.
	 */
	if (iVersion > 4) {
		if (!prmSpecified(msr->prm,"bBinary")) 
			FDL_read(fdl,"bBinary",&msr->param.bBinary);
		}
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
	if (iVersion > 4) {
	    if (!prmSpecified(msr->prm,"dEtaCourant"))
		FDL_read(fdl,"dEtaCourant",&msr->param.dEtaCourant);
	    }
	if (!prmSpecified(msr->prm,"bEpsAccStep"))
		FDL_read(fdl,"bEpsAccStep",&msr->param.bEpsAccStep);
	if (!prmSpecified(msr->prm,"bNonSymp"))
		FDL_read(fdl,"bNonSymp",&msr->param.bNonSymp);
	if (!prmSpecified(msr->prm,"iMaxRung"))
		FDL_read(fdl,"iMaxRung",&msr->param.iMaxRung);
	if (!prmSpecified(msr->prm,"dEwCut"))
		FDL_read(fdl,"dEwCut",&msr->param.dEwCut);
	if (!prmSpecified(msr->prm,"dEwhCut"))
		FDL_read(fdl,"dEwhCut",&msr->param.dEwhCut);
	FDL_read(fdl,"dPeriod",&msr->param.dPeriod);
	if (iVersion > 3) {
		FDL_read(fdl,"dxPeriod",&msr->param.dxPeriod);
		FDL_read(fdl,"dyPeriod",&msr->param.dyPeriod);
		FDL_read(fdl,"dzPeriod",&msr->param.dzPeriod);
		}
	else {
		msr->param.dxPeriod = msr->param.dPeriod;
		msr->param.dyPeriod = msr->param.dPeriod;
		msr->param.dzPeriod = msr->param.dPeriod;
		}
	FDL_read(fdl,"dHubble0",&msr->param.csm->dHubble0);
	FDL_read(fdl,"dOmega0",&msr->param.csm->dOmega0);
	if (iVersion > 4) {
	    FDL_read(fdl,"dLambda",&msr->param.csm->dLambda);
	    FDL_read(fdl,"dOmegaRad",&msr->param.csm->dOmegaRad);
	    FDL_read(fdl,"dQuintess",&msr->param.csm->dQuintess);
	    }
	if (iVersion > 3) {
		if (!prmSpecified(msr->prm,"dTheta"))
			FDL_read(fdl,"dTheta",&msr->param.dTheta);
		if (!prmSpecified(msr->prm,"dTheta2"))
			FDL_read(fdl,"dTheta2",&msr->param.dTheta2);
		}
	else {
		if (!prmSpecified(msr->prm,"dTheta"))
			msr->param.dTheta = msr->dCrit;
		if (!prmSpecified(msr->prm,"dTheta2"))
			msr->param.dTheta2 = msr->dCrit;
		}
	if (iVersion > 5) {
		if (!prmSpecified(msr->prm,"bDoGravity"))
			FDL_read(fdl,"bDoGravity",&msr->param.bDoGravity);
		if (!prmSpecified(msr->prm,"bDoGas"))
			FDL_read(fdl,"bDoGas",&msr->param.bDoGas);
		if (!prmSpecified(msr->prm,"dEtaCourant"))
			FDL_read(fdl,"dEtaCourant",&msr->param.dEtaCourant);
		if (!prmSpecified(msr->prm,"iStartStep"))
			FDL_read(fdl,"iStartStep",&msr->param.iStartStep);
		if (!prmSpecified(msr->prm,"dFracNoDomainDecomp"))
			FDL_read(fdl,"dFracNoDomainDecomp",&msr->param.dFracNoDomainDecomp);
#ifndef NOCOOLING
#if defined(COOLING_COSMO) && defined(GASOLINE)
		if (!prmSpecified(msr->prm,"dMassFracHelium"))
			FDL_read(fdl,"dMassFracHelium",&msr->param.CoolParam.dMassFracHelium);
#endif
#endif
		}
	if (iVersion > 3) {
	    FDL_read(fdl, "max_order", &msr->nMaxOrder);
	    FDL_read(fdl, "max_order_gas", &msr->nMaxOrderGas);
	    FDL_read(fdl, "max_order_dark", &msr->nMaxOrderDark);
	    }
	else {
	    msr->nMaxOrder = msr->N - 1;
	    msr->nMaxOrderGas = -1;
	    msr->nMaxOrderDark = msr->nMaxOrder;
	    }
	if (iVersion > 4) {
		FDL_read(fdl,"bFandG",&msr->param.bFandG);
		FDL_read(fdl,"bHeliocentric",&msr->param.bHeliocentric);
		FDL_read(fdl,"dCentMass",&msr->param.dCentMass);
		FDL_read(fdl,"bLogHalo",&msr->param.bLogHalo);
		FDL_read(fdl,"bHernquistSpheroid",&msr->param.bHernquistSpheroid);
		FDL_read(fdl,"bMiyamotoDisk",&msr->param.bMiyamotoDisk);
		}
	else {
		msr->param.bFandG = 0;
		msr->param.bHeliocentric = 0;
		msr->param.dCentMass = 0.0;
		msr->param.bLogHalo = 0;
		msr->param.bHernquistSpheroid = 0;
		msr->param.bMiyamotoDisk = 0;
		}

	/*
	 * Check if redshift file is present, and if so reread it --JPG
	 */
	/* Store old values in temporary space */
	pdOutTimeTmp = malloc(msr->nMaxOuts*sizeof(double));
	nMaxOutsTmp = msr->nMaxOuts;
	nOutsTmp = msr->nOuts;
	for (i=0;i<msr->nOuts;++i) {
	    pdOutTimeTmp[i] = msr->pdOutTime[i];
	}
	/* Calculate initial time to give to msrReadOuts*/
	msrReadOuts(msr,dTime - (msr->param.nSteps*msr->param.dDelta));
	if (msr->nOuts == 0) { /* No redshift file...use old values */
	    free(msr->pdOutTime);
	    msr->pdOutTime = pdOutTimeTmp;
	    msr->nMaxOuts = nMaxOutsTmp;
	    msr->nOuts = nOutsTmp;
	}

	if (msr->param.bVDetails)
		printf("Reading checkpoint file...\nN:%d Time:%g\n",msr->N,dTime);
	in.nFileStart = 0;
	in.nFileEnd = msr->N - 1;
	in.nDark = msr->nDark;
	in.nGas = msr->nGas;
	in.nStar = msr->nStar;
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
	in.fPeriod[0] = msr->param.dxPeriod;
	in.fPeriod[1] = msr->param.dyPeriod;
	in.fPeriod[2] = msr->param.dzPeriod;
	in.iVersion = iVersion;
	in.iOffset = FDL_offset(fdl,"particle_array");
	FDL_finish(fdl);
	if(msr->param.bParaRead)
	    pstReadCheck(msr->pst,&in,sizeof(in),NULL,NULL);
	else
	    msrOneNodeReadCheck(msr,&in);
	if (msr->param.bVDetails)
		puts("Checkpoint file has been successfully read.");
	inset.nGas = msr->nGas;
	inset.nDark = msr->nDark;
	inset.nStar = msr->nStar;
	inset.nMaxOrderGas = msr->nMaxOrderGas;
	inset.nMaxOrderDark = msr->nMaxOrderDark;
	inset.nMaxOrder = msr->nMaxOrder;
	pstSetNParts(msr->pst,&inset,sizeof(inset),NULL,NULL);
	intype.nSuperCool = msr->param.nSuperCool;
	pstSetParticleTypes(msr->pst,&intype,sizeof(intype),NULL,NULL);
	
	i = msrCountType(msr, TYPE_GAS, TYPE_GAS);
	assert(i == msr->nGas);
	i = msrCountType(msr, TYPE_DARK, TYPE_DARK);
	assert(i == msr->nDark);
	i = msrCountType(msr, TYPE_STAR, TYPE_STAR);
	assert(i == msr->nStar);
	
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
	_msrMakePath(plcl->pszDataPath,in->achOutFile,achOutFile);

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
	char achOutFile[PST_FILENAME_SIZE],achTmp[PST_FILENAME_SIZE];
	int i;
	LCL *plcl = msr->pst->plcl;
	FDL_CTX *fdl;
	char *pszFdl;
	int iVersion,iNotCorrupt;
	static int first = 1;
	
	/*
	 ** Add Data Subpath for local and non-local names.
	 */
	if(first) {
	    sprintf(achTmp,"%s.chk0",msr->param.achOutName);
	    first = 0;
	} else {
		sprintf(achTmp,"%s.chk1",msr->param.achOutName);
		first = 1;
		}
	_msrMakePath(msr->param.achDataSubPath,achTmp,in.achOutFile);
	/*
	 ** Add local Data Path.
	 */
	_msrMakePath(plcl->pszDataPath,in.achOutFile,achOutFile);
	pszFdl = getenv("PKDGRAV_CHECKPOINT_FDL");
	if (pszFdl == NULL) {
          fprintf(stderr,"PKDGRAV_CHECKPOINT_FDL environment variable not set: no Checkpoint written\n");
          return;
        }	  
	fdl = FDL_create(achOutFile,pszFdl);
	iVersion = CHECKPOINT_VERSION;
	FDL_write(fdl,"version",&iVersion);
	iNotCorrupt = 0;
	FDL_write(fdl,"not_corrupt_flag",&iNotCorrupt);
	/*
	 ** Checkpoint header.
	 */
	FDL_write(fdl,"number_of_particles",&msr->N);
	FDL_write(fdl,"number_of_gas_particles",&msr->nGas);
	FDL_write(fdl,"number_of_dark_particles",&msr->nDark);
	FDL_write(fdl,"number_of_star_particles",&msr->nStar);
	FDL_write(fdl, "max_order", &msr->nMaxOrder);
	FDL_write(fdl, "max_order_gas", &msr->nMaxOrderGas);
	FDL_write(fdl, "max_order_dark", &msr->nMaxOrderDark);
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
	FDL_write(fdl,"bComove",&msr->param.csm->bComove);
	FDL_write(fdl,"bParaRead",&msr->param.bParaRead);
	FDL_write(fdl,"bParaWrite",&msr->param.bParaWrite);
	FDL_write(fdl,"bCannonical",&msr->param.bCannonical);
	FDL_write(fdl,"bStandard",&msr->param.bStandard);
	FDL_write(fdl,"bKDK",&msr->param.bKDK);
	FDL_write(fdl,"bBinary",&msr->param.bBinary);
	FDL_write(fdl,"bDoGravity",&msr->param.bDoGravity);
	FDL_write(fdl,"bDoGas",&msr->param.bDoGas);
	FDL_write(fdl,"bFandG",&msr->param.bFandG);
	FDL_write(fdl,"bHeliocentric",&msr->param.bHeliocentric);
	FDL_write(fdl,"bLogHalo",&msr->param.bLogHalo);
	FDL_write(fdl,"bHernquistSpheroid",&msr->param.bHernquistSpheroid);
	FDL_write(fdl,"bMiyamotoDisk",&msr->param.bMiyamotoDisk);
	FDL_write(fdl,"nBucket",&msr->param.nBucket);
	FDL_write(fdl,"iOutInterval",&msr->param.iOutInterval);
	FDL_write(fdl,"iLogInterval",&msr->param.iLogInterval);
	FDL_write(fdl,"iCheckInterval",&msr->param.iCheckInterval);
	FDL_write(fdl,"iExpOrder",&msr->param.iOrder);
	FDL_write(fdl,"iEwOrder",&msr->param.iEwOrder);
	FDL_write(fdl,"nReplicas",&msr->param.nReplicas);
	FDL_write(fdl,"iStartStep",&msr->param.iStartStep);
	FDL_write(fdl,"nSteps",&msr->param.nSteps);
	FDL_write(fdl,"dExtraStore",&msr->param.dExtraStore);
	FDL_write(fdl,"dDelta",&msr->param.dDelta);
	FDL_write(fdl,"dEta",&msr->param.dEta);
	FDL_write(fdl,"dEtaCourant",&msr->param.dEtaCourant);
	FDL_write(fdl,"bEpsAccStep",&msr->param.bEpsAccStep);
	FDL_write(fdl,"bNonSymp",&msr->param.bNonSymp);
	FDL_write(fdl,"iMaxRung",&msr->param.iMaxRung);
	FDL_write(fdl,"dEwCut",&msr->param.dEwCut);
	FDL_write(fdl,"dEwhCut",&msr->param.dEwhCut);
	FDL_write(fdl,"dPeriod",&msr->param.dPeriod);
	FDL_write(fdl,"dxPeriod",&msr->param.dxPeriod);
	FDL_write(fdl,"dyPeriod",&msr->param.dyPeriod);
	FDL_write(fdl,"dzPeriod",&msr->param.dzPeriod);
	FDL_write(fdl,"dHubble0",&msr->param.csm->dHubble0);
	FDL_write(fdl,"dOmega0",&msr->param.csm->dOmega0);
	FDL_write(fdl,"dLambda",&msr->param.csm->dLambda);
	FDL_write(fdl,"dOmegaRad",&msr->param.csm->dOmegaRad);
	FDL_write(fdl,"dQuintess",&msr->param.csm->dQuintess);
	FDL_write(fdl,"dTheta",&msr->param.dTheta);
	FDL_write(fdl,"dTheta2",&msr->param.dTheta2);
	FDL_write(fdl,"dCentMass",&msr->param.dCentMass);
#if defined(GASOLINE) && !defined(NOCOOLING) && defined(COOLING_COSMO)
	FDL_write(fdl,"dMassFracHelium",&msr->param.CoolParam.dMassFracHelium);
#else
	{
	FLOAT dummy = 0.75; /* Nasty! */
	FDL_write(fdl,"dMassFracHelium",&dummy);
	}
#endif
	FDL_write(fdl,"dFracNoDomainDecomp",&msr->param.dFracNoDomainDecomp);
	if (msr->param.bVDetails)
		printf("Writing checkpoint file...\nTime:%g\n",dTime);
	/*
	 ** Do a parallel or serial write to the output file.
	 */
	msrCalcWriteStart(msr);
	in.iOffset = FDL_offset(fdl,"particle_array");
	if(msr->param.bParaWrite)
	    pstWriteCheck(msr->pst,&in,sizeof(in),NULL,NULL);
	else
	    msrOneNodeWriteCheck(msr, &in);
	if (msr->param.bVDetails)
		puts("Checkpoint file has been successfully written.");
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
#ifndef COLLISIONS
		if (msr->param.bVWarnings)
			printf("WARNING: Could not open redshift input file:%s\n",achFile);
#endif
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
			msr->pdOutTime[i] = csmExp2Time(msr->param.csm,a);
			break;
		case 'a':
			ret = sscanf(&achIn[1],"%lf",&a);
			if (ret != 1) goto NoMoreOuts;
			msr->pdOutTime[i] = csmExp2Time(msr->param.csm,a);
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
			msr->pdOutTime[i] = csmExp2Time(msr->param.csm,a);
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
	return(msr->param.csm->bComove);
	}


int msrKDK(MSR msr)
{
	return(msr->param.bCannonical && msr->param.bKDK);
	}


int msrDoSun(MSR msr)
{
	if (msr->param.bFandG) return(0);
	else return(1);
	}


double msrSoft(MSR msr)
{
	return(msr->param.dSoft);
	}


void msrSwitchTheta(MSR msr,double dTime)
{
	double a;

	a = csmTime2Exp(msr->param.csm,dTime);
	if (a >= msr->param.daSwitchTheta) msr->dCrit = msr->param.dTheta2; 
	}


#define SUPPRESSMASSCHECKREPORT      
double msrMassCheck(MSR msr,double dMass,char *pszWhere)
{
	struct outMassCheck out;
	
#ifdef GROWMASS
	out.dMass = 0.0;
#ifndef SUPPRESSMASSCHECKREPORT      
	if (msr->param.bVDetails) puts("doing mass check...");
#endif
	pstMassCheck(msr->pst,NULL,0,&out,NULL);
	if (dMass < 0.0) return(out.dMass);
	else if (fabs(dMass - out.dMass) > 1e-12*dMass) {
		printf("ERROR: Mass not conserved (%s): %.15e != %.15e!\n",
			   pszWhere,dMass,out.dMass);
		}	
#endif
	return(out.dMass);
	}

void msrMassMetalsEnergyCheck(MSR msr,double *dTotMass, 
        double *dTotMetals, double *dTotEnergy,char *pszWhere)
{
	struct outMassMetalsEnergyCheck out;
	
#ifdef GROWMASS
	out.dTotMass = 0.0;
	out.dTotMetals = 0.0;
	out.dTotEnergy = 0.0;
#else
#ifndef SUPPRESSMASSCHECKREPORT      
	if (msr->param.bVDetails) puts("doing mass check...");
#endif
	pstMassMetalsEnergyCheck(msr->pst,NULL,0,&out,NULL);
	if (*dTotMass < 0.0) {}
	else {
		if ( fabs(*dTotMass - out.dTotMass) > 1e-12*(*dTotMass) ) {
			printf("ERROR: Mass not conserved (%s): %.15e != %.15e!\n",
				   pszWhere,*dTotMass,out.dTotMass);
			}
            /* Commented out because metals are almost conserved.  Error
             * comes because some of the gas mass gets converted into stars
             * so conservation error is reported as a net loss in metals
             * Remaining metals are in stars
			 if ( fabs(*dTotMetals - out.dTotMetals) > 1e-12*(*dTotMetals) ) {
			 printf("ERROR: Metal mass not conserved (%s): %.15e != %.15e!\n",
			 pszWhere,*dTotMetals,out.dTotMetals);
			 }*/
#ifdef STARFORM
		if ( fabs(*dTotEnergy - out.dTotEnergy*msr->param.dDeltaStarForm) > 1e-12*(*dTotEnergy) ) {
			printf("ERROR: SN Energy not conserved (%s): %.15e != %.15e!\n",
				   pszWhere,*dTotEnergy,out.dTotEnergy*msr->param.dDeltaStarForm);
			}
#endif
		}
#endif
	*dTotMass = out.dTotMass;
	*dTotMetals = out.dTotMetals;
	*dTotEnergy = out.dTotEnergy;
	return;
	}

void
msrInitStep(MSR msr)
{
    struct inSetRung in;

    in.iRung = msr->param.iMaxRung - 1;
    pstSetRung(msr->pst, &in, sizeof(in), NULL, NULL);
    msr->iCurrMaxRung = in.iRung;
    }


void
msrSetRung(MSR msr, int iRung)
{
    struct inSetRung in;

    in.iRung = iRung;
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


int msrCurrMaxRungInclDF(MSR msr)
{
	int iRung = msr->iCurrMaxRung;
	if (msr->bDumpFrame && iRung < msr->df->iMaxRung) iRung = msr->df->iMaxRung;
    return iRung;
    }


double msrEta(MSR msr)
{
    return msr->param.dEta;
    }

double msrEtaCourant(MSR msr)
{
    return msr->param.dEtaCourant;
    }


void msrBallMax(MSR msr, int iRung, int bGreater)
{
    struct inBallMax in;

    in.iRung = iRung;
    in.bGreater = bGreater;
    in.dhFac = 1.0+msr->param.ddHonHLimit;
    pstBallMax(msr->pst, &in, sizeof(in), NULL, NULL);
    }

/*
 * bGreater = 1 => activate all particles at this rung and greater.
 */
void msrActiveRung(MSR msr, int iRung, int bGreater)
{
    struct inActiveRung in;

    in.iRung = iRung;
    in.bGreater = bGreater;
    pstActiveRung(msr->pst, &in, sizeof(in), &(msr->nActive), NULL);
    }

void msrActiveTypeOrder(MSR msr, unsigned int iTestMask )
{
    struct inActiveTypeOrder in;
    int nActive;

    in.iTestMask = iTestMask;
    pstActiveTypeOrder(msr->pst,&in,sizeof(in),&nActive,NULL);

    if (iTestMask & TYPE_ACTIVE)       msr->nActive      = nActive;
    if (iTestMask & TYPE_TREEACTIVE)   msr->nTreeActive   = nActive;
    if (iTestMask & TYPE_SMOOTHACTIVE) msr->nSmoothActive = nActive;
    }

void msrActiveOrder(MSR msr)
{
    pstActiveOrder(msr->pst,NULL,0,&(msr->nActive),NULL);
    }

void msrActiveExactType(MSR msr, unsigned int iFilterMask, unsigned int iTestMask, unsigned int iSetMask) 
{
    struct inActiveType in;
    int nActive;

    in.iFilterMask = iFilterMask;
    in.iTestMask = iTestMask;
    in.iSetMask = iSetMask;

    pstActiveExactType(msr->pst,&in,sizeof(in),&nActive,NULL);

    if (iSetMask & TYPE_ACTIVE      ) msr->nActive       = nActive;
    if (iSetMask & TYPE_TREEACTIVE  ) msr->nTreeActive   = nActive;
    if (iSetMask & TYPE_SMOOTHACTIVE) msr->nSmoothActive = nActive;
    }

void msrActiveType(MSR msr, unsigned int iTestMask, unsigned int iSetMask) 
{
    struct inActiveType in;
    int nActive;

    in.iTestMask = iTestMask;
    in.iSetMask = iSetMask;

    pstActiveType(msr->pst,&in,sizeof(in),&nActive,NULL);

    if (iSetMask & TYPE_ACTIVE      ) msr->nActive       = nActive;
    if (iSetMask & TYPE_TREEACTIVE  ) msr->nTreeActive   = nActive;
    if (iSetMask & TYPE_SMOOTHACTIVE) msr->nSmoothActive = nActive;
    }

void msrSetType(MSR msr, unsigned int iTestMask, unsigned int iSetMask) 
{
    struct inActiveType in;
    int nActive;

    in.iTestMask = iTestMask;
    in.iSetMask = iSetMask;

    pstSetType(msr->pst,&in,sizeof(in),&nActive,NULL);
    }

void msrResetType(MSR msr, unsigned int iTestMask, unsigned int iSetMask) 
{
    struct inActiveType in;
    int nActive;

    in.iTestMask = iTestMask;
    in.iSetMask = iSetMask;

    pstResetType(msr->pst,&in,sizeof(in),&nActive,NULL);

    if (msr->param.bVDetails) printf("nResetType: %d\n",nActive);
    }

int msrCountType(MSR msr, unsigned int iFilterMask, unsigned int iTestMask) 
{
    struct inActiveType in;
    int nActive;

    in.iFilterMask = iFilterMask;
    in.iTestMask = iTestMask;

    pstCountType(msr->pst,&in,sizeof(in),&nActive,NULL);

    return nActive;
    }

void msrActiveMaskRung(MSR msr, unsigned int iSetMask, int iRung, int bGreater) 
{
    struct inActiveType in;
    int nActive;

    in.iTestMask = (~0);
    in.iSetMask = iSetMask;
    in.iRung = iRung;
    in.bGreater = bGreater;

    pstActiveMaskRung(msr->pst,&in,sizeof(in),&nActive,NULL);

    if (iSetMask & TYPE_ACTIVE      ) msr->nActive       = nActive;
    if (iSetMask & TYPE_TREEACTIVE  ) msr->nTreeActive   = nActive;
    if (iSetMask & TYPE_SMOOTHACTIVE) msr->nSmoothActive = nActive;
    }

void msrActiveTypeRung(MSR msr, unsigned int iTestMask, unsigned int iSetMask, int iRung, int bGreater) 
{
    struct inActiveType in;
    int nActive;

    in.iTestMask = iTestMask;
    in.iSetMask = iSetMask;
    in.iRung = iRung;
    in.bGreater = bGreater;

    pstActiveTypeRung(msr->pst,&in,sizeof(in),&nActive,NULL);

    if (iSetMask & TYPE_ACTIVE      ) msr->nActive       = nActive;
    if (iSetMask & TYPE_TREEACTIVE  ) msr->nTreeActive   = nActive;
    if (iSetMask & TYPE_SMOOTHACTIVE) msr->nSmoothActive = nActive;
    }

int msrCurrRung(MSR msr, int iRung)
{
    struct inCurrRung in;
    struct outCurrRung out;

    in.iRung = iRung;
    pstCurrRung(msr->pst, &in, sizeof(in), &out, NULL);
    return out.iCurrent;
    }

void
msrGravStep(MSR msr, double dTime)
{
    struct inGravStep in;
    double a;

    if (msrComove(msr)) {
        a = csmTime2Exp(msr->param.csm,dTime);
        in.dEta = msrEta(msr)*pow(a,1.5);
        }
    else {
        in.dEta = msrEta(msr);
        }
    pstGravStep(msr->pst,&in,sizeof(in),NULL,NULL);
    }

void
msrAccelStep(MSR msr,double dTime)
{
    struct inAccelStep in;
    double a;

    in.dEta = msrEta(msr);
    a = csmTime2Exp(msr->param.csm,dTime);
    if (msr->param.bCannonical) {
		in.dVelFac = 1.0/(a*a);
	}
    else {
		in.dVelFac = 1.0;
	}
    in.dAccFac = 1.0/(a*a*a);
    in.bDoGravity = msrDoGravity(msr);
    in.bEpsAcc = msr->param.bEpsAccStep;
    in.bSqrtPhi = msr->param.bSqrtPhiStep;
    in.dhMinOverSoft = (msr->param.bLowerSoundSpeed ? msr->param.dhMinOverSoft : 0);
    pstAccelStep(msr->pst,&in,sizeof(in),NULL,NULL);
    }

void
msrDensityStep(MSR msr,double dTime)
{
    struct inDensityStep in;
    double expand;

    if (msr->param.bVDetails) printf("Calculating Rung Densities...\n");
    msrSmooth(msr,dTime,SMX_DENSITY,0);
    in.dEta = msrEta(msr);
    expand = csmTime2Exp(msr->param.csm,dTime);
    in.dRhoFac = 1.0/(expand*expand*expand);
    pstDensityStep(msr->pst,&in,sizeof(in),NULL,NULL);
    }

void
msrInitDt(MSR msr)
{
    struct inInitDt in;
    
    in.dDelta = msrDelta(msr);
    pstInitDt(msr->pst,&in,sizeof(in),NULL,NULL);
    }

void msrDtToRung(MSR msr, int iRung, double dDelta, int bAll)
{
    struct inDtToRung in;
    struct outDtToRung out;

    in.iRung = iRung;
    in.dDelta = dDelta;
    in.iMaxRung = msrMaxRung(msr);
    in.bAll = bAll;

    pstDtToRung(msr->pst, &in, sizeof(in), &out, NULL);

    if(out.iMaxRungIdeal > msrMaxRung(msr)) {
	fprintf(stderr, "WARNING, TIMESTEPS TOO LARGE: nMaxRung (%d) is greater than ideal rung (%d)\n", 
		msrMaxRung(msr), out.iMaxRungIdeal);
	}
    if (out.nMaxRung <= msr->param.nTruncateRung && out.iMaxRung > iRung) {
	if (msr->param.bVDetails)
	    printf("n_CurrMaxRung = %d  (iCurrMaxRung = %d):  Promoting particles to iCurrMaxrung = %d\n",
		   out.nMaxRung,out.iMaxRung,out.iMaxRung-1);

	in.iMaxRung = out.iMaxRung; /* Note this is the forbidden rung so no -1 here */
	pstDtToRung(msr->pst, &in, sizeof(in), &out, NULL);
	}

    msr->iCurrMaxRung = out.iMaxRung;
    }

/* Not fully tested: */
void msrTopStepSym(MSR msr, double dStep, double dTime, double dDelta, 
				   int iRung, double *pdActiveSum)
{
    double dMass = -1.0;
    int iSec;
    double dWMax, dIMax, dEMax;
	int nActive;

	if(msrCurrMaxRung(msr) >= iRung) { /* particles to be kicked? */
	    if(iRung < msrMaxRung(msr)-1) {
			if (msr->param.bVDetails) printf("Adjust, iRung: %d\n", iRung);
			msrActiveRung(msr, iRung, 1);
			msrDrift(msr,dTime,0.5*dDelta);
			dTime += 0.5*dDelta;
			msrInitDt(msr);
			if (msr->param.bGravStep || msr->param.bAccelStep) {
			    msrInitAccel(msr);
			    msrDomainDecomp(msr,iRung,1);
			    msrActiveRung(msr,iRung,1);
			    msrUpdateSoft(msr,dTime);
			    msrBuildTree(msr,0,dMass,0);
			    msrGravity(msr,dStep,msrDoSun(msr),&iSec,&dWMax,&dIMax,&dEMax,&nActive);
				if (msr->param.bGravStep) {
					msrGravStep(msr,dTime);
					}
				if (msr->param.bAccelStep) {
					msrAccelStep(msr,dTime);
					}
			    }
			if (msr->param.bDensityStep) {
			    msrDomainDecomp(msr,iRung,1);
			    msrActiveRung(msr,iRung,1);
			    msrBuildTree(msr,0,dMass,1);
			    msrDensityStep(msr,dTime);
			    }
			msrDtToRung(msr,iRung,dDelta,0);
			msrDrift(msr,dTime,-0.5*dDelta);
			dTime -= 0.5*dDelta;
			}
		/*
		 ** Actual Stepping.
		 */
		msrTopStepSym(msr, dStep, dTime, 0.5*dDelta,iRung+1,pdActiveSum);
		dStep += 1.0/(2 << iRung);
		msrActiveRung(msr, iRung, 0);
		msrInitAccel(msr);
#ifdef GASOLINE
		if(msrSphCurrRung(msr, iRung, 0)) {
			if (msr->param.bVDetails) printf("SPH, iRung: %d\n", iRung);
           	        msrActiveTypeRung(msr, TYPE_GAS, TYPE_ACTIVE, iRung, 0 );
                        msrDomainDecomp(msr,iRung,0);
           	        msrActiveTypeRung(msr, TYPE_GAS, TYPE_ACTIVE, iRung, 0 );
           	        msrActiveType(msr, TYPE_GAS, TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE );
			msrBuildTree(msr,1,-1.0,1);
			msrSmooth(msr,dTime,SMX_DENSITY,1);
			msrGetGasPressure(msr);
			if (msrDoGas(msr)) {
			  if (msr->param.bBulkViscosity) {
  			    msrReSmooth(msr,dTime,SMX_DIVVORT,1);
			    msrActiveRung(msr, iRung, 0);
			    msrReSmooth(msr,dTime,SMX_HKPRESSURETERMS,1);
			    }
			  else {
                            if (msr->param.bViscosityLimiter
				|| msr->param.bShockTracker
				|| msr->param.bStarForm) 
			      msrReSmooth(msr, dTime, SMX_DIVVORT, 1);
			    msrUpdateShockTracker(msr, dDelta);
  			    msrSphViscosityLimiter(msr, dTime);
			    msrActiveRung(msr, iRung, 0);
			    msrReSmooth(msr,dTime,SMX_SPHPRESSURETERMS,1);
                            }
			  }
			}
#endif
		if(msrCurrRung(msr, iRung)) {
		    if(msrDoGravity(msr)) {
   		        if (msr->param.bVDetails) printf("Gravity, iRung: %d\n", iRung);
                        msrDomainDecomp(msr,iRung,0);
			msrActiveRung(msr, iRung, 0);
			msrUpdateSoft(msr,dTime);
			msrBuildTree(msr,0,dMass,0);
			msrGravity(msr,dStep,msrDoSun(msr),&iSec,&dWMax,&dIMax,&dEMax,&nActive);
			*pdActiveSum += (double)nActive/msr->N;
			}
		    }
		if (msr->param.bVDetails) printf("Kick, iRung: %d\n", iRung);
		msrKickDKD(msr, dTime, dDelta);
		msrRungStats(msr);
		msrTopStepSym(msr,dStep,dTime+0.5*dDelta,0.5*dDelta,iRung+1,
					  pdActiveSum);
		}
	else {    
		if (msr->param.bVDetails) printf("Drift, iRung: %d\n",iRung-1);
		msrDrift(msr,dTime,dDelta);
		}
	}


void msrRungStats(MSR msr)
{
	if (msr->param.bVRungStat) {
		struct inRungStats in;
		struct outRungStats out;
		int i;

		printf("Rung distribution: (");
		for (i=0;i<msr->param.iMaxRung;++i) {
			if (i != 0) printf(",");
			in.iRung = i;
			pstRungStats(msr->pst,&in,sizeof(in),&out,NULL);
			msr->nRung[i] = out.nParticles;
			printf("%d",out.nParticles);
			}
		printf(")\n");
		}
	}

void msrTopStepNS(MSR msr, double dStep, double dTime, double dDelta, int
				  iRung, int iAdjust, double *pdActiveSum)
{
    double dMass = -1.0;
    int iSec;
    double dWMax,dIMax,dEMax;
	int nActive;

	if(msrCurrMaxRung(msr) >= iRung) { /* particles to be kicked? */
		if(iAdjust && (iRung < msrMaxRung(msr)-1)) {
			if (msr->param.bVDetails) printf("Adjust, iRung: %d\n", iRung);
			msrActiveRung(msr,iRung,1);
 			msrActiveType(msr,TYPE_ALL,TYPE_SMOOTHACTIVE|TYPE_TREEACTIVE);
			msrInitDt(msr);
			if (msr->param.bGravStep) {
				msrGravStep(msr,dTime);
				}
			if (msr->param.bAccelStep) {
			    msrAccelStep(msr,dTime);
				}
			if (msr->param.bDensityStep) {
				msrDomainDecomp(msr,iRung,1);
			    msrActiveRung(msr,iRung,1);
			    msrBuildTree(msr,0,dMass,1);
			    msrDensityStep(msr,dTime);
			    }
#ifdef GASOLINE
			if (msr->param.bSphStep) {
				msrSphStep(msr,dTime);
				}
#endif
			msrDtToRung(msr,iRung,dDelta,1);
			}
		/*
		 ** Actual Stepping.
		 */
		msrTopStepNS(msr,dStep,dTime,0.5*dDelta,iRung+1,0,pdActiveSum);
		dStep += 1.0/(2 << iRung);
		msrActiveRung(msr,iRung,0);
		msrDomainDecomp(msr,iRung,0);
		msrActiveRung(msr,iRung,0);
		msrInitAccel(msr);
#ifdef GASOLINE
		if(msrSphCurrRung(msr, iRung, 0)) {
			if (msr->param.bVDetails) printf("SPH, iRung: %d\n", iRung);
			msrActiveTypeRung(msr,TYPE_GAS,TYPE_ACTIVE,iRung,0);
			msrActiveType(msr,TYPE_GAS,TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
			msrBuildTree(msr,1,-1.0,1);
			msrSmooth(msr,dTime,SMX_DENSITY,1);
			msrGetGasPressure(msr);
			if (msrDoGas(msr)) {
 			   if (msr->param.bBulkViscosity) {
  			       msrReSmooth(msr,dTime,SMX_DIVVORT,1);
			       msrReSmooth(msr,dTime,SMX_HKPRESSURETERMS,1);
			       } 
			   else {
			       if (msr->param.bViscosityLimiter
				   || msr->param.bShockTracker
				   || msr->param.bStarForm)
				   msrReSmooth(msr,dTime,SMX_DIVVORT,1);
			       msrUpdateShockTracker(msr, dDelta);
			       msrSphViscosityLimiter(msr, dTime);
			       msrReSmooth(msr,dTime,SMX_SPHPRESSURETERMS,1);
			       }
			   }
		}
#endif
		if(msrCurrRung(msr, iRung)) {
		    if(msrDoGravity(msr)) {
				if (msr->param.bVDetails) printf("Gravity, iRung: %d\n", iRung);
				msrActiveType(msr,TYPE_ALL,TYPE_TREEACTIVE);
				msrUpdateSoft(msr,dTime);
				msrBuildTree(msr,0,dMass,0);
				msrGravity(msr,dStep,msrDoSun(msr),&iSec,&dWMax,&dIMax,&dEMax,
						   &nActive);
				*pdActiveSum += (double)nActive/msr->N;
				}
		    }
		if (msr->param.bVDetails) printf("Kick, iRung: %d\n", iRung);
		msrKickDKD(msr,dTime,dDelta);
		msrTopStepNS(msr,dStep,dTime+0.5*dDelta,0.5*dDelta,iRung+1,1,
					 pdActiveSum);
		}
	else {    
		if (msr->param.bVDetails) printf("Drift, iRung: %d\n", iRung-1);
		msrDrift(msr,dTime,dDelta);
		}
	}

void msrTopStepDKD(MSR msr, double dStep, double dTime, double dDelta, 
				double *pdMultiEff)
{
	int iRung = 0;

#ifdef COLLISIONS
	assert(0); /* DKD multi-stepping unsupported for COLLISIONS */
#endif

	*pdMultiEff = 0.0;
	if(msr->param.bNonSymp)
		msrTopStepNS(msr,dStep,dTime,dDelta,iRung,1,pdMultiEff);
	else
		msrTopStepSym(msr,dStep,dTime,dDelta,iRung,pdMultiEff);

	if (msr->param.bVStep)
		printf("Multistep Efficiency (average number of microsteps per step):%f\n",
			   *pdMultiEff);
	}

void msrTopStepKDK(MSR msr,
				   double dStep,	/* Current step */
				   double dTime,	/* Current time */
				   double dDelta,	/* Time step */
				   int iRung,		/* Rung level */
				   int iKickRung,	/* Gravity on all rungs from iRung
									   to iKickRung */
				   int iAdjust,		/* Do an adjust? */
				   double *pdActiveSum,
				   double *pdWMax,
				   double *pdIMax,
				   double *pdEMax,
				   int *piSec)
{
    double dMass = -1.0;
    int nActive;

    if(iAdjust && (iRung < msrMaxRung(msr)-1)) {
		if (msr->param.bVDetails) printf("Adjust, iRung: %d\n",iRung);
		msrActiveRung(msr, iRung, 1);
		msrActiveType(msr,TYPE_ALL,TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
		msrInitDt(msr);
		if (msr->param.bGravStep) {
			msrGravStep(msr,dTime);
			}
		if (msr->param.bAccelStep) {
		    msrAccelStep(msr,dTime);
			}
		if (msr->param.bDensityStep) {
			msrDomainDecomp(msr,iRung,1);
		    msrActiveRung(msr,iRung,1);
		    msrBuildTree(msr,0,dMass,1);
		    msrDensityStep(msr,dTime);
		    }
		if (msr->param.bDeltaAccelStep) {
		  /* Ensure we have a Density tree ready to use */
			if (!msr->param.bDeltaAccelStepGasTree || 
				msr->iTreeType != MSR_TREE_DENSITY) {
#ifdef DELTAACCELACTIVE
			    msrActiveRung(msr,iRung,1);
			    msrActiveType(msr,TYPE_ACTIVE,TYPE_TREEACTIVE);
		        msrBuildTree(msr,0,-1.0,1);
#else
			    msrActiveType(msr,TYPE_ALL,TYPE_TREEACTIVE);
			    msrBuildTree(msr,0,-1.0,1);
#endif
			    }
		    msrActiveRung(msr,iRung,1);
			msrActiveType(msr,TYPE_ACTIVE,TYPE_SMOOTHACTIVE);
			msrResetType(msr,TYPE_GAS,TYPE_SMOOTHDONE);
			/* This smooth sets dt directly -- hardwired coefficient */
			msrSmooth(msr,dTime,SMX_DELTAACCEL,0);
		    }

#ifdef GASOLINE
		if (msr->param.bSphStep) {
			msrSphStep(msr,dTime);
			}
#endif
		msrDtToRung(msr,iRung,dDelta,1);
		if (iRung == 0) {
		  /*
		  msrReorder(msr);
		  msrOutArray(msr,"test.dt",OUT_DT_ARRAY);
		  msrActiveOrder(msr);
		  */
		  msrRungStats(msr);
		  }
		}
    if (msr->param.bVDetails) printf("Kick, iRung: %d\n",iRung);
    msrActiveRung(msr,iRung,0);
    msrActiveType(msr,TYPE_ALL,TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE );
#ifdef GASOLINE
    msrUpdateuDot(msr,dTime,0.5*dDelta,1);
#endif
    msrKickKDKOpen(msr,dTime,0.5*dDelta);
    if (msrCurrMaxRungInclDF(msr) > iRung) {
		/*
		 ** Recurse.
		 */
		msrTopStepKDK(msr,dStep,dTime,0.5*dDelta,iRung+1,iRung+1,0,
					  pdActiveSum,pdWMax,pdIMax,pdEMax,piSec);
		dStep += 1.0/(2 << iRung);
		dTime += 0.5*dDelta;
		msrActiveRung(msr,iRung,0);
#ifdef GASOLINE
		msrUpdateuDot(msr,dTime,0.5*dDelta,0); /* Need forward uDot for Upreds */
#endif
		msrTopStepKDK(msr,dStep,dTime,0.5*dDelta,iRung+1,iKickRung,1,
					  pdActiveSum,pdWMax,pdIMax,pdEMax,piSec);
		}
    else {
		/* This Drifts everybody */
		if (msr->param.bVDetails) printf("Drift, iRung: %d\n", iRung);
#ifdef GASOLINE
		msrDrift(msr,dTime,0.5*dDelta);
		dTime += 0.5*dDelta;
		dStep += 1.0/(2 << iRung);
		msrActiveRung(msr,iRung,0);
		msrUpdateuDot(msr,dTime,0.5*dDelta,0); /* Need forward uDot for uPred */
		msrDrift(msr,dTime,0.5*dDelta);
		dTime += 0.5*dDelta;
		dStep += 1.0/(2 << iRung);
#else
		msrDrift(msr,dTime,dDelta);
		dTime += dDelta;
		dStep += 1.0/(1 << iRung);
#endif
#ifdef SIMPLESF
		msrSimpleStarForm(msr, dTime, dDelta);
#endif
#ifdef STARFORM
                /* only form stars at user defined intervals */
                if ( iKickRung <= msr->param.iStarFormRung )
                    msrFormStars(msr, dTime, max(dDelta,msr->param.dDeltaStarForm));
#endif
                /* only accrete onto sinks at user defined intervals */
                if ( iKickRung <= msr->param.iSinkRung )
                    msrDoSinks(msr);
		/* 
		 ** Dump Frame
		 */
		if (msr->param.dDumpFrameTime > 0 && dTime >= msr->df->dTime)
			msrDumpFrame( msr, dTime, dStep );
		else if (msr->param.dDumpFrameStep > 0 && dStep >= msr->df->dStep) 
			msrDumpFrame( msr, dTime, dStep );

		/* 
		 ** Calculate Forces (if required)
		 */
                msrActiveType(msr,TYPE_ALL,TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE );
		msrActiveMaskRung(msr,TYPE_ACTIVE,iKickRung,1);
		if (msr->nActive) {
			msrDomainDecomp(msr,iKickRung,1);
			msrInitAccel(msr);

			if (msr->param.bVStep) printf("Forces, Step:%f nActive %i\n",dStep,msr->nActive);
			if(msrDoGravity(msr)) {
				if (msr->param.bDoSelfGravity) {
					msrActiveRung(msr,iKickRung,1);
					msrUpdateSoft(msr,dTime);
					msrActiveType(msr,TYPE_ALL,TYPE_TREEACTIVE);
					if (msr->param.bVDetails)
						printf("Gravity, iRung: %d to %d\n", iRung, iKickRung);
					msrBuildTree(msr,0,dMass,0);
					}
				msrGravity(msr,dStep,msrDoSun(msr),piSec,pdWMax,pdIMax,pdEMax,&nActive);
				*pdActiveSum += (double)nActive/msr->N;
				}
			
#ifdef GASOLINE
			if (msr->param.bVDetails)
				printf("SPH, iRung: %d to %d\n",iRung,iKickRung);
			msrActiveTypeRung(msr,TYPE_GAS,TYPE_ACTIVE,iKickRung,1);
			msrActiveType(msr,TYPE_GAS,TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
			printf("nActive %d nTreeActive %d nSmoothActive %d\n",msr->nActive,
				   msr->nTreeActive,msr->nSmoothActive);
			if(msrDoGas(msr) && msrSphCurrRung(msr,iKickRung,1)) {
				msrBuildTree(msr,1,-1.0,1);
 				msrActiveTypeRung(msr,TYPE_GAS,TYPE_ACTIVE,iKickRung,1);
				if (msr->param.bFastGas && msr->nActive < msr->nGas*msr->param.dFracFastGas) {
			        msrResetType(msr,TYPE_GAS,TYPE_SMOOTHDONE|TYPE_NbrOfACTIVE|TYPE_Scatter|TYPE_DensZeroed );
					msrActiveType(msr,TYPE_ACTIVE,TYPE_SMOOTHACTIVE|TYPE_DensACTIVE );
					if (msr->param.bVDetails)
						printf("Dens Active Particles: %d\n",msr->nSmoothActive );
					/* Density for Actives and mark Gather neighbours */
					msrSmooth(msr,dTime,SMX_MARKDENSITY,1); 
					/* mark Scatter Neighbours */
					msrMarkSmooth(msr,dTime,1,TYPE_Scatter); 
					/* They need density too... */
					msrActiveType(msr,TYPE_ACTIVE|TYPE_NbrOfACTIVE|TYPE_Scatter, TYPE_DensACTIVE );
					/* ...but don't redo Actives in smooth*/
					msrActiveExactType(msr,TYPE_DensACTIVE|TYPE_ACTIVE, 
									   TYPE_DensACTIVE,TYPE_SMOOTHACTIVE);
					/* Density for Neighbours and mark scatter neighbours of actives */
					msrSmooth(msr,dTime,SMX_MARKIIDENSITY,1);
					/* mark Scatter Neighbours of Neighbours */
					msrMarkSmooth(msr,dTime,1,TYPE_Scatter);
					/* Scatter Density contribution from those particles... */
					msrActiveExactType(msr,TYPE_Scatter|TYPE_DensACTIVE, 
									   TYPE_Scatter,TYPE_SMOOTHACTIVE);
					msrSmooth(msr,dTime,SMX_MARKIIDENSITY,1); 
					
					/* We want direct neighbours of Actives only */
					msrActiveType(msr,TYPE_NbrOfACTIVE,TYPE_SMOOTHACTIVE);
					if (msr->param.bVDetails)
						printf("Density Zeroed: %d ",
							   msrCountType(msr,TYPE_DensZeroed,TYPE_DensZeroed));
					if (msr->param.bVDetails)
						printf("Neighbours: %d ",
							   msrCountType(msr,TYPE_ACTIVE|TYPE_NbrOfACTIVE,TYPE_NbrOfACTIVE));
			        }
				else {
					msrResetType(msr,TYPE_GAS,TYPE_SMOOTHDONE|TYPE_NbrOfACTIVE );
					msrActiveType(msr,TYPE_ACTIVE,TYPE_SMOOTHACTIVE|TYPE_DensACTIVE );
					if (msr->param.bVDetails)
						printf("Dens Active Particles: %d\n",msr->nSmoothActive );
					msrActiveType(msr,TYPE_GAS,TYPE_SMOOTHACTIVE );
					
					msrSmooth(msr,dTime,SMX_MARKDENSITY,1);
					
					if (msr->param.bVDetails)
						printf("Neighbours: %d ",
							   msrCountType(msr,TYPE_NbrOfACTIVE,TYPE_NbrOfACTIVE ) );
					msrActiveType(msr,TYPE_NbrOfACTIVE,TYPE_SMOOTHACTIVE);
			        }

				if (msr->param.bVDetails)
					printf("Smooth Active Particles: %d\n",msr->nSmoothActive);

				if (msr->param.bViscosityLimiter
					|| msr->param.bBulkViscosity
					|| msr->param.bShockTracker
					|| msr->param.bStarForm) {
					msrReSmooth(msr,dTime,SMX_DIVVORT,1);
					}
				msrSphViscosityLimiter(msr, dTime);
				
				msrGetGasPressure(msr);
				
				if (msr->param.bShockTracker) { 
			        msrReSmooth(msr,dTime,SMX_SPHPRESSURE,1);
					msrUpdateShockTracker(msr, dDelta);
			        if (msr->param.bBulkViscosity) 
						msrReSmooth(msr,dTime,SMX_HKVISCOSITY,1);
					else
						msrReSmooth(msr,dTime,SMX_SPHVISCOSITY,1);
			        }
				else {
			        if (msr->param.bBulkViscosity) 
						msrReSmooth(msr,dTime,SMX_HKPRESSURETERMS,1);     
					else
						msrReSmooth(msr,dTime,SMX_SPHPRESSURETERMS,1); 
			        }

				msrBallMax(msr,iKickRung,1);
				}
#endif /* GASOLINE */
			}

		}
    if (msr->param.bVDetails) printf("Kick, iRung: %d\n",iRung);
    msrActiveRung(msr,iRung,0);
#ifdef GASOLINE
    msrUpdateuDot(msr,dTime,0.5*dDelta,1);
#endif
    msrKickKDKClose(msr,dTime,0.5*dDelta);

	}

int
msrMaxOrder(MSR msr)
{
    return msr->nMaxOrder;
    }

void
msrAddDelParticles(MSR msr)
{
    struct outColNParts *pColNParts;
    int *pNewOrder;
    struct inSetNParts in;
    struct inSetParticleTypes intype;
    int iOut;
    int i;
    
    if (msr->param.bVDetails) printf("Changing Particle number\n");
    pColNParts = malloc(msr->nThreads*sizeof(*pColNParts));
    pstColNParts(msr->pst, NULL, 0, pColNParts, &iOut);
    /*
     * Assign starting numbers for new particles in each processor.
     */
    pNewOrder = malloc(msr->nThreads*sizeof(*pNewOrder));
    for(i=0;i<msr->nThreads;i++) {
		/*
		 * Detect any changes in particle number, and force a tree
		 * build.
		 */
		if (pColNParts[i].nNew != 0 || pColNParts[i].nDeltaGas != 0 ||
			pColNParts[i].nDeltaDark != 0 || pColNParts[i].nDeltaStar != 0)
			msr->iTreeType = MSR_TREE_NONE;
		pNewOrder[i] = msr->nMaxOrder + 1;
		msr->nMaxOrder += pColNParts[i].nNew;
		msr->nGas += pColNParts[i].nDeltaGas;
		msr->nDark += pColNParts[i].nDeltaDark;
		msr->nStar += pColNParts[i].nDeltaStar;
		}
    msr->N = msr->nGas + msr->nDark + msr->nStar;
#ifndef GASOLINE
    msr->nMaxOrderDark = msr->nMaxOrder;
#endif

    pstNewOrder(msr->pst,pNewOrder,sizeof(*pNewOrder)*msr->nThreads,NULL,NULL);

    if (msr->param.bVDetails)
	printf("New numbers of particles: %d gas %d dark %d star\n",
	       msr->nGas, msr->nDark, msr->nStar);
    
    in.nGas = msr->nGas;
    in.nDark = msr->nDark;
    in.nStar = msr->nStar;
    in.nMaxOrderGas = msr->nMaxOrderGas;
    in.nMaxOrderDark = msr->nMaxOrderDark;
    in.nMaxOrder = msr->nMaxOrder;
    pstSetNParts(msr->pst,&in,sizeof(in),NULL,NULL);
    intype.nSuperCool = msr->param.nSuperCool;
	/* This shouldn't really be necessary -- it is undesirable to do a fix-up like this */
    pstSetParticleTypes(msr->pst,&intype,sizeof(intype),NULL,NULL); 

    i = msrCountType(msr, TYPE_GAS, TYPE_GAS);
    assert(i == msr->nGas);
    i = msrCountType(msr, TYPE_DARK, TYPE_DARK);
    assert(i == msr->nDark);
    i = msrCountType(msr, TYPE_STAR, TYPE_STAR);
    assert(i == msr->nStar);

    free(pNewOrder);
    free(pColNParts);
    }

void
msrDoSinks(MSR msr)
{
	double sec,sec1,dsec,dMass;
	int nAccreted;

    if(msr->param.bDoSinks == 0 || msr->nStar == 0) return;
	    
	sec = msrTime();

    dMass = msrMassCheck(msr, -1, "Accrete onto Sinks: Initial Value");

    if (msr->iTreeType != MSR_TREE_DENSITY) {
	    msrActiveType(msr,TYPE_GAS,TYPE_TREEACTIVE);
		msrBuildTree(msr,0,-1.0,1);
	    }
	msrResetType(msr,TYPE_STAR,TYPE_SMOOTHDONE);
	msrActiveType(msr,TYPE_STAR,TYPE_ACTIVE|TYPE_SMOOTHACTIVE);
	/* Using smooth we kill only up to nSmooth particles at once */
	msrSmooth(msr,0.0,SMX_SINKACCRETE,1);

	nAccreted = msr->nGas;

    msrMassCheck(msr, dMass, "Accrete onto Sinks: before particle adjustment");

    msrAddDelParticles(msr);
    msrMassCheck(msr, dMass, "Accrete onto Sinks: after particle adjustment");

	nAccreted -= msr->nGas;

	sec1 = msrTime();
	dsec = sec1 - sec;
	printf("Sinks Done (%d accreted) Calculated, Wallclock: %f secs\n\n",nAccreted,dsec);

    }

int msrDoDensity(MSR msr)
{
	return(msr->param.bDoDensity);
	}

int msrDoGravity(MSR msr)
{
	return(msr->param.bDoGravity);
	}

int msrDoGas(MSR msr)
{
	return(msr->param.bDoGas);
	}

void msrInitAccel(MSR msr)
{
	pstInitAccel(msr->pst,NULL,0,NULL,NULL);
	}

void msrInitTimeSteps(MSR msr,double dTime,double dDelta) 
{
	double dMass = -1.0;

	msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
	msrInitDt(msr);
	if (msr->param.bGravStep) {
		msrGravStep(msr,dTime);
		}
	if (msr->param.bAccelStep) {
		msrAccelStep(msr,dTime);
		}
	if (msr->param.bDensityStep) {
		msrDomainDecomp(msr,0,1);
		msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
		msrBuildTree(msr,0,dMass,1);
		msrDensityStep(msr,dTime);
		}
#ifdef GASOLINE
	if (msr->param.bSphStep) {
		msrSphStep(msr,dTime);
		}
#endif
	msrDtToRung(msr,0,dDelta,1);
	msrRungStats(msr);
	}

#ifdef GASOLINE

void msrGetGasPressure(MSR msr)
{
	struct inGetGasPressure in;
  
	in.iGasModel = (enum GasModel) msr->param.iGasModel;

	switch (in.iGasModel) {

	case GASMODEL_ADIABATIC:
	case GASMODEL_ISOTHERMAL:
	case GASMODEL_COOLING:
		in.gamma = msr->param.dConstGamma;
		in.gammam1 = in.gamma-1;
		break;
	case GASMODEL_GLASS:
#ifdef GLASS
		in.dGlassPoverRhoL = msr->param.dGlassPoverRhoL;
		in.dGlassPoverRhoR = msr->param.dGlassPoverRhoR;
		in.dGlassxL = msr->param.dGlassxL;
		in.dGlassxR = msr->param.dGlassxR;
		in.dxBoundL = -0.5*msr->param.dxPeriod;
		in.dxBoundR = +0.5*msr->param.dxPeriod;
#else
		assert(0);
#endif
		break;
		}

	pstGetGasPressure(msr->pst,&in,sizeof(in),NULL,NULL);

	if (msr->param.bLowerSoundSpeed) msrLowerSoundSpeed(msr);
	}

void msrLowerSoundSpeed(MSR msr)
{
	struct inLowerSoundSpeed in;
  
	in.dhMinOverSoft = msr->param.dhMinOverSoft;

	pstLowerSoundSpeed(msr->pst,&in,sizeof(in),NULL,NULL);
	}


void msrUpdateuDot(MSR msr,double dTime,double dDelta,int bUpdateY)
{
	struct inUpdateuDot in;
	struct outUpdateuDot out;
	double a;
	
	in.duDelta = dDelta;
	dTime += dDelta/2.0;
	a = csmTime2Exp(msr->param.csm,dTime);
	in.z = 1/a - 1;
	in.dTime = dTime;
	in.iGasModel = msr->param.iGasModel;
	in.bUpdateY = bUpdateY;

	pstUpdateuDot(msr->pst,&in,sizeof(in),&out,NULL);

	printf("UpdateUdot: Avg Wallclock %f, Max Wallclock %f\n",
	       out.SumTime/out.nSum,out.MaxTime);
	}

void msrUpdateShockTracker(MSR msr,double dDelta)
{
        struct inUpdateShockTracker in;

        if (!msr->param.bShockTracker) return;

	in.dDelta = dDelta;
	in.dShockTrackerA = msr->param.dShockTrackerA;
	in.dShockTrackerB = msr->param.dShockTrackerB;

	pstUpdateShockTracker(msr->pst,&in,sizeof(struct inUpdateShockTracker),NULL,NULL);

	printf("UpdateShockTracker Done.\n");
	}

void msrInitSph(MSR msr,double dTime)
{
#ifndef NOCOOLING
	struct inInitEnergy in;
#endif
	double a;

	msrActiveType(msr,TYPE_GAS,TYPE_ACTIVE|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);

	msrBuildTree(msr,1,-1.0,1);
	msrSmooth(msr,dTime,SMX_DENSITY,1);
	msrBallMax(msr, 0, 1);

#ifndef NOCOOLING
	switch (msr->param.iGasModel) {
	case GASMODEL_COOLING:
	    if(msr->param.bRestart)
		break;		/* Already OK from checkpoint */
		/*
		* Get a consistent initial state where energy is consistent with 
		* the initial density and input temperature and the ionization
		* fraction is the consistent equilibrium state.
		**/
		in.dTuFac = msr->param.dGasConst/(msr->param.dConstGamma - 1)/
			msr->param.dMeanMolWeight;
		a = csmTime2Exp(msr->param.csm,dTime);
		in.z = 1/a - 1;
		in.dTime = dTime;
		pstInitEnergy(msr->pst, &in, sizeof(in), NULL, NULL);
		break;
		}
#endif

	if (msrDoGas(msr)) {
	        if (msr->param.bViscosityLimiter || msr->param.bBulkViscosity
		    || msr->param.bStarForm) {
		        msrReSmooth(msr,dTime,SMX_DIVVORT,1);
			}
		msrSphViscosityLimiter(msr, dTime);

		msrGetGasPressure(msr);
			
		if (msr->param.bShockTracker) { 
			msrReSmooth(msr,dTime,SMX_SPHPRESSURE,1);
			msrUpdateShockTracker(msr, 0.0);
			if (msr->param.bBulkViscosity) 
			        msrReSmooth(msr,dTime,SMX_HKVISCOSITY,1);
			else
			        msrReSmooth(msr,dTime,SMX_SPHVISCOSITY,1);
		        }
		else {
			if (msr->param.bBulkViscosity) 
				msrReSmooth(msr,dTime,SMX_HKPRESSURETERMS,1);     
			else
				msrReSmooth(msr,dTime,SMX_SPHPRESSURETERMS,1); 
		        }

   	        msrUpdateuDot(msr,dTime,0.5*msr->param.dDelta,0);
		}

	}

#ifndef NOCOOLING
void msrInitCooling(MSR msr)
{
	int cntTable;
	int nTableRows, nTableColumns;
	void *dTableData = NULL;
	char TableFileSuffix[20];
	struct inInitCooling in;

	in.dGmPerCcUnit = msr->param.dGmPerCcUnit;
	in.dComovingGmPerCcUnit = msr->param.dComovingGmPerCcUnit;
	in.dErgPerGmUnit = msr->param.dErgPerGmUnit;
	in.dSecUnit = msr->param.dSecUnit;
	in.dKpcUnit = msr->param.dKpcUnit;
	in.z = 60.0; /*dummy value*/
	in.dTime = 0.0; /* dummy value */
	in.CoolParam = msr->param.CoolParam;

	pstInitCooling(msr->pst,&in,sizeof(struct inInitCooling),NULL,NULL);

	/* Read in tables from files as necessary */
	cntTable = 0;
	for (;;) {
		CoolTableReadInfo( &msr->param.CoolParam, cntTable, &nTableColumns, TableFileSuffix );
		if (!nTableColumns) break;

		cntTable++;
		nTableRows = msrReadASCII(msr, TableFileSuffix, nTableColumns, NULL);
		if (nTableRows) {
			assert( sizeof(double)*nTableRows*nTableColumns <= CL_NMAXBYTETABLE );
			dTableData = malloc(sizeof(double)*nTableRows*nTableColumns);
			assert( dTableData != NULL );
			nTableRows = msrReadASCII(msr,TableFileSuffix, 7, dTableData );
			
			pstCoolTableRead(msr->pst,dTableData,sizeof(double)*nTableRows*nTableColumns,NULL,NULL);
			}
		}
	}
#endif

int msrSphCurrRung(MSR msr, int iRung, int bGreater)
{
    struct inSphCurrRung in;
    struct outSphCurrRung out;

    in.iRung = iRung;
    in.bGreater = bGreater;
    pstSphCurrRung(msr->pst,&in,sizeof(in),&out,NULL);
    return out.iCurrent;
    }

void msrSphStep(MSR msr, double dTime)
{
    struct inSphStep in;
    
    if (!msrDoGas(msr)) return;

    in.dCosmoFac = csmTime2Exp(msr->param.csm,dTime);
    in.dEtaCourant = msrEtaCourant(msr);
    in.dEtauDot = msr->param.dEtauDot;
    in.bViscosityLimitdt = msr->param.bViscosityLimitdt;
    pstSphStep(msr->pst,&in,sizeof(in),NULL,NULL);
    }

void msrSphViscosityLimiter(MSR msr, double dTime)
{
    struct inSphViscosityLimiter in;
    
    in.bOn = msr->param.bViscosityLimiter;
    in.bShockTracker = msr->param.bShockTracker;

    pstSphViscosityLimiter(msr->pst,&in,sizeof(in),NULL,NULL);
    }

#endif /* GASOLINE */

int msrDumpFrameInit(MSR msr, double dTime, double dStep, int bRestart) {
	LCL *plcl = &msr->lcl;
	char achFile[160];
	
	if (msr->param.dDumpFrameStep > 0 || msr->param.dDumpFrameTime > 0) {
		msr->bDumpFrame = 1;
		/*
		 ** Add Data Subpath for local and non-local names.
		 */
		achFile[0] = 0;
		sprintf(achFile,"%s/%s.director",msr->param.achDataSubPath,
				msr->param.achOutName);
		
		dfInitialize( &msr->df, dTime, msr->param.dDumpFrameTime, dStep, 
					 msr->param.dDumpFrameStep, msr->param.dDelta, 
					 msr->param.iMaxRung, msr->param.bVDetails,
					 achFile );

		/* Read in photogenic particle list */
		if (msr->df->bGetPhotogenic) {
		  achFile[0] = 0;
		  sprintf(achFile,"%s/%s.photogenic",msr->param.achDataSubPath,
				msr->param.achOutName);
		  msrSetTypeFromFile( msr, achFile, TYPE_PHOTOGENIC );
		}

		if(!bRestart)
			msrDumpFrame( msr, dTime, dStep );
                return 1;
		} else { return 0; }
	}

void msrDumpFrame(MSR msr, double dTime, double dStep)
{
	double sec,dsec1,dsec2,dExp;

	sec = msrTime();

	if (msr->df->iDimension == DF_3D) {
#ifdef VOXEL
		/* 3D Voxel Projection */
		struct inDumpVoxel in;
		assert(0);

		dfSetupVoxel( msr->df, dTime, dStep, &in );

		pstDumpVoxel(msr->pst, &in, sizeof(struct inDumpVoxel), NULL, NULL );
		dsec1 = msrTime() - sec;
		
		dfFinishVoxel( msr->df, dTime, dStep, &in );
		
		dsec2 = msrTime() - sec;
		
		printf("DF Dumped Voxel %i at %g (Wallclock: Render %f tot %f secs).\n",
			   msr->df->nFrame-1,dTime,dsec1,dsec2);
#endif
		}
	else {
		/* 2D Projection */
	    struct inDumpFrame in;
		void *Image; 
		int nImage;
		double com[12];

		if (msr->df->bGetCentreOfMass) {
		  pstCOM(msr->pst, NULL, 0, &com[0], NULL);
		  }

		if (msr->df->bGetPhotogenic) {
		  int type = TYPE_PHOTOGENIC;
		  pstCOMByType(msr->pst, &type, sizeof(int), &com[0], NULL);
		  }

		if (msr->df->bGetOldestStar) {
		  pstOldestStar(msr->pst, NULL, 0, &com[0], NULL);
		  }

		dExp = csmTime2Exp(msr->param.csm,dTime);
		dfSetupFrame( msr->df, dTime, dStep, dExp, &com[0], &in );

		Image = dfAllocateImage( in.nxPix, in.nyPix );
		
		pstDumpFrame(msr->pst, &in, sizeof(struct inDumpFrame), Image, &nImage );
		dsec1 = msrTime() - sec;
		
		dfFinishFrame( msr->df, dTime, dStep, &in, Image );
		
		dsec2 = msrTime() - sec;
		
		printf("DF Dumped Image %i at %g (Wallclock: Render %f tot %f secs).\n",
			   msr->df->nFrame-1,dTime,dsec1,dsec2);
		
		dfFreeImage( Image );
		}
	}

void msrFormStars(MSR msr, double dTime, double dDelta)
{
#ifdef STARFORM
    struct inFormStars in;
    struct outFormStars outFS;
    struct inFeedback inFB;
    struct outFeedback outFB;
    double dTotMass = -1.0, dTotMetals = -1.0, dTotEnergy = -1.0;
    double dTotSNEnergy = 0.0;
    int i;
    int iDum;

	double sec,dsec;

	sec = msrTime();

    if(msr->param.bStarForm == 0)
		return;
    
    in.dTime = dTime;
    msr->param.stfm->dDeltaT = dDelta;
    in.stfm = *msr->param.stfm;
    inFB.dTime = dTime;
    inFB.dDelta = dDelta;
    inFB.fb  = *msr->param.fb;
    inFB.sn  = *msr->param.sn;
    
    if (msr->param.bVDetails) printf("Form Stars ... ");

    msrMassMetalsEnergyCheck(msr, &dTotMass, &dTotMetals, &dTotEnergy, "Form Stars");
    
    msrActiveType(msr,TYPE_GAS,TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE|TYPE_ACTIVE);
    msrDomainDecomp(msr, 0, 1);
    msrBuildTree(msr,1,dTotMass,1);
    pstFormStars(msr->pst, &in, sizeof(in), &outFS, NULL);
    if (msr->param.bVDetails)
		printf("%d Stars formed with mass %g, %d gas deleted\n",
			   outFS.nFormed, outFS.dMassFormed, outFS.nDeleted);
    /* there are two gas particle deletion criteria:
	   
       1) in pstFormStars: gas particles with mass less than
       stfm->dMinGasMass are marked for deletion
	   
       2) in DeleteGas (see smoothfcn.c): gas particles with 
       mass less than dMinMassFrac of the average mass of neighbouring
       gas particles are also marked for deletion 
	   
       - eh, Feb 7/01*/

    /*
     * Find low mass gas particles and mark them for deletion.
	 * For better numerical treatment
     */
    if (msr->param.bVDetails) printf("Delete Gas ...\n");
    msrSmooth(msr, dTime, SMX_DELETE_GAS, 0);
    
    /*
     * Record star formation events XXX - not done.
     * NB.  At the moment each star is a star formation event.
     */
      
    /*
     * Distribute mass, and metals of deleted gas particles.
     */
    if (msr->param.bVDetails) printf("Distribute Deleted Gas ...\n");
    msrActiveType(msr, TYPE_DELETED, TYPE_SMOOTHACTIVE);
    msrSmooth(msr, dTime, SMX_DIST_DELETED_GAS, 1);
    /*
     * adjust particle numbers
     */
    msrAddDelParticles(msr);
    msrMassCheck(msr, dTotMass, "Form stars: after particle adjustment");

	dsec = msrTime() - sec;
	printf("Star Formation Calculated, Wallclock: %f secs\n\n",dsec);

    /*
     * Calculate energy of SN for any stars for the next timestep.  This
     * requires looking at past star forming events.  Also calculate
     * mass loss.
     */
    if(msr->param.bFeedBack) {
		if (msr->param.bVDetails) printf("Calculate Feedback ...\n");
		sec = msrTime();
		pstFeedback(msr->pst, &inFB, sizeof(inFB),
					&outFB, &iDum);
		if(msr->param.bVDetails) {
			printf("Feedback totals: mass, energy, metalicity\n");
			for(i = 0; i < FB_NFEEDBACKS; i++){
				printf("feedback %d: %g %g %g\n", i,
					   outFB.fbTotals[i].dMassLoss,
					   outFB.fbTotals[i].dEnergy,
					   outFB.fbTotals[i].dMassLoss != 0.0 ?
					   outFB.fbTotals[i].dMetals
					   /outFB.fbTotals[i].dMassLoss : 0.0);
                                dTotMetals += outFB.fbTotals[i].dMetals;
                                dTotSNEnergy += outFB.fbTotals[i].dEnergy;
                                }
			}


		/*
		 * spread mass lost from SN, (along with energy and metals)
		 * to neighboring gas particles.
		 */
		if (msr->param.bVDetails) printf("Distribute SN Energy ...\n");
		msrActiveType(msr, TYPE_GAS, TYPE_ACTIVE|TYPE_TREEACTIVE);
		msrBuildTree(msr,1,-1.0,1);

		msrResetType(msr, TYPE_STAR, TYPE_SMOOTHDONE);
		msrActiveType(msr, TYPE_STAR, TYPE_SMOOTHACTIVE);
		assert(msr->nSmoothActive == msr->nStar);
		msrSmooth(msr, dTime, SMX_DIST_SN_ENERGY, 1);
		msrMassMetalsEnergyCheck(msr, &dTotMass, &dTotMetals, 
                            &dTotSNEnergy, "Form stars: after feedback");

		dsec = msrTime() - sec;
		printf("Feedback Calculated, Wallclock: %f secs\n\n",dsec);
		}

#endif
    }

void msrSimpleStarForm(MSR msr, double dTime, double dDelta)
{
/* Note: Must be called with an SPH tree built and available */
#ifdef SIMPLESF
    struct inSimpleStarForm in;
    struct outSimpleStarForm out;
	double a,d1,d2;

    double dMass = -1.0;
	double sec,sec1,dsec;

    if(msr->param.bStarForm == 0) return;
    
	sec = msrTime();

    a = csmTime2Exp(msr->param.csm,dTime);

	/* Convert input parameters to code units */
    in.dRateCoeff = msr->param.SSF_dEfficiency*sqrt(4.*M_PI/pow(a,3)); /* G=1 */
	in.dTMax = msr->param.SSF_dTMax;
	d1 = msr->param.SSF_dComovingDenMin;
	d2 = msr->param.SSF_dPhysDenMin/msr->param.dGmPerCcUnit*pow(a,3);
	in.dDenMin = (d1>d2 ? d1 : d2);
	in.dDelta = dDelta;

	in.dTime = dTime;
	in.dInitStarMass = msr->param.SSF_dInitStarMass;
	in.dESNPerStarMass = msr->param.SSF_dESNPerStarMass/msr->param.dErgPerGmUnit;
#define SECONDSPERYEAR   31557600.
	in.dtCoolingShutoff = msr->param.SSF_dtCoolingShutoff*SECONDSPERYEAR/msr->param.dSecUnit;
	in.bdivv = msr->param.SSF_bdivv;
    
    if (msr->param.bVDetails) printf("Simple Star Form ... ");

    dMass = msrMassCheck(msr, -1.0, "Form Stars");
	msrActiveType(msr,0,TYPE_SMOOTHACTIVE|TYPE_SMOOTHDONE|TYPE_NbrOfACTIVE);
	/* New stars will be set to TYPE_SMOOTHACTIVE when created */
	
    pstSimpleStarForm(msr->pst, &in, sizeof(in), &out, NULL);
    if (msr->param.bVDetails)
		printf("%d Stars formed with mass %g, %d gas deleted\n",
			   out.nFormed, out.dMassFormed, out.nDeleted);

    /*
     * adjust particle numbers
     */
    msrAddDelParticles(msr);
    msrMassCheck(msr, dMass, "Form stars: after particle adjustment");

	sec1 = msrTime();
	dsec = sec1 - sec;
	printf("Star Formation Calculated, Wallclock: %f secs\n\n",dsec);

	if (msr->param.bFeedBack && out.nFormed) {
		/* Build a tree to distribute energy from SN, if nFormed > 0  */

		/* Any new stars have been set to SMOOTHACTIVE */
		msrActiveType(msr,TYPE_GAS,TYPE_TREEACTIVE|TYPE_ACTIVE);
		msrBuildTree(msr,1,dMass,1);
		msrSmooth(msr, dTime, SMX_SIMPLESF_FEEDBACK, 1);
		dsec = msrTime() - sec1;
		printf("Feedback Calculated, Wallclock: %f secs\n\n",dsec);
		}
#endif
    }

#ifdef GLASS

void msrInitGlass(MSR msr)
{
    struct inRandomVelocities in;
    
    in.dMaxVelocityL = msr->param.dGlassVL; 
    in.dMaxVelocityR = msr->param.dGlassVR; 
    pstRandomVelocities(msr->pst,&in,sizeof(in),NULL,NULL);
    }

#endif

#ifdef COLLISIONS

int
msrNumRejects(MSR msr)
{
	struct outNumRejects out;

	pstNumRejects(msr->pst,NULL,0,&out,NULL);
	return out.nRej;
	}

void
msrFindRejects(MSR msr)
{
	/*
	 ** Checks initial conditions for particles with overlapping physical
	 ** or Hill spheres. The latter case only makes sense for particles
	 ** orbiting a massive central body, like the Sun, and is controlled by
	 ** the value of msr->dCentMass. Rejects are written to REJECTS_FILE
	 ** (cf. ssdefs.h). This procedure is intended to be called iteratively
	 ** from an external initial-conditions program.
	 */

	int nRej = 0;

	if (msr->param.bVStart)	puts("Checking for rejected ICs...");
	msrDomainDecomp(msr,0,1);
	msrActiveType(msr,TYPE_ALL,TYPE_SMOOTHACTIVE|TYPE_TREEACTIVE);
	msrBuildTree(msr,0,-1.0,1);
	msrSmooth(msr,0.0,SMX_REJECTS,1); /* 1=use combiner cache */
	nRej = msrNumRejects(msr);
	if (nRej) {
		printf("%i reject%s found!\n",nRej,(nRej==1?"":"s"));
		msrReorder(msr);
		msrOutArray(msr,REJECTS_FILE,OUT_REJECTS_ARRAY);
		_msrExit(msr,1);
		}
	else {
		puts("No rejects found.");
		_msrExit(msr,0);
		}
	assert(0); /* unreachable statement */
	}

void
msrOneNodeReadSS(MSR msr,struct inReadSS *in)
{
    int i,id;
    int *nParts;
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

    pstOneNodeReadInit(msr->pst,in,sizeof(*in),nParts,&nid);
    assert(nid == msr->nThreads*sizeof(*nParts));
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
	_msrMakePath(plcl->pszDataPath,in->achInFile,achInFile);

    nStart = nParts[0];
	assert(msr->pMap[0] == 0);
    for (i=1;i<msr->nThreads;++i) {
		id = msr->pMap[i];
		/* 
		 * Read particles into the local storage.
		 */
		assert(plcl->pkd->nStore >= nParts[id]);
		pkdReadSS(plcl->pkd,achInFile,nStart,nParts[id]);
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
    pkdReadSS(plcl->pkd,achInFile,0,nParts[0]);
    }

double
msrReadSS(MSR msr)
{
	SSIO ssio;
	SSHEAD head;
	struct inReadSS in;
	struct inSetParticleTypes intype;
	char achInFile[PST_FILENAME_SIZE];
	LCL *plcl = msr->pst->plcl;
	double dTime;

	if (msr->param.achInFile[0]) {
		/*
		 ** Add Data Subpath for local and non-local names.
		 */
		_msrMakePath(msr->param.achDataSubPath,msr->param.achInFile,in.achInFile);
		/*
		 ** Add local Data Path.
		 */
		_msrMakePath(plcl->pszDataPath,in.achInFile,achInFile);

		if (ssioOpen(achInFile,&ssio,SSIO_READ)) {
			printf("Could not open InFile:%s\n",achInFile);
			_msrExit(msr,1);
			}
		}
	else {
		printf("No input file specified\n");
		_msrExit(msr,1);
		}

	/* Read header */

	if (ssioHead(&ssio,&head)) {
		printf("Could not read header of InFile:%s\n",achInFile);
		_msrExit(msr,1);
		}
	if (ssioClose(&ssio)) {
		printf("Could not close InFile:%s\n",achInFile);
		_msrExit(msr,1);
		}

	msr->N = msr->nDark = head.n_data;
	msr->nGas = msr->nStar = 0;
	msr->nMaxOrder = msr->N - 1;
	msr->nMaxOrderGas = msr->nGas - 1; /* always -1 */
	msr->nMaxOrderDark = msr->nDark - 1;

	dTime = head.time;
	if (msr->param.bVStart) {
		double tTo;
		printf("Input file...N=%i,Time=%g\n",msr->N,dTime);
		tTo = dTime + msr->param.nSteps*msr->param.dDelta;
		printf("Simulation to Time:%g\n",tTo);
		}

	in.nFileStart = 0;
	in.nFileEnd = msr->N - 1;
	in.nDark = msr->nDark;
	in.nGas = msr->nGas;	/* always zero */
	in.nStar = msr->nStar;	/* always zero */
	in.iOrder = msr->param.iOrder;
	/*
	 ** Since pstReadSS causes the allocation of the local particle
	 ** store, we need to tell it the percentage of extra storage it
	 ** should allocate for load balancing differences in the number of
	 ** particles.
	 */
	in.fExtraStore = msr->param.dExtraStore;

	in.fPeriod[0] = msr->param.dxPeriod;
	in.fPeriod[1] = msr->param.dyPeriod;
	in.fPeriod[2] = msr->param.dzPeriod;

	if (msr->param.bParaRead)
	    pstReadSS(msr->pst,&in,sizeof(in),NULL,NULL);
	else
	    msrOneNodeReadSS(msr,&in);
	if (msr->param.bVDetails) puts("Input file successfully read.");
	/*
	 ** Set particle ACTIVE flags to correspond to appropriate type.
	 */
	intype.nSuperCool = msr->param.nSuperCool;
	assert(intype.nSuperCool == 0); /* better be zero... */
	pstSetParticleTypes(msr->pst,&intype,sizeof(intype),NULL,NULL);
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

void
msrOneNodeWriteSS(MSR msr,struct inWriteSS *in)
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
	_msrMakePath(plcl->pszDataPath,in->achOutFile,achOutFile);

    /* 
     * First write our own particles.
     */
    pkdWriteSS(plcl->pkd,achOutFile,plcl->nWriteStart);
    nStart = plcl->pkd->nLocal;
	assert(msr->pMap[0] == 0);
    for (i=1;i<msr->nThreads;++i) {
		id = msr->pMap[i];
		/* 
		 * Swap particles with the remote processor.
		 */
		inswap = 0;
		mdlReqService(pst0->mdl,id,PST_SWAPALL,&inswap,sizeof(inswap));
		pkdSwapAll(plcl->pkd,id);
		mdlGetReply(pst0->mdl,id,NULL,NULL);
		/* 
		 * Write the swapped particles.
		 */
		pkdWriteSS(plcl->pkd,achOutFile,nStart);
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

void
msrWriteSS(MSR msr,char *pszFileName,double dTime)
{
	SSIO ssio;
	SSHEAD head;
	struct inWriteSS in;
	char achOutFile[PST_FILENAME_SIZE];
	LCL *plcl = msr->pst->plcl;

	/*
	 ** Calculate where each processor should start writing.
	 ** This sets plcl->nWriteStart.
	 */
	msrCalcWriteStart(msr);
	/*
	 ** Add Data Subpath for local and non-local names.
	 */
	_msrMakePath(msr->param.achDataSubPath,pszFileName,in.achOutFile);
	/*
	 ** Add local Data Path.
	 */
	_msrMakePath(plcl->pszDataPath,in.achOutFile,achOutFile);
	
	if (ssioOpen(achOutFile,&ssio,SSIO_WRITE)) {
		printf("Could not open OutFile:%s\n",achOutFile);
		_msrExit(msr,1);
		}

	/* Write header */

	head.time = dTime;
	head.n_data = msr->N;
	head.pad = -1;

	if (ssioHead(&ssio,&head)) {
		printf("Could not write header of OutFile:%s\n",achOutFile);
		_msrExit(msr,1);
		}
	if (ssioClose(&ssio)) {
		printf("Could not close OutFile:%s\n",achOutFile);
		_msrExit(msr,1);
		}

	if(msr->param.bParaWrite)
	    pstWriteSS(msr->pst,&in,sizeof(in),NULL,NULL);
	else
		msrOneNodeWriteSS(msr,&in);

	if (msr->param.bVDetails) puts("Output file successfully written.");
	}

void
msrPlanetsKDK(MSR msr,double dStep,double dTime,double dDelta,double *pdWMax,
			  double *pdIMax,double *pdEMax,int *piSec)
{
	struct inKick in;
	struct outKick out;

	msrActiveRung(msr,0,1); /* just in case */
	in.dvFacOne = 1.0;
	in.dvFacTwo = 0.5*dDelta;
    if (msr->param.bVDetails) printf("Planets Kick\n");
	pstKick(msr->pst,&in,sizeof(in),&out,NULL);
	printf("Kick: Avg Wallclock %f, Max Wallclock %f\n",
	       out.SumTime/out.nSum,out.MaxTime);
	if (msr->param.bVDetails) printf("Planets Drift\n");
	msrPlanetsDrift(msr,dStep,dTime,dDelta);
	dTime += 0.5*dDelta; /* not used */
	dStep += 1.0;
	msrActiveRung(msr,0,1);
	msrInitAccel(msr);
	if(msrDoGravity(msr)) {
		int nDum;
		if (msr->param.bVDetails) printf("Planets Gravity\n");
		msrDomainDecomp(msr,0,1);
		msrActiveRung(msr,0,1);
		msrUpdateSoft(msr,dTime);
		msrBuildTree(msr,0,-1.0,0);
		msrGravity(msr,dStep,msrDoSun(msr),piSec,pdWMax,pdIMax,pdEMax,&nDum);
		}
    if (msr->param.bVDetails) printf("Planets Kick\n");
	pstKick(msr->pst,&in,sizeof(in),&out,NULL);
	printf("Kick: Avg Wallclock %f, Max Wallclock %f\n",
	       out.SumTime/out.nSum,out.MaxTime);
	}

void
msrPlanetsDrift(MSR msr,double dStep,double dTime,double dDelta)
{
	struct inDrift in;
	double dSmall,dSubTime,dNext;

	msrMarkEncounters(msr,0); /* initialize */

	dSmall = msr->param.dSmallStep;
	dSubTime = 0;

	do {
		if (dSmall > 0)
			msrNextEncounter(msr,dSubTime,dDelta,&dNext);
		else
			dNext = dDelta;
		if (dNext > 0) { /* Kepler drift to next encounter or end of step */
			in.dDelta = dNext;
			in.bPeriodic = 0; /* (in.fCenter[] unused) */
			in.bFandG = 1;/*DEBUG redundant: iDriftType already KEPLER*/
			in.fCentMass = msr->param.dCentMass;
			pstDrift(msr->pst,&in,sizeof(in),NULL,NULL);
			msr->iTreeType = MSR_TREE_NONE;
			dTime += dNext;
			dSubTime += dNext;
			}
		if (dSubTime < dDelta) { /* Handle "normal-step" encounters */
			msrMarkEncounters(msr,dSubTime + dSmall);
			msr->param.bFandG = 0; /* force direct Sun term included */
			msrLinearKDK(msr,dStep + dSubTime/dDelta,dTime,dSmall);
			msr->param.bFandG = 1; /* restore flag */
			dTime += dSmall;
			dSubTime += dSmall;
			}
		} while (dSubTime < dDelta);
	}

void
msrNextEncounter(MSR msr,double dStart,double dEnd,double *dNext)
{
	struct outNextEncounter out;
	double sec = msrTime();

	if (msr->param.bVStep)
		printf("Start encounter search (dStart=%g,dEnd=%g)...\n",dStart,dEnd);

	/* construct q-Q tree -- must do this every time */

	/*
	 ** walk q-Q tree to get neighbour lists. each list is pruned via the
	 ** "lambda-a" phase space criterion, then further refined by solving
	 ** for intersecting ellipsoidal toroids. finally, a minimum encounter
	 ** time is assigned to each particle and the global minimum is returned.
	 ** dStart is checked to see whether this is an update or the first call
	 ** this step.
	 */

	/* new smooth operation goes here */

	/*** IF DSTART > 0, WE WANT TO DO AN UPDATE! ***/
	/* ie. RECONSIDER ALL PAIRS WITH 0 < dtCol < dStart */

	/*DEBUG FOLLOWING WON'T WORK!!! INIT'N MUST OCCUR AT PST LEVEL!!*/
	out.dt = dEnd - dStart; /* Kepler drift if no more encounters */
	pstNextEncounter(msr->pst,NULL,0,&out,NULL);

	/* dNext could be zero, indicating still working on NORMAL steps */

	*dNext = out.dt;

	if (msr->param.bVStep)
		printf("Encounter search completed, time = %g sec\n",msrTime() - sec);
	}

void
msrMarkEncounters(MSR msr,double dtMax)
{
	struct inMarkEncounters in;

	/* initializes if dTMax <= 0 */
	/* SET iDriftType to NORMAL and set ACTIVE; otherwise reset ACTIVE */
	/* note msrMarkEncounters() is necessary because msrFindEncounter()
	   must be called *before* the Kepler drift (to see how long to drift
	   for) and all particles must be updated at that time */
	in.dt = dtMax;
	pstMarkEncounters(msr->pst,&in,sizeof(in),NULL,NULL);
	}

void
msrLinearKDK(MSR msr,double dStep,double dTime,double dDelta)
{
    if (msr->param.bVDetails) printf("Linear Kick\n");
    msrKickKDKOpen(msr,dTime,0.5*dDelta);
	if (msr->param.bVDetails) printf("Linear Drift\n");
	msrDrift(msr,dTime,dDelta);
	dTime += 0.5*dDelta;
	dStep += 1.0;
	msrInitAccel(msr);
	if(msrDoGravity(msr)) {
		double dDum;
		int iDum;
		if (msr->param.bVDetails) printf("Linear Gravity\n");
		msrDomainDecomp(msr,0,1);
		msrUpdateSoft(msr,dTime);
		msrBuildTree(msr,0,-1.0,0);
		msrGravity(msr,dStep,msrDoSun(msr),&iDum,&dDum,&dDum,&dDum,&iDum);
		}
    if (msr->param.bVDetails) printf("Linear Kick\n");
    msrKickKDKClose(msr,dTime,0.5*dDelta);
	}

static char *
_msrParticleLabel(MSR msr,int iColor)
{
	/* For use with msrDoCollisions() only */

#ifdef SAND_PILE
	if (iColor < 0) {
		static char ach[256];
		WALLS *w = &msr->param.CP.walls;
		int wall;
		if (iColor < -w->nWalls) {
			int endpt;
			wall = (-iColor - w->nWalls - 1)/2;
			endpt = (-iColor - w->nWalls - 1)%2;
			(void) sprintf(ach,"WALL %i ENDPT %i",wall,endpt);
			ach[255] = '\0'; /* paranoid */
			return ach;
			}
		wall = -iColor - 1;
		(void) sprintf(ach,"WALL %i",wall);
		ach[255] = '\0';
		return ach;
		}
#endif

	switch (iColor) {
	case SUN:
		return "SUN";
	case JUPITER:
		return "JUPITER";
	case SATURN:
		return "SATURN";
	case URANUS:
		return "URANUS";
	case NEPTUNE:
		return "NEPTUNE";
	case PLANETESIMAL:
		return "PLANETESIMAL";
	default:
		return "UNKNOWN";
		}
	}

void
msrDoCollisions(MSR msr,double dTime,double dDelta)
{
	/*
	 ** Performs smooth operation to determine if a collision occurs
	 ** in the next interval. If so, the collision is processed and
	 ** collision flags are updated so that subsequent searches in the
	 ** interval can use resmooth over far fewer particles. This
	 ** continues until no further collisions occur in the interval.
	 */

	struct inSmooth smooth;
	struct outNextCollision next;
	struct inGetColliderInfo inGet;
	struct outGetColliderInfo outGet;
	struct inDoCollision inDo;
	struct outDoCollision outDo;
	struct inResetColliders reset;
	COLLIDER *c1 = &inDo.Collider1,*c2 = &inDo.Collider2,*c;
	double sec;
	unsigned int nCol=0,nMis=0,nMrg=0,nBnc=0,nFrg=0;
	int first_pass;

#ifdef SAND_PILE
	if (msr->param.nSmooth < 1) return; /* might hit a wall */
#else
	if (msr->param.nSmooth < 2) return; /* won't find any colliders */
#endif
	if (msr->param.bVStep)
		printf("Start collision search (dTime=%e,dDelta=%e)...\n",
			   dTime,dDelta);
	sec = msrTime();
	msrActiveType(msr,TYPE_ALL,TYPE_ALLACTIVE);
	smooth.nSmooth = msr->param.nSmooth;
	smooth.bPeriodic = msr->param.bPeriodic;
	smooth.bSymmetric = 0;
	smooth.iSmoothType = SMX_COLLISION;
	smooth.dfBall2OverSoft2 = 0.0; /* No softening limit */
	smooth.smf.dTime = dTime;
	smooth.smf.dStart = 0.0;
	smooth.smf.dEnd = dDelta;
	smooth.smf.dCollapseLimit = msr->param.CP.dCollapseLimit;
#ifdef SLIDING_PATCH
	smooth.smf.fLx = msr->param.dxPeriod;
	smooth.smf.dOrbFreq = msr->param.dOrbFreq;
#endif
#ifdef SAND_PILE
	smooth.smf.walls = msr->param.CP.walls; /* structure copy */
#endif
	inDo.bPeriodic = smooth.bPeriodic;
#ifdef SLIDING_PATCH
	inDo.dOrbFreq = smooth.smf.dOrbFreq;
	inDo.dTime = smooth.smf.dTime;
#endif
	first_pass = 1;
	do {
		if (first_pass) {
			if (msr->iTreeType != MSR_TREE_DENSITY) msrBuildTree(msr,0,-1.0,1);
			pstSmooth(msr->pst,&smooth,sizeof(smooth),NULL,NULL);
			first_pass = 0;
			}
		else {
			assert(msr->iTreeType == MSR_TREE_DENSITY);
			/* following assumes inSmooth and inReSmooth are identical... */
			/*DEBUG pstReSmooth() has caused problems in parallel (MDL_CACHE_LINE):
			  see DCR's e-mails in pkd folder around end of Oct '01.*/
			pstReSmooth(msr->pst,&smooth,sizeof(smooth),NULL,NULL);
			}
		pstNextCollision(msr->pst,NULL,0,&next,NULL);
		/*
		 ** The following assert ensures that no two collisions occur at
		 ** precisely the same instant. Physically this is possible but
		 ** in our case it's probably a sign of trouble, so it's not
		 ** allowed for now. Note CheckForCollision() asserts this too,
		 ** for a more limited case. CAVEAT: the tests rely on *exact*
		 ** double-precision equality which is almost impossible to achieve
		 ** -- the test should really be fabs(dt - dt0) < PRECISION, but
		 ** this is expensive. Consequently, a simultaneous collision may
		 ** be missed, resulting in an overlap during the next step!
		 */
#ifndef FIX_COLLAPSE
		assert(next.dt > smooth.smf.dStart);
#endif
		/* process the collision */
		if (COLLISION(next.dt)) {
			assert(next.iOrder1 >= 0);
#ifndef SAND_PILE
			assert(next.iOrder2 >= 0);
#endif
			inDo.dt = next.dt;
			inGet.iOrder = next.iOrder1;
			pstGetColliderInfo(msr->pst,&inGet,sizeof(inGet),&outGet,NULL);
			*c1 = outGet.Collider; /* struct copy */
			assert(c1->id.iOrder == inGet.iOrder);
#ifdef SAND_PILE
			if (next.iOrder2 < 0) { /* wall or endpoint collision */
				c2->id.iOrder = next.iOrder2;
				if (msr->param.iCollLogOption == COLL_LOG_VERBOSE) {
					/* the rest is for logging purposes only */
					c2->id.iPid = c2->id.iIndex = -1;
					c2->fMass = DBL_MAX;
					c2->fRadius = 0.0;
					c2->r[0] = c2->r[1] = c2->r[2] =
						c2->v[0] = c2->v[1] = c2->v[2] =
							c2->w[0] = c2->w[1] = c2->w[2] = 0.0;
					c2->dt = next.dt;
					c2->iColor = next.iOrder2;
					}
				}
			else {
#endif
				inGet.iOrder = next.iOrder2;
				pstGetColliderInfo(msr->pst,&inGet,sizeof(inGet),&outGet,NULL);
				*c2 = outGet.Collider;
				assert(c2->id.iOrder == inGet.iOrder);
#ifdef SAND_PILE
				}
#endif
			inDo.CP = msr->param.CP;
#ifdef AGGS
			if (COLLIDER_IS_AGG(c1))
				msrAggsAdvanceCollider(msr,c1,next.dt);
			if (COLLIDER_IS_AGG(c2))
				msrAggsAdvanceCollider(msr,c2,next.dt);
			inDo.CP.iAggNewID = msr->iAggNewID; /*DEBUG temporary cheat*/
#endif
			pstDoCollision(msr->pst,&inDo,sizeof(inDo),&outDo,NULL);
			msr->dTcoll += outDo.dT; /* account for kinetic energy loss */
			++nCol;
			switch (outDo.iOutcome) {
#ifdef FIX_COLLAPSE
			case MISS:
				++nMis;
				--nCol;
				break;
#endif
			case MERGE:
				++nMrg;
				break;
			case BOUNCE:
				++nBnc;
				break;
			case FRAG:
				++nFrg;
				break;
			default:
				assert(0); /* unknown outcome */
				}
			/*
			 ** If there's a merger, the deleted particle has ACTIVE set
			 ** to zero, so we can still use the old tree and ReSmooth().
			 ** Deleted particles are cleaned up outside this loop.
			 */
#ifdef IGNORE_FOR_NOW/*DEBUG*/
			if (outDo.iOutcome & FRAG) {
				/* note this sets msr->iTreeType to MSR_TREE_NONE */
				msrAddDelParticles(msr);
				/* have to rebuild tree here and do new smooth... */
				/* also would have to reactivate all particles */
				}
#endif
			reset.iOrder1 = c1->id.iOrder;
			reset.iOrder2 = c2->id.iOrder;
			pstResetColliders(msr->pst,&reset,sizeof(reset),NULL,NULL);
			smooth.smf.dStart = next.dt;
#ifdef AGGS
			/* do this AFTER reset because update may set SMOOTHACTIVE on other particles */
			msrAggsCollisionUpdate(msr,c1,c2,outDo.iOutcome,outDo.Out,outDo.nOut);
#endif
			switch (msr->param.iCollLogOption) { /* log collision if requested */
			case COLL_LOG_NONE:
				break;
			case COLL_LOG_VERBOSE:
				{
				FILE *fp;
				int i;

				fp = fopen(msr->param.achCollLog,"a");
				assert(fp);
				for (i=0;i<3;i++) {
					c1->r[i] += c1->v[i]*next.dt;
					c2->r[i] += c2->v[i]*next.dt;
					}
				fprintf(fp,"%s-%s COLLISION:T=%e\n",
						_msrParticleLabel(msr,c1->iColor),
						_msrParticleLabel(msr,c2->iColor),dTime + next.dt);
				fprintf(fp,"***1:p=%i,o=%i,i=%i,oi=%i,M=%e,R=%e,dt=%e,rung=%i,"
						"r=(%e,%e,%e),v=(%e,%e,%e),w=(%e,%e,%e)\n",
						c1->id.iPid,c1->id.iOrder,c1->id.iIndex,c1->id.iOrgIdx,
						c1->fMass,c1->fRadius,c1->dt,c1->iRung,
						c1->r[0],c1->r[1],c1->r[2],
						c1->v[0],c1->v[1],c1->v[2],
						c1->w[0],c1->w[1],c1->w[2]);
				fprintf(fp,"***2:p=%i,o=%i,i=%i,oi=%i,M=%e,R=%e,dt=%e,rung=%i,"
						"r=(%e,%e,%e),v=(%e,%e,%e),w=(%e,%e,%e)\n",
						c2->id.iPid,c2->id.iOrder,c2->id.iIndex,c2->id.iOrgIdx,
						c2->fMass,c2->fRadius,c2->dt,c2->iRung,
						c2->r[0],c2->r[1],c2->r[2],
						c2->v[0],c2->v[1],c2->v[2],
						c2->w[0],c2->w[1],c2->w[2]);
				fprintf(fp,"***OUTCOME=%s dT=%e\n",
#ifdef FIX_COLLAPSE
						outDo.iOutcome == MISS ? "MISS" :
#endif
						outDo.iOutcome == MERGE ? "MERGE" :
						outDo.iOutcome == BOUNCE ? "BOUNCE" :
						outDo.iOutcome == FRAG ? "FRAG" : "UNKNOWN",outDo.dT);
				for (i=0;i<outDo.nOut;++i) {
					c = &outDo.Out[i];
					fprintf(fp,"***out%i:p=%i,o=%i,i=%i,oi=%i,M=%e,R=%e,rung=%i,"
							"r=(%e,%e,%e),v=(%e,%e,%e),w=(%e,%e,%e)\n",i,
							c->id.iPid,c->id.iOrder,c->id.iIndex,c->id.iOrgIdx,
							c->fMass,c->fRadius,c->iRung,
							c->r[0],c->r[1],c->r[2],
							c->v[0],c->v[1],c->v[2],
							c->w[0],c->w[1],c->w[2]);
					}
				fclose(fp);
				break;
				}
			case COLL_LOG_TERSE:
				{
				/*
				 ** FORMAT: For each event, time (double), collider 1 iOrder
				 ** (int), collider 2 iOrder (int), number of post-collision
				 ** particles (int), iOrder for each of these (n * int).
				 */

				FILE *fp;
				XDR xdrs;
				double dDum;
				int i;

				if (outDo.iOutcome != MERGE && outDo.iOutcome != FRAG)
					break; /* only care when particle indices change */
				fp = fopen(msr->param.achCollLog,"a");
				assert(fp);
				xdrstdio_create(&xdrs,fp,XDR_ENCODE);
				dDum = dTime + next.dt;
				(void) xdr_double(&xdrs,&dDum);
				(void) xdr_int(&xdrs,&c1->id.iOrgIdx); /*DEBUG not iOrder!*/
				(void) xdr_int(&xdrs,&c2->id.iOrgIdx);
				(void) xdr_int(&xdrs,&outDo.nOut);
				for (i=0;i<outDo.nOut;i++)
					(void) xdr_int(&xdrs,&outDo.Out[i].id.iOrgIdx);
				xdr_destroy(&xdrs);
				(void) fclose(fp);
				break;
				}
			default:
				assert(1); /* invalid collision log option */
				} /* logging */
			} /* if collision */
		} while (COLLISION(next.dt) && smooth.smf.dStart < smooth.smf.dEnd);
	msrAddDelParticles(msr); /* clean up any deletions */
	if (msr->param.nSmooth > msr->N) {
		msr->param.nSmooth = msr->N;
		if (msr->param.bVWarnings)
			printf("WARNING: nSmooth reduced to %i\n",msr->N);
		}
	if (msr->param.bVStep) {
		double dsec = msrTime() - sec;
		printf("%i collision%s: %i miss%s, %i merger%s, %i bounce%s, %i frag%s\n",
			   nCol,nCol==1?"":"s",nMis,nMis==1?"":"es",nMrg,nMrg==1?"":"s",
			   nBnc,nBnc==1?"":"s",nFrg,nFrg==1?"":"s");
		printf("Collision search completed, time = %g sec\n",dsec);
		}
	}

#ifdef OLD_KEPLER
void
msrBuildQQTree(MSR msr,int bActiveOnly,double dMass)
{
	struct inBuildTree in;
	struct outBuildTree out;
	struct inColCells inc;
	KDN *pkdn;
	double sec,dsec;
	int iDum,nCell;

	if (msr->param.bVDetails) printf("Domain Decomposition...\n");
	sec = msrTime();
	pstQQDomainDecomp(msr->pst,NULL,0,NULL,NULL);
	msrMassCheck(msr,dMass,"After pstDomainDecomp in msrBuildTree");
	dsec = msrTime() - sec;
	if (msr->param.bVDetails) {
		printf("Domain Decomposition complete, Wallclock: %f secs\n\n",dsec);
		}
	if (msr->param.bVDetails) printf("Building local trees...\n");
	/*
	 ** First make sure the particles are in Active/Inactive order.
	 */
	msrActiveOrder(msr);
	in.nBucket = msr->param.nBucket;
	msr->bGravityTree = 0;
	msr->iTreeType = MSR_TREE_QQ;
	in.bActiveOnly = bActiveOnly;
	sec = msrTime();
	pstQQBuildTree(msr->pst,&in,sizeof(in),&out,&iDum);
	msrMassCheck(msr,dMass,"After pstBuildTree in msrBuildQQ");
	dsec = msrTime() - sec;
	if (msr->param.bVDetails) {
		printf("Tree built, Wallclock: %f secs\n\n",dsec);
		}
	nCell = 1<<(1+(int)ceil(log((double)msr->nThreads)/log(2.0)));
	pkdn = malloc(nCell*sizeof(KDN));
	assert(pkdn != NULL);
	inc.iCell = ROOT;
	inc.nCell = nCell;
	pstColCells(msr->pst,&inc,sizeof(inc),pkdn,NULL);
	pstDistribCells(msr->pst,pkdn,nCell*sizeof(KDN),NULL,NULL);
	free(pkdn);
	msrMassCheck(msr,dMass,"After pstDistribCells in msrBuildQQ");
	}
#endif

#endif /* COLLISIONS */

#ifdef AGGS

void msrAggsToBodyAxes(MSR msr,int iAggIdx,Aggregate *agg)
{
	/* requires call to msrAggsGetAxes() first */

	struct inAggsToBodyAxes in;

	in.iAggIdx = iAggIdx;
	matrixTranspose(agg->lambda,in.spaceToBody);
	pstAggsToBodyAxes(msr->pst,&in,sizeof(in),NULL,NULL);
	}

void msrAggsGetAxes(MSR msr,int iAggIdx,Aggregate *agg)
{
	/* requires call to msrAggsGetCOM() first */

	struct inAggsGetAxes in;
	struct outAggsGetAxes out;
	Matrix spaceToBody;
	Vector h;
	int k;

	/* compute inertia tensor I and angular momentum L */
	in.iAggIdx = iAggIdx;
	vectorCopy(agg->r_com,in.r_com);
	vectorCopy(agg->v_com,in.v_com);
	pstAggsGetAxes(msr->pst,&in,sizeof(in),&out,NULL);

	/*
	 ** The inertia tensor computed is only valid for the diagonal and
	 ** upper triangle.  Use symmetry to fill out matrix for jacobi().
	 */
	out.I[1][0] = out.I[0][1]; /* row-column notation */
	out.I[2][0] = out.I[0][2];
	out.I[2][1] = out.I[1][2];

	/* compute moments and principal axes (lambda) */
	jacobi(out.I,agg->moments,agg->lambda);

	/* transform L to body frame (not stored) */
	matrixTranspose(agg->lambda,spaceToBody);
	matrixTransform(spaceToBody,out.L,h);

	/* compute spin vector omega = I^{-1} h */
	for (k=0;k<3;k++) /* shortcut for diagonal matrix in body frame */
		agg->omega[k] = h[k]/agg->moments[k];
	}

void msrAggsGetCOM(MSR msr,int iAggIdx,Aggregate *agg)
{
	struct inAggsGetCOM in;
	struct outAggsGetCOM out;
	int k;

	in.iAggIdx = iAggIdx;
	pstAggsGetCOM(msr->pst,&in,sizeof(in),&out,NULL);
	assert(out.m > 0.0);

	for (k=0;k<3;k++) {
		agg->r_com[k] = out.mr[k]/out.m;
		agg->v_com[k] = out.mv[k]/out.m;
		}
	}

void msrAggsUpdate(MSR msr,int iAggID)
{
	Aggregate *agg;

	assert(iAggID >= 0 && iAggID < msr->nAggs);
	agg = &msr->pAggs[iAggID];
	assert(agg->bAssigned);

	/* get center of mass info, needed before axes */
	msrAggsGetCOM(msr,iAggID,agg);
	/* get moments, principal axes, and spin vector */
	msrAggsGetAxes(msr,iAggID,agg);
	/* get particle positions in body frame */
	msrAggsToBodyAxes(msr,iAggID,agg);
	}

void msrAggsGetNewID(MSR msr)
{
	int i;

	/* find next unassigned aggregate; grow buffer if needed */

	for (i=0;i<msr->nAggs;i++)
		if (!msr->pAggs[i].bAssigned) {
			msr->iAggNewID = i;
			return;
			}

	msr->iAggNewID = msr->nAggs;
	msr->nAggs <<= 1; /* double buffer size */
	if (msr->param.bVDetails)
		(void) printf("Doubling aggregate buffer to %i.\n",msr->nAggs);
	msr->pAggs = (Aggregate *) realloc(msr->pAggs,msr->nAggs*sizeof(Aggregate));
	assert(msr->pAggs);
	for (i=msr->iAggNewID;i<msr->nAggs;i++)
		msr->pAggs[i].bAssigned = 0;
	}

void msrAggsFind(MSR msr)
{
	/*
	 ** Finds and initializes aggregates. Must be called immediately
	 ** after initial conditions are loaded.
	 */

	struct outAggsFind outFind;
	struct inAggsConfirm inConfirm;
	struct outAggsConfirm outConfirm;
	Aggregate *agg;
	int i,nAssigned = 0;

	if (msr->param.bVDetails)
		(void) printf("Allocating space for aggregates...\n");

	/*
	 ** We don't know in advance how many aggregates there will be.
	 ** Worse, the aggregate indices could be in any order, and may
	 ** not even be continguous. Strategy: loop through all particles
	 ** to find the largest aggregate index (stored as -1 - org_idx)
	 ** and assume storage for that many aggregates (+ 1) is required.
	 ** Some slots may be left empty though.
	 */

	pstAggsFind(msr->pst,NULL,0,&outFind,NULL);
	if (outFind.iMaxIdx == -1) {
		/* no aggregates found -- make some space anyway */
		msr->nAggs = AGGS_INIT_BUFF_SIZE;
		}
	else {
		/* round up buffer size to next power of 2 */
		msr->nAggs = 1 << (int) (log((outFind.iMaxIdx + 1) << 1)/M_LN2);
		if (msr->nAggs < AGGS_INIT_BUFF_SIZE)
			msr->nAggs = AGGS_INIT_BUFF_SIZE;
		}
	assert(msr->nAggs > 1);

	msr->pAggs = (Aggregate *) malloc(msr->nAggs*sizeof(Aggregate));
	assert(msr->pAggs);

	if (msr->param.bVDetails)
		(void) printf("Space for %i aggregates allocated.\n",msr->nAggs);

	for (i=0;i<msr->nAggs;i++) {
		agg = &msr->pAggs[i];
		agg->bAssigned = 0;
		/* make sure at least one particle belongs to each aggregate */
		if (i <= outFind.iMaxIdx) {
			inConfirm.iAggIdx = i;
			pstAggsConfirm(msr->pst,&inConfirm,sizeof(inConfirm),&outConfirm,NULL);
			agg->bAssigned = outConfirm.bAssigned;
			}
		if (agg->bAssigned) {
			++nAssigned;
			msrAggsUpdate(msr,i);
			}
		}

	msrAggsGetNewID(msr);

	if (msr->param.bVDetails)
		(void) printf("%i aggregate%s initialized.\n",nAssigned,
					  nAssigned == 1 ? "" : "s");

	if (msr->param.bVWarnings && msr->nAggs > AGGS_INIT_BUFF_SIZE &&
		nAssigned < (msr->nAggs >> 2))
		(void) fprintf(stderr,
					   "WARNING: Inefficient use of aggregate buffer\n");
	}

void msrAggsKick(MSR msr,double dt)
{
	/* dt should be 0.5*dDelta */

	struct inAggsSetSpaceVel in;
	Aggregate *agg;
	int i,k;

	for (i=0;i<msr->nAggs;i++) {
		agg = &msr->pAggs[i];
		if (agg->bAssigned) {
			for (k=0;k<3;k++)
				agg->v_com[k] += agg->a_com[k]*dt;
			in.iAggIdx = i;
			vectorCopy(agg->v_com,in.v_com);
			vectorCopy(agg->omega,in.omega);
			matrixCopy(agg->lambda,in.lambda);
			pstAggsSetSpaceVel(msr->pst,&in,sizeof(in),NULL,NULL);
			}
		}
	}

void msrAggsAdvance(MSR msr,Aggregate *agg,double dToTime)
{
	/* Integrate to get new lambda, omega, and COM position */

	FLOAT vars[12];
	FLOAT time = 0.0;			/* irrelevant time */
	double dt;
	int k;

	dt = dToTime - agg->dLastUpdate;

	if (dt == 0.0) return; /* nothing to be done */

	/* Make sure time step is reasonable */

/*DEBUG!! fix collapse uses negative step...
	assert(dt > 0.0);
*/

/*DEBUG ???
	assert(dt < 0.01*2*M_PI/sqrt(agg->omega[0]*agg->omega[0] +
								 agg->omega[1]*agg->omega[1] +
								 agg->omega[2]*agg->omega[2]));
*/

	/*
	 ** Build vars.  To make things easier in aggsEulerDerivs(), we
	 ** copy the lambda data by *columns* (i.e. by principal axes).
	 */
	for (k=0;k<3;k++) {
		vars[k] = agg->omega[k];
		vars[k + 3] = agg->lambda[k][0]; /* q1 */
		vars[k + 6] = agg->lambda[k][1]; /* q2 */
		vars[k + 9] = agg->lambda[k][2]; /* q3 */
		}

	/*
	 ** Take an rk4 step.  Note that the integration time is
	 ** largely irrevelant, hence we pass 0.  The actual time
	 ** in general terms does not matter, due to the structure
	 ** of the Euler equations.
	 **DEBUG note agg passed for moments, right?
	 */
	aggsRungeStep(dt,time,vars,12,agg,aggsEulerDerivs,&time,vars);

	/* Unpack vars */
	for (k=0;k<3;k++) {
		agg->omega[k] = vars[k];
		agg->lambda[k][0] = vars[k + 3];
		agg->lambda[k][1] = vars[k + 6];
		agg->lambda[k][2] = vars[k + 9];
		}
 
	/* COM "drift" */
	for (k=0;k<3;k++)
		agg->r_com[k] += agg->v_com[k]*dt;
	}

void msrAggsAdvanceOpen(MSR msr)
{
	Aggregate *agg;
	int i;

	for (i=0;i<msr->nAggs;i++) {
		agg = &msr->pAggs[i];
		if (agg->bAssigned)
			agg->dLastUpdate = 0.0;
		}
	}

void msrAggsAdvanceClose(MSR msr,double dt)
{
	/* dt should be dDelta */

	struct inAggsSetSpacePos in;
	Aggregate *agg;
	int i;

	for (i=0;i<msr->nAggs;i++) {
		agg = &msr->pAggs[i];
		if (agg->bAssigned) {
			msrAggsAdvance(msr,agg,dt);
			/*DEBUG following could be separate msr routine...*/
			in.iAggIdx = i;
			vectorCopy(agg->r_com,in.r_com);
			matrixCopy(agg->lambda,in.lambda);
			pstAggsSetSpacePos(msr->pst,&in,sizeof(in),NULL,NULL);
			}
		}
	}

void msrAggsAdvanceCollider(MSR msr,COLLIDER *c,double dToTime)
{
	Aggregate *agg = &msr->pAggs[COLLIDER_AGG_ID(c)];

	msrAggsAdvance(msr,agg,dToTime);

	vectorCopy(agg->r_com,c->agg.r_com);
	vectorCopy(agg->v_com,c->agg.v_com);
	vectorCopy(agg->omega,c->agg.omega);
	vectorCopy(agg->moments,c->agg.moments);
	matrixCopy(agg->lambda,c->agg.lambda);
	}

void msrAggsCollisionUpdate(MSR msr,COLLIDER *c1,COLLIDER *c2,
							int iOutcome,COLLIDER c[],int nOut)
{
	struct inAggsMerge inMerge;
	int iAggID;

	assert(iOutcome == MERGE && nOut == 1); /* for now */

	/* currently output not used */

	if (iOutcome == MERGE) {
		if (COLLIDER_IS_AGG(c1) && COLLIDER_IS_AGG(c2)) {
			assert(COLLIDER_AGG_ID(c1) != COLLIDER_AGG_ID(c2));
			if (COLLIDER_AGG_ID(c1) < COLLIDER_AGG_ID(c2)) {
				inMerge.iOldID = COLLIDER_AGG_ID(c2);
				inMerge.iNewID = COLLIDER_AGG_ID(c1);
				}
			else { /* (COLLIDER_AGG_ID(c2) < COLLIDER_AGG_ID(c1)) */
				inMerge.iOldID = COLLIDER_AGG_ID(c1);
				inMerge.iNewID = COLLIDER_AGG_ID(c2);
				}
			pstAggsMerge(msr->pst,&inMerge,sizeof(inMerge),NULL,NULL);
			msr->pAggs[inMerge.iOldID].bAssigned = 0;
			iAggID = inMerge.iNewID;
			}
		else if (COLLIDER_IS_AGG(c1) && !COLLIDER_IS_AGG(c2)) {
			iAggID = COLLIDER_AGG_ID(c1);
			}
		else if (COLLIDER_IS_AGG(c2) && !COLLIDER_IS_AGG(c1)) {
			iAggID = COLLIDER_AGG_ID(c2);
			}
		else { /* (!COLLIDER_IS_AGG(c1) && !COLLIDER_IS_AGG(c2)) */
			iAggID = msr->iAggNewID;
			msr->pAggs[iAggID].bAssigned = 1;
			msrAggsGetNewID(msr);
			}
		msrAggsUpdate(msr,iAggID);
		}
	}

void msrAggsGravity(MSR msr)
{
	struct inAggsAccel inAccel;
	struct outAggsAccel outAccel;
	struct inAggsTorque inTorque;
	struct outAggsTorque outTorque;
	Matrix spaceToBody;
	Aggregate *agg;
	int i,k;

	for (i=0;i<msr->nAggs;i++) {
		agg = &msr->pAggs[i];
		if (agg->bAssigned) {
			vectorZero(agg->a_com);
			vectorZero(agg->torque);
			if (msr->param.bDoGravity) {
				/* get acceleration of center of mass */
				inAccel.iAggIdx = i;
				pstAggsAccel(msr->pst,&inAccel,sizeof(inAccel),&outAccel,NULL);
				assert(outAccel.m > 0.0);
				for (k=0;k<3;k++) agg->a_com[k] = outAccel.ma[k]/outAccel.m;
				if (msr->param.bDoSelfGravity) {
					/* get torque around center of mass */
					inTorque.iAggIdx = i;
					for (k=0;k<3;k++) {
						inTorque.r_com[k] = agg->r_com[k];
						inTorque.a_com[k] = agg->a_com[k];
						}
					pstAggsTorque(msr->pst,&inTorque,sizeof(inTorque),&outTorque,NULL);
					/* transform torque vector to body frame */
					matrixTranspose(agg->lambda,spaceToBody);
					matrixTransform(spaceToBody,outTorque.torque,agg->torque);
					}
				}
			}
		}
	}

void msrAggsActivate(MSR msr)
{
	pstAggsActivate(msr->pst,NULL,0,NULL,NULL);
	}

void msrAggsDeactivate(MSR msr)
{
	pstAggsDeactivate(msr->pst,NULL,0,NULL,NULL);
	}

#endif /* AGGS */


/* Note if dDataOut is NULL it just counts the number of valid input lines */
int msrReadASCII(MSR msr, char *extension, int nDataPerLine, double *dDataOut)
{
	char achFile[PST_FILENAME_SIZE];
	char ach[PST_FILENAME_SIZE];
	LCL *plcl = &msr->lcl;
	FILE *fp;
	int i,ret;
	char achIn[160];
	double *dData;

	if (dDataOut == NULL) 
		dData = malloc(sizeof(double)*nDataPerLine);
	else
		dData = dDataOut;
	
	assert(nDataPerLine > 0 && nDataPerLine <= 10);
	/*
	 ** Add Data Subpath for local and non-local names.
	 */
	achFile[0] = 0;
	sprintf(achFile,"%s/%s.%s",msr->param.achDataSubPath,
			msr->param.achOutName, extension);
	/*
	 ** Add local Data Path.
	 */
	if (plcl->pszDataPath) {
		strcpy(ach,achFile);
		sprintf(achFile,"%s/%s",plcl->pszDataPath,ach);
		}
	fp = fopen(achFile,"r");
	if (!fp) {
		if (msr->param.bVWarnings)
			printf("WARNING: Could not open .%s input file:%s\n",extension,achFile);
		return 0;
		}

	i = 0;
	while (1) {
		if (!fgets(achIn,160,fp)) goto Done;
		switch (nDataPerLine) {
		case 1:
			ret = sscanf(achIn,"%lf",dData); 
			break;
		case 2:
			ret = sscanf(achIn,"%lf %lf",dData,dData+1); 
			break;
		case 3:
			ret = sscanf(achIn,"%lf %lf %lf",dData,dData+1,dData+2); 
			break;
		case 4:
			ret = sscanf(achIn,"%lf %lf %lf %lf",dData,dData+1,dData+2,dData+3); 
			break;
		case 5:
			ret = sscanf(achIn,"%lf %lf %lf %lf %lf",dData,dData+1,dData+2,dData+3,dData+4); 
			break;
		case 6:
			ret = sscanf(achIn,"%lf %lf %lf %lf %lf %lf",dData,dData+1,dData+2,dData+3,dData+4,dData+5); 
			break;
		case 7:
			ret = sscanf(achIn,"%lf %lf %lf %lf %lf %lf %lf",
						 dData,dData+1,dData+2,dData+3,dData+4,dData+5,dData+6); 
			break;
		case 8:
			ret = sscanf(achIn,"%lf %lf %lf %lf %lf %lf %lf %lf",
						 dData,dData+1,dData+2,dData+3,dData+4,dData+5,dData+6,dData+7); 
			break;
		case 9:
			ret = sscanf(achIn,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",
						 dData,dData+1,dData+2,dData+3,dData+4,dData+5,dData+6,dData+7,dData+8); 
			break;
		case 10:
			ret = sscanf(achIn,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
						 dData,dData+1,dData+2,dData+3,dData+4,dData+5,dData+6,dData+7,dData+8,dData+9); 
			break;
		default:
			ret = EOF;
			assert(0);
			}
		if (ret != nDataPerLine) goto Done;
		++i;
		if (dDataOut != NULL) dData += nDataPerLine;
		}
 Done:
	fclose(fp);
	if (dDataOut != NULL && msr->param.bVDetails) printf("Read %i lines from %s\n",i,achFile);
	if (dDataOut == NULL) free(dData);
	return i;
	}

int msrSetTypeFromFile(MSR msr, char *file, int iSetMask)
{
	FILE *fp;
	struct inSetTypeFromFile in;
	struct outSetTypeFromFile out;

	in.iSetMask = iSetMask;
	in.biGasOrder = 1; /* Set from Parent Gas iOrder for stars? */
	assert(strlen(file) < PST_SETTYPEFROMFILEMAXLEN);
	strcpy( &(in.file[0]), file );

	fp = fopen( file, "r" );
	if (!fp) {
	  fprintf(stderr,"ERROR: Could not open iOrder list file:%s\n",file);
	  assert(0);
	  }
	fclose(fp);

	pstSetTypeFromFile(msr->pst,&in,sizeof(in),&out,NULL);
	
	if (msr->param.bVDetails) printf("%d iOrder numbers read.  %d direct and %d Gas Parent iOrder photogenic particles selected.",out.niOrder,out.nSet,out.nSetiGasOrder);
	
	return out.nSet+out.nSetiGasOrder;
	}
