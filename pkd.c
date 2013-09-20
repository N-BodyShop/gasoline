#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <sys/time.h>

#ifdef CRAY_XT3
#include "../xdr/types.h"
#include "../xdr/xdr.h"
#else
#include <rpc/types.h>
#include <rpc/xdr.h>
#endif

#include "pkd.h"
#include "ewald.h"
#include "grav.h"
#include "walk.h"
#include "opentype.h"
#include "mdl.h"
#include "tipsydefs.h"

#ifdef COLLISIONS
#include "ssdefs.h" /* in turn includes ssio.h */
#include "collision.h"
#endif

#ifdef AGGS
#include "aggs.h"
#endif
#include "debug.h"

#ifdef SLIDING_PATCH
#include <sys/types.h> /* for getpid() */
#include <unistd.h> /* ditto */
#include <time.h> /* for time() */
#ifndef SGN
int SGN(double x);
#define SGN(x) ((x) < 0.0 ? (-1) : (x) > 0.0 ? 1 : 0)
#endif
#endif

#ifdef _LARGE_FILES
#define fseek fseeko
#endif


double pkdGetTimer(PKD pkd,int iTimer)
{
	return(pkd->ti[iTimer].sec);
	}

double pkdGetSystemTimer(PKD pkd,int iTimer)
{
	return(pkd->ti[iTimer].system_sec);
	}

double pkdGetWallClockTimer(PKD pkd,int iTimer)
{
	return(pkd->ti[iTimer].wallclock_sec);
	}


void pkdClearTimer(PKD pkd,int iTimer)
{
	int i;

	if (iTimer >= 0) {
		pkd->ti[iTimer].sec = 0.0;
		pkd->ti[iTimer].system_sec = 0.0;
		pkd->ti[iTimer].wallclock_sec = 0.0;
		pkd->ti[iTimer].iActive = 0;
		}
	else {
		for (i=0;i<MAX_TIMERS;++i) {
			pkd->ti[i].sec = 0.0;
			pkd->ti[i].system_sec = 0.0;
			pkd->ti[i].wallclock_sec = 0.0;
			pkd->ti[i].iActive = 0;
			}
		}
	}


void pkdStartTimer(PKD pkd,int iTimer)
{
	struct timezone tz;
	struct timeval tv;
	tz.tz_minuteswest = 0;
	tz.tz_dsttime = 0;

	pkd->ti[iTimer].iActive++;

	if (pkd->ti[iTimer].iActive == 1) {
		pkd->ti[iTimer].stamp = mdlCpuTimer(pkd->mdl);
		gettimeofday(&tv,&tz);
		pkd->ti[iTimer].wallclock_stamp = tv.tv_sec + 1e-6*(double) tv.tv_usec;
#ifndef _CRAYMPP
		{
	    struct rusage ru;
	    
	    getrusage(0,&ru);
	    pkd->ti[iTimer].system_stamp = (double)ru.ru_stime.tv_sec + 1e-6*(double)ru.ru_stime.tv_usec;
		}
#endif
		}
	}


void pkdStopTimer(PKD pkd,int iTimer)
{
	double sec;
	struct timeval tv;
	struct timezone tz;

	sec = -pkd->ti[iTimer].stamp;
	pkd->ti[iTimer].stamp = mdlCpuTimer(pkd->mdl);
	sec += pkd->ti[iTimer].stamp;
	if (sec < 0.0) sec = 0.0;
	pkd->ti[iTimer].sec += sec;

	sec = -pkd->ti[iTimer].wallclock_stamp;
	tz.tz_minuteswest = 0;
	tz.tz_dsttime = 0;
	gettimeofday( &tv, &tz );
	pkd->ti[iTimer].wallclock_stamp = tv.tv_sec + 1e-6*(double)tv.tv_usec;
	sec += pkd->ti[iTimer].wallclock_stamp;
	if (sec < 0.0) sec = 0.0;
	pkd->ti[iTimer].wallclock_sec += sec;

#ifndef _CRAYMPP
	{
	struct rusage ru;

	sec = -pkd->ti[iTimer].system_stamp;
	getrusage(0,&ru);
	pkd->ti[iTimer].system_stamp = ((double)ru.ru_stime.tv_sec + 1e-6*(double)ru.ru_stime.tv_usec);
	sec += pkd->ti[iTimer].system_stamp;
	if (sec < 0.0) sec = 0.0;
	pkd->ti[iTimer].system_sec += sec;
	}
#endif
	pkd->ti[iTimer].iActive--;
	}


void pkdInitialize(PKD *ppkd,MDL mdl,int iOrder,int nStore,int nLvl,
				   FLOAT *fPeriod,FLOAT dxInflow,FLOAT dxOutflow,int nDark,int nGas,int nStar)
{
	PKD pkd;
	int j;
	
	pkd = (PKD)malloc(sizeof(struct pkdContext));
	mdlassert(mdl,pkd != NULL);
	pkd->mdl = mdl;
	pkd->iOrder = iOrder;
	pkd->idSelf = mdlSelf(mdl);
	pkd->nThreads = mdlThreads(mdl);
	pkd->nStore = nStore;
	pkd->nLocal = 0;
	pkd->nDark = nDark;
	pkd->nGas = nGas;
	pkd->nStar = nStar;
	pkd->nMaxOrderGas = NIORDERGASBUFFER + nGas - 1; 
	pkd->nMaxOrderDark = NIORDERGASBUFFER + nGas + nDark - 1;
	pkd->nRejects = 0;
	for (j=0;j<3;++j) {
		pkd->fPeriod[j] = fPeriod[j];
		}
	pkd->dxInflow = dxInflow;
	pkd->dxOutflow = dxOutflow;
	pkd->sinkLog.nLog = 0; /* Nothing in log, even if not needed */
	/*
	 ** Allocate the main particle store.
	 ** Need to use mdlMalloc() since the particles will need to be
	 ** visible to all other processors thru mdlAquire() later on.
	 **
	 ** We need one EXTRA storage location at the very end to use for 
	 ** calculating acceleration on arbitrary positions in space, for example
	 ** determining the force on the sun. The easiest way to do this is to
	 ** allocate one hidden particle, which won't interfere with the rest of
	 ** the code (hopefully). pkd->pStore[pkd->nStore] is this particle.
	 */
	pkd->pStore = mdlMalloc(pkd->mdl,(nStore+1)*sizeof(PARTICLE));
	mdlassert(pkd->mdl,pkd->pStore != NULL);
	pkd->kdNodes = NULL;
	pkd->piLeaf = NULL;
	pkd->kdTop = NULL;
	/*
	 ** Allocate initial interaction lists
	 */
	pkd->nMaxPart = 500;
	pkd->nMaxCellSoft = 500;
	pkd->nMaxCellNewt = 500;
	pkd->nSqrtTmp = 500;
	pkd->ilp = malloc(pkd->nMaxPart*sizeof(ILP));
	mdlassert(pkd->mdl,pkd->ilp != NULL);
	pkd->ilcs = malloc(pkd->nMaxCellSoft*sizeof(ILCS));
	mdlassert(pkd->mdl,pkd->ilcs != NULL);
	pkd->ilcn = malloc(pkd->nMaxCellNewt*sizeof(ILCN));
	mdlassert(pkd->mdl,pkd->ilcn != NULL);
	pkd->sqrttmp = malloc(pkd->nSqrtTmp*sizeof(double));
	mdlassert(pkd->mdl,pkd->sqrttmp != NULL);
	pkd->d2a = malloc(pkd->nSqrtTmp*sizeof(double));
	mdlassert(pkd->mdl,pkd->d2a != NULL);
	/*
	 ** Ewald stuff!
	 */
	pkd->nMaxEwhLoop = 100;
	pkd->ewt = malloc(pkd->nMaxEwhLoop*sizeof(EWT));
	mdlassert(pkd->mdl,pkd->ewt != NULL);
	/*
	 * Cooling
	 */
#ifdef GASOLINE
#ifndef NOCOOLING
	pkd->Cool = CoolInit();
#endif	
#endif	

	*ppkd = pkd;


	}


void pkdFinish(PKD pkd)
{
	if (pkd->kdNodes) {
		/*
		 ** Close caching space and free up nodes.
		 */
		mdlFinishCache(pkd->mdl,CID_CELL);
		mdlFree(pkd->mdl,pkd->kdNodes);
		}
	if (pkd->kdTop) free(pkd->kdTop);
	if (pkd->piLeaf) free(pkd->piLeaf);
	free(pkd->ilp);
	free(pkd->ilcs);
	free(pkd->ilcn);
	free(pkd->sqrttmp);
	free(pkd->d2a);
	free(pkd->ewt);
#if defined(GASOLINE) && !defined(NOCOOLING)
	CoolFinalize(pkd->Cool);
#endif
	mdlFree(pkd->mdl,pkd->pStore);
	free(pkd);
	}


void pkdSeek(PKD pkd,FILE *fp,int nStart,int bStandard)
{
	long lStart;

	/*
	 ** Seek according to true XDR size structures when bStandard is true.
	 ** This may be a bit dicey, but it should work as long
	 ** as no one changes the tipsy binary format!
	 */
	if (bStandard) lStart = 32;
	else lStart = sizeof(struct dump);
	if (nStart > pkd->nGas) {
		if (bStandard) lStart += pkd->nGas*48L;
		else lStart += pkd->nGas*sizeof(struct gas_particle);
		nStart -= pkd->nGas;
		if (nStart > pkd->nDark) {
			if (bStandard) lStart += pkd->nDark*36L;
			else lStart += pkd->nDark*sizeof(struct dark_particle);
			nStart -= pkd->nDark;
			if (bStandard) lStart += nStart*44L;
			else lStart += nStart*sizeof(struct star_particle);
			}
		else {
			if (bStandard) lStart += nStart*36L;
			else lStart += nStart*sizeof(struct dark_particle);
			}
		}
	else {
		if (bStandard) lStart += nStart*48L;
		else lStart += nStart*sizeof(struct gas_particle);
		}
	assert(fseek(fp,lStart,0) == 0);
	}


void pkdGenericSeek(PKD pkd,FILE *fp,long nStart,int iHeader,int iElement)
{
  long lStart;
  
  lStart = iHeader + nStart*((long)iElement);
  fseek(fp,lStart,0);
}


void pkdReadTipsy(PKD pkd,char *pszFileName,int nStart,int nLocal,
		  int bStandard,int iReadIOrder,double dvFac,double dTuFac)
    {
    FILE *fp,*fpiord = NULL;
#if defined(SIMPLESF) || defined(STARFORM)
    FILE *fpmStar = NULL, *fptCoolAgain = NULL;
#endif
    int i,j, iSetMask;
    PARTICLE *p;
    struct dark_particle dp;
    struct gas_particle gp;
    struct star_particle sp;
    float fTmp;
    
    pkd->nLocal = nLocal;
    pkd->nActive = nLocal;
    /*
    ** General initialization.
    */
    for (i=0;i<nLocal;++i) {
        p = &pkd->pStore[i];
        TYPEClear(p);
        p->iRung = 0;
        p->fWeight = 1.0;
        p->fDensity = 0.0;
#ifdef DRHODT
        p->fDensity_t = 0.0;
        p->fDensity_PdV = 0.0;
        p->fDensity_PdVcorr = 0.0;
#endif
        p->fBall2 = 0.0;
        p->fBallMax = 0.0;
#ifdef GASOLINE
        p->dt = FLT_MAX;
        p->dtNew = FLT_MAX;
        p->u = 0.0;
        p->uPred = 0.0;
#ifdef MASSNONCOOL
        p->fMassNoncool = 0;
#endif
#ifdef UNONCOOL
        p->uNoncool = 0.;
        p->uNoncoolPred = 0.;
        p->uNoncoolDot = 0.;
#endif
#ifdef STARSINK
        SINK_Lx(p) = 0.0;
        SINK_Ly(p) = 0.0;
        SINK_Lz(p) = 0.0;
#endif
#ifdef STARFORM
        p->uDotFB = 0.0;
        p->fNSN = 0.0;
        p->fNSNtot = 0.0;
        p->fMOxygenOut = 0.0;
        p->fMIronOut = 0.0;
        p->fMFracOxygen = 0.0;
        p->fMFracIron = 0.0;
        p->fTimeCoolIsOffUntil = 0.0;
#endif
#ifdef SIMPLESF
        p->fMassStar = 0;
#endif
#ifdef SHOCKTRACK
        p->ShockTracker = 0.0;
#endif
#ifdef VARALPHA
/*
  p->alpha = 1.0;
  p->alphaPred = 1.0;
*/
        p->alpha = ALPHAMIN;
        p->alphaPred = ALPHAMIN;
        p->divv = 0;
#ifdef DODVDS
        p->dvds = 0;
#endif
#endif
#ifndef NOCOOLING		
        /* Place holders -- later fixed in pkdInitEnergy */
        CoolDefaultParticleData( &p->CoolParticle );
#endif
        p->c = 0.0;
        p->fMetals = 0.0;
        p->fTimeForm = 1e37;
#ifdef SINKING
        p->fTrueMass = 0;
        p->fSinkingTime = 1e37;  
        p->iSinkingOnto = -1;
#endif
#endif
#ifdef NEED_VPRED
        for (j=0;j<3;++j) {
            p->vPred[j] = 0.0;
            }
#endif
        }
    /*
    ** Seek past the header and up to nStart.
    */
    fp = fopen(pszFileName,"r");
    mdlassert(pkd->mdl,fp != NULL);
    /*
    ** Seek to right place in file.
    */
    pkdSeek(pkd,fp,nStart,bStandard);
    
    /* Open iOrder file if requested */
    if (iReadIOrder) {
        char atmp[512];
        sprintf(atmp,"%s.iord",pszFileName);
        fpiord = fopen(atmp,"r");   
        mdlassert(pkd->mdl,fpiord != NULL);
        /*
        ** Seek to right place in file
        */
        switch(iReadIOrder) {
        case 1:
            pkdGenericSeek(pkd,fpiord,nStart,sizeof(int),sizeof(int));
            break;
        case 2:
            pkdGenericSeek(pkd,fpiord,nStart,sizeof(int),sizeof(long long));
            if (bStandard) assert(sizeof(long)==sizeof(long long));
            break;
        case 3:
            pkdGenericSeek(pkd,fpiord,nStart,sizeof(int),sizeof(pkd->pStore[0].iOrder));
            if (bStandard) assert(sizeof(int)==sizeof(pkd->pStore[0].iOrder));
            break;
        default:
            fprintf(stderr,"Don't understand iOrder format: %d\n",iReadIOrder);
            mdlassert(pkd->mdl,0);
            }
        }
  
#ifdef STARFORM
        {
        char atmp[512];
        sprintf(atmp,"%s.coolontime",pszFileName);
        fptCoolAgain = fopen(atmp,"r");
        if (fptCoolAgain!=NULL) {
            /*
            ** Seek to right place in file
            */
            pkdGenericSeek(pkd,fptCoolAgain,nStart,sizeof(int),sizeof(float));
        }
        else {
            if(pkd->idSelf == 0)
                fprintf(stderr, "Could not open %s,  skipped.\n",atmp);
        }
        }
#endif

#ifdef SIMPLESF
        {
        char atmp[512];
        sprintf(atmp,"%s.tCoolAgain",pszFileName);
        fptCoolAgain = fopen(atmp,"r");
        if (fptCoolAgain!=NULL) {
            /*
            ** Seek to right place in file
            */
            pkdGenericSeek(pkd,fptCoolAgain,nStart,sizeof(int),sizeof(float));
        }
        else {
            fprintf(stderr, "Could not open %s,  skipped.\n",atmp);
        }
        
        sprintf(atmp,"%s.mStar",pszFileName);
        fpmStar = fopen(atmp,"r");
        if (fpmStar!=NULL) {
            /*
            ** Seek to right place in file
            */
            pkdGenericSeek(pkd,fpmStar,nStart,sizeof(int),sizeof(float));
        }
        else {
            fprintf(stderr, "Could not open %s,  skipped.\n",atmp);
        }
        }
#endif

	/*
	 ** Read Stuff!
	 */
        
	if (bStandard) {
		FLOAT vTemp;
		long LongTmp;
		int IntTmp;
		XDR xdrs,xdrstoc,xdrsiord;
		xdrstdio_create(&xdrs,fp,XDR_DECODE);
#ifdef STARFORM
		xdrstdio_create(&xdrstoc,fptCoolAgain,XDR_DECODE);
#endif
		if (iReadIOrder) xdrstdio_create(&xdrsiord,fpiord,XDR_DECODE);
        
		for (i=0;i<nLocal;++i) {
			p = &pkd->pStore[i];
			p->iOrder = nStart + i; /* temporary */
#if NIORDERGASBUFFER
            if (p->iOrder >= pkd->nGas) p->iOrder += NIORDERGASBUFFER;
#endif            
			if (pkdIsGasByOrder(pkd,p)) {
                iSetMask = TYPE_GAS; /* saves identity based on Tipsy file in case iOrder changed */
				xdr_float(&xdrs,&fTmp);
				p->fMass = fTmp;
				assert(p->fMass > 0.0);
#ifdef SINKING
				p->fTrueMass = fTmp;
#endif
				for (j=0;j<3;++j) {
					xdr_float(&xdrs,&fTmp);
					p->r[j] = fTmp;
					}
				for (j=0;j<3;++j) {
					xdr_float(&xdrs,&fTmp);
					vTemp = fTmp;
					p->v[j] = dvFac*vTemp;			
#ifdef NEED_VPRED
					p->vPred[j] = dvFac*vTemp;
#endif
                    }
#ifdef GASOLINE
				xdr_float(&xdrs,&fTmp);
				p->fDensity = fTmp;
				/*
				** Convert Temperature to Thermal energy.
				*/
				xdr_float(&xdrs,&fTmp);
				vTemp = fTmp;
#ifdef FBPARTICLETSTART
                if (fTmp > FBPARTICLETSTART) TYPESet(p, TYPE_FEEDBACK);
#endif
				p->u = dTuFac*vTemp;
				p->uPred = dTuFac*vTemp;
// Special purpose hack for testing noncooling
#ifdef UNONCOOLINIT
				p->uNoncool = p->u;
				p->uNoncoolPred = p->u;
                p->u = 0;
                p->uPred = 0;
#endif
#ifdef COOLDEBUG
				assert(p->u >= 0.0);
				assert(p->uPred >= 0.0);
#endif
				xdr_float(&xdrs,&fTmp);
				p->fSoft = fTmp;
#ifdef CHANGESOFT
				p->fSoft0 = fTmp;
#endif
				xdr_float(&xdrs,&fTmp);
				p->fMetals = fTmp;
#ifdef DIFFUSION
				p->fMetalsPred = fTmp;
#ifdef MASSDIFF
				p->fMass0 = p->fMass;
#endif
#endif				
#ifdef STARFORM
				/* O and Fe ratio based on Asplund et al 2009 */
				if (p->fMetals && !p->fMFracOxygen && 
				    !p->fMFracIron) {
                    p->fMFracOxygen = 0.43 * p->fMetals;
                    p->fMFracIron = 0.098 * p->fMetals;
                    }
#endif
#else /* now not GASOLINE */
				xdr_float(&xdrs,&fTmp); /* Dens */
				xdr_float(&xdrs,&fTmp); /* T */
				xdr_float(&xdrs,&fTmp); /* eps */
				p->fSoft = fTmp;
#ifdef CHANGESOFT
				p->fSoft0 = fTmp;
#endif
				xdr_float(&xdrs,&fTmp); /* metals */
#endif
				xdr_float(&xdrs,&fTmp);
				p->fPot = fTmp;
                }
			else if (pkdIsDarkByOrder(pkd,p)) {
                iSetMask = TYPE_DARK;
                xdr_float(&xdrs,&fTmp);
                p->fMass = fTmp;
                assert(p->fMass >= 0.0);
#ifdef SINKING
                p->fTrueMass = fTmp;
#endif
                for (j=0;j<3;++j) {
                    xdr_float(&xdrs,&fTmp);
                    p->r[j] = fTmp;
                    }
                for (j=0;j<3;++j) {
                    xdr_float(&xdrs,&fTmp);
                    vTemp = fTmp;
                    p->v[j] = dvFac*vTemp;			
                    }
                xdr_float(&xdrs,&fTmp);
                p->fSoft = fTmp;
#ifdef CHANGESOFT				
                p->fSoft0 = fTmp;
#endif
                xdr_float(&xdrs,&fTmp);
                p->fPot = fTmp;
                }
            else if (pkdIsStarByOrder(pkd,p)) {
                iSetMask = TYPE_STAR;
                xdr_float(&xdrs,&fTmp);
                p->fMass = fTmp;
#ifdef STARFORM
                p->fMassForm = fTmp;
#endif
                assert(p->fMass >= 0.0);
#ifdef SINKING
                p->fTrueMass = fTmp;
#endif
                for (j=0;j<3;++j) {
                    xdr_float(&xdrs,&fTmp);
                    p->r[j] = fTmp;
                    }
                for (j=0;j<3;++j) {
                    xdr_float(&xdrs,&fTmp);
                    vTemp = fTmp;
                    p->v[j] = dvFac*vTemp;			
                    }
#ifdef GASOLINE
                xdr_float(&xdrs,&fTmp);
                p->fMetals = fTmp;
#ifdef DIFFUSION
                p->fMetalsPred = fTmp;
#ifdef MASSDIFF
                p->fMass0 = p->fMass;
#endif
#endif				
#ifdef STARFORM
                /* O and Fe ratio based on Asplund et al 2009 */
                if (p->fMetals && !p->fMFracOxygen && 
                    !p->fMFracIron) {
                    p->fMFracOxygen = 0.43 * p->fMetals;
                    p->fMFracIron = 0.098 * p->fMetals;
                    }
#endif
                xdr_float(&xdrs,&fTmp);
                p->fTimeForm = fTmp;
#else /* not GASOLINE */
                xdr_float(&xdrs,&fTmp);
                xdr_float(&xdrs,&fTmp);
#endif
                xdr_float(&xdrs,&fTmp);
                p->fSoft = fTmp;
#ifdef CHANGESOFT
                p->fSoft0 = fTmp;
#endif
                xdr_float(&xdrs,&fTmp);
                p->fPot = fTmp;
                }
            else mdlassert(pkd->mdl,0); /* unrecognized type */
            /* particle read from tipsy */
        
#ifdef STARFORM
            if (fptCoolAgain!=NULL) {
                xdr_float(&xdrstoc,&fTmp);
                if (pkdIsGasByOrder(pkd,p)) p->fTimeCoolIsOffUntil = fTmp;
                }
#endif
            
#ifdef SIMPLESF
            if (fptCoolAgain!=NULL) {
                fread(&fTmp,sizeof(float),1,fptCoolAgain);
                if (pkdIsGasByOrder(pkd,p)) p->fTimeForm = fTmp;
                }
            if (fpmStar!=NULL) {
                fread(&fTmp,sizeof(float),1,fpmStar);
                if (pkdIsGasByOrder(pkd,p)) p->fMassStar = fTmp;
                }
#endif
#ifdef INFLOWOUTFLOW
            if (p->r[0] < pkd->dxInflow) { TYPESet(p,TYPE_INFLOW); assert(p->v[0] > 0); }
            if (p->r[0] > pkd->dxOutflow) { TYPESet(p,TYPE_OUTFLOW); assert(p->v[0] > 0); }
#endif
            TYPESet(p,iSetMask);
            /* Read iOrder last so byOrder Types not messed up */
            switch (iReadIOrder) {
            case 0:
                break;
            case 1:
                xdr_int(&xdrsiord,&IntTmp);
//			    fread(&IntTmp,sizeof(IntTmp),1,fpiord);
                p->iOrder = IntTmp;
                break;
            case 2:
                xdr_long(&xdrsiord,&LongTmp);
//			    fread(&LongTmp,sizeof(LongTmp),1,fpiord);
                p->iOrder = LongTmp;
                break;
            case 3:
                xdr_int(&xdrsiord,&IntTmp);
// see assert above -- I have to assume iOrder is int
                p->iOrder = IntTmp;
                break;
                }
            }
        xdr_destroy(&xdrs);
#ifdef STARFORM
        xdr_destroy(&xdrstoc);
#endif
        if (iReadIOrder) xdr_destroy(&xdrsiord);
        } /* standard read done */
    
    /* native format read */
	else { 
		long long LongTmp;
		int IntTmp;
		for (i=0;i<nLocal;++i) {
			p = &pkd->pStore[i];
			p->iOrder = nStart + i;
#if NIORDERGASBUFFER
            if (p->iOrder >= pkd->nGas) p->iOrder += NIORDERGASBUFFER;
#endif            
			if (pkdIsGasByOrder(pkd,p)) {
				iSetMask = TYPE_GAS;
				fread(&gp,sizeof(struct gas_particle),1,fp);
				for (j=0;j<3;++j) {
					p->r[j] = gp.pos[j];
					p->v[j] = dvFac*gp.vel[j];
#ifdef NEED_VPRED
                    p->vPred[j] = dvFac*gp.vel[j];
#endif
                    }
                p->fMass = gp.mass;
                assert(p->fMass >= 0.0);
#ifdef SINKING
                p->fTrueMass = gp.mass;
#endif
                p->fSoft = gp.hsmooth;
#ifdef CHANGESOFT
                p->fSoft0 = gp.hsmooth;
#endif
                p->fPot = gp.phi;
#ifdef GASOLINE
                p->fDensity = gp.rho;
                p->u = dTuFac*gp.temp;
                p->uPred = dTuFac*gp.temp;
#ifdef COOLDEBUG
                assert(p->u >= 0.0);
                assert(p->uPred >= 0.0);
#endif
                p->fMetals = gp.metals;
#ifdef DIFFUSION
                p->fMetalsPred = gp.metals;
#ifdef MASSDIFF
                p->fMass0 = p->fMass;
#endif
#endif				
#ifdef STARFORM
                /* O and Fe ratio based on Asplund et al 2009 */
                if (p->fMetals && !p->fMFracOxygen && !p->fMFracIron) {
                    p->fMFracOxygen = 0.43 * p->fMetals;
                    p->fMFracIron = 0.098 * p->fMetals;
                    }
#endif
#endif
				}
			else if (pkdIsDarkByOrder(pkd,p)) {
                iSetMask = TYPE_DARK;
                fread(&dp,sizeof(struct dark_particle),1,fp);
                for (j=0;j<3;++j) {
                    p->r[j] = dp.pos[j];
                    p->v[j] = dvFac*dp.vel[j];
                    }
                p->fMass = dp.mass;
                assert(p->fMass >= 0.0);
#ifdef SINKING
                p->fTrueMass = dp.mass;
#endif
                p->fSoft = dp.eps;
#ifdef CHANGESOFT
                p->fSoft0 = dp.eps;
#endif
                p->fPot = dp.phi;
                }
            else if (pkdIsStarByOrder(pkd,p)) {
                iSetMask = TYPE_STAR;
                fread(&sp,sizeof(struct star_particle),1,fp);
                for (j=0;j<3;++j) {
                    p->r[j] = sp.pos[j];
                    p->v[j] = dvFac*sp.vel[j];
                    }
                p->fMass = sp.mass;
#ifdef STARFORM
                p->fMassForm = sp.mass;
#endif
                assert(p->fMass >= 0.0);
#ifdef SINKING
                p->fTrueMass = sp.mass;
#endif
                p->fSoft = sp.eps;
#ifdef CHANGESOFT
                p->fSoft0 = sp.eps;
#endif
                p->fPot = sp.phi;
#ifdef GASOLINE
                p->fMetals = sp.metals;
                p->fTimeForm = sp.tform;		
#ifdef STARFORM
                /* O and Fe ratio based on Asplund et al 2009 */
                if (p->fMetals && !p->fMFracOxygen
                    && !p->fMFracIron) {
                    p->fMFracOxygen = 0.43 * p->fMetals;
                    p->fMFracIron = 0.098 * p->fMetals;
                    }
#endif
#endif
                }
            else mdlassert(pkd->mdl,0); /* unrecognized particle type */
            
            /* tipsy particle read done */
        
#if defined(SIMPLESF) || defined(STARFORM)
            if (fptCoolAgain!=NULL) {
                fread(&fTmp,sizeof(float),1,fptCoolAgain);
#if defined(STARFORM)
                if (pkdIsGasByOrder(pkd,p)) p->fTimeCoolIsOffUntil = fTmp;
#else
                if (pkdIsGasByOrder(pkd,p)) p->fTimeForm = fTmp;
#endif
                }
            if (fpmStar!=NULL) {
                fread(&fTmp,sizeof(float),1,fpmStar);
                if (pkdIsGasByOrder(pkd,p)) p->fMassForm = fTmp;
                }
#endif
#ifdef INFLOWOUTFLOW
            if (p->r[0] < pkd->dxInflow) { TYPESet(p,TYPE_INFLOW); assert(p->v[0] > 0); }
            if (p->r[0] > pkd->dxOutflow) { TYPESet(p,TYPE_OUTFLOW); assert(p->v[0] > 0); }
#endif
            TYPESet(p,iSetMask); /* needed to get max order info */
			/* Read iOrder last so Types not messed up */
			switch (iReadIOrder) {
			case 0:
#if !(NIORDERGASBUFFER)
				p->iOrder = nStart + i; /* This should be redundant */
#endif
				break;
			case 1:
				fread(&IntTmp,sizeof(IntTmp),1,fpiord);
				p->iOrder = IntTmp;
				break;
			case 2:
				fread(&LongTmp,sizeof(LongTmp),1,fpiord);
				p->iOrder = LongTmp;
				break;
			case 3:
				fread(&p->iOrder,sizeof(p->iOrder),1,fpiord);
				break;
				}
		    }
        } /* native read done */

	if (fpiord!=NULL) fclose(fpiord);
	fclose(fp);
    }


void pkdCalcBound(PKD pkd,BND *pbnd,BND *pbndActive,BND *pbndTreeActive, BND *pbndBall,BNDDT *pbndDt)
{
        /* Faster by assuming active order */

	int i,j;
	FLOAT fBall,r;

	/*
	 ** Initialize the bounds to 0 at the beginning
	 */
	for (j=0;j<3;++j) {
		pbnd->fMin[j] = FLOAT_MAXVAL;
		pbnd->fMax[j] = -FLOAT_MAXVAL;
		pbndActive->fMin[j] = FLOAT_MAXVAL;
		pbndActive->fMax[j] = -FLOAT_MAXVAL;
		pbndBall->fMin[j] = FLOAT_MAXVAL;
		pbndBall->fMax[j] = -FLOAT_MAXVAL;
		}
	/*
	 ** Calculate Active Bounds assume TreeActive is the same as Active (Gnah!)
	 */
	for (i=0;i<pkd->nActive;++i) {
	        fBall = pkd->pStore[i].fBallMax;
	        for (j=0;j<3;++j) {
		        r = pkd->pStore[i].r[j];
			if (r < pbnd->fMin[j]) pbnd->fMin[j] = r;
			if (r > pbnd->fMax[j]) pbnd->fMax[j] = r;
	                if (r < pbndActive->fMin[j]) pbndActive->fMin[j] = r;
			if (r > pbndActive->fMax[j]) pbndActive->fMax[j] = r;
			if (r-fBall < pbndBall->fMin[j]) pbndBall->fMin[j] = r-fBall;
			if (r+fBall > pbndBall->fMax[j]) pbndBall->fMax[j] = r+fBall;
		        }
		}
	/*
	 ** Calculate Local Bounds.
	 */
	for (i=pkd->nActive;i<pkd->nLocal;++i) {
		for (j=0;j<3;++j) {
		        r = pkd->pStore[i].r[j];
			if (r < pbnd->fMin[j]) pbnd->fMin[j] = r;
			if (r > pbnd->fMax[j]) pbnd->fMax[j] = r;
			}
		}

	for (j=0;j<3;++j) {
		pbndTreeActive->fMin[j] = pbndActive->fMin[j];
		pbndTreeActive->fMax[j] = pbndActive->fMax[j];
		}

#ifdef LONGRANGESTEP
	if (pbndDt != NULL) {
	    DIAGDIST2(pbndDt->drMax2,pbnd->fMin,pbnd->fMax);
	    pbndDt->cMax = -FLOAT_MAXVAL;
	    for (j=0;j<3;++j) {
		pbndDt->vMin[j] = FLOAT_MAXVAL;
		pbndDt->vMax[j] = -FLOAT_MAXVAL;
		}
	    for (i=0;i<pkd->nLocal;++i) {
		double v2;
		if (TYPETest(&(pkd->pStore[i]),TYPE_GAS)) {
		    if (pkd->pStore[i].c > pbndDt->cMax)
			pbndDt->cMax = pkd->pStore[i].c;
		    for (j=0;j<3;++j) {
			if (pkd->pStore[i].vPred[j] < pbndDt->vMin[j]) 
			    pbndDt->vMin[j] = pkd->pStore[i].vPred[j];
			if (pkd->pStore[i].vPred[j] > pbndDt->vMax[j])
			    pbndDt->vMax[j] = pkd->pStore[i].vPred[j];
			}
		    }
		}
	    }
#endif
	}


void pkdCalcBound_old(PKD pkd,BND *pbnd,BND *pbndActive,BND *pbndTreeActive, BND *pbndBall)
{
	int i,j;
	FLOAT fBall;

	/*
	 ** Initialize the bounds to 0 at the beginning
	 */
	for (j=0;j<3;++j) {
		pbnd->fMin[j] = FLOAT_MAXVAL;
		pbnd->fMax[j] = -FLOAT_MAXVAL;
		pbndActive->fMin[j] = FLOAT_MAXVAL;
		pbndActive->fMax[j] = -FLOAT_MAXVAL;
		pbndTreeActive->fMin[j] = FLOAT_MAXVAL;
		pbndTreeActive->fMax[j] = -FLOAT_MAXVAL;
		pbndBall->fMin[j] = FLOAT_MAXVAL;
		pbndBall->fMax[j] = -FLOAT_MAXVAL;
		}
	/*
	 ** Calculate Local Bounds.
	 */
	for (i=0;i<pkd->nLocal;++i) {
		for (j=0;j<3;++j) {
			if (pkd->pStore[i].r[j] < pbnd->fMin[j]) 
				pbnd->fMin[j] = pkd->pStore[i].r[j];
			if (pkd->pStore[i].r[j] > pbnd->fMax[j])
				pbnd->fMax[j] = pkd->pStore[i].r[j];
			}
		}
	/*
	 ** Calculate Active Bounds.
	 */
	for (i=0;i<pkd->nLocal;++i) {
		if (TYPEQueryACTIVE(&(pkd->pStore[i]))) {
			for (j=0;j<3;++j) {
				if (pkd->pStore[i].r[j] < pbndActive->fMin[j]) 
					pbndActive->fMin[j] = pkd->pStore[i].r[j];
				if (pkd->pStore[i].r[j] > pbndActive->fMax[j])
					pbndActive->fMax[j] = pkd->pStore[i].r[j];
				}
			}
		}
	/*
	 ** Calculate TreeActive Bounds.
	 */
	for (i=0;i<pkd->nLocal;++i) {
		if (TYPEQueryTREEACTIVE(&(pkd->pStore[i]))) {
			for (j=0;j<3;++j) {
				if (pkd->pStore[i].r[j] < pbndTreeActive->fMin[j]) 
					pbndTreeActive->fMin[j] = pkd->pStore[i].r[j];
				if (pkd->pStore[i].r[j] > pbndTreeActive->fMax[j])
					pbndTreeActive->fMax[j] = pkd->pStore[i].r[j];
				}
			}
		}
	/*
	 ** Calculate fBall Bounds on TreeActive Particles
	 */
	for (i=0;i<pkd->nLocal;++i) {
		if (TYPEQueryTREEACTIVE(&(pkd->pStore[i]))) {
		        fBall = pkd->pStore[i].fBallMax;
			for (j=0;j<3;++j) {
				if (pkd->pStore[i].r[j]-fBall < pbndBall->fMin[j]) 
					pbndBall->fMin[j] = pkd->pStore[i].r[j]-fBall;
				if (pkd->pStore[i].r[j]+fBall > pbndBall->fMax[j])
					pbndBall->fMax[j] = pkd->pStore[i].r[j]+fBall;
				}
			}
		}
    
    }


void pkdGasWeight(PKD pkd)
{
    PARTICLE *p;
    int i;

    for(i=0;i<pkdLocal(pkd);++i) {
		p = &pkd->pStore[i];
		if (TYPETest( p, TYPE_GAS )) p->fWeight=1.0;
		else p->fWeight=0.0;
		}
    }

void pkdRungDDWeight(PKD pkd, int iMaxRung, double dWeight)
{
    PARTICLE *p;
    int i;
    float fRungWeight[50],sum;

    mdlassert(pkd->mdl,iMaxRung<50);
    fRungWeight[0]=1.0;
    sum=1.0;
    for (i=1;i<=iMaxRung;i++) {
                sum*=2.0;
                fRungWeight[i] = dWeight* sum + (1-dWeight);
    }
  
    for(i=0;i<pkdLocal(pkd);++i) {
		p = &pkd->pStore[i];
		p->fWeight *= fRungWeight[p->iRung];
		}
    }

/*
 ** Partition particles between iFrom and iTo into those < fSplit and
 ** those >= to fSplit.  Find number and weight in each partition.
 */
int pkdWeight(PKD pkd,int d,FLOAT fSplit,int iSplitSide,int iFrom,int iTo,
			  int *pnLow,int *pnHigh,FLOAT *pfLow,FLOAT *pfHigh)
{
	int i,iPart;
	FLOAT fLower,fUpper;

	/*
	 ** First partition the memory about fSplit for particles iFrom to iTo.
	 */
	if (iSplitSide) {
		iPart = pkdLowerPart(pkd,d,fSplit,iFrom,iTo);
		*pnLow = pkdLocal(pkd)-iPart;
		*pnHigh = iPart;
		}
	else {
		iPart = pkdUpperPart(pkd,d,fSplit,iFrom,iTo);
		*pnLow = iPart;
		*pnHigh = pkdLocal(pkd)-iPart;
		}
	/*
	 ** Calculate the lower weight and upper weight BETWEEN the particles
	 ** iFrom to iTo!
	 */
	fLower = 0.0;
	for (i=iFrom;i<iPart;++i) {
		fLower += pkd->pStore[i].fWeight;
		}
	fUpper = 0.0;
	for (i=iPart;i<=iTo;++i) {
		fUpper += pkd->pStore[i].fWeight;
		}
	if (iSplitSide) {
		*pfLow = fUpper;
		*pfHigh = fLower;
		}
	else {
		*pfLow = fLower;
		*pfHigh = fUpper;
		}
	return(iPart);
	}


/*
 ** Partition particles between iFrom and iTo into those < fSplit and
 ** those >= to fSplit.  Find number and weight in each partition.
 */
int pkdWeightWrap(PKD pkd,int d,FLOAT fSplit,FLOAT fSplit2, int iSplitSide,int iFrom,int iTo,
			  int *pnLow,int *pnHigh,FLOAT *pfLow,FLOAT *pfHigh)
{
	int i,iPart;
	FLOAT fLower,fUpper;

	/*
	 ** First partition the memory about fSplit for particles iFrom to iTo.
	 */
	if (!iSplitSide) {
		iPart = pkdLowerPartWrap(pkd,d,fSplit,fSplit2,iFrom,iTo);
		*pnLow = iPart;
		*pnHigh = pkdLocal(pkd)-iPart;
		}
	else {
		iPart = pkdUpperPartWrap(pkd,d,fSplit,fSplit2,iFrom,iTo);
		*pnHigh = iPart;
		*pnLow = pkdLocal(pkd)-iPart;
		}
	/*
	 ** Calculate the lower weight and upper weight BETWEEN the particles
	 ** iFrom to iTo!
	 */
	/* Not needed */
	fLower = 0.0;
	for (i=iFrom;i<iPart;++i) {
		fLower += pkd->pStore[i].fWeight;
		}
	fUpper = 0.0;
	for (i=iPart;i<=iTo;++i) {
		fUpper += pkd->pStore[i].fWeight;
		}
	if (!iSplitSide) {
		*pfLow = fLower;
		*pfHigh = fUpper;
		}
	else {
		*pfLow = fUpper;
		*pfHigh = fLower;
		}

	return(iPart);
	}


int pkdOrdWeight(PKD pkd,int iOrdSplit,int iSplitSide,int iFrom,int iTo,
				 int *pnLow,int *pnHigh)
{
	int iPart;
	
	/*
	 ** First partition the memory about fSplit for particles iFrom to iTo.
	 */
	if (iSplitSide) {
		iPart = pkdLowerOrdPart(pkd,iOrdSplit,iFrom,iTo);
		*pnLow = pkdLocal(pkd)-iPart;
		*pnHigh = iPart;
		}
	else {
		iPart = pkdUpperOrdPart(pkd,iOrdSplit,iFrom,iTo);
		*pnLow = iPart;
		*pnHigh = pkdLocal(pkd)-iPart;
		}
	return(iPart);
	}


int pkdLowerPart(PKD pkd,int d,FLOAT fSplit,int i,int j)
{
	PARTICLE pTemp;
	/*
	int i0=i,j0=j,iold;
	*/
#ifdef OLD_KEPLER
	mdlassert(pkd->mdl,d < 5);
	if (d > 2) return pkdLowerQQPart(pkd,d,fSplit,i,j);
#else
	mdlassert(pkd->mdl,d < 3);
#endif

	if (i > j) goto done1;
	i--;
	j++;
	while (++i<j && pkd->pStore[i].r[d] >= fSplit);
        while (i<--j && pkd->pStore[j].r[d] < fSplit);
	if (i>=j) goto done1;
	pTemp = pkd->pStore[i];
	pkd->pStore[i] = pkd->pStore[j];
	pkd->pStore[j] = pTemp;

	while (1) {
	     while (pkd->pStore[++i].r[d] >= fSplit);
	     while (pkd->pStore[--j].r[d] < fSplit);
	     if (i > j) goto done1;
	     pTemp = pkd->pStore[i];
	     pkd->pStore[i] = pkd->pStore[j];
	     pkd->pStore[j] = pTemp;
             }
 done1:
	/*     iold=i;
     i=i0;
     j=j0;
	
	if (i > j) goto done;
    while (1) {
        while (pkd->pStore[i].r[d] >= fSplit)
            if (++i > j) goto done;
        while (pkd->pStore[j].r[d] < fSplit)
            if (i > --j) goto done;
		pTemp = pkd->pStore[i];
		pkd->pStore[i] = pkd->pStore[j];
		pkd->pStore[j] = pTemp;
        }
	
	
 done:
    mdlassert(pkd->mdl,i==iold);
	*/
    return(i);
	}


int pkdUpperPart(PKD pkd,int d,FLOAT fSplit,int i,int j)
{
	PARTICLE pTemp;
	/*	
	int i0=i,j0=j,iold;
	*/

#ifdef OLD_KEPLER
	mdlassert(pkd->mdl,d < 5);
	if (d > 2) return pkdUpperQQPart(pkd,d,fSplit,i,j);
#else
	mdlassert(pkd->mdl,d < 3);
#endif

	if (i > j) goto done1;
	i--;
	j++;
	while (++i<j && pkd->pStore[i].r[d] < fSplit);
	while (i<--j && pkd->pStore[j].r[d] >= fSplit);
	if (i>=j) goto done1;
	pTemp = pkd->pStore[i];
	pkd->pStore[i] = pkd->pStore[j];
	pkd->pStore[j] = pTemp;

	while (1) {
	     while (pkd->pStore[++i].r[d] < fSplit);
	     while (pkd->pStore[--j].r[d] >= fSplit);
	     if (i > j) goto done1;
	     pTemp = pkd->pStore[i];
	     pkd->pStore[i] = pkd->pStore[j];
	     pkd->pStore[j] = pTemp;
             }
 done1:
	/*
     iold=i;
     i=i0;
     j=j0;
	
	if (i > j) goto done;
    while (1) {
        while (pkd->pStore[i].r[d] < fSplit)
            if (++i > j) goto done;
        while (pkd->pStore[j].r[d] >= fSplit)
            if (i > --j) goto done;

		pTemp = pkd->pStore[i];
		pkd->pStore[i] = pkd->pStore[j];
		pkd->pStore[j] = pTemp;
        }
 done:
    mdlassert(pkd->mdl,i==iold);
	*/
    return(i);
	}


int pkdLowerPartWrap(PKD pkd,int d,FLOAT fSplit1,FLOAT fSplit2,int i,int j)
{
	PARTICLE pTemp;

#ifdef OLD_KEPLER
	mdlassert(pkd->mdl,d < 5);
	if (d > 2) return pkdLowerQQPart(pkd,d,fSplit1,i,j);
#else
	mdlassert(pkd->mdl,d < 3);
#endif

	if (fSplit1 > fSplit2) {
	     if (i > j) goto done1;
	     i--;
	     j++;
	     while (++i<j && (pkd->pStore[i].r[d] < fSplit2 || pkd->pStore[i].r[d] >= fSplit1));
	     while (i<--j && (pkd->pStore[j].r[d] >= fSplit2 && pkd->pStore[j].r[d] < fSplit1));
	     if (i>=j) goto done1;
	     pTemp = pkd->pStore[i];
	     pkd->pStore[i] = pkd->pStore[j];
	     pkd->pStore[j] = pTemp;
	     
	     while (1) {
	       while (pkd->pStore[++i].r[d] < fSplit2 || pkd->pStore[i].r[d] >= fSplit1);
	       while (pkd->pStore[--j].r[d] >= fSplit2 && pkd->pStore[j].r[d] < fSplit1);
	       if (i > j) goto done1;
	       pTemp = pkd->pStore[i];
	       pkd->pStore[i] = pkd->pStore[j];
	       pkd->pStore[j] = pTemp;
               }
	     }
	else {
	     if (i > j) goto done1;
	     i--;
	     j++;
	     while (++i<j && (pkd->pStore[i].r[d] < fSplit2 && pkd->pStore[i].r[d] >= fSplit1));
	     while (i<--j && (pkd->pStore[j].r[d] >= fSplit2 || pkd->pStore[j].r[d] < fSplit1));
	     if (i>=j) goto done1;
	     pTemp = pkd->pStore[i];
	     pkd->pStore[i] = pkd->pStore[j];
	     pkd->pStore[j] = pTemp;
	     
	     while (1) {
	       while (pkd->pStore[++i].r[d] < fSplit2 && pkd->pStore[i].r[d] >= fSplit1);
	       while (pkd->pStore[--j].r[d] >= fSplit2 || pkd->pStore[j].r[d] < fSplit1);
	       if (i > j) goto done1;
	       pTemp = pkd->pStore[i];
	       pkd->pStore[i] = pkd->pStore[j];
	       pkd->pStore[j] = pTemp;
               }
	     }
	

 done1:
	return(i);
        }

int pkdUpperPartWrap(PKD pkd,int d,FLOAT fSplit1,FLOAT fSplit2,int i,int j)
{
	PARTICLE pTemp;

#ifdef OLD_KEPLER
	mdlassert(pkd->mdl,d < 5);
	if (d > 2) return pkdUpperQQPart(pkd,d,fSplit1,i,j);
#else
	mdlassert(pkd->mdl,d < 3);
#endif

	if (fSplit1 > fSplit2) {
	     if (i > j) goto done1;
	     i--;
	     j++;
	     while (++i<j && (pkd->pStore[i].r[d] >= fSplit2 && pkd->pStore[i].r[d] < fSplit1));
	     while (i<--j && (pkd->pStore[j].r[d] < fSplit2 || pkd->pStore[j].r[d] >= fSplit1));
	     if (i>=j) goto done1;
	     pTemp = pkd->pStore[i];
	     pkd->pStore[i] = pkd->pStore[j];
	     pkd->pStore[j] = pTemp;
	     
	     while (1) {
	       while (pkd->pStore[++i].r[d] >= fSplit2 && pkd->pStore[i].r[d] < fSplit1);
	       while (pkd->pStore[--j].r[d] < fSplit2 || pkd->pStore[j].r[d] >= fSplit1);
	       if (i > j) goto done1;
	       pTemp = pkd->pStore[i];
	       pkd->pStore[i] = pkd->pStore[j];
	       pkd->pStore[j] = pTemp;
               }
	     }
	else {
	     if (i > j) goto done1;
	     i--;
	     j++;
	     while (++i<j && (pkd->pStore[i].r[d] >= fSplit2 || pkd->pStore[i].r[d] < fSplit1));
	     while (i<--j && (pkd->pStore[j].r[d] < fSplit2 && pkd->pStore[j].r[d] >= fSplit1));
	     if (i>=j) goto done1;
	     pTemp = pkd->pStore[i];
	     pkd->pStore[i] = pkd->pStore[j];
	     pkd->pStore[j] = pTemp;
	     
	     while (1) {
	       while (pkd->pStore[++i].r[d] >= fSplit2 || pkd->pStore[i].r[d] < fSplit1);
	       while (pkd->pStore[--j].r[d] < fSplit2 && pkd->pStore[j].r[d] >= fSplit1);
	       if (i > j) goto done1;
	       pTemp = pkd->pStore[i];
	       pkd->pStore[i] = pkd->pStore[j];
	       pkd->pStore[j] = pTemp;
               }
	     }
	

 done1:
	return(i);
        }



int pkdLowerOrdPart(PKD pkd,int nOrdSplit,int i,int j)
{
	PARTICLE pTemp;

	if (i > j) goto done;
    while (1) {
        while (pkd->pStore[i].iOrder >= nOrdSplit)
            if (++i > j) goto done;
        while (pkd->pStore[j].iOrder < nOrdSplit)
            if (i > --j) goto done;
		pTemp = pkd->pStore[i];
		pkd->pStore[i] = pkd->pStore[j];
		pkd->pStore[j] = pTemp;
        }
 done:
    return(i);
	}


int pkdUpperOrdPart(PKD pkd,int nOrdSplit,int i,int j)
{
	PARTICLE pTemp;

	if (i > j) goto done;
    while (1) {
        while (pkd->pStore[i].iOrder < nOrdSplit)
            if (++i > j) goto done;
        while (pkd->pStore[j].iOrder >= nOrdSplit)
            if (i > --j) goto done;
		pTemp = pkd->pStore[i];
		pkd->pStore[i] = pkd->pStore[j];
		pkd->pStore[j] = pTemp;
        }
 done:
    return(i);
	}


int pkdActiveTypeOrder(PKD pkd, unsigned int iTestMask)
{
	PARTICLE pTemp;
	int i=0;
	int j=pkdLocal(pkd)-1;
	/*
	int i0=i,j0=j,iold;
	*/

	if (i > j) goto done1;
	i--;
	j++;
	while (++i<j && TYPETest(&(pkd->pStore[i]), iTestMask ));
        while (i<--j && !TYPETest(&(pkd->pStore[j]), iTestMask ));
	if (i>=j) goto done1;
	pTemp = pkd->pStore[i];
	pkd->pStore[i] = pkd->pStore[j];
	pkd->pStore[j] = pTemp;

	while (1) {
	     while (TYPETest(&(pkd->pStore[++i]), iTestMask ));
	     while (!TYPETest(&(pkd->pStore[--j]), iTestMask ));
	     if (i > j) goto done1;
	     pTemp = pkd->pStore[i];
	     pkd->pStore[i] = pkd->pStore[j];
	     pkd->pStore[j] = pTemp;
             }
 done1:
	/*	iold=i;
	i=i0;
	j=j0;
	
        if (i > j) goto done;
    while (1) {
        while (TYPETest(&(pkd->pStore[i]), iTestMask ))
            if (++i > j) goto done;
        while (!TYPETest(&(pkd->pStore[j]), iTestMask ))
            if (i > --j) goto done;
                pTemp = pkd->pStore[i];
                pkd->pStore[i] = pkd->pStore[j];
                pkd->pStore[j] = pTemp;
        }
	
 done:
	assert (i==iold);
	*/
	if ( iTestMask & TYPE_ACTIVE )       pkd->nActive = i;
	if ( iTestMask & TYPE_TREEACTIVE )   pkd->nTreeActive = i;
	if ( iTestMask & TYPE_SMOOTHACTIVE ) pkd->nSmoothActive = i;

	return (i);
	}


int pkdActiveOrder(PKD pkd)
{
	PARTICLE pTemp;
	int i=0;
	int j=pkdLocal(pkd)-1;
	/*
	int i0=i,j0=j,iold;
	*/
	
	if (i > j) goto done1;
	i--;
	j++;
	while (++i<j && TYPEQueryACTIVE(&(pkd->pStore[i])));
        while (i<--j && !TYPEQueryACTIVE(&(pkd->pStore[j])));
	if (i>=j) goto done1;
	pTemp = pkd->pStore[i];
	pkd->pStore[i] = pkd->pStore[j];
	pkd->pStore[j] = pTemp;

	while (1) {
	     while (TYPEQueryACTIVE(&(pkd->pStore[++i])));
	     while (!TYPEQueryACTIVE(&(pkd->pStore[--j])));
	     if (i > j) goto done1;
	     pTemp = pkd->pStore[i];
	     pkd->pStore[i] = pkd->pStore[j];
	     pkd->pStore[j] = pTemp;
             }
	
 done1:
	/*
	iold=i;
	i=i0;
	j=j0;

        if (i > j) goto done;
    while (1) {
        while (TYPEQueryACTIVE(&(pkd->pStore[i])))
            if (++i > j) goto done;
        while (!TYPEQueryACTIVE(&(pkd->pStore[j])))
            if (i > --j) goto done;
                pTemp = pkd->pStore[i];
                pkd->pStore[i] = pkd->pStore[j];
                pkd->pStore[j] = pTemp;
        }
	
 done:
        assert (i==iold);
	*/
	return (pkd->nActive = i);
	}


int pkdColRejects_Active_Inactive(PKD pkd,int d,FLOAT fSplit,FLOAT fSplitInactive,
				  int iSplitSide)
{
	PARTICLE pTemp;
	int nSplit,nSplitInactive,iRejects,i,j;

	mdlassert(pkd->mdl,pkd->nRejects == 0);
	if (iSplitSide) {
		nSplit = pkdLowerPart(pkd,d,fSplit,0,pkdActive(pkd)-1);
		}
	else {
		nSplit = pkdUpperPart(pkd,d,fSplit,0,pkdActive(pkd)-1);
		}
	if (iSplitSide) {
		nSplitInactive = pkdLowerPart(pkd,d,fSplitInactive,
									  pkdActive(pkd),pkdLocal(pkd)-1);
		}
	else {
		nSplitInactive = pkdUpperPart(pkd,d,fSplitInactive,
									  pkdActive(pkd),pkdLocal(pkd)-1);
		}
	/*
	for(i = 0; i < nSplit; ++i)
	    mdlassert(pkd->mdl,TYPEQueryACTIVE(&(pkd->pStore[i])));
	for(i = pkdActive(pkd); i < nSplitInactive; ++i)
	    mdlassert(pkd->mdl,!TYPEQueryACTIVE(&(pkd->pStore[i])));
	*/

	nSplitInactive -= pkdActive(pkd);
	/*
	 ** Now do some fancy rearrangement.
	 */
	i = nSplit;
	j = nSplit+nSplitInactive;
	while (j < pkdActive(pkd) + nSplitInactive) {
		pTemp = pkd->pStore[i];
		pkd->pStore[i] = pkd->pStore[j];
		pkd->pStore[j] = pTemp;
		++i; ++j;
		}
	pkd->nRejects = pkdLocal(pkd) - nSplit - nSplitInactive;
	iRejects = pkdFreeStore(pkd) - pkd->nRejects;
	pkd->nActive = nSplit;
	pkd->nLocal = nSplit + nSplitInactive;
	/*
	 ** Move rejects to High memory.
	 */
	for (i=pkd->nRejects-1;i>=0;--i)
		pkd->pStore[iRejects+i] = pkd->pStore[pkd->nLocal+i];
	return(pkd->nRejects);
	}


int pkdColRejects(PKD pkd,int d,FLOAT fSplit,FLOAT fSplitInactive,
				  int iSplitSide)
{
	int nSplit,iRejects,i;

	mdlassert(pkd->mdl,pkd->nRejects == 0);
	if (!iSplitSide) {
		nSplit = pkdLowerPartWrap(pkd,d,fSplitInactive,fSplit,0,pkdLocal(pkd)-1);
		}
	else {
		nSplit = pkdUpperPartWrap(pkd,d,fSplitInactive,fSplit,0,pkdLocal(pkd)-1);
		}

	pkd->nRejects = pkdLocal(pkd) - nSplit;
	iRejects = pkdFreeStore(pkd) - pkd->nRejects;
	pkd->nLocal = nSplit;
	/*
	 ** Move rejects to High memory.
	 */
	for (i=pkd->nRejects-1;i>=0;--i)
		pkd->pStore[iRejects+i] = pkd->pStore[pkd->nLocal+i];
	return(pkd->nRejects);
	}


int pkdSwapRejects(PKD pkd,int idSwap)
{
	size_t nBuf;
	size_t nOutBytes,nSndBytes,nRcvBytes;

	if (idSwap != -1) {
		nBuf = (pkdSwapSpace(pkd))*sizeof(PARTICLE);
		nOutBytes = pkd->nRejects*sizeof(PARTICLE);
		mdlassert(pkd->mdl,pkdLocal(pkd) + pkd->nRejects <= pkdFreeStore(pkd));
		mdlSwap(pkd->mdl,idSwap,nBuf,&pkd->pStore[pkdLocal(pkd)],
				nOutBytes,&nSndBytes,&nRcvBytes);
		pkd->nLocal += nRcvBytes/sizeof(PARTICLE);
		pkd->nRejects -= nSndBytes/sizeof(PARTICLE);
		}
	return(pkd->nRejects);
	}

void pkdSwapAll(PKD pkd, int idSwap)
{
    size_t nBuf;
    size_t nOutBytes,nSndBytes,nRcvBytes;
    int i;
    int iBuf;
    int bGood;
    
    /*
     ** Move particles to High memory.
     */
    iBuf = pkdSwapSpace(pkd);
    for (i=pkdLocal(pkd)-1;i>=0;--i)
	pkd->pStore[iBuf+i] = pkd->pStore[i];

    nBuf = pkdFreeStore(pkd)*sizeof(PARTICLE);
    nOutBytes = pkdLocal(pkd)*sizeof(PARTICLE);
    bGood = mdlSwap(pkd->mdl,idSwap,nBuf,&pkd->pStore[0], nOutBytes,
			&nSndBytes, &nRcvBytes);
    mdlassert(pkd->mdl,bGood);
    mdlassert(pkd->mdl,nSndBytes/sizeof(PARTICLE) == pkdLocal(pkd));
    pkd->nLocal = nRcvBytes/sizeof(PARTICLE);
    }

int pkdSwapSpace(PKD pkd)
{
	return(pkdFreeStore(pkd) - pkdLocal(pkd));
	}


int pkdFreeStore(PKD pkd)
{
	return(pkd->nStore);
	}

int pkdActive(PKD pkd)
{
	return(pkd->nActive);
	}

int pkdTreeActive(PKD pkd)
{
	return(pkd->nTreeActive);
	}

int pkdSmoothActive(PKD pkd)
{
	return(pkd->nSmoothActive);
	}

int pkdInactive(PKD pkd)
{
	return(pkd->nLocal - pkd->nActive);
	}

int pkdTreeInactive(PKD pkd)
{
	return(pkd->nLocal - pkd->nTreeActive);
	}

int pkdLocal(PKD pkd)
{
	return(pkd->nLocal);
	}

void pkdTotals(PKD pkd, int *nDark, int *nGas, int *nStar)
{
    int i,nd,ng,ns;
    nd=ng=ns=0;
    for (i=0;i<pkd->nLocal;++i) {
        if (pkdIsDark(pkd,&pkd->pStore[i])) nd++;
        if (pkdIsGas(pkd,&pkd->pStore[i])) ng++;
        if (pkdIsStar(pkd,&pkd->pStore[i])) ns++;
        }
    *nDark=nd;
    *nGas=ng;
    *nStar=ns;
    }

int pkdNodes(PKD pkd)
{
	return(pkd->nNodes);
	}


void pkdDomainColor(PKD pkd)
{
#ifdef COLORCODE
	int i;
	
	for (i=0;i<pkd->nLocal;++i) {
		pkd->pStore[i].fColor = (FLOAT)pkd->idSelf;
		}
#endif
	}


int pkdColOrdRejects(PKD pkd,int nOrdSplit,int iSplitSide)
{
	int nSplit,iRejects,i;

	if (iSplitSide) nSplit = pkdLowerOrdPart(pkd,nOrdSplit,0,pkdLocal(pkd)-1);
	else nSplit = pkdUpperOrdPart(pkd,nOrdSplit,0,pkdLocal(pkd)-1);
	pkd->nRejects = pkdLocal(pkd) - nSplit;
	iRejects = pkdFreeStore(pkd) - pkd->nRejects;
	pkd->nLocal = nSplit;
	/*
	 ** Move rejects to High memory.
	 */
	for (i=pkd->nRejects-1;i>=0;--i)
		pkd->pStore[iRejects+i] = pkd->pStore[pkd->nLocal+i];
	return(pkd->nRejects);
	}


int cmpParticles(const void *pva,const void *pvb)
{
	PARTICLE *pa = (PARTICLE *)pva;
	PARTICLE *pb = (PARTICLE *)pvb;

	return(pa->iOrder - pb->iOrder);
	}


void pkdLocalOrder(PKD pkd)
{
	qsort(pkd->pStore,pkdLocal(pkd),sizeof(PARTICLE),cmpParticles);
	}


void pkdWriteTipsy(PKD pkd,char *pszFileName,int nStart,
		   int bStandard,double dvFac,double duTFac,int iGasModel)
{
  PARTICLE *p;
  FILE *fp;
  int i,j;
  struct dark_particle dp;
  struct gas_particle gp;
  struct star_particle sp;
  int nout;
  float fTmp;
  
  /*
  ** Seek past the header and up to nStart.
  */
  fp = fopen(pszFileName,"r+");
  mdlassert(pkd->mdl,fp != NULL);
  pkdSeek(pkd,fp,nStart,bStandard);
  if (bStandard) {
    FLOAT vTemp;
    XDR xdrs;
    /* 
    ** Write Stuff!
    */
    xdrstdio_create(&xdrs,fp,XDR_ENCODE);
    for (i=0;i<pkdLocal(pkd);++i) {
      p = &pkd->pStore[i];
      if (pkdIsDark(pkd,p)) {
	fTmp = p->fMass;
	xdr_float(&xdrs,&fTmp);
	for (j=0;j<3;++j) {
	  fTmp = p->r[j];
	  xdr_float(&xdrs,&fTmp);
	}
	for (j=0;j<3;++j) {
	  vTemp = dvFac*p->v[j];			
	  fTmp = vTemp;
	  xdr_float(&xdrs,&fTmp);
	}
#ifdef CHANGESOFT
	fTmp = p->fSoft0;
#else
	fTmp = p->fSoft;
#endif
	xdr_float(&xdrs,&fTmp);
	fTmp = p->fPot;
	xdr_float(&xdrs,&fTmp);
      }
      else if (pkdIsGas(pkd,p)) {
#ifdef SINKING
	if (TYPETest(p,TYPE_SINK)) 
	    fTmp = p->fTrueMass;
	else
#endif
	    fTmp = p->fMass;
	xdr_float(&xdrs,&fTmp);
	for (j=0;j<3;++j) {
	  fTmp = p->r[j];
	  xdr_float(&xdrs,&fTmp);
	}
#ifdef SINKINGOUTVPRED
	if (TYPETest(p,TYPE_SINKING)) {
#ifdef SINKDBG
	  if (p->iOrder == 55) printf("SINKINGKICKOUT %d, %g %g  %g %g  %g %g\n",p->iOrder,p->vPred[0],p->v[0],p->vPred[1],p->v[1],p->vPred[2],p->v[2]);
#endif
	  for (j=0;j<3;++j) {
	    vTemp = dvFac*p->vPred[j];			
	    fTmp = vTemp;
	    xdr_float(&xdrs,&fTmp);
	  }
	}
	else 
#endif
	  {
	    for (j=0;j<3;++j) {
	      vTemp = dvFac*p->v[j];			
	      fTmp = vTemp;
	      xdr_float(&xdrs,&fTmp);
	    }
	  }
	fTmp = p->fDensity;
	xdr_float(&xdrs,&fTmp);
#ifdef GASOLINE
	/*
	** Convert thermal energy to tempertature.
	*/
	switch (iGasModel) {
	case GASMODEL_COOLING:
#ifndef NOCOOLING
	    vTemp = CoolCodeEnergyToTemperature( pkd->Cool, &p->CoolParticle, p->u, p->fMetals );
#else
	    mdlassert(pkd->mdl,0);
#endif
	    break;
	default:
	    vTemp = duTFac*p->u;
	    }
	fTmp = vTemp;
	xdr_float(&xdrs,&fTmp);
	/* fTmp = sqrt(0.25*p->fBall2);  Write softening in tipsy outputs */
#ifdef CHANGESOFT
	fTmp = p->fSoft0;
#else
	fTmp = p->fSoft;
#endif
#ifdef DRHODT 
	/* Horrible hack -- overwrite soft output */
#ifdef DRHODTDIVOUT
	fTmp = p->fDivv_PdV;
#else
	fTmp = p->fDensity_PdV;
#endif
#endif
	xdr_float(&xdrs,&fTmp);
#ifdef SINKING
	if (TYPETest( p, TYPE_SINKING)) {
	    xdr_int(&xdrs,&p->iSinkingOnto); /* output sinkingonto integer in place of metals */
	    }
	else if (TYPETest( p, TYPE_SINK)) {
	    fTmp = p->fMass;
	    xdr_float(&xdrs,&fTmp);
	    }
	else {
	    fTmp = p->fMetals;
	    xdr_float(&xdrs,&fTmp);
	    }
#else
	fTmp = p->fMetals;
#ifdef DRHODT 
	/* Horrible hack -- overwrite metals output */
#ifdef DRHODTDIVOUT
	fTmp = p->fDivv_t;
#else
	fTmp = p->fDensity_t;
#endif
#endif
	xdr_float(&xdrs,&fTmp);
#endif
#else /* not gasoline */
	fTmp = 0.0;
	xdr_float(&xdrs,&fTmp);
#ifdef CHANGESOFT
	fTmp = p->fSoft0;
#else
	fTmp = p->fSoft;
#endif
	xdr_float(&xdrs,&fTmp);
	fTmp = 0.0;
	xdr_float(&xdrs,&fTmp);
#endif
	fTmp = p->fPot;
#ifdef DRHODT 
	/* Horrible hack -- overwrite pot output */
#ifdef DRHODTDIVOUT
	fTmp = p->fDivv_PdVcorr;
#else
	fTmp = p->fDensity_PdVcorr;
#endif
#endif
	xdr_float(&xdrs,&fTmp);
      }
      else if (pkdIsStar(pkd,p)) {
#ifdef SINKING
	if (TYPETest(p,TYPE_SINK)) 
	    fTmp = p->fTrueMass;
	else
#endif
	    fTmp = p->fMass;
	xdr_float(&xdrs,&fTmp);
	for (j=0;j<3;++j) {
	  fTmp = p->r[j];
	  xdr_float(&xdrs,&fTmp);
	}
	for (j=0;j<3;++j) {
	  vTemp = dvFac*p->v[j];			
	  fTmp = vTemp;
	  xdr_float(&xdrs,&fTmp);
	}
#ifdef GASOLINE
#ifdef SINKING
	if (TYPETest( p, TYPE_SINK)) 
	    fTmp = p->fMass;
	else
#endif
	    fTmp = p->fMetals;
	xdr_float(&xdrs,&fTmp);
	fTmp = p->fTimeForm;
	xdr_float(&xdrs,&fTmp);
#else
	fTmp = 0.0;
	xdr_float(&xdrs,&fTmp);
	xdr_float(&xdrs,&fTmp);			
#endif
#ifdef CHANGESOFT
	fTmp = p->fSoft0;
#else
	fTmp = p->fSoft;
#endif
	xdr_float(&xdrs,&fTmp);
	fTmp = p->fPot;
	xdr_float(&xdrs,&fTmp);
      }
      else mdlassert(pkd->mdl,0);
    }
    xdr_destroy(&xdrs);
  }
  else {
    /* 
    ** Write Stuff!
    */
    for (i=0;i<pkdLocal(pkd);++i) {
      p = &pkd->pStore[i];
      if (pkdIsDark(pkd,p)) {
	for (j=0;j<3;++j) {
	  dp.pos[j] = p->r[j];
	  dp.vel[j] = dvFac*p->v[j];
	}
	dp.mass = p->fMass;
#ifdef CHANGESOFT
	dp.eps = p->fSoft0;
#else
	dp.eps = p->fSoft;
#endif
	dp.phi = p->fPot;
	nout = fwrite(&dp,sizeof(struct dark_particle),1,fp);
	mdlassert(pkd->mdl,nout == 1);
      }
      else if (pkdIsGas(pkd,p)) {
	for (j=0;j<3;++j) {
	  gp.pos[j] = p->r[j];
	  gp.vel[j] = dvFac*p->v[j];
	}
#ifdef SINKING
	if (TYPETest(p,TYPE_SINK)) 
	    gp.mass = p->fTrueMass;
	else
#endif
	    gp.mass = p->fMass;
#ifdef CHANGESOFT
	gp.hsmooth = p->fSoft0;
#else
	gp.hsmooth = p->fSoft;
#endif
	gp.phi = p->fPot;
	gp.rho = p->fDensity;
#ifdef GASOLINE
	gp.temp = duTFac*p->u;
#ifdef SINKING
	if (TYPETest( p, TYPE_SINK)) 
	    gp.metals = p->fMass;
	else 
	    *((int *) (&gp.metals)) = p->iSinkingOnto;
#else
	gp.metals = p->fMetals;
#endif
#else
	gp.temp = 0.0;
	gp.metals = 0.0;
#endif
	nout = fwrite(&gp,sizeof(struct gas_particle),1,fp);
	mdlassert(pkd->mdl,nout == 1);
      }
      else if (pkdIsStar(pkd,p)) {
	for (j=0;j<3;++j) {
	  sp.pos[j] = p->r[j];
	  sp.vel[j] = dvFac*p->v[j];
	}
#ifdef SINKING
	if (TYPETest(p,TYPE_SINK)) 
	    sp.mass = p->fTrueMass;
	else
#endif
	    sp.mass = p->fMass;
#ifdef CHANGESOFT
	sp.eps = p->fSoft0;
#else
	sp.eps = p->fSoft;
#endif
	sp.phi = p->fPot;
#ifdef GASOLINE
#ifdef SINKING
	if (TYPETest(p,TYPE_SINK)) 
	    sp.metals = p->fMass;
	else
#endif
	    sp.metals = p->fMetals;

	sp.tform = p->fTimeForm;
#else
	sp.metals = 0.0;
	sp.tform = 0.0;
#endif
	nout = fwrite(&sp,sizeof(struct star_particle),1,fp);
	mdlassert(pkd->mdl,nout == 1);
      }
      else mdlassert(pkd->mdl,0);
    }
  }
  nout = fclose(fp);
  mdlassert(pkd->mdl,nout == 0);
}

void pkdTreeZip(PKD pkd,char *pszFileName, double *dmin, double *dmax) 
    {
    char file[500];
    FILE *fp;
    TZX *tz;
    int nBits_nPerBucket = 4;
    int nBits_Position = 6;
    int nBitsMaxPrecision_PerDirection = TZNBITS_IN_KEY/3;
    int nBitsMinPrecision_PerDirection = 0;

    PARTICLE *p;
    int i;

    printf("id %d: TreeZip: Building tree %d\n",pkd->idSelf,pkd->nLocal);
    if (!pkd->nLocal) return;

    tz = tzInit( &dmin[0], &dmax[0], nBits_nPerBucket, nBits_Position, nBitsMaxPrecision_PerDirection, nBitsMinPrecision_PerDirection );
    for (i=0; i<pkd->nLocal; i++) {
	p = &pkd->pStore[i];
        tzAddPos( tz, &(p->r[0]), p->iOrder );
	}

    sprintf(file,"%s/p.%d",pszFileName,pkd->idSelf);
    fp = fopen(file,"a");
    printf("id %d: TreeZip: Opening %s\n",pkd->idSelf,file);
    assert(fp!=NULL);
    tzOutputFile( tz, fp );
    tzWriteHeader( tz );
    tzWriteTreeZip( tz ) ;
    fclose(fp);

    printf("id %d: Leftovers:  nParticle %i, nBucket %i, nNode %i\n",pkd->idSelf,tz->nParticle,tz->nBucket,tz->nNode);
    printf("id %d: nBits Particle %i, Total %i, Written %i (%i)\n",pkd->idSelf,tz->nParticlebits,tz->nTotalbits,tz->nWritebits,tz->nWritebits/8);

    tzFinalize( tz );

    }

void pkdOutputBlackHoles(PKD pkd,char *pszFileName, double dvFac)
{
#ifdef STARFORM
        PARTICLE *p;
        FILE *fp;
        int i;
        int nout;

        fp = fopen(pszFileName,"a");
        mdlassert(pkd->mdl,fp != NULL);
                /*
                 ** Write Stuff!
                 */
        for (i=0; i<pkd->nLocal; i++) {
            p = &pkd->pStore[i];
	    if(p->fTimeForm < 0.0)
		fprintf(fp, "%i %g %g %g %g %g %g %g %g\n", p->iOrder,
			p->fMass, p->r[0], p->r[1], p->r[2], dvFac*p->v[0],
			dvFac*p->v[1], dvFac*p->v[2], p->fPot);
            }
        nout = fclose(fp);
        mdlassert(pkd->mdl,nout == 0);
#endif
        }

void pkdSetSoft(PKD pkd,double dSoft)
{
	PARTICLE *p;
	int i,n;

	p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;++i) {
#ifdef CHANGESOFT
		p[i].fSoft0 = dSoft;
#endif
		p[i].fSoft = dSoft;
		}
	}

#ifdef CHANGESOFT
void pkdPhysicalSoft(PKD pkd,double dSoftMax,double dFac,int bSoftMaxMul)
{
	PARTICLE *p;
	int i,n;

	p = pkd->pStore;
	n = pkdLocal(pkd);
	
	assert(dFac > 0.0);
	if (bSoftMaxMul) {
	        for (i=0;i<n;++i) { 
		  /*		  assert(p[i].fSoft0 > 0.);
				  printf("Compton:  %f\n",p[i].fSoft0);*/
                  
		        assert(p[i].fSoft0 > 0.0);
	                p[i].fSoft = p[i].fSoft0*dFac;
			assert(p[i].fSoft > 0.0);
		        }
	        }
	else {
 	        assert(dSoftMax > 0.0);
	        for (i=0;i<n;++i) {
	                p[i].fSoft = p[i].fSoft0*dFac;
	                if (p[i].fSoft > dSoftMax) p[i].fSoft = dSoftMax;
		        }
	        }
	}

void pkdPreVariableSoft(PKD pkd,int iVariableSoftType)
{
	PARTICLE *p;
	int i,n;

	p = pkd->pStore;
	n = pkdLocal(pkd);

#ifndef DENSSOFT
	for (i=0;i<n;++i) {
	        if (TYPETest(&(p[i]),iVariableSoftType) && TYPEQueryACTIVE(&(p[i]))) {
			  p[i].fSoft = p[i].fBall2;
			  }
	        }
#endif
	}

void pkdPostVariableSoft(PKD pkd,double dSoftMax,int bSoftMaxMul,int iVariableSoftType)
{
	PARTICLE *p;
	int i,n;
	double dTmp;

	p = pkd->pStore;
	n = pkdLocal(pkd);
	
	if (bSoftMaxMul) {
	        for (i=0;i<n;++i) {
	                if (TYPETest(&(p[i]),iVariableSoftType) && TYPEQueryACTIVE(&(p[i]))) {
#ifdef DENSSOFT
	                          p[i].fSoft = pow((p[i].fMass*1.90986/p[i].fDensity),.3333333333);
#else
		                  dTmp = sqrt(p[i].fBall2*.25);
	                          p[i].fBall2 = p[i].fSoft;
	                          p[i].fSoft = (dTmp <= p[i].fSoft0*dSoftMax ? dTmp : p[i].fSoft0*dSoftMax);
#endif
                                  }
		        }
	        }
	else {
	        for (i=0;i<n;++i) {
#ifdef CHECKSOFT			  
			  if (p[i].iOrder == CHECKSOFT) fprintf(stderr,"Particle %i: %g %g %g %i %g\n",p[i].iOrder,p[i].fDensity,sqrt(p[i].fBall2*.25),p[i].fSoft,p[i].iActive,dSoftMax);
#endif
			  if (TYPETest(&(p[i]),iVariableSoftType) && TYPEQueryACTIVE(&(p[i]))) {
#ifdef DENSSOFT
				  p[i].fSoft = pow((p[i].fMass*1.90986/p[i].fDensity),.3333333333);
#else
				  dTmp = sqrt(p[i].fBall2*.25);
				  p[i].fBall2 = p[i].fSoft;
				  p[i].fSoft = (dTmp <= dSoftMax ? dTmp : dSoftMax);
#endif
			          }
#ifdef CHECKSOFT			  
			  if (p[i].iOrder == CHECKSOFT) fprintf(stderr,"Particle %iA: %g %g %g %i %g\n",p[i].iOrder,p[i].fDensity,sqrt(p[i].fBall2*.25),p[i].fSoft,p[i].iActive,dSoftMax);
#endif
		        }
	        }
	}
#endif

void pkdCombine(KDN *p1,KDN *p2,KDN *pOut)
{
	int j;

	/*
	 ** Combine the bounds.
	 */
	for (j=0;j<3;++j) {
		if (p2->bnd.fMin[j] < p1->bnd.fMin[j])
			pOut->bnd.fMin[j] = p2->bnd.fMin[j];
		else
			pOut->bnd.fMin[j] = p1->bnd.fMin[j];
		if (p2->bnd.fMax[j] > p1->bnd.fMax[j])
			pOut->bnd.fMax[j] = p2->bnd.fMax[j];
		else
			pOut->bnd.fMax[j] = p1->bnd.fMax[j];

		if (p2->bndBall.fMin[j] < p1->bndBall.fMin[j])
			pOut->bndBall.fMin[j] = p2->bndBall.fMin[j];
		else
			pOut->bndBall.fMin[j] = p1->bndBall.fMin[j];
		if (p2->bndBall.fMax[j] > p1->bndBall.fMax[j])
			pOut->bndBall.fMax[j] = p2->bndBall.fMax[j];
		else
			pOut->bndBall.fMax[j] = p1->bndBall.fMax[j];

#ifdef LONGRANGESTEP
		if (p2->bndDt.vMin[j] < p1->bndDt.vMin[j])
			pOut->bndDt.vMin[j] = p2->bndDt.vMin[j];
		else
			pOut->bndDt.vMin[j] = p1->bndDt.vMin[j];
		if (p2->bndDt.vMax[j] > p1->bndDt.vMax[j])
			pOut->bndDt.vMax[j] = p2->bndDt.vMax[j];
		else
			pOut->bndDt.vMax[j] = p1->bndDt.vMax[j];
#endif
		}

#ifdef LONGRANGESTEP
	DIAGDIST2(pOut->bndDt.drMax2,pOut->bnd.fMin,pOut->bnd.fMax);
	if (p2->bndDt.cMax > p1->bndDt.cMax)
	    pOut->bndDt.cMax = p2->bndDt.cMax;
	else
	    pOut->bndDt.cMax = p1->bndDt.cMax;
#endif
	/*
	 ** Find the center of mass and mass weighted softening.
	 */
        pOut->fMass = p1->fMass + p2->fMass;
	pOut->fSoft = p1->fMass*p1->fSoft + p2->fMass*p2->fSoft;
	for (j=0;j<3;++j) {
		pOut->r[j] = p1->fMass*p1->r[j] + p2->fMass*p2->r[j];
		}
	if (pOut->fMass > 0) {
		pOut->fSoft /= pOut->fMass;
		for (j=0;j<3;++j) {
			pOut->r[j] /= pOut->fMass;
			}
		}
	}


void pkdCalcCell(PKD pkd,KDN *pkdn,FLOAT *rcm,int iOrder,
				 struct pkdCalcCellStruct *pcc)
{
	int pj;
	int bSinkCell=0;
	double m,dx,dy,dz,d2,d1;
	struct pkdCalcCellStruct cc;
#ifdef  RADIATIVEBOX
	double fLum = 0.0,dxg,dyg,dzg,d2g,mg;
	int j;
#endif

	/*
	 ** Initialize moments.
	 ** Initialize various B numbers.
	 */
	switch (iOrder) {	
	case 4:
		cc.Hxxxx = 0.0;
		cc.Hxyyy = 0.0;
		cc.Hxxxy = 0.0;
		cc.Hyyyy = 0.0;
		cc.Hxxxz = 0.0;
		cc.Hyyyz = 0.0;
		cc.Hxxyy = 0.0;
		cc.Hxxyz = 0.0;
		cc.Hxyyz = 0.0;
		cc.Hxxzz = 0.0;
		cc.Hxyzz = 0.0;
		cc.Hxzzz = 0.0;
		cc.Hyyzz = 0.0;
		cc.Hyzzz = 0.0;
		cc.Hzzzz = 0.0;
		cc.B6 = 0.0;
	case 3:
		cc.Oxxx = 0.0;
		cc.Oxyy = 0.0;
		cc.Oxxy = 0.0;
		cc.Oyyy = 0.0;
		cc.Oxxz = 0.0;
		cc.Oyyz = 0.0;
		cc.Oxyz = 0.0;
		cc.Oxzz = 0.0;
		cc.Oyzz = 0.0;
		cc.Ozzz = 0.0;
		cc.B5 = 0.0;
	default:
		cc.Qxx = 0.0;
		cc.Qyy = 0.0;
		cc.Qzz = 0.0;
		cc.Qxy = 0.0;
		cc.Qxz = 0.0;
		cc.Qyz = 0.0;
		cc.B2 = 0.0;
		cc.B3 = 0.0;
		cc.B4 = 0.0;
		}
	cc.Bmax = 0.0;
#ifdef  RADIATIVEBOX
	cc.fLW = 0.0; /* total LW radiation in box*/
	cc.gmass = 0.0; /* total mass of gas*/
	cc.gmom = 0.0; /* moment of gas*/
	for (j=0;j<3;++j) cc.cLumLW[j] = 0; /* luminosity weighted center of LW radiation */
#endif
	/*
	 ** Calculate moments and B numbers about center-of-mass.
	 */
	if (pkdn == NULL) pkdn = &pkd->kdNodes[pkd->iRoot];
	for (pj=pkdn->pLower;pj<=pkdn->pUpper;++pj) {
		m = pkd->pStore[pj].fMass;
		dx = pkd->pStore[pj].r[0] - rcm[0];
		dy = pkd->pStore[pj].r[1] - rcm[1];
		dz = pkd->pStore[pj].r[2] - rcm[2];
		d2 = dx*dx + dy*dy + dz*dz;
		d1 = sqrt(d2);
		if (d1 > cc.Bmax) cc.Bmax = d1;
#ifdef STARSINK
		if (TYPETest(&(pkd->pStore[pj]),TYPE_STAR)) bSinkCell=1;
#endif
		cc.B2 += m*d2;
		cc.B3 += m*d2*d1;
		cc.B4 += m*d2*d2;
		cc.B5 += m*d2*d2*d1;
		cc.B6 += m*d2*d2*d2;
#ifdef COMPLETE_LOCAL
		d2 = 0.0;
#endif
#ifdef  RADIATIVEBOX
		if (TYPETest(&(pkd->pStore[pj]),TYPE_GAS)) { /*find the moment of the gas -- this will be used to determine average distance of the gas from the average source of radiation for the purposes of combining boxes*/
		  mg = pkd->pStore[pj].fMass;
		  dxg = pkd->pStore[pj].r[0] - rcm[0];
		  dyg = pkd->pStore[pj].r[1] - rcm[1];
		  dzg = pkd->pStore[pj].r[2] - rcm[2];
		  d2g = dxg*dxg + dyg*dyg + dzg*dzg;
		  cc.gmom += mg*d2g;
		  cc.gmass += mg;
		}
		
		if (TYPETest(&(pkd->pStore[pj]),TYPE_STAR)) {
		  fLum = CoolLymanWerner(pkd->Cool, pkd->pStore[pj].fMassForm, pkd->pStore[pj].CoolParticle.dLymanWerner);
		  pkd->pStore[pj].CoolParticle.dLymanWerner = fLum;
		  for (j=0;j<3;++j) cc.cLumLW[j] += fLum*pkd->pStore[pj].r[j]; /*average location of LW source*/
		  cc.fLW += fLum; /*Total LW radiation*/
		}
 /* Sums the LW flux emitted from the star particles in the bucket*/
#endif
		switch (iOrder) {
		case 4:
			/*
			 ** Calculate reduced hexadecapole moment...
			 */
			cc.Hxxxx += m*(dx*dx*dx*dx - 6.0/7.0*d2*(dx*dx - 0.1*d2));
			cc.Hxyyy += m*(dx*dy*dy*dy - 3.0/7.0*d2*dx*dy);
			cc.Hxxxy += m*(dx*dx*dx*dy - 3.0/7.0*d2*dx*dy);
			cc.Hyyyy += m*(dy*dy*dy*dy - 6.0/7.0*d2*(dy*dy - 0.1*d2));
			cc.Hxxxz += m*(dx*dx*dx*dz - 3.0/7.0*d2*dx*dz);
			cc.Hyyyz += m*(dy*dy*dy*dz - 3.0/7.0*d2*dy*dz);
			cc.Hxxyy += m*(dx*dx*dy*dy - 1.0/7.0*d2*(dx*dx + dy*dy - 0.2*d2));
			cc.Hxxyz += m*(dx*dx*dy*dz - 1.0/7.0*d2*dy*dz);
			cc.Hxyyz += m*(dx*dy*dy*dz - 1.0/7.0*d2*dx*dz);
			cc.Hxxzz += m*(dx*dx*dz*dz - 1.0/7.0*d2*(dx*dx + dz*dz - 0.2*d2));
			cc.Hxyzz += m*(dx*dy*dz*dz - 1.0/7.0*d2*dx*dy);
			cc.Hxzzz += m*(dx*dz*dz*dz - 3.0/7.0*d2*dx*dz);
			cc.Hyyzz += m*(dy*dy*dz*dz - 1.0/7.0*d2*(dy*dy + dz*dz - 0.2*d2));
			cc.Hyzzz += m*(dy*dz*dz*dz - 3.0/7.0*d2*dy*dz);
			cc.Hzzzz += m*(dz*dz*dz*dz - 6.0/7.0*d2*(dz*dz - 0.1*d2));
		case 3:
			/*
			 ** Calculate reduced octopole moment...
			 */
			cc.Oxxx += m*(dx*dx*dx - 0.6*d2*dx);
			cc.Oxyy += m*(dx*dy*dy - 0.2*d2*dx);
			cc.Oxxy += m*(dx*dx*dy - 0.2*d2*dy);
			cc.Oyyy += m*(dy*dy*dy - 0.6*d2*dy);
			cc.Oxxz += m*(dx*dx*dz - 0.2*d2*dz);
			cc.Oyyz += m*(dy*dy*dz - 0.2*d2*dz);
			cc.Oxyz += m*dx*dy*dz;
			cc.Oxzz += m*(dx*dz*dz - 0.2*d2*dx);
			cc.Oyzz += m*(dy*dz*dz - 0.2*d2*dy);
			cc.Ozzz += m*(dz*dz*dz - 0.6*d2*dz);;
		default:
			/*
			 ** Calculate quadrupole moment...
			 */
			cc.Qxx += m*dx*dx;
			cc.Qyy += m*dy*dy;
			cc.Qzz += m*dz*dz;
			cc.Qxy += m*dx*dy;
			cc.Qxz += m*dx*dz;
			cc.Qyz += m*dy*dz;
			}
		}
#ifdef STARSINK
	if (bSinkCell) {
	    cc.Bmax*=2;
	    bSinkCell=0;
	    }
#endif
#ifdef  RADIATIVEBOX
	if (cc.gmass > 0) cc.gmom /= cc.gmass; /*Finish determining the moment of the gas*/
	if (cc.fLW > 0) {
	  for (j=0;j<3;++j) {
	    cc.cLumLW[j] /= cc.fLW; /*Finish average center of LW radiation*/
	  }
	}
#endif
	*pcc = cc;
	}


double fcnAbsMono(KDN *pkdn,double r)
{
	double t;

	t = r*(r - pkdn->mom.Bmax);
	t *= t;
	t = (3.0*pkdn->mom.B2 - 2.0*pkdn->mom.B3/r)/t;
	return t;
	}


double fcnAbsQuad(KDN *pkdn,double r)
{
	double t;

	t = r*(r - pkdn->mom.Bmax);
	t *= r*t;
	t = (4.0*pkdn->mom.B3 - 3.0*pkdn->mom.B4/r)/t;
	return t;
	}


double fcnAbsOct(KDN *pkdn,double r)
{
	double t;

	t = r*r*(r - pkdn->mom.Bmax);
	t *= t;
	t = (5.0*pkdn->mom.B4 - 4.0*pkdn->mom.B5/r)/t;
	return t;
	}


double fcnAbsHex(KDN *pkdn,double r)
{
	double t;

	t = r*r*(r - pkdn->mom.Bmax);
	t *= r*t;
	t = (6.0*pkdn->mom.B5 - 5.0*pkdn->mom.B6/r)/t;
	return t;
	}


double dRootBracket(KDN *pkdn,double dErrBnd,double (*fcn)(KDN *,double))
{
	double dLower,dUpper,dMid;
	double dLowerErr,dErrDif,dErrCrit;
	int iter = 0;

	dErrCrit = 1e-6*dErrBnd;
	dLower = (1.0 + 1e-6)*pkdn->mom.Bmax;
	dLowerErr = (*fcn)(pkdn,dLower);
	if (dLowerErr < dErrBnd) {
/*
		printf("iter:%d %g\n",iter,dLower);
*/
		return dLower;
		}
	/*
	 ** Hunt phase for upper bound.
	 */
	dUpper = 2*dLower;
	while (1) {
		dErrDif = dErrBnd - (*fcn)(pkdn,dUpper);
		if (dErrDif > dErrCrit) break;
		dUpper = 2*dUpper;
		}
	/*
	 ** Bracket root.
	 */
	while (1) {
		dMid = 0.5*(dLower + dUpper);
		dErrDif = dErrBnd - (*fcn)(pkdn,dMid);
		if (dErrDif < 0) {
			dLower = dMid;
			}
		else {
			dUpper = dMid;
			if (dErrDif < dErrCrit) break;
			}
		if (++iter > 32) break;
		}
/*
	printf("iter:%d %g\n",iter,dMid);
*/
	return dMid;
	}


double pkdCalcOpen(KDN *pkdn,int iOpenType,double dCrit,int iOrder)
{
	double dOpen=0;

	if (iOpenType == OPEN_ABSPAR) {
		switch (iOrder) {
		case 1:
			dOpen = dRootBracket(pkdn,dCrit,fcnAbsMono);
			break;
		case 2:
			dOpen = dRootBracket(pkdn,dCrit,fcnAbsQuad);
			break;
		case 3:
			dOpen = dRootBracket(pkdn,dCrit,fcnAbsOct);
			break;
		case 4:
			dOpen = dRootBracket(pkdn,dCrit,fcnAbsHex);
			break;
		default:
			assert(0);
			}
		}
	else if (iOpenType == OPEN_JOSH) {
		/*
		 ** Set opening criterion to an approximation of Josh's theta.
		 */
		dOpen = 2/sqrt(3.0)*pkdn->mom.Bmax/dCrit;
		if (dOpen < pkdn->mom.Bmax) dOpen = pkdn->mom.Bmax;
		}
	else {
		/*
		 ** Set opening criterion to the minimal, i.e., Bmax.
		 */
		dOpen = pkdn->mom.Bmax;
		}
	return(dOpen);
	}


void pkdUpPass(PKD pkd,int iCell,int iOpenType,double dCrit,
			   int iOrder, int bGravity)
{
	KDN *c;
	PARTICLE *p;
	int l,u,pj,j;
	double dOpen;

	c = pkd->kdNodes;
	p = pkd->pStore;
	l = c[iCell].pLower;
	u = c[iCell].pUpper;
	if (c[iCell].iDim >= 0) {
		pkdUpPass(pkd,LOWER(iCell),iOpenType,dCrit,iOrder, bGravity);
		pkdUpPass(pkd,UPPER(iCell),iOpenType,dCrit,iOrder, bGravity);
		pkdCombine(&c[LOWER(iCell)],&c[UPPER(iCell)],&c[iCell]);
		}
	else {
		c[iCell].fMass = 0.0;
		c[iCell].fSoft = 0.0;
		for (j=0;j<3;++j) {
			/*
			 ** We initialize the bounds to these extreme crazy values so
			 ** that any particle will set the bounds correctly, if there 
			 ** are no particles in the loop then these remain the bounds.
			 */
			c[iCell].bnd.fMin[j] = FLOAT_MAXVAL;
			c[iCell].bnd.fMax[j] = -FLOAT_MAXVAL;
			c[iCell].bndBall.fMin[j] = FLOAT_MAXVAL;
			c[iCell].bndBall.fMax[j] = -FLOAT_MAXVAL;
#ifdef LONGRANGESTEP
			c[iCell].bndDt.vMin[j] = FLOAT_MAXVAL;
			c[iCell].bndDt.vMax[j] = -FLOAT_MAXVAL;
#endif
			c[iCell].r[j] = 0.0;
			}
#ifdef LONGRANGESTEP
		c[iCell].bndDt.cMax = -FLOAT_MAXVAL;
#endif
		for (pj=l;pj<=u;++pj) {
#ifdef LONGRANGESTEP
		        if (p[pj].c > c[iCell].bndDt.cMax)
			        c[iCell].bndDt.cMax = p[pj].c;
#endif
			for (j=0;j<3;++j) {
				if (p[pj].r[j] < c[iCell].bnd.fMin[j])
					c[iCell].bnd.fMin[j] = p[pj].r[j];
				if (p[pj].r[j] > c[iCell].bnd.fMax[j])
					c[iCell].bnd.fMax[j] = p[pj].r[j];

				if (p[pj].r[j]-p[pj].fBallMax < c[iCell].bndBall.fMin[j])
					c[iCell].bndBall.fMin[j] = p[pj].r[j]-p[pj].fBallMax;
				if (p[pj].r[j]+p[pj].fBallMax > c[iCell].bndBall.fMax[j])
					c[iCell].bndBall.fMax[j] = p[pj].r[j]+p[pj].fBallMax;
#ifdef LONGRANGESTEP
				if (p[pj].vPred[j] < c[iCell].bndDt.vMin[j])
					c[iCell].bndDt.vMin[j] = p[pj].vPred[j];
				if (p[pj].vPred[j] > c[iCell].bndDt.vMax[j])
					c[iCell].bndDt.vMax[j] = p[pj].vPred[j];
#endif

				}
			/*
			 ** Find center of mass and total mass and mass weighted softening.
			 */
			c[iCell].fMass += p[pj].fMass;
			c[iCell].fSoft += p[pj].fMass*p[pj].fSoft;
			for (j=0;j<3;++j) {
				c[iCell].r[j] += p[pj].fMass*p[pj].r[j];
				}
		        }
#ifdef LONGRANGESTEP
		assert( !isnan(c[iCell].bndDt.cMax) );
		DIAGDIST2(c[iCell].bndDt.drMax2,c[iCell].bnd.fMin,c[iCell].bnd.fMax);
#endif
		if (c[iCell].fMass > 0) {
			for (j=0;j<3;++j) {
				c[iCell].r[j] /= c[iCell].fMass;
				}
			c[iCell].fSoft /= c[iCell].fMass;
			}
		}
	/*
	 ** Calculate multipole moments.
	 */
	if(bGravity) {
	    pkdCalcCell(pkd,&c[iCell],c[iCell].r,iOrder,&c[iCell].mom);
	    dOpen = pkdCalcOpen(&c[iCell],iOpenType,dCrit,iOrder);
	    c[iCell].fOpen2 = dOpen*dOpen;
	    }
	}


/*
 ** JST's Select
 */
void pkdSelect(PKD pkd,int d,int k,int l,int r)
{
	PARTICLE *p,t;
	FLOAT v;
	int i,j;

	p = pkd->pStore;
	while (r > l) {
		v = p[k].r[d];
		t = p[r];
		p[r] = p[k];
		p[k] = t;
		i = l - 1;
		j = r;
		while (1) {
			while (i < j) if (p[++i].r[d] >= v) break;
			while (i < j) if (p[--j].r[d] <= v) break;
			t = p[i];
			p[i] = p[j];
			p[j] = t;
			if (j <= i) break;
			}
		p[j] = p[i];
		p[i] = p[r];
		p[r] = t;
		if (i >= k) r = i - 1;
		if (i <= k) l = i + 1;
		}
	}


int NumBinaryNodes(PKD pkd,int nBucket,int pLower,int pUpper)
{
	BND bnd;
	int i,j,m,d,l,u;
	FLOAT fSplit;
	int bGoodBounds;	/* is the cell a finite size? */
	PARTICLE *pStore = pkd->pStore;

	if (pLower > pUpper) return(0);
	else {
		/*
		 ** We need to find the bounding box.
		 */
		for (j=0;j<3;++j) {
			bnd.fMin[j] = pStore[pLower].r[j];
			bnd.fMax[j] = pStore[pLower].r[j];
			}
		for (i=pLower+1;i<=pUpper;++i) {
			for (j=0;j<3;++j) {
				if (pStore[i].r[j] < bnd.fMin[j]) 
					bnd.fMin[j] = pStore[i].r[j];
				else if (pStore[i].r[j] > bnd.fMax[j])
					bnd.fMax[j] = pStore[i].r[j];
				}
			}	
		bGoodBounds = 0;
		for (j=0;j<3;++j) {
			if (bnd.fMax[j] > bnd.fMin[j])
			    bGoodBounds = 1;
			}
		if ((pUpper-pLower+1 > nBucket) && bGoodBounds) {
			/*
			 ** Now we need to determine the longest axis.
			 */
			d = 0;
			for (j=1;j<3;++j) {
				if (bnd.fMax[j]-bnd.fMin[j] > bnd.fMax[d]-bnd.fMin[d]) {
					d = j;
					}
				}
			/*
			 ** Now we do the split.
			 */
			fSplit = 0.5*(bnd.fMin[d]+bnd.fMax[d]);
			m = pkdUpperPart(pkd,d,fSplit,pLower,pUpper);
			/*
			 ** Recursive node count.
			 */
			l = NumBinaryNodes(pkd,nBucket,pLower,m-1);
			u = NumBinaryNodes(pkd,nBucket,m,pUpper);
			/*
			 ** Careful, this assert only applies when we are doing the
			 ** squeezing!
			 */
			mdlassert(pkd->mdl,l > 0 && u > 0);
			return(l+u+1);
			}
		else {
			return(1);
			}
		}
	}


int BuildBinary(PKD pkd,int nBucket,int pLower,int pUpper,int iOpenType,
		double dCrit,int iOrder, int bGravity)
{
	KDN *pkdn;
	int i,j,m,d,c;
	FLOAT fm;
	double dOpen;
	int bGoodBounds;	/* Is the cell a finite size? */
	BND bnd, bndBall;
	PARTICLE *pStore;

	if (pLower > pUpper) return(-1);
	else {
		/*
		 ** Grab a cell from the cell storage!
		 */
		mdlassert(pkd->mdl,pkd->iFreeCell < pkd->nNodes);
		c = pkd->iFreeCell++;
		pkdn = &pkd->kdNodes[c];
		pkdn->pLower = pLower;
		pkdn->pUpper = pUpper;
		pStore = pkd->pStore;
		/*
		 ** We need to find the bounding box.
		 */
		for (j=0;j<3;++j) {
			bnd.fMin[j] = pStore[pLower].r[j];
			bnd.fMax[j] = pStore[pLower].r[j];
			bndBall.fMin[j] = pStore[pLower].r[j]-pStore[pLower].fBallMax;
			bndBall.fMax[j] = pStore[pLower].r[j]+pStore[pLower].fBallMax;
			}
		for (i=pLower+1;i<=pUpper;++i) {
			for (j=0;j<3;++j) {
				if (pStore[i].r[j] < bnd.fMin[j]) 
					bnd.fMin[j] = pStore[i].r[j];
				else if (pStore[i].r[j] > bnd.fMax[j])
					bnd.fMax[j] = pStore[i].r[j];
				if (pStore[i].r[j]-pStore[i].fBallMax < bndBall.fMin[j]) 
					bndBall.fMin[j] = pStore[i].r[j]-pStore[i].fBallMax;
				if (pStore[i].r[j]+pStore[i].fBallMax > bndBall.fMax[j])
					bndBall.fMax[j] = pStore[i].r[j]+pStore[i].fBallMax;
				}
			}	
		bGoodBounds = 0;
		for (j=0;j<3;++j) {
			if (bnd.fMax[j] > bnd.fMin[j])
			    bGoodBounds = 1;
			}

		pkdn->bnd = bnd;
		pkdn->bndBall = bndBall;
		
		if ((pUpper-pLower+1 > nBucket) && bGoodBounds) {
			/*
			 ** Now we need to determine the longest axis.
			 */
			d = 0;
			for (j=1;j<3;++j) {
				if (bnd.fMax[j]-bnd.fMin[j] > 
					bnd.fMax[d]-bnd.fMin[d]) {
					d = j;
					}
				}
			pkdn->iDim = d;
			/*
			 ** Now we do the split.
			 */
			pkdn->fSplit = 0.5*(bnd.fMin[d]+bnd.fMax[d]);
			m = pkdUpperPart(pkd,d,pkdn->fSplit,pLower,pUpper);
			pkdn->iLower = BuildBinary(pkd,nBucket,pLower,m-1,
						   iOpenType, dCrit,iOrder,
						   bGravity);
			pkdn->iUpper = BuildBinary(pkd,nBucket,m,pUpper,
						   iOpenType, dCrit,iOrder,
						   bGravity);
			/*
			 ** Careful, this assert only applies when we are doing the
			 ** squeezing!
			 */
			mdlassert(pkd->mdl,pkdn->iLower != -1 && pkdn->iUpper != -1);
			/*
			 ** Now calculate the mass, center of mass and mass weighted
			 ** softening radius.
			 */
			pkdn->fMass = 0.0;
			pkdn->fSoft = 0.0;
			for (j=0;j<3;++j) {
				pkdn->r[j] = 0.0;
				}
			if (pkdn->iLower != -1) {
				fm = pkd->kdNodes[pkdn->iLower].fMass;
				pkdn->fMass += fm;
				pkdn->fSoft += fm*pkd->kdNodes[pkdn->iLower].fSoft;
				for (j=0;j<3;++j) {
					pkdn->r[j] += fm*pkd->kdNodes[pkdn->iLower].r[j];
					}
				}
			if (pkdn->iUpper != -1) {
				fm = pkd->kdNodes[pkdn->iUpper].fMass;
				pkdn->fMass += fm;
				pkdn->fSoft += fm*pkd->kdNodes[pkdn->iUpper].fSoft;
				for (j=0;j<3;++j) {
					pkdn->r[j] += fm*pkd->kdNodes[pkdn->iUpper].r[j];
					}
				}
			if (pkdn->fMass > 0) {
				pkdn->fSoft /= pkdn->fMass;
				for (j=0;j<3;++j) {
					pkdn->r[j] /= pkdn->fMass;
					}
				}
			}
		else {
			pkdn->iDim = -1; /* it is a bucket! */
			pkdn->fSplit = 0.0;
			pkdn->iLower = -1;
			pkdn->iUpper = -1;
			/*
			 ** Calculate the bucket quantities.
			 */
			pkdn->fMass = 0.0;
			pkdn->fSoft = 0.0;
			for (j=0;j<3;++j) {
				pkdn->r[j] = 0.0;
				}
			for (i=pkdn->pLower;i<=pkdn->pUpper;++i) {
				fm = pStore[i].fMass;
				pkdn->fMass += fm;
				pkdn->fSoft += fm*pStore[i].fSoft;
				for (j=0;j<3;++j) {
					pkdn->r[j] += fm*pStore[i].r[j];
					}
				}
			if (pkdn->fMass > 0) {
				pkdn->fSoft /= pkdn->fMass;
				for (j=0;j<3;++j) {
					pkdn->r[j] /= pkdn->fMass;
					}
				}
			}
		/*
		 ** Calculate multipole moments.
		 */
		if(bGravity) {
		    pkdCalcCell(pkd,pkdn,pkdn->r,iOrder,&pkdn->mom);
		    dOpen = pkdCalcOpen(pkdn,iOpenType,dCrit,iOrder);
		    pkdn->fOpen2 = dOpen*dOpen;
		    }
		return(c);
		}
	}


void pkdThreadTree(PKD pkd,int iCell,int iNext)
{
	int l,u;

	if (iCell == -1) return;
	else if (pkd->kdNodes[iCell].iDim != -1) {
		l = pkd->kdNodes[iCell].iLower;
		u = pkd->kdNodes[iCell].iUpper;
		if (u == -1) {
			mdlassert(pkd->mdl,l != -1);
			pkdThreadTree(pkd,l,iNext);
			}
		else if (l == -1) {
			mdlassert(pkd->mdl,u != -1);
			pkdThreadTree(pkd,u,iNext);
			/* 
			 ** It is convenient to change the "down" pointer in this case.
			 */
			pkd->kdNodes[iCell].iLower = u;
			}
		else {
			pkdThreadTree(pkd,l,u);
			pkdThreadTree(pkd,u,iNext);
			}
		pkd->kdNodes[iCell].iUpper = iNext;
		}
	else {
		pkd->kdNodes[iCell].iLower = -1;	/* Just make sure! */
		pkd->kdNodes[iCell].iUpper = iNext;
		}
	}

/*
 * Builds a binary tree by splitting the largest spatial axis in half
 * The bounds for each subcell are recalculated (squeezed) to be exact
 *
 */
void pkdBuildBinary(PKD pkd,int nBucket,int iOpenType,double dCrit,
		    int iOrder,int bTreeActiveOnly,int bGravity, KDN *pRoot)
{
	int bEmpty = 0;

	/*
	 ** Make sure the particles are in Active/Inactive order.
	 */
	pkdActiveTypeOrder(pkd, TYPE_TREEACTIVE );
	if (pkd->kdNodes) {
		/*
		 ** Close caching space and free up nodes.
		 */
		mdlFinishCache(pkd->mdl,CID_CELL);
		mdlFree(pkd->mdl,pkd->kdNodes);
		}
	/*
	 ** First problem is to figure out how many cells we need to
	 ** allocate. We need to do this in advance because we need to
	 ** synchronize to allocate the memory with mdlMalloc().
	 */
	if (bTreeActiveOnly) {
		pkd->nNodes = NumBinaryNodes(pkd,nBucket,0,pkd->nTreeActive-1);
		}
	else {
		pkd->nNodes = NumBinaryNodes(pkd,nBucket,0,pkd->nLocal-1);
		}
	/*
	 ** We need at least one particle per processor.
	mdlassert(pkd->mdl,pkd->nNodes > 0);
	 ** Let's relax this --trq
	 */
	if(pkd->nNodes == 0) {
	    pkd->nNodes = 1;
	    bEmpty = 1;
	}

	/*
	 ** Need to allocate a special extra cell that we will use to calculate
	 ** the acceleration on an arbitrary point in space.
	 */
	pkd->kdNodes = mdlMalloc(pkd->mdl,(pkd->nNodes + 1)*sizeof(KDN));
	mdlassert(pkd->mdl,pkd->kdNodes != NULL);
	/*
	 ** Now we really build the tree.
	 */
	pkd->iFreeCell = 0;
	if (bTreeActiveOnly) {
		pkd->iRoot = BuildBinary(pkd,nBucket,0,pkd->nTreeActive-1,
					 iOpenType,dCrit,iOrder, bGravity);
		}
	else {
		pkd->iRoot = BuildBinary(pkd,nBucket,0,pkd->nLocal-1,
					 iOpenType,dCrit,iOrder, bGravity);
		}
	mdlassert(pkd->mdl,bEmpty || pkd->iFreeCell == pkd->nNodes);
	if(bEmpty) {
	    KDN *pkdn;
	    int j;
	    
	    /*
	     * Set up an empty bucket.
	     */
	    pkd->iRoot = 0;
	    pkdn = &pkd->kdNodes[0];
	    pkdn->pLower = 0;
	    pkdn->pUpper = -1;
	    for (j=0;j<3;++j) {
	        pkdn->bnd.fMin[j] = FLOAT_MAXVAL;
		pkdn->bnd.fMax[j] = -FLOAT_MAXVAL;
		pkdn->bndBall.fMin[j] = FLOAT_MAXVAL;
		pkdn->bndBall.fMax[j] = -FLOAT_MAXVAL;
	        }
#ifdef LONGRANGESTEP
	    pkdn->bndDt.cMax = -FLOAT_MAXVAL;
#endif
	    pkdn->iDim = -1; /* it is a bucket! */
	    pkdn->fSplit = 0.0;
	    pkdn->iLower = -1;
	    pkdn->iUpper = -1;
	    pkdn->fMass = 0.0;
	    pkdn->fSoft = 0.0;
	    for (j=0;j<3;++j) {
	        pkdn->r[j] = 0.0;
		}
	    if(bGravity) {
	        pkdCalcCell(pkd,pkdn,pkdn->r,iOrder,&pkdn->mom);
		pkdn->fOpen2 = 0.0;
	      }
	    }
	/*
	 ** Thread the tree.
	 */
	pkdThreadTree(pkd,pkd->iRoot,-1);
	*pRoot = pkd->kdNodes[pkd->iRoot];
	/*
	 ** Finally activate a read only cache for remote access.
	 */
	mdlROcache(pkd->mdl,CID_CELL,pkd->kdNodes,sizeof(KDN),pkdNodes(pkd));
	}


/*
 * Builds a binary tree by splitting on the largest spatial axis
 * so as to equal divide the particle by number
 * After the tree is built the bounds for each subcell are later recalculated (squeezed) to be exact
 *
 */
void pkdBuildLocal(PKD pkd,int nBucket,int iOpenType,double dCrit,
		   int iOrder,int bTreeActiveOnly,int bGravity, KDN *pRoot)
{
	int l,n,i,d,m,j,diff;
	KDN *c;
	char ach[256];
	BND bndDum;
	int bGoodBounds;

	/*
	 ** Make sure the particles are in Active/Inactive order.
	 */
	pkdActiveTypeOrder(pkd, TYPE_TREEACTIVE );
	pkd->nBucket = nBucket;
	if (bTreeActiveOnly) n = pkd->nTreeActive;
	else n = pkd->nLocal;
	pkd->nLevels = 1;
	l = 1;
	while (n > nBucket) {
		n = n>>1;
		l = l<<1;
		++pkd->nLevels;
		}
	pkd->nSplit = l;
	pkd->nNodes = l<<1;
	if (pkd->kdNodes) {
		/*
		 ** Close caching space and free up nodes.
		 */
		mdlFinishCache(pkd->mdl,CID_CELL);
		mdlFree(pkd->mdl,pkd->kdNodes);
		}
	/*
	 ** Need to allocate a special extra cell that we will use to calculate
	 ** the acceleration on an arbitrary point in space.
	 */
	pkd->kdNodes = mdlMalloc(pkd->mdl,(pkd->nNodes + 1)*sizeof(KDN));
	mdlassert(pkd->mdl,pkd->kdNodes != NULL);
	pkd->iFreeCell = pkd->nNodes;
	sprintf(ach,"nNodes:%d nSplit:%d nLevels:%d nBucket:%d\n",
			pkd->nNodes,pkd->nSplit,pkd->nLevels,nBucket);
	mdlDiag(pkd->mdl,ach);
	/*
	 ** Set up ROOT node
	 */
	c = pkd->kdNodes;
	pkd->iRoot = ROOT;
	c[pkd->iRoot].pLower = 0;
	if (bTreeActiveOnly) {
		c[pkd->iRoot].pUpper = pkd->nTreeActive-1;
		}
	else {
		c[pkd->iRoot].pUpper = pkd->nLocal-1;
		}
	/*
	 ** determine the local bound of the particles.
	 */
	if (bTreeActiveOnly) {
	    pkdCalcBound(pkd,&bndDum,&bndDum,&c[pkd->iRoot].bnd,&bndDum,NULL);
		}
	else {
	    pkdCalcBound(pkd,&c[pkd->iRoot].bnd,&bndDum,&bndDum,&bndDum,NULL);
		}
	i = pkd->iRoot;
	while (1) {
		bGoodBounds = 0;
		for (j=0;j<3;++j) {
			if (c[i].bnd.fMax[j] > c[i].bnd.fMin[j])
				bGoodBounds = 1;
			}
		if (i < pkd->nSplit && (c[i].pUpper - c[i].pLower) > 0 && bGoodBounds) {
			d = 0;
			for (j=1;j<3;++j) {
				if (c[i].bnd.fMax[j]-c[i].bnd.fMin[j] > 
					c[i].bnd.fMax[d]-c[i].bnd.fMin[d]) d = j;
				}
			c[i].iDim = d;
			m = (c[i].pLower + c[i].pUpper)/2;
			pkdSelect(pkd,d,m,c[i].pLower,c[i].pUpper);
			c[i].fSplit = pkd->pStore[m].r[d];
			c[i].iLower = LOWER(i);
			c[i].iUpper = UPPER(i);
			c[LOWER(i)].bnd = c[i].bnd;
			c[LOWER(i)].bnd.fMax[d] = c[i].fSplit;
			/* Unknown here; fix in pkdUpPass
			c[LOWER(i)].bndBall = c[i].bndBall; */
			c[LOWER(i)].pLower = c[i].pLower;
			c[LOWER(i)].pUpper = m;
			c[UPPER(i)].bnd = c[i].bnd;
			c[UPPER(i)].bnd.fMin[d] = c[i].fSplit;
			/* Unknown here; fix in pkdUpPass
			c[UPPER(i)].bndBall = c[i].bndBall; */
			c[UPPER(i)].pLower = m+1;
			c[UPPER(i)].pUpper = c[i].pUpper;
			diff = (m-c[i].pLower+1)-(c[i].pUpper-m);
			mdlassert(pkd->mdl,diff == 0 || diff == 1);
			i = LOWER(i);
			}
		else {
			c[i].iDim = -1;		/* to indicate a bucket! */
			c[i].iLower = -1;
			c[i].iUpper = -1;
			SETNEXT(i);
			if (i == pkd->iRoot) break;
			}
		}
	pkdUpPass(pkd,pkd->iRoot,iOpenType,dCrit,iOrder, bGravity); 
	/*
	 ** Thread the tree.
	 */
	pkdThreadTree(pkd,pkd->iRoot,-1);
	*pRoot = c[pkd->iRoot];
	/*
	 ** Finally activate a read only cache for remote access.
	 */
	mdlROcache(pkd->mdl,CID_CELL,pkd->kdNodes,sizeof(KDN),pkdNodes(pkd));
	}


void pkdBucketWeight(PKD pkd,int iBucket,FLOAT fWeight)
{
	KDN *pbuc;
	int pj;
	
	pbuc = &pkd->kdNodes[iBucket];
	for (pj=pbuc->pLower;pj<=pbuc->pUpper;++pj) {
		if (TYPEQueryACTIVE(&(pkd->pStore[pj])))
			pkd->pStore[pj].fWeight = fWeight;
		}
	}


void pkdColorCell(PKD pkd,int iCell,FLOAT fColor)
{
#ifdef COLORCODE
	KDN *pkdn;
	int pj;
	
	pkdn = &pkd->kdNodes[iCell];
	for (pj=pkdn->pLower;pj<=pkdn->pUpper;++pj) {
		pkd->pStore[pj].fColor = fColor;
		}
#endif	
	}

void
pkdGravAll(PKD pkd,int nReps,int bPeriodic,int iOrder,int bEwald,int iEwOrder,
	   double fEwCut,double fEwhCut, int bComove, double dRhoFac,
	   int bDoSun,double dSunSoft, double *aSun,int *nActive,
		   double *pdPartSum,double *pdCellSum,double *pdSoftSum,CASTAT *pcs,
		   double *pdFlop)
{
	KDN *c = pkd->kdNodes;
	int iCell,n;
	FLOAT fWeight;
	double dFlopI, dFlopE;
	int i,j;
	BND bndActive;
	BND bndTmp;
	char achDiag[256];
	const int bDetailTimer = 0; /* Change to 1 to get detailed timing on
				     * interacts, walks, and Ewald
				     */

	*pdFlop = 0.0;
	pkdClearTimer(pkd,1);
	pkdClearTimer(pkd,2);
	pkdClearTimer(pkd,3);
	/*
	 ** Set up Ewald tables and stuff.
	 */
	if (bPeriodic && bEwald) {
	    pkdEwaldInit(pkd,fEwhCut,iEwOrder);	/* ignored in Flop count! */
	}
	/*
	 ** Start particle caching space (cell cache already active).
	 */
	mdlROcache(pkd->mdl,CID_PARTICLE,pkd->pStore,sizeof(PARTICLE),
			   pkdLocal(pkd));
	/*
	 ** Walk over the local buckets!
	 */
	*nActive = 0;
	*pdPartSum = 0.0;
	*pdCellSum = 0.0;
	*pdSoftSum = 0.0;
	iCell = pkd->iRoot;
	if(!bDetailTimer) pkdStartTimer(pkd,2);
	while (iCell != -1) {
		if (c[iCell].iLower != -1) {
			iCell = c[iCell].iLower;
			continue;
			}
		n = 0;
		for(j = 0; j < 3; j++) {
		    bndActive.fMin[j] = FLOAT_MAXVAL;
		    bndActive.fMax[j] = -FLOAT_MAXVAL;
		    }
		for (i=c[iCell].pLower;i<=c[iCell].pUpper;++i) {
			if (TYPEQueryACTIVE(&(pkd->pStore[i]))) {
			    ++n;
			    for (j=0;j<3;++j) {
					if (pkd->pStore[i].r[j] < bndActive.fMin[j]) 
						bndActive.fMin[j] = pkd->pStore[i].r[j];
					if (pkd->pStore[i].r[j] > bndActive.fMax[j])
						bndActive.fMax[j] = pkd->pStore[i].r[j];
					}
			    }
			}
		if (n > 0) {
			/* 
			 * set bounds to bounds of active particles.
			 */
			bndTmp = c[iCell].bnd;
			c[iCell].bnd = bndActive;
			/*
			 ** Calculate gravity on this bucket.
			 */
			if(bDetailTimer) pkdStartTimer(pkd,1);
			pkdBucketWalk(pkd,iCell,nReps,iOrder); /* ignored in Flop count! */
			if(bDetailTimer) pkdStopTimer(pkd,1);
			c[iCell].bnd = bndTmp;
			*nActive += n;
			*pdPartSum += n*pkd->nPart + 
				n*(2*(c[iCell].pUpper-c[iCell].pLower) - n + 1)/2;
			*pdCellSum += n*pkd->nCellNewt;
			*pdSoftSum += n*pkd->nCellSoft;
			if(bDetailTimer) pkdStartTimer(pkd,2);
			dFlopI = pkdBucketInteract(pkd,iCell,iOrder);
			*pdFlop += dFlopI;
			if(bDetailTimer) pkdStopTimer(pkd,2);
			/*
			 * Now do Ewald part.
			 */
			if (bPeriodic && bEwald) {
				if(bDetailTimer) pkdStartTimer(pkd,3);
				dFlopE = pkdBucketEwald(pkd,iCell,nReps,fEwCut,iEwOrder);
				*pdFlop += dFlopE;
				if(bDetailTimer) pkdStopTimer(pkd,3);
				}
			else {
			    dFlopE = 0.0;
			    }
			fWeight = dFlopI + dFlopE;
			if(bComove && !bPeriodic) {
			    /* 
			     * Add gravity from the rest of the
			     * Universe.  This adds force and
			     * potential from a uniform (negative!)
			     * density sphere.
			     */
			    KDN *pkdn = &pkd->kdNodes[iCell];
			    int nP = pkdn->pUpper - pkdn->pLower + 1;
			    PARTICLE *p = &pkd->pStore[pkdn->pLower];
			    int jP;
			    
			    for(jP=0;jP<nP;++jP) {
				int k;
				double r2;
				
				if (!TYPEQueryACTIVE(&(p[jP]))) continue;
				r2 = 0.0;
				for(k = 0; k < 3; k++) {
				    r2 += p[jP].r[k]*p[jP].r[k];
				    p[jP].a[k] += dRhoFac*p[jP].r[k];
				    }
				p[jP].fPot -= 0.5*dRhoFac*r2;
				}
			    }
			/*
			 ** pkdBucketWeight, only updates the weights of the active
			 ** particles. Although this isn't really a requirement it
			 ** might be a good idea, if weights correspond to different 
			 ** tasks at different times.
			 */
			pkdBucketWeight(pkd,iCell,fWeight);
			}
		iCell = c[iCell].iUpper;
		}
	if(!bDetailTimer) pkdStopTimer(pkd,2);
	if (bDoSun) {
		const double dTinyBox = 1e-14;
		int iDummy = pkd->nStore;
		/*
		 ** Calculate the indirect interaction for solar system problems.
		 ** we need a "dummy" particle that is pointed to by the cell,
		 ** pkd->iFreeCell which also needs its bounds set correctly for 
		 ** this to work.
		 ** Don't allow periodic BCs at the moment.
		 */
		mdlassert(pkd->mdl,nReps == 0);
		mdlassert(pkd->mdl,bPeriodic == 0);
		for (j=0;j<3;++j) {
			pkd->pStore[iDummy].r[j] = 0.0;
			pkd->pStore[iDummy].a[j] = 0.0;
			c[pkd->iFreeCell].bnd.fMin[j] = -dTinyBox;
			c[pkd->iFreeCell].bnd.fMax[j] = dTinyBox;
			}
		TYPESet(&(pkd->pStore[iDummy]),TYPE_ACTIVE);
		pkd->pStore[iDummy].fPot = 0;
		pkd->pStore[iDummy].fSoft = dSunSoft;
		c[pkd->iFreeCell].fSoft = dSunSoft;
		c[pkd->iFreeCell].pLower = iDummy;
		c[pkd->iFreeCell].pUpper = iDummy;
		c[pkd->iFreeCell].iDim = -1;	/* It is really a bucket! */

		pkdBucketWalk(pkd,pkd->iFreeCell,nReps,iOrder);
		pkdBucketInteract(pkd,pkd->iFreeCell,iOrder);

		/*
		 ** Now we should have the indirect acceleration contained in the 
		 ** particle pointed to by iDummy. Now we have to bring this acceleration
		 ** back up to the master level. We can then take care of it with an
		 ** external potential call.
		 */
		for (j=0;j<3;++j) {
			aSun[j] = pkd->pStore[iDummy].a[j];
			}
		}
	/*
	 ** Get caching statistics.
	 */
	pcs->dcNumAccess = mdlNumAccess(pkd->mdl,CID_CELL);
	pcs->dcMissRatio = mdlMissRatio(pkd->mdl,CID_CELL);
	pcs->dcCollRatio = mdlCollRatio(pkd->mdl,CID_CELL);
	pcs->dcMinRatio = mdlMinRatio(pkd->mdl,CID_CELL);
	pcs->dpNumAccess = mdlNumAccess(pkd->mdl,CID_PARTICLE);
	pcs->dpMissRatio = mdlMissRatio(pkd->mdl,CID_PARTICLE);
	pcs->dpCollRatio = mdlCollRatio(pkd->mdl,CID_PARTICLE);
	pcs->dpMinRatio = mdlMinRatio(pkd->mdl,CID_PARTICLE);
	/*
	 ** Stop particle caching space.
	 */
	mdlFinishCache(pkd->mdl,CID_PARTICLE);
	sprintf(achDiag, "nMaxPart: %d, nMaxSoftCell: %d, nMaxNewtCell: %d\n",
		pkd->nMaxPart, pkd->nMaxCellSoft, pkd->nMaxCellNewt);
	mdlDiag(pkd->mdl, achDiag);
	}

#ifdef  RADIATIVEBOX
/* Walk down the tree within a node setting the LW flux in each cell equal to the maximum of the typical flux in that cell of the typical flux in its parent cell, CC*/
void pkdLocalFinishLWTree(PKD pkd, int iCell,
			  double fPrevLW, /* Total LW luminosity in parent cell*/
			  FLOAT *PrevcLumLW /* Luminosity-weighted center in parent cell*/
			  )        
{
	PARTICLE *p;
	KDN *pkdn,*pkdnL;
	int nPart,pj,iRight,iLeft,j;
	double fDistanceCell2 = 0, fDistance2 = 0, fPrevAveLW, fDistance2min, smooth2;
	nPart = pkd->nPart;
	p = pkd->pStore;
	pkdn = &pkd->kdNodes[iCell];

	/*	Determine flux at center of mass of cell from average source of in the parent cell */
	for (j = 0; j<3; ++j){
	  fDistanceCell2 =+ (PrevcLumLW[j] - pkdn->r[j])*(PrevcLumLW[j] - pkdn->r[j]);
	}
	if (fDistanceCell2 < pkdn->mom.gmom) fDistanceCell2 = pkdn->mom.gmom; /* At minimum, this average distance between gas mass center of child cell and flux center of parent cell should be the moment of the gas in the child cell*/
	fPrevAveLW = fPrevLW/fDistanceCell2; /* Calculated typical flux with radiative source being in the parent cell*/
	if (pkdn->mom.gmom == 0 || fPrevAveLW > pkdn->mom.fLW/pkdn->mom.gmom){
	  /* If using the radiation from the parent cell would typically provide more flux, set total luminosity and average source position equal to that in parent cell */
	  pkdn->mom.fLW = fPrevLW;
	  for (j=0;j<3;++j) {
	    pkdn->mom.cLumLW[j] = PrevcLumLW[j];
	  }
	} 
	if (pkdn->iLower != -1) {
	  pkdnL =  &pkd->kdNodes[pkdn->iLower];
	  iLeft = pkdn->iLower;
	  pkdLocalFinishLWTree(pkd,pkdn->iLower,pkdn->mom.fLW,pkdn->mom.cLumLW);
	  if (pkdnL->iUpper != -1){
	    iRight = pkdnL->iUpper;
	    pkdLocalFinishLWTree(pkd,pkdnL->iUpper,pkdn->mom.fLW,pkdn->mom.cLumLW);	      
	  }
	}
	else {
	  /* If you are in a bucket, calculate the flux at each gas particle from the total luminosity at the average position of the source */
	  for (pj=pkdn->pLower;pj<=pkdn->pUpper;++pj) {
	    if (TYPETest(&(p[pj]),TYPE_GAS))
	      {
		fDistance2 = 0;
		for (j=0;j<3;++j) {
		  fDistance2 += pow((p[pj].r[j] - pkdn->mom.cLumLW[j]),2);
		}		
		fDistance2min = p[pj].fSoft*p[pj].fSoft*0.25;
		if (fDistance2 < fDistance2min) {
		  fDistance2 = fDistance2min;
		}
		assert(fDistance2 > 0);
		p[pj].CoolParticle.dLymanWerner =  pkdn->mom.fLW/(4.0*M_PI*fDistance2);
	      }	    
	  }
	}
}

/* Walk down the top level of the tree setting the LW flux in each cell equal to the maximum of the typical flux in that cell of the typical flux in its parent cell, CC*/
void pkdFinishLWTree(PKD pkd)
{
       KDN *pbuc,*pkdn,*pkdnL, *pkdnR;
       int iCell,iLeft,iRight,id,idLL,idLU,idRL,idRU,j;
       double fDistanceCell2, fPrevAveLW;

	/*
	** Walk the top tree first, finding local trees to
	** continue walking.
	*/
	iCell = ROOT;
	pkdn = &pkd->kdTop[iCell];
	while (1) {
	  id = pkd->kdTop[iCell].pLower;
	  if (id == pkd->idSelf) {
	    /* Start parallel tree traversal*/
	    pkdLocalFinishLWTree(pkd,pkd->iRoot,pkdn->mom.fLW,pkdn->mom.cLumLW);
	    SETNEXT(iCell);
	      if (iCell == ROOT)   break;
	  }
	  else if (id >= 0) {
	    /*Ignore, other processors will take care of */
	    SETNEXT(iCell);
	      if (iCell == ROOT) break;
	   }
	  else {	
	    pkdn = &pkd->kdTop[iCell];   
	    iLeft = LOWER(iCell);
	    iRight = UPPER(iCell);

	    pkdnL = &pkd->kdTop[iLeft];
	    idLL = pkd->kdTop[iLeft].pLower;
	    idLU = pkd->kdTop[iLeft].pUpper;
	    /*	Determine flux at center of mass of cell from average source of in the parent cell */
	    fDistanceCell2 = 0;
	    for (j = 0; j<3; ++j){
	      fDistanceCell2 =+ (pkdn->mom.cLumLW[j] - pkdnL->r[j])*(pkdn->mom.cLumLW[j] - pkdnL->r[j]);
		}
	    if (fDistanceCell2 < pkdnL->mom.gmom) fDistanceCell2 = pkdnL->mom.gmom;/* At minimum, this average distance between gas mass center of child cell and flux center of parent cell should be the moment of the gas in the child cell*/
	    fPrevAveLW = pkdn->mom.fLW/fDistanceCell2;
	    if (pkdnL->mom.gmom == 0 || fPrevAveLW > pkdnL->mom.fLW/pkdnL->mom.gmom) {
	      /* If using the radiation from the parent cell would typically provide more flux, set total luminosity and average source position equal to that in parent cell */
	      pkdnL->mom.fLW = pkdn->mom.fLW;
	      for (j=0;j<3;++j) {
		pkdnL->mom.cLumLW[j] = pkdn->mom.cLumLW[j];
	      }	      
	    }
	    /* As above for right node*/
	    pkdnR = &pkd->kdTop[iRight];
	    idRL = pkd->kdTop[iRight].pLower;
	    idRU = pkd->kdTop[iRight].pUpper;
	    fDistanceCell2 = 0;
	    for (j = 0; j<3; ++j){
	      fDistanceCell2 =+ (pkdn->mom.cLumLW[j] - pkdnR->r[j])*(pkdn->mom.cLumLW[j] - pkdnR->r[j]);
		}
	    if (fDistanceCell2 < pkdnR->mom.gmom) fDistanceCell2 = pkdnR->mom.gmom;
	    fPrevAveLW = pkdn->mom.fLW/fDistanceCell2;
	    if (pkdnR->mom.gmom == 0 || fPrevAveLW > pkdnR->mom.fLW/pkdnR->mom.gmom){
	      pkdnR->mom.fLW = pkdn->mom.fLW;
	      for (j=0;j<3;++j) {
		pkdnR->mom.cLumLW[j] = pkdn->mom.cLumLW[j];
	      }	      
	    }
	    iCell = LOWER(iCell);
	    if (iCell == ROOT) break; /* If you are at the top of the tree, (make like a tree and) leave*/
	  }
	}
}
#endif

void pkdSunIndirect(PKD pkd,double *aSun,int bDoSun,double dSunMass,double dSunSoft)
{
	PARTICLE *p;
	double t,idt2,a,b;
	int i,j,n;

	p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;++i) {
		if (TYPEQueryACTIVE(&(p[i]))) {
			t = 0;
			for (j=0;j<3;++j) t += p[i].r[j]*p[i].r[j];
			if(t == 0)
				a = b = 0;
			else
				SPLINE(t, (p[i].fSoft + dSunSoft), a, b);
			p[i].fPot -= dSunMass*a;
			idt2 = (p[i].fMass + dSunMass)*b;
			if (idt2 > p[i].dtGrav) p[i].dtGrav = idt2;
			/*
			 ** The way that aSun[j] gets added in here is confusing, but this
			 ** is the correct way! (aSun[] is the acceleration on the Sun).
			 */
			if (bDoSun) {
				b *= dSunMass;
				for (j=0;j<3;++j) {					
					p[i].a[j] -= (aSun[j] + p[i].r[j]*b);
					}				
				}
			else {
				b *= p[i].fMass;
				for (j=0;j<3;++j) {
					p[i].a[j] -= (aSun[j] - p[i].r[j]*b);
					}
				}
			}
		}
	}


void pkdLogHalo(PKD pkd, double Vcirc, double Eps, double Flat)
{
	PARTICLE *p;
	int i,n;

	/* Note: Old: did not include 0.5 in potential -- see B&T 2ed, pg 75 
	const double Vcirc = 194.9848486; in (kpc)/(4.691822922e16s) 128 km/s x 1.5233; (x sqrt(2) now?) 
	const double Eps = 12; kpc 
	const double Flat = 1.0; flattening */
	p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;++i) {
		if (TYPEQueryACTIVE(&(p[i]))) {
			double x = p[i].r[0];
			double y = p[i].r[1];
			double z = p[i].r[2];
			/*
			 **	Do the logarithmic halo potential.
			 */
			double r2 = (x*x + y*y + z*z/(Flat*Flat));
			double C = Vcirc*Vcirc/(r2 + Eps*Eps);	
			p[i].a[0] -= x*C;
			p[i].a[1] -= y*C;
			p[i].a[2] -= z*C/(Flat*Flat);
			p[i].fPot += 0.5*Vcirc*Vcirc*log(r2 + Eps*Eps);
			}
		}
	}


void pkdHernquistSpheroid(PKD pkd)
{
	PARTICLE *p;
	int i,n;

	const double M_s = 3.4e5;	/* in 10^5 M_sun */
	const double c = 0.7;		/* in kpc */
	p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;++i) {
		if (TYPEQueryACTIVE(&(p[i]))) {
			double x = p[i].r[0];
			double y = p[i].r[1];
			double z = p[i].r[2];
			/*
			 **	Do the spheroid potential
			 */
			double r = sqrt(x*x + y*y + z*z);
			double A = 1.0/(r + c);	
			p[i].a[0] -= M_s*A*A*x/r;
			p[i].a[1] -= M_s*A*A*y/r;
			p[i].a[2] -= M_s*A*A*z/r;
			p[i].fPot -= M_s/(r + c);
			}
		}
	}


void pkdNFWSpheroid(PKD pkd, double M_200, double r_200, double c, double dSoft)
{
  PARTICLE *p;
  int i,n;
  /* M_200, r_200, c, and dSoft are now set in param file as
   * dNFWm200, dNFWr200, dNFWconc, and dNFWsoft */
  const double G = 1;
  /* Assuming G = 1 (this sets a timescale) */
  /* TimeUnit = sqrt(kpc^3/(G*1e12 Msun)) = 1.1285945e+09 yr */
  /* Vunit = 2073.8081 km/s */
  
  const double eps = c*dSoft/r_200; 
  
  /* r=r_200 (x=1, cx=c), M=M_200 */
  const double M_const = M_200/
    ( (1./3.)*eps*eps/(1+eps)/(1+eps) 
      - 1/(1+eps) + log((1+c)/(1+eps)) + 1/(1+c) );
  p = pkd->pStore;
  n = pkdLocal(pkd);
  for (i=0;i<n;++i) {
    if (TYPEQueryACTIVE(&(p[i]))) {
      double x = p[i].r[0];
      double y = p[i].r[1];
      double z = p[i].r[2];
      /*
      **	Do the spheroid potential
      */
      double r = sqrt(x*x + y*y + z*z);
      double A, fPot;
      double cx = r*(c/r_200);
      
      if (cx < eps) {
	fPot = G*M_const*c/r_200* ((1./6.)*( cx*cx - eps*eps )/
				   (eps*(1+eps)*(1+eps)) 
				   - (1./3.)*eps/(1+eps)/(1+eps) - 1/(1+eps));
	
	A = G*M_const*(1./3.)/(eps*(1+eps)*(1+eps))
	  *(c*c*c)/(r_200*r_200*r_200);
      }			  
      else {
	fPot = G*M_const*c/r_200 * 
	  (( -(1./3.) *eps*eps/(1+eps)/(1+eps) - eps/(1+eps) 
	     - log((1+cx)/(1+eps)) ) /cx);
	
	A = G*M_const* 
	  ( (1./3.)*eps*eps/(1+eps)/(1+eps) 
	    - 1/(1+eps) + log((1+cx)/(1+eps)) + 1/(1+cx) )
	  /(r*r*r);
      }
      
      /*fprintf(stderr,"%i: %f %f %f  %f %f %f %f\n",p[i].iOrder,x,y,z,r,fPot,-A*r,A); */
      p[i].a[0] -= A*x;
      p[i].a[1] -= A*y;
      p[i].a[2] -= A*z;
      p[i].fPot += fPot;
    }
  }
}

/** JH Feb 4 2004 
elliptical potential using Hernquist model for luminous matter
and Hernquist or NFW model for dark matter */ 

void pkdElliptical(PKD pkd, int bEllipticalDarkNFW)
{
	PARTICLE *p;
	int i,n;
/* if bEllipticalDarkNFW=0 in param file uses hernquist dark matter otherwise uses NFW dark matter */
	const double M_l =1.45e11/6.45e9 ;	/* mass of luminous matter in 6.45e9 M_sun */
	const double a = 2.856/3.5;     /* hernquist scale length for light matter in 6 kpc */
	const double pi=  3.1415; 
	const double M_d= 7.2e11/6.45e9;     /* mass of dark matter in 6.45e9 M_sun */
        const double d= 31.7/3.5;     /* hernquist scale length for dark matter in 6 kpc */
        const double r_s=14.4/3.5 ;     /* NFW scale radius in 6 kpc  */
        const double rho_s= (0.010e9/6.45e9)*3.5*3.5*3.5 ;   /* NFW density at scale radius M_sun/kpc^3 converted to dMsolUnits*/
        const double eps1= 1e-5 ;   /* small softening for hernquist model */
        const double dsoft=1e-5 ;
	p = pkd->pStore;
	n = pkdLocal(pkd);
	/*	printf("bEllipticalDarkNFW= %d", bEllipticalDarkNFW); */
	for (i=0;i<n;++i) {
		if (TYPEQueryACTIVE(&(p[i]))) {
			double x = p[i].r[0];
			double y = p[i].r[1];
			double z = p[i].r[2];
			/*
			 **	Do the spheroid potential
	                 */
			double r = sqrt(x*x + y*y + z*z);
			double A = 1.0/(r + a)  ;	/*light matter component */
			double Ah=  1.0/(r+ d) ;        /*hernquist dark matter component */
			double Anfw= 4*pi*rho_s*r_s*r_s*r_s ;     /* nfw dark matter component */
			if (!bEllipticalDarkNFW) {
			  p[i].a[0] -= M_l*A*A*x/(r + eps1) + M_d*Ah*Ah*x/(r+ eps1);
			  p[i].a[1] -= M_l*A*A*y/(r +eps1)+ M_d*Ah*Ah*y/(r + eps1);
			  p[i].a[2] -= M_l*A*A*z/(r + eps1)+ M_d*Ah*Ah*z/(r +eps1);
			  p[i].fPot -= M_l*A + M_d*Ah ;
			} 
			else {
			      if (r>=dsoft) {
                              p[i].a[0]-= M_l*A*A*x/(r+ eps1) +Anfw*(r_s/(r+r_s)/(r*r) +log(r+r_s)/(r*r) -1/(r*r) -log(r_s)/(r*r))*x/r ; 
			    p[i].a[1]-= M_l*A*A*y/(r+eps1)  +Anfw*(r_s/(r+r_s)/(r*r) +log(r+r_s)/(r*r) -1/(r*r) -log(r_s)/(r*r))*y/r ;   
			    p[i].a[2]-= M_l*A*A*z/(r+ eps1) +Anfw*(r_s/(r+r_s)/(r*r) +log(r+r_s)/(r*r) -1/(r*r) -log(r_s)/(r*r))*z/r ;   
			    p[i].fPot-= M_l*A -Anfw*(-log(r+r_s) + log(r_s))/r;
			    }
			      else {
			    p[i].a[0]-= M_l*A*A*x/(r+ eps1) +Anfw*(r_s/(dsoft+r_s)/(dsoft*dsoft) +log(dsoft+r_s)/(dsoft*dsoft) -1/(dsoft*dsoft) -log(r_s)/(dsoft*dsoft))*x/dsoft ; 
			    p[i].a[1]-= M_l*A*A*y/(r+eps1)  +Anfw*(r_s/(dsoft+r_s)/(dsoft*dsoft) +log(dsoft+r_s)/(dsoft*dsoft) -1/(dsoft*dsoft) -log(r_s)/(dsoft*dsoft))*y/dsoft ;   
			    p[i].a[2]-= M_l*A*A*z/(r+ eps1) +Anfw*(r_s/(dsoft+r_s)/(dsoft*dsoft) +log(dsoft+r_s)/(dsoft*dsoft) -1/(dsoft*dsoft) -log(r_s)/(dsoft*dsoft))*z/dsoft ;   
			    p[i].fPot-= M_l*A -Anfw*(-log(dsoft+r_s) + log(r_s))/dsoft; 
			    }
			  }
		}
	}
}


void pkdHomogSpheroid(PKD pkd)
{
	PARTICLE *p;
	int i,n;

#ifdef SINKING
	const double M_s = 1;	
	const double r_s = 1;
#else
	const double M_s = 5;	
	const double r_s = 10;
#endif
	p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;++i) {
		if (TYPEQueryACTIVE(&(p[i]))) {
			double x = p[i].r[0];
			double y = p[i].r[1];
			double z = p[i].r[2];
			/*
			 **	Do the spheroid potential
			 */
			double r = sqrt(x*x + y*y + z*z);
			double Mr,A;
			Mr = (r < r_s ? pow(r/r_s,3.)*M_s : M_s);
			if (r > 0) {
			    A = Mr/(r*r*r);
			    p[i].a[0] -= A*x;
			    p[i].a[1] -= A*y;
			    p[i].a[2] -= A*z;
			    p[i].fPot += 2*( r < r_s ? 0.5*(Mr/r-3*M_s/r_s) : -M_s/r );
			    }
			}
		}
	}

void pkdGalaxyDiskVerticalPotentialForce(PKD pkd, double Vc, double R)
{
		  /*
			  -  This is the external disk potential that is used together with Chris 
			  -  Gatopolous' Enzo initial conditions for a disk slice.	The initial 
			  -  values Chris used for Vc and R were 220 km/s and 6 kpc respectively.
			  -  */
		 PARTICLE *p;
		int i,n;
		p = pkd->pStore;
		n = pkdLocal(pkd);
		for (i=0;i<n;++i) 
			{
				if (TYPEQueryACTIVE(&(p[i]))) 
					{
						  double z = p[i].r[2];
						  double g = Vc*Vc*z/(R*R+z*z);
						  p[i].a[2] -= g;
						  p[i].fPot += g*z;
					}
				}
}
void pkdBodyForce(PKD pkd, double dConst)
{
	PARTICLE *p;
	int i,n;

	p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;++i) {
		if (TYPEQueryACTIVE(&(p[i]))) {
#ifdef GLASS
		    double ir2= dConst/(p[i].r[0]*p[i].r[0]+p[i].r[1]*p[i].r[1]+p[i].r[2]*p[i].r[2]+p[i].fSoft*p[i].fSoft);
		    p[i].a[0] -= p[i].r[0]*ir2;
		    p[i].a[1] -= p[i].r[1]*ir2;
		    p[i].a[2] -= p[i].r[2]*ir2;
#else
			if (p[i].r[2]>0) {
				p[i].a[2] -= dConst;
				p[i].fPot += dConst*p[i].r[2];
				}
			else {
				p[i].a[2] += dConst;
				p[i].fPot -= dConst*p[i].r[2];
				}
#endif
			}
		}
	}

void pkdMiyamotoDisk(PKD pkd)
{
	PARTICLE *p;
	int i,n;

	const double M_d = 1.0e6;	/* in 10^5 M_sun */
	const double aa = 6.5;		/* in kpc */
	const double b = 0.26;		/* in kpc */
	p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;++i) {
		if (TYPEQueryACTIVE(&(p[i]))) {
			double x = p[i].r[0];
			double y = p[i].r[1];
			double z = p[i].r[2];
			/*
			 **	Do the Miyamoto Disk potential
			 */
			double U = sqrt(z*z + b*b);
			double W = aa + U;
			double A = sqrt(x*x + y*y + W*W);
			p[i].a[0] -= M_d*x/(A*A*A);
			p[i].a[1] -= M_d*y/(A*A*A);
			p[i].a[2] -= M_d*z*W/(A*A*A*U);
			p[i].fPot -= M_d/A;
			}
		}
	}

/*
 ** This is a copy of the Miyamoto function above, but change it as you
 ** like.
 */
void pkdTimeVarying(PKD pkd,double dTime)
{
	const double M_d0 = 1.0e6;	/* in 10^5 M_sun */
	const double aa = 6.5;		/* in kpc */
	const double b = 0.26;		/* in kpc */
	PARTICLE *p;
	int i,n;
	double M_d;
	
	p = pkd->pStore;
	n = pkdLocal(pkd);
	M_d = M_d0 * dTime;
	for (i=0;i<n;++i) {
		if (TYPEQueryACTIVE(&(p[i]))) {
			double x = p[i].r[0];
			double y = p[i].r[1];
			double z = p[i].r[2];
			/*
			 **	Do the Miyamoto Disk potential
			 */
			double U = sqrt(z*z + b*b);
			double W = aa + U;
			double A = sqrt(x*x + y*y + W*W);
			p[i].a[0] -= M_d*x/(A*A*A);
			p[i].a[1] -= M_d*y/(A*A*A);
			p[i].a[2] -= M_d*z*W/(A*A*A*U);
			p[i].fPot -= M_d/A;
			}
		}
	}


/*DEBUG #define SPINUP*/

#ifdef ROT_FRAME
void
pkdRotFrame(PKD pkd,double dOmega,double dOmegaDot)
{
	/* WARNING: p[i].fPot not updated */

#ifdef SPINUP
	PARTICLE *p = pkd->pStore;
	int i;
	for (i=0;i<pkd->nLocal;i++) {
		p[i].a[0] -= dOmegaDot*p[i].r[1];
		p[i].a[1] += dOmegaDot*p[i].r[0];
		}
#else
	/* note Omega & dOmega/dt are assumed to be in the +z direction */

	PARTICLE *p;
	double w2,w2r[2],wxv[2],dwxr[2];
	int i,k;

	p = pkd->pStore;
	w2 = dOmega*dOmega;
	for (i=0;i<pkd->nLocal;i++) {
		if (!TYPEQueryACTIVE(&p[i])) continue;
		w2r[0] = -w2*p[i].r[0];
		w2r[1] = -w2*p[i].r[1];
		wxv[0] = -dOmega*p[i].vPred[1];
		wxv[1] =  dOmega*p[i].vPred[0];
		dwxr[0] = -dOmegaDot*p[i].r[1];
		dwxr[1] =  dOmegaDot*p[i].r[0];
		for (k=0;k<2;k++)
			p[i].a[k] -= (w2r[k] + 2*wxv[k] + dwxr[k]);
		}
#endif
	}
#endif

void pkdCalcEandL(PKD pkd,double *T,double *U,double *Eth,double L[])
{
	/* L is calculated with respect to the origin (0,0,0) */

	PARTICLE *p;
	FLOAT rx,ry,rz,vx,vy,vz;
	int i,n;

#ifdef COLLISIONS
	FLOAT wx,wy,wz,moi;
#endif

	p = pkd->pStore;
	n = pkdLocal(pkd);
	*T = 0.0;
	*U = 0.0;
	*Eth = 0.0;
	L[0] = L[1] = L[2] = 0;
	for (i=0;i<n;++i) {
		rx = p[i].r[0]; ry = p[i].r[1]; rz = p[i].r[2];
		vx = p[i].v[0]; vy = p[i].v[1]; vz = p[i].v[2];
		*T += 0.5*p[i].fMass*(vx*vx + vy*vy + vz*vz);
		*U += 0.5*p[i].fMass*p[i].fPot;
#ifdef GASOLINE
		if (pkdIsGas(pkd,&p[i]))
			*Eth += p[i].fMass*p[i].u
#ifdef MASSNONCOOL
                + p[i].fMassNoncool*(p[i].uNoncool-p[i].u)
#else
#ifdef UNONCOOL
                + p[i].fMass*p[i].uNoncool
#endif
#endif
                ;
#endif
		L[0] += p[i].fMass*(ry*vz - rz*vy);
		L[1] += p[i].fMass*(rz*vx - rx*vz);
		L[2] += p[i].fMass*(rx*vy - ry*vx);
#ifdef COLLISIONS
		wx = p[i].w[0]; wy = p[i].w[1]; wz = p[i].w[2];
		moi = 1.6*p[i].fMass*p[i].fSoft*p[i].fSoft; /* 2/5 MR^2: unf. sphere */
		*T += 0.5*moi*(wx*wx + wy*wy + wz*wz);
		L[0] += moi*wx;
		L[1] += moi*wy;
		L[2] += moi*wz;
#endif
		}
	}


void pkdCalcEandLExt(PKD pkd,double *dMass,double dSumMR[],double dSumMV[],
					 double *dPot)
{
	PARTICLE *p;
	FLOAT m,r2;
	int i,k,n;

	/* Currently this is for the heliocentric reference frame only */

	p = pkd->pStore;
	n = pkdLocal(pkd);
	*dMass = *dPot = 0;
	for (k=0;k<3;k++) dSumMR[k] = dSumMV[k] = 0;
	for (i=0;i<n;i++) {
		*dMass += (m = p[i].fMass);
		r2 = 0;
		for (k=0;k<3;k++) {
			dSumMR[k] += m*p[i].r[k];
			dSumMV[k] += m*p[i].v[k];
			r2 += p[i].r[k]*p[i].r[k];
			}
		if (r2 > 0) *dPot += m/sqrt(r2);
		}
	}

void pkdMassInR(PKD pkd, double R, double *pdMass, FLOAT *com)
{
    PARTICLE *p;
    int i;
    int nLocal;
    double dMass = 0.0;
    FLOAT fCom[3] = {0.0,0.0,0.0};

    p = pkd->pStore;
    nLocal = pkdLocal(pkd);
    for(i = 0; i < nLocal; i++) {
	double r2 = 0.0;
	int k;
	
	for(k = 0; k < 3; k++) {
	    r2 += p[i].r[k]*p[i].r[k];
	    }
	
	if(r2 <= R*R) {
	    dMass += p[i].fMass;
	    for(k = 0; k < 3; k++)
		fCom[k] += p[i].fMass*p[i].r[k];
	    }
	}
    *pdMass = dMass;
    for(i = 0; i < 3; i++)
	com[i] = fCom[i];
    }

/* 
 * Use the f and g functions to advance an unperturbed orbit.
 */
void fg(MDL mdl,double mu,FLOAT *x,FLOAT *v,double dt) {
	double f,g,fd,gd;			/* Gauss's f, g, fdot and gdot */
	double r,vsq;
	double u;					/* r v cos(phi) */
	double a;					/* semi-major axis */
	double e;					/* eccentricity */
	double ec,es;				/* e cos(E), e sin(E) */
	double en;					/* mean motion */
	double nf;					/* orbital frequency */
	double dec;					/* delta E */
	double dm;					/* delta mean anomoly */
	double lo = -2*M_PI;
	double up = 2*M_PI;
	double w;					/* function to zero */
	double wp;					/* first derivative */
	double wpp;					/* second derivative */
	double wppp;				/* third derivative */
	double dx,s,c,sh;
	double next;
	double converge;			/* converge criterion */
	int iter,i,j;
	const double DOUBLE_EPS = 1.2e-16;
	const int MAX_ITER = 256; /*DEBUG originally 32*/

	/* 
	 * Evaluate some orbital quantites.
	 */
	r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
	vsq = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
	u = x[0]*v[0] + x[1]*v[1] + x[2]*v[2];
	a = 1/(2/r-vsq/mu);
	en = sqrt(mu/(a*a*a));
	ec = 1-r/a;
	es = u/(en*a*a);
	e = sqrt(ec*ec + es*es);
#if (0)
	/*
	 ** Only for GR!
	 */
	dt *= 1 - 1.5*mu/(csq*a);
#endif
	nf = en/(2*M_PI);
	j = (int)(nf*dt);
	dt -= j/nf;					/* reduce to single orbit */
	dm = en*dt - es;
	if ((es*cos(dm)+ec*sin(dm)) > 0) dec = dm + 0.85*e;
	else dec = dm - 0.85*e;
	dm = en*dt;					/* reset mean anomoly */
	converge = fabs(dm*DOUBLE_EPS);
	/*
	 * First double precision iterations for dec solution.
	 * This solves Kepler's equation in difference form:
	 * dM = dE - e cos(E_0) sin(dE) + e sin(E_0)(1 - cos(dE))
	 * Note that (1 - cos(dE)) == 2*sin(dE/2)*sin(dE/2).  The latter is
	 * probably more accurate for small dE.
	 */
	for(iter=1;iter<=MAX_ITER;iter++) {
		s = sin(dec);
		c = cos(dec);
		sh = sin(0.5*dec);
		w = dec - ec*s + es*2*sh*sh - dm;
		if (w > 0) up = dec;
		else lo = dec;
		wp = 1 - ec*c + es*s;
		wpp = ec*s + es*c;
		wppp = ec*c - es*s;
		dx = -w/wp;
		dx = -w/(wp + 0.5*dx*wpp);
		dx = -w/(wp + 0.5*dx*wpp + dx*dx*wppp/6);
		next = dec + dx;
		if (fabs(dx) <= converge) break;
		if (next > lo && next < up) dec = next;
		else dec = 0.5*(lo + up);
		if (dec==lo || dec==up) break;
		}
	if (iter>MAX_ITER) {
		char ach[256];

		sprintf(ach,"dec soln failed, dconverge: %g dx: %g dec: %g\n",
				converge, dx, dec);
		mdlDiag(mdl,ach);
		mdlassert(mdl, 0);
		}
	/*
	 * Update the orbit.
	 * See Danby, eq. 6.8.11 and 6.8.12, for expressions for f, g, fdot,
	 * and gdot.
	 * JS: changed f and gd expressions to take advantage of the half angle 
	 * formula given above.
	 */
	f = 1 - (a/r)*2*sh*sh;
	g = dt + (s-dec)/en;
	fd = -(a/(r*wp))*en*s;
	gd = 1 - 2*sh*sh/wp;
	for (i=0;i<3;i++) {
		s = f*x[i]+g*v[i];
		v[i] = fd*x[i]+gd*v[i];
		x[i] = s;
		}
	}

void
pkdDrift(PKD pkd,double dDelta,FLOAT fCenter[3],int bPeriodic,int bInflowOutflow,int bFandG,
		 FLOAT fCentMass, double dTime)
{
	PARTICLE *p;
	int i,j,n;
#ifdef SLIDING_PATCH
	FLOAT fShear;
#endif
	int bInBox = 1;
	FLOAT lfCenter[3];

	mdlDiag(pkd->mdl, "Into pkddrift\n");
	for(j = 0; j < 3; j++)
		lfCenter[j] = fCenter[j];
	n = pkdLocal(pkd);
#ifdef OLD_KEPLER
	if (p->iDriftType == KEPLER) /*DEBUG a bit ugly...*/
#else
	if (bFandG)
#endif
		for (i=0;i<n;++i) {
			p = &pkd->pStore[i];
			fg(pkd->mdl,fCentMass + p->fMass,p->r,p->v,dDelta);
		}
	else {
		p = pkd->pStore;
		for (i=0;i<n;++i,++p) {
#ifdef AGGS
			if (IS_AGG(p))
				continue; /* skip aggregate particles; handled in msrAggsAdvance() */
#endif
#ifdef SLIDING_PATCH
			fShear = 0.0;
			p->bAzWrap = 0; /* reset azimuthal wrap flag */
#endif
#ifdef SINKING
			if (TYPETest( p, TYPE_SINKING)) {
			    FLOAT r0 = p->rSinking0Mag;
			    FLOAT r1 = r0 + p->vSinkingr0*(dTime-p->fSinkingTime);
			    FLOAT r2 = r1 + p->vSinkingr0*dDelta;
			    FLOAT thfac, th1, th2, costh1, sinth1, costh2, sinth2, sqr01, sqr02;
			    if (r2 < 0.1*r0) r2 = 0.1*r0; /* HACK */
			    if (r1 < 0.1*r0) r1 = 0.1*r0; /* HACK */
			    thfac = p->vSinkingTang0Mag*2/(p->vSinkingr0);
			    sqr01 = sqrt(r0/r1);
			    sqr02 = sqrt(r0/r2);
			    th1 = thfac*(1-sqr01);
			    th2 = thfac*(1-sqr02);
			    costh1 = cos(th1);
			    sinth1 = sin(th1);
			    costh2 = cos(th2);
			    sinth2 = sin(th2);
/*#define SINKFREEZE*/
#ifndef SINKFREEZE
			    for (j=0;j<3;j++) {
				p->r[j] += (r2*costh2-r1*costh1)*p->rSinking0Unit[j]+(r2*sinth2-r1*sinth1)*p->vSinkingTang0Unit[j];
/* Do vpred adjustment here -- no need to do it in kickvpred */
				p->vPred[j] += p->vSinkingTang0Mag*(
				    (-sinth2*sqr02+sinth1*sqr01)*p->rSinking0Unit[j]
				    +(costh2*sqr02-costh1*sqr01)*p->vSinkingTang0Unit[j])
				    +p->vSinkingr0*((costh2-costh1)*p->rSinking0Unit[j]
						    +(sinth2-sinth1)*p->vSinkingTang0Unit[j]);
				}
#endif
#ifdef SINKDBG
			    if (p->iOrder == 55) {
				printf("SINKINGDRIFT %d: %g %g %lf %lf, %g %g\n",p->iOrder,th1,th2,r1,r2,dTime,dTime+dDelta);
				}
#endif
			    }
		    
#endif
#ifdef INFLOWOUTFLOW
			if (bInflowOutflow) {
/*			    if (p->r[0] < pkd->dxInflow) TYPESet(p,TYPE_INFLOW);*/
			    if (p->r[0] > pkd->dxOutflow) TYPESet(p,TYPE_OUTFLOW); 
			    }
#endif

#ifdef GASOLINE
#ifdef DRHODT
			if(pkdIsGas(pkd, p)) {
			    p->fDensity *= exp(-p->fDivv_t*dDelta); // Predictor for density
			    if(dDelta > 0.0)
			       assert(p->fDensity > 0.0);
			    }
			p->fDensity_t *= exp(-p->fDivv_t*dDelta);
			p->fDensity_PdV *= exp(-p->fDivv_PdV*dDelta);
			p->fDensity_PdVcorr *= exp(-p->fDivv_PdVcorr*dDelta);
#endif
#endif
			for (j=0;j<3;++j) {
			        p->r[j] += dDelta*p->v[j];
				if (bPeriodic) {
					if (p->r[j] >= lfCenter[j] + 0.5*pkd->fPeriod[j]) {
						p->r[j] -= pkd->fPeriod[j];
#ifdef SLIDING_PATCH
						if (j == 0) { /* radial wrap */
						  fShear = 1.5*pkd->PP->dOrbFreq*pkd->PP->dWidth; /* (dWidth is same as fPeriod[0]) */
						  p->r[1] += SHEAR(-1,pkd->dTime + dDelta,pkd->PP);
						}
						/* apply randomization only if there was no radial wrap (fShear=0) */
						if (pkd->PP->bRandAzWrap == 1 && j == 1 && fShear == 0.0) /* azimuthal wrap */
						  p->bAzWrap = 1;
#endif
						}
					if (p->r[j] < lfCenter[j] - 0.5*pkd->fPeriod[j]) {
						p->r[j] += pkd->fPeriod[j];
#ifdef SLIDING_PATCH
						if (j == 0) {
						  fShear = - 1.5*pkd->PP->dOrbFreq*pkd->PP->dWidth;
						  p->r[1] += SHEAR(1,pkd->dTime + dDelta,pkd->PP);
						}
						/* apply randomization only if there was no radial wrap (fShear=0) */
						if (pkd->PP->bRandAzWrap == 1 && j == 1 && fShear == 0.0)
						  p->bAzWrap = 1;
#endif
						}
					bInBox = bInBox && (p->r[j] >= lfCenter[j]-0.5*pkd->fPeriod[j]);
					bInBox = bInBox && (p->r[j] <  lfCenter[j]+0.5*pkd->fPeriod[j]);
					assert(bInBox);
					}
			    }
#ifdef INFLOWOUTFLOW
			if (bInflowOutflow) {
			    if (p->r[0] > pkd->dxInflow && TYPETest(p,TYPE_INFLOW)) {  
#if (0)
                                /* duplicate -- assumes constant inflow properties: rho, v, T, etc... */ 
				PARTICLE pInflow;
				pInflow = *p;
				pInflow.r[0] -= (pkd->dxInflow+pkd->fPeriod[0]*.5); /* create new Inflow particle */
				pkdNewParticle(pkd, pInflow);
#endif
				TYPEReset(p, TYPE_INFLOW); /* Original particle joins regular particles now */
				}
			    if (p->r[0] < pkd->dxOutflow && TYPETest(p,TYPE_OUTFLOW)) { /* delete */ 
				pkdDeleteParticle(pkd, p);
				}
			    }
#endif

#ifdef SINKING
#ifdef SINKDBG
			if (p->iOrder == 55) {
			    printf("SINKDRIFT %g, %g %g %g, %g %g %g: %d %g\n",dTime+dDelta,p->v[0],p->v[1],p->v[2],p->r[0],p->r[1],p->r[2],p->iOrder,p->fMass);
			    }
#endif
#endif
#ifdef SLIDING_PATCH
			p->v[1] += fShear;
			p->dPy -= fShear/3.0; /* Angular momentum is
						 also changed. */
#endif
			mdlassert(pkd->mdl, bInBox);
			}
		}
	mdlDiag(pkd->mdl, "Out of pkddrift\n");
	}

void pkdKick(PKD pkd, double dvFacOne, double dvFacTwo, double dvPredFacOne,
	     double dvPredFacTwo, double duDelta, double duPredDelta, int iGasModel,
    double z, double duDotLimit, double dTimeEnd, UNCC uncc )
{
	PARTICLE *p;
	int i,j,n;
#ifdef GLASS 
#define dvFacOneSTD     1.0
#define dvPredFacOneSTD 1.0
#else
#define dvFacOneSTD     dvFacOne
#define dvPredFacOneSTD dvPredFacOne
#endif
	pkdClearTimer(pkd,1);
	pkdStartTimer(pkd,1);

	p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;++i,++p) {
#ifdef AGGS
		if (IS_AGG(p)) /* skip aggregate particles; handled in msrAggsKick() */
			continue;
#endif
		if (TYPEQueryACTIVE(p)) {
#ifdef GASOLINE
			if (pkdIsGas(pkd, p)) {
#ifdef SPH1D
                p->a[1]=0; p->a[2]=0; 
#endif                      
#ifdef GLASSZ
                p->a[0]=0; p->a[1]=0;
#endif
#ifdef ACCZERO
                p->a[0]=0; p->a[1]=0; p->a[2]=0;
#endif
				for (j=0;j<3;++j) {
					p->vPred[j] = p->v[j]*dvPredFacOne + p->a[j]*dvPredFacTwo;
					p->v[j] = p->v[j]*dvFacOne + p->a[j]*dvFacTwo;
				    }
#ifdef VARALPHA
				    {
				    double dalphadt;
				    dalphadt = - p->alphaPred*0.1*p->c/sqrt(0.25*p->fBall2);
#ifdef DODVDS
				    if (p->dvds < 0) dalphadt -= p->dvds*1.5;
#else
				    if (p->divv < 0) dalphadt -= p->divv;
#endif
				    p->alphaPred = p->alpha + dalphadt*duPredDelta;
				    if (p->alphaPred < ALPHAMIN) p->alphaPred = ALPHAMIN;
				    p->alpha = p->alpha + dalphadt*duDelta;
				    if (p->alpha < ALPHAMIN) p->alpha = ALPHAMIN;
				    }
#endif
#ifdef SINKING
				if (TYPETest( p, TYPE_SINKING)) {
				    FLOAT r0 = p->rSinking0Mag;
				    /* For kick std dTime is midpoint of kick period, put dTime to end for this exact calc */
				    FLOAT r2 = r0 + p->vSinkingr0*(dTimeEnd-p->fSinkingTime);
				    FLOAT thfac, sqr02, th2, costh2, sinth2;
				    if (r2 < 0.1*r0) r2 = 0.1*r0; /* HACK */
				    thfac = p->vSinkingTang0Mag*2/(p->vSinkingr0);
				    sqr02 = sqrt(r0/r2);
				    th2 = thfac*(1-sqr02);
				    costh2 = cos(th2);
				    sinth2 = sin(th2);
				    /* v does not include motion around sink -- add it back */
				    for (j=0;j<3;j++) {
					p->vPred[j] += p->vSinkingTang0Mag*sqr02*(-sinth2*p->rSinking0Unit[j]+costh2*p->vSinkingTang0Unit[j])+p->vSinkingr0*(costh2*p->rSinking0Unit[j]+sinth2*p->vSinkingTang0Unit[j]);
					}
				    }
#ifdef SINKDBG
				if (p->iOrder == 55) printf("SINKINGKICK %d %g, %g %g  %g %g  %g %g\n",p->iOrder,dTimeEnd,p->vPred[0],p->v[0],p->vPred[1],p->v[1],p->vPred[2],p->v[2]);
#endif
#endif /* SINKING */
				if (iGasModel != GASMODEL_ISOTHERMAL && iGasModel != GASMODEL_GLASS) {
#ifndef NOCOOLING				
                    p->uPred = p->u + p->uDot*duPredDelta;
                    p->u = p->u + p->uDot*duDelta;
                    if (p->u < 0) {
                        FLOAT uold = p->u - p->uDot*duDelta;
#ifdef FBPARTICLE
                        fprintf(stderr,"FBP Negative! %d: %g %g %g %g\n",p->iOrder,uold,p->u,p->uDot,duDelta);
#endif
                        p->uPred = uold*exp(p->uDot*duPredDelta/uold);
                        p->u = uold*exp(p->uDot*duDelta/uold);
                        }
#ifdef UNONCOOL
                    p->uNoncoolPred = p->uNoncool + p->uNoncoolDot*duPredDelta;
                    p->uNoncool = p->uNoncool + p->uNoncoolDot*duDelta;
                    if (p->uNoncoolPred < 0) p->uNoncoolPred = 0;
                    if (p->uNoncool < 0) p->uNoncool = 0;
#ifdef MASSNONCOOL
					FLOAT upnc52, up52, fMassFlux;
					FLOAT ph = sqrt(0.25*p->fBall2);
				   upnc52 = pow(p->uNoncoolPred, 2.5);
				   up52 = pow(p->uPred, 2.5);
				   FLOAT fFactor = duPredDelta*uncc.gpc.dEvapCoeffCode*ph*ph*3.1415;
				   fMassFlux = fFactor*(upnc52-up52);
				   dbgprint("EVAPINTERNAL: %d %e %e %e %e %e %e %e %e\n", 
						   p->iOrder, duDelta, duPredDelta, fMassFlux, ph, p->fMass-p->fMassNoncool, p->fMassNoncool, p->uPred, p->uNoncoolPred);
				   if(fMassFlux > 0) { // Make sure that the flow is in the right direction
					   // If all the mass becomes hot, switch to being single-phase
					   if(fMassFlux > (p->fMass-p->fMassNoncool)) {
						   p->uPred = (p->uPred*p->fMass + p->uNoncoolPred*p->fMassNoncool)/(p->fMass+p->fMassNoncool);
						   p->u = (p->u*p->fMass + p->uNoncool*p->fMassNoncool)/(p->fMass+p->fMassNoncool);
						   p->uDot = (p->uDot*p->fMass + p->uNoncoolDot*p->fMassNoncool)/(p->fMass+p->fMassNoncool);
						   p->fMassNoncool = 0;
						   p->uNoncool = 0;
						   p->uNoncoolDot = 0;
						   p->uNoncoolPred = 0;
					   }
					   else {
						   p->uNoncoolPred = (p->uPred*fMassFlux + p->uNoncoolPred*p->fMassNoncool)/(fMassFlux+p->fMassNoncool);
						   p->uNoncool = (p->u*fMassFlux + p->uNoncool*p->fMassNoncool)/(fMassFlux+p->fMassNoncool);
						   p->fMassNoncool += fMassFlux;
						   assert(p->fMassNoncool >= 0);
						   assert(p->uPred > 0);
						   assert(p->uNoncoolPred > 0);
						   assert(p->u > 0);
						   assert(p->uNoncool > 0);
					   }
				   }
				   else if (p->uPred > p->uNoncoolPred) { // No sense in keeping the noncooling mass around if it is much colder than the regular mass
						   p->uPred = (p->uPred*p->fMass + p->uNoncoolPred*p->fMassNoncool)/(p->fMass+p->fMassNoncool);
						   p->u = (p->u*p->fMass + p->uNoncool*p->fMassNoncool)/(p->fMass+p->fMassNoncool);
						   p->uDot = (p->uDot*p->fMass + p->uNoncoolDot*p->fMassNoncool)/(p->fMass+p->fMassNoncool);
						   p->fMassNoncool = 0;
						   p->uNoncool = 0;
						   p->uNoncoolDot = 0;
						   p->uNoncoolPred = 0;
				   }
#endif
#endif
#else /* NOCOOLING */
                    p->uPred = p->u + UDOT_HYDRO(p)*duPredDelta;
                    p->u = p->u + UDOT_HYDRO(p)*duDelta;
#endif /* NOCOOLING */
#if defined(PRES_HK) || defined(PRES_MONAGHAN) 
                    if (p->uPred < 0) p->uPred = 0;
                    if (p->u < 0) p->u = 0;
#endif /* PRES_HK  PRES_MONAGHAN */
                    }
#ifdef DIFFUSION
				p->fMetalsPred = p->fMetals + p->fMetalsDot*duPredDelta;
				p->fMetals = p->fMetals + p->fMetalsDot*duDelta;
#ifdef MASSDIFF
				p->fMass = p->fMass0 + p->fMassDot*duPredDelta;
				p->fMass0 = p->fMass0 + p->fMassDot*duDelta;
#endif
#ifdef STARFORM
				p->fMFracOxygenPred = p->fMFracOxygen + p->fMFracOxygenDot*duPredDelta;
				p->fMFracOxygen = p->fMFracOxygen + p->fMFracOxygenDot*duDelta;
				p->fMFracIronPred = p->fMFracIron + p->fMFracIronDot*duPredDelta;
				p->fMFracIron = p->fMFracIron + p->fMFracIronDot*duDelta;
#endif /* STARFORM */
#endif /* DIFFUSION */
			    }
			else 
#endif /* GASOLINE */
			    { /* Not gas or not -DGASOLINE */
			    for (j=0;j<3;++j) {
#if defined(NEED_VPRED) && !defined(GASOLINE)
				p->vPred[j] = p->v[j]*dvPredFacOneSTD + p->a[j]*dvPredFacTwo;
#endif /* NEED_VPRED */
				p->v[j] = p->v[j]*dvFacOneSTD + p->a[j]*dvFacTwo;
			        }
			    }
		    }
	    }

	pkdStopTimer(pkd,1);
	mdlDiag(pkd->mdl, "Done pkdkick\n");
	}

void pkdEmergencyAdjust(PKD pkd, int iRung, int iMaxRung, double dDelta, double dDeltaThresh, int *pnUn, int *piMaxRungIdeal, int *pnMaxRung, int *piMaxRungOut)
    {
#ifdef GASOLINE
	PARTICLE *p;
	int i,j,n;
    int iTempRung;
    int iMaxRungOut=0;
    int iMaxRungIdeal=0;
    int nMaxRung=0;
    int nUn=0;
 
    assert(dDeltaThresh < dDelta);
	p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;++i,++p) {
		if (pkdIsGas(pkd, p)) { // Only test non-active gas particles 
            if (p->iRung < iRung && p->dtNew < dDeltaThresh) {
                /* Emergency Adjust */
                nUn++;

                p->dt = p->dtNew;
				iTempRung = pkdOneParticleDtToRung( iRung,dDelta,p->dt );
                assert(iTempRung > iRung);
                
				if(iTempRung >= iMaxRungIdeal)
					iMaxRungIdeal = iTempRung+1;
				if(iTempRung >= iMaxRung) {
					iTempRung = iMaxRung-1;
                    }
				p->iRung = iTempRung;
                
                if(p->iRung > iMaxRungOut) {
                    iMaxRungOut = p->iRung;
                    nMaxRung = 1;
                    }
                else if (p->iRung == iMaxRungOut) 
                    nMaxRung ++;
                
                /* UnKick -- revert to predicted values -- low order, non symplectic :( */
				for (j=0;j<3;++j) {
					p->v[j] = p->vPred[j];
				    }
#ifndef NOCOOLING				
                p->u = p->uPred;
#ifdef UNONCOOL
                p->uNoncool = p->uNoncoolPred;
#endif
#endif
#ifdef DIFFUSION
				p->fMetals = p->fMetalsPred;
#ifdef STARFORM
				p->fMFracOxygen = p->fMFracOxygenPred;
				p->fMFracIron = p->fMFracIronPred;
#endif /* STARFORM */
#endif /* DIFFUSION */
                }
            }
        }
    *pnUn = nUn;
    *piMaxRungIdeal = iMaxRungIdeal;
    *pnMaxRung = nMaxRung;
    *piMaxRungOut = iMaxRungOut;
#endif /* GASOLINE */
    }

void pkdKickPatch(PKD pkd, double dvFacOne, double dvFacTwo,
		  double dOrbFreq, /* Orbital Frequency of Patch */
		  int bOpen)	/* Is this an opening Kick? */

{
#ifdef SLIDING_PATCH
	PARTICLE *p;
	int i,j,n;

	pkdClearTimer(pkd,1);
	pkdStartTimer(pkd,1);

	p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;++i,++p) {
#ifdef AGGS
		if (IS_AGG(p)) /* skip aggregate particles; handled in msrAggsKick() */
			continue;
#endif
		if (TYPEQueryACTIVE(p)) {
		    if(!bOpen) {
			/* perform Cross Hamiltonian */
			p->v[0] += 2.0*dvFacTwo*dOrbFreq*p->dPy;
			p->v[1] = p->dPy - 2*dOrbFreq*p->r[0];
			}
		    for (j=0;j<3;++j) {
			p->v[j] = p->v[j]*dvFacOne + p->a[j]*dvFacTwo;
			}
		    if(bOpen) {
			p->dPy = p->v[1] + 2.0*dOrbFreq*p->r[0];
			/* perform Cross Hamiltonian */
			/* note that the y velocity is being set up
			 * so that the subsequent "Drift" will advance
			 * the positions all the way to the end of the
			 * closing cross Hamiltonian.
			 */
			p->v[0] += 2.0*dvFacTwo*dOrbFreq*p->dPy;
			/* The terms are: the drift (py), the cross
			 * term at the beginning of the drift, and the
			 * cross term at the end of the drift.
			 */
			p->v[1] = p->dPy - dOrbFreq*p->r[0]
			    - dOrbFreq*(p->r[0] + 2.0*dvFacTwo*p->v[0]);
			}
		    }
	    }

	pkdStopTimer(pkd,1);
	mdlDiag(pkd->mdl, "Done pkdkick\n");
#endif
	}

void pkdGravInflow(PKD pkd, double r) {

	PARTICLE *p;
	int i,n;

	p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;++i) {
		if (TYPEQueryACTIVE(&(p[i]))) {
		    p[i].a[0] -= 0.1; /* Hardwired for test case */
		    p[i].fPot += 0.1*p[i].r[0];
		}
	    }
    }

void pkdCreateInflow(PKD pkd, int Ny, int iGasModel, double dTuFac, double pmass, double x, double vx, double density, double temp, double metals, double eps, double dt, int iRung)
    {

    PARTICLE p;
    int j,k;
    int start,end;
    double yy,zz,dx;
    
    dx = pkd->fPeriod[1]/Ny;

    start = pkd->idSelf*Ny/pkd->nThreads;
    end = (pkd->idSelf+1)*Ny/pkd->nThreads;

    for (j=start;j<end;j++) {
	yy = pkd->fPeriod[1]*((1.0*(j+0.5))/Ny-0.5);
	for (k=0;k<Ny;k++) {
	    zz = pkd->fPeriod[2]*((1.0*(k+(j&1)*0.5+0.25))/Ny-0.5);
	    TYPEClear(&p);
	    TYPESet(&p,TYPE_GAS);
	    TYPESet(&p,TYPE_INFLOW);
	    p.fMass = pmass;
	    p.r[0] = x;
	    p.r[1] = yy;
	    p.r[2] = zz;
	    p.v[0] = vx;
	    p.v[1] = 0;
	    p.v[2] = 0;
#ifdef NEED_VPRED
	    p.vPred[0] = vx;
	    p.vPred[1] = 0;
	    p.vPred[2] = 0;
#endif
	    p.fSoft = eps;
#ifdef CHANGESOFT				
	    p.fSoft0 = eps;
#endif
	    p.fDensity = density;

	    p.iRung = iRung;
	    p.dt = dt;
	    p.fWeight = 1.0;
	    p.fBall2 = 4*dx*dx;
	    p.fBallMax = p.fBall2;
#ifdef GASOLINE
#ifndef NOCOOLING
	    if (iGasModel == GASMODEL_COOLING) { 
#ifdef COOLDEBUG
	        pkd->Cool->p = &p; /* Send in particle pointer only for temporary debug */
#endif
		CoolInitEnergyAndParticleData( pkd->Cool, &p.CoolParticle, &p.u, density, temp, p.fMetals );
		p.uPred = p.u;
		}
	    else
#endif
		{
		p.u = dTuFac*temp;
		p.uPred = dTuFac*temp;
		}
	    p.c = sqrt(5./3.)*(5./3.-1)*p.uPred; /* Hack */
	    p.uDotPdV = 0;
	    p.uDotAV = 0;
	    p.uDotDiff = 0;
#ifndef NOCOOLING
	    p.uDot = 0;
#endif
	    p.fMetals = metals;
#ifdef DIFFUSION
	    p.fMetalsPred = metals;
#endif				
#ifdef STARFORM
	    p.uDotFB = 0.0;
	    p.fNSN = 0.0;
	    p.fNSNtot = 0.0;
	    p.fMOxygenOut = 0.0;
	    p.fMIronOut = 0.0;
	    p.fMFracOxygen = 0.0;
	    p.fMFracIron = 0.0;
	    p.fTimeCoolIsOffUntil = 0.0;
#endif
#ifdef SIMPLESF
	    p.fMassStar = 0;
#endif
#ifdef SHOCKTRACK
	    p.ShockTracker = 0.0;
#endif
#ifdef VARALPHA
	    p.alpha = ALPHAMIN;
	    p.alphaPred = ALPHAMIN;
	    p.divv = 0;
#ifdef DODVDS
	    p.dvds = 0;
#endif
#endif
	    p.c = 0.0;
	    p.fTimeForm = 0.0;
#endif


	    pkdNewParticle(pkd, p);
	    }
	}
    }

void pkdReadCheck(PKD pkd,char *pszFileName,int iVersion,int iOffset,
				  int nStart,int nLocal)
{
	FILE *fp;
	CHKPART cp;
	long lStart;
	int i,j;

	pkd->nLocal = nLocal;
	/*
	 ** Seek past the header and up to nStart.
	 */
	fp = fopen(pszFileName,"r");
	mdlassert(pkd->mdl,fp != NULL);
	lStart = iOffset+nStart*sizeof(CHKPART);
	fseek(fp,lStart,SEEK_SET);
	/*
	 ** Read Stuff!
	 */
	for (i=0;i<nLocal;++i) {
		PARTICLE *p = &pkd->pStore[i];
	    
		fread(&cp,sizeof(CHKPART),1,fp);
		p->iOrder = cp.iOrder;
		p->iActive = cp.iActive;
		p->fMass = cp.fMass;
		assert(p->fMass >= 0.0);
#ifdef CHANGESOFT
		p->fSoft0 = cp.fSoft;
#endif
		p->fSoft = cp.fSoft;
        p->dt = FLT_MAX;
        p->dtNew = FLT_MAX;
		for (j=0;j<3;++j) {
			p->r[j] = cp.r[j];
			p->v[j] = cp.v[j];
#ifdef NEED_VPRED
			p->vPred[j] = cp.v[j];
#endif
			}
#ifdef INFLOWOUTFLOW
		if (p->r[0] < pkd->dxInflow) TYPESet(p,TYPE_INFLOW);
		if (p->r[0] > pkd->dxOutflow) TYPESet(p,TYPE_OUTFLOW);
#endif
#ifdef GASOLINE
		p->u = cp.u;
		p->uPred = cp.u;
#ifdef MASSNONCOOL
        p->fMassNoncool = cp.fMassNoncool;
#endif
#ifdef UNONCOOL
#ifdef UNONCOOLMERGE
        p->u += cp.uNoncool;
		p->uPred += cp.uNoncool;
        p->uNoncool = 0;
#else
		p->uNoncool = cp.uNoncool;
#endif
		p->uNoncoolPred = p->uNoncool;
		assert(p->uNoncool >= 0);
        p->uNoncoolDot = 0;
#endif
#ifdef STARSINK
		SINK_Lx(p) = cp.Lx;
		SINK_Ly(p) = cp.Ly;
		SINK_Lz(p) = cp.Lz;
#endif

#ifdef VARALPHA
		p->alpha = cp.alpha;
		p->alphaPred = cp.alpha;
#endif
#ifdef COOLDEBUG
		if (p->iOrder == 842079) fprintf(stderr,"Particle %i in pStore[%i]\n",p->iOrder,(int) (p-pkd->pStore));
#endif
#ifdef STARFORM 
		p->uDotFB = 0.0;
		p->fTimeForm = cp.fTimeForm;
		p->fMassForm = cp.fMassForm;
		p->iGasOrder = cp.iGasOrder;
        p->fTimeCoolIsOffUntil = cp.fTimeCoolIsOffUntil;
        p->fNSN = 0.0;
        p->fNSNtot = 0.0;
        p->fMFracOxygen = cp.fMFracOxygen;
        p->fMFracIron = cp.fMFracIron;
#ifdef DIFFUSION
        p->fMFracOxygenPred = cp.fMFracOxygen;
        p->fMFracIronPred = cp.fMFracIron;
#endif
#endif
#ifdef SIMPLESF
		p->fMassStar = cp.fMassStar;
		p->fTimeForm = cp.fTimeForm;
		for (j=0;j<3;++j) {
			p->rForm[j] = cp.rForm[j];
			p->vForm[j] = cp.vForm[j];
			}
		p->fDensity = cp.fDenForm;
		p->iGasOrder = cp.iGasOrder;
#endif
		p->fMetals = cp.fMetals;
#ifdef DIFFUSION
		p->fMetalsPred = cp.fMetals;

#ifdef MASSDIFF
		p->fMass0 = cp.fMass;
#endif
#endif
#ifndef NOCOOLING		
		p->CoolParticle = cp.CoolParticle;
#endif
#ifdef SINKING
		p->rSinking0Unit[0] = cp.rSinking0Unit[0];
		p->rSinking0Unit[1] = cp.rSinking0Unit[1];
		p->rSinking0Unit[2] = cp.rSinking0Unit[2];
		p->rSinking0Mag = cp.rSinking0Mag;
		p->vSinkingTang0Unit[0] = cp.vSinkingTang0Unit[0];
		p->vSinkingTang0Unit[1] = cp.vSinkingTang0Unit[1];
		p->vSinkingTang0Unit[2] = cp.vSinkingTang0Unit[2];
		p->vSinkingTang0Mag = cp.vSinkingTang0Mag;
		p->vSinkingr0 = cp.vSinkingr0;
		p->fSinkingTime = cp.fSinkingTime;  
		p->iSinkingOnto = cp.iSinkingOnto;
		p->fTrueMass = cp.fTrueMass;
#endif
#endif
#ifdef COLLISIONS
		p->iOrgIdx = cp.iOrgIdx;
		for (j=0;j<3;++j)
			p->w[j] = cp.w[j];
		p->iColor = cp.iColor;
#endif
//		TYPEClear(p);
        TYPEReset(p,TYPE_ACTIVE|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE|TYPE_SMOOTHDONE|TYPE_DensACTIVE|TYPE_NbrOfACTIVE|TYPE_Scatter|TYPE_DensZeroed);
		TYPESet(p,TYPE_ACTIVE);
		p->iRung = 0;
		p->fWeight = 1.0;
		p->fDensity = 0.0;
		p->fBall2 = 0.0;
		p->fBallMax = 0.0;
		}
	fclose(fp);	
	}


void pkdWriteCheck(PKD pkd,char *pszFileName,int iOffset,int nStart)
{
	FILE *fp;
	CHKPART cp;
	PARTICLE *p;
	long lStart;
	int i,j,nLocal;
	int nout;

	/*
	 ** Seek past the header and up to nStart.
	 */
	fp = fopen(pszFileName,"r+");
	mdlassert(pkd->mdl,fp != NULL);
	lStart = iOffset+nStart*sizeof(CHKPART);
	nout = fseek(fp,lStart,0);
	mdlassert(pkd->mdl,nout == 0);
	/* 
	 ** Write Stuff!
	 */
	nLocal = pkdLocal(pkd);
	for (i=0;i<nLocal;++i) {
        p = &pkd->pStore[i];
		cp.iOrder = p->iOrder;
		cp.iActive = p->iActive;
		cp.fMass = p->fMass;
#ifdef CHANGESOFT
		cp.fSoft = p->fSoft0;
#else
		cp.fSoft = p->fSoft;
#endif
		for (j=0;j<3;++j) {
			cp.r[j] = p->r[j];
			cp.v[j] = p->v[j];
			}
#ifdef GASOLINE
		cp.u = p->u;
#ifdef MASSNONCOOL
		cp.fMassNoncool = p->fMassNoncool;
#endif
#ifdef UNONCOOL
		cp.uNoncool = p->uNoncool;
#endif
#ifdef STARSINK
		cp.Lx = SINK_Lx(p);
		cp.Ly = SINK_Ly(p);
		cp.Lz = SINK_Lz(p);
#endif
#ifdef VARALPHA
		cp.alpha = p->alpha;
#endif
		cp.fMetals = p->fMetals;
#ifndef NOCOOLING		
		cp.CoolParticle = p->CoolParticle;
#endif
#ifdef STARFORM
		cp.fTimeForm = p->fTimeForm;
		cp.fMassForm = p->fMassForm;
		cp.iGasOrder = p->iGasOrder;
        cp.fTimeCoolIsOffUntil = p->fTimeCoolIsOffUntil;
        cp.fMFracOxygen = p->fMFracOxygen;
        cp.fMFracIron = p->fMFracIron;
#endif
#ifdef SIMPLESF
		cp.fMassStar = p->fMassStar;
		cp.fTimeForm = p->fTimeForm;
		for (j=0;j<3;++j) {
			cp.rForm[j] = p->rForm[j];
			cp.vForm[j] = p->vForm[j];
			}
		cp.fDenForm = p->fDensity;
		cp.iGasOrder = p->iGasOrder;
#endif
#ifdef SINKING
		cp.rSinking0Unit[0] = p->rSinking0Unit[0];
		cp.rSinking0Unit[1] = p->rSinking0Unit[1];
		cp.rSinking0Unit[2] = p->rSinking0Unit[2];
		cp.rSinking0Mag = p->rSinking0Mag;
		cp.vSinkingTang0Unit[0] = p->vSinkingTang0Unit[0];
		cp.vSinkingTang0Unit[1] = p->vSinkingTang0Unit[1];
		cp.vSinkingTang0Unit[2] = p->vSinkingTang0Unit[2];
		cp.vSinkingTang0Mag = p->vSinkingTang0Mag;
		cp.vSinkingr0 = p->vSinkingr0;
		cp.fSinkingTime = p->fSinkingTime;  
		cp.iSinkingOnto = p->iSinkingOnto;
		cp.fTrueMass = p->fTrueMass;
#endif
#endif
#ifdef COLLISIONS
		cp.iOrgIdx = p->iOrgIdx;
		for (j=0;j<3;++j)
			cp.w[j] = p->w[j];
		cp.iColor = p->iColor;
#endif /* COLLISIONS */
		nout = fwrite(&cp,sizeof(CHKPART),1,fp);
		mdlassert(pkd->mdl,nout == 1);
		}
	nout = fclose(fp);
	mdlassert(pkd->mdl,nout == 0);
	}


void pkdDistribCells(PKD pkd,int nCell,KDN *pkdn)
{
	int i;

	if (pkd->kdTop != NULL) free(pkd->kdTop);
	if (pkd->piLeaf != NULL) free(pkd->piLeaf);
	pkd->kdTop = malloc(nCell*sizeof(KDN));
	mdlassert(pkd->mdl,pkd->kdTop != NULL);
	pkd->piLeaf = malloc(pkd->nThreads*sizeof(int));
	mdlassert(pkd->mdl,pkd->piLeaf != NULL);
	for (i=1;i<nCell;++i) {
		if (pkdn[i].pUpper) {
			pkd->kdTop[i] = pkdn[i];
			if (pkdn[i].pLower >= 0) pkd->piLeaf[pkdn[i].pLower] = i;
			}
		}
	}


void pkdCalcRoot(PKD pkd,struct ilCellNewt *pcc)
{
	KDN *pkdn;
	int pj;
	double m,dx,dy,dz;
	double d2;

	/*
	 ** Initialize moments.
	 */
	pcc->xxxx = 0.0;
	pcc->xyyy = 0.0;
	pcc->xxxy = 0.0;
	pcc->yyyy = 0.0;
	pcc->xxxz = 0.0;
	pcc->yyyz = 0.0;
	pcc->xxyy = 0.0;
	pcc->xxyz = 0.0;
	pcc->xyyz = 0.0;
	pcc->xxzz = 0.0;
	pcc->xyzz = 0.0;
	pcc->xzzz = 0.0;
	pcc->yyzz = 0.0;
	pcc->yzzz = 0.0;
	pcc->zzzz = 0.0;
	pcc->xxx = 0.0;
	pcc->xyy = 0.0;
	pcc->xxy = 0.0;
	pcc->yyy = 0.0;
	pcc->xxz = 0.0;
	pcc->yyz = 0.0;
	pcc->xyz = 0.0;
	pcc->xzz = 0.0;
	pcc->yzz = 0.0;
	pcc->zzz = 0.0;
	/*
	 ** Calculate moments.
	 */
	pkdn = &pkd->kdNodes[pkd->iRoot];
	for (pj=pkdn->pLower;pj<=pkdn->pUpper;++pj) {
		m = pkd->pStore[pj].fMass;
		dx = pkd->pStore[pj].r[0] - pkdn->r[0];
		dy = pkd->pStore[pj].r[1] - pkdn->r[1];
		dz = pkd->pStore[pj].r[2] - pkdn->r[2];
#ifdef REDUCED_EWALD
		d2 = dx*dx + dy*dy + dz*dz;
#else
		d2 = 0.0;
#endif
		pcc->xxxx += m*(dx*dx*dx*dx - 6.0/7.0*d2*(dx*dx - 0.1*d2));
		pcc->xyyy += m*(dx*dy*dy*dy - 3.0/7.0*d2*dx*dy);
		pcc->xxxy += m*(dx*dx*dx*dy - 3.0/7.0*d2*dx*dy);
		pcc->yyyy += m*(dy*dy*dy*dy - 6.0/7.0*d2*(dy*dy - 0.1*d2));
		pcc->xxxz += m*(dx*dx*dx*dz - 3.0/7.0*d2*dx*dz);
		pcc->yyyz += m*(dy*dy*dy*dz - 3.0/7.0*d2*dy*dz);
		pcc->xxyy += m*(dx*dx*dy*dy - 1.0/7.0*d2*(dx*dx + dy*dy - 0.2*d2));
		pcc->xxyz += m*(dx*dx*dy*dz - 1.0/7.0*d2*dy*dz);
		pcc->xyyz += m*(dx*dy*dy*dz - 1.0/7.0*d2*dx*dz);
		pcc->xxzz += m*(dx*dx*dz*dz - 1.0/7.0*d2*(dx*dx + dz*dz - 0.2*d2));
		pcc->xyzz += m*(dx*dy*dz*dz - 1.0/7.0*d2*dx*dy);
		pcc->xzzz += m*(dx*dz*dz*dz - 3.0/7.0*d2*dx*dz);
		pcc->yyzz += m*(dy*dy*dz*dz - 1.0/7.0*d2*(dy*dy + dz*dz - 0.2*d2));
		pcc->yzzz += m*(dy*dz*dz*dz - 3.0/7.0*d2*dy*dz);
		pcc->zzzz += m*(dz*dz*dz*dz - 6.0/7.0*d2*(dz*dz - 0.1*d2));
		/*
		 ** Calculate octopole moment...
		 */
		pcc->xxx += m*(dx*dx*dx - 0.6*d2*dx);
		pcc->xyy += m*(dx*dy*dy - 0.2*d2*dx);
		pcc->xxy += m*(dx*dx*dy - 0.2*d2*dy);
		pcc->yyy += m*(dy*dy*dy - 0.6*d2*dy);
		pcc->xxz += m*(dx*dx*dz - 0.2*d2*dz);
		pcc->yyz += m*(dy*dy*dz - 0.2*d2*dz);
		pcc->xyz += m*dx*dy*dz;
		pcc->xzz += m*(dx*dz*dz - 0.2*d2*dx);
		pcc->yzz += m*(dy*dz*dz - 0.2*d2*dy);
		pcc->zzz += m*(dz*dz*dz - 0.6*d2*dz);
		}
	}


void pkdDistribRoot(PKD pkd,struct ilCellNewt *pcc)
{
	KDN *pkdn;
	double tr;

	pkd->ilcnRoot = *pcc;
	/*
	 ** Must set the quadrupole, mass and cm.
	 */
	pkdn = &pkd->kdTop[ROOT];
#ifdef REDUCED_EWALD
	tr = pkdn->mom.Qxx + pkdn->mom.Qyy + pkdn->mom.Qzz;
#else
	tr = 0.0;
#endif
	pkd->ilcnRoot.m = pkdn->fMass;
	pkd->ilcnRoot.x = pkdn->r[0];
	pkd->ilcnRoot.y = pkdn->r[1];
	pkd->ilcnRoot.z = pkdn->r[2];
	pkd->ilcnRoot.xx = pkdn->mom.Qxx - tr/3.0;
	pkd->ilcnRoot.xy = pkdn->mom.Qxy;
	pkd->ilcnRoot.xz = pkdn->mom.Qxz;
	pkd->ilcnRoot.yy = pkdn->mom.Qyy - tr/3.0;
	pkd->ilcnRoot.yz = pkdn->mom.Qyz;
	pkd->ilcnRoot.zz = pkdn->mom.Qzz - tr/3.0;
	}


double pkdMassCheck(PKD pkd) 
{
	double dMass=0.0;
	int i,iRej;
	int nLocal = pkdLocal(pkd);

	for (i=0;i<nLocal;++i) {
		dMass += pkd->pStore[i].fMass;
		}
	iRej = pkdFreeStore(pkd) - pkd->nRejects;
	for (i=0;i<pkd->nRejects;++i) {
		dMass += pkd->pStore[iRej+i].fMass;
		}
	return(dMass);
	}

void pkdMassMetalsEnergyCheck(PKD pkd, double *dTotMass, double *dTotMetals, 
                    double *dTotOx, double *dTotFe, double *dTotEnergy) 
{
	int i;
	*dTotMass=0.0;
	*dTotMetals=0.0;
	*dTotOx=0.0;
	*dTotFe=0.0;
	*dTotEnergy=0.0;

	for (i=0;i<pkdLocal(pkd);++i) {
		*dTotMass += pkd->pStore[i].fMass;
#ifdef GASOLINE 
                *dTotMetals += pkd->pStore[i].fMass*pkd->pStore[i].fMetals;
#ifdef STARFORM
                *dTotOx += pkd->pStore[i].fMass*pkd->pStore[i].fMFracOxygen;
                *dTotFe += pkd->pStore[i].fMass*pkd->pStore[i].fMFracIron;
                if ( TYPETest(&pkd->pStore[i], TYPE_GAS) ){
		  *dTotEnergy += pkd->pStore[i].fMass*pkd->pStore[i].uDotFB;
                    }
#endif
#endif
		}
	}


void
pkdSetRung(PKD pkd, int iRung)
{
    int i;
    
    for(i=0;i<pkdLocal(pkd);++i) {
		pkd->pStore[i].iRung = iRung;
		}
    }

void
pkdBallMax(PKD pkd, int iRung, int bGreater, double dhFac)
{
    int i;
    
    for (i=0;i<pkdLocal(pkd);++i) {
		if (TYPETest(&(pkd->pStore[i]),TYPE_GAS) &&
			(pkd->pStore[i].iRung == iRung ||
			 (bGreater && pkd->pStore[i].iRung > iRung))) {
			pkd->pStore[i].fBallMax = sqrt(pkd->pStore[i].fBall2)*dhFac;
			}
		}
    return;
    }

int
pkdActiveRung(PKD pkd, int iRung, int bGreater)
{
    int i;
    int nActive;
    char out[128];
    int nLocal = pkdLocal(pkd);
    
    nActive = 0;
    for (i=0;i<nLocal;++i) {
		if(pkd->pStore[i].iRung == iRung ||
		   (bGreater && pkd->pStore[i].iRung > iRung)) {
			TYPESet(&(pkd->pStore[i]),TYPE_ACTIVE);
			++nActive;
			}
		else
			TYPEReset(&(pkd->pStore[i]),TYPE_ACTIVE);
		}
    sprintf(out, "nActive: %d\n", nActive);
    mdlDiag(pkd->mdl, out);
    pkd->nActive = nActive;
    return nActive;
    }

int
pkdCurrRung(PKD pkd, int iRung)
{
    int i;
    int iCurrent;
    
    iCurrent = 0;
    for(i=0;i<pkdLocal(pkd);++i) {
		if(pkd->pStore[i].iRung == iRung) {
			iCurrent = 1;
			break;
			}
		}
    return iCurrent;
    }

void
pkdGravStep(PKD pkd,double dEta)
{
	double dt;
    int i;

    for (i=0;i<pkdLocal(pkd);i++) {
		if (TYPEQueryACTIVE(&(pkd->pStore[i]))) {
			mdlassert(pkd->mdl,pkd->pStore[i].dtGrav > 0.0);
			dt = dEta/sqrt(pkd->pStore[i].dtGrav);
			if (dt < pkd->pStore[i].dt)
				pkd->pStore[i].dt = dt;
			}
		}
    }

void
pkdAccelStep(PKD pkd,double dEta,double dVelFac,double dAccFac,int bDoGravity,
             int bEpsAcc,int bSqrtPhi,double dhMinOverSoft)
{
    int i;
    double vel;
    double acc;
    int j;
    double dT;

    for (i=0;i<pkdLocal(pkd);++i) {
		if (TYPEQueryACTIVE(&(pkd->pStore[i]))) {
			vel = 0;
			acc = 0;
			for (j=0;j<3;j++) {
				vel += pkd->pStore[i].v[j]*pkd->pStore[i].v[j];
				acc += pkd->pStore[i].a[j]*pkd->pStore[i].a[j];
				}
			mdlassert(pkd->mdl,vel >= 0.0);
			vel = sqrt(vel)*dVelFac;
			acc = sqrt(acc)*dAccFac;
			dT = FLOAT_MAXVAL;
			if (bEpsAcc && acc>0) {
#ifdef GASOLINE
#ifdef EPSACCH
			    if (pkdIsGas(pkd, &(pkd->pStore[i]))) {
                    dT = dEta*sqrt(sqrt(0.25*pkd->pStore[i].fBall2)/acc);
                    }		        
#else
			    if (pkdIsGas(pkd, &(pkd->pStore[i])) && dhMinOverSoft < 1 && pkd->pStore[i].fBall2<4.0*pkd->pStore[i].fSoft*pkd->pStore[i].fSoft) {
			        if (pkd->pStore[i].fBall2 > 4.0*dhMinOverSoft*dhMinOverSoft
                        *pkd->pStore[i].fSoft*pkd->pStore[i].fSoft) 
                        dT = dEta*sqrt(sqrt(0.25*pkd->pStore[i].fBall2)/acc);
			        else 
                        dT = dEta*sqrt((dhMinOverSoft*pkd->pStore[i].fSoft)/acc);
			        }
#endif 
                else 
#endif /* GASOLINE */
			        dT = dEta*sqrt(pkd->pStore[i].fSoft/acc);
                }
			if (bSqrtPhi && acc>0) {
				/*
				 ** NOTE: The factor of 3.5 keeps this criterion in sync
				 ** with DensityStep. The nominal value of dEta for both
				 ** cases is then 0.02-0.03.
				 */
				double dtemp =
					dEta*3.5*sqrt(dAccFac*fabs(pkd->pStore[i].fPot))/acc;
				if (dtemp < dT)
					dT = dtemp;
				}
			if (dT < pkd->pStore[i].dt)
				pkd->pStore[i].dt = dT;
			}
		}
    }

void
pkdDensityStep(PKD pkd,double dEta,double dRhoFac)
{
    int i;
    double dT;
    
    for (i=0;i<pkdLocal(pkd);++i) {
		if (TYPEQueryACTIVE(&(pkd->pStore[i]))) {
			dT = dEta/sqrt(pkd->pStore[i].fDensity*dRhoFac);
			if (dT < pkd->pStore[i].dt)
				pkd->pStore[i].dt = dT;
			}
		}
    }


int
pkdOneParticleDtToRung( int iRung,double dDelta,double dt)
    {
    int iSteps,iTempRung;
    double dSteps = dDelta/dt;
    
    assert(dSteps < 2.1e9); /* avoid integer overflow */
    
    iSteps = floor(dSteps);
    /* insure that integer boundary goes
       to the lower rung. */
    if(fmod(dDelta,dt) == 0.0) iSteps--;
    
    iTempRung = iRung;
    if(iSteps < 0)
        iSteps = 0;
    
    while(iSteps) {
        ++iTempRung;
        iSteps >>= 1;
        }
    
    return iTempRung;
    }


int
pkdDtToRung(PKD pkd,int iRung,double dDelta,int iMaxRung,
    int bAll, /* 0 => symplectic case */
    int bDiagExceed, /* Diagnostics -- only in iMaxRung not reduced */
    int *pnMaxRung,	/* number of particles on MaxRung */
    int *piMaxRungIdeal)  /* preferred max rung */
    {

    int i;
    int iMaxRungOut;
    int iTempRung;
    int nMaxRung,nExceed;
    int iMaxRungIdeal;
    int bDiag;

    nExceed = (bDiagExceed ? 100 : 0);
    
    iMaxRungOut = 0;
    iMaxRungIdeal = 0;
    nMaxRung = 0;
    bDiag = 0;
    for(i=0;i<pkdLocal(pkd);++i) {
		if(pkd->pStore[i].iRung >= iRung) {
			mdlassert(pkd->mdl,TYPEQueryACTIVE(&(pkd->pStore[i])));
			if(bAll) {          /* Assign all rungs at iRung and above */
                assert(pkd->pStore[i].fSoft > 0.0);
                assert(pkd->pStore[i].dt > 0.0);
				iTempRung = pkdOneParticleDtToRung( iRung,dDelta,pkd->pStore[i].dt);

				if(iTempRung >= iMaxRungIdeal)
					iMaxRungIdeal = iTempRung+1;
				if(iTempRung >= iMaxRung) {
                    if (nExceed) {
                        nExceed--;
                        bDiag = 1;
                        }
					iTempRung = iMaxRung-1;
                    }
                
#ifdef GASOLINE
                if (bDiag) {
                    bDiag = 0;
                    if (pkdIsGas(pkd, &pkd->pStore[i])) {
                        PARTICLE *p = &pkd->pStore[i];
                        double ph = sqrt(p->fBall2*0.25);
#ifndef NOCOOLING
                        double pTemp = CoolCodeEnergyToTemperature( pkd->Cool, &p->CoolParticle, p->u, p->fMetals );
#ifdef UNONCOOL
                        double pTempTot = CoolCodeEnergyToTemperature( pkd->Cool, &p->CoolParticle, p->u+p->uNoncool, p->fMetals );
#else
                        double pTempTot = p->u;
#endif
#else
                        double pTemp = p->u, pTempTot = -1;
#endif
                        double h1=0,h2=0;
#ifdef STARFORM
                        h1 = p->uDotFB;
#endif
#ifndef NOCOOLING
                        h2 = p->uDot;
#endif
                        fprintf(stderr,"p %d exceeds maxrung: %g %g %g  T %g %g udot %g %g h %g %g dt %g %g %g %g %g %g %g\n",p->iOrder,p->fDensity,p->c,sqrt(p->v[0]*p->v[0]+p->v[1]*p->v[1]+p->v[2]*p->v[2]),pTemp,pTempTot,h1,h2,ph,p->fSoft,p->dt,0.4*(ph/(p->c + 0.6*(p->c))),0.4*(ph/(p->c + 0.6*(p->c + 2*p->mumax))),0.2*sqrt(ph/sqrt(p->a[0]*p->a[0]+p->a[1]*p->a[1]+p->a[2]*p->a[2])),0.25*p->u/UDOT_HYDRO(p),1/2.8*ph*ph/(DIFFRATE(p)),p->dtOld);
                        }
                    }
#endif

				pkd->pStore[i].iRung = iTempRung;
				pkd->pStore[i].dtOld = pkd->pStore[i].dt;
				}
			else {
				if(dDelta <= pkd->pStore[i].dt) {
					pkd->pStore[i].iRung = iRung;
					}
				else {
					pkd->pStore[i].iRung = iRung+1;
					}
                }
			}
		if(pkd->pStore[i].iRung > iMaxRungOut) {
			iMaxRungOut = pkd->pStore[i].iRung;
		        nMaxRung = 1;
                        }
		else if (pkd->pStore[i].iRung == iMaxRungOut) 
		        nMaxRung ++;
		        
		}

    *pnMaxRung = nMaxRung;
    *piMaxRungIdeal = iMaxRungIdeal;
    return iMaxRungOut;
    }

void
pkdInitDt(PKD pkd,double dDelta)
{
    int i;
    
    for(i=0;i<pkdLocal(pkd);++i) {
		if(TYPEQueryACTIVE(&(pkd->pStore[i])))
			pkd->pStore[i].dt = dDelta;
		}
    }
    
int pkdRungParticles(PKD pkd,int iRung)
{
	int i,n;

	n = 0;
	for (i=0;i<pkdLocal(pkd);++i) {
		if (pkd->pStore[i].iRung == iRung) ++n;
		}
	return n;
	}

void 
pkdMoveParticle(PKD pkd, double *xcenter,double *xoffset,int iOrder)
{
    int i,j;
    
    for (i=0;i<pkdLocal(pkd);i++) {
	if (pkd->pStore[i].iOrder == iOrder) {
	    
#ifdef SLIDING_PATCH
	    pkd->pStore[i].v[1]+=1.5*pkd->PP->dOrbFreq*pkd->pStore[i].r[0];
#endif	
	    for (j=0;j<3;j++)
		pkd->pStore[i].r[j]=xcenter[j]+xoffset[j];
#ifdef SLIDING_PATCH
	    pkd->pStore[i].v[1]-=1.5*pkd->PP->dOrbFreq*pkd->pStore[i].r[0];
#endif
	    }
	}
    }


void
pkdUnDeleteParticle(PKD pkd, PARTICLE *p)
{
    assert(TYPETest(p, TYPE_DELETED)); 

    TYPEReset(p, TYPE_DELETED); 
    p->iOrder = -2 - p->iOrder; 
    }

void
pkdDeleteParticle(PKD pkd, PARTICLE *p)
{
    p->iOrder = -2 - p->iOrder; 
    TYPEClearACTIVE(p); 
    TYPESet(p, TYPE_DELETED); 
    }

void
pkdNewParticle(PKD pkd, PARTICLE p)
{
    mdlassert(pkd->mdl,pkd->nLocal < pkd->nStore);
    pkd->pStore[pkd->nLocal] = p;
    pkd->pStore[pkd->nLocal].iOrder = -1;
    pkd->nLocal++;
    }

void
pkdColNParts(PKD pkd, int *pnNew, int *nAddGas, int *nAddDark,
	     int *nAddStar, int *nDelGas, int *nDelDark, int *nDelStar)
{
    int pi, pj;
    int nNew = 0;
    int naddGas = 0;
    int naddDark = 0;
    int naddStar = 0;
    int ndelGas = 0;
    int ndelDark = 0;
    int ndelStar = 0;
    int newnLocal;
    PARTICLE *p;
    
    newnLocal = pkdLocal(pkd);
    for(pi = 0, pj = 0; pi < pkdLocal(pkd); pi++) {
	if(pj < pi)
	    pkd->pStore[pj] = pkd->pStore[pi];
	p = &pkd->pStore[pi];
	if(p->iOrder == -1) {
	    ++pj;
	    ++nNew;
#if NIORDERGASBUFFER
	    if (TYPETest(p,TYPE_GAS)) 
		++naddGas;
	    else if (TYPETest(p,TYPE_STAR)) 
		++naddStar;
	    else if (TYPETest(p,TYPE_DARK)) 
		++naddDark;
	    else assert(0);
#else
#ifdef GASOLINE
	    ++naddStar;
#else
	    ++naddDark;
#endif
#endif
	    if(TYPEQueryACTIVE(p))
		++pkd->nActive;
	    continue;
	    }
	else if(p->iOrder < -1){
	    --newnLocal;
#if NIORDERGASBUFFER
	    p->iOrder = 2000000000;
	    if (TYPETest(p,TYPE_GAS)) 
		++ndelGas;
	    else if (TYPETest(p,TYPE_STAR)) 
		++ndelStar;
	    else if (TYPETest(p,TYPE_DARK)) 
		++ndelDark;
	    else assert(0);
#else
	    p->iOrder = -2 - p->iOrder;
	    if(pkdIsGas(pkd, p))
		++ndelGas;
	    else if(pkdIsDark(pkd, p))
		++ndelDark;
	    else if(pkdIsStar(pkd, p))
		++ndelStar;
	    else
		mdlassert(pkd->mdl,0);
#endif
	    if(TYPEQueryACTIVE(p))
		--pkd->nActive;
	    }
	else {
	    ++pj;
	    }
	}

    *pnNew = nNew;
    *nAddGas = naddGas;
    *nAddDark = naddDark;
    *nAddStar = naddStar;
    *nDelGas = ndelGas;
    *nDelDark = ndelDark;
    *nDelStar = ndelStar;
    pkd->nLocal = newnLocal;
    }

void
pkdNewOrder(PKD pkd,int nStartGas, int nStartDark, int nStartStar) 
/* How does it feel? */
{
    int pi;
    PARTICLE *p;

    /* Sink Log to be updated? */
    if (pkd->sinkLog.nLog) {
	SINKEVENT *pSE;
	int iOrdSink = nStartStar;
	int iForm, iLog;
	iForm = pkd->sinkLog.nFormOrdered;
	iLog = pkd->sinkLog.nLogOrdered;
	pSE = &(pkd->sinkLog.SinkEventTab[iLog]);
	for(pi=0;pi<pkdLocal(pkd);pi++) {
	    p = &(pkd->pStore[pi]);
	    if(p->iOrder == -1) {
		if (pkdIsStar(pkd, p)) {
		    while (pSE->iOrdSink != -1) { /* Find next sink creation accretion in log */
			printf("%d New order Non-sink form accrete -- skipped %d %d %d\n",pkd->idSelf,iLog,pSE->iOrdSink,pSE->iOrdVictim);
			pSE++;
			iLog++;
			assert(iLog < pkd->sinkLog.nLog);
			}
		    while (pSE->iOrdSink == -1) { /* Count through sink creation accretions in log */
			printf("%d New order Sink form accrete -- skipped %d %d %d %d\n",pkd->idSelf,iLog,pSE->iOrdSink,pSE->iOrdVictim,pkd->sinkLog.nAccrete);
			pSE++;
			iLog++;
			}
		    printf("%d New order Sink form -- %d %d %d %d %d\n",pkd->idSelf,iLog,pSE->iOrdSink,pSE->iOrdVictim,iForm,pkd->sinkLog.nForm);
		    /* Must be followed by sink creation event in log */
		    assert(iLog < pkd->sinkLog.nLog);
		    assert(pSE->iOrdSink == iForm);
		    assert(pSE->iOrdVictim == -1);
		    pSE->iOrdSink = iOrdSink;
		    pSE++;
		    iForm++; iLog++;
		    pkd->sinkLog.nLogOrdered = iLog;

		    iOrdSink++; /* Just counted here -- set below */
		    pkd->sinkLog.nFormOrdered++;
		    }
		}
	    }
	printf("%d New order sink ordering -- %d %d %d\n",pkd->idSelf,iForm,pkd->sinkLog.nFormOrdered,pkd->sinkLog.nForm);
	assert(iForm == pkd->sinkLog.nForm);
	}
    /* END sink Log */

    for(pi=0;pi<pkdLocal(pkd);pi++) {
	p = &(pkd->pStore[pi]);
	if(p->iOrder == -1) {
	    if (pkdIsGas(pkd, p)) 
		p->iOrder = nStartGas++;
	    else if (pkdIsDark(pkd, p)) 
		p->iOrder = nStartDark++;
	    else {
#ifdef STARFORM
            /* Also record iOrder in the starLog table. */
            pkd->starLog.seTab[pkd->starLog.nOrdered].iOrdStar = nStartStar;
            pkd->starLog.nOrdered++;
            assert(pkd->starLog.nOrdered <= pkd->starLog.nLog);
#endif
            p->iOrder = nStartStar++;
            }
	    }
	}
    }

void
pkdGetNParts(PKD pkd, struct outGetNParts *out )
{
    int pi;
    int n;
    int nGas;
    int nDark;
    int nStar;
	int iMaxOrderGas;
	int iMaxOrderDark;
	int iMaxOrderStar;
    PARTICLE *p;
    
    n = 0;
    nGas = 0;
    nDark = 0;
    nStar = 0;
    iMaxOrderGas = -1;
    iMaxOrderDark = -1;
    iMaxOrderStar = -1;
    for(pi = 0; pi < pkdLocal(pkd); pi++) {
		p = &pkd->pStore[pi];
		n++;
	    if(pkdIsGas(pkd, p)) {
			++nGas;
			if (p->iOrder > iMaxOrderGas) iMaxOrderGas = p->iOrder;
			}
	    else if(pkdIsDark(pkd, p)) {
			++nDark;
			if (p->iOrder > iMaxOrderDark) iMaxOrderDark = p->iOrder;
			}
	    else if(pkdIsStar(pkd, p)) {
			++nStar;
			if (p->iOrder > iMaxOrderStar) iMaxOrderStar = p->iOrder;
			}
		}

	out->n  = n;
    out->nGas = nGas;
    out->nDark = nDark;
    out->nStar = nStar;
    out->iMaxOrderGas = iMaxOrderGas;
    out->iMaxOrderDark = iMaxOrderDark;
    out->iMaxOrderStar = iMaxOrderStar;
    }

void
pkdSetNParts(PKD pkd,int nGas,int nDark,int nStar,int nMaxOrder, int nMaxOrderGas,
	     int nMaxOrderDark)
{
    pkd->nGas = nGas;
    pkd->nDark = nDark;
    pkd->nStar = nStar;
    pkd->nMaxOrder = nMaxOrder; 
    pkd->nMaxOrderGas = nMaxOrderGas;
    pkd->nMaxOrderDark = nMaxOrderDark;
}

void pkdCoolVelocity(PKD pkd,int nSuperCool,double dCoolFac,
					 double dCoolDens,double dCoolMaxDens)
{
	PARTICLE *p;
	int i,j,n;

	p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;++i) {
		if (p[i].iOrder < nSuperCool && p[i].fDensity >= dCoolDens &&
			p[i].fDensity < dCoolMaxDens) {
			for (j=0;j<3;++j) {
#ifdef SUPERCOOL
				p[i].v[j] = p[i].vMean[j] + (p[i].v[j]-p[i].vMean[j])*dCoolFac;
#else
				p[i].v[j] *= dCoolFac;
#endif
				}
			}
		}
	}

int pkdActiveExactType(PKD pkd, unsigned int iFilterMask, unsigned int iTestMask, unsigned int iSetMask)
{
    PARTICLE *p;
    int i, nActive = 0;
    int nLocal = pkdLocal(pkd);

    for(i=0;i<nLocal;++i) { 
		p = &pkd->pStore[i];
		/* DEBUG: Paranoia check */
		mdlassert(pkd->mdl,TYPETest(p,TYPE_ALL));
		if (TYPEFilter(p,iFilterMask,iTestMask)) {
			TYPESet(p,iSetMask);
			nActive++;
			}
		else {
			TYPEReset(p,iSetMask);
			}
		}
    if (iSetMask & TYPE_ACTIVE) pkd->nActive = nActive;
    if (iSetMask & TYPE_TREEACTIVE) pkd->nTreeActive = nActive;
    if (iSetMask & TYPE_SMOOTHACTIVE) pkd->nSmoothActive = nActive;
    return nActive;
    }

int pkdSetType(PKD pkd, unsigned int iTestMask, unsigned int iSetMask)
{
    PARTICLE *p;
    int i, nActive = 0;

    for(i=0;i<pkdLocal(pkd);++i) {
		p = &pkd->pStore[i];
#ifdef COOLDEBUG
		if (p->iOrder == 842079) fprintf(stderr,"Particle %i in pStore[%i]\n",p->iOrder,(int) (p-pkd->pStore));
		assert(p->u >= 0.0);
		assert(p->uPred >= 0.0);
#endif
		/* DEBUG: Paranoia check */
		mdlassert(pkd->mdl,TYPETest(p,TYPE_ALL));
		if (TYPETest(p,iTestMask)) {
			TYPESet(p,iSetMask);
			nActive++;
			}
		}
    /*
	   Need to fix this up:
	   if (iSetMask & TYPE_ACTIVE) pkd->nActive = nActive;
	   if (iSetMask & TYPE_TREEACTIVE) pkd->nTreeActive = nActive;
	   if (iSetMask & TYPE_SMOOTHACTIVE) pkd->nSmoothActive = nActive;
	   */
    return nActive;
    }

int pkdResetType(PKD pkd, unsigned int iTestMask, unsigned int iSetMask)
{
    PARTICLE *p;
    int i, nActive = 0;

    for(i=0;i<pkdLocal(pkd);++i) {
		p = &pkd->pStore[i];
		/* DEBUG: Paranoia check */
		mdlassert(pkd->mdl,TYPETest(p,TYPE_ALL));
		if (TYPETest(p,iTestMask)) {
			TYPEReset(p,iSetMask);
			nActive++;
			}
		}
    /*
	   if (iSetMask & TYPE_ACTIVE) pkd->nActive = nActive;
	   if (iSetMask & TYPE_TREEACTIVE) pkd->nTreeActive = nActive;
	   if (iSetMask & TYPE_SMOOTHACTIVE) pkd->nSmoothActive = nActive;
	   */
    return nActive;
    }

int pkdCountType(PKD pkd, unsigned int iFilterMask, unsigned int iTestMask)
{
    PARTICLE *p;
    int i, nActive = 0;
    int nLocal = pkdLocal(pkd);

    for(i=0;i<nLocal;++i) {
		p = &pkd->pStore[i];
		if (TYPEFilter(p,iFilterMask,iTestMask)) {
			nActive++;
			}
		}
    {
	char debug[100];
	sprintf(debug, "Filter %d:%d, Counted: %d\n",iFilterMask,iTestMask,nActive);
	mdlDiag(pkd->mdl,debug);
	}
    return nActive;
    }

int pkdActiveType(PKD pkd, unsigned int iTestMask, unsigned int iSetMask)
{
    PARTICLE *p;
    int i, nActive = 0;
    int nLocal = pkdLocal(pkd);

    for(i=0;i<nLocal;++i) {
		p = &pkd->pStore[i];
		/* DEBUG: Paranoia check */
		/* mdlassert(pkd->mdl,TYPETest(p,TYPE_ALL)); */
		if (TYPETest(p,iTestMask)) {
			TYPESet(p,iSetMask);
			nActive++;
			}
		else {
			TYPEReset(p,iSetMask);
			}
		}
    if (iSetMask & TYPE_ACTIVE      ) pkd->nActive       = nActive;
    if (iSetMask & TYPE_TREEACTIVE  ) pkd->nTreeActive   = nActive;
    if (iSetMask & TYPE_SMOOTHACTIVE) pkd->nSmoothActive = nActive;
    return nActive;
    }

int
pkdActiveMaskRung(PKD pkd, unsigned iSetMask, int iRung, int bGreater)
{
    PARTICLE *p;
    int i;
    int nActive;
    char out[128];
    
    nActive = 0;
    for(i=0;i<pkdLocal(pkd);++i) {
        p = &pkd->pStore[i];
		if(p->iRung == iRung || (bGreater && p->iRung > iRung)) {
			TYPESet(p,iSetMask);
			++nActive;
			}
		else
			TYPEReset( p, iSetMask );
		}
    sprintf(out,"nActive: %d\n",nActive);
    mdlDiag(pkd->mdl,out);

    if ( iSetMask & TYPE_ACTIVE      ) pkd->nActive       = nActive;
    if ( iSetMask & TYPE_TREEACTIVE  ) pkd->nTreeActive   = nActive;
    if ( iSetMask & TYPE_SMOOTHACTIVE) pkd->nSmoothActive = nActive;
    return nActive;
    }

int
pkdActiveTypeRung(PKD pkd, unsigned iTestMask, unsigned iSetMask, int iRung, int bGreater)
{
    PARTICLE *p;
    int i;
    int nActive;
    char out[128];
    int nLocal = pkdLocal(pkd);
    
    nActive = 0;
    for(i=0;i<nLocal;++i) {
        p = &pkd->pStore[i];
        /* DEBUG: Paranoia check */
        mdlassert(pkd->mdl,TYPETest(p,TYPE_ALL));
		if(TYPETest(p,iTestMask) && 
           (p->iRung == iRung || (bGreater && p->iRung > iRung))) {
			TYPESet(p,iSetMask);
			++nActive;
			}
		else
			TYPEReset( p, iSetMask );
		}
    sprintf(out,"nActive: %d\n",nActive);
    mdlDiag(pkd->mdl,out);

    if ( iSetMask & TYPE_ACTIVE      ) pkd->nActive       = nActive;
    if ( iSetMask & TYPE_TREEACTIVE  ) pkd->nTreeActive   = nActive;
    if ( iSetMask & TYPE_SMOOTHACTIVE) pkd->nSmoothActive = nActive;
    return nActive;
    }

struct SortStruct {
  int iOrder;
  int iStore;
};

int CompSortStruct(const void * a, const void * b) {
  return ( ( ((struct SortStruct *) a)->iOrder < ((struct SortStruct *) b)->iOrder ? -1 : 1 ) );
}

void pkdSetTypeFromFileSweep(PKD pkd, int iSetMask, char *file, 
	   struct SortStruct *ss, int nss, int *pniOrder, int *pnSet) {
  int niOrder = 0, nSet = 0;
  int iOrder, iOrderOld, nRet;
  FILE *fp;
  int iss;

  fp = fopen( file, "r" );
  assert( fp != NULL );

  iss = 0;
  iOrderOld = -1;
  while ( (nRet=fscanf( fp, "%d\n", &iOrder )) == 1 ) {
	niOrder++;
	assert( iOrder > iOrderOld );
	iOrderOld = iOrder;
	while (ss[iss].iOrder < iOrder) {
	  iss++;
	  if (iss >= nss) goto DoneSS;
	}
	if (iOrder == ss[iss].iOrder) {
	  TYPESet(&(pkd->pStore[ss[iss].iStore]),iSetMask);
	  nSet++;
	}
  }
 
DoneSS:
  /* Finish reading file to verify consistency across processors */
  while ( (nRet=fscanf( fp, "%d\n", &iOrder )) == 1 ) {
	niOrder++;
	assert( iOrder > iOrderOld );
	iOrderOld = iOrder;
	}
  fclose(fp);

  *pniOrder += niOrder;
  *pnSet += nSet;

  return;
}


int pkdSetTypeFromFile(PKD pkd, int iSetMask, int biGasOrder, char *file, int *pniOrder, int *pnSet, int *pnSetiGasOrder)
{
  struct SortStruct *ss;
  int i,nss;

  *pniOrder = 0;
  *pnSet = 0;
  *pnSetiGasOrder = 0;

  nss = pkdLocal(pkd);
  ss = malloc(sizeof(*ss)*nss);
  assert( ss != NULL );

  for(i=0;i<pkdLocal(pkd);++i) {
	ss[i].iOrder = 	pkd->pStore[i].iOrder;
	ss[i].iStore = i;
  }

  qsort( ss, nss, sizeof(*ss), CompSortStruct );

  pkdSetTypeFromFileSweep(pkd, iSetMask, file, ss, nss, pniOrder, pnSet);

#if defined(STARFORM) || defined(SIMPLESF)
  if (biGasOrder) {
	for(i=0;i<pkdLocal(pkd);++i) {
	  ss[i].iOrder =   ( TYPETest(&pkd->pStore[i], TYPE_STAR) ? 
						 pkd->pStore[i].iGasOrder : -1 );
	  ss[i].iStore = i;
	}
	
	qsort( ss, nss, sizeof(*ss), CompSortStruct );
	
	pkdSetTypeFromFileSweep(pkd, iSetMask, file, ss, nss, pniOrder, pnSetiGasOrder);
  }
#endif

  free( ss );

  return *pnSet+*pnSetiGasOrder;
}


void pkdSetParticleTypes(PKD pkd, int nSuperCool)
{
    PARTICLE *p;
    int i, iSetMask;
    
    for(i=0;i<pkdLocal(pkd);++i) {
		p = &pkd->pStore[i];
		iSetMask = 0;
		if (pkdIsGasByOrder(pkd,p)) iSetMask |= TYPE_GAS;
		if (pkdIsDarkByOrder(pkd,p)) iSetMask |= TYPE_DARK;
		if (pkdIsStarByOrder(pkd,p)) iSetMask |= TYPE_STAR;
		if (p->iOrder < nSuperCool) iSetMask |= TYPE_SUPERCOOL;
#ifdef STARSINK
		if (iSetMask & TYPE_STAR) {
		    iSetMask |= TYPE_SINK;
#ifdef SINKING
		    if (p->iSinkingOnto < 0) p->iSinkingOnto=0;
		    
#endif
		    }
#endif
#ifdef SINKING
		if ((iSetMask&TYPE_GAS) && p->iSinkingOnto >= 0) iSetMask |= TYPE_SINKING;
#endif
		TYPESet(p,iSetMask);
		}
    } 

int
pkdSoughtParticleList(PKD pkd, int iTypeSought, int nMax, int *n, struct SoughtParticle *sp)
{
    int i, nFound;
	PARTICLE *p = pkd->pStore;

	nFound = 0;
    for(i=0;i<pkdLocal(pkd);++i) {
 	  if (TYPETest((p+i),iTypeSought)) {
		/* Don't halt on overflow -- just stop copying */
		if (nFound < nMax) {
		  sp[nFound].iOrder = p[i].iOrder;
		  sp[nFound].iActive = p[i].iActive;
		  sp[nFound].x = p[i].r[0];
		  sp[nFound].y = p[i].r[1];
		  sp[nFound].z = p[i].r[2];
		  sp[nFound].m = p[i].fMass;
		}
		nFound ++;
	  }
	}
	*n = nFound;
	for (i=0;i<nFound;i++) {
	  printf("pkd sub star %d: %g %g %g  %g\n",i,sp[i].x,sp[i].y,sp[i].z,sp[i].m);
	}

	return (nFound);
}

void
pkdCoolUsingParticleList(PKD pkd, int nList, struct SoughtParticle *l)
    {
#ifndef NOCOOLING
    int i,j;
	double r2,r2min,dx,mass_min;
	PARTICLE *p = pkd->pStore;

	assert(nList > 0);
    for(i=0;i<pkdLocal(pkd);++i) {
 	  if (TYPEQueryACTIVE(&p[i]) && TYPETest(&p[i],TYPE_GAS)) {
		dx = p[i].r[0]-l[0].x;
		r2 = dx*dx;
		dx = p[i].r[1]-l[0].y;
		r2 += dx*dx;
		dx = p[i].r[2]-l[0].z;
		r2 += dx*dx;
		r2min = r2;
		mass_min = l[0].m;
		for (j=1;j<nList;j++) {
		  dx = p[i].r[0]-l[j].x;
		  r2 = dx*dx;
		  dx = p[i].r[1]-l[j].y;
		  r2 += dx*dx;
		  dx = p[i].r[2]-l[j].z;
		  r2 += dx*dx;
		  if (r2 < r2min) { 
		      r2min = r2;
		      mass_min = l[0].m;
		      }
		}
#ifdef COOLING_DISK
		p[i].CoolParticle.r = sqrt(r2min);
		p[i].CoolParticle.StarMass = mass_min;
#endif
	  }
	}
#endif
}

void pkdGrowMass(PKD pkd,int nGrowMass, int iGrowType, double dDeltaM, double dMinM, double dMaxM)
{
    int i;
    PARTICLE *p;

    for(i=0;i<pkdLocal(pkd);++i) {
	p = pkd->pStore + i;
	if (TYPETest(p, iGrowType) && p->iOrder < nGrowMass) {
	    p->fMass += dDeltaM;
	    if (p->fMass < dMinM) p->fMass = dMinM;
	    else if (p->fMass > dMaxM) p->fMass = dMaxM;
	    }
	}
    }

void pkdInitAccel(PKD pkd)
{
    int i,j;
    
    for(i=0;i<pkdLocal(pkd);++i) {
		if (TYPEQueryACTIVE(&(pkd->pStore[i]))) {
  		        pkd->pStore[i].fPot = 0;
			for (j=0;j<3;++j) {
				pkd->pStore[i].a[j] = 0;
#if defined(GASOLINE) && defined(SHOCKTRACK)
				pkd->pStore[i].aPres[j] = 0;
#endif
				pkd->pStore[i].dtGrav = 0;
				}
			}
		}
    }

void pkdModifyAccel(PKD pkd, double a)
{
    int i,j;
    PARTICLE *p;
    FLOAT dxSmooth, xSmooth, factor;
    dxSmooth = (pkd->dxInflow+pkd->fPeriod[0]*.5);
    xSmooth = pkd->dxInflow + dxSmooth;

    for(i=0;i<pkdLocal(pkd);++i) {
      p = &(pkd->pStore[i]);
      if (p->r[0] < pkd->dxInflow) {
#ifdef GASOLINE
/*	    factor = 0;
	    if (p->r[0] > pkd->dxInflow && p->r[0] < xSmooth) factor = 1-(p->r[0]-pkd->dxInflow)/dxSmooth;*/
	    pkd->pStore[i].uDotPdV = 0;
	    pkd->pStore[i].uDotAV = 0;
	    pkd->pStore[i].uDotDiff = 0;
#endif
	    pkd->pStore[i].a[0] = a; 
	    pkd->pStore[i].a[1] = 0; 
	    pkd->pStore[i].a[2] = 0; 
	    }
	else if (p->r[0] > pkd->dxOutflow) {
#ifdef GASOLINE
/*	    factor = 0;
	    if (p->r[0] > pkd->dxInflow && p->r[0] < xSmooth) factor = 1-(p->r[0]-pkd->dxInflow)/dxSmooth;*/
	    pkd->pStore[i].uDotPdV = 0;
	    pkd->pStore[i].uDotAV = 0;
	    pkd->pStore[i].uDotDiff = 0;
#endif
	    pkd->pStore[i].a[0] = a; 
	    pkd->pStore[i].a[1] = 0; 
	    pkd->pStore[i].a[2] = 0; 
	    }
	else if (TYPEQueryACTIVE(&(pkd->pStore[i]))) {
	    pkd->pStore[i].a[0] += a; /* Note: Added, not set in this case */
	    }
	}
}

int pkdIsGasByOrder(PKD pkd,PARTICLE *p) {
	if (p->iOrder <= pkd->nMaxOrderGas) return 1;
	else return 0;
	}

int pkdIsDarkByOrder(PKD pkd,PARTICLE *p) {
	if (p->iOrder > pkd->nMaxOrderGas && p->iOrder <= pkd->nMaxOrderDark)
	    return 1;
	else return 0;
	}

int pkdIsStarByOrder(PKD pkd,PARTICLE *p) {
	if (p->iOrder > pkd->nMaxOrderDark) return 1;
	else return 0;
	}

#ifdef GASOLINE

#ifdef UNONCOOL
double pkduNoncoolConvRate(PKD pkd, UNCC uncc, FLOAT fBall2, double uNoncoolPred, double uPred) 
    {
    double rate,ueff;
    if (uncc.dNoncoolConvRate > 0) return uncc.dNoncoolConvRate;
    ueff = uNoncoolPred+uPred;
    if (ueff < uncc.dNoncoolConvUMin) ueff = uncc.dNoncoolConvUMin;
    rate = uncc.dNoncoolConvRateMul*sqrt((ueff)/(fBall2*0.25));
    if (rate > uncc.dNoncoolConvRateMax) return uncc.dNoncoolConvRateMax;
    return rate;
    }
#endif

void pkdUpdateuDot(PKD pkd, double duDelta, double dTime, double z, UNCC uncc, int iGasModel, int bUpdateState )
    {
#ifndef NOCOOLING	
    PARTICLE *p;
	int i,n;
	int bCool = 0;
	COOL *cl = NULL;
	COOLPARTICLE cp;
	double E,dt = 0,dtUse,uDotSansCooling;
#ifdef COOLING_MOLECULARH
	double correL = 1.0; /* Correlation length used when calculating shielding*/
#endif
#endif
    pkdClearTimer(pkd,1);
    pkdStartTimer(pkd,1);
    
#ifndef NOCOOLING
    switch (iGasModel) {
    case GASMODEL_COOLING:
        bCool = 1;
        cl = pkd->Cool;
        CoolSetTime( cl, dTime, z  );
        dt = CoolCodeTimeToSeconds( cl, duDelta );
        break;
        }
    
    p = pkd->pStore;
    n = pkdLocal(pkd);

    for (i=0;i<n;++i,++p) {
        if(TYPEFilter(p,TYPE_GAS|TYPE_ACTIVE,TYPE_GAS|TYPE_ACTIVE)) {
            double fDensity;

            double PoverRho,PoverRhoGas,PoverRhoNoncool,PoverRhoFloorJeans,cGas;
            double uNoncoolDotConv=0, uNoncoolPredTmp;
            double uNoncoolDotFB=0, uDotFBThermal=0, uDotPdVNJ, uncPdVFrac;
            int bUpdateStd=1;

#if defined(STARFORM)
#ifdef UNONCOOL
            uNoncoolDotFB = p->uDotFB; //Note: FB Added to uNoncool OR u (not both)
#else
            uDotFBThermal = p->uDotFB;
#endif
#endif /* STARFORM */

            pkdGasPressureParticle(pkd, &uncc.gpc, p, &PoverRhoFloorJeans, &PoverRhoNoncool, &PoverRhoGas, &cGas );
            PoverRho = PoverRhoGas + PoverRhoNoncool ;
            if (PoverRho < PoverRhoFloorJeans) PoverRho = PoverRhoFloorJeans;

#ifdef MASSNONCOOL        

            uncPdVFrac = p->fMassNoncool*p->uNoncoolPred;
            uncPdVFrac = uncPdVFrac/(uncPdVFrac+(p->fMass-p->fMassNoncool)*p->uPred);
            uDotPdVNJ = p->uDotPdV*(PoverRhoNoncool+PoverRhoGas)/(PONRHOFLOOR + PoverRho); /* remove JeansFloor */
            
            bUpdateStd = (p->fMassNoncool < 0.9*p->fMass);

            if (p->fMassNoncool > 0) {
                assert(p->uNoncool >= 0);
                
                uDotSansCooling = (uDotPdVNJ+p->uDotAV)*uncPdVFrac // Fraction of PdV related to uNoncool 
                    + p->uNoncoolDotDiff + uNoncoolDotFB;
				if ( bCool && p->uNoncool > 0) {
                /*if ( 0) {//I am paranoid about letting the hot phase cool, turned off by default*/
                    cp = p->CoolParticle;
                    E = p->uNoncool;
					dtUse = dt;
                    fDensity = PoverRhoGas/(uncc.gpc.gammam1*p->uNoncool); /* Density of bubble part of particle */
                    CoolIntegrateEnergyCode(cl, &cp, &E, uDotSansCooling, fDensity, p->fMetals, p->r, dtUse);
                    p->uNoncoolDot = (E - p->uNoncool)/duDelta;
                    if (bUpdateState && !bUpdateStd) p->CoolParticle = cp;
                    }
                else 
                    p->uNoncoolDot = uDotSansCooling;
                }
            else {
                p->uNoncoolDot = 0;
				uDotFBThermal = p->uDotFB;
			}

            assert(p->uPred > 0);
            fDensity = PoverRhoGas/(uncc.gpc.gammam1*p->uPred); /* Density of non-bubble part of particle */

            uDotSansCooling = (p->uDotPdV+p->uDotAV)*(1-uncPdVFrac) // Fraction of PdV related to u thermal
                    + p->uDotDiff;
#else /* !MASSNONCOOL */

            fDensity = p->fDensity; /* Density for cooling */
#ifdef DENSITYU
            if (p->fDensityU < fDensity) fDensity = p->fDensityU;
#endif
#ifdef UNONCOOL
            /* 2nd order estimator for Conv -- note that PdV etc ...should already be 2nd order via Leap Frog */
            uNoncoolPredTmp = p->uNoncool+p->uNoncoolDot*duDelta*0.5;
#ifdef UNONCOOLCONV
            uNoncoolDotConv = uNoncoolPredTmp*
                pkduNoncoolConvRate(pkd,uncc,p->fBall2,uNoncoolPredTmp,p->u+p->uDot*duDelta*0.5); 
#endif
            
            p->uNoncoolDot = p->uDotPdV*PoverRhoNoncool/(PONRHOFLOOR + PoverRho) // Fraction of PdV related to uNoncool 
                - uNoncoolDotConv + uNoncoolDotFB + p->uNoncoolDotDiff;
#endif    
            uDotSansCooling = p->uDotPdV*PoverRhoGas/(PONRHOFLOOR + PoverRho) // Fraction of PdV related to u thermal
                + p->uDotAV                                                   // Only u thermal energy gets shock heating
                + uNoncoolDotConv + uDotFBThermal + p->uDotDiff;
#endif /* !MASSNONCOOL */

            if ( bCool ) {
                cp = p->CoolParticle;
                E = p->u;
#ifdef COOLDEBUG
                cl->p = p; /* Send in particle pointer only for temporary debug */
#endif
                dtUse = dt;
#ifdef STARFORM
                if ( dTime < p->fTimeCoolIsOffUntil) {
                    dtUse = -dt;
                    p->uDot = uDotSansCooling;
                    }
#endif
#ifdef SIMPLESF 
                if (dTime < p->fTimeForm) {
                    dtUse = -dt;
                    if (p->uDot < UDOTHYDRO(p)) p->uDot = UDOTHYDRO(p);
                    }
#endif
                
#ifdef COOLING_MOLECULARH
#ifdef NEWSHEAR
				/***** particle diffusion method ******/
                if (p->diff != 0) 
                    correL = p->c * (0.25*p->fBall2)/p->diff; 
                if (correL > sqrt(0.25*p->fBall2) || p->diff == 0) 
                    correL = sqrt(0.25*p->fBall2); /*minimum is particle smoothing*/
#else /*NEWSHEAR*/
#ifdef PARTSHEAR
                /***** particle shear method********/ 
                if (p->curlv[0]*p->curlv[0] + p->curlv[1]*p->curlv[1] + p->curlv[2]*p->curlv[2] != 0) 
                    correL = p->c/sqrt(p->curlv[0]*p->curlv[0] + p->curlv[1]*p->curlv[1] + p->curlv[2]*p->curlv[2]);				
#else /*PARTSHEAR*/
                /*#ifdef COLUMNLENGTH */ /*Made using the smoothing length the default, as it has been used that way in all production runs to Jun 4th, 2012, CC*/
                /***** From particle smoothing.  This works best for determining the correlation length.  CC 7/20/11 ******/
                correL = sqrt(0.25*p->fBall2);
                /*#endif COLUMNLENGTH*/
#endif /*PARTSHEAR*/
#endif /*NEWSHEAR*/
                CoolIntegrateEnergyCode(cl, &cp, &E, uDotSansCooling, fDensity, p->fMetals, p->r, dtUse, correL); /* If doing H2, send the correlation length to calculate the shielding*/
#else /* !COOLING_MOLECULARH */
                
#ifdef COOLING_BOLEY
                cp.mrho = pow(p->fMass/p->fDensity, 1./3.);
#endif
                CoolIntegrateEnergyCode(cl, &cp, &E, uDotSansCooling, fDensity, p->fMetals, p->r, dtUse);
#endif /* !COOLING_MOLECULARH */
                
                mdlassert(pkd->mdl,E > 0);
                
                if (dtUse > 0 || uDotSansCooling*duDelta + p->u < 0) p->uDot = (E - p->u)/duDelta;
                if (bUpdateState && bUpdateStd) p->CoolParticle = cp;
                }
            else { 
                p->uDot = uDotSansCooling;
                }
            }
        }
#endif /*NOT NOCOOLING*/
    pkdStopTimer(pkd,1);
    }

void pkdUpdateShockTracker(PKD pkd, double dDelta, double dShockTrackerA, double dShockTrackerB )
{
#ifdef SHOCKTRACK
	PARTICLE *p;
	int i,n;
	double conh,factor;
	double dv2,a2,ap2,adotap=1;

        if (dShockTrackerB == 0) printf("Doing cheap shock tracking\n");
	p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;++i,++p) {
		if(TYPEFilter(p,TYPE_GAS|TYPE_ACTIVE,TYPE_GAS|TYPE_ACTIVE)) {
		        p->ShockTracker = 0;

			a2 = ((p->a[0]*p->a[0])+(p->a[1]*p->a[1])+(p->a[2]*p->a[2]));
			ap2 = ((p->aPres[0]*p->aPres[0])+(p->aPres[1]*p->aPres[1])+(p->aPres[2]*p->aPres[2]));
			adotap = ((p->a[0]*p->aPres[0])+(p->a[1]*p->aPres[1])+(p->a[2]*p->aPres[2]));

#if 0
			if (!(i%2000)) {
			  /*
printf("aP %d: %g %g %g %g %g %g\n",i,p->a[0],p->a[1],p->a[2],p->aPres[0],p->aPres[1],p->aPres[2]);
			  */
printf("r %g PdV %g a %g %g %g SW %g %g %g\n",sqrt(p->r[0]*p->r[0]+p->r[1]*p->r[1]+p->r[2]*p->r[2]),
       UDOT_HYDRO(p),((p->a[0]*p->aPres[0])+(p->a[1]*p->aPres[1])+(p->a[2]*p->aPres[2])),sqrt(a2),sqrt(ap2),
       sqrt(a2*0.25*p->fBall2)/(p->c*p->c),
       (p->a[0]*p->gradrho[0]+p->a[1]*p->gradrho[1]+p->a[2]*p->gradrho[2])/
       (p->gradrho[0]*p->gradrho[0]+p->gradrho[1]*p->gradrho[1]+p->gradrho[2]*p->gradrho[2]+
	0.4/p->fBall2)/(p->c*p->c),p->BalsaraSwitch );
			}
#endif

			/* Rarefaction or Gravity dominated compression */
			if ( UDOT_HYDRO(p) < 0 || adotap < 0.5*a2) p->ShockTracker = 0;
			else {
			 if (dShockTrackerB == 0) {
			  if (a2 < ap2) 
			    dv2 = sqrt(a2*0.25*p->fBall2);
			  else
			    dv2 = sqrt(ap2*0.25*p->fBall2);
			 }
			 else {
			  if (a2 < ap2) 
			    dv2 = -p->fDensity*(p->a[0]*p->gradrho[0]+p->a[1]*p->gradrho[1]+p->a[2]*p->gradrho[2])/
			      (p->gradrho[0]*p->gradrho[0]+p->gradrho[1]*p->gradrho[1]+p->gradrho[2]*p->gradrho[2]+
			       p->fDensity*p->fDensity*dShockTrackerB/p->fBall2);
			  else
			    dv2 = -p->fDensity*(p->aPres[0]*p->gradrho[0]+p->aPres[1]*p->gradrho[1]+p->aPres[2]*p->gradrho[2])/
			      (p->gradrho[0]*p->gradrho[0]+p->gradrho[1]*p->gradrho[1]+p->gradrho[2]*p->gradrho[2]+
			       p->fDensity*p->fDensity*dShockTrackerB/p->fBall2);

			  }
			 /*
			  if (dv2 < dShockTrackerA*p->c*p->c) p->ShockTracker = 0;
			 */
			 p->ShockTracker = dv2/(dShockTrackerA*p->c*p->c);
			 if (p->ShockTracker > 1) p->ShockTracker = 1;
			 }
			}
		}

#endif
        }

double pkdPoverRhoFloorJeansParticle(PKD pkd, double dResolveJeans, PARTICLE *p) 
    {
    double e2,l2;
    /*
     * Add pressure floor to keep Jeans Mass
     * resolved.  In comparison with Agertz et
     * al. 2009, dResolveJeans should be 3.0:
     * P_min = 3*G*max(h,eps)^2*rho^2
     * Note that G = 1 in our code
     */
#ifdef JEANSSOFTONLY
    l2 = p->fSoft*p->fSoft;
#else
    l2 = 0.25*p->fBall2; 
#ifdef JEANSSOFT
    e2 = p->fSoft*p->fSoft; 
    if (l2 < e2) l2 = e2; /* Jeans scale can't be smaller than softening */
#endif
#endif
    return l2*dResolveJeans*p->fDensity;
    }

void pkdGasPressureParticle(PKD pkd, struct GasPressureContext *pgpc, PARTICLE *p, 
    double *pPoverRhoFloorJeans, double *pPoverRhoNoncool, double *pPoverRhoGas, double *pcGas ) 
    {
#ifdef PCONST
    p->u = PCONST/(pgpc->gammam1*p->fDensity);
    p->uPred = p->u;
#endif
#ifdef MASSNONCOOL
    {
    double frac = p->fMassNoncool/p->fMass;
    /* Note: assuming that P/rho = (gamma-1) u (e.g. cooling_metal)
       some non-standard cooling may assume otherwise */
    *pPoverRhoGas = pgpc->gammam1*(p->uNoncoolPred*frac + p->uPred*(1-frac));
    *pPoverRhoNoncool = 0;
    *pcGas = sqrt(pgpc->gamma*(*pPoverRhoGas));
    }
#else /* !MASSNONCOOL */
#ifndef NOCOOLING
    if (pgpc->iGasModel == 2) {
        COOL *cl = pkd->Cool;

        CoolCodePressureOnDensitySoundSpeed( cl, &p->CoolParticle, p->uPred, p->fDensity, pgpc->gamma, pgpc->gammam1, pPoverRhoGas, pcGas );
        }
    else
#endif
    {
    *pPoverRhoGas = pgpc->gammam1*p->uPred;
    *pcGas = sqrt(pgpc->gamma*(*pPoverRhoGas));
    }

#ifdef UNONCOOL
    *pPoverRhoNoncool = (GAMMA_NONCOOL-1)*p->uNoncoolPred;
#else
    *pPoverRhoNoncool = 0;
#endif
#endif /* !MASSNONCOOL */

    *pPoverRhoFloorJeans = pkdPoverRhoFloorJeansParticle(pkd, pgpc->dResolveJeans, p);
    }

void  pkdSetThermalCond(PKD pkd, struct GasPressureContext *pgpc, PARTICLE *p) 
    {
#ifdef THERMALCOND
    double fThermalCond = pgpc->dThermalCondCoeffCode*pow(p->uPred,2.5); /* flux = coeff grad u   coeff ~ flux x h/u */ 
    double Tp = CoolCodeEnergyToTemperature(pkd->Cool, &p->CoolParticle, p->uPred, p->fMetals );
	if (Tp < pgpc->dEvapMinTemp) fThermalCond = 0;
    double fThermalCond2 = pgpc->dThermalCond2CoeffCode*pow(p->uPred,0.5);
	if (Tp < pgpc->dEvapMinTemp) fThermalCond2 = 0;
    double fSat = p->fDensity*p->c*p->fThermalLength; /* Max flux x L/u */
    double fThermalCondSat = pgpc->dThermalCondSatCoeff*fSat;
    double fThermalCond2Sat = pgpc->dThermalCond2SatCoeff*fSat;

//    printf("Saturated %d %g %g %g %g %g %g %g %g %g\n",p->iOrder,p->r[0],p->fDensity,p->uPred/4.80258,fThermalCond,fThermalCond2,fThermalCondSat,fThermalCond2Sat,p->fThermalLength,sqrt(p->fBall2*0.25));
    p->fThermalCond = (fThermalCond < fThermalCondSat ? fThermalCond : fThermalCondSat) +
        (fThermalCond2 < fThermalCond2Sat ? fThermalCond2 : fThermalCond2Sat);
#endif
    }


/* Note: Uses uPred */
void pkdGasPressure(PKD pkd, struct GasPressureContext *pgpc)
{
    PARTICLE *p;
    int i;

    p = pkd->pStore;
    for(i=0;i<pkdLocal(pkd);++i,++p) {
		if (pkdIsGas(pkd,p)) {
            double PoverRho,PoverRhoGas,PoverRhoNoncool,PoverRhoMinJeans,cGas,PoverRhoJeans;

            pkdGasPressureParticle(pkd, pgpc, p, &PoverRhoMinJeans, &PoverRhoNoncool, &PoverRhoGas, &cGas ); 
            PoverRho = PoverRhoGas + PoverRhoNoncool;
            PoverRhoJeans = (PoverRho < PoverRhoMinJeans ? PoverRhoMinJeans - PoverRho : 0);
            PoverRho += PoverRhoJeans;
            p->PoverRho2 = PoverRho/p->fDensity;
            p->c = sqrt(cGas*cGas+GAMMA_NONCOOL*PoverRhoNoncool+GAMMA_JEANS*PoverRhoJeans);
#ifdef THERMALCOND
            pkdSetThermalCond(pkd,pgpc,p);
#endif
#ifdef DTADJUST
                {
                double uTotDot, dt;
                uTotDot = p->uDot;
#ifdef UNONCOOL
                uTotDot += p->uNoncoolDot;
#endif
                if (uTotDot > 0) dt = pgpc->dtFacCourant*sqrt(p->fBall2*0.25/(4*(p->c*p->c+GAMMA_NONCOOL*uTotDot*p->dt)));
                else dt = pgpc->dtFacCourant*sqrt(p->fBall2*0.25)/(2*(p->c));
                if (dt < p->dt) p->dt = dt; // Update to scare the neighbours
                }
#endif
            }
#ifdef DEBUG
        if (pkdIsGas(pkd,p) && (p->iOrder % 1000)==0) {
			printf("Pressure %i: %i %i %f %f %f  %f %f %f %f %f\n",
			       p->iOrder,TYPEQueryACTIVE(p),TYPEQueryTREEACTIVE(p),
			       p->r[0],p->r[1],p->r[2],sqrt(0.25*p->fBall2),p->fDensity,p->uPred,
			       p->PoverRho2*p->fDensity*p->fDensity,p->c);
			}
#endif            
		}
    }

void pkdGetDensityU(PKD pkd, double uMin)
{
#ifdef DENSITYU
    PARTICLE *p;
    int i;

    p = pkd->pStore;
    for(i=0;i<pkdLocal(pkd);++i,++p) {
		if (pkdIsGas(pkd,p) && TYPEQueryACTIVE(p)) {
   			p->fDensityU /= (p->uPred + uMin);
		}
	}
#endif
}

void pkdLowerSoundSpeed(PKD pkd, double dhMinOverSoft)
{
    PARTICLE *p;
    double dfBall2MinOverSoft,dfBall2Min,ratio;
    int i;

    dfBall2MinOverSoft = 4*dhMinOverSoft*dhMinOverSoft;

    p = pkd->pStore;
    for(i=0;i<pkdLocal(pkd);++i,++p) {
		if (pkdIsGas(pkd,p)) {
		        dfBall2Min = dfBall2MinOverSoft*p->fSoft*p->fSoft;
			if (p->fBall2 < dfBall2Min) {
			  ratio = p->fBall2/dfBall2Min;
			  p->PoverRho2 *= ratio;
			  p->c *= sqrt(ratio);
#ifdef DEBUG
		if ((p->iOrder % 100)==0) {
			printf("Pressure %i: %g %g, %g %g, %g %g\n",
			       p->iOrder,
			       sqrt(p->fBall2),sqrt(dfBall2Min),
			       p->PoverRho2,p->PoverRho2/ratio,
			       p->c,p->c/sqrt(ratio));
			}
#endif            
			  }
			}
		}
    }

void pkdInitEnergy(PKD pkd, double dTuFac, double z, double dTime )
    {
    PARTICLE *p;
    int i;
#ifndef NOCOOLING
    COOL *cl;
    double T,E;
    
    cl = pkd->Cool;
    CoolSetTime( cl, dTime, z  );
#endif
    
    p = pkd->pStore;
    for(i=0;i<pkdLocal(pkd);++i,++p) {
        if (TYPEQueryTREEACTIVE(p) && pkdIsGas(pkd,p)) {
#ifndef NOCOOLING
			T = p->u / dTuFac;
#ifdef COOLDEBUG
			cl->p = p; /* Send in particle pointer only for temporary debug */
#endif
			CoolInitEnergyAndParticleData( cl, &p->CoolParticle, &E, p->fDensity, T, p->fMetals );
			p->u = E;
#endif
            p->uPred = p->u;
#ifdef DEBUG
            if ((p->iOrder % 1000)==0) {
                printf("InitEnergy %i: %f %g   %f %f %f %g\n",
                    p->iOrder,T,p->u * cl->dErgPerGmUnit,
                    p->CoolParticle.HI,p->CoolParticle.HeI,p->CoolParticle.HeII,p->fDensity*cl->dComovingGmPerCcUnit);
                } 
#endif            
            }
        }
    }

#ifdef GLASS
/* Currently wired to have no more than two regions with different
   Pressures (densities) split by x=0 with a linear connection */
void pkdGlassGasPressure(PKD pkd, void *vin)
{
    PARTICLE *p;
    double PoverRho,xx,nsp=2.5;
    int i;
    struct GlassData *in = vin;

    for(i=0;i<pkdLocal(pkd);++i) {
                p = &pkd->pStore[i];
                if (TYPEQueryTREEACTIVE(p)) {
#define NEWGLASS
#ifdef NEWGLASS
		    PoverRho = in->dGlassPoverRhoL*(p->r[0]*p->r[0]+p->r[1]*p->r[1]+p->r[2]*p->r[2])+1/p->fDensity;
#else
    		        if (p->r[0] < -nsp*in->dGlassxL) {
			  if (p->r[0] > in->dxBoundL + nsp*in->dGlassxL)
			       PoverRho=in->dGlassPoverRhoL;  
			  else {
			       xx = ( p->r[0] - in->dxBoundL + nsp*in->dGlassxR )
				 / ( nsp*in->dGlassxL+nsp*in->dGlassxR );
			       xx = xx*xx* ( -2*xx + 3 );
		               PoverRho = in->dGlassPoverRhoR + xx
				 *(in->dGlassPoverRhoL - in->dGlassPoverRhoR);
              		       }
			  }
                        else if (p->r[0] > nsp*in->dGlassxR) {
			  if (p->r[0] < in->dxBoundR - nsp*in->dGlassxR)
			       PoverRho=in->dGlassPoverRhoR;  
			  else {
			       xx = ( p->r[0] - in->dxBoundR + nsp*in->dGlassxR )
				 / ( nsp*in->dGlassxL+nsp*in->dGlassxR );
			       xx = xx*xx* ( -2*xx + 3 );
		               PoverRho = in->dGlassPoverRhoR + xx
				 *(in->dGlassPoverRhoL - in->dGlassPoverRhoR);
			       }
			  }
			else {
			       xx = ( p->r[0] + nsp*in->dGlassxL )
				 / ( nsp*in->dGlassxL+nsp*in->dGlassxR );
			       xx = xx*xx* ( -2*xx + 3 );
		               PoverRho = in->dGlassPoverRhoL + xx
				 *(in->dGlassPoverRhoR - in->dGlassPoverRhoL);
			       }
#endif

			p->u = PoverRho;
			p->uPred = PoverRho;
			p->PoverRho2 = PoverRho/p->fDensity;
   			p->c = sqrt(in->dGamma*PoverRho);
		        }
#if (0)
		if (pkdIsGas(pkd,p) && (p->iOrder % 1000)==0) 
		        printf("Glass P %i: %i %i %f %f %f  %f %f %f %f %f\n",
			       p->iOrder,TYPEQueryACTIVE(p),TYPEQueryTREEACTIVE(p),
			       p->r[0],p->r[1],p->r[2],sqrt(0.25*p->fBall2),p->fDensity,p->uPred,p->PoverRho2*p->fDensity*p->fDensity,p->c);
#endif
                }
    }

#endif

void pkdKickRhopred(PKD pkd, double dHubbFac, double dDelta)
{
	PARTICLE *p;
	int i,n;

#ifdef DEBUG
	printf("pkdKickRhopred: %g %g\n",dHubbFac,dDelta);
#endif
	p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;++i) {
#ifdef DEBUG
		if (pkdIsGas(pkd,(p+i)) && ((p+i)->iOrder % 3000)==0) {
            printf("Rhopreding %i: %i %i %f %f %f %f   %f %f %f %f\n",
                (p+i)->iOrder,TYPEQueryACTIVE(p+i),
                TYPEQueryTREEACTIVE(p+i),
                sqrt(0.25*(p+i)->fBall2),(p+i)->fDensity,(p+i)->u,(p+i)->uPred,
                (p+i)->a[0],(p+i)->a[1],(p+i)->a[2],UDOT_HYDRO(p+i));
			}
#endif
		if(TYPEFilter( &p[i], TYPE_GAS|TYPE_ACTIVE, TYPE_GAS )) {
			p[i].fDensity = p[i].fDensity*(1 + dDelta*(dHubbFac - p[i].divv));
			}
		}
	}

int pkdSphCurrRung(PKD pkd, int iRung, int bGreater)
{
    int i;
    int iCurrent;
    
    iCurrent = 0;
    for(i = 0; i < pkdLocal(pkd); ++i) {
		if(pkdIsGas(pkd, &pkd->pStore[i]) &&
		   (pkd->pStore[i].iRung == iRung ||
			(pkd->pStore[i].iRung > iRung && bGreater))) {
			iCurrent = 1;
			break;
			}
		}
    return iCurrent;
    }

double pkdDtFacCourant( double dEtaCourant, double dCosmoFac ) {
    return dEtaCourant*dCosmoFac*2/1.6;
    }

/* DTTEST should be defined to the dt value that going under will trigger the print */
#ifdef DTTEST
#define DTSAVE(_val,_label)  { \
    if ((_label)[0]=='0') dTnSave=0; \
    dTSave[dTnSave]=_val; \
    strncpy(&dTLabel[dTnSave][0], _label, 3); \
    dTLabel[dTnSave][3]='\0'; \
    dTnSave++; \
    }
#else
#define DTSAVE(_val,_label)
#endif

void
pkdSphStep(PKD pkd, double dCosmoFac, double dEtaCourant, double dEtauDot, double dDiffCoeff, double dEtaDiffusion, double dResolveJeans, int bViscosityLimitdt, double *pdtMinGas)
    {
    int i;
    PARTICLE *p;    
    double dT,dTC,dTu,dTD,ph;
#ifdef DTTEST
    double dTSave[12],dtCut=DTTEST;
    char dTLabel[12][4];
    int dTnSave=0,nFail=0;
#endif

    *pdtMinGas = DBL_MAX;
    for(i=0;i<pkdLocal(pkd);++i) {
        p = &pkd->pStore[i];
        if(pkdIsGas(pkd, p)) {

            if (TYPEQueryACTIVE(p)) {
                ph = sqrt(0.25*p->fBall2);
#ifdef SINKING
                if (TYPETest( p, TYPE_SINKING)) {
                    p->dt = FLT_MAX; /* reset later to min gas step */
                    continue;
                    }
#endif
                /*
                 * Courant condition goes here.
                 */
                DTSAVE(p->dt,"0IN");
#if defined(DRHODT) || defined(DTADJUST)
                dT = p->dtNew; /* Start with estimate from prior Sph force calculations */ 
                DTSAVE(dT,"SPH");
                p->dtNew = FLT_MAX;
#else /* If not doing DTADJUST */
                if (p->mumax>0.0) {
                    if (bViscosityLimitdt) 
                        dTC = dEtaCourant*dCosmoFac*(ph/(p->c + 0.6*(p->c + 2*p->BalsaraSwitch*p->mumax)));
                    else
                        dTC = dEtaCourant*dCosmoFac*(ph/(p->c + 0.6*(p->c + 2*p->mumax)));
                    }
                else {
#if defined(PRES_HK) || defined(PRES_MONAGHAN) || defined(SIMPLESF)
                    dTC = dEtaCourant*dCosmoFac*(ph/(1.6*p->c+(10./FLT_MAX)));
#else
                    dTC = dEtaCourant*dCosmoFac*(ph/(1.6*p->c));
#endif
                    }
                DTSAVE(dTC,"COU");
                dT = dTC;
#endif

#ifdef DTADJUST
                    {
                    double uTotDot, dtExtrap;
#ifdef MASSNONCOOL
					double x = p->fMassNoncool/p->fMass;
					uTotDot = p->uNoncoolDot*(1-x)+p->uDot*x;
#else
                    uTotDot = p->uDot;
#ifdef UNONCOOL
                    uTotDot += p->uNoncoolDot;
#endif /* MASSNONCOOL */
#endif
                    if (uTotDot > 0) {
                        dtExtrap = pkdDtFacCourant(dEtaCourant,dCosmoFac)
                            *sqrt(p->fBall2*0.25/(4*(p->c*p->c+GAMMA_NONCOOL*uTotDot*p->dt)));
                        DTSAVE(dtExtrap,"UEX");
                        if (dtExtrap < dT) dT = dtExtrap; 
                        }
                    }
#endif
                if (dEtauDot > 0.0 && p->uDotPdV < 0.0) { /* Prevent rapid adiabatic cooling */
                    double PoverRhoFloorJeans=pkdPoverRhoFloorJeansParticle(pkd, dResolveJeans, p);
                    double uEff = PONRHOFLOOR+PoverRhoFloorJeans/(GAMMA_JEANS-1)+p->u;
                            
                    assert(p->u > 0.0);
#ifdef UNONCOOL
                    uEff += p->uNoncool;
#endif
                    dTu = dEtauDot*uEff/fabs(p->uDotPdV);
                    DTSAVE(dTu,"PDV");
                    if (dTu < dT) dT = dTu;
                    }
#ifdef DIFFUSION
#ifdef THERMALCOND
            /* h^2/(2.77Q) Linear stability from Brookshaw */
                if (p->fThermalCond > 0 || (p->diff > 0 && dDiffCoeff > 0)) {
                    dTD = dEtaDiffusion*ph*ph*dCosmoFac*dCosmoFac
                        /(dDiffCoeff*p->diff + ph/p->fThermalLength*p->fThermalCond/p->fDensity);  
                    DTSAVE(dTD,"DIF");
                    if (dTD < dT) dT = dTD;
                    }
#else
            /* h^2/(2.77Q) Linear stability from Brookshaw */
                if (p->diff > 0 && dDiffCoeff > 0) {
                    dTD = dEtaDiffusion*ph*ph*dCosmoFac*dCosmoFac
                        /(dDiffCoeff*p->diff);  
                    DTSAVE(dTD,"DIF");
                    if (dTD < dT) dT = dTD;
                    }
#endif
#endif
#ifdef DTTEST                
                if ((dT < dtCut && nFail < 10)) {
                    int j;
                    double T;
                    nFail++;
                    dtCut=dT; /* Try to print extreme low dt's */
                    fprintf(stderr,"dt problem %d: dt %g < %g",p->iOrder,dT,DTTEST);
                    for (j=0;j<dTnSave;j++) fprintf(stderr,", %3s %g",&dTLabel[j][0],dTSave[j]);
                    fprintf(stderr,"\n");
#ifndef NOCOOLING
                    T = CoolCodeEnergyToTemperature( pkd->Cool, &p->CoolParticle, p->uPred, p->fMetals );
#endif
#ifdef THERMALCOND
                    fprintf(stderr,"u %g T %g %g c %g h %g divv %g rho %g Z %g dtdiff %g %g\n",p->uPred,p->uPred/4802.57,T,p->c,sqrt(0.25*p->fBall2),p->divv,p->fDensity,p->fMetals,p->uPred/(fabs(p->uDotDiff)+1e-20),p->fThermalCond);
#else
                    fprintf(stderr,"u %g T %g %g c %g h %g divv %g rho %g Z %g dtdiff %g \n",p->uPred,p->uPred/4802.57,T,p->c,sqrt(0.25*p->fBall2),p->divv,p->fDensity,p->fMetals,p->uPred/(fabs(p->uDotDiff)+1e-20));
#endif
                    }
#endif
                if(dT < p->dt) p->dt = dT;
                }
            /* This code relies on SPH step being done last -- not good */
	        if (p->dt < *pdtMinGas) { *pdtMinGas = p->dt; }
            }
#ifdef SINKING
        else if (TYPETest( p, TYPE_SINK )) {
            if (p->dt < *pdtMinGas) { *pdtMinGas = p->dt; }
            }
#endif
        }
    }

void
pkdSinkStep(PKD pkd, double dtMax)
{
    int i;
    PARTICLE *p;    

    for(i=0;i<pkdLocal(pkd);++i) {
        p = &pkd->pStore[i];
        if(TYPETest( p, TYPE_SINK|TYPE_SINKING ) && TYPEQueryACTIVE(p)) {
	    if (dtMax < p->dt) p->dt = dtMax;
	    }
	}
}

void
pkdSetSphStep(PKD pkd, double dt)
{
    int i;
    PARTICLE *p;    
	int n=0;

    for(i=0;i<pkdLocal(pkd);++i) {
        p = &pkd->pStore[i];
        if(TYPETest( p, TYPE_GAS )) p->dt = dt;
	}
}

void 
pkdSphViscosityLimiter(PKD pkd, int bOn, int bShockTracker)
{
    int i;
    PARTICLE *p;    

    if (bOn) {
        for(i=0;i<pkdLocal(pkd);++i) {
			p = &pkd->pStore[i];
                        /* Only set values for particles with fresh curlv, divv from smooth */
			if(pkdIsGas(pkd, p) && TYPETest(p,TYPE_ACTIVE)) {
				if (p->divv!=0.0) {         	 
					p->BalsaraSwitch = fabs(p->divv)/
						(fabs(p->divv)+sqrt(p->curlv[0]*p->curlv[0]+
											p->curlv[1]*p->curlv[1]+
											p->curlv[2]*p->curlv[2]));
					}
				else { 
					p->BalsaraSwitch = 0;
					}
				}
			}
        }
    else {
        for(i=0;i<pkdLocal(pkd);++i) {
			p = &pkd->pStore[i];
			if(pkdIsGas(pkd, p)) {
			        p->BalsaraSwitch = 1;
				}
			}
        }
    }

void pkdDensCheck(PKD pkd, int iRung, int bGreater, int iMeasure, void *data) {
#if (0)
    int i;
    struct {     
		double dMaxDensError;
		double dAvgDensError;
		int nError;
		int nTotal;
    } *tmp=data;
    double error;

    char ach[256];

    tmp->dMaxDensError=0;
    tmp->dAvgDensError=0;
    tmp->nError=0;
    tmp->nTotal=0;

    if (!iMeasure) {
        for(i = 0; i < pkdLocal(pkd); ++i) {
			if(TYPETest(&(pkd->pStore[i]),TYPE_GAS) &&
			   (pkd->pStore[i].iRung == iRung ||
				(bGreater && pkd->pStore[i].iRung > iRung))) {
				if (pkd->pStore[i].fDensity == 0) {
					sprintf(ach, "dens zero i: %d dens %g iAct %d\n",
							pkd->pStore[i].iOrder,pkd->pStore[i].fDensity,
							pkd->pStore[i].iActive);
					mdlDiag(pkd->mdl, ach);
					}
				pkd->pStore[i].fDensSave = pkd->pStore[i].fDensity;
				}
			}
		return;
        }

    for(i=0;i<pkdLocal(pkd);++i) {
		if(TYPETest(&(pkd->pStore[i]),TYPE_GAS) &&
		   (pkd->pStore[i].iRung == iRung ||
			(bGreater && pkd->pStore[i].iRung > iRung))) {
			error = abs((pkd->pStore[i].fDensSave - pkd->pStore[i].fDensity)/pkd->pStore[i].fDensity);
			tmp->dAvgDensError += error;
			tmp->nTotal++;
			if (error>tmp->dMaxDensError) tmp->dMaxDensError=error;
			if (error>1e-5) {
				tmp->nError++;
				sprintf(ach, "dens error i: %d save %g dens %g  iAct %d\n",
						pkd->pStore[i].iOrder,pkd->pStore[i].fDensSave,
						pkd->pStore[i].fDensity,pkd->pStore[i].iActive);
				mdlDiag(pkd->mdl, ach);
				}
			}
		}

    tmp->dAvgDensError/=tmp->nTotal; 
    return;
#endif
	}

#endif /* GASOLINE */

#ifdef GLASS

void
pkdRandomVelocities(PKD pkd, double dMaxVL, double dMaxVR)
{
    int i,j;
    PARTICLE *p;  
    double v;

    for(i=0;i<pkdLocal(pkd);++i) {
        p = &pkd->pStore[i];
        if (p->r[0]<0.0) {
			for (j=0;j<3;j++) {
				p->v[j]+=(v=((random()%20001)/10000.0-1.0)*dMaxVL);
				p->vPred[j]+=v;
				}
			}
        else {
			for (j=0;j<3;j++) {
				p->v[j]+=(v=((random()%20001)/10000.0-1.0)*dMaxVR);
				p->vPred[j]+=v;
				}
			}
		}
    }

#endif /* GLASS */

#ifdef COLLISIONS

int
pkdNumRejects(PKD pkd)
{
	int i,nRej = 0;

	for (i=0;i<pkd->nLocal;i++)
		if (REJECT(&pkd->pStore[i])) ++nRej;

	return nRej;
	}

void
pkdReadSS(PKD pkd,char *pszFileName,int nStart,int nLocal)
{
	SSIO ssio;
	SSDATA data;
	PARTICLE *p;
	int i,j, iSetMask;

	pkd->nLocal = nLocal;
	pkd->nActive = nLocal;
	/*
	 ** General initialization (modeled after pkdReadTipsy()).
	 */
	for (i=0;i<nLocal;++i) {
		p = &pkd->pStore[i];
		TYPEClear(p);
		p->iRung = 0;
		p->fWeight = 1.0;
		p->fDensity = 0.0;
		p->fBall2 = 0.0;
		p->fBallMax = 0.0;
		p->iDriftType = NORMAL; /*DEBUG must initialize b/c pkdDrift()...*/
#ifdef SAND_PILE
		p->bStuck = 0;
#endif
		}
	/*
	 ** Seek past the header and up to nStart.
	 */
	if (ssioOpen(pszFileName,&ssio,SSIO_READ))
		mdlassert(pkd->mdl,0); /* unable to open ss file */
	if (ssioSetPos(&ssio,SSHEAD_SIZE + nStart*SSDATA_SIZE))
		mdlassert(pkd->mdl,0); /* unable to seek in ss file */
	/*
	 ** Read Stuff!
	 */
	for (i=0;i<nLocal;++i) {
		p = &pkd->pStore[i];
		p->iOrder = nStart + i;
		if (!pkdIsDarkByOrder(pkd,p)) /* determined by p->iOrder */
			mdlassert(pkd->mdl,0); /* only dark particles allowed in ss file */
		iSetMask = TYPE_DARK;
		if (ssioData(&ssio,&data))
			mdlassert(pkd->mdl,0); /* error during read in ss file */

		/* This is changed because the data.org_idx appears to always be zero. DCR made a change that has probably not been incorporated yet. This will work as long as pkdReadSS is only called from msrOneNodeReadSS. This is the original line of code:
		p->iOrgIdx = data.org_idx;
		*/
		p->iOrgIdx = p->iOrder;
		p->fMass = data.mass;
		p->fSoft = SOFT_FROM_SSDATA(&data); /* half physical radius */
#ifdef CHANGESOFT 
 		p->fSoft0 = p->fSoft;
#endif
		for (j=0;j<3;++j) p->r[j] = data.pos[j];
		for (j=0;j<3;++j) p->v[j] = data.vel[j];
		for (j=0;j<3;++j) p->w[j] = data.spin[j];
		p->iColor = data.color;
#ifdef NEED_VPRED
		for (j=0;j<3;++j) p->vPred[j] = p->v[j];
#endif
		TYPESet(p,iSetMask);
		}
	if (ssioClose(&ssio))
		mdlassert(pkd->mdl,0); /* unable to close ss file */
	}

void
pkdWriteSS(PKD pkd,char *pszFileName,int nStart)
{
	SSIO ssio;
	SSDATA data;
	PARTICLE *p;
	int i,j;

	/*
	 ** Seek past the header and up to nStart.
	 */
	if (ssioOpen(pszFileName,&ssio,SSIO_UPDATE))
		mdlassert(pkd->mdl,0); /* unable to open ss file */
	if (ssioSetPos(&ssio,SSHEAD_SIZE + nStart*SSDATA_SIZE))
		mdlassert(pkd->mdl,0); /* unable to seek in ss file */
	/* 
	 ** Write Stuff!
	 */
	for (i=0;i<pkdLocal(pkd);++i) {
		p = &pkd->pStore[i];
		if (!pkdIsDark(pkd,p))
			mdlassert(pkd->mdl,0); /* only dark particles allowed in ss file */
		data.org_idx = p->iOrgIdx;
		data.mass = p->fMass;
		data.radius = RADIUS(p); /*DEBUG note CHANGESOFT not supported (what is it, anyway?)*/
		for (j=0;j<3;++j) data.pos[j]  = p->r[j];
		for (j=0;j<3;++j) data.vel[j]  = p->v[j];
		for (j=0;j<3;++j) data.spin[j] = p->w[j];
		data.color = p->iColor;
		if (ssioData(&ssio,&data))
			mdlassert(pkd->mdl,0); /* unable to write in ss file */
		}
	if (ssioClose(&ssio))
		mdlassert(pkd->mdl,0); /* unable to close ss file */
	}

void
pkdKickUnifGrav(PKD pkd,double dvx,double dvy,double dvz)
{
	PARTICLE *p;
	int i,n;

	p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;i++)
#ifdef SAND_PILE
		if (TYPEQueryACTIVE(&p[i]) && !p[i].bStuck) {
#else
		if (TYPEQueryACTIVE(&p[i])) {
#endif
			p[i].v[0] += dvx;
			p[i].v[1] += dvy;
			p[i].v[2] += dvz;
			}
	}

void
pkdNextEncounter(PKD pkd,double *dt)
{
	PARTICLE *p;
	int i;

	for (i=0;i<pkdLocal(pkd);++i) {
		p = &pkd->pStore[i];
		if (p->dtCol < *dt) *dt = p->dtCol;
		}
	}

void
pkdMarkEncounters(PKD pkd,double dt)
{
	PARTICLE *p;
	int i;

	if (dt <= 0) { /* initialize */
		for (i=0;i<pkdLocal(pkd);++i) {
			p = &pkd->pStore[i];
			TYPESet(p,TYPE_ACTIVE);/*DEBUG redundant?*/
			p->iDriftType = KEPLER;
			p->dtCol = DBL_MAX;
/*DEBUG old neighbour stuff		p->idNbr.iPid = p->idNbr.iIndex =
				p->idNbr.iOrder = -1; */
			}
		}
	else {
		for (i=0;i<pkdLocal(pkd);++i) {
			p = &pkd->pStore[i];
			if (p->dtCol < dt) {
				TYPESet(p,TYPE_ACTIVE);
				p->iDriftType = NORMAL;
				}
			else {
				TYPEReset(p,TYPE_ACTIVE);
				p->iDriftType = KEPLER;
				}
			}
		}
	}

#endif /* COLLISIONS */

#ifdef SLIDING_PATCH

void pkdPatch(PKD pkd)
{
  PARTICLE *p;
  PATCH_PARAMS *PP;
  int i;

  PP = pkd->PP;

  for (i=0;i<pkd->nLocal;i++) {
	p = &pkd->pStore[i];
	if (!TYPEQueryACTIVE(p)) continue;
	/*
	** Apply the velocity independendent part of Hill's Equations
	** (Kick2 in the symplectic integration paper)
	*/
	p->a[0] -= PP->dOrbFreq*PP->dOrbFreq*p->r[0];
	p->a[2] -= PP->dOrbFreqZ2*p->r[2];
  }

  if (PP->bExtPert) {
	double rPert[3];      /* position of perturber relative to saturn */
	double rpp[3];        /* position of patch relative to perturber */
	double R[3];          /* position of a patch particle relative to perturber (goal) */
	double Rmag,Rmag3inv; /* magnitude of R, and inverse cube */

	/*
	** NOTE: positive x: away from saturn, y: direction of orbital velocity, z: RH rule
	** NOTE: All vectors are calculated in the *patch* frame (rotating wrt saturn)
	*/

	/* calculate position of perturber relative to saturn, in the patch frame,
	** by using a rotation matrix:
    **       cos(omega_patch * t)        sin(omega_patch * t)
    **      -sin(omega_patch * t)        cos(omega_patch * t)
	*/

	rPert[0] = (PP->dPertOrbDist*cos((PP->dPertOrbFreq*pkd->dTime) + PP->dPertPhase)*cos(PP->dOrbFreq*pkd->dTime))
	  + (PP->dPertOrbDist*sin((PP->dPertOrbFreq*pkd->dTime) + PP->dPertPhase)*sin(PP->dOrbFreq*pkd->dTime));

	rPert[1] = (PP->dPertOrbDist*cos((PP->dPertOrbFreq*pkd->dTime) + PP->dPertPhase)*(-sin(PP->dOrbFreq*pkd->dTime)))
	  + (PP->dPertOrbDist*sin((PP->dPertOrbFreq*pkd->dTime) + PP->dPertPhase)*cos(PP->dOrbFreq*pkd->dTime));

	rPert[2] = PP->dPertMaxZ*cos((PP->dPertOrbFreqZ*pkd->dTime) + PP->dPertPhaseZ);

	/* calculate vector pointing from perturber (center) to patch center */
	rpp[0] = PP->dOrbDist - rPert[0];
	rpp[1] = -rPert[1];
	rpp[2] = -rPert[2];

	for (i=0;i<pkd->nLocal;i++) {
	  p = &pkd->pStore[i];
	  if (!TYPEQueryACTIVE(p)) continue;

	  /* calculate force/mass of an external perturber on sliding-patch particles */

	  /* calculate position vector */
	  R[0] = rpp[0] + p->r[0];
	  R[1] = rpp[1] + p->r[1];
	  R[2] = rpp[2] + p->r[2];

	  Rmag = sqrt(R[0]*R[0] + R[1]*R[1] + R[2]*R[2]);
	  Rmag3inv = 1.0/(Rmag*Rmag*Rmag);

	  p->a[0] -= PP->dPertMass*R[0]*Rmag3inv;
	  p->a[1] -= PP->dPertMass*R[1]*Rmag3inv;
	  p->a[2] -= PP->dPertMass*R[2]*Rmag3inv;
	}
  }
}

double _pkdUniRan(void)
{
  return (double) random()/RAND_MAX;
}

double _pkdGauRan(void)
{
  /* generates Gaussian deviate, mean 0, std dev 1 */
  /* based on NRiC 7.2 gasdev() */

  static int iNeedDev = 1;
  static double dSecondDev;
  double x,y,r2,f;

  if (iNeedDev) {
	do {
	  x = 2.0*_pkdUniRan() - 1.0;
	  y = 2.0*_pkdUniRan() - 1.0;
	  r2 = x*x + y*y;
	} while (r2 >= 1.0 || r2 == 0.0);
	f = sqrt(-2.0*log(r2)/r2);
	dSecondDev = x*f;
	iNeedDev = 0;
	return y*f;
  } else {
	iNeedDev = 1;
	return dSecondDev;
  }
}

int pkdRandAzWrap(PKD pkd)
{
  /*DEBUG CURRENTLY IGNORING SPIN*/
  /*
  ** This routine probably only valid for low surface densities
  ** that do not lead to strong wake/aggregate formation.
  */

  const double Sqrt2OverPi = M_2_SQRTPI/M_SQRT2; /* sqrt(2/PI) */

  static int bFirstCall = 1;

  PATCH_PARAMS *PP;
  PARTICLE *p;
  double amp,phase;
  int i,iSide,nRand = 0;

  assert(pkd->nLocal > 0);

  if (bFirstCall) {
	/*
	** Seeding strategy: since nodes in a cluster may be well
	** sychronized, we mod the seed (PID) by the iOrder of the first
	** particle on the processor (plus the number of local particles
	** in case the iOrder is zero).
	*/
	unsigned int seed = (int) time(NULL) % getpid() + getppid();

	(void) printf("pkdRandAzWrap(): processor %i random seed = %u\n",pkd->idSelf,seed);
	srandom(seed);
	bFirstCall = 0;
  }

  PP = pkd->PP;

  for (i=0;i<pkd->nLocal;i++) {
	p = &pkd->pStore[i];
	TYPESet(p,TYPE_TREEACTIVE); /* all particles go into density tree */
	TYPEReset(p,TYPE_SMOOTHACTIVE); /* not all are centers for overlap searches */
	if (p->bAzWrap) {
	  /* assign uniform random radial position depending on side */
	  iSide = SGN(p->r[0]); /* -1 ==> left side of patch; +1 ==> right side */
	  if (iSide == 0) iSide = 2*(random()%2) - 1; /* handle rare case */
	  /* flip side if equilibrium stream on opposite side */
	  if ((iSide == -1 && PP->iStripOption == STRIP_RIGHT_ONLY) ||
		  (iSide ==  1 && PP->iStripOption == STRIP_LEFT_ONLY))
		iSide *= -1;
	  p->r[0] = iSide*(PP->dStripInner + (PP->dStripOuter - PP->dStripInner)*_pkdUniRan());
	  /* assign random velocity relative to shear */
	  p->v[0] = _pkdGauRan()*PP->dVelDispX*sqrt(PP->dAvgMass/p->fMass);
	  p->v[1] = _pkdGauRan()*PP->dVelDispY*sqrt(PP->dAvgMass/p->fMass);
	  p->v[1] -= 1.5*PP->dOrbFreq*p->r[0]; /* add shear */
	  /* deal with the third dimension */
	  amp = PP->dAvgVertAmp*Sqrt2OverPi*sqrt(-2.0*log(_pkdUniRan()))*sqrt(PP->dAvgMass/p->fMass); /* Rayleigh deviate */
	  phase = 2.0*M_PI*_pkdUniRan();
	  p->r[2] = amp*cos(phase);
	  p->v[2] = - amp*sin(phase);
	  /* set flag to do an overlap search around this particle */
	  TYPESet(p,TYPE_SMOOTHACTIVE);
	  ++nRand;
	  /*DEBUG*/
	  if ((iSide == -1 && p->r[1] > -0.4*PP->dLength) ||
		  (iSide == 1 && p->r[1] < 0.4*PP->dLength))
		(void) fprintf(stderr,"DEBUG iOrd=%i iSide=%i y=%g (Ly=%g)\n",p->iOrder,iSide,p->r[1],PP->dLength);
	}
  }

  return nRand;
}

#endif /* SLIDING_PATCH */

#ifdef SIMPLE_GAS_DRAG
void
_pkdSimpleGasVel(double r[],int iFlowOpt,double dTime,double u[])
{
	/* From Paolo Tanga's e-mail Feb 19, 2001: M, f_1 and f_2 are
	   parameters (double precision) that can be considered internal
	   to _pkdGasVel, since they define the streamfunction. It define
	   the attached velocity field. f_1 and f_2 are scaling parameters
	   that allow to control the size of the vortex and its
	   width/length ratio.  Note that I have centered a vortex in 0,0,
	   which should be good in the local coordinates of a single
	   patch. Anyway, the given function is not suitable for the
	   simulation of an entire disk...! */

	switch (iFlowOpt) {
	case 1: { /* PATCH_FLOW */
		const double M	= 1.0e-7;
		const double f1	= 1.0e7;
		const double f2	= 1.0e7;
#ifndef SLIDING_PATCH
		assert(0); /* this formula only valid in patch frame */
#endif

		u[0] = M*cos(f2*r[1] - M_PI_2)*cos(f1*r[0] - M_PI);
		u[1] = M*sin(f2*r[1] - M_PI_2)*sin(f1*r[0] - M_PI);
		u[2] = 0.0;
		break;
		}
	default:
		assert(0); /* no other cases defined */
		}
}

void pkdSimpleGasDrag(PKD pkd,int iFlowOpt,int bEpstein,double dGamma,
					  double dTime)
{
	PARTICLE *p;
	double u[3],dCoef;
	int i;

	for (i=0;i<pkd->nLocal;i++) {
		p = &pkd->pStore[i];
		if (!TYPEQueryACTIVE(p)) continue;
		_pkdSimpleGasVel(p->r,iFlowOpt,dTime,u); /* get gas velocity */
		if (bEpstein)
			dCoef = -dGamma/RADIUS(p); /* i.e. -gamma/R */
		else
			dCoef = -dGamma/(RADIUS(p)*RADIUS(p)); /* i.e. -gamma/R^2 */
		/* note: following assumes dust density << gas density */
		p->a[0] += dCoef*(p->vPred[0] - u[0]);
		p->a[1] += dCoef*(p->vPred[1] - u[1]);
		p->a[2] += dCoef*(p->vPred[2] - u[2]);
		}
	}
#endif /* SIMPLE_GAS_DRAG */

#ifdef NEED_VPRED
#ifdef GASOLINE
void
pkdKickVpred(PKD pkd,double dvFacOne,double dvFacTwo,double duDelta,
    int iGasModel,double z,double duDotLimit, double dTimeEnd,UNCC uncc)
{
	PARTICLE *p;
	int i,j,n;

	mdlDiag(pkd->mdl, "Into Vpred\n");

	pkdClearTimer(pkd,1);
	pkdStartTimer(pkd,1);

	p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;++i,++p) {
		if (pkdIsGas(pkd,p)) {
#ifdef GLASSZ
		        p->a[0]=0; p->a[1]=0;
#endif
			for (j=0;j<3;++j) {
				p->vPred[j] = p->vPred[j]*dvFacOne + p->a[j]*dvFacTwo;
				}
#ifdef VARALPHA
			    {
			    double dalphadt;
			    dalphadt = - p->alphaPred*0.1*p->c/sqrt(0.25*p->fBall2);
#ifdef DODVDS
			    if (p->dvds < 0) dalphadt -= p->dvds*1.5;
#else
			    if (p->divv < 0) dalphadt -= p->divv;
#endif
			    p->alphaPred = p->alpha + dalphadt*duDelta;
			    if (p->alphaPred < ALPHAMIN) p->alphaPred = ALPHAMIN;
			    }
#endif
#ifdef SINKING
#if (0)
/* Vpred routines can be done in drift so I will do them there -- more efficient 
   Perhaps pkdkickvpred should be removed completely */
			if (TYPETest( p, TYPE_SINKING)) {
			    FLOAT r0 = p->rSinking0Mag;
			    /* For kick std dTime is midpoint of kick period, put dTime to end for this exact calc */
			    FLOAT r2 = r0 + p->vSinkingr0*(dTimeEnd-p->fSinkingTime);
			    FLOAT thfac, sqr02, th2, costh2, sinth2;
			    if (r2 < 0.1*r0) r2 = 0.1*r0; /* HACK */
			    thfac = p->vSinkingTang0Mag*2/(p->vSinkingr0);
			    sqr02 = sqrt(r0/r2);
			    th2 = thfac*(1-sqr02);
			    costh2 = cos(th2);
			    sinth2 = sin(th2);
			    /* v does not include motion around sink -- add it back */
			    for (j=0;j<3;j++) {
				p->vPred[j] += p->vSinkingTang0Mag*sqr02*(-sinth2*p->rSinking0Unit[j]+costh2*p->vSinkingTang0Unit[j])+p->vSinkingr0*(costh2*p->rSinking0Unit[j]+sinth2*p->vSinkingTang0Unit[j]);
				}
			    }
#endif /* (0) */
#ifdef SINKDBG
			if (p->iOrder == 55) printf("SINKINGKICKVPRED %d %g, %g %g  %g %g  %g %g\n",p->iOrder,dTimeEnd,p->vPred[0],p->v[0],p->vPred[1],p->v[1],p->vPred[2],p->v[2]);
#endif
#endif
			if (iGasModel != GASMODEL_ISOTHERMAL && iGasModel != GASMODEL_GLASS) {
#ifndef NOCOOLING
#ifdef COOLDEBUG
				if (p->uPred+p->uDot*duDelta < 0) 
					fprintf(stderr,"upred error %i: %g %g %g -> %g %i\n",p->iOrder,p->uPred,p->uDot,duDelta,p->uPred + p->uDot*duDelta,p->iRung);
#endif
			  p->uPred = p->uPred + p->uDot*duDelta;
			  if (p->uPred < 0) {
			      FLOAT uold = p->uPred - p->uDot*duDelta;
#ifdef FBPARTICLE
                  fprintf(stderr,"FBP Negative! %d: %g %g %g %g\n",p->iOrder,uold,p->u,p->uDot,duDelta);
#endif
			      p->uPred = uold*exp(p->uDot*duDelta/uold);
			      }
#ifdef UNONCOOL
			  p->uNoncoolPred = p->uNoncoolPred + p->uNoncoolDot*duDelta;
			  if (p->uNoncoolPred < 0) p->uNoncoolPred = 0;
#endif /* UNONCOOL */
#else /* NOCOOLING is defined: */
              p->uPred = p->uPred + UDOT_HYDRO(p)*duDelta;
#endif
#if defined(PRES_HK) || defined(PRES_MONAGHAN) || defined(SIMPLESF)
			  if (p->uPred < 0) p->uPred = 0;
#endif
			  mdlassert(pkd->mdl,p->uPred >= 0.0);
			  }
#ifdef DIFFUSION
			p->fMetalsPred = p->fMetalsPred + p->fMetalsDot*duDelta;
#ifdef MASSDIFF
			p->fMass = p->fMass + p->fMassDot*duDelta;
#endif
#ifdef STARFORM
			p->fMFracOxygenPred = p->fMFracOxygenPred + p->fMFracOxygenDot*duDelta;
			p->fMFracIronPred = p->fMFracIronPred + p->fMFracIronDot*duDelta;
#endif /* STARFORM */
#endif /* DIFFUSION */
		    }
		}

	pkdStopTimer(pkd,1);
	mdlDiag(pkd->mdl, "Done Vpred\n");
	}
#else /* NOT GASOLINE */
void
pkdKickVpred(PKD pkd,double dvFacOne,double dvFacTwo)
{
	PARTICLE *p;
	int i,j,n;

	p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;i++)
#ifdef GLASS
		    if (!pkdIsGas(pkd,&p[i])) {
			for (j=0;j<3;j++)
			    p[i].vPred[j] = p[i].vPred[j] + p[i].a[j]*dvFacTwo;
			}
		    else
#endif	
			{
			for (j=0;j<3;j++)
			    p[i].vPred[j] = p[i].vPred[j]*dvFacOne + p[i].a[j]*dvFacTwo;
			}
	}
#endif
#endif /* NEED_VPRED */

void
pkdCOM(PKD pkd, double *com)
{
    int i;
    int nLocal = pkdLocal(pkd);
	double m;
    
	com[0] = 0;
	com[1] = 0;
	com[2] = 0;
	com[3] = 0;
	com[4] = 0;
	com[5] = 0;
	com[6] = 0;
	com[7] = 0;
	com[8] = 0;
	com[9] = 0;
	com[10] = 0;
	com[11] = 0;

    for (i=0;i<nLocal;++i) {
	  m = pkd->pStore[i].fMass;
	  if ( TYPETest(&pkd->pStore[i], TYPE_GAS) ) {
		com[0] += m*pkd->pStore[i].r[0];
		com[1] += m*pkd->pStore[i].r[1];
		com[2] += m*pkd->pStore[i].r[2];
		com[3] += m;
		}
	  else if ( TYPETest(&pkd->pStore[i], TYPE_DARK) ) {
		com[4] += m*pkd->pStore[i].r[0];
		com[5] += m*pkd->pStore[i].r[1];
		com[6] += m*pkd->pStore[i].r[2];
		com[7] += m;
		}
	  else if ( TYPETest(&pkd->pStore[i], TYPE_STAR) ) {
		com[8] += m*pkd->pStore[i].r[0];
		com[9] += m*pkd->pStore[i].r[1];
		com[10] += m*pkd->pStore[i].r[2];
		com[11] += m;
		}
  	  }
    }

void
pkdCOMByType(PKD pkd, int type, double *com)
{
    int i;
    int nLocal = pkdLocal(pkd);
	double m;
    
	com[0] = 0;
	com[1] = 0;
	com[2] = 0;
	com[3] = 0;

    for (i=0;i<nLocal;++i) {
	  if ( TYPETest(&pkd->pStore[i], type) ) {
		m = pkd->pStore[i].fMass;
		com[0] += m*pkd->pStore[i].r[0];
		com[1] += m*pkd->pStore[i].r[1];
		com[2] += m*pkd->pStore[i].r[2];
		com[3] += m;
		}
  	  }
    }

void
pkdOldestStar(PKD pkd, double *com)
{
#if defined(STARFORM) || defined (SIMPLESF)
    int nLocal = pkdLocal(pkd);
    int i;
#endif

	com[0] = 0;
	com[1] = 0;
	com[2] = 0;
	com[3] = FLT_MAX;

#if defined(STARFORM) || defined (SIMPLESF)
    for (i=0;i<nLocal;++i) {
	  if ( TYPETest(&pkd->pStore[i], TYPE_STAR) && pkd->pStore[i].fTimeForm < com[3]) {
		com[0] = pkd->pStore[i].r[0];
		com[1] = pkd->pStore[i].r[1];
		com[2] = pkd->pStore[i].r[2];
		com[3] = pkd->pStore[i].fTimeForm;
		}
  	  }
#endif
   }

int pkdSetSink(PKD pkd, double dSinkMassMin)
{
#ifdef GASOLINE
    PARTICLE *p;
    int i,nSink = 0;
    int nLocal = pkdLocal(pkd);

    for(i=0;i<nLocal;++i) { 
		p = &pkd->pStore[i];
#ifdef STARSINK
                if ((TYPETest(p,TYPE_STAR))) {
#else
#if (0)
		if ((TYPETest(p,TYPE_STAR) && p->fMass >= dSinkMassMin) || p->fMetals < 0) {
#else 
		if ((TYPETest(p,TYPE_STAR) && p->fTimeForm < 0)) {
#endif
#endif
		    TYPESet(p,TYPE_SINK);
		    nSink++;
		    }
		}

    return nSink;
#else
    assert(0);	/* GASOLINE needed for sink particles to work */
    return -1; /* to keep compiler happy */
#endif
    }

void pkdFormSinks(PKD pkd, int bJeans, double dJConst2, int bDensity, double dDensityCut, double dTime, int iKickRung, int bSimple, int *nCandidates, double *pJvalmin)
{
#ifdef GASOLINE
    int i;
    PARTICLE *p;
    int n = pkdLocal(pkd);
    double Jval, Jvalmin = FLT_MAX;
#endif
    
    *nCandidates = 0;
    
#ifdef GASOLINE
    for(i = 0; i < n; ++i) {
        p = &pkd->pStore[i];
    
        if(TYPETest( p, TYPE_GAS ) && p->iRung >= iKickRung) {
#ifdef SINKING
	    if (TYPETest( p, TYPE_SINKING )) continue;
#endif
/* Jeans Mass compared to nJeans particle masses */
	    Jval =  dJConst2*(p->c*p->c*p->c*p->c*p->c*p->c)/(p->fMass*p->fMass*p->fDensity);
	    if (Jval < Jvalmin) Jvalmin = Jval;
	    if ((bJeans && Jval < 1) ||
		(bDensity && p->fDensity >= dDensityCut)) {
		(*nCandidates)++;
		if (bSimple) {
		    TYPESet(p, TYPE_SINK); /* Is now a sink! */
		    p->fMetals = -dTime;
		    }
		else {
		    TYPESet(p, TYPE_SMOOTHACTIVE); /* Is now a candidate */
		    }
		}
	    }
	}
    *pJvalmin = Jvalmin;
#endif
}

void pkdSinkLogInit(PKD pkd)
{
    SINKLOG *pSinkLog = &pkd->sinkLog;
    
    pSinkLog->nLog = 0;
    pSinkLog->nMaxLog = 1000;	/* inital size of buffer */
    pSinkLog->nLogOrdered = 0;
    pSinkLog->nFormOrdered = 0;
    pSinkLog->nForm = 0;
    pSinkLog->nAccrete = 0;
    pSinkLog->nMerge = 0;
    pSinkLog->SinkEventTab = malloc(pSinkLog->nMaxLog*sizeof(SINKEVENT));
}

void pkdSinkLogFlush(PKD pkd, char *pszFileName)
{
    FILE *fp;
    int iLog;
    XDR xdrs;
    
    if(pkd->sinkLog.nLog == 0)
	return;
    
    assert(pkd->sinkLog.nForm == pkd->sinkLog.nFormOrdered);
    
    fp = fopen(pszFileName, "a");
    assert(fp != NULL);
    xdrstdio_create(&xdrs,fp,XDR_ENCODE);
    /* Note: Could treat different events differently 
       -- Accretion events do not have useful info except iOrd of Victim 
          e.g. if iOrdSink == -1 (sink formation accrete), -2 (std accrete) */
    for(iLog = 0; iLog < pkd->sinkLog.nLog; iLog++){
	SINKEVENT *pSEv = &(pkd->sinkLog.SinkEventTab[iLog]); 
	xdr_int(&xdrs, &(pSEv->iOrdSink));
	xdr_int(&xdrs, &(pSEv->iOrdVictim));
	xdr_double(&xdrs, &(pSEv->time));
	xdr_double(&xdrs, &(pSEv->mass));
	xdr_double(&xdrs, &(pSEv->r[0]));
	xdr_double(&xdrs, &(pSEv->r[1]));
	xdr_double(&xdrs, &(pSEv->r[2]));
	xdr_double(&xdrs, &(pSEv->v[0]));
	xdr_double(&xdrs, &(pSEv->v[1]));
	xdr_double(&xdrs, &(pSEv->v[2]));
	xdr_double(&xdrs, &(pSEv->L[0]));
	xdr_double(&xdrs, &(pSEv->L[1]));
	xdr_double(&xdrs, &(pSEv->L[2]));
	}
    xdr_destroy(&xdrs);
    fclose(fp);
    pkd->sinkLog.nLog = 0;
    pkd->sinkLog.nForm = 0;
    pkd->sinkLog.nAccrete = 0;
    pkd->sinkLog.nMerge = 0;
    pkd->sinkLog.nLogOrdered = 0;
    pkd->sinkLog.nFormOrdered = 0;
    }
    


