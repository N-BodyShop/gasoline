#include <stdlib.h>
#include <stddef.h>
#include <malloc.h>
#include <math.h>
#include <assert.h>
#include <sys/time.h>

#include <rpc/types.h>
#include <rpc/xdr.h>

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
#include "aggs.h" /*DEBUG include this from collision.h?...*/
#endif

#ifdef _LARGE_FILES
#define fseek fseeko
#endif

/*
 BalsaraSwitch hackery
*/
/*
#define VISCSWITCHCUT 0.6
*/

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
				   FLOAT *fPeriod,int nDark,int nGas,int nStar)
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
	pkd->nMaxOrderGas = nGas - 1;
	pkd->nMaxOrderDark = nGas + nDark - 1;
	pkd->nRejects = 0;
	for (j=0;j<3;++j) {
		pkd->fPeriod[j] = fPeriod[j];
		}
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
	pkd->cl = clInit();
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
	clFinalize(pkd->cl);
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


void pkdGenericSeek(PKD pkd,FILE *fp,int nStart,int iHeader,int iElement)
{
	long lStart;

	lStart = iHeader + nStart*iElement;
	fseek(fp,lStart,0);
	}


void pkdReadTipsy(PKD pkd,char *pszFileName,int nStart,int nLocal,
				  int bStandard,int iReadIOrder,double dvFac,double dTuFac)
{
	FILE *fp,*fp2 = NULL;
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
		p->fBall2 = 0.0;
		p->fBallMax = 0.0;
#ifdef GASOLINE
		p->u = 0.0;
		p->uPred = 0.0;
#ifdef SUPERNOVA
                p->uSN = 0.0;
	        p->PdVSN = 0.0;
#endif
#ifdef STARFORM
		p->fESNrate = 0.0;
#endif
#ifdef SHOCKTRACK
		p->ShockTracker = 0.0;
#endif
#ifndef NOCOOLING		
		/* Place holders -- later fixed in pkdInitEnergy */
		p->Y_HI = 0.75;
		p->Y_HeI = 0.0625;
		p->Y_HeII = 0.0;
#endif
		p->c = 0.0;
		p->fMetals = 0.0;
		p->fTimeForm = 0.0;
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
		char atmp[160];
		sprintf(atmp,"%s.iord",pszFileName);
		fp2 = fopen(atmp,"r");
		mdlassert(pkd->mdl,fp2 != NULL);
		/*
		 ** Seek to right place in file, assumes 4 byte iOrders
		 */
		switch(iReadIOrder) {
		case 1:
			pkdGenericSeek(pkd,fp2,nStart,sizeof(int),sizeof(int));
			break;
		case 2:
			pkdGenericSeek(pkd,fp2,nStart,sizeof(int),sizeof(long long));
			break;
		case 3:
			pkdGenericSeek(pkd,fp2,nStart,sizeof(int),sizeof(pkd->pStore[0].iOrder));
			break;
		default:
			fprintf(stderr,"Don't understand iOrder format: %d\n",iReadIOrder);
			mdlassert(pkd->mdl,0);
			}
		}

	/*
	 ** Read Stuff!
	 */
	if (bStandard) {
		FLOAT vTemp;
		long long LongTmp;
		int IntTmp;
		XDR xdrs;
		xdrstdio_create(&xdrs,fp,XDR_DECODE);
		for (i=0;i<nLocal;++i) {
			p = &pkd->pStore[i];
			p->iOrder = nStart + i; /* temporary */
			if (pkdIsDarkByOrder(pkd,p)) {
				iSetMask = TYPE_DARK;
				xdr_float(&xdrs,&fTmp);
				p->fMass = fTmp;
				assert(p->fMass >= 0);
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
			else if (pkdIsGasByOrder(pkd,p)) {
				iSetMask = TYPE_GAS;
				xdr_float(&xdrs,&fTmp);
				p->fMass = fTmp;
				assert(p->fMass > 0);
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
				p->u = dTuFac*vTemp;
				p->uPred = dTuFac*vTemp;
#ifdef COOLDEBUG
				if (p->iOrder == 842079) fprintf(stderr,"Particle %i in pStore[%i]\n",p->iOrder,(int) (p-pkd->pStore));
				assert(p->u >= 0);
				assert(p->uPred >= 0);
#endif
				xdr_float(&xdrs,&fTmp);
				p->fSoft = fTmp;
#ifdef CHANGESOFT
				p->fSoft0 = fTmp;
#endif
				xdr_float(&xdrs,&fTmp);
				p->fMetals = fTmp;
#else
				xdr_float(&xdrs,&fTmp);
				xdr_float(&xdrs,&fTmp);
				xdr_float(&xdrs,&fTmp);
				p->fSoft = fTmp;
#ifdef CHANGESOFT
				p->fSoft0 = fTmp;
#endif
				xdr_float(&xdrs,&fTmp);
#endif
				xdr_float(&xdrs,&fTmp);
				p->fPot = fTmp;
				}
			else if (pkdIsStarByOrder(pkd,p)) {
				iSetMask = TYPE_STAR;
				xdr_float(&xdrs,&fTmp);
				p->fMass = fTmp;
				assert(p->fMass >= 0);
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
				xdr_float(&xdrs,&fTmp);
				p->fTimeForm = fTmp;
#else
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
			else mdlassert(pkd->mdl,0);

			TYPESet(p,iSetMask);
			switch (iReadIOrder) {
			case 0:
				break;
			case 1:
				fread(&IntTmp,sizeof(IntTmp),1,fp2);
				p->iOrder = IntTmp;
				break;
			case 2:
				fread(&LongTmp,sizeof(LongTmp),1,fp2);
				p->iOrder = LongTmp;
				break;
			case 3:
				fread(&p->iOrder,sizeof(p->iOrder),1,fp2);
				break;
				}
			}
		xdr_destroy(&xdrs);
		}
	else {
		long long LongTmp;
		int IntTmp;
		for (i=0;i<nLocal;++i) {
			p = &pkd->pStore[i];
			p->iOrder = nStart + i;

			if (pkdIsDarkByOrder(pkd,p)) {
				iSetMask = TYPE_DARK;
				fread(&dp,sizeof(struct dark_particle),1,fp);
				for (j=0;j<3;++j) {
					p->r[j] = dp.pos[j];
					p->v[j] = dvFac*dp.vel[j];
					}
				p->fMass = dp.mass;
				assert(p->fMass >= 0);
				p->fSoft = dp.eps;
#ifdef CHANGESOFT
				p->fSoft0 = dp.eps;
#endif
				p->fPot = dp.phi;
				}
			else if (pkdIsGasByOrder(pkd,p)) {
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
				assert(p->fMass >= 0);
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
				if (p->iOrder == 842079) fprintf(stderr,"Particle %i in pStore[%i]\n",p->iOrder,(int) (p-pkd->pStore));
				assert(p->u >= 0);
				assert(p->uPred >= 0);
#endif
				p->fMetals = gp.metals;
#endif
				}
			else if (pkdIsStarByOrder(pkd,p)) {
				iSetMask = TYPE_STAR;
				fread(&sp,sizeof(struct star_particle),1,fp);
				for (j=0;j<3;++j) {
					p->r[j] = sp.pos[j];
					p->v[j] = dvFac*sp.vel[j];
					}
				p->fMass = sp.mass;
				assert(p->fMass >= 0);
				p->fSoft = sp.eps;
#ifdef CHANGESOFT
				p->fSoft0 = sp.eps;
#endif
				p->fPot = sp.phi;
#ifdef GASOLINE
				p->fMetals = sp.metals;
				p->fTimeForm = sp.tform;		
#endif
				}
			else mdlassert(pkd->mdl,0);

			TYPESet(p,iSetMask); /* needed to get max order info */
			switch (iReadIOrder) {
			case 0:
				p->iOrder = nStart + i;
				break;
			case 1:
				fread(&IntTmp,sizeof(IntTmp),1,fp2);
				p->iOrder = IntTmp;
				break;
			case 2:
				fread(&LongTmp,sizeof(LongTmp),1,fp2);
				p->iOrder = LongTmp;
				break;
			case 3:
				fread(&p->iOrder,sizeof(p->iOrder),1,fp2);
				break;
				}
			}
		}
	if (fp2!=NULL) fclose(fp2);
	fclose(fp);
	}


void pkdCalcBound(PKD pkd,BND *pbnd,BND *pbndActive,BND *pbndTreeActive, BND *pbndBall)
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
	int nBuf;
	int nOutBytes,nSndBytes,nRcvBytes;

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
    int nBuf;
    int nOutBytes,nSndBytes,nRcvBytes;
    int i;
    int iBuf;
    
    /*
     ** Move particles to High memory.
     */
    iBuf = pkdSwapSpace(pkd);
    for (i=pkdLocal(pkd)-1;i>=0;--i)
	pkd->pStore[iBuf+i] = pkd->pStore[i];

    nBuf = pkdFreeStore(pkd)*sizeof(PARTICLE);
    nOutBytes = pkdLocal(pkd)*sizeof(PARTICLE);
    mdlSwap(pkd->mdl,idSwap,nBuf,&pkd->pStore[0], nOutBytes,
			&nSndBytes, &nRcvBytes);
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
	fp = fopen(pszFileName,"a");
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
				fTmp = p->fDensity;
				xdr_float(&xdrs,&fTmp);
#ifdef GASOLINE
				/*
				 ** Convert thermal energy to tempertature.
				 */
				switch (iGasModel) {
				case GASMODEL_COOLING:
				case GASMODEL_COOLING_NONEQM:
#ifndef NOCOOLING
					vTemp = clTemperature(2*(pkd->cl)->Y_H - p->Y_HI + 
										  3*(pkd->cl)->Y_He - 2*p->Y_HeI - p->Y_HeII,
										  p->u*(pkd->cl)->dErgPerGmUnit);
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
				xdr_float(&xdrs,&fTmp);
#ifdef DEBUG
				/* Store divv in metals for now */
				fTmp = p->divv;
#else
				fTmp = p->fMetals;
#endif
				xdr_float(&xdrs,&fTmp);
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
				xdr_float(&xdrs,&fTmp);
				}
			else if (pkdIsStar(pkd,p)) {
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
				gp.metals = p->fMetals;
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
				sp.mass = p->fMass;
#ifdef CHANGESOFT
				sp.eps = p->fSoft0;
#else
				sp.eps = p->fSoft;
#endif
				sp.phi = p->fPot;
#ifdef GASOLINE
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


void pkdSetSoft(PKD pkd,double dSoft)
{
	PARTICLE *p;
	int i,n;

	p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;++i) {
#ifdef CHANGESOFT
		p[i].fSoft0 = dSoft;
#else
		p[i].fSoft = dSoft;
#endif
		}
	}

#ifdef CHANGESOFT
void pkdPhysicalSoft(PKD pkd,double dSoftMax,double dFac,int bSoftMaxMul)
{
	PARTICLE *p;
	int i,n;

	p = pkd->pStore;
	n = pkdLocal(pkd);
	
	assert(dFac > 0);
	if (bSoftMaxMul) {
	        for (i=0;i<n;++i) {
		        assert(p[i].fSoft0 > 0);
	                p[i].fSoft = p[i].fSoft0*dFac;
			assert(p[i].fSoft > 0);
		        }
	        }
	else {
 	        assert(dSoftMax > 0);
	        for (i=0;i<n;++i) {
	                p[i].fSoft = p[i].fSoft0*dFac;
	                if (p[i].fSoft > dSoftMax) p[i].fSoft = dSoftMax;
		        }
	        }
	}

void pkdPreVariableSoft(PKD pkd)
{
	PARTICLE *p;
	int i,n;

	p = pkd->pStore;
	n = pkdLocal(pkd);
	
	for (i=0;i<n;++i) {
	        if (TYPEQueryACTIVE(&(p[i]))) p[i].fSoft = p[i].fBall2;
	        }
	}

void pkdPostVariableSoft(PKD pkd,double dSoftMax,int bSoftMaxMul)
{
	PARTICLE *p;
	int i,n;
	double dTmp;

	p = pkd->pStore;
	n = pkdLocal(pkd);
	
	if (bSoftMaxMul) {
	        for (i=0;i<n;++i) {
	                if (TYPEQueryACTIVE(&(p[i]))) {
		                  dTmp = sqrt(p[i].fBall2*.25);
	                          p[i].fBall2 = p[i].fSoft;
	                          p[i].fSoft = (dTmp <= p[i].fSoft0*dSoftMax ? dTmp : p[i].fSoft0*dSoftMax);
                                  }
		        }
	        }
	else {
	        for (i=0;i<n;++i) {
	                if (TYPEQueryACTIVE(&(p[i]))) {
		                  dTmp = sqrt(p[i].fBall2*.25);
	                          p[i].fBall2 = p[i].fSoft;
	                          p[i].fSoft = (dTmp <= dSoftMax ? dTmp : dSoftMax);
                                  }
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
		}
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
	double m,dx,dy,dz,d2,d1;

	/*
	 ** Initialize moments.
	 ** Initialize various B numbers.
	 */
	switch (iOrder) {	
	case 4:
		pcc->Hxxxx = 0.0;
		pcc->Hxyyy = 0.0;
		pcc->Hxxxy = 0.0;
		pcc->Hyyyy = 0.0;
		pcc->Hxxxz = 0.0;
		pcc->Hyyyz = 0.0;
		pcc->Hxxyy = 0.0;
		pcc->Hxxyz = 0.0;
		pcc->Hxyyz = 0.0;
		pcc->Hxxzz = 0.0;
		pcc->Hxyzz = 0.0;
		pcc->Hxzzz = 0.0;
		pcc->Hyyzz = 0.0;
		pcc->Hyzzz = 0.0;
		pcc->Hzzzz = 0.0;
		pcc->B6 = 0.0;
	case 3:
		pcc->Oxxx = 0.0;
		pcc->Oxyy = 0.0;
		pcc->Oxxy = 0.0;
		pcc->Oyyy = 0.0;
		pcc->Oxxz = 0.0;
		pcc->Oyyz = 0.0;
		pcc->Oxyz = 0.0;
		pcc->Oxzz = 0.0;
		pcc->Oyzz = 0.0;
		pcc->Ozzz = 0.0;
		pcc->B5 = 0.0;
	default:
		pcc->Qxx = 0.0;
		pcc->Qyy = 0.0;
		pcc->Qzz = 0.0;
		pcc->Qxy = 0.0;
		pcc->Qxz = 0.0;
		pcc->Qyz = 0.0;
		pcc->B2 = 0.0;
		pcc->B3 = 0.0;
		pcc->B4 = 0.0;
		}
	pcc->Bmax = 0.0;
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
		if (d1 > pcc->Bmax) pcc->Bmax = d1;
		pcc->B2 += m*d2;
		pcc->B3 += m*d2*d1;
		pcc->B4 += m*d2*d2;
		pcc->B5 += m*d2*d2*d1;
		pcc->B6 += m*d2*d2*d2;
#ifdef COMPLETE_LOCAL
		d2 = 0.0;
#endif
		switch (iOrder) {
		case 4:
			/*
			 ** Calculate reduced hexadecapole moment...
			 */
			pcc->Hxxxx += m*(dx*dx*dx*dx - 6.0/7.0*d2*(dx*dx - 0.1*d2));
			pcc->Hxyyy += m*(dx*dy*dy*dy - 3.0/7.0*d2*dx*dy);
			pcc->Hxxxy += m*(dx*dx*dx*dy - 3.0/7.0*d2*dx*dy);
			pcc->Hyyyy += m*(dy*dy*dy*dy - 6.0/7.0*d2*(dy*dy - 0.1*d2));
			pcc->Hxxxz += m*(dx*dx*dx*dz - 3.0/7.0*d2*dx*dz);
			pcc->Hyyyz += m*(dy*dy*dy*dz - 3.0/7.0*d2*dy*dz);
			pcc->Hxxyy += m*(dx*dx*dy*dy - 1.0/7.0*d2*(dx*dx + dy*dy - 0.2*d2));
			pcc->Hxxyz += m*(dx*dx*dy*dz - 1.0/7.0*d2*dy*dz);
			pcc->Hxyyz += m*(dx*dy*dy*dz - 1.0/7.0*d2*dx*dz);
			pcc->Hxxzz += m*(dx*dx*dz*dz - 1.0/7.0*d2*(dx*dx + dz*dz - 0.2*d2));
			pcc->Hxyzz += m*(dx*dy*dz*dz - 1.0/7.0*d2*dx*dy);
			pcc->Hxzzz += m*(dx*dz*dz*dz - 3.0/7.0*d2*dx*dz);
			pcc->Hyyzz += m*(dy*dy*dz*dz - 1.0/7.0*d2*(dy*dy + dz*dz - 0.2*d2));
			pcc->Hyzzz += m*(dy*dz*dz*dz - 3.0/7.0*d2*dy*dz);
			pcc->Hzzzz += m*(dz*dz*dz*dz - 6.0/7.0*d2*(dz*dz - 0.1*d2));
		case 3:
			/*
			 ** Calculate reduced octopole moment...
			 */
			pcc->Oxxx += m*(dx*dx*dx - 0.6*d2*dx);
			pcc->Oxyy += m*(dx*dy*dy - 0.2*d2*dx);
			pcc->Oxxy += m*(dx*dx*dy - 0.2*d2*dy);
			pcc->Oyyy += m*(dy*dy*dy - 0.6*d2*dy);
			pcc->Oxxz += m*(dx*dx*dz - 0.2*d2*dz);
			pcc->Oyyz += m*(dy*dy*dz - 0.2*d2*dz);
			pcc->Oxyz += m*dx*dy*dz;
			pcc->Oxzz += m*(dx*dz*dz - 0.2*d2*dx);
			pcc->Oyzz += m*(dy*dz*dz - 0.2*d2*dy);
			pcc->Ozzz += m*(dz*dz*dz - 0.6*d2*dz);;
		default:
			/*
			 ** Calculate quadrupole moment...
			 */
			pcc->Qxx += m*dx*dx;
			pcc->Qyy += m*dy*dy;
			pcc->Qzz += m*dz*dz;
			pcc->Qxy += m*dx*dy;
			pcc->Qxz += m*dx*dz;
			pcc->Qyz += m*dy*dz;
			}
		}
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
			c[iCell].r[j] = 0.0;
			}
		for (pj=l;pj<=u;++pj) {
			for (j=0;j<3;++j) {
				if (p[pj].r[j] < c[iCell].bnd.fMin[j])
					c[iCell].bnd.fMin[j] = p[pj].r[j];
				if (p[pj].r[j] > c[iCell].bnd.fMax[j])
					c[iCell].bnd.fMax[j] = p[pj].r[j];

				if (p[pj].r[j]-p[pj].fBallMax < c[iCell].bndBall.fMin[j])
					c[iCell].bndBall.fMin[j] = p[pj].r[j]-p[pj].fBallMax;
				if (p[pj].r[j]+p[pj].fBallMax > c[iCell].bndBall.fMax[j])
					c[iCell].bndBall.fMax[j] = p[pj].r[j]+p[pj].fBallMax;
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

	if (pLower > pUpper) return(0);
	else {
		/*
		 ** We need to find the bounding box.
		 */
		for (j=0;j<3;++j) {
			bnd.fMin[j] = pkd->pStore[pLower].r[j];
			bnd.fMax[j] = pkd->pStore[pLower].r[j];
			}
		for (i=pLower+1;i<=pUpper;++i) {
			for (j=0;j<3;++j) {
				if (pkd->pStore[i].r[j] < bnd.fMin[j]) 
					bnd.fMin[j] = pkd->pStore[i].r[j];
				else if (pkd->pStore[i].r[j] > bnd.fMax[j])
					bnd.fMax[j] = pkd->pStore[i].r[j];
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
		/*
		 ** We need to find the bounding box.
		 */
		for (j=0;j<3;++j) {
			pkdn->bnd.fMin[j] = pkd->pStore[pLower].r[j];
			pkdn->bnd.fMax[j] = pkd->pStore[pLower].r[j];
			pkdn->bndBall.fMin[j] = pkd->pStore[pLower].r[j]-pkd->pStore[pLower].fBallMax;
			pkdn->bndBall.fMax[j] = pkd->pStore[pLower].r[j]+pkd->pStore[pLower].fBallMax;
			}
		for (i=pLower+1;i<=pUpper;++i) {
			for (j=0;j<3;++j) {
				if (pkd->pStore[i].r[j] < pkdn->bnd.fMin[j]) 
					pkdn->bnd.fMin[j] = pkd->pStore[i].r[j];
				else if (pkd->pStore[i].r[j] > pkdn->bnd.fMax[j])
					pkdn->bnd.fMax[j] = pkd->pStore[i].r[j];
				if (pkd->pStore[i].r[j]-pkd->pStore[i].fBallMax < pkdn->bndBall.fMin[j]) 
					pkdn->bndBall.fMin[j] = pkd->pStore[i].r[j]-pkd->pStore[i].fBallMax;
				if (pkd->pStore[i].r[j]+pkd->pStore[i].fBallMax > pkdn->bndBall.fMax[j])
					pkdn->bndBall.fMax[j] = pkd->pStore[i].r[j]+pkd->pStore[i].fBallMax;
				}
			}	
		bGoodBounds = 0;
		for (j=0;j<3;++j) {
			if (pkdn->bnd.fMax[j] > pkdn->bnd.fMin[j])
			    bGoodBounds = 1;
			}
		if ((pUpper-pLower+1 > nBucket) && bGoodBounds) {
			/*
			 ** Now we need to determine the longest axis.
			 */
			d = 0;
			for (j=1;j<3;++j) {
				if (pkdn->bnd.fMax[j]-pkdn->bnd.fMin[j] > 
					pkdn->bnd.fMax[d]-pkdn->bnd.fMin[d]) {
					d = j;
					}
				}
			pkdn->iDim = d;
			/*
			 ** Now we do the split.
			 */
			pkdn->fSplit = 0.5*(pkdn->bnd.fMin[d]+pkdn->bnd.fMax[d]);
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
				fm = pkd->pStore[i].fMass;
				pkdn->fMass += fm;
				pkdn->fSoft += fm*pkd->pStore[i].fSoft;
				for (j=0;j<3;++j) {
					pkdn->r[j] += fm*pkd->pStore[i].r[j];
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
	pkdActiveTypeOrder(pkd, TYPE_ACTIVE|TYPE_TREEACTIVE );
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
	    printf("id:%d has an empty local tree\n",pkd->idSelf);
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
	pkdActiveTypeOrder(pkd, TYPE_ACTIVE|TYPE_TREEACTIVE );
	pkd->nBucket = nBucket;
	if (bTreeActiveOnly) n = pkd->nTreeActive;
	else n = pkd->nLocal;
	if (n == 0) {
		printf("id:%d has an empty local tree\n",pkd->idSelf);
		}
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
		pkdCalcBound(pkd,&bndDum,&bndDum,&c[pkd->iRoot].bnd,&bndDum);
		}
	else {
		pkdCalcBound(pkd,&c[pkd->iRoot].bnd,&bndDum,&bndDum,&bndDum);
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
			c[LOWER(i)].bndBall = c[i].bndBall; /* too big: fixed in pkdUpPass */
			c[LOWER(i)].pLower = c[i].pLower;
			c[LOWER(i)].pUpper = m;
			c[UPPER(i)].bnd = c[i].bnd;
			c[UPPER(i)].bnd.fMin[d] = c[i].fSplit;
			c[UPPER(i)].bndBall = c[i].bndBall; /* too big: fixed in pkdUpPass */
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
		   double fEwCut,double fEwhCut,int bDoSun,double *aSun,int *nActive,
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
			pkdStartTimer(pkd,1);
			pkdBucketWalk(pkd,iCell,nReps,iOrder); /* ignored in Flop count! */
			pkdStopTimer(pkd,1);
			c[iCell].bnd = bndTmp;
			*nActive += n;
			*pdPartSum += n*pkd->nPart + 
				n*(2*(c[iCell].pUpper-c[iCell].pLower) - n + 1)/2;
			*pdCellSum += n*pkd->nCellNewt;
			*pdSoftSum += n*pkd->nCellSoft;
			pkdStartTimer(pkd,2);
			dFlopI = pkdBucketInteract(pkd,iCell,iOrder);
			*pdFlop += dFlopI;
			pkdStopTimer(pkd,2);
			/*
			 * Now do Ewald part.
			 */
			if (bPeriodic && bEwald) {
				pkdStartTimer(pkd,3);
				dFlopE = pkdBucketEwald(pkd,iCell,nReps,fEwCut,iEwOrder);
				*pdFlop += dFlopE;
				pkdStopTimer(pkd,3);
				}
			else {
			    dFlopE = 0.0;
			    }
			fWeight = dFlopI + dFlopE;
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
			pkd->pStore[iDummy].fSoft = 0.0;
			c[pkd->iFreeCell].bnd.fMin[j] = -dTinyBox;
			c[pkd->iFreeCell].bnd.fMax[j] = dTinyBox;
			}
		TYPESet(&(pkd->pStore[iDummy]),TYPE_ACTIVE);
		pkd->pStore[iDummy].fPot = 0;
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


void pkdSunIndirect(PKD pkd,double *aSun,int bDoSun,double dSunMass)
{
	PARTICLE *p;
	double t,idt2;
	int i,j,n;

	p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;++i) {
		if (TYPEQueryACTIVE(&(p[i]))) {
			t = 0;
			for (j=0;j<3;++j) t += p[i].r[j]*p[i].r[j];
			t = (t == 0 ? 0 : 1/sqrt(t)); /* gravity at origin = zero */
			p[i].fPot -= dSunMass*t;
			t = t*t*t;
			idt2 = (p[i].fMass + dSunMass)*t;
			if (idt2 > p[i].dtGrav) p[i].dtGrav = idt2;
			/*
			 ** The way that aSun[j] gets added in here is confusing, but this
			 ** is the correct way! (aSun[] is the acceleration on the Sun).
			 */
			if (bDoSun) {
				t *= dSunMass;
				for (j=0;j<3;++j) {					
					p[i].a[j] -= (aSun[j] + p[i].r[j]*t);
					}				
				}
			else {
				t *= p[i].fMass;
				for (j=0;j<3;++j) {
					p[i].a[j] -= (aSun[j] - p[i].r[j]*t);
					}
				}
			}
		}
	}


void pkdLogHalo(PKD pkd)
{
	PARTICLE *p;
	int i,n;

	/* v in (kpc)/(4.691822922e16s) 128 km/s x 1.5233 */
	const double v = 194.9848486;
	const double d = 12.0;	/* in kpc */
	const double fl = 1.0;	/* flattening */
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
			double r = sqrt(x*x + y*y + z*z/(fl*fl));
			double C = 1.0/(r*r + d*d);	
			p[i].a[0] -= 2*v*v*x*C;
			p[i].a[1] -= 2*v*v*y*C;
			p[i].a[2] -= 2*v*v*z*C/(fl*fl);
			p[i].fPot += v*v*log(r*r + d*d);
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


void pkdNFWSpheroid(PKD pkd)
{
	PARTICLE *p;
	int i,n;

	const double M_200 = 2.;	/* Units 1e12 Solar masses */
	const double r_200 = 200;       /* Units kpc */
	const double G = 1;
	/* Assuming G = 1 (this sets a timescale) */
	/* TimeUnit = sqrt(kpc^3/(G*1e12 Msun)) = 1.1285945e+09 yr */
        /* Vunit = 2073.8081 km/s */

        const double c = 11; /* NFW concentration (cf. Lucio) */
	const double dSoft = 0.5; /* kpc */
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
			
			if (r < eps) {
			  fPot = G*M_const*c/r_200* ((1./6.)*( cx*cx - eps*eps )/(eps*(1+eps)*(1+eps)) 
			    - (1./3.)*eps/(1+eps)/(1+eps) - 1/(1+eps));

			  A = G*M_const* 		
			    (1./3.)/(eps*(1+eps)*(1+eps))
			    *r_200*r_200*r_200/(c*c*c);
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


void pkdHomogSpheroid(PKD pkd)
{
	PARTICLE *p;
	int i,n;

	const double M_s = 5;	
	const double r_s = 10;
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
			A = Mr/(r*r*r);
			p[i].a[0] -= A*x;
			p[i].a[1] -= A*y;
			p[i].a[2] -= A*z;
			p[i].fPot += 2*( r < r_s ? 0.5*(Mr/r-3*M_s/r_s) : -M_s/r );
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
			*Eth += p[i].fMass*p[i].u;
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
		exit(1);
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
pkdDrift(PKD pkd,double dDelta,FLOAT fCenter[3],int bPeriodic,int bFandG,
		 FLOAT fCentMass)
{
	PARTICLE *p;
	int i,j,n;
#ifdef SLIDING_PATCH
	FLOAT fShear;
#endif

	mdlDiag(pkd->mdl, "Into pkddrift\n");
	n = pkdLocal(pkd);
	for (i=0;i<n;++i) {
		p = &pkd->pStore[i];
#ifdef AGGS
		if (IS_AGG(p))
			continue;
#endif
#ifdef OLD_KEPLER
		if (p->iDriftType == KEPLER) /*DEBUG a bit ugly...*/
#else
		if (bFandG)
#endif
			fg(pkd->mdl,fCentMass + p->fMass,p->r,p->v,dDelta);
		else {
#ifdef SLIDING_PATCH
			fShear = 0;
#endif
			for (j=0;j<3;++j) {
				p->r[j] += dDelta*p->v[j];
				if (bPeriodic) {
					if (p->r[j] >= fCenter[j] + 0.5*pkd->fPeriod[j]) {
						p->r[j] -= pkd->fPeriod[j];
#ifdef SLIDING_PATCH
						if (j == 0) {
							fShear = 1.5*pkd->dOrbFreq*pkd->fPeriod[0];
							p->r[1] += SHEAR(-1,pkd->fPeriod[0],
											   pkd->fPeriod[1],pkd->dOrbFreq,
											   pkd->dTime + dDelta);
							}
#endif
						}
					if (p->r[j] < fCenter[j] - 0.5*pkd->fPeriod[j]) {
						p->r[j] += pkd->fPeriod[j];
#ifdef SLIDING_PATCH
						if (j == 0) {
							fShear = - 1.5*pkd->dOrbFreq*pkd->fPeriod[0];
							p->r[1] += SHEAR(1,pkd->fPeriod[0],
											   pkd->fPeriod[1],pkd->dOrbFreq,
											   pkd->dTime + dDelta);
							}
#endif
						}
					mdlassert(pkd->mdl,p->r[j] >= fCenter[j]-0.5*pkd->fPeriod[j]);
					mdlassert(pkd->mdl,p->r[j] <  fCenter[j]+0.5*pkd->fPeriod[j]);
					}
				}
#ifdef SLIDING_PATCH
			p->v[1] += fShear;
#endif
			}
		}
	mdlDiag(pkd->mdl, "Out of pkddrift\n");
	}


void pkdKick(PKD pkd, double dvFacOne, double dvFacTwo, double dvPredFacOne,
	     double dvPredFacTwo, double duDelta, double duPredDelta, int iGasModel,
	     double z, double duDotLimit )
{
	PARTICLE *p;
	int i,j,n;

	pkdClearTimer(pkd,1);
	pkdStartTimer(pkd,1);

	p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;++i,++p) {
		if (TYPEQueryACTIVE(p)) {
#ifdef NEED_VPRED
#ifdef GASOLINE
			if (pkdIsGas(pkd, p)) {
				for (j=0;j<3;++j) {
					p->vPred[j] = p->v[j]*dvPredFacOne + p->a[j]*dvPredFacTwo;
				}
				if (iGasModel != GASMODEL_ISOTHERMAL) {
#ifndef NOCOOLING				
				  p->uPred = p->u + p->uDot*duPredDelta;
				  p->u = p->u + p->uDot*duDelta;
#else
#ifdef COOLDEBUG
				if (p->u+p->uDot*duPredDelta < 0) {
					fprintf(stderr,"upred from u error %i: %g %g %g -> %g %i\n",p->iOrder,p->u,p->uDot,duPredDelta,p->uPred + p->uDot*duPredDelta,p->iRung);
					assert(0);
					}
				if (p->u+p->uDot*duDelta < 0) {
					fprintf(stderr,"u error %i: %g %g %g -> %g %i\n",p->iOrder,p->u,p->uDot,duDelta,p->uPred + p->uDot*duDelta,p->iRung);
					assert(0);
					}
#endif
				  p->uPred = p->u + p->PdV*duPredDelta;
				  p->u = p->u + p->PdV*duDelta;
#endif
#if defined(PRES_HK) || defined(PRES_MONAGHAN) 
				  if (p->uPred < 0) p->uPred = 0;
				  if (p->u < 0) p->u = 0;
#endif
#ifdef SUPERNOVA
				  p->uSN += p->PdVSN*duDelta;	  
#endif
				  }
				}
#else
			for (j=0;j<3;++j) {
				p->vPred[j] = p->v[j]*dvPredFacOne + p->a[j]*dvPredFacTwo;
				}
#endif
#endif /* NEED_VPRED */
			for (j=0;j<3;++j) {
				p->v[j] = p->v[j]*dvFacOne + p->a[j]*dvFacTwo;
				}
			}
	    }

	pkdStopTimer(pkd,1);
	mdlDiag(pkd->mdl, "Done pkdkick\n");
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
		p->fMass = cp.fMass;
		assert(p->fMass >= 0);
#ifdef CHANGESOFT
		p->fSoft0 = cp.fSoft;
#else
		p->fSoft = cp.fSoft;
#endif
		for (j=0;j<3;++j) {
			p->r[j] = cp.r[j];
			p->v[j] = cp.v[j];
#ifdef NEED_VPRED
			p->vPred[j] = cp.v[j];
#endif
			}
#ifdef GASOLINE
		p->u = cp.u;
		p->uPred = cp.u;
#ifdef COOLDEBUG
		if (p->iOrder == 842079) fprintf(stderr,"Particle %i in pStore[%i]\n",p->iOrder,(int) (p-pkd->pStore));
		assert(p->u >= 0);
#endif
#ifdef SUPERNOVA
                p->uSN = 0.0;
                p->PdVSN = 0.0;
#endif
#ifdef STARFORM
		p->fESNrate = 0.0;
		p->fTimeForm = cp.fTimeForm;
		for (j=0;j<3;++j) {
			p->rForm[j] = cp.rForm[j];
			p->vForm[j] = cp.vForm[j];
			}
		p->fDensity = cp.fDenForm;
		p->iGasOrder = cp.iGasOrder;
#endif
		p->fMetals = cp.fMetals;
#ifndef NOCOOLING		
		p->Y_HI = cp.Y_HI;
		p->Y_HeI = cp.Y_HeI;
		p->Y_HeII = cp.Y_HeII;
#endif
#endif
#ifdef COLLISIONS
		p->iOrgIdx = cp.iOrgIdx;
		for (j=0;j<3;++j)
			p->w[j] = cp.w[j];
		p->iColor = cp.iColor;
#endif
		TYPEClear(p);
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
	long lStart;
	int i,j,nLocal;
	int nout;

	/*
	 ** Seek past the header and up to nStart.
	 */
	fp = fopen(pszFileName,"r+");
	mdlassert(pkd->mdl,fp != NULL);
	lStart = iOffset+nStart*sizeof(CHKPART);
	fseek(fp,lStart,0);
	/* 
	 ** Write Stuff!
	 */
	nLocal = pkdLocal(pkd);
	for (i=0;i<nLocal;++i) {
		cp.iOrder = pkd->pStore[i].iOrder;
		cp.fMass = pkd->pStore[i].fMass;
#ifdef CHANGESOFT
		cp.fSoft = pkd->pStore[i].fSoft0;
#else
		cp.fSoft = pkd->pStore[i].fSoft;
#endif
		for (j=0;j<3;++j) {
			cp.r[j] = pkd->pStore[i].r[j];
			cp.v[j] = pkd->pStore[i].v[j];
			}
#ifdef GASOLINE
		cp.u = pkd->pStore[i].u;
		cp.fMetals = pkd->pStore[i].fMetals;
#ifndef NOCOOLING		
		cp.Y_HI = pkd->pStore[i].Y_HI;
		cp.Y_HeI = pkd->pStore[i].Y_HeI;
		cp.Y_HeII = pkd->pStore[i].Y_HeII;
#endif
#ifdef STARFORM
		cp.fTimeForm = pkd->pStore[i].fTimeForm;
		for (j=0;j<3;++j) {
			cp.rForm[j] = pkd->pStore[i].rForm[j];
			cp.vForm[j] = pkd->pStore[i].vForm[j];
			}
		cp.fDenForm = pkd->pStore[i].fDensity;
		cp.iGasOrder = pkd->pStore[i].iGasOrder;
#endif
#endif
#ifdef COLLISIONS
		cp.iOrgIdx = pkd->pStore[i].iOrgIdx;
		for (j=0;j<3;++j)
			cp.w[j] = pkd->pStore[i].w[j];
		cp.iColor = pkd->pStore[i].iColor;
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

	for (i=0;i<pkdLocal(pkd);++i) {
		dMass += pkd->pStore[i].fMass;
		}
	iRej = pkdFreeStore(pkd) - pkd->nRejects;
	for (i=0;i<pkd->nRejects;++i) {
		dMass += pkd->pStore[iRej+i].fMass;
		}
	return(dMass);
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
    
    nActive = 0;
    for (i=0;i<pkdLocal(pkd);++i) {
#ifdef COOLDEBUG
		if (pkd->pStore[i].iOrder == 842079) fprintf(stderr,"Particle %i in pStore[%i] (Rung)\n",pkd->pStore[i].iOrder, i);
		assert(pkd->pStore[i].u >= 0);
		assert(pkd->pStore[i].uPred >= 0);
#endif
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
			mdlassert(pkd->mdl,pkd->pStore[i].dtGrav);
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
			mdlassert(pkd->mdl,vel >= 0);
			vel = sqrt(vel)*dVelFac;
			acc = sqrt(acc)*dAccFac;
			dT = FLOAT_MAXVAL;
			if (bEpsAcc && acc>0) {
#ifdef GASOLINE
			     if (pkdIsGas(pkd, &(pkd->pStore[i])) && dhMinOverSoft < 1 && pkd->pStore[i].fBall2<4.0*pkd->pStore[i].fSoft*pkd->pStore[i].fSoft) {
			        if (pkd->pStore[i].fBall2 > 4.0*dhMinOverSoft
				    *pkd->pStore[i].fSoft*pkd->pStore[i].fSoft) 
				   dT = dEta*sqrt(sqrt(0.25*pkd->pStore[i].fBall2)/acc);
			        else 
				   dT = dEta*sqrt((dhMinOverSoft*pkd->pStore[i].fSoft)/acc);
			        }
		             else 
#endif
			        dT = dEta*sqrt(pkd->pStore[i].fSoft/acc);
			        }
			if (bSqrtPhi) {
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
pkdDtToRung(PKD pkd,int iRung,double dDelta,int iMaxRung,int bAll, int *pnMaxRung)
{
    int i;
    int iMaxRungOut;
    int iTempRung;
    int iSteps;
    int nMaxRung;
    
    iMaxRungOut = 0;
    nMaxRung = 0;
    for(i=0;i<pkdLocal(pkd);++i) {
		if(pkd->pStore[i].iRung >= iRung) {
			mdlassert(pkd->mdl,TYPEQueryACTIVE(&(pkd->pStore[i])));
			if(bAll) {          /* Assign all rungs at iRung and above */
			        assert(pkd->pStore[i].fSoft > 0);
			        assert(pkd->pStore[i].dt > 0);
				iSteps = floor(dDelta/pkd->pStore[i].dt);
				/* insure that integer boundary goes
				   to the lower rung. */
				if(fmod(dDelta,pkd->pStore[i].dt) == 0.0)
				    iSteps--;
				iTempRung = iRung;
				if(iSteps < 0)
				    iSteps = 0;
				while(iSteps) {
					++iTempRung;
					iSteps >>= 1;
					}
				if(iTempRung >= iMaxRung)
					iTempRung = iMaxRung-1;
				pkd->pStore[i].iRung = iTempRung;
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
pkdColNParts(PKD pkd, int *pnNew, int *nDeltaGas, int *nDeltaDark,
	     int *nDeltaStar)
{
    int pi, pj;
    int nNew;
    int ndGas;
    int ndDark;
    int ndStar;
    int newnLocal;
    PARTICLE *p;
    
    nNew = 0;
    ndGas = 0;
    ndDark = 0;
    ndStar = 0;
    newnLocal = pkdLocal(pkd);
    for(pi = 0, pj = 0; pi < pkdLocal(pkd); pi++) {
	if(pj < pi)
	    pkd->pStore[pj] = pkd->pStore[pi];
	p = &pkd->pStore[pi];
	if(p->iOrder == -1) {
	    ++pj;
	    ++nNew;
#ifdef GASOLINE
	    ++ndStar;
#else
	    ++ndDark;
#endif
	    if(TYPEQueryACTIVE(p))
		++pkd->nActive;
	    continue;
	    }
	else if(p->iOrder < -1){
	    --newnLocal;
	    p->iOrder = -2 - p->iOrder;
	    if(pkdIsGas(pkd, p))
		--ndGas;
	    else if(pkdIsDark(pkd, p))
		--ndDark;
	    else if(pkdIsStar(pkd, p))
		--ndStar;
	    else
		mdlassert(pkd->mdl,0);
	    if(TYPEQueryACTIVE(p))
		--pkd->nActive;
	    }
	else {
	    ++pj;
	    }
	}

    *pnNew = nNew;
    *nDeltaGas = ndGas;
    *nDeltaDark = ndDark;
    *nDeltaStar = ndStar;
    pkd->nLocal = newnLocal;
    }

void
pkdNewOrder(PKD pkd,int nStart)
{
    int pi;
    
    for(pi=0;pi<pkdLocal(pkd);pi++) {
		if(pkd->pStore[pi].iOrder == -1) {
			pkd->pStore[pi].iOrder = nStart++;
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
/*  pkd->nMaxOrder = nMaxOrder; */
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

    for(i=0;i<pkdLocal(pkd);++i) { 
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
		assert(p->u >= 0);
		assert(p->uPred >= 0);
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

    for(i=0;i<pkdLocal(pkd);++i) {
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

    for(i=0;i<pkdLocal(pkd);++i) {
		p = &pkd->pStore[i];
#ifdef COOLDEBUG
		if (p->iOrder == 842079) fprintf(stderr,"Particle %i in pStore[%i] Active Type\n",p->iOrder,(int) (p-pkd->pStore));
		assert(p->u >= 0);
		assert(p->uPred >= 0);
#endif
		/* DEBUG: Paranoia check */
		mdlassert(pkd->mdl,TYPETest(p,TYPE_ALL));
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
    
    nActive = 0;
    for(i=0;i<pkdLocal(pkd);++i) {
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

		TYPESet(p,iSetMask);
		}
    }

void pkdGrowMass(PKD pkd,int nGrowMass, double dDeltaM)
{
#ifdef GROWMASS
    int i;

    for(i=0;i<pkdLocal(pkd);++i) {
		if (pkd->pStore[i].iOrder < nGrowMass) {
			pkd->pStore[i].fMass += dDeltaM;
			}
		}
#endif
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
#ifdef SUPERNOVA

struct outCountSupernova pkdCountSupernova(PKD pkd, double dMetal, double dRhoCut, double dTMin, double dTMax,
					   double duTFac,int iGasModel)
{
    PARTICLE *p;
    int i;
    struct outCountSupernova ret;
    double Temp;

    ret.dMassMetalRhoCut=0;
    ret.dMassMetalTotal=0;
    ret.dMassNonMetalRhoCut=0;
    ret.dMassNonMetalTotal=0;
    ret.dMassTotal=0;

    p = pkd->pStore;

    if (dTMin>0 || dTMax < 1e20) {
      for(i=0;i<pkdLocal(pkd);++i,++p) {
	ret.dMassTotal += p->fMass;
	if (pkdIsGas(pkd,p)) {
	  switch (iGasModel) {
	  case GASMODEL_COOLING:
	  case GASMODEL_COOLING_NONEQM:
#ifndef NOCOOLING
	    Temp = clTemperature(2*(pkd->cl)->Y_H - p->Y_HI + 
			    3*(pkd->cl)->Y_He - 2*p->Y_HeI - p->Y_HeII,
			    p->u*(pkd->cl)->dErgPerGmUnit);
#endif
	    break;
	  default:
	    Temp = duTFac*p->u;
	  }
	  if (p->fMetals > dMetal) {
	    ret.dMassMetalTotal += p->fMass;
	    if (p->fDensity > dRhoCut && Temp > dTMin && Temp <= dTMax)
	      ret.dMassMetalRhoCut += p->fMass;
	  }
	  else {
	    ret.dMassNonMetalTotal += p->fMass;
	    if (p->fDensity > dRhoCut && Temp > dTMin && Temp <= dTMax)
	      ret.dMassNonMetalRhoCut += p->fMass;
	  }
	}
      }
    }
    else {
      for(i=0;i<pkdLocal(pkd);++i,++p) {
	ret.dMassTotal += p->fMass;
	if (pkdIsGas(pkd,p)) {
	  if (p->fMetals > dMetal) {
	    ret.dMassMetalTotal += p->fMass;
	    if (p->fDensity > dRhoCut)
	      ret.dMassMetalRhoCut += p->fMass;
	  }
	  else {
	    ret.dMassNonMetalTotal += p->fMass;
	    if (p->fDensity > dRhoCut)
	      ret.dMassNonMetalRhoCut += p->fMass;
	  }
	}
      }
    }
    return ret;
}

void pkdAddSupernova(PKD pkd, double dMetal, double dRhoCut, double dTMin, double dTMax, 
		     double duTFac,int iGasModel, double dPdVMetal, double dPdVNonMetal)
{
    PARTICLE *p;
    int i;
    double Temp;

    p = pkd->pStore;

    if (dTMin>0 || dTMax < 1e20) {
      for(i=0;i<pkdLocal(pkd);++i,++p) {
	if (TYPEQueryACTIVE(p) && pkdIsGas(pkd,p)) {
	  switch (iGasModel) {
	  case GASMODEL_COOLING:
	  case GASMODEL_COOLING_NONEQM:
#ifndef NOCOOLING
	    Temp = clTemperature(2*(pkd->cl)->Y_H - p->Y_HI + 
			    3*(pkd->cl)->Y_He - 2*p->Y_HeI - p->Y_HeII,
			    p->u*(pkd->cl)->dErgPerGmUnit);
#endif
	    break;
	  default:
	    Temp = duTFac*p->u;
	  }
	  if (p->fMetals > dMetal) {
	    if (p->fDensity > dRhoCut && Temp > dTMin && Temp <= dTMax) {
	      if ((p->iOrder%1000)==0) printf("Cluster: %i fDensity %g u %g PdV %g SN_PdV %g\n",p->iOrder,p->fDensity,p->u,p->PdV, dPdVMetal);
	      p->PdV += dPdVMetal;
	      p->PdVSN = dPdVMetal;
	    }
	  }
	  else {
	    if (p->fDensity > dRhoCut && Temp > dTMin && Temp <= dTMax) {
	      if ((p->iOrder%1000)==0) printf("Non-Cluster: %i fDensity %g u %g PdV %g SN_PdV %g\n",p->iOrder,p->fDensity,p->u,p->PdV, dPdVNonMetal);
	      p->PdV += dPdVNonMetal;
	      p->PdVSN = dPdVNonMetal;
	    }
	  }
	}
      }
    }
    else {
      for(i=0;i<pkdLocal(pkd);++i,++p) {
	if (TYPEQueryACTIVE(p) && pkdIsGas(pkd,p)) {
	  if (p->fMetals > dMetal) {
	    if (p->fDensity > dRhoCut) {
	      if ((p->iOrder%1000)==0) printf("Cluster: %i fDensity %g u %g PdV %g SN_PdV %g\n",p->iOrder,p->fDensity,p->u,p->PdV, dPdVMetal);
	      p->PdV += dPdVMetal;
	      p->PdVSN = dPdVMetal;
	    }
	  }
	  else {
	    if (p->fDensity > dRhoCut) {
	      if ((p->iOrder%1000)==0) printf("Non-Cluster: %i fDensity %g u %g PdV %g SN_PdV %g\n",p->iOrder,p->fDensity,p->u,p->PdV, dPdVNonMetal);
	      p->PdV += dPdVNonMetal;
	      p->PdVSN = dPdVNonMetal;
	    }
	  }
	}
      }
    }
}


#endif

void pkdUpdateuDot(PKD pkd, double duDelta, double dTime, double z, int iGasModel, int bUpdateY )
{
	PARTICLE *p;
	int i,n;
	int bCool = 0;
	CL *cl = NULL;
	PERBARYON Y;
	double E,dt = 0;

	pkdClearTimer(pkd,1);
	pkdStartTimer(pkd,1);

#ifndef NOCOOLING
	switch (iGasModel) {
	case GASMODEL_COOLING:
	case GASMODEL_COOLING_NONEQM:
		bCool = 1;
		cl = pkd->cl;
		clRatesRedshift( cl, z, dTime );
		dt = duDelta * cl->dSecUnit;
		break;
		}

	p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;++i,++p) {
		if(TYPEFilter(p,TYPE_GAS|TYPE_ACTIVE,TYPE_GAS|TYPE_ACTIVE)) {
			if (bCool) {
				pkdPARTICLE2PERBARYON(&Y,p,cl->Y_H,cl->Y_He);
				E = p->u*cl->dErgPerGmUnit;
#ifdef STARFORM
				cl->p = p;
				clIntegrateEnergy(cl, &Y, &E, (p->fESNrate + p->PdV)
								  *cl->dErgPerGmPerSecUnit, 
								  p->fDensity*cl->dComovingGmPerCcUnit, dt);
#else
				clIntegrateEnergy(cl, &Y, &E, p->PdV*cl->dErgPerGmPerSecUnit, 
								  p->fDensity*cl->dComovingGmPerCcUnit, dt);
#endif

				p->uDot = (E*cl->diErgPerGmUnit - p->u)/duDelta;
#ifdef COOLDEBUG
				if (E*cl->diErgPerGmUnit<1e-6*p->u || p->iOrder == 784461 || p->iOrder == 602270 || p->iOrder == 96299 || p->iOrder == 722701) 
					fprintf(stderr,"udot error? %i: %g %g %g -> %g %i\n",p->iOrder,p->uPred,p->uDot,duDelta,p->u + p->uDot*duDelta,p->iRung);
#endif
				if (bUpdateY) pkdPERBARYON2PARTICLE(&Y, p);
				}
			else { 
#ifdef STARFORM
				p->uDot = p->PdV + p->fESNrate;
#else
				p->uDot = p->PdV;
#endif
				}
			}
		}
#endif
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
       p->PdV,((p->a[0]*p->aPres[0])+(p->a[1]*p->aPres[1])+(p->a[2]*p->aPres[2])),sqrt(a2),sqrt(ap2),
       sqrt(a2*0.25*p->fBall2)/(p->c*p->c),
       (p->a[0]*p->gradrho[0]+p->a[1]*p->gradrho[1]+p->a[2]*p->gradrho[2])/
       (p->gradrho[0]*p->gradrho[0]+p->gradrho[1]*p->gradrho[1]+p->gradrho[2]*p->gradrho[2]+
	0.4/p->fBall2)/(p->c*p->c),p->BalsaraSwitch );
			}
#endif

			/* Rarefaction or Gravity dominated compression */
			if ( p->PdV < 0 || adotap < 0.5*a2) p->ShockTracker = 0;
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

/* Note: Uses uPred */
void pkdAdiabaticGasPressure(PKD pkd, double gammam1, double gamma)
{
    PARTICLE *p;
    double PoverRho;
    int i;

    p = pkd->pStore;
    for(i=0;i<pkdLocal(pkd);++i,++p) {
		if (pkdIsGas(pkd,p)) {
			PoverRho = gammam1*p->uPred;
			p->PoverRho2 = PoverRho/p->fDensity;
   			p->c = sqrt(gamma*PoverRho);
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


void pkdLowerSoundSpeed(PKD pkd, double dhMinOverSoft)
{
    PARTICLE *p;
    double PoverRho;
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
    CL *cl;
    PERBARYON Y;
    RATE r;
    double T;

#ifndef NOCOOLING
    cl = pkd->cl;
    clRatesRedshift( cl, z, dTime );
#endif

    p = pkd->pStore;
    for(i=0;i<pkdLocal(pkd);++i,++p) {
		if (TYPEQueryTREEACTIVE(p) && pkdIsGas(pkd,p)) {
#ifndef NOCOOLING
			T = p->u / dTuFac;
			p->Y_HI = cl->Y_H;
			p->Y_HeI = cl->Y_He;
			p->Y_HeII = 0.0;
			pkdPARTICLE2PERBARYON(&Y,p,cl->Y_H,cl->Y_He);
			clRates(cl,&r,T);
			clAbunds(cl,&Y,&r,p->fDensity*cl->dComovingGmPerCcUnit);
			pkdPERBARYON2PARTICLE(&Y,p);
			p->u = clThermalEnergy(Y.Total,T)*cl->diErgPerGmUnit;
#endif
			p->uPred = p->u;
#ifdef DEBUG
			if ((p->iOrder % 1000)==0) {
				printf("InitEnergy %i: %f %g   %f %f %f %g\n",
					   p->iOrder,T,p->u * cl->dErgPerGmUnit,
					   Y.HI,Y.HeI,Y.HeII,p->fDensity*cl->dComovingGmPerCcUnit);
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
    struct inGetGasPressure *in = vin;

    for(i=0;i<pkdLocal(pkd);++i) {
                p = &pkd->pStore[i];
                if (TYPEQueryTREEACTIVE(p)) {
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

			p->u = PoverRho;
			p->uPred = PoverRho;
			p->PoverRho2 = PoverRho/p->fDensity;
   			p->c = sqrt(in->gamma*PoverRho);
		        }
#ifdef DEBUG
		if (pkdIsGas(pkd,p) && (p->iOrder % 1000)==0) 
		        printf("Glass P %i: %i %i %f %f %f  %f %f %f %f %f\n",
			       p->iOrder,TYPEQueryACTIVE(p),TYPEQueryTREEACTIVE(p),
			       p->r[0],p->r[1],p->r[2],sqrt(0.25*p->fBall2),p->fDensity,p->uPred,p->PoverRho2*p->fDensity*p->fDensity,p->c);
#endif
                }
    }

#endif

/*
void pkdKickVpred(PKD pkd, double dvFacOne, double dvFacTwo, double duDelta,int iGasModel,
		  double z, double duDotLimit)
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
			for (j=0;j<3;++j) {
				p->vPred[j] = p->vPred[j]*dvFacOne + p->a[j]*dvFacTwo;
				}
			if (iGasModel != GASMODEL_ISOTHERMAL) {
#ifndef NOCOOLING
			  p->uPred = p->uPred + p->uDot*duDelta;
#else
			  p->uPred = p->uPred + p->PdV*duDelta;
#endif
#if defined(PRES_HK) || defined(PRES_MONAGHAN) 
                          if (p->uPred < 0) p->uPred = 0;
#endif
			  mdlassert(pkd->mdl,p->uPred > 0);
                          }
			}
		}

	pkdStopTimer(pkd,1);
	mdlDiag(pkd->mdl, "Done Vpred\n");
	}
*/

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
			       (p+i)->a[0],(p+i)->a[1],(p+i)->a[2],(p+i)->PdV);
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

void
pkdSphStep(PKD pkd, double dCosmoFac, double dEtaCourant, double dEtauDot, int bViscosityLimitdt)
{
    int i;
    PARTICLE *p;    
    double dT,dTu;

    for(i=0;i<pkdLocal(pkd);++i) {
        p = &pkd->pStore[i];
        if(pkdIsGas(pkd, p) && TYPEQueryACTIVE(p)) {
			/*
			 * Courant condition goes here.
			 */
	       if (p->mumax>0.0) {
				if (bViscosityLimitdt) 
				  dT = dEtaCourant*dCosmoFac*(sqrt(0.25*p->fBall2)/(p->c + 0.6*(p->c + 2*p->BalsaraSwitch*p->mumax)));
				else
				  dT = dEtaCourant*dCosmoFac*(sqrt(0.25*p->fBall2)/(p->c + 0.6*(p->c + 2*p->mumax)));
	                   }
	       else
			   dT = dEtaCourant*dCosmoFac*(sqrt(0.25*p->fBall2)/(1.6*p->c));

	       if (dEtauDot > 0.0 && p->PdV < 0.0) { /* Prevent rapid adiabatic cooling */
			    assert(p->u > 0);
				dTu = dEtauDot*p->u/fabs(p->PdV);
				if (dTu < dT) 
					dT = dTu;
				}
		if(dT < p->dt)
				p->dt = dT;
			}
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

			if(pkdIsGas(pkd, p)) {
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
/* Currently wired to have no more than two regions with different
   Pressures (densities) split by x=0 with a linear connection */
void pkdGlassGasPressure(PKD pkd, void *vin)
{
    PARTICLE *p;
    double PoverRho,xx,nsp=2.5;
    int i;
    struct inGetGasPressure *in = vin;

    for(i=0;i<pkdLocal(pkd);++i) {
		p = &pkd->pStore[i];
		if (TYPEQueryTREEACTIVE(p)) {
			if (p->r[0] < -nsp*in->dGlassxL) {
				if (p->r[0] > in->dxBoundL + nsp*in->dGlassxL)
					PoverRho=in->dGlassPoverRhoL;  
				else {
					xx =	 (p->r[0] - in->dxBoundL + nsp*in->dGlassxR)/
						(nsp*in->dGlassxL+nsp*in->dGlassxR);
			       xx = xx*xx*(-2*xx + 3);
					PoverRho = in->dGlassPoverRhoR +
						xx*(in->dGlassPoverRhoL - in->dGlassPoverRhoR);
					}
				}
			else if (p->r[0] > nsp*in->dGlassxR) {
				if (p->r[0] < in->dxBoundR - nsp*in->dGlassxR)
					PoverRho=in->dGlassPoverRhoR;  
				else {
					xx = (p->r[0] - in->dxBoundR + nsp*in->dGlassxR)/
						(nsp*in->dGlassxL+nsp*in->dGlassxR);
					xx = xx*xx*(-2*xx + 3);
					PoverRho = in->dGlassPoverRhoR +
						xx*(in->dGlassPoverRhoL - in->dGlassPoverRhoR);
					}
				}
			else {
				xx = (p->r[0] + nsp*in->dGlassxL)/
					(nsp*in->dGlassxL+nsp*in->dGlassxR);
				xx = xx*xx* ( -2*xx + 3 );
				PoverRho = in->dGlassPoverRhoL +
					xx*(in->dGlassPoverRhoR - in->dGlassPoverRhoL);
				}

			p->u = PoverRho;
			p->uPred = PoverRho;
			p->PoverRho2 = PoverRho/p->fDensity;
   			p->c = sqrt(in->gamma*PoverRho);
			}
#ifdef DEBUG
		if (pkdIsGas(pkd,p) && (p->iOrder % 1000)==0) 
			printf("Glass P %i: %i %i %f %f %f  %f %f %f %f %f\n",
			       p->iOrder,TYPEQueryACTIVE(p),TYPEQueryTREEACTIVE(p),
			       p->r[0],p->r[1],p->r[2],sqrt(0.25*p->fBall2),p->fDensity,p->uPred,p->PoverRho2*p->fDensity*p->fDensity,p->c);
#endif
		}
    }

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
		p->iOrgIdx = data.org_idx;
		p->fMass = data.mass;
		p->fSoft = 0.5*data.radius;
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
#ifdef CHANGESOFT
		data.radius = 2*p->fSoft0;
#else
		data.radius = 2*p->fSoft;
#endif
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
/*DEBUG			p->idNbr.iPid = p->idNbr.iIndex =
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

void
pkdPatch(PKD pkd,double dOrbFreqZ2)
{
	PARTICLE *p;
	int i;

	for (i=0;i<pkd->nLocal;i++) {
		p = &pkd->pStore[i];
		if (!TYPEQueryACTIVE(p)) continue;
		/*
		 ** Apply Hill's Equations, using *predicted* velocities...
		 */
		p->a[0] += pkd->dOrbFreq*(2*p->vPred[1] + 3*pkd->dOrbFreq*p->r[0]);
		p->a[1] -= 2*pkd->dOrbFreq*p->vPred[0];
		p->a[2] -= dOrbFreqZ2*p->r[2];
		}
	}

#endif /* SLIDING_PATCH */

#ifdef SIMPLE_GAS_DRAG
void
_pkdSimpleGasVel(double r[],int iFlowOpt,double u[])
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

void
pkdSimpleGasDrag(PKD pkd,int iFlowOpt,int bEpstein,double dGamma,
				 double dOmegaZ)
{
	PARTICLE *p = pkd->pStore;
	double u[3],dCoef;
	int i;

	for (i=0;i<pkd->nLocal;i++) {
		if (!TYPEQueryACTIVE(&p[i])) continue;
		_pkdSimpleGasVel(p[i].r,iFlowOpt,u); /* get gas velocity */
		if (bEpstein)
			dCoef = -dGamma/(2*p[i].fSoft); /* i.e. -gamma/R */
		else
			dCoef = -dGamma/(4*p[i].fSoft*p[i].fSoft); /* i.e. -gamma/R^2 */
		/* note: following assumes dust density << gas density */
		p[i].a[0] += dCoef*(p[i].vPred[0] - u[0]);
		p[i].a[1] += dCoef*(p[i].vPred[1] - u[1]);
		p[i].a[2] += dCoef*(p[i].vPred[2] - u[2]);
		}
	}
#endif /* SIMPLE_GAS_DRAG */

#ifdef NEED_VPRED
#ifdef GASOLINE
void
pkdKickVpred(PKD pkd,double dvFacOne,double dvFacTwo,double duDelta,
			 int iGasModel,double z,double duDotLimit)
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
			for (j=0;j<3;++j) {
				p->vPred[j] = p->vPred[j]*dvFacOne + p->a[j]*dvFacTwo;
				}
			if (iGasModel != GASMODEL_ISOTHERMAL) {
#ifndef NOCOOLING
#ifdef COOLDEBUG
				if (p->uPred+p->uDot*duDelta < 0) 
					fprintf(stderr,"upred error %i: %g %g %g -> %g %i\n",p->iOrder,p->uPred,p->uDot,duDelta,p->uPred + p->uDot*duDelta,p->iRung);
#endif
			  p->uPred = p->uPred + p->uDot*duDelta;
#else
			  p->uPred = p->uPred + p->PdV*duDelta;
#endif
#if defined(PRES_HK) || defined(PRES_MONAGHAN) 
			  if (p->uPred < 0) p->uPred = 0;
#endif
			  mdlassert(pkd->mdl,p->uPred >= 0);
			  }
			}
		}

	pkdStopTimer(pkd,1);
	mdlDiag(pkd->mdl, "Done Vpred\n");
	}
#else
void
pkdKickVpred(PKD pkd,double dvFacOne,double dvFacTwo)
{
	PARTICLE *p;
	int i,j,n;

	p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;i++)
		for (j=0;j<3;j++)
			p[i].vPred[j] = p[i].vPred[j]*dvFacOne + p[i].a[j]*dvFacTwo;
	}
#endif
#endif /* NEED_VPRED */
