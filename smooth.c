
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <limits.h>
#include <assert.h>
#include "smooth.h"
#include "pkd.h"
#include "smoothfcn.h"

/* Mark is for indicating which local particles are already in the PQ 
   This is needed because nbr search uses multiple approaches that re-find same nbr 
   INQUEUE does same thing for non-local nbrs using a hash approach */
#ifdef MARK
#define smTestMARK(pi__) smx->piMark[pi__]
#define smSetMARK(pi__) smx->piMark[pi__]=1
#define smResetMARK(pi__) smx->piMark[pi__]=0
#else
#define smTestMARK(pi__) TYPETest( (&p[pi__]), TYPE_MARK ) 
#define smSetMARK(pi__) TYPESet( (&p[pi__]), TYPE_MARK )
#define smResetMARK(pi__) TYPEReset( (&p[pi__]), TYPE_MARK )
#endif 

//#define NSI_DEBUG(xxx) xxx
#define NSI_DEBUG(xxx) 

void smInitList(SMX smx, int nSmooth) {
    /*
    ** Allocate Nearest-Neighbor List.
    */
    smx->nListSize = (smx->nSmooth > 1 ? smx->nSmooth : 2);
    smx->nnList = (NN *)malloc(smx->nListSize*sizeof(NN));
    assert(smx->nnList != NULL);
    smx->pbRelease = (int *)malloc(smx->nListSize*sizeof(int));
    assert(smx->pbRelease != NULL);
    /*
    ** Allocate priority queue.
    */
    smx->pq = (PQ *)malloc(nSmooth*sizeof(PQ));
    assert(smx->pq != NULL);
    PQ_INIT(smx->pq,nSmooth);
    /*
    ** Set up "in priority queue" hash table stuff.
    */
    smx->nHash = (int)ceil(nSmooth/PQ_LOAD_FACTOR);
    smx->pqHash = (PQ **)malloc(smx->nHash*sizeof(PQ *));
    assert(smx->pqHash != NULL);
    }


void smFinishList(SMX smx) {
    free(smx->pqHash);
    free(smx->pq);
    free(smx->nnList);
    free(smx->pbRelease);
    }


int smInitialize(SMX *psmx,PKD pkd,SMF *smf,int nSmooth,int bPeriodic,
    int bSymmetric,int iSmoothType,int bSmooth, double dfBall2OverSoft2 )
    {
    SMX smx;
    void (*initParticle)(void *) = NULL;
    void (*initTreeParticle)(void *) = NULL;
    void (*init)(void *) = NULL;
    void (*comb)(void *,void *) = NULL;
    int pi;
    int nTree;

    smx = (SMX)malloc(sizeof(struct smContext));
    assert(smx != NULL);
    smx->pkd = smf->pkd = pkd;
    smx->nSmooth = nSmooth;
    smx->nSmoothMax = 10*nSmooth;
    smx->fBall2InnerFrac = 0.5;
    smx->nSmoothInner = pow(smx->fBall2InnerFrac,1.5)*smx->nSmooth;
    smx->bPeriodic = bPeriodic;
    if(smx->bPeriodic) {
        smx->lx = smx->pkd->fPeriod[0];
        smx->ly = smx->pkd->fPeriod[1];
        smx->lz = smx->pkd->fPeriod[2];
        }
    else {
        smx->lx = FLOAT_MAXVAL;
        smx->ly = FLOAT_MAXVAL;
        smx->lz = FLOAT_MAXVAL;
        }
    smx->dfBall2OverSoft2 = dfBall2OverSoft2;
    smx->iLowhFix = ( dfBall2OverSoft2 > 0.0 ? LOWHFIX_HOVERSOFT : LOWHFIX_NONE );
    smx->bUseBallMax = 1;
#ifdef SLIDING_PATCH
    smx->dTime = smf->dTime;
    smx->PP = &smf->PP;
#endif

    switch (iSmoothType) {
    case SMX_NULL:
        smx->fcnSmooth = NullSmooth;
        initParticle = NULL; /* Original Particle */
        initTreeParticle = NULL; /* Original Particle */
        init = NULL; /* Cached copies */
        comb = NULL;
        smx->fcnPost = NULL;
        smx->bUseBallMax = 0;
        break;
    case SMX_DENSITY:
        smx->fcnSmooth = bSymmetric?DensitySym:Density;
        initParticle = initDensity; /* Original Particle */
        initTreeParticle = NULL; /* Original Particle */
        init = initDensity; /* Cached copies */
        comb = combDensity;
        smx->fcnPost = NULL;
        break;
    case SMX_DENSITYTMP:
        smx->fcnSmooth = bSymmetric?DensityTmpSym:DensityTmp;
        initParticle = initDensityTmp; /* Original Particle */
        initTreeParticle = NULL; /* Original Particle */
        init = initDensityTmp; /* Cached copies */
        comb = combDensityTmp;
        smx->fcnPost = NULL;
        break;
    case SMX_MARKDENSITY:
        smx->fcnSmooth = bSymmetric?MarkDensitySym:MarkDensity;
        initParticle = initParticleMarkDensity; /* Original Particle */
        initTreeParticle = NULL; /* Original Particle */
        init = initMarkDensity; /* Cached copies */
        comb = combMarkDensity;
        smx->fcnPost = NULL;
        break;
    case SMX_MARKIIDENSITY:
        smx->fcnSmooth = bSymmetric?MarkIIDensitySym:MarkIIDensity;
        initParticle = initParticleMarkIIDensity; /* Original Particle */
        initTreeParticle = NULL; /* Original Particle */
        init = initMarkIIDensity; /* Cached copies */
        comb = combMarkIIDensity;
        smx->fcnPost = NULL;
        break;
    case SMX_MARK:
        smx->fcnSmooth = NULL;
        initParticle = NULL;
        initTreeParticle = NULL;
        init = initMark;
        comb = combMark;
        smx->fcnPost = NULL;
        break;
    case SMX_DT:
        smx->fcnSmooth = NULL;
        initParticle = NULL;
        initTreeParticle = NULL;
        init = NULL;
        comb = NULL;
        smx->fcnPost = NULL;
        break;
    case SMX_DELTAACCEL:
        smx->fcnSmooth = DeltaAccel;
        initParticle = NULL; /* Original Particle */
        initTreeParticle = NULL; /* Original Particle */
        init = NULL; /* Cached copies */
        comb = combDeltaAccel;
        smx->fcnPost = NULL;
        break;
    case SMX_SINKACCRETETEST:
        smx->fcnSmooth = SinkAccreteTest;
        initParticle = initSinkAccreteTest; /* Original Particle */
        initTreeParticle = initSinkAccreteTest; /* Original Particle */
        init = initSinkAccreteTest; /* Cached copies */
        comb = combSinkAccreteTest;
        smx->fcnPost = NULL;
        smx->iLowhFix = LOWHFIX_SINKRADIUS;
        smx->bUseBallMax = 0;
        break;
    case SMX_SINKACCRETE:
        smx->fcnSmooth = SinkAccrete;
        initParticle = NULL; /* Original Particle */
        initTreeParticle = NULL; /* Original Particle */
        init = initSinkAccrete; /* Cached copies */
        comb = combSinkAccrete;
        smx->fcnPost = NULL;
        smx->iLowhFix = LOWHFIX_SINKRADIUS;
        smx->bUseBallMax = 0;
        break;
    case SMX_SINKINGAVERAGE:
        smx->fcnSmooth = SinkingAverage;
        initParticle = NULL; /* Original Particle */
        initTreeParticle = NULL; /* Original Particle */
        init = NULL; /* Cached copies */
        comb = NULL;
        smx->fcnPost = NULL;
        break;
    case SMX_SINKINGFORCESHARE:
        smx->fcnSmooth = SinkingForceShare;
        initParticle = NULL; /* Original Particle */
        initTreeParticle = NULL; /* Original Particle */
        init = initSinkingForceShare; /* Cached copies */
        comb = combSinkingForceShare;
        smx->fcnPost = NULL;
        smx->iLowhFix = LOWHFIX_SINKRADIUS_BUFF;
        smx->bUseBallMax = 0;
        break;
    case SMX_BHDENSITY:
        smx->fcnSmooth = BHSinkDensity;
        initParticle = initBHSinkDensity; /* Original Particle */
        initTreeParticle = initTreeParticleBHSinkDensity; /* Original Particle */
        init = initBHSinkDensity; /* Cached copies */
        comb = combBHSinkDensity;
        smx->fcnPost = NULL;
        break;
    case SMX_BHSINKACCRETE:
        smx->fcnSmooth = BHSinkAccrete;
        initParticle = NULL; /* Original Particle */
        initTreeParticle = initTreeParticleBHSinkAccrete; /* Original Particle */
        init = initBHSinkAccrete; /* Cached copies */
        comb = combBHSinkAccrete;
        smx->fcnPost = postBHSinkAccrete;
        break;
    case SMX_BHSINKIDENTIFY:
        smx->fcnSmooth = BHSinkIdentify;
        initParticle = NULL; /* Original Particle */
        initTreeParticle = NULL; /* Original Particle */
        init = initBHSinkIdentify; /* Cached copies */
        comb = combBHSinkIdentify;
        smx->fcnPost = NULL;
        break;
    case SMX_BHSINKMERGE:
        smx->fcnSmooth = BHSinkMerge;
        initParticle = NULL; /* Original Particle */
        initTreeParticle = NULL; /* Original Particle */
        init = initBHSinkMerge; /* Cached copies */
        comb = combBHSinkMerge;
        smx->fcnPost = NULL;
        break;
    case SMX_SINKFORMTEST:
        smx->fcnSmooth = SinkFormTest;
        initParticle = initSinkFormTest; /* Original Particle */
        initTreeParticle = initSinkFormTest; /* Original Particle */
        init = initSinkFormTest; /* Cached copies */
        comb = combSinkFormTest;
        smx->fcnPost = NULL;
        smx->iLowhFix = LOWHFIX_SINKRADIUS;
        smx->bUseBallMax = 0;
        break;
    case SMX_SINKFORM:
        smx->fcnSmooth = SinkForm;
        initParticle = NULL; /* Original Particle */
        initTreeParticle = NULL; /* Original Particle */
        init = initSinkForm; /* Cached copies */
        comb = combSinkForm;
        smx->fcnPost = NULL;
        smx->iLowhFix = LOWHFIX_SINKRADIUS;
        smx->bUseBallMax = 0;
        break;
    case SMX_SINKMERGETEST:
        smx->fcnSmooth = SinkMergeTest;
        initParticle = initSinkMergeTest; /* Original Particle */
        initTreeParticle = initSinkMergeTest; /* Original Particle */
        init = initSinkMergeTest; /* Cached copies */
        comb = combSinkMergeTest;
        smx->fcnPost = NULL;
        break;
    case SMX_SINKMERGE:
        smx->fcnSmooth = SinkMerge;
        initParticle = NULL; /* Original Particle */
        initTreeParticle = NULL; /* Original Particle */
        init = initSinkMerge; /* Cached copies */
        comb = combSinkMerge;
        smx->fcnPost = NULL;
        break;
#ifdef SUPERCOOL
    case SMX_MEANVEL:
        smx->fcnSmooth = bSymmetric?MeanVelSym:MeanVel;
        initParticle = initMeanVel;
        initTreeParticle = NULL;
        init = initMeanVel;
        comb = combMeanVel;
        smx->fcnPost = NULL;
        break;
#endif
#ifdef GASOLINE
    case SMX_SPHPRESSURETERMS:
        smx->fcnSmooth = bSymmetric?SphPressureTermsSym:SphPressureTerms;
        initParticle = initSphPressureTermsParticle; /* Original Particle */
        initTreeParticle = NULL;
        init = initSphPressureTerms; /* Cached copies */
        comb = combSphPressureTerms;
        smx->fcnPost = NULL;
        break;
    case SMX_DENDVDX:
        smx->fcnSmooth = DenDVDX;
        initParticle = NULL;
        initTreeParticle = NULL;
        init = initDenDVDX;
        comb = combDenDVDX;
        smx->fcnPost = postDenDVDX;
        break;
    case SMX_SURFACENORMAL:
        printf("Surface Normal\n");
        smx->fcnSmooth = SurfaceNormal;
        initParticle = NULL;
        initTreeParticle = NULL;
        init = initSurfaceNormal;
        comb = combSurfaceNormal;
        smx->fcnPost = NULL;
        break;
    case SMX_SURFACEAREA:
        printf("Surface Area\n");
        smx->fcnSmooth = SurfaceArea;
        initParticle = NULL;
        initTreeParticle = NULL;
        init = initSurfaceArea;
        comb = combSurfaceArea;
        smx->fcnPost = NULL;
        break;
    case SMX_SMOOTHBSW:
        smx->fcnSmooth = SmoothBSw;
        initParticle = NULL;
        initTreeParticle = NULL;
        init = initSmoothBSw;
        comb = combSmoothBSw;
        smx->fcnPost = NULL;
        break;
    case SMX_DIVVORT:
        printf("Div vort\n");
        smx->fcnSmooth = bSymmetric?DivVortSym:DivVort;
        initParticle = initDivVort;
        initTreeParticle = NULL;
        init = initDivVort;
        comb = combDivVort;
        smx->fcnPost = NULL;
        break;
    case SMX_SHOCKTRACK:
        smx->fcnSmooth = bSymmetric?ShockTrackSym:ShockTrack;
        initParticle = initShockTrack;
        initTreeParticle = NULL;
        init = initShockTrack;
        comb = combShockTrack;
        smx->fcnPost = NULL;
        break;
    case SMX_HKPRESSURETERMS:
        smx->fcnSmooth = bSymmetric?HKPressureTermsSym:HKPressureTerms;
        initParticle = initHKPressureTermsParticle; /* Original Particle */
        initTreeParticle = NULL;
        init = initHKPressureTerms; /* Cached copies */
        comb = combHKPressureTerms;
        smx->fcnPost = NULL;
        break;
    case SMX_SPHPRESSURE:
        smx->fcnSmooth = bSymmetric?SphPressureSym:SphPressure;
        initParticle = initSphPressureParticle; /* Original Particle */
        init = initSphPressure; /* Cached copies */
        comb = combSphPressure;
        smx->fcnPost = postSphPressure;
        break;
    case SMX_SPHVISCOSITY:
        smx->fcnSmooth = bSymmetric?SphViscositySym:SphViscosity;
        initParticle = initSphViscosityParticle; /* Original Particle */
        initTreeParticle = NULL;
        init = initSphViscosity; /* Cached copies */
        comb = combSphViscosity;
        smx->fcnPost = NULL;
        break;
    case SMX_HKVISCOSITY:
        smx->fcnSmooth = bSymmetric?HKViscositySym:HKViscosity;
        initParticle = initHKViscosityParticle; /* Original Particle */
        initTreeParticle = NULL;
        init = initHKViscosity; /* Cached copies */
        comb = combHKViscosity;
        smx->fcnPost = NULL;
        break;
#ifdef STARFORM
    case SMX_STARCLUSTERFORM:
        smx->fcnSmooth = StarClusterForm;
        initParticle = NULL; /* Original Particle */
        initTreeParticle = NULL; /* Original Particle */
        init = initStarClusterForm; /* Cached copies */
        comb = combStarClusterForm;
        smx->fcnPost = NULL;
        smx->iLowhFix = LOWHFIX_SINKRADIUS;
        smx->bUseBallMax = 0;
        break;
    case SMX_DELETE_GAS:
        assert(bSymmetric == 0);
        smx->fcnSmooth = DeleteGas;
        initParticle = NULL;
        initTreeParticle = NULL;
        init = NULL;
        comb = NULL;
        smx->fcnPost = NULL;
        smx->bUseBallMax = 0;
        break;
    case SMX_DIST_DELETED_GAS:
        assert(bSymmetric != 0);
        smx->fcnSmooth = DistDeletedGas;
        initParticle = NULL;
        initTreeParticle = NULL;
        init = initDistDeletedGas;
        comb = combDistDeletedGas;
        smx->fcnPost = NULL;
        smx->bUseBallMax = 0;
        break;
    case SMX_PROMOTE_TO_HOT_GAS:
        assert(bSymmetric != 0);
        smx->fcnSmooth = PromoteToHotGas;
        initParticle = initPromoteToHotGas;
        initTreeParticle = initPromoteToHotGas;
        init = initPromoteToHotGas;
        comb = combPromoteToHotGas;
        smx->fcnPost = NULL;
        smx->bUseBallMax = 0;
        break;
    case SMX_SHARE_WITH_HOT_GAS:
        assert(bSymmetric != 0);
        smx->fcnSmooth = ShareWithHotGas;
        initParticle = NULL;
        initTreeParticle = NULL;
        init = initShareWithHotGas;
        comb = combShareWithHotGas;
        smx->fcnPost = NULL;
        smx->bUseBallMax = 0;
        break;
    case SMX_DIST_FB_ENERGY:
        assert(bSymmetric != 0);
        smx->fcnSmooth = DistFBEnergy;
        initParticle = NULL;
        initTreeParticle = initTreeParticleDistFBEnergy;
        init = initDistFBEnergy;
        comb = combDistFBEnergy;
        smx->fcnPost = postDistFBEnergy;
        smx->bUseBallMax = 0;
        break;
#endif

#endif         
#ifdef COLLISIONS
    case SMX_REJECTS:
        assert(bSymmetric != 0);
        smx->fcnSmooth = FindRejects;
        initParticle = initFindRejects;
        initTreeParticle = NULL;
        init = initFindRejects;
        comb = combFindRejects;
        smx->fcnPost = NULL;
        smx->bUseBallMax = 0;
        break;
    case SMX_COLLISION:
        assert(bSymmetric == 0);
        smx->fcnSmooth = CheckForCollision;
        initParticle = NULL;
        initTreeParticle = NULL;
        init = NULL;
        comb = NULL;
        smx->fcnPost = NULL;
        smx->bUseBallMax = 0;
        break;
    case SMX_FINDBINARY:
        smx->fcnSmooth = FindBinary;
        initParticle = NULL;
        init = NULL;
        comb = NULL;
        smx->fcnPost = NULL;
        break;
#endif /* COLLISIONS */
#ifdef SLIDING_PATCH
    case SMX_FIND_OVERLAPS:
        assert(bSymmetric != 0);
        smx->fcnSmooth = FindOverlaps;
        initParticle = initFindOverlaps;
        initTreeParticle = NULL;
        init = initFindOverlaps;
        comb = combFindOverlaps;
        smx->fcnPost = NULL;
        break;
#endif /* SLIDING_PATCH */
    default:
        assert(0);
        }
    /*
    ** Initialize the ACTIVE particles in the tree.
    ** There are other particles in the tree -- just not active.
    */
    nTree = pkd->kdNodes[pkd->iRoot].pUpper + 1;
    for (pi=0;pi<nTree;++pi) {
        if (TYPEQuerySMOOTHACTIVE(&(pkd->pStore[pi]))) {
            TYPEReset( &(pkd->pStore[pi]),TYPE_SMOOTHDONE );
            /*                  if (bSmooth) pkd->pStore[pi].fBall2 = -1.0;*/
            if (initParticle != NULL) {
                initParticle(&pkd->pStore[pi]);
                }
            }
        else if (TYPEQueryTREEACTIVE(&(pkd->pStore[pi]))) {
            /*TYPEReset( &(pkd->pStore[pi]),TYPE_SMOOTHDONE );*/
            /*                  if (bSmooth) pkd->pStore[pi].fBall2 = -1.0;*/
            if (initTreeParticle != NULL) {
                initTreeParticle(&pkd->pStore[pi]);
                }
            }
        }
    /*
    ** Start particle caching space (cell cache is already active).
    */
    if (bSymmetric) {
        assert(init != NULL);
        assert(comb != NULL);
        mdlCOcache(pkd->mdl,CID_PARTICLE,pkd->pStore,sizeof(PARTICLE),
            nTree,init,comb);
        }
    else {
        mdlROcache(pkd->mdl,CID_PARTICLE,pkd->pStore,sizeof(PARTICLE),
            nTree);
        }
    /*
    ** Allocate mark array.
    */
#ifdef MARK
    if(pkdLocal(pkd) == 0)
        smx->piMark = (int *)malloc(sizeof(int));
    else
        smx->piMark = (int *)malloc(pkdLocal(pkd)*sizeof(int));
    assert(smx->piMark != NULL);
#endif
    smInitList(smx,nSmooth);

    *psmx = smx;        
    return(1);
    }


void smFinish(SMX smx,SMF *smf, CASTAT *pcs)
    {
    PKD pkd = smx->pkd;
    int pi;

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
    mdlFinishCache(smx->pkd->mdl,CID_PARTICLE);
    /*
    ** Now do any post calculations, these ususlly involve some sort of
    ** normalizations of the smoothed quantities, usually division by
    ** the local density! Do NOT put kernel normalizations in here as
    ** these do not depend purely on local properties in the case of
    ** "Gather-Scatter" kernel.
    */
    if (smx->fcnPost != NULL) {
        for (pi=0;pi<pkd->nTreeActive;++pi) {
            smx->fcnPost(&pkd->pStore[pi],smf); 
            }
        }
    /*
    ** Free up context storage.
    */
#ifdef MARK
    free(smx->piMark);
#endif
    smFinishList(smx);
    free(smx);
    }

void smLargefBallCheck(SMX smx,PARTICLE *p,FLOAT lx, FLOAT ly, FLOAT lz) {
#if (defined(SLIDING_PATCH) && INTERNAL_WARNINGS)
    /*
    ** For periodic boundary conditions, make sure search ball has not
    ** exceeded half the spatial period. If it has, this probably means
    ** nSmooth is too large. HOWEVER, in the case, e.g., of a periodic
    ** patch with no z boundary, a particle could legitimately get far
    ** enough from the plane to generate an error. In this case only a
    ** warning is issued (once per processor, if warnings are enabled).
    **
    ** DCR 3/15/05: now want to relax the x boundary condition too!
    ** (to look at thin rings)...  so let's just make it general.
    */
    
    if ((lx < FLOAT_MAXVAL && p->fBall2 >= 0.25*lx*lx) ||
        (ly < FLOAT_MAXVAL && p->fBall2 >= 0.25*ly*ly) ||
        (lz < FLOAT_MAXVAL && p->fBall2 >= 0.25*lz*lz)) {
        static int bGiveWarning = 1;
        if (bGiveWarning) {
            (void) fprintf(stderr,"WARNING: Large search ball (iOrder = %i)..."
                "lx = %g ly = %g lz = %g fBall = %g x = %g y = %g z = %g\n",
                p->iOrder,lx,ly,lz,sqrt(p->fBall2),p->r[0],p->r[1],p->r[2]);
#if (INTERNAL_WARNINGS_ONCE)
            bGiveWarning = 0;
#endif
            }
        }
#else /* SLIDING_PATCH && INTERNAL_WARNINGS */
#endif
    }


PQ *smBallSearch(SMX smx,PQ *pq,FLOAT *ri,int *cpStart)
    {
    MDL mdl = smx->pkd->mdl;
    KDN *c = smx->pkd->kdNodes;
    PARTICLE *p = smx->pkd->pStore;
    int cell,cp,ct,pj,pUpper,idSelf,bPeriodic;
    FLOAT fBall2,fDist2,dx,dy,dz,lx,ly,lz,sx,sy,sz,x,y,z;

    x = ri[0];    y = ri[1];    z = ri[2];
    lx = smx->lx; ly = smx->ly; lz = smx->lz;
    idSelf = smx->pkd->idSelf;
    fBall2 = pq->fKey;
    cell = ROOT;
    bPeriodic = 0;
    /*
    ** First find the "local" Bucket.
    ** This could be the closest to ri[3]! 
    ** Warning: in parallel ri[3] SHOULD be contained in the LOCAL DOMAIN!
    */
    while (c[cell].iDim >= 0) {
        if (ri[c[cell].iDim] < c[cell].fSplit) cell = LOWER(cell);
        else cell = UPPER(cell);
        }
    pUpper = c[cell].pUpper;
    for (pj=c[cell].pLower;pj<=pUpper;++pj) {
        dx = x - p[pj].r[0];
        dy = y - p[pj].r[1];
        dz = z - p[pj].r[2];
        fDist2 = dx*dx + dy*dy + dz*dz;
        if (fDist2 < fBall2) {
            *cpStart = cell;
            if (smTestMARK(pj)) continue; /* Don't include particle in queue twice */
            if (pq->id == idSelf) smResetMARK(pq->p); /* unmark local particle as in queue */
            else {
                mdlRelease(mdl,CID_PARTICLE,pq->pPart); /* pq->pPart is non-local -- release */
                PQ_HASHDEL(smx->pqHash,smx->nHash,pq); /* remove from pqHash that marks non-local nbrs */
                pq->id = idSelf;
                }
            smSetMARK(pj); /* Mark new local particle as in the queue */
            pq->fKey = fDist2;
            pq->dx = dx;
            pq->dy = dy;
            pq->dz = dz;
            pq->p = pj;
            pq->pPart = &p[pj];
            pq->ax = 0.0;
            pq->ay = 0.0;
            pq->az = 0.0;
            PQ_REPLACE(pq);
            fBall2 = pq->fKey;
            }
        }
    while (cell != ROOT) {
        cp = SIBLING(cell);
        ct = cp;
        SETNEXT(ct);
        while (1) {
            INTERSECT(&c[cp],fBall2,lx,ly,lz,x,y,z,sx,sy,sz,bPeriodic,
                GetNext_1);
            if (c[cp].iDim >= 0) {
                cp = LOWER(cp);
                continue;
                }
            else {
                pUpper = c[cp].pUpper;
                for (pj=c[cp].pLower;pj<=pUpper;++pj) {
                    dx = sx - p[pj].r[0];
                    dy = sy - p[pj].r[1];
                    dz = sz - p[pj].r[2];
                    fDist2 = dx*dx + dy*dy + dz*dz;
                    if (fDist2 < fBall2) {
                        *cpStart = PARENT(cell);
                        if (smTestMARK(pj)) continue;
                        if (pq->id == idSelf) smResetMARK(pq->p);
                        else {
                            mdlRelease(mdl,CID_PARTICLE,pq->pPart);
                            PQ_HASHDEL(smx->pqHash,smx->nHash,pq);
                            pq->id = idSelf;
                            }
                        smSetMARK(pj);
                        pq->fKey = fDist2;
                        pq->dx = dx;
                        pq->dy = dy;
                        pq->dz = dz;
                        pq->p = pj;
                        pq->pPart = &p[pj];
                        pq->ax = sx - x;
                        pq->ay = sy - y;
                        pq->az = sz - z;
                        PQ_REPLACE(pq);
                        fBall2 = pq->fKey;
                        }
                    }
                }
        GetNext_1:
            SETNEXT(cp);
            if (cp == ct) break;
            }
        cell = PARENT(cell);
        }
    if (bPeriodic) *cpStart = 0;
    return(pq);
    }


PQ *smBallSearchNP(SMX smx,PQ *pq,FLOAT *ri,int *cpStart)
    {
    MDL mdl = smx->pkd->mdl;
    KDN *c = smx->pkd->kdNodes;
    PARTICLE *p = smx->pkd->pStore;
    int cell,cp,ct,pj,pUpper,idSelf;
    FLOAT fBall2,fDist2,dx,dy,dz,x,y,z;

    x = ri[0];    y = ri[1];    z = ri[2];
    idSelf = smx->pkd->idSelf;
    fBall2 = pq->fKey;
    cell = ROOT;
    /*
    ** First find the "local" Bucket.
    ** This could mearly be the closest to ri[3]! 
    ** Warning: in parallel ri[3] SHOULD be contained in the LOCAL DOMAIN!
    */
    while (c[cell].iDim >= 0) {
        if (ri[c[cell].iDim] < c[cell].fSplit) cell = LOWER(cell);
        else cell = UPPER(cell);
        }
    pUpper = c[cell].pUpper;
    for (pj=c[cell].pLower;pj<=pUpper;++pj) {
        dx = x - p[pj].r[0];
        dy = y - p[pj].r[1];
        dz = z - p[pj].r[2];
        fDist2 = dx*dx + dy*dy + dz*dz;
        if (fDist2 < fBall2) {
            *cpStart = cell;
            if (smTestMARK(pj)) continue;
            if (pq->id == idSelf) smResetMARK(pq->p);
            else {
                mdlRelease(mdl,CID_PARTICLE,pq->pPart);
                PQ_HASHDEL(smx->pqHash,smx->nHash,pq);
                pq->id = idSelf;
                }
            smSetMARK(pj);
            pq->fKey = fDist2;
            pq->dx = dx;
            pq->dy = dy;
            pq->dz = dz;
            pq->p = pj;
            pq->pPart = &p[pj];
            PQ_REPLACE(pq);
            fBall2 = pq->fKey;
            }
        }
    while (cell != ROOT) {
        cp = SIBLING(cell);
        ct = cp;
        SETNEXT(ct);
        while (1) {
            INTERSECTNP(&c[cp],fBall2,x,y,z,GetNext_1);
            if (c[cp].iDim >= 0) {
                cp = LOWER(cp);
                continue;
                }
            else {
                pUpper = c[cp].pUpper;
                for (pj=c[cp].pLower;pj<=pUpper;++pj) {
                    dx = x - p[pj].r[0];
                    dy = y - p[pj].r[1];
                    dz = z - p[pj].r[2];
                    fDist2 = dx*dx + dy*dy + dz*dz;
                    if (fDist2 < fBall2) {
                        *cpStart = PARENT(cell);
                        if (smTestMARK(pj)) continue;
                        if (pq->id == idSelf) smResetMARK(pq->p);
                        else {
                            mdlRelease(mdl,CID_PARTICLE,pq->pPart);
                            PQ_HASHDEL(smx->pqHash,smx->nHash,pq);
                            pq->id = idSelf;
                            }
                        smSetMARK(pj);
                        pq->fKey = fDist2;
                        pq->dx = dx;
                        pq->dy = dy;
                        pq->dz = dz;
                        pq->p = pj;
                        pq->pPart = &p[pj];
                        PQ_REPLACE(pq);
                        fBall2 = pq->fKey;
                        }
                    }
                }
        GetNext_1:
            SETNEXT(cp);
            if (cp == ct) break;
            }
        cell = PARENT(cell);
        }
    return(pq);
    }


/*
 * I don't increase the PQ size here: I assume it isn't in use if we
 * need this function.
 */
void smGrowList(SMX smx)
    {
    smx->nListSize *= 1.5;
    
    smx->nnList = (NN *) realloc(smx->nnList,smx->nListSize*sizeof(NN));
    assert(smx->nnList != NULL);
    smx->pbRelease =
        (int *) realloc(smx->pbRelease,smx->nListSize*sizeof(int));
    assert(smx->pbRelease != NULL);
    }

int smBallGather(SMX smx,FLOAT fBall2,FLOAT *ri)
    {
    KDN *c = smx->pkd->kdNodes;
    PARTICLE *p = smx->pkd->pStore;
    int pj,nCnt,cp,pUpper;
    FLOAT dx,dy,dz,x,y,z,lx,ly,lz,sx,sy,sz,fDist2;
    int iDum;

    x = ri[0];    y = ri[1];    z = ri[2];
    lx = smx->lx; ly = smx->ly; lz = smx->lz;
    nCnt = 0;
    cp = ROOT;
    if(c[cp].pLower > c[cp].pUpper) {
        return(0);              /* Empty Tree */
        }
            
    while (1) {
        INTERSECT(&c[cp],fBall2,lx,ly,lz,x,y,z,sx,sy,sz,iDum,GetNextCell);
        /*
        ** We have an intersection to test.
        */
        if (c[cp].iDim >= 0) {
            cp = LOWER(cp);
            continue;
            }
        else {
            pUpper = c[cp].pUpper;
            for (pj=c[cp].pLower;pj<=pUpper;++pj) {
                dx = sx - p[pj].r[0];
                dy = sy - p[pj].r[1];
                dz = sz - p[pj].r[2];
                fDist2 = dx*dx + dy*dy + dz*dz;
                if (fDist2 <= fBall2) {
                    if(nCnt >= smx->nListSize)
                        smGrowList(smx);
                    smx->nnList[nCnt].fDist2 = fDist2;
                    smx->nnList[nCnt].dx = dx;
                    smx->nnList[nCnt].dy = dy;
                    smx->nnList[nCnt].dz = dz;
                    smx->nnList[nCnt].pPart = &p[pj];
                    smx->nnList[nCnt].iIndex = pj;
                    smx->nnList[nCnt].iPid = smx->pkd->idSelf;
                    smx->pbRelease[nCnt++] = 0;
                    }
                }
            }
    GetNextCell:
        SETNEXT(cp);
        if (cp == ROOT) break;
        }
    assert(nCnt <= smx->nListSize);
    return(nCnt);
    }


int smBallGatherNP(SMX smx,FLOAT fBall2,FLOAT *ri,int cp)
    {
    KDN *c = smx->pkd->kdNodes;
    PARTICLE *p = smx->pkd->pStore;
    int pj,nCnt,pUpper;
    FLOAT dx,dy,dz,x,y,z,fDist2;

    x = ri[0];    y = ri[1];    z = ri[2];
    nCnt = 0;
    if(c[cp].pLower > c[cp].pUpper) {
        return(0);              /* Empty Tree */
        }
    while (1) {
        INTERSECTNP(&c[cp],fBall2,x,y,z,GetNextCell);
        /*
        ** We have an intersection to test.
        */
        if (c[cp].iDim >= 0) {
            cp = LOWER(cp);
            continue;
            }
        else {
            pUpper = c[cp].pUpper;
            for (pj=c[cp].pLower;pj<=pUpper;++pj) {
                dx = x - p[pj].r[0];
                dy = y - p[pj].r[1];
                dz = z - p[pj].r[2];
                fDist2 = dx*dx + dy*dy + dz*dz;
                if (fDist2 <= fBall2) {
                    if(nCnt >= smx->nListSize)
                        smGrowList(smx);
                    smx->nnList[nCnt].fDist2 = fDist2;
                    smx->nnList[nCnt].dx = dx;
                    smx->nnList[nCnt].dy = dy;
                    smx->nnList[nCnt].dz = dz;
                    smx->nnList[nCnt].pPart = &p[pj];
                    smx->nnList[nCnt].iIndex = pj;
                    smx->nnList[nCnt].iPid = smx->pkd->idSelf;
                    smx->pbRelease[nCnt++] = 0;
                    }
                }
            }
    GetNextCell:
        SETNEXT(cp);
        if (cp == ROOT) break;
        }
    assert(nCnt <= smx->nListSize);
    return(nCnt);
    }


void smBallScatter(SMX smx,FLOAT *ri,int iMarkType)
    {
    KDN *c = smx->pkd->kdNodes;
    PARTICLE *p = smx->pkd->pStore;
    int pj,cp,pUpper;
    FLOAT dx,dy,dz,x,y,z,lx,ly,lz,sx,sy,sz,fDist2;
    int iDum;

    /*mdlDiag(smx->pkd->mdl, "smBallScatter:Start\n" );*/
    x = ri[0];    y = ri[1];    z = ri[2];
    lx = smx->lx; ly = smx->ly; lz = smx->lz;

    cp = ROOT;
    while (1) {
        /*mdlDiag(smx->pkd->mdl, "smINTERSECTSCATTER:Before\n" );*/
        INTERSECTSCATTER(&c[cp],lx,ly,lz,x,y,z,sx,sy,sz,iDum,GetNextCell);
        /*mdlDiag(smx->pkd->mdl, "smINTERSECTSCATTER:After\n" );*/
        /*
        ** We have an intersection to test.
        */
        if (c[cp].iDim >= 0) {
            cp = LOWER(cp);
            continue;
            }
        else {
            pUpper = c[cp].pUpper;
            for (pj=c[cp].pLower;pj<=pUpper;++pj) {
                dx = sx - p[pj].r[0];
                dy = sy - p[pj].r[1];
                dz = sz - p[pj].r[2];
                fDist2 = dx*dx + dy*dy + dz*dz;
                if (fDist2 <= p[pj].fBallMax*p[pj].fBallMax) {
                    TYPESet(&(p[pj]),iMarkType);
                    }
                }
            }
    GetNextCell:
        SETNEXT(cp);
        if (cp == ROOT) break;
        }

    /*mdlDiag(smx->pkd->mdl, "smBallScatter:End\n" );*/
    }


void smBallScatterNP(SMX smx,FLOAT *ri,int iMarkType,int cp)
    {
    KDN *c = smx->pkd->kdNodes;
    PARTICLE *p = smx->pkd->pStore;
    int pj,pUpper;
    FLOAT dx,dy,dz,x,y,z,fDist2;

    x = ri[0];    y = ri[1];    z = ri[2];

    while (1) {
        INTERSECTSCATTERNP(&c[cp],x,y,z,GetNextCell);
        /*
        ** We have an intersection to test.
        */
        if (c[cp].iDim >= 0) {
            cp = LOWER(cp);
            continue;
            }
        else {
            pUpper = c[cp].pUpper;
            for (pj=c[cp].pLower;pj<=pUpper;++pj) {
                dx = x - p[pj].r[0];
                dy = y - p[pj].r[1];
                dz = z - p[pj].r[2];
                fDist2 = dx*dx + dy*dy + dz*dz;
                if (fDist2 <= p[pj].fBallMax*p[pj].fBallMax) {
                    TYPESet(&(p[pj]),iMarkType);
                    }
                }
            }
    GetNextCell:
        SETNEXT(cp);
        if (cp == ROOT) break;
        }
    }

#define INVDTOPENANGLE2 1

void smDtBallNP(SMX smx,FLOAT *ri,FLOAT *vi,FLOAT *pdt2,FLOAT fBall2,int cp)
    { 
#ifdef GASOLINE
    KDN *c = smx->pkd->kdNodes,*pkdn;
    PARTICLE *p = smx->pkd->pStore;
    int pj,pUpper;
    FLOAT dt2,dtEst2,dx,dy,dz,dvx,dvy,dvz,x,y,z,vx,vy,vz,fDist2,vdotr;

    x = ri[0];     y = ri[1];    z = ri[2];
    vx = vi[0];    vy = vi[1];   vz = vi[2];
    dt2 = (*pdt2);

    while (1) {
        pkdn = &c[cp];
        DTINTERSECTNP(pkdn,dt2,fDist2,x,y,z,vx,vy,vz,GetNextCell);
        /*
        ** We have an intersection to test.
        */
        if (fDist2 > INVDTOPENANGLE2*pkdn->bndDt.drMax2) {
            /* Particle Cell - passed Opening Angle Test */
            DTESTIMATOR(&c[cp],dtEst2,x,y,z,vx,vy,vz);
            dtEst2 = dtEst2*dtEst2;
            if (dtEst2 < dt2) dt2 = dtEst2;
            goto GetNextCell;
            }
        if (c[cp].iDim >= 0) {
            cp = LOWER(cp);
            continue;
            }
        else {
            pUpper = c[cp].pUpper;
            for (pj=c[cp].pLower;pj<=pUpper;++pj) {
                dx = x - p[pj].r[0];
                dy = y - p[pj].r[1];
                dz = z - p[pj].r[2];
                fDist2 = (dx*dx + dy*dy + dz*dz); 
                if (fDist2 < fBall2) continue;
                dvx = vx - p[pj].vPred[0];
                dvy = vy - p[pj].vPred[1];
                dvz = vz - p[pj].vPred[2];
                if ((vdotr = dvx*dx + dvy*dy + dvz*dz) < 0) {
                    double fDist = sqrt(fDist2);
                    dtEst2 = fDist/(p[pj].c-vdotr/fDist);
                    dtEst2 = dtEst2*dtEst2;
                    if (dtEst2 < dt2) dt2 = dtEst2;
                    }
                else {
                    dtEst2 = fDist2/(p[pj].c*p[pj].c);
                    if (dtEst2 < dt2) dt2 = dtEst2;
                    }
                }
            }
    GetNextCell:
        SETNEXT(cp);
        if (cp == ROOT) break;
        }
    *pdt2 = dt2;
#endif
    }


void smDtBall(SMX smx,FLOAT *ri,FLOAT *vi,FLOAT *pdt2,FLOAT fBall2)
    {
#ifdef GASOLINE
    KDN *c = smx->pkd->kdNodes, *pkdn;
    PARTICLE *p = smx->pkd->pStore;
    int pj,cp,pUpper;
    FLOAT dt2,dtEst2,dx,dy,dz,dvx,dvy,dvz;
    FLOAT x,y,z,lx,ly,lz,sx,sy,sz,vx,vy,vz,fDist2,vdotr;
    int iDum;

    /*mdlDiag(smx->pkd->mdl, "smBallScatter:Start\n" );*/
    x = ri[0];
    y = ri[1];
    z = ri[2];
    vx = vi[0];
    vy = vi[1];
    vz = vi[2];
    if (smx->bPeriodic) {
        lx = smx->pkd->fPeriod[0];
        ly = smx->pkd->fPeriod[1];
        lz = smx->pkd->fPeriod[2];
        }
    else {
        lx = FLOAT_MAXVAL;
        ly = FLOAT_MAXVAL;
        lz = FLOAT_MAXVAL;
        }
    dt2 = (*pdt2);

    cp = ROOT;
    while (1) {
        pkdn = &c[cp];
        DTINTERSECT(pkdn,dt2,fDist2,lx,ly,lz,x,y,z,sx,sy,sz,iDum,vx,vy,vz,GetNextCell);
        /*
        ** We have an intersection to test.
        */
        if (fDist2 > INVDTOPENANGLE2*pkdn->bndDt.drMax2) {
            /* Particle Cell - passed Opening Angle Test */
            DTESTIMATOR(&c[cp],dtEst2,sx,sy,sz,vx,vy,vz);
            dtEst2 = dtEst2*dtEst2;
            if (dtEst2 < dt2) dt2 = dtEst2;
            goto GetNextCell;
            }
        if (c[cp].iDim >= 0) {
            cp = LOWER(cp);
//	    if ((x+.5)*(x+.5)+y*y+z*z < 0.033*0.033) printf("DT OPENCHILD %d\n",cp);
            continue;
            }
        else {
            pUpper = c[cp].pUpper;
            for (pj=c[cp].pLower;pj<=pUpper;++pj) {
                dx = sx - p[pj].r[0];
                dy = sy - p[pj].r[1];
                dz = sz - p[pj].r[2];
                fDist2 = (dx*dx + dy*dy + dz*dz);
                if (fDist2 < fBall2) continue;
                dvx = vx - p[pj].vPred[0];
                dvy = vy - p[pj].vPred[1];
                dvz = vz - p[pj].vPred[2];
                if ((vdotr = dvx*dx + dvy*dy + dvz*dz) < 0) {
                    double fDist = sqrt(fDist2);
                    dtEst2 = fDist/(p[pj].c-vdotr/fDist);
                    dtEst2 = dtEst2*dtEst2;
                    if (dtEst2 < dt2) dt2 = dtEst2;
                    }
                else {
                    dtEst2 = fDist2/(p[pj].c*p[pj].c);
                    if (dtEst2 < dt2) dt2 = dtEst2;
                    }
//		if ((x+.5)*(x+.5)+y*y+z*z < 0.033*0.033) printf("DT PART: %d %f  %f %f %f   %f %f %f   %f %f %f\n",pj,sqrt(dt2),sqrt(fDist2),p[pj].c,vdotr/sqrt(fDist2),p[pj].r[0],p[pj].r[1],p[pj].r[2],x,y,z);
                }
            }
    GetNextCell:
//	if ((x+.5)*(x+.5)+y*y+z*z < 0.033*0.033 && !bGood) printf("FAIL\n");
        SETNEXT(cp);
        if (cp == ROOT) break;
        }
    *pdt2 = dt2;
    /*mdlDiag(smx->pkd->mdl, "smBallScatter:End\n" );*/
#endif
    }

int smResmoothParticle(SMX smx,SMF *smf,int pi,FLOAT fBall2,int cp,FLOAT lx,FLOAT ly,FLOAT lz) {
    PKD pkd = smx->pkd;
    MDL mdl = smx->pkd->mdl;
    PARTICLE *p = smx->pkd->pStore;
    PARTICLE *pPart;
    KDN *pkdn;
    int pj,nCnt,id,i;
    FLOAT x,y,z,sx,sy,sz,dx,dy,dz,fDist2;
    int iDum;

    if (cp) {
        nCnt = smBallGatherNP(smx,fBall2,p[pi].r,cp);
        smx->fcnSmooth(&p[pi],nCnt,smx->nnList,smf);
        /*
        ** Try a cache check to improve responsiveness.
        */
        mdlCacheCheck(mdl);
        }
    else {
        if (smx->bPeriodic) {
            nCnt = smBallGather(smx,fBall2,p[pi].r);
            }
        else {
            nCnt = smBallGatherNP(smx,fBall2,p[pi].r,ROOT);
            }
        /*
        ** Start non-local search.
        */
        x = p[pi].r[0];
        y = p[pi].r[1];
        z = p[pi].r[2];
        cp = ROOT;
        id = -1;    /* We are in the LTT now ! */
        while (1) {
            if (id == pkd->idSelf) goto SkipLocal;
            if (id >= 0) pkdn = mdlAquire(mdl,CID_CELL,cp,id);
            else pkdn = &pkd->kdTop[cp];
            if (pkdn->pUpper < 0) goto GetNextCell;
            INTERSECT(pkdn,fBall2,lx,ly,lz,x,y,z,sx,sy,sz,iDum,GetNextCell);
            if (pkdn->iDim >= 0 || id == -1) {
                if (id >= 0) mdlRelease(mdl,CID_CELL,pkdn);
                pkdLower(pkd,cp,id);
                continue;
                }
            for (pj=pkdn->pLower;pj<=pkdn->pUpper;++pj) {
                pPart = mdlAquire(mdl,CID_PARTICLE,pj,id);
                dx = sx - pPart->r[0];
                dy = sy - pPart->r[1];
                dz = sz - pPart->r[2];
                fDist2 = dx*dx + dy*dy + dz*dz;
                if (fDist2 <= fBall2) {
                    if(nCnt >= smx->nListSize) smGrowList(smx);
                    smx->nnList[nCnt].fDist2 = fDist2;
                    smx->nnList[nCnt].dx = dx;
                    smx->nnList[nCnt].dy = dy;
                    smx->nnList[nCnt].dz = dz;
                    smx->nnList[nCnt].pPart = pPart;
                    smx->nnList[nCnt].iIndex = pj;
                    smx->nnList[nCnt].iPid = id;
                    smx->pbRelease[nCnt++] = 1;
                    continue;
                    }
                mdlRelease(mdl,CID_PARTICLE,pPart);
                }
        GetNextCell:
            if (id >= 0) mdlRelease(mdl,CID_CELL,pkdn);
        SkipLocal:
            pkdNext(pkd,cp,id);
            if (pkdIsRoot(cp,id)) break;
            }
        smx->fcnSmooth(&p[pi],nCnt,smx->nnList,smf);
        /*
        ** Release non-local particle pointers.
        */
        for (i=0;i<nCnt;++i) {
            if (smx->pbRelease[i]) {
                mdlRelease(mdl,CID_PARTICLE,smx->nnList[i].pPart);
                }
            }
        }

    return nCnt;
    }

PQ *smLoadPQ(SMX smx, int pi, int nSmooth, int nTree, FLOAT x, FLOAT y, FLOAT z,  FLOAT lx, FLOAT ly, FLOAT lz) {
    PKD pkd = smx->pkd;
    MDL mdl = smx->pkd->mdl;
    KDN *c = pkd->kdNodes;
    PARTICLE *p = pkd->pStore;
    KDN *pkdn;
    PARTICLE *pPart;
    int i,pj;
    int cell,id;
    FLOAT dx,dy,dz,sx,sy,sz;
    PQ *pq,*pqi;
    int iDum;
    int nQueue,iLoad;
    KDN dpkdn; /* dummy node */

    cell = ROOT;
    while (c[cell].iDim >= 0) {
        if (p[pi].r[c[cell].iDim] < c[cell].fSplit) cell = LOWER(cell);
        else cell = UPPER(cell);
        }
    /*
    ** Add local stuff to the prioq.
    */
    iLoad = c[cell].pLower;
    nQueue = 0;
    if (iLoad > nTree - nSmooth) iLoad = nTree - nSmooth;
    /*
    ** Have to add non-local particles to the queue?
    */
    for (id=0;id<mdlThreads(pkd->mdl) && (iLoad < 0);++id) {
        if (id == pkd->idSelf) continue;
        pkdn = mdlAquire(mdl,CID_CELL,ROOT,id);
        for (pj=pkdn->pLower;pj<=pkdn->pUpper && (iLoad < 0);++pj) {
            pqi = &smx->pq[nQueue];
            pPart = mdlAquire(mdl,CID_PARTICLE,pj,id);
            dx = x - pPart->r[0];
            dy = y - pPart->r[1];
            dz = z - pPart->r[2];
            pqi->ax = 0.0;
            pqi->ay = 0.0;
            pqi->az = 0.0;
            /*
            ** If nSmooth ~ N, may need to load the queue with ghosts...
            */
            for (i=0;i<3;i++) dpkdn.bnd.fMin[i] = dpkdn.bnd.fMax[i] = pPart->r[i];

            INTERSECT(&dpkdn,FLOAT_MAXVAL,lx,ly,lz,x,y,z,sx,sy,sz,iDum,
                dumlabel1);
            dx = sx - pPart->r[0];
            dy = sy - pPart->r[1];
            dz = sz - pPart->r[2];
            pqi->ax = sx - x;
            pqi->ay = sy - y;
            pqi->az = sz - z;
        dumlabel1:
            pqi->fKey = dx*dx + dy*dy + dz*dz;
            pqi->dx = dx;
            pqi->dy = dy;
            pqi->dz = dz;
            pqi->p = pj;
            pqi->id = id;
            pqi->pPart = pPart;
            PQ_HASHADD(smx->pqHash,smx->nHash,pqi);
            /*
            ** Note: in this case we DO NOT want to release pPart!
            */
            ++iLoad;
            ++nQueue;
            }
        mdlRelease(mdl,CID_CELL,pkdn);
        }
    assert(iLoad >= 0);
    for (;nQueue<nSmooth;++nQueue,++iLoad) {
        pqi = &smx->pq[nQueue];
        smSetMARK(iLoad);
        dx = x - p[iLoad].r[0];
        dy = y - p[iLoad].r[1];
        dz = z - p[iLoad].r[2];
        pqi->ax = 0.0;
        pqi->ay = 0.0;
        pqi->az = 0.0;
        for (i=0;i<3;i++)
            dpkdn.bnd.fMin[i] = dpkdn.bnd.fMax[i] = p[iLoad].r[i];
        INTERSECT(&dpkdn,FLOAT_MAXVAL,lx,ly,lz,x,y,z,sx,sy,sz,iDum,dumlabel2);
        dx = sx - p[iLoad].r[0];
        dy = sy - p[iLoad].r[1];
        dz = sz - p[iLoad].r[2];
        pqi->ax = sx - x;
        pqi->ay = sy - y;
        pqi->az = sz - z;
    dumlabel2:
        pqi->fKey = dx*dx + dy*dy + dz*dz;
        pqi->dx = dx;
        pqi->dy = dy;
        pqi->dz = dz;
        pqi->p = iLoad;
        pqi->id = pkd->idSelf;
        pqi->pPart = &p[iLoad];
        }
    PQ_BUILD(smx->pq,nSmooth,pq);
    return (pq);
    }

PQ *smRecentrePQ(SMX smx, PQ *pqNext, int nSmooth, FLOAT x, FLOAT y, FLOAT z) {
    int i;
    PQ *pq,*pqi;
    FLOAT dx,dy,dz;

    for (i=0,pqi=smx->pq;i<nSmooth;++i,++pqi) {
        if (pqi == pqNext) continue;
        pqi->ax -= pqNext->ax;
        pqi->ay -= pqNext->ay;
        pqi->az -= pqNext->az;
        dx = x - pqi->pPart->r[0] + pqi->ax;
        dy = y - pqi->pPart->r[1] + pqi->ay;
        dz = z - pqi->pPart->r[2] + pqi->az;
        pqi->fKey = dx*dx + dy*dy + dz*dz;
        pqi->dx = dx;
        pqi->dy = dy;
        pqi->dz = dz;
        }
    pqNext->fKey = 0.0;
    pqNext->ax = 0.0;
    pqNext->ay = 0.0;
    pqNext->az = 0.0;
    pqNext->dx = 0.0;
    pqNext->dy = 0.0;
    pqNext->dz = 0.0;
    PQ_BUILD(smx->pq,nSmooth,pq);

    return (pq);
    }

void smLoadNNFromPQ(NN *nnList,PQ *pqi) {
    nnList->iPid = pqi->id;
    nnList->iIndex = pqi->p;
    nnList->pPart = pqi->pPart;
    nnList->fDist2 = pqi->fKey;
    nnList->dx = pqi->dx;
    nnList->dy = pqi->dy;
    nnList->dz = pqi->dz;
    }


PQ *smBallSearchAll(SMX smx, PQ *pq, int pi, FLOAT x, FLOAT y, FLOAT z,  FLOAT lx, FLOAT ly, FLOAT lz) {
    PKD pkd = smx->pkd;
    MDL mdl = smx->pkd->mdl;
    PARTICLE *p = pkd->pStore;
    KDN *pkdn;
    PARTICLE *pPart;
    int pj;
    int cell,idcell,cp,id,ct,idct;
    FLOAT fBall2,fDist2,dx,dy,dz,sx,sy,sz;
    int iDum;

    if (smx->bPeriodic) {
        pq = smBallSearch(smx,pq,p[pi].r,&p[pi].cpStart);
        }
    else {
        pq = smBallSearchNP(smx,pq,p[pi].r,&p[pi].cpStart);
        }
    /*
    ** Start non-local search.
    */
    fBall2 = pq->fKey;
    idcell = -1;        /* We are in the LTT now ! */
    cell = pkd->piLeaf[pkd->idSelf];
    while (!pkdIsRoot(cell,idcell)) {
        cp = SIBLING(cell);
        id = idcell;
        ct = cp;
        idct = id;
        pkdNext(pkd,ct,idct);
        /*
        ** Check for any cache work.
        */
        while (1) {
            if (id >= 0) pkdn = mdlAquire(mdl,CID_CELL,cp,id);
            else pkdn = &pkd->kdTop[cp];
            if (pkdn->pUpper < 0) goto GetNext_2;
            INTERSECT(pkdn,fBall2,lx,ly,lz,x,y,z,sx,sy,sz,iDum,GetNext_2);
            if (pkdn->iDim >= 0 || id == -1) {
                if (id >= 0) mdlRelease(mdl,CID_CELL,pkdn);
                pkdLower(pkd,cp,id);
                continue;
                }
            for (pj=pkdn->pLower;pj<=pkdn->pUpper;++pj) {
                pPart = mdlAquire(mdl,CID_PARTICLE,pj,id);
                dx = sx - pPart->r[0];
                dy = sy - pPart->r[1];
                dz = sz - pPart->r[2];
                fDist2 = dx*dx + dy*dy + dz*dz;
                if (fDist2 < fBall2) {
                    p[pi].cpStart = 0;
                    PQ_INQUEUE(smx->pqHash,smx->nHash,pj,id,NextParticle); /* Check if already in queue, jump to next particle if true */
                    if (pq->id == pkd->idSelf) smResetMARK(pq->p); /* unmark furthest particle as in queue (local) */
                    else {
                        mdlRelease(mdl,CID_PARTICLE,pq->pPart); /*  release furthest particle (non local) */
                        PQ_HASHDEL(smx->pqHash,smx->nHash,pq); /* remove furthest particle from pqHash that marks non-local nbrs */
                        }
                    pq->fKey = fDist2;
                    pq->dx = dx;
                    pq->dy = dy;
                    pq->dz = dz;
                    pq->p = pj;
                    pq->id = id;
                    pq->pPart = pPart;
                    pq->ax = sx - x;
                    pq->ay = sy - y;
                    pq->az = sz - z;
                    PQ_HASHADD(smx->pqHash,smx->nHash,pq); /* add to pqHash that marks non-local nbrs */ 
                    PQ_REPLACE(pq);
                    fBall2 = pq->fKey;
                    /*
                    ** Note: in this case we DO NOT want to release pPart!
                    */
                    /*
                    ** Setting p[pi].cpStart to zero means we have to search 
                    ** non-local tree when doing a BallGather (smReSmooth).
                    */
                    continue;
                    }
            NextParticle:
                mdlRelease(mdl,CID_PARTICLE,pPart); // Release non-neighbours only
                }
        GetNext_2:
            if (id >= 0) mdlRelease(mdl,CID_CELL,pkdn); // Can always release cells
            pkdNext(pkd,cp,id);
            if (cp == ct && id == idct) break;
            }
        pkdParent(pkd,cell,idcell);
        }
    return(pq);
    }
    
void smSmooth(SMX smx,SMF *smf)
    {
    PKD pkd = smx->pkd;
    MDL mdl = smx->pkd->mdl;
    KDN *c = pkd->kdNodes;
    PARTICLE *p = pkd->pStore;
    int nSmooth,i,j,pi,pNext,nCnt;
    int cp;
    FLOAT fBall2,x,y,z,lx,ly,lz,r2Next;
    PQ *pq,*pqi,*pqNext;
    int iDum;
    int nTree,nLocal;
    int nSmoothed = 0,nSmoothedInner = 0,nSmoothedFixh = 0;
    int nSmoothInnerFail=0;
#ifdef NSMOOTHINNER
    int  nInner,nSmoothInnerCut = smx->nSmoothInner*0.8;
#endif
    nSmooth = smx->nSmooth;
    nTree = c[pkd->iRoot].pUpper + 1;
    nLocal = pkd->nLocal;
    lx = smx->lx; ly = smx->ly; lz = smx->lz;
    /*
    ** Clear Mark array and pqHash.
    */
    for (pi=0;pi<pkd->nLocal;++pi) {
        smResetMARK(pi);
#ifdef NSMOOTHINNER
        TYPEReset(&p[pi],TYPE_RESMOOTHINNER);
#endif
        }
    for (j=0;j<smx->nHash;++j) smx->pqHash[j] = NULL;

    for (pNext=0;pNext<nLocal;++pNext) { 
        /* Check every particle -- some done out of order by the snake */
        if (!TYPEFilter(&(pkd->pStore[pNext]),TYPE_SMOOTHACTIVE|TYPE_SMOOTHDONE,TYPE_SMOOTHACTIVE)) continue;

        pi = pNext;
        x = p[pi].r[0];    y = p[pi].r[1];    z = p[pi].r[2];
        /* Build PQ of nSmooth particles -- not a ball yet */
        pq = smLoadPQ(smx,pi,nSmooth,nTree,x,y,z,lx,ly,lz);
        
        for (;;) { /* Snake Loop */
            /*
            ** Priority Queue must be built. 'pi' must be defined.
            ** Squeeze PQ to actual ball 
            */
            pq = smBallSearchAll(smx,pq,pi,x,y,z,lx,ly,lz);
            /*
            ** Create Nearest-Neighbor List and try to pick next particle.
            */
            pqNext = NULL;
            r2Next = 2.0*pq->fKey;   /* arbitrarily bigger than pq->fKey! */
            fBall2 = pq->fKey;
            nCnt = 0;
            /* 
            ** If desired reject fBall2 below a minimum 
            ** We may get many more than nSmooth neighbours here -- resort to a ReSmooth
            */
            if (smx->iLowhFix && 
                    ((smx->iLowhFix==LOWHFIX_HOVERSOFT && fBall2 < smx->dfBall2OverSoft2*p[pi].fSoft*p[pi].fSoft) ||
                    (smx->iLowhFix==LOWHFIX_SINKRADIUS && fBall2 < smf->dSinkRadius*smf->dSinkRadius) ||
                    (smx->iLowhFix==LOWHFIX_SINKRADIUS_BUFF && fBall2 < smf->dSinkRadius*smf->dSinkRadius*1.1)) ) {
                /* We Resmooth for this guy later */
                p[pi].fBall2 = -1.0; /* any value < 0 will do -- see code after "DoneSmooth:" below */
                TYPESet(&p[pi],TYPE_SMOOTHDONE);

                /* Get the next in line */
                for (i=0,pqi=smx->pq;i<nSmooth;++i,++pqi) {
                    if (pqi->id != pkd->idSelf) continue;
                    if (!TYPEFilter(&(pkd->pStore[pqi->p]),TYPE_SMOOTHACTIVE|TYPE_SMOOTHDONE,TYPE_SMOOTHACTIVE)) continue;
                    if (pqi->fKey < r2Next) { /* prefer candidate closest to old centre */
                        pqNext = pqi;
                        r2Next = pqNext->fKey;
                        }
                    }
                }
            else {
                /* Limit fBall2 growth to help stability and neighbour finding */
                /*if (smx->bUseBallMax && p[pi].fBallMax > 0.0 && fBall2 > p[pi].fBallMax*p[pi].fBallMax)*/
                    /*fBall2=p[pi].fBallMax*p[pi].fBallMax;*/

                p[pi].fBall2 = fBall2;
                TYPESet(&p[pi],TYPE_SMOOTHDONE);

#ifdef NSMOOTHINNER
                nInner = 0;
#endif
                for (i=0,pqi=smx->pq;i<nSmooth;++i,++pqi) {
                    /*
                    ** There are cases like in collisions where we do want to include the
                    ** single nearest neighbor. So include the most distant particle.
                    **
                    if (pqi == pq) continue;
                    */
                    /*
                    ** Move relevant data into Nearest Neighbor array.
                    */
                    if (pqi->fKey <= fBall2) {
#ifdef NSMOOTHINNER
                        if (pqi->fKey <= fBall2*smx->fBall2InnerFrac) nInner++;
#endif
                        smLoadNNFromPQ(&smx->nnList[nCnt],pqi);
                        ++nCnt;
                        } 
                    if (pqi->id != pkd->idSelf) continue;
                    if (!TYPEFilter(&(pkd->pStore[pqi->p]),TYPE_SMOOTHACTIVE|TYPE_SMOOTHDONE,TYPE_SMOOTHACTIVE)) continue;
                    if (pqi->fKey < r2Next) {
                        pqNext = pqi;
                        r2Next = pqNext->fKey;
                        }
                    }

#ifdef NSMOOTHINNER
                if (nInner < nSmoothInnerCut) { /* Defer smooth */
                    ISORT *isort;
                    isort = (ISORT *) malloc(sizeof(ISORT)*nCnt);
                    for (i=0;i<nCnt;++i) {
                        isort[i].r2 = smx->nnList[i].fDist2;
                        }
                    qsort( isort, nCnt, sizeof(ISORT), CompISORT );
                    if(nCnt > smx->nSmoothInner) { 
                        NSI_DEBUG(fprintf(stderr,"%d: nSmooth Inner Fail %d : nSm %d %d %d  fBall: %f %f\n",pkd->idSelf,p[pi].iOrder,nInner,nSmoothInnerCut,smx->nSmoothInner,sqrt(p[pi].fBall2),sqrt(isort[smx->nSmoothInner-1].r2/smx->fBall2InnerFrac)););
                        p[pi].fBall2 = -isort[smx->nSmoothInner-1].r2/smx->fBall2InnerFrac; 
                    }
                    else {
                        NSI_DEBUG(fprintf(stderr,"%d: nSmooth Inner Fail %d : nSm %d %d %d  fBall: %f %f\n",pkd->idSelf,p[pi].iOrder,nInner,nSmoothInnerCut,nCnt,sqrt(p[pi].fBall2),sqrt(isort[nCnt-1].r2/smx->fBall2InnerFrac)););
                        p[pi].fBall2 = -isort[nCnt-1].r2/smx->fBall2InnerFrac; 
                    }
                    TYPESet(&p[pi],TYPE_RESMOOTHINNER);
                    nSmoothInnerFail++;
                    free(isort);
                    }
                else  
#endif 
                    {
                    smx->fcnSmooth(&p[pi],nCnt,smx->nnList,smf);
                    nSmoothed++;
                    smLargefBallCheck(smx,&p[pi],lx,ly,lz);
                    }
                }

            /*
            ** Need to do a CACHE recombine (ie. finish up Smooth)
            ** to deliver info to particles on other processors.
            */
            /*
            ** Try a cache check to improve responsiveness.
            */
            mdlCacheCheck(mdl);

            if (!pqNext) break; /* end Snake */
            /*
            ** Recalculate the priority queue using the previous particles.
            ** "THE SNAKE"
            */
            pi = pqNext->p;
            x = p[pi].r[0];    y = p[pi].r[1];    z = p[pi].r[2];
            pq = smRecentrePQ(smx,pqNext,nSmooth,x,y,z);
            }
        /*
        ** No next candidate in current nbr list 
        ** Clean up Marks, pqHash, and release all acquired pointers!
        ** Start PQ from scratch again
        */
        for (i=0,pqi=smx->pq;i<nSmooth;++i,++pqi) {
            if (pqi->id == pkd->idSelf) smResetMARK(pqi->p);
            else {
                mdlRelease(mdl,CID_PARTICLE,pqi->pPart);
                PQ_HASHDEL(smx->pqHash,smx->nHash,pqi);
                }
            }
        }
    /* 
    ** Regular smooth complete -- standard PQ no-longer needed
    ** fix failed smooth candidates 
    */
#ifdef NSMOOTHINNER
    if (nSmoothInnerFail) {
        smFinishList(smx);
        nSmooth = smx->nSmooth = smx->nSmoothMax;
        smInitList(smx,nSmooth);
        /* Clear Mark array and pqHash */
        for (pi=0;pi<pkd->nLocal;++pi) {
            smResetMARK(pi);
            if (TYPETest(&p[pi],TYPE_RESMOOTHINNER)) TYPEReset( &p[pi],TYPE_SMOOTHDONE );
            }
        for (j=0;j<smx->nHash;++j) smx->pqHash[j] = NULL;
        }
#endif

    for (pNext=0;pNext<nLocal;++pNext) {
        if (!TYPETest(&(p[pNext]),TYPE_SMOOTHACTIVE) || p[pNext].fBall2 >= 0.0) continue;
        pi = pNext;
#ifdef NSMOOTHINNER
        if (TYPETest(&p[pi],TYPE_RESMOOTHINNER)) {
            /* Cap max nbrs:   original nSmooth < nCnt <= nSmoothMax */
            x = p[pi].r[0];        y = p[pi].r[1];         z = p[pi].r[2];
            pq = smLoadPQ(smx,pi,nSmooth,nTree,x,y,z,lx,ly,lz);
            for (;;) { /* SNAKE loop */
                pq = smBallSearchAll(smx,pq,pi,x,y,z,lx,ly,lz);
                pqNext = NULL;
                r2Next = 2.0*pq->fKey;  
                fBall2 = pq->fKey;
                p[pi].fBall2 = fabs(p[pi].fBall2);
                if (fBall2 > p[pi].fBall2) fBall2 = p[pi].fBall2;
                else p[pi].fBall2 = fBall2;
                TYPESet(&p[pi],TYPE_SMOOTHDONE);
                
                nCnt = 0;
                NSI_DEBUG(nInner = 0;)
                for (i=0,pqi=smx->pq;i<nSmooth;++i,++pqi) {
                    NSI_DEBUG(if (pqi->fKey <= fBall2*smx->fBall2InnerFrac) nInner++;);
                    if (pqi->fKey <= fBall2) {
                        smLoadNNFromPQ(&smx->nnList[nCnt],pqi);
                        ++nCnt;
                        } 
                    if (pqi->id != pkd->idSelf) continue;
                    if (!TYPEFilter(&(pkd->pStore[pqi->p]),TYPE_RESMOOTHINNER|TYPE_SMOOTHDONE,TYPE_RESMOOTHINNER)) continue;
                    if (pqi->fKey < r2Next) {
                        pqNext = pqi;
                        r2Next = pqNext->fKey;
                        }
                    }
                NSI_DEBUG(fprintf(stderr,"%d nSmooth Inner Redo %d : nSm %d %d %d  fBall: %f\n",pkd->idSelf,p[pi].iOrder,nInner,smx->nSmoothInner,nCnt,sqrt(p[pi].fBall2)););


                smx->fcnSmooth(&p[pi],nCnt,smx->nnList,smf);
                nSmoothed++;
                nSmoothedInner++;
                smLargefBallCheck(smx,&p[pi],lx,ly,lz);
                mdlCacheCheck(mdl);
                if (!pqNext) break;
                pi = pqNext->p;
                x = p[pi].r[0];    y = p[pi].r[1];    z = p[pi].r[2];
                pq = smRecentrePQ(smx,pqNext,nSmooth,x,y,z);
                }
            /* Snake end: Clean up marking/release cache */
            for (i=0,pqi=smx->pq;i<nSmooth;++i,++pqi) {
                if (pqi->id == pkd->idSelf) smResetMARK(pqi->p);
                else {
                    mdlRelease(mdl,CID_PARTICLE,pqi->pPart);
                    PQ_HASHDEL(smx->pqHash,smx->nHash,pqi);
                    }
                }
            }
        else
#endif 
            {
            switch (smx->iLowhFix) {
            case LOWHFIX_HOVERSOFT:
                fBall2 = smx->dfBall2OverSoft2*p[pi].fSoft*p[pi].fSoft;
                break;
            case LOWHFIX_SINKRADIUS:
                fBall2 = smf->dSinkRadius*smf->dSinkRadius;
                break;
            case LOWHFIX_SINKRADIUS_BUFF:
                fBall2 = smf->dSinkRadius*smf->dSinkRadius*1.1;
                break;
            default:
                fprintf(stderr,"Illegal value for iLowhFix %d in smooth\n",smx->iLowhFix);
                assert(0);
                }
            /* Do a Ball Gather at r^2 = fBall2 */
            p[pi].fBall2 = fBall2;
            p[pi].cpStart = cp = 0;
            smResmoothParticle(smx,smf,pi,fBall2,cp,lx,ly,lz);
            nSmoothed++;
            nSmoothedFixh++;
            }
        }
    assert(nSmoothedInner == nSmoothInnerFail);
    smf->nSmoothed = nSmoothed;
    smf->nSmoothedInner = nSmoothedInner;
    smf->nSmoothedFixh = nSmoothedFixh;
    }


void smReSmooth(SMX smx,SMF *smf)
    {
    PKD pkd = smx->pkd;
    PARTICLE *p = smx->pkd->pStore;
    int pi,cp;
    FLOAT lx,ly,lz,fBall2;
    int nTree;
        
    if (smx->bPeriodic) {
        lx = pkd->fPeriod[0];
        ly = pkd->fPeriod[1];
        lz = pkd->fPeriod[2];
        }
    else {
        lx = FLOAT_MAXVAL;
        ly = FLOAT_MAXVAL;
        lz = FLOAT_MAXVAL;
        }

    nTree = pkd->kdNodes[pkd->iRoot].pUpper + 1;
    for (pi=0;pi<nTree;++pi) {
        if (!TYPEFilter(&(p[pi]),TYPE_SMOOTHACTIVE|TYPE_SMOOTHDONE, TYPE_SMOOTHACTIVE)) continue;
        /*
        ** Do a Ball Gather at the radius of the most distant particle
        ** which smSmooth sets in p[pi].fBall2.
        */
        fBall2 = p[pi].fBall2;
        cp = p[pi].cpStart;

        /* Resmooth for one particle.  Handles all cases: local, periodic and non-local */
        smResmoothParticle(smx,smf,pi,fBall2,cp,lx,ly,lz);
        TYPESet(&p[pi],TYPE_SMOOTHDONE);
        }
    }

void smMarkSmooth(SMX smx,SMF *smf,int iMarkType)
    {
    PKD pkd = smx->pkd;
    MDL mdl = smx->pkd->mdl;
    PARTICLE *p = smx->pkd->pStore;
    PARTICLE *pPart;
    KDN *pkdn;
    int pi,pj,cp,id;
    FLOAT x,y,z,lx,ly,lz,sx,sy,sz,dx,dy,dz,fDist2;
    int iDum;
    int nTree;
    char DiagStr[100];
        
    mdlDiag(smx->pkd->mdl,"smMarkSmooth:Start\n" );
    lx = smx->lx; ly = smx->ly; lz = smx->lz;
    nTree = pkd->kdNodes[pkd->iRoot].pUpper + 1;
    for (pi=0;pi<nTree;++pi) {
        if (!TYPEFilter(&(p[pi]),TYPE_SMOOTHACTIVE|TYPE_SMOOTHDONE,
                TYPE_SMOOTHACTIVE)) continue;
        /*              sprintf( DiagStr, "smMarkSmooth: Particle %d\n", p[pi].iOrder);*/
        mdlDiag(smx->pkd->mdl, DiagStr );
        /*
        ** Do a Ball Scatter to this particle
        */
        TYPESet(&p[pi],TYPE_SMOOTHDONE);

        if (smx->bPeriodic) {
            smBallScatter(smx,p[pi].r,iMarkType);
            }
        else {
            smBallScatterNP(smx,p[pi].r,iMarkType,ROOT);
            }
        /*
        ** Start non-local search.
        */
        x = p[pi].r[0];
        y = p[pi].r[1];
        z = p[pi].r[2];
        cp = ROOT;
        id = -1;        /* We are in the LTT now ! */
        while (1) {
            if (id == pkd->idSelf) goto SkipLocal;
            if (id >= 0) pkdn = mdlAquire(mdl,CID_CELL,cp,id);
            else pkdn = &pkd->kdTop[cp];
            if (pkdn->pUpper < 0) goto GetNextCell;
            INTERSECTSCATTER(pkdn,lx,ly,lz,x,y,z,sx,sy,sz,iDum,GetNextCell);
            if (pkdn->iDim >= 0 || id == -1) {
                if (id >= 0) mdlRelease(mdl,CID_CELL,pkdn);
                pkdLower(pkd,cp,id);
                continue;
                }
            for (pj=pkdn->pLower;pj<=pkdn->pUpper;++pj) {
                pPart = mdlAquire(mdl,CID_PARTICLE,pj,id);
                dx = sx - pPart->r[0];
                dy = sy - pPart->r[1];
                dz = sz - pPart->r[2];
                fDist2 = dx*dx + dy*dy + dz*dz;
                if (fDist2 <= pPart->fBallMax*pPart->fBallMax) 
                    TYPESet(pPart, iMarkType);  
                mdlRelease(mdl,CID_PARTICLE,pPart);
                }
        GetNextCell:
            if (id >= 0) mdlRelease(mdl,CID_CELL,pkdn);
        SkipLocal:
            pkdNext(pkd,cp,id);
            if (pkdIsRoot(cp,id)) break;
            }
        }
    }


void smDtSmooth(SMX smx,SMF *smf)
    {
#ifdef GASOLINE
    PKD pkd = smx->pkd;
    MDL mdl = smx->pkd->mdl;
    PARTICLE *p = smx->pkd->pStore;
    PARTICLE *pPart;
    KDN *pkdn;
    int pi,pj,cp,id;
    FLOAT dt2,dtEst2,dtNew,dtMin,dx,dy,dz,dvx,dvy,dvz;
    FLOAT x,y,z,lx,ly,lz,sx,sy,sz,vx,vy,vz,fDist2,vdotr,fBall2;
    int iDum;
    int nTree;
    char DiagStr[100];
        
    mdlDiag(smx->pkd->mdl,"smDtSmooth:Start\n" );
    if (smx->bPeriodic) {
        lx = pkd->fPeriod[0];
        ly = pkd->fPeriod[1];
        lz = pkd->fPeriod[2];
        }
    else {
        lx = FLOAT_MAXVAL;
        ly = FLOAT_MAXVAL;
        lz = FLOAT_MAXVAL;
        }
    dtMin = FLOAT_MAXVAL;

    nTree = pkd->kdNodes[pkd->iRoot].pUpper + 1;
    for (pi=0;pi<nTree;++pi) {
        if (!TYPEFilter(&(p[pi]),TYPE_SMOOTHACTIVE|TYPE_SMOOTHDONE,
                TYPE_SMOOTHACTIVE)) continue;
        /*              sprintf( DiagStr, "smDtSmooth: Particle %d\n", p[pi].iOrder);*/
        mdlDiag(smx->pkd->mdl, DiagStr );
        /*
        ** Do a Ball Gather to this particle
        */
        TYPESet(&p[pi],TYPE_SMOOTHDONE);
        dt2 = p[pi].dt;
        dt2 /= smf->dEtaCourantLong; /* factor of (1/eta) corrected later */
        dt2 = dt2*dt2;
        fBall2 = p[pi].fBall2;
        /* Hopefully sph basic dt small so most intersects fail */

        if (smx->bPeriodic) {
            smDtBall(smx,p[pi].r,p[pi].v,&dt2,fBall2);
            }
        else {
            smDtBallNP(smx,p[pi].r,p[pi].v,&dt2,fBall2,ROOT);
            }
        /*
        ** Start non-local search.
        */
        x = p[pi].r[0];
        y = p[pi].r[1];
        z = p[pi].r[2];
//	if ((x+.5)*(x+.5)+y*y+z*z < 0.033*0.033) printf("DT for: %d %f  %f %f %f  \n",pi,p[pi].c,x,y,z);
        vx = p[pi].vPred[0];
        vy = p[pi].vPred[1];
        vz = p[pi].vPred[2];
        cp = ROOT;
        id = -1;        /* We are in the LTT now ! */
        while (1) {
            if (id == pkd->idSelf) goto SkipLocal;
            if (id >= 0) pkdn = mdlAquire(mdl,CID_CELL,cp,id);
            else pkdn = &pkd->kdTop[cp];
            if (pkdn->pUpper < 0) goto GetNextCell;
            DTINTERSECT(pkdn,dt2,fDist2,lx,ly,lz,x,y,z,sx,sy,sz,iDum,vx,vy,vz,GetNextCell);
            if (fDist2 > INVDTOPENANGLE2*pkdn->bndDt.drMax2) {
                /* Particle Cell - passed Opening Angle Test */
                DTESTIMATOR(pkdn,dtEst2,sx,sy,sz,vx,vy,vz);
                dtEst2 = dtEst2*dtEst2;
                if (dtEst2 < dt2) dt2 = dtEst2;
                goto GetNextCell;
                }
            if (pkdn->iDim >= 0 || id == -1) {
                if (id >= 0) mdlRelease(mdl,CID_CELL,pkdn);
                pkdLower(pkd,cp,id);
                continue;
                }
            for (pj=pkdn->pLower;pj<=pkdn->pUpper;++pj) {
                pPart = mdlAquire(mdl,CID_PARTICLE,pj,id);
                dx = sx - pPart->r[0];
                dy = sy - pPart->r[1];
                dz = sz - pPart->r[2];
                fDist2 = (dx*dx + dy*dy + dz*dz);
                if (fDist2 < fBall2) continue;
                dvx = vx - pPart->vPred[0];
                dvy = vy - pPart->vPred[1];
                dvz = vz - pPart->vPred[2];
                if ((vdotr = dvx*dx + dvy*dy + dvz*dz) < 0) {
                    double fDist = sqrt(fDist2);
                    dtEst2 = fDist/(pPart->c-vdotr/fDist);
                    dtEst2 = dtEst2*dtEst2;
                    if (dtEst2 < dt2) dt2 = dtEst2;
                    }
                else {
                    dtEst2 = fDist2/(pPart->c*pPart->c);
                    if (dtEst2 < dt2) dt2 = dtEst2;
                    }

                mdlRelease(mdl,CID_PARTICLE,pPart);
                }
        GetNextCell:
            if (id >= 0) mdlRelease(mdl,CID_CELL,pkdn);
        SkipLocal:
            pkdNext(pkd,cp,id);
            if (pkdIsRoot(cp,id)) break;
            }
        dtNew = smf->dEtaCourantLong*sqrt(dt2); /* Multiply by eta */
        if (dtNew < dtMin) dtMin = dtNew;
        p[pi].dt = dtNew;
        }
    smf->dtMin = dtMin;
#endif
    }
