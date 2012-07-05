#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "ewald.h"
#include "pkd.h"
#include "meval.h"
#include "qeval.h"

#ifdef __crayx1
#include "erf.c"
#endif

int pkdBucketEwald(PKD pkd,int iBucket,int nReps,double fEwCut,int iOrder)
{
	KDN *pkdn;
	PARTICLE *p;
	struct ilCellNewt mom;
	int i,j,n,ix,iy,iz,nEwReps,bInHole,bInHolex,bInHolexy;
	double L,alpha,alpha2,alphan,k1,ka;
	double fEwCut2;
	double fPot,ax,ay,az;
	double dx,dy,dz,dxo,dyo,dzo,r2,r,dir,dir2,a;
	double gam[6];
	double hdotx,s,c;
	int nFlop;
	int nActive = 0;
	int nLoop = 0;
#ifdef REDUCED_EWALD
	int nMultiFlop[5] = QEVAL_FLOP;
#else
	int nMultiFlop[5] = MEVAL_FLOP;
#endif
	
#ifdef __crayx1
	/* Optimization for vector processors */
	int tmp_n;
	int ixStride, iyStride;
#endif


	if (!iOrder) return 0;
	mom = pkd->ilcnRoot;
	pkdn = &pkd->kdNodes[iBucket];
	n = pkdn->pUpper - pkdn->pLower + 1;
	p = &pkd->pStore[pkdn->pLower];
	nEwReps = ceil(fEwCut);
	L = pkd->fPeriod[0];
	fEwCut2 = fEwCut*fEwCut*L*L;
	nEwReps = nEwReps > nReps ? nEwReps : nReps;
	alpha = 2.0/L;
	alpha2 = alpha*alpha;
	k1 = M_PI/(alpha2*L*L*L);
	ka = 2.0*alpha/sqrt(M_PI);
	for(j=0;j<n;++j) {
		if (!TYPEQueryACTIVE(&(p[j]))) continue;
		fPot = mom.m*k1;
		ax = 0.0;
		ay = 0.0;
		az = 0.0;
		dx = p[j].r[0] - mom.x;
		dy = p[j].r[1] - mom.y;
		dz = p[j].r[2] - mom.z;

#ifndef __crayx1
		for (ix=-nEwReps;ix<=nEwReps;++ix) {
			bInHolex = (ix >= -nReps && ix <= nReps);
			dxo = dx + ix*L;
			for(iy=-nEwReps;iy<=nEwReps;++iy) {
				bInHolexy = (bInHolex && iy >= -nReps && iy <= nReps);
				dyo = dy + iy*L;
				for(iz=-nEwReps;iz<=nEwReps;++iz) {
#else
		/* Optimized for vector machines */

		ixStride = (2*nEwReps + 1)*(2*nEwReps + 1);
		iyStride = 2*nEwReps + 1;
		for(tmp_n = 0; tmp_n < ixStride*iyStride; tmp_n++) {
			ix = tmp_n/ixStride - nEwReps;
			iy = (tmp_n%ixStride)/iyStride - nEwReps;
			iz = tmp_n%iyStride - nEwReps;
			bInHolex = (ix >= -nReps && ix <= nReps);
			dxo = dx + ix*L;
				bInHolexy = (bInHolex && iy >= -nReps && iy <= nReps);
				dyo = dy + iy*L;

#endif
					bInHole = (bInHolexy && iz >= -nReps && iz <= nReps);
					/*
					 ** Scoring for Ewald inner stuff = (+,*)
					 **		Visible ops 		= (8,28)
					 **		sqrt, 1/sqrt est. 	= (6,11)
					 **		exp est.			= (6,11)  same as sqrt.
					 **		erf/erfc est.		= (12,22) twice a sqrt.	
					 **		Subtotal			= (32,72) = 104
					 **		Total				= 104 + nMultiFlop[iOrder]
					 */
					dzo = dz + iz*L;
					r2 = dxo*dxo + dyo*dyo + dzo*dzo;
					if (r2 > fEwCut2 && !bInHole) continue;
					if (r2 < 3.0e-3*L*L) {
					/*
					 * For small r, series expand about
					 * the origin to avoid errors caused
					 * by cancellation of large terms.
					 */
					  alphan = ka;
					  r2 *= alpha2;
					  gam[0] = alphan*(r2/3 - 1);
					  alphan *= 2*alpha2;
					  gam[1] = alphan*(r2/5 - 1.0/3.0);
					  alphan *= 2*alpha2;
					  gam[2] = alphan*(r2/7 - 1.0/5.0);
					  alphan *= 2*alpha2;
					  gam[3] = alphan*(r2/9 - 1.0/7.0);
					  alphan *= 2*alpha2;
					  gam[4] = alphan*(r2/11 - 1.0/9.0);
					  alphan *= 2*alpha2;
					  gam[5] = alphan*(r2/13 - 1.0/11.0);
					  }
					else {
					    r = sqrt(r2);
					    dir = 1/r;
					    dir2 = dir*dir;
					    a = exp(-r2*alpha2);
					    a *= ka*dir2;
					    if (bInHole) gam[0] = -erf(alpha*r);
					    else gam[0] = erfc(alpha*r);
					    gam[0] *= dir;
					    gam[1] = gam[0]*dir2 + a;
					    alphan = 2*alpha2;
					    gam[2] = 3*gam[1]*dir2 + alphan*a;
					    alphan *= 2*alpha2;
					    gam[3] = 5*gam[2]*dir2 + alphan*a;
					    alphan *= 2*alpha2;
					    gam[4] = 7*gam[3]*dir2 + alphan*a;
					    alphan *= 2*alpha2;
					    gam[5] = 9*gam[4]*dir2 + alphan*a;
					    }
#ifdef REDUCED_EWALD
					QEVAL(iOrder,mom,gam,dxo,dyo,dzo,ax,ay,az,fPot);
#else
					MEVAL(iOrder,mom,gam,dxo,dyo,dzo,ax,ay,az,fPot);
#endif
					++nLoop;

#ifndef __crayx1
					}
				}
#endif
			}
		/*
		 ** Try a cache check to improve responsiveness.
		 */
		mdlCacheCheck(pkd->mdl);
		/*
		 ** Scoring for the h-loop (+,*)
		 ** 	Without trig = (10,14)
		 **	    Trig est.	 = 2*(6,11)  same as 1/sqrt scoring.
		 **		Total        = (22,36)
		 **					 = 58
		 */
		for (i=0;i<pkd->nEwhLoop;++i) {
			hdotx = pkd->ewt[i].hx*dx + pkd->ewt[i].hy*dy + pkd->ewt[i].hz*dz;
			c = cos(hdotx);
			s = sin(hdotx);
			fPot += pkd->ewt[i].hCfac*c + pkd->ewt[i].hSfac*s;
			ax += pkd->ewt[i].hx*(pkd->ewt[i].hCfac*s - pkd->ewt[i].hSfac*c);
			ay += pkd->ewt[i].hy*(pkd->ewt[i].hCfac*s - pkd->ewt[i].hSfac*c);
			az += pkd->ewt[i].hz*(pkd->ewt[i].hCfac*s - pkd->ewt[i].hSfac*c);
			}
		p[j].fPot += fPot;
		p[j].a[0] += ax;
		p[j].a[1] += ay;
		p[j].a[2] += az;
		++nActive;
		/*
		 ** Try a cache check to improve responsiveness.
		 */
		mdlCacheCheck(pkd->mdl);
	    }
	nFlop = nLoop*(104 + nMultiFlop[iOrder]) + 
		nActive*pkd->nEwhLoop*58;
	return(nFlop);
	}



void pkdEwaldInit(PKD pkd,double fhCut,int iOrder)
{
	struct ilCellNewt mom;
	int i,hReps,hx,hy,hz,h2;
	double alpha,k4,L;
	double gam[6],mfacc,mfacs;
	double ax,ay,az;

	mom = pkd->ilcnRoot;
	/*
	 ** Now setup stuff for the h-loop.
	 */
	hReps = ceil(fhCut);
	L = pkd->fPeriod[0];
	alpha = 2.0/L;
	k4 = M_PI*M_PI/(alpha*alpha*L*L);
	i = 0;
	for (hx=-hReps;hx<=hReps;++hx) {
		for (hy=-hReps;hy<=hReps;++hy) {
			for (hz=-hReps;hz<=hReps;++hz) {
				h2 = hx*hx + hy*hy + hz*hz;
				if (h2 == 0) continue;
				if (h2 > fhCut*fhCut) continue;
				if (i == pkd->nMaxEwhLoop) {
					pkd->nMaxEwhLoop *= 2;
					pkd->ewt = (EWT *) realloc(pkd->ewt,pkd->nMaxEwhLoop*sizeof(EWT));
					assert(pkd->ewt != NULL);
					}
				gam[0] = exp(-k4*h2)/(M_PI*h2*L);
				gam[1] = 2*M_PI/L*gam[0];
				gam[2] = -2*M_PI/L*gam[1];
				gam[3] = 2*M_PI/L*gam[2];
				gam[4] = -2*M_PI/L*gam[3];
				gam[5] = 2*M_PI/L*gam[4];
				gam[1] = 0.0;
				gam[3] = 0.0;
				gam[5] = 0.0;
				ax = 0.0;
				ay = 0.0;
				az = 0.0;
				mfacc = 0.0;
				QEVAL(iOrder,mom,gam,hx,hy,hz,ax,ay,az,mfacc);
				gam[0] = exp(-k4*h2)/(M_PI*h2*L);
				gam[1] = 2*M_PI/L*gam[0];
				gam[2] = -2*M_PI/L*gam[1];
				gam[3] = 2*M_PI/L*gam[2];
				gam[4] = -2*M_PI/L*gam[3];
				gam[5] = 2*M_PI/L*gam[4];
				gam[0] = 0.0;
				gam[2] = 0.0;
				gam[4] = 0.0;
				ax = 0.0;
				ay = 0.0;
				az = 0.0;
				mfacs = 0.0;
				QEVAL(iOrder,mom,gam,hx,hy,hz,ax,ay,az,mfacs);
				pkd->ewt[i].hx = 2*M_PI/L*hx;
				pkd->ewt[i].hy = 2*M_PI/L*hy;
				pkd->ewt[i].hz = 2*M_PI/L*hz;
				pkd->ewt[i].hCfac = mfacc;
				pkd->ewt[i].hSfac = mfacs;
				++i;
				}
			}
		}
	pkd->nEwhLoop = i;
	}
