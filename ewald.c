#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <assert.h>
#include "ewald.h"
#include "pkd.h"


void pkdBucketEwald(PKD pkd,int iBucket,int nReps,float fEwCut)
{
	KDN *pkdn,*pbox;
	PARTICLE *p;
	int i,j,n,ix,iy,iz,nEwReps,bInHolex,bInHolexy,bInHole;
	float L,alpha,alpha2,alpha4,k1,k6;
	float fEwCut2;
	float Mass,Qxx,Qyy,Qzz,Qxy,Qxz,Qyz,Qtr;
	float fPot,ax,ay,az;
	float dx,dy,dz,dxo,dyo,dzo,r2,r,dir,dir2,aq,bq,dq,eq,eralph;
	float qirx,qiry,qirz,qir,qir3;
	float hdotx,s;
	
	pkdn = &pkd->kdNodes[iBucket];
	pbox = &pkd->kdTop[ROOT];
	n = pkdn->pUpper - pkdn->pLower + 1;
	p = &pkd->pStore[pkdn->pLower];
	nEwReps = ceil(fEwCut);
	L = pkd->fPeriod[0];
	fEwCut2 = fEwCut*fEwCut*L*L;
	nEwReps = nEwReps > nReps ? nEwReps : nReps;
	alpha = 2.0/L;
	alpha2 = alpha*alpha;
	alpha4 = alpha2*alpha2;
	k1 = M_PI/(alpha2*L*L*L);
	k6 = 2.0*alpha/sqrt(M_PI);
	Mass = pbox->fMass;
	Qxx = pbox->mom.Qxx;
	Qyy = pbox->mom.Qyy;
	Qzz = pbox->mom.Qzz;
	Qxy = pbox->mom.Qxy;
	Qxz = pbox->mom.Qxz;
	Qyz = pbox->mom.Qyz;
	Qtr = 0.5*(Qxx + Qyy + Qzz);
	for(j=0;j<n;++j) {
		fPot = Mass*k1;
		ax = 0.0;
		ay = 0.0;
		az = 0.0;
		dx = pbox->r[0] - p[j].r[0];
		dy = pbox->r[1] - p[j].r[1];
		dz = pbox->r[2] - p[j].r[2];
		for (ix=-nEwReps;ix<=nEwReps;++ix) {
			bInHolex = (ix >= -nReps && ix <= nReps);
			dxo = dx + ix*L;
			for(iy=-nEwReps;iy<=nEwReps;++iy) {
				bInHolexy = (bInHolex && iy >= -nReps && iy <= nReps);
				dyo = dy + iy*L;
				for(iz=-nEwReps;iz<=nEwReps;++iz) {
					bInHole = (bInHolexy && iz >= -nReps && iz <= nReps);
					dzo = dz + iz*L;
					r2 = dxo*dxo + dyo*dyo + dzo*dzo;
					if (r2 == 0.0) continue;
					if (r2 > fEwCut2 && !bInHole) continue;
					r = sqrt(r2);
					dir = 1/r;
					dir2 = dir*dir;
					aq = exp(-r2*alpha2);
#if 0
					/*
					 ** Using this approximation seems to cause problems 
					 ** when large large numbers of particles are involved.
					 ** On the 1.35 million particle run we had to use 
					 ** the actual functions for erf and erfc, and used 
					 ** the GNU versions of these on the KSR due to a bug
					 ** in the KSR math library.
					 */
					t = 1/(1 + 0.3275911*alpha*r);
					eralph = aq*t*(0.254829592 - t*(0.284496736 - t*(1.421413741 - t*(1.453152027 - t*1.061405429))));
					if (bInHole) eralph -= 1.0;
#endif
					if (bInHole) eralph = -erf(alpha*r);
					else eralph = erfc(alpha*r);
					eralph *= dir;
					aq *= k6*dir2;
					bq = aq + eralph*dir2;
					dq = 2*alpha2*aq + 3*bq*dir2;
					eq = 4*alpha4*aq + 5*dq*dir2;
					qirx = Qxx*dxo + Qxy*dyo + Qxz*dzo;
					qiry = Qxy*dxo + Qyy*dyo + Qyz*dzo;
					qirz = Qxz*dxo + Qyz*dyo + Qzz*dzo;
					qir = 0.5*(qirx*dxo + qiry*dyo + qirz*dzo);
					qir3 = Mass*bq + eq*qir - dq*Qtr;
					fPot -= Mass*eralph + dq*qir - bq*Qtr;
					ax += qir3*dxo - dq*qirx;
					ay += qir3*dyo - dq*qiry;
					az += qir3*dzo - dq*qirz;
					}
				}
			}
		for (i=0;i<pkd->nEwhLoop;++i) {
			hdotx = pkd->ewt[i].hx*dx + pkd->ewt[i].hy*dy + pkd->ewt[i].hz*dz;
			s = sin(hdotx);
			fPot -= pkd->ewt[i].k2*cos(hdotx);
			ax += pkd->ewt[i].k3x*s;
			ay += pkd->ewt[i].k3y*s;
			az += pkd->ewt[i].k3z*s;
			}
		p[j].fPot += fPot;
		p[j].a[0] += ax;
		p[j].a[1] += ay;
		p[j].a[2] += az;
	    }
	}



void pkdEwaldInit(PKD pkd,float fhCut)
{
	KDN *pbox;
	int i,hReps,hx,hy,hz,h2;
	float alpha,k2,k3,k4,k5,k7,exph2,qfac,dir3,L;

	hReps = ceil(fhCut);
	L = pkd->fPeriod[0];
	alpha = 2.0/L;
	k2 = 1/(M_PI*L);
	k3 = 2/(L*L);
	k4 = M_PI*M_PI/(alpha*alpha*L*L);
	k5 = 2.0*M_PI/L;
	k7 = 4*M_PI*M_PI/(L*L);
	i = 0;
	pbox = &pkd->kdTop[ROOT];
	for (hx=-hReps;hx<=hReps;++hx) {
		for (hy=-hReps;hy<=hReps;++hy) {
			for (hz=-hReps;hz<=hReps;++hz) {
				h2 = hx*hx + hy*hy + hz*hz;
				if (h2 == 0) continue;
				if (h2 > fhCut*fhCut) continue;
				if (i == pkd->nMaxEwhLoop) {
					pkd->nMaxEwhLoop *= 2;
					pkd->ewt = realloc(pkd->ewt,pkd->nMaxEwhLoop*sizeof(EWT));
					assert(pkd->ewt != NULL);
					}
				pkd->ewt[i].hx = k5*hx;
				pkd->ewt[i].hy = k5*hy;
				pkd->ewt[i].hz = k5*hz;
				exph2 = exp(-k4*h2)/h2;
				qfac = 0.5*(pbox->mom.Qxx*hx*hx + pbox->mom.Qyy*hy*hy + 
							pbox->mom.Qzz*hz*hz) + pbox->mom.Qxy*hx*hy +
								pbox->mom.Qxz*hx*hz + pbox->mom.Qyz*hy*hz;
				dir3 = k3*exph2*(pbox->fMass - k7*qfac);
				pkd->ewt[i].k2 = k2*exph2*(pbox->fMass - k7*qfac);
				pkd->ewt[i].k3x = hx*dir3;
				pkd->ewt[i].k3y = hy*dir3;
				pkd->ewt[i].k3z = hz*dir3;
				++i;
				}
			}
		}
	pkd->nEwhLoop = i;
	}



