#ifdef COLLISIONS

#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "pkd.h"
#include "collision.h" /*DEBUG not the nicest way of doing this*/

#ifdef SPECIAL_PARTICLES

static void
oblate(const FLOAT r0[3],double dRadEq,double J2,double J4,FLOAT fMass,
	   const FLOAT r[3],FLOAT a[3])
{
	/*
	 ** NOTE: assumes particle mass negligible compared to special particle,
	 ** so center of mass of system is center of special particle.
	 */

	double R2,R4,dx,dy,dz,r2,ir2,ir,ir3,ir5,ir7,z2r2,cz1,cz2,cz3,j2,j4,cxy,cz;

	R2 = dRadEq*dRadEq;
	R4 = R2*R2;
	dx = r0[0] - r[0];
	dy = r0[1] - r[1];
	dz = r0[2] - r[2];
	r2 = dx*dx + dy*dy + dz*dz;
	assert(r2 > 0);
	ir2 = 1/r2;
	ir = sqrt(ir2);
	ir3 = ir*ir2;
	ir5 = ir3*ir2;
	ir7 = ir5*ir2;
	z2r2 = dz*dz*ir2; /* sin^2(phi) */
	cz1 = 5*z2r2;
	cz2 = (21*z2r2 - 14)*z2r2 + 1;
	cz3 = (63*z2r2 - 70)*z2r2 + 15;
	j2 = 1.5*fMass*J2*R2*ir5;
	j4 = 0.625*fMass*J4*R4*ir7;
	cxy = j2*(1 - cz1) + 3*j4*cz2;
	cz = j2*(3 - cz1) + j4*cz3;
	a[0] += cxy*dx;
	a[1] += cxy*dy;
	a[2] += cz*dz;
	}

void
pkdGetSpecialParticles(PKD pkd,int nSpecial,int iId[],double dCentMass,
					   SPECIAL_PARTICLE_INFO sInfo[])
{
	/* retrieves current info for "special" particles */

	PARTICLE *p;
	int i,j,k;

	for (i=0;i<nSpecial;i++) {
		if (iId[i] == -1) {
			/* super special case: reference frame itself has special gravity */
			/* NOTE: all processors will do this -- should be OK...? */
			sInfo[i].iOrder = -2; /*DEBUG use symbols? -1 reserved... */
			sInfo[i].fMass = dCentMass;
			for (k=0;k<3;k++) sInfo[i].r[k] = 0;
			}
		else {
			for (j=0;j<pkdLocal(pkd);j++) {
				p = &pkd->pStore[j];
				if (p->iOrgIdx == iId[i]) {
					assert(sInfo[i].iOrder == -1); /* only one match allowed */
					sInfo[i].iOrder = p->iOrder;
					sInfo[i].fMass = p->fMass;
					for (k=0;k<3;k++) sInfo[i].r[k] = p->r[k];
					}
				}
			}
		}
	}

void
pkdDoSpecialParticles(PKD pkd,int nSpecial,int bNonInertial,
					  SPECIAL_PARTICLE_DATA sData[],
					  SPECIAL_PARTICLE_INFO sInfo[],FLOAT aFrame[])
{
	/* applies effects of "special" particles on other particles */

	PARTICLE *p;
	int i,j;

	for (i=0;i<nSpecial;i++) {
		for (j=0;j<pkdLocal(pkd);j++) {
			p = &pkd->pStore[j];
			if (p->iOrder == sInfo[i].iOrder)
				continue;
			if (sData[i].iType & SPECIAL_OBLATE) {
				oblate(p->r,sData[i].oblate.dRadEq,sData[i].oblate.J2,
					   sData[i].oblate.J4,sInfo[i].fMass,sInfo[i].r,p->a);
				}
			if (sData[i].iType & SPECIAL_GR) {
				assert(0); /* not implemented yet */
				}
			}
		if (bNonInertial) {
			FLOAT r[3]={0,0,0};
			int k;
			/* get acceleration on frame */
			if (sInfo[i].iOrder < 0)
				continue; /* skip if frame center is also special */
			for (k=0;k<3;k++) aFrame[k] = 0;
			if (sData[i].iType & SPECIAL_OBLATE) {
				oblate(r,sData[i].oblate.dRadEq,sData[i].oblate.J2,
					   sData[i].oblate.J4,sInfo[i].fMass,sInfo[i].r,aFrame);
				}
			if (sData[i].iType & SPECIAL_GR) {
				assert(0); /* not implemented yet */
				}
			}
		}
	}

#endif /* SPECIAL_PARTICLES */

#endif /* COLLISION */
