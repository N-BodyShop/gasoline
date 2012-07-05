#ifdef COLLISIONS

#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "pkd.h"
#include "collision.h" /*DEBUG not the nicest way of doing this*/

#ifdef SPECIAL_PARTICLES

static void
do_oblate(const FLOAT r0[3],struct oblate_s *oblate,
		  SPECIAL_PARTICLE_INFO *sInfo,FLOAT a[3])
{
	/*DEBUG: Kevin has a newer version that includes the effect of
	  body obliquity -- it's still in development. The version here
	  ignores obliquity. -- DCR 4/3/03*/

	/*
	 ** NOTE: assumes particle mass negligible compared to special particle,
	 ** so center of mass of system is center of special particle.
	 ** (Cf. Eqn (5-9) in Lang 1986).
	 */

	double R2,R4,dx,dy,dz,r2,ir2,ir,ir3,ir5,ir7,z2r2,cz1,cz2,cz3,j2,j4,cxy,cz;

	R2 = oblate->dRadEq*oblate->dRadEq;
	R4 = R2*R2;
	dx = r0[0] - sInfo->r[0];
	dy = r0[1] - sInfo->r[1];
	dz = r0[2] - sInfo->r[2];
	r2 = dx*dx + dy*dy + dz*dz;
	assert(r2 > 0.0);
	ir2 = 1/r2;
	ir = sqrt(ir2);
	ir3 = ir*ir2;
	ir5 = ir3*ir2;
	ir7 = ir5*ir2;
	z2r2 = dz*dz*ir2; /* sin^2(phi) */
	cz1 = 5*z2r2;
	cz2 = (21*z2r2 - 14)*z2r2 + 1;
	cz3 = (63*z2r2 - 70)*z2r2 + 15;
	j2 = 1.5*sInfo->fMass*oblate->J2*R2*ir5;
	j4 = 0.625*sInfo->fMass*oblate->J4*R4*ir7;
	cxy = j2*(cz1 - 1) + 3*j4*cz2;
	cz = j2*(cz1 - 3) + j4*cz3;
	a[0] += cxy*dx;
	a[1] += cxy*dy;
	a[2] += cz*dz;
	}

#ifdef OLD_COMET_FORCE
static void
do_force(struct force_s *force,SPECIAL_PARTICLE_INFO *sInfo,FLOAT r[3],FLOAT a[3])
{
#ifdef OLD_VERSION
	/*DEBUG for now, direct force on a line 135 deg cw from position vector*/
	/*      ROTATION ASSUMED TO BE IN XY PLANE*/
	double fx,fy,f;
	fx = -sInfo->r[0] + sInfo->r[1];
	fy = -sInfo->r[0] - sInfo->r[1];
	f = sqrt(fx*fx + fy*fy);
	fx *= force->dMag/f;
	fy *= force->dMag/f;
	a[0] += fx/sInfo->fMass;
	a[1] += fy/sInfo->fMass;
#endif
	/*DEBUG */
	double x,y,t;
	x = sInfo->r[0] - r[0]; /* assume xy plane; r[0,1] = com of rubble pile */
	y = sInfo->r[1] - r[1];
	t = acos(x/sqrt(x*x + y*y))*180.0/M_PI;
	if (y >= 0.0 && t <= 120.0) /* up to 120 deg counterclockwise from +x-axis */
		a[1] -= force->dMag/sInfo->fMass;
	else if (y <= 0.0 && t >= 60.0) /* up to 120 deg counterclockwise from -x-axis for symmetry (no net drift) */
		a[1] += force->dMag/sInfo->fMass;
	}
#endif

static void
do_noghostpert(const FLOAT r0[3],SPECIAL_MASTER_INFO *mInfo,
	       SPECIAL_PARTICLE_INFO *sInfo,double **rGhosts,int nGhosts,
	       FLOAT a[3])
{
  int i,j;
  double r[3],rmag,mir3;

  for (i=0;i<nGhosts;i++) {
    /* Calculate distance */
    rmag=0;
    for (j=0;j<3;j++) {
      r[j]=r0[j]-rGhosts[i][j];
      rmag+=r[j]*r[j];
    }
    rmag=sqrt(rmag);

    /* Subtract gravity */
    mir3=sInfo->fMass/(rmag*rmag*rmag);
    for (j=0;j<3;j++)
      a[j] += mir3*r[j]; /* Force is negative, must add */
  }
}

static void
getghostcoords(SPECIAL_MASTER_INFO *mInfo,SPECIAL_PARTICLE_INFO *sInfo,double **rGhosts,int nGhosts) {
  int i;
  int ix,iy,iz,nx,ny,nz;


  nx=(mInfo->dxPeriod < FLOAT_MAXVAL ? mInfo->nReplicas : 0);
  ny=(mInfo->dyPeriod < FLOAT_MAXVAL ? mInfo->nReplicas : 0);
  nz=(mInfo->dzPeriod < FLOAT_MAXVAL ? mInfo->nReplicas : 0);

  i=0;
  ix=-nx;
  iy=-ny;
  iz=-nz;

  do {
    if (ix != 0 || iy != 0 || iz != 0) {
      rGhosts[i][0]=sInfo->r[0]+ix*mInfo->dxPeriod; 
      rGhosts[i][1]=sInfo->r[1]+iy*mInfo->dyPeriod; 
      rGhosts[i][2]=sInfo->r[2]+iz*mInfo->dzPeriod;
#ifdef SLIDING_PATCH
    if (mInfo->dxPeriod < FLOAT_MAXVAL)
      rGhosts[i][1] -= 1.5*mInfo->dOmega*mInfo->dTime*mInfo->dxPeriod;
#endif /* SLIDING_PATCH */
    } else {
      i--; /* reduce i as 0,0,0 shouldn't be counted. */
    }

    /* This is actually elegant. */
    if (nz > 0) {
      iz++;
      if (iz > nz) {
	iz=-nz;
	if (ny > 0) {
	  iy++;
	  if (iy > ny) {
	    iy=-ny;
	    if (nx > 0) {
	      ix++;
	    }
	  }
	} else if (nx > 0) {
	  ix++;
	}
      }
    } else if (ny > 0) {
      iy++;
      if (iy > ny) {
	iy=-ny;
	if (nx > 0) {
	  ix++;
	}
      }
    } else if (nx > 0) {
      ix++;
    }
  } while (++i<nGhosts);
}

void
pkdGetSpecialParticles(PKD pkd,int nSpecial,int iId[],SPECIAL_MASTER_INFO *mInfo,
					   SPECIAL_PARTICLE_INFO *sInfo)
{
	/* retrieves current info for "special" particles */

	PARTICLE *p;
	int i,j,k;

	for (i=0;i<nSpecial;i++) {
		if (iId[i] == -1) {
			/* super special case: reference frame itself has special gravity */
			/* NOTE: all processors will do this -- should be OK...? */
			sInfo[i].iOrder = -2; /*DEBUG use symbols? -1 reserved... */
			sInfo[i].fMass = mInfo->dCentMass;
			sInfo[i].fRadius = 0.0; /* undefined */
			for (k=0;k<3;k++) sInfo[i].r[k] = 0;
			}
		else {
			for (j=0;j<pkdLocal(pkd);j++) {
				p = &pkd->pStore[j];
				if (p->iOrgIdx == iId[i]) {
					assert(sInfo[i].iOrder == -1); /* only one match allowed */
					sInfo[i].iOrder = p->iOrder;
					sInfo[i].fMass = p->fMass;
					sInfo[i].fRadius = RADIUS(p);
					for (k=0;k<3;k++) sInfo[i].r[k] = p->r[k];
					}
				}
			}
		}
	}

void
pkdDoSpecialParticles(PKD pkd,int nSpecial,SPECIAL_MASTER_INFO *mInfo,
					  SPECIAL_PARTICLE_DATA sData[],
					  SPECIAL_PARTICLE_INFO sInfo[],FLOAT aFrame[])
{
	/* applies effects of "special" particles on other particles */

	PARTICLE *p;
	int i,j;
	int nGhosts=0;
	double **rGhosts=NULL;
#ifdef OLD_COMET_STUFF
	/*THIS STUFF #ifdef'D OUT TO SIMPLIFY MERGING OF RORY'S VERSION
	  WITH UMD SVN VERSION 5/23/06*/
	int nLocal=pkdLocal(pkd);

	/*DEBUG some (all?) of this stuff is for the comet spin-up problem*/
	/*DEBUG this first stuff is for do_force(): need to know COM position of rubble pile -- won't work in parallel!!*/
	FLOAT r[3],rc[3];
	int n,nn;
	do {
		r[0] = r[1] = r[2] = 0.0;
		for (i=n=0;i<nLocal;i++) {
			p = &pkd->pStore[i];
			if (p->iColor != 2) {
				r[0] += p->r[0]; /*DEBUG assumes equal-mass particles...*/
				r[1] += p->r[1];
				r[2] += p->r[2];
				++n;
			}
		}
		assert(n > 0);
		rc[0] = r[0]/n;
		rc[1] = r[1]/n;
		rc[2] = r[2]/n;
		/* loop again to remove any particles far from main clump */
		for (i=nn=0;i<nLocal;i++) {
			p = &pkd->pStore[i];
			if (p->iColor == 2) continue;
			if (sqrt((p->r[0] - rc[0])*(p->r[0] - rc[0]) +
					 (p->r[1] - rc[1])*(p->r[1] - rc[1]) +
					 (p->r[2] - rc[2])*(p->r[2] - rc[2])) > 4.0e-8) { /*DEBUG 8 km */
				if (p->iColor != 3) {
					(void) fprintf(stderr,"special particle %i escaped\n",p->iOrder);
/*DEBUG					assert(0);*/
				}
				p->iColor = 2;
				continue;
			}
			++nn;
		}
		assert(nn > 0);
	} while (nn < n);

	assert(nn == n);
#endif /*OLD_COMET_STUFF*/

	for (i=0;i<nSpecial;i++) {
	  if (sData[i].iType & SPECIAL_NOGHOSTPERT) {
	    nGhosts=(mInfo->dxPeriod < FLOAT_MAXVAL ? 2*mInfo->nReplicas+1 : 1)*
	      (mInfo->dyPeriod < FLOAT_MAXVAL ? 2*mInfo->nReplicas+1 : 1)*
	       (mInfo->dzPeriod < FLOAT_MAXVAL ? 2*mInfo->nReplicas+1 : 1)-1;
	    if (nGhosts > 0) {
	      rGhosts=malloc(nGhosts*sizeof(double*));
	      for (j=0;j<nGhosts;j++) {
		rGhosts[j]=malloc(3*sizeof(double));
	      }
	    getghostcoords(mInfo,&sInfo[i],rGhosts,nGhosts);
	    }
	  }
	  for (j=0;j<pkdLocal(pkd);j++) {
	    p = &pkd->pStore[j];
#ifdef OLD_COMET_FORCE
			if (p->iOrder == sInfo[i].iOrder) {
				if (sData[i].iType & SPECIAL_FORCE) {
					/* special case: force applies ONLY to special particle */
					double z = sInfo[i].r[2] - rc[2];
					if (fabs(z) > 3*sInfo[i].fRadius) {
						(void) fprintf(stderr,"%i: Strong off-axis torque\n",p->iOrder);
						sInfo[i].iOrder = -1;
						p->iColor = 16;
						continue;
						}
					do_force(&sData[i].force,&sInfo[i],rc,p->a);
					}
				else
					continue;
				}
#else
			continue;
#endif
			if (sData[i].iType & SPECIAL_OBLATE) {
				do_oblate(p->r,&sData[i].oblate,&sInfo[i],p->a);
				}
			if (sData[i].iType & SPECIAL_GR) {
				assert(0); /* not implemented yet */
				}
			if (sData[i].iType & SPECIAL_NOGHOSTPERT && nGhosts > 0) {
			  do_noghostpert(p->r,mInfo,&sInfo[i],rGhosts,nGhosts,p->a);
			        }
		        }
		if (mInfo->bNonInertial) {
			FLOAT r[3]={0,0,0};
			int k;
			/* get acceleration on frame */
			if (sInfo[i].iOrder < 0)
				continue; /* skip if frame center is also special */
			for (k=0;k<3;k++) aFrame[k] = 0;
			if (sData[i].iType & SPECIAL_OBLATE) {
				do_oblate(r,&sData[i].oblate,&sInfo[i],aFrame);
				}
			if (sData[i].iType & SPECIAL_GR) {
				assert(0); /* not implemented yet */
				}
			if (sData[i].iType & SPECIAL_NOGHOSTPERT) {
			        assert(0); /* Already in patch frame */
			        }
			}
		}
	  if (rGhosts != NULL) {
	    free(rGhosts);
	  }
	}

#endif /* SPECIAL_PARTICLES */

#endif /* COLLISION */
