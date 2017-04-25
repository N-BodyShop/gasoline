
#ifdef COLLISIONS

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "collision.h"

#ifdef RUBBLE_ZML
#include "rubble.h"
#endif

#define BOUNCE_OK 0
#define NEAR_MISS 1

void pkdNextCollision(PKD pkd,double *dt,int *iOrder1,int *iOrder2)
{
	/*
	 ** Returns time and iOrder of particles for earliest predicted
	 ** collision on this processor. Argument initialization must
	 ** occur in calling routine.
	 */

	PARTICLE *p;
	int i,nLocal = pkdLocal(pkd);

	for (i=0;i<nLocal;i++) {
		p = &pkd->pStore[i];
		if (!TYPEQueryACTIVE(p)) continue; /* skip over inactive particles */
		if (p->iOrder < 0) continue; /* skip over deleted particles */
		if (p->dtCol < *dt) {
			*dt = p->dtCol;
			*iOrder1 = p->iOrder;
			*iOrder2 = p->iOrderCol;
			}
		}
	}

double LastKickTime(int iRung,double dBaseTime,double dTimeNow) 
{
	/*
	 ** Determines the last time a particle received a force update.
	 */

	double dRungTime;
	int nRungSteps;

	dRungTime = dBaseTime/(1<<iRung);
	nRungSteps = (int)(dTimeNow/dRungTime);

	return dTimeNow - (dRungTime*nRungSteps);
}

void MergerReorder(PKD pkd,const COLLIDER *c1,const COLLIDER *c2,const COLLIDER *c,
				   COLLIDER *cOut,double dt,const COLLISION_PARAMS *CP,
				   int bDiagInfo,const FLOAT fOffset[],
#ifdef RUBBLE_ZML
				   double dMassInDust,
#endif
#ifdef SLIDING_PATCH
				   FLOAT fShear,
#endif
				   int bPeriodic)
{
  /*
  ** NOTE:
  **   c1 = (input) pointer to 1st particle's pre-collision info
  **   c2 = (input) pointer to 2nd particle's pre-collision info
  **    c = (input) pointer to array of post-collision particle
  **        info, adjusted to start of step
  ** cOut = (output) pointer to array of post-collision particle
  **        info, intended for collision log
  */

  const PARTICLE_ID *pMrg=NULL,*pDel=NULL,*pOth=NULL;

  assert(c1 != NULL);
  assert(c2 != NULL);
  assert(c != NULL);
  assert(cOut != NULL);
#ifdef RUBBLE_ZML
  assert(CP != NULL);  /* Only used for the Rubble */
#endif

  if (c1->id.iPid == pkd->idSelf) { /* local particle */
	/*
	** Keep this particle if it has larger mass, or, if both
	** particles are the same mass, keep this particle if it
	** has a lower original index, or, if both particles have
	** same original index, keep this particle if it has a
	** lower processor number, or, if both particles are on
	** same processor, keep particle with lower iOrder.
	*/
	if (c1->fMass > c2->fMass ||
		(c1->fMass == c2->fMass &&
		 (c1->id.iOrgIdx < c2->id.iOrgIdx ||
		  (c1->id.iOrgIdx == c2->id.iOrgIdx &&
		   (c1->id.iPid < c2->id.iPid ||
			(c1->id.iPid == c2->id.iPid &&
			 c1->id.iOrder < c2->id.iOrder))))))
	  pMrg = &c1->id; /* new centre-of-mass particle */
	else
	  pDel = &c1->id; /* this particle is removed */
	pOth = &c2->id;
  }
  if (c2->id.iPid == pkd->idSelf) { /* same thing for other particle */
	if (c2->fMass > c1->fMass ||
		(c2->fMass == c1->fMass &&
		 (c2->id.iOrgIdx < c1->id.iOrgIdx ||
		  (c2->id.iOrgIdx == c1->id.iOrgIdx &&
		   (c2->id.iPid < c1->id.iPid ||
			(c2->id.iPid == c1->id.iPid &&
			 c2->id.iOrder < c1->id.iOrder))))))
	  pMrg = &c2->id;
	else
	  pDel = &c2->id;
	pOth = &c1->id;
  }
  /* 
  ** Store or delete particle as appropriate, and make
  ** sure master knows ID of surviving particle for
  ** tracking purposes.  INT_MAX is used to indicate that
  ** the other collider is now gone.
  */
  if (pMrg != NULL) {
	PutColliderInfo(pkd, c,INT_MAX,&pkd->pStore[pMrg->iIndex],dt);
	if (bDiagInfo) cOut->id = *pMrg; /* struct copy */
#ifdef RUBBLE_ZML
	{
	  /*DEBUG -- CHECK FOR ABNORMAL DENSITY*/
	  double rho;
	  rho = c->fMass/(4.0/3.0*M_PI*CUBE(c->fRadius));
	  if (rho > 4.0e6 && pkd->pStore[pMrg->iIndex].iColor == PLANETESIMAL) /* f=1 */
		(void) printf("BAD MERGED DENSITY: iOrder=%i iOrgIdx=%i rho=%g\n",pMrg->iOrder,pMrg->iOrgIdx,rho);
	}
	if (CP->iRubForcedOutcome == RUB_FORCED_NONE) {
	  /*
	  ** Store any dust created by an interpolated collision into
	  ** particle dust parameter.
	  **DEBUG Should this be a +=?
	  */
	  pkd->pStore[pMrg->iIndex].dDustMass = dMassInDust;
	}
#endif
  }
  if (pDel != NULL) {
	pkdDeleteParticle(pkd,&pkd->pStore[pDel->iIndex]);
	if (bDiagInfo && pMrg == NULL) cOut->id = *pOth; /* may need merger info */
  }
}

void SetMergerRung(const COLLIDER *c1,const COLLIDER *c2,COLLIDER *c,
				   double dBaseStep,double dTimeNow,int iTime0)
{
  double lastkick1,lastkick2;
  FLOAT m1,m2,m;
  int k;

  m1 = c1->fMass;
  m2 = c2->fMass;
  m = m1 + m2;

  if (c1->iRung!=c2->iRung) {
    lastkick1=LastKickTime(c1->iRung,dBaseStep,(dTimeNow-iTime0));
    lastkick2=LastKickTime(c2->iRung,dBaseStep,(dTimeNow-iTime0));
    if (m1 > m2) {
      c->iRung=c1->iRung;
      if (lastkick1>lastkick2) {
	for (k=0;k<3;k++) {
	  c->v[k]+=(m1/m)*c1->a[k]*(lastkick1-lastkick2);
	}
      } else {
	for (k=0;k<3;k++) {
	  c->v[k]+=(m2/m)*c2->a[k]*(lastkick2-lastkick1);
	}
      }
    } else { /* c2->fMass > c1->fMass */
      c->iRung=c2->iRung;
      if (lastkick1>lastkick2) {
	for (k=0;k<3;k++) {
	  c->v[k]+=(m1/m)*c1->a[k]*(lastkick1-lastkick2);
	}
      } else {
	for (k=0;k<3;k++) {
	  c->v[k]+=(m2/m)*c2->a[k]*(lastkick2-lastkick1);
	}
      }
    }
  } else { /* c1->iRung == c2->iRung */
    c->iRung=c1->iRung;
  }
  /* This was the old way of fixing the resultant particle's rung
   * c->iRung = (c2->fMass > c1->fMass ? c2->iRung : c1->iRung); */
}

void CalcOffset(PKD pkd,COLLIDER *c1,COLLIDER *c2, FLOAT fOffset[]
#ifdef SLIDING_PATCH
				,FLOAT *fShear
#endif
				)
{
  int i;

  for (i=0;i<3;i++) {
    if (c2->r[i] - c1->r[i] >= 0.5*pkd->fPeriod[i]) {
      fOffset[i] -= pkd->fPeriod[i];
      c2->r[i] -= pkd->fPeriod[i];
    }
    else if (c2->r[i] - c1->r[i] < - 0.5*pkd->fPeriod[i]) {
      fOffset[i] += pkd->fPeriod[i];
      c2->r[i] += pkd->fPeriod[i];
    }
#ifdef SLIDING_PATCH
    if (i == 0) {
      if (fOffset[0] < 0) {
	*fShear = 1.5*pkd->PP->dOrbFreq*pkd->fPeriod[0];
	fOffset[1] = SHEAR(-1,pkd->dTime,pkd->PP);
      }
      else if (fOffset[0] > 0) {
	*fShear = - 1.5*pkd->PP->dOrbFreq*pkd->fPeriod[0];
	fOffset[1] = SHEAR(1,pkd->dTime,pkd->PP);
      }
      c2->r[1] += fOffset[1];
      c2->v[1] += *fShear;
    }
#endif
  }
}

void pkdFindTightestBinary(PKD pkd,double *dBindEn,int *iOrder1,int *iOrder2,int *n)
{
  /* Find the local binary with the highest binding energy. */

  PARTICLE *p;
  int i,nLocal = pkdLocal(pkd);

  for (i=0;i<nLocal;i++) {
    p = &pkd->pStore[i];
    if (!TYPEQueryACTIVE(p)) continue;
    if (p->iOrder<0) continue;
    if (p->dtCol < *dBindEn) { /* This particle is more tightly bound */
      *n=2;
      *dBindEn = p->dtCol;
      *iOrder1 = p->iOrder;
      *iOrder2 = p->iOrderCol;
    }
  }
}

void pkdMergeBinary(PKD pkd,const COLLIDER *pc1,const COLLIDER *pc2,COLLIDER *cOut,
					int bPeriodic,double dBaseStep,double dTimeNow,int iTime0,
					double dDensity,int *bool)
{
  COLLIDER c1,c2,c;
  FLOAT rc1[3],rc2[3],vc1[3],vc2[3],ac1[3],ac2[3];
  FLOAT m1,m2,m;
  int k,bDiagInfo;
  FLOAT rsq1=0,rsq2=0,vsq1=0,vsq2=0,r2,r1,ang_mom,ke;
  const double dDenFac = 4.0/3*M_PI;  
  FLOAT fOffset[3] = {0,0,0};
#ifdef SLIDING_PATCH
  FLOAT fShear=0.0;
#endif

  c1 = *pc1;
  c2 = *pc2;

  /* First find center of mass position and velocity */

  m1=c1.fMass;
  m2=c2.fMass;
  m=m1+m2;
  r1=c1.fRadius;
  r2=c2.fRadius;

  bDiagInfo = (c1.id.iPid == pkd->idSelf); 
 
  if (bPeriodic) 
    CalcOffset(pkd,&c1,&c2,fOffset
#ifdef SLIDING_PATCH
	       ,&fShear
#endif
	       );

  for (k=0;k<3;k++) {
    c.r[k] = (m1*c1.r[k] + m2*c2.r[k])/m;
    c.v[k] = (m1*c1.v[k] + m2*c2.v[k])/m;
    c.a[k] = (m1*c1.a[k] + m2*c2.a[k])/m;
    c.w[k] = (m1*c1.w[k] + m2*c2.w[k])/m;/*DEBUG can't do this with spins, can you?--DCR*/
    rc1[k] = c1.r[k] - c.r[k];
    rc2[k] = c2.r[k] - c.r[k];
    vc1[k] = c1.v[k] - c.v[k];
    vc2[k] = c2.v[k] - c.v[k];
    ac1[k] = c1.a[k] - c.a[k];
    ac2[k] = c2.a[k] - c.a[k];
    rsq1+=rc1[k]*rc1[k];
    vsq1+=vc1[k]*vc1[k];
    rsq2+=rc2[k]*rc2[k];
    vsq2+=vc2[k]*vc2[k];
  }
  c.fMass=m;
  if (dDensity)
    c.fRadius = pow(m/(dDenFac*dDensity),1.0/3);
  else
    c.fRadius = (m2 > m1 ? r2*pow(m/m2,1.0/3) : r1*pow(m/m1,1.0/3));

  /* This much angular momentum and kinetic energy are lost. *//*DEBUG why is ang mom lost?--DCR*/
  ang_mom = m1*sqrt(rsq1)*sqrt(vsq1) + m2*sqrt(rsq2)*sqrt(vsq1);
  ke=0.5*(m1*vsq1 + m2*vsq2);

  /* For now both colliders are the same type */
  c.iColor = c1.iColor;

  SetMergerRung(&c1,&c2,&c,dBaseStep,dTimeNow,iTime0);
  MergerReorder(pkd,&c1,&c2,(const COLLIDER *) &c,&c,-1,NULL/*DEBUG--hack*/,bDiagInfo,fOffset,
#ifdef RUBBLE_ZML
				0.0/*DEBUG--hack*/,
#endif
#ifdef SLIDING_PATCH 
		fShear,
#endif
		bPeriodic);

  *cOut = c;

  /* For parallel purposes, let's now set the boolean to 1 */
  *bool=1;
}

void pkdGetColliderInfo(PKD pkd,int iOrder,COLLIDER *c)
{
	/*
	 ** Returns collider info for particle with matching iOrder.
	 */

	PARTICLE *p;
	int i,j,nLocal = pkdLocal(pkd);

	for (i=0;i<nLocal;i++) {
		p = &pkd->pStore[i];
		if (p->iOrder == iOrder) {
			assert(TYPEQueryACTIVE(p));
			c->id.iPid = pkd->idSelf;
			c->id.iOrder = iOrder;
			c->id.iIndex = i;
			c->id.iOrgIdx = p->iOrgIdx;
			c->fMass = p->fMass;
			c->fRadius = RADIUS(p); /* twice softening radius */
			for (j=0;j<3;j++) {
				c->r[j] = p->r[j];
				c->v[j] = p->v[j];
				c->w[j] = p->w[j];
				c->a[j] = p->a[j];
				}
			c->iColor = p->iColor;
			c->dt = p->dt;
			c->iRung = p->iRung;
			c->bTinyStep = p->bTinyStep;
			return;
			}
		}
	}

void PutColliderInfo(PKD pkd, const COLLIDER *c,int iOrder2,PARTICLE *p,
		     double dt)
{
	/*
	 ** Stores collider info in particle structure (except id & color).
	 ** Also, dt is stored in dtPrevCol for inelastic collapse checks.
	 **
	 ** dt is set to -1 in a binary merge.
	 **
	 ** NOTE: Because colliding particles have their positions traced back to
	 ** the start of the step using their NEW velocities, it is possible for
	 ** a neighbour to no longer lie inside a previous collider's search ball.
	 ** To fix this, the ball radius is expanded by a conservative estimate of
	 ** the displacement amount (using the "Manhattan metric").  This is
	 ** needed so we can use resmooth instead of smooth to save time after
	 ** a collision occurs during the search interval.
	 */

	int i;
	double r;

	p->fMass = c->fMass;
	p->fSoft = SOFT(c); /* half particle radius */
	r = fabs(c->r[0] - p->r[0]) +
		fabs(c->r[1] - p->r[1]) +
		fabs(c->r[2] - p->r[2]); /* measured AFTER collider traced back */
	for (i=0;i<3;i++) {
		p->r[i] = c->r[i];
		p->v[i] = c->v[i];
		p->w[i] = c->w[i];
		p->a[i] = c->a[i];
#ifdef NEED_VPRED
		p->vPred[i] = c->v[i] - dt*p->a[i];
#endif
		}
#ifdef SLIDING_PATCH
	p->dPy = p->v[1] + 2.0*pkd->PP->dOrbFreq*p->r[0];
#endif
	p->iRung = c->iRung;
	p->fBall2 += 2*sqrt(p->fBall2)*r + r*r;
	p->dtPrevCol = dt; /* N.B. this is set to -1 in a binary merge */
	p->iPrevCol = iOrder2; /* stored to avoid false collisions */
}

void pkdSetBall(PKD pkd, double dDelta, double fac)
{
  PARTICLE *p;
  double vsq;
  int i,j,k,nLocal = pkdLocal(pkd);

  for (i=0;i<nLocal;i++) {
    p = &pkd->pStore[i];
    vsq=0;
    for (j=0;j<3;j++) 
      vsq+=p->v[j]*p->v[j];
    /* search radius = (fac*dDelta*v + 4*softlength)^2 */
    p->fBall2=fac*fac*dDelta*dDelta*vsq + 4.0*fac*dDelta*sqrt(vsq)*p->fSoft + 16.0*p->fSoft*p->fSoft;
    if (pkd->fPeriod[0]<FLOAT_MAXVAL || pkd->fPeriod[1]<FLOAT_MAXVAL ||
	pkd->fPeriod[2]<FLOAT_MAXVAL) {
      for (k=0;k<3;k++) 
	assert(p->fBall2 < 0.25*pkd->fPeriod[k]*pkd->fPeriod[k]);
    }
    p->cpStart = 0;
  }
}

void pkdMerge(PKD pkd,const COLLIDER *c1,const COLLIDER *c2,double dDensity,
			  COLLIDER **cOut,int *pnOut,int iTime0,double dBaseStep,
			  double dTimeNow)
{
	/*
	 ** Merges colliders into new centre-of-mass particle (cOut[0])
	 ** with radius determined by the supplied fixed bulk density.
	 ** If dDensity is zero, the radius of the merged particle is
	 ** instead set by the bulk density of the largest collider.
	 */

	const double dDenFac = 4.0/3*M_PI;

	COLLIDER *c;
	FLOAT com_pos[3],com_vel[3],rc1[3],rc2[3],vc1[3],vc2[3],ang_mom[3];
	FLOAT com_acc[3],ac1[3],ac2[3];
	FLOAT m1,m2,m,r1,r2,r,i1,i2,i;
	int k;

	m1 = c1->fMass;
	m2 = c2->fMass;
	m = m1 + m2;
	r1 = c1->fRadius;
	r2 = c2->fRadius;
	if (dDensity)
		r = pow(m/(dDenFac*dDensity),1.0/3);
	else
		r = pow(r1*r1*r1 + r2*r2*r2,1.0/3); /* conserves volume *//*DEBUG used to be r = (m2 > m1 ? r2*pow(m/m2,1.0/3) : r1*pow(m/m1,1.0/3));*/
	i1 = 0.4*m1*r1*r1;
	i2 = 0.4*m2*r2*r2;
	i = 0.4*m*r*r;

	for (k=0;k<3;k++) {
		com_pos[k] = (m1*c1->r[k] + m2*c2->r[k])/m;
		rc1[k] = c1->r[k] - com_pos[k];
		rc2[k] = c2->r[k] - com_pos[k];
		com_vel[k] = (m1*c1->v[k] + m2*c2->v[k])/m;
		vc1[k] = c1->v[k] - com_vel[k];
		vc2[k] = c2->v[k] - com_vel[k];
		com_acc[k] = (m1*c1->a[k] + m2*c2->a[k])/m;
		ac1[k] = c1->a[k] - com_acc[k];
		ac2[k] = c2->a[k] - com_acc[k];
	}

	ang_mom[0] = m1*(rc1[1]*vc1[2] - rc1[2]*vc1[1]) + i1*c1->w[0] +
		         m2*(rc2[1]*vc2[2] - rc2[2]*vc2[1]) + i2*c2->w[0];
	ang_mom[1] = m1*(rc1[2]*vc1[0] - rc1[0]*vc1[2]) + i1*c1->w[1] +
				 m2*(rc2[2]*vc2[0] - rc2[0]*vc2[2]) + i2*c2->w[1];
	ang_mom[2] = m1*(rc1[0]*vc1[1] - rc1[1]*vc1[0]) + i1*c1->w[2] +
				 m2*(rc2[0]*vc2[1] - rc2[1]*vc2[0]) + i2*c2->w[2];

	*pnOut = 1;
	*cOut = (COLLIDER *) malloc(*pnOut*sizeof(COLLIDER));
	assert(*cOut != NULL);

	c = &(*cOut)[0];

	/* Note: id info set in pkdDoCollision(), used only for log */

	c->fMass = m;
	c->fRadius = r;

	for (k=0;k<3;k++) {
	  c->r[k] = com_pos[k];
	  c->v[k] = com_vel[k];
	  c->w[k] = ang_mom[k]/i;
	  c->a[k] = com_acc[k];
	}

	/* Set merger's timestep to iRung of largest mass. */
	/* We have to do some messy logic to correct for the possibility
	   of the two colliders being on different rungs, and hence having
	   received their last kick at different times. */
	SetMergerRung(c1,c2,c,dBaseStep,dTimeNow,iTime0);
}

int pkdBounce(PKD pkd,const COLLIDER *c1,const COLLIDER *c2,
			  double dEpsN,double dEpsT,int bFixCollapse,
			  COLLIDER **cOut,int *pnOut)
{
	/* Bounces colliders, preserving particle order */

	COLLIDER *co1,*co2;
	FLOAT n[3],s1[3],s2[3],v[3],s[3],u[3],un[3],ut[3],p[3],q[3];
	FLOAT m1,m2,m,r1,r2,i1,i2,mu,alpha,beta;
	FLOAT a,b,c,d;
	int i;

/*DEBUG verbose EpsN, EpsT output
	(void) printf("e_n = %g e_t = %g\n",dEpsN,dEpsT);
*/

	*pnOut = 2;
	*cOut = (COLLIDER *) malloc(*pnOut*sizeof(COLLIDER));
	assert(*cOut != NULL);

	(*cOut)[0] = *c1;
	(*cOut)[1] = *c2;

	m1 = c1->fMass;
	m2 = c2->fMass;
	m = m1 + m2;
	r1 = c1->fRadius;
	r2 = c2->fRadius;
	i1 = 0.4*m1*r1*r1;
	i2 = 0.4*m2*r2*r2;
	mu = m1*m2/m;
	alpha = 2.5*(1/m1 + 1/m2);
	beta = 1/(1 + alpha*mu);

	a = 0;
	for (i=0;i<3;i++) {
		n[i] = (c2->r[i] - c1->r[i]);
		a += n[i]*n[i];
		}
	a = 1/sqrt(a);
	for (i=0;i<3;i++)
		n[i] *= a;

	s1[0] = r1*(c1->w[1]*n[2] - c1->w[2]*n[1]);
	s1[1] = r1*(c1->w[2]*n[0] - c1->w[0]*n[2]);
	s1[2] = r1*(c1->w[0]*n[1] - c1->w[1]*n[0]);

	s2[0] = r2*(c2->w[2]*n[1] - c2->w[1]*n[2]);
	s2[1] = r2*(c2->w[0]*n[2] - c2->w[2]*n[0]);
	s2[2] = r2*(c2->w[1]*n[0] - c2->w[0]*n[1]);

	for (i=0;i<3;i++) {
		v[i] = c2->v[i] - c1->v[i];
		s[i] = s2[i] - s1[i];
		}

	for (i=0;i<3;i++)
		u[i] = v[i] + s[i];

	a = u[0]*n[0] + u[1]*n[1] + u[2]*n[2];
	if (a >= 0) {
		static int bGiveWarning = 1;
		if (bGiveWarning) {
			(void) fprintf(stderr,"WARNING: %i & %i -- near miss? (a = %g)\n",
						   c1->id.iOrder,c2->id.iOrder,a);
			bGiveWarning = 0;
			}

		if (bFixCollapse)
			return NEAR_MISS; /* particles remain unchanged */
		else
			assert(0); /* near miss not allowed */
		}
	for (i=0;i<3;i++) {
		un[i] = a*n[i];
		ut[i] = u[i] - un[i];
		}

#ifdef SAND_PILE
	if (pkd->pStore[c1->id.iIndex].bStuck) {
		mu = m2;
		beta = 2.0/7;
		}
	else if (pkd->pStore[c2->id.iIndex].bStuck) {
		mu = m1;
		beta = 2.0/7;
		}
#endif

	a = (1 + dEpsN);
	b = beta*(1 - dEpsT);
	for (i=0;i<3;i++)
		p[i] = a*un[i] + b*ut[i];

	a = mu*b;
	q[0] = a*(n[1]*u[2] - n[2]*u[1]);
	q[1] = a*(n[2]*u[0] - n[0]*u[2]);
	q[2] = a*(n[0]*u[1] - n[1]*u[0]);

	a =   m2/m;
	b = - m1/m;
	c = r1/i1;
	d = r2/i2;

#ifdef SAND_PILE
	if (pkd->pStore[c1->id.iIndex].bStuck) {
		a = 0;
		b = -1;
		c = 0;
		d = r2/i2;
		}
	else if (pkd->pStore[c2->id.iIndex].bStuck) {
		a = 1;
		b = 0;
		c = r1/i1;
		d = 0;
		}
#endif

	co1 = &(*cOut)[0];
	co2 = &(*cOut)[1];

	for (i=0;i<3;i++) {
		co1->v[i] += a*p[i];
		co2->v[i] += b*p[i];
		co1->w[i] += c*q[i];
		co2->w[i] += d*q[i];
		}

/*DEBUG -- dT check
	{
	double dT;
	dT =
		- mu*(v[0]*p[0] + v[1]*p[1] + v[2]*p[2]) +
			0.5*mu*(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]) +
				(r1*c1->w[0] + r2*c2->w[0])*q[0] +
					(r1*c1->w[1] + r2*c2->w[1])*q[1] +
						(r1*c1->w[2] + r2*c2->w[2])*q[2] +
							0.5*alpha*(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]);
	(void) printf("COLLISION %i & %i e_n=%f e_t=%f dT %e\n",
				  c1->id.iOrder,c2->id.iOrder,dEpsN,dEpsT,dT);
	}
*/

	return BOUNCE_OK;
	}

void pkdFrag(PKD pkd,const COLLIDER *c1,const COLLIDER *c2,
			 COLLIDER **cOut,int *pnOut)
{
	/*
	 ** Fragments colliders into at most MAX_NUM_FRAG pieces.
	 ** Not implemented yet. Development notes:
	 ** - be sure new particles are ACTIVE
	 ** - may need to assert(*pnOut <= MAX_NUM_FRAG)
	 ** - remember to set id info for logging purposes
	 ** - and decide on starting iRung values
	 ** - and assign colors
	 */

	assert(0);
	*cOut = NULL;
	*pnOut = 0;
	}

void pkdDoCollision(PKD pkd,double dt,const COLLIDER *pc1,const COLLIDER *pc2,
					int bPeriodic,int iTime0,double dBaseStep,double dTimeNow,
					const COLLISION_PARAMS *CP,int *piOutcome,double *dT,
					COLLIDER *cOut,int *pnOut)
{
	/*
	 ** Processes collision by advancing particle coordinates to impact
	 ** time, determining the collision outcome, and moving the resulting
	 ** particles back to their (new) start-of-step positions. Diagnostic
	 ** info is returned if the first collider resides on the processor
	 ** (the second collider may not be local, in which case this routine
	 ** is called twice in parallel, once for each particle). Local copies
	 ** of the collider info are used because the particle coordinates need
	 ** to be advanced to the moment of collision. With periodic boundary
	 ** conditions, the second collider may be offset by one period as well.
	 ** The offset is removed before calling PutColliderInfo().
	 */

	COLLIDER c1,c2,*c = NULL;
	FLOAT fOffset[3] = {0.0,0.0,0.0};
	double v2,ve2,dImpactEnergy;
	int bDiagInfo,iOutcome,i,j,n; /*DEBUG would bReportInfo be better?*/
#ifdef SLIDING_PATCH
	FLOAT fShear = 0.0,ct,st,wx,wy;
#endif

/*DEBUG verbose collision output
	(void) printf("COLLISION %i & %i (dt = %.16e)\n",
				  pc1->id.iOrder,pc2->id.iOrder,dt);
*/
#ifdef RUBBLE_ZML
	double dMassInDust;
#endif

	/* Get local copies of collider info */

	c1 = *pc1; /* struct copy */
	c2 = *pc2;

	bDiagInfo = (c1.id.iPid == pkd->idSelf);

	if (bDiagInfo && dT) *dT = 0;

#ifdef SAND_PILE
#ifdef TUMBLER
	if (c2.id.iOrder < 0) { /* wall collision */
		const WALLS *w = &CP->walls;
		WALL mywall = w->wall[-c2.id.iOrder - 1];

		assert(c1.id.iPid == pkd->idSelf);
		assert(pkd->pStore[c1.id.iIndex].bStuck == 0);
		for (i=0;i<3;i++) c1.r[i] += c1.v[i]*dt;
		if (mywall.dEpsN == 0) { /* sticky wall! */
			for (i=0;i<3;i++) c1.v[i] = c1.w[i] = 0;
			pkd->pStore[c1.id.iIndex].bStuck = 1;
			*piOutcome = MERGE;
			}
		else { /* bounce */
			double m1,r1,i1,vWall,ndotr=0,rprime[3],rlength=0;
			double n[3],s1[3],s2[3],v2[3],u[3],udotn=0,un[3],ut[3],w[3];
			double dEpsN,dEpsT;
			m1 = c1.fMass;
			r1 = c1.fRadius;
			i1 = 0.4*m1*r1*r1;
			vWall = (mywall.hotParam)*CP->dSlideLimit*(1. - mywall.dEpsN)/(1. + mywall.dEpsN);

			for (i=0;i<3;i++) {
				ndotr += mywall.n[i] * c1.r[i];
				}
			for (i=0;i<3;i++) {
				rprime[i] = c1.r[i] - ndotr * mywall.n[i];
				rlength += rprime[i]*rprime[i];
				}
			rlength = sqrt(rlength);

			if (mywall.type == 1) {	/* cylindrical wall */
				for (i=0;i<3;i++) {
					n[i] = rprime[i] / rlength;
					}
				s2[0] = mywall.n[1]*n[2]-mywall.n[2]*n[1];
				s2[1] = mywall.n[2]*n[0]-mywall.n[0]*n[2];
				s2[2] = mywall.n[0]*n[1]-mywall.n[1]*n[0];
				for (i=0;i<3;i++) {
					s2[i] *= mywall.radius*mywall.omega;
					v2[i] = -vWall*n[i];
					}
				}
			else {				/* infinite flat wall */
				for (i=0;i<3;i++) {
					n[i] = (ndotr < mywall.ndotp) ? mywall.n[i] : -mywall.n[i];
					}

				s2[0] = mywall.n[1]*rprime[2]-mywall.n[2]*rprime[1];
				s2[1] = mywall.n[2]*rprime[0]-mywall.n[0]*rprime[2];
				s2[2] = mywall.n[0]*rprime[1]-mywall.n[1]*rprime[0];
				for (i=0;i<3;i++) {
					s2[i] *= mywall.omega;
					v2[i] = -vWall*n[i];
					}
				}

			s1[0] = r1 * (c1.w[1]*n[2] - c1.w[2]*n[1]);
			s1[1] = r1 * (c1.w[2]*n[0] - c1.w[0]*n[2]);
			s1[2] = r1 * (c1.w[0]*n[1] - c1.w[1]*n[0]);
			for (i=0;i<3;i++) {
				u[i] = v2[i] - c1.v[i] + s2[i] - s1[i];
				udotn += u[i] * n[i];
				}
			assert(udotn < 0.0);
			for (i=0;i<3;i++) {
				un[i] = udotn * n[i];
				ut[i] = u[i] - un[i];
				}

			if (c1.bTinyStep) {
				dEpsN = CP->dCollapseEpsN;
				dEpsT = CP->dCollapseEpsT;
				}
			else if (CP->iSlideOption == MaxTrv &&
					 fabs(c1.v[0]*n[0] + c1.v[1]*n[1] + c1.v[2]*n[2]) < CP->dSlideLimit) {
				dEpsN = CP->dSlideEpsN;
				dEpsT = CP->dSlideEpsT;
				}
			else {
				dEpsN = mywall.dEpsN;
				dEpsT = mywall.dEpsT;
				}

			w[0] = n[1]*u[2] - n[2]*u[1];
			w[1] = n[2]*u[0] - n[0]*u[2];
			w[2] = n[0]*u[1] - n[1]*u[0];
			for (i=0;i<3;i++) {
				c1.v[i] += (1 + dEpsN)*un[i] + (2.0/7)*(1 - dEpsT)*ut[i];
				c1.w[i] += (2.0/7)*(m1/i1)*(1 - dEpsT)*r1*w[i];
				}
			for (i=0;i<3;i++) c1.r[i] -= c1.v[i]*dt;
			*piOutcome = BOUNCE;
			}
		PutColliderInfo(pkd, &c1,c2.id.iOrder,&pkd->pStore[c1.id.iIndex],dt);
		cOut[0] = c1;
		cOut[1] = c2;
		*pnOut = 2;
		return;
		}
#else /* TUMBLER */
	if (c2.id.iOrder < 0) { /* wall or endpoint collision */
		const WALLS *w = &CP->walls;
		int wall,endpt;
		assert(c1.id.iPid == pkd->idSelf);
		assert(pkd->pStore[c1.id.iIndex].bStuck == 0);
		if (c2.id.iOrder < -w->nWalls) { /* endpoint */
			wall = (-c2.id.iOrder - w->nWalls - 1)/2;
			endpt = (-c2.id.iOrder - w->nWalls - 1)%2;
			}
		else { /* wall */
			wall = -c2.id.iOrder - 1;
			endpt = -1; /* dummy statement */
			}
		for (i=0;i<3;i++) c1.r[i] += c1.v[i]*dt;
		if (w->wall[wall].dEpsN == 0) { /* sticky wall! */
			for (i=0;i<3;i++) c1.v[i] = c1.w[i] = 0;
			pkd->pStore[c1.id.iIndex].bStuck = 1;
			*piOutcome = MERGE;
			}
		else { /* bounce */
			double m1,r1,i1,nx=0,nz=0,a,s1x,s1z,ux,uz,unx,unz,utx,utz;
			double dEpsN,dEpsT;
			m1 = c1.fMass;
			r1 = c1.fRadius;
			i1 = 0.4*m1*r1*r1;
			if (c2.id.iOrder < -w->nWalls) { /* endpoint */
				switch (endpt) {
				case 0:
					nx = w->wall[wall].x1 - c1.r[0];
					nz = w->wall[wall].z1 - c1.r[2];
					break;
				case 1:
					nx = w->wall[wall].x2 - c1.r[0];
					nz = w->wall[wall].z2 - c1.r[2];
					break;
				default:
					assert(0);
					}
				}
			else { /* wall */
				nx = w->wall[wall].z2 - w->wall[wall].z1;
				nz = w->wall[wall].x1 - w->wall[wall].x2;
				if ((w->wall[wall].x1 - c1.r[0])*nx +
					(w->wall[wall].z1 - c1.r[2])*nz < 0) {
					nx = -nx;
					nz = -nz;
					}
				}
			a = 1/sqrt(nx*nx + nz*nz);
			nx *= a;
			nz *= a;
			s1x =  r1*c1.w[1]*nz;
			s1z = -r1*c1.w[1]*nx;
			ux = -c1.v[0] - s1x;
			uz = -c1.v[2] - s1z;
			a = ux*nx + uz*nz;
			assert(a < 0.0);
			unx = a*nx;
			unz = a*nz;
			utx = ux - unx;
			utz = uz - unz;
			if (c1.bTinyStep) {
				dEpsN = CP->dCollapseEpsN;
				dEpsT = CP->dCollapseEpsT;
				}
			else if (CP->iSlideOption == MaxTrv &&
					 fabs(c1.v[0]*nx + c1.v[2]*nz) < CP->dSlideLimit) {
				dEpsN = CP->dSlideEpsN;
				dEpsT = CP->dSlideEpsT;
				}
			else {
				dEpsN = w->wall[wall].dEpsN;
				dEpsT = w->wall[wall].dEpsT;
				}
			c1.v[0] += (1 + dEpsN)*unx + (2.0/7)*(1 - dEpsT)*utx;
			c1.v[2] += (1 + dEpsN)*unz + (2.0/7)*(1 - dEpsT)*utz;
			c1.w[1] += (2.0/7)*(m1/i1)*(1 - dEpsT)*r1*(ux*nz - uz*nx);
			for (i=0;i<3;i++) c1.r[i] -= c1.v[i]*dt;
			*piOutcome = BOUNCE;
			}
		PutColliderInfo(pkd, &c1,c2.id.iOrder,&pkd->pStore[c1.id.iIndex],dt);
		cOut[0] = c1;
		cOut[1] = c2;
		*pnOut = 2;
		return;
		}
#endif /* !TUMBLER */
#endif /* SAND_PILE */

	/*
	 ** Determine whether second particle is supposed to be an image,
	 ** and adjust accordingly. This is assumed to be the case if the
	 ** particle positions are separated by more than half a box length.
	 ** This may fail for small box sizes or fast particles!
	 */

	if (bPeriodic) {
	  CalcOffset(pkd,&c1,&c2,fOffset
#ifdef SLIDING_PATCH
		     ,&fShear
#endif
		     );
		}

	/* Advance coordinates to impact time */

	for (i=0;i<3;i++) {
		c1.r[i] += c1.v[i]*dt;
		c2.r[i] += c2.v[i]*dt;
	}

#ifdef SLIDING_PATCH
	/* For sliding patch, need to transform spins to patch frame as well */

	ct = cos(pkd->PP->dOrbFreq*(pkd->dTime + dt));
	st = sin(pkd->PP->dOrbFreq*(pkd->dTime + dt));
	wx =  ct*c1.w[0] + st*c1.w[1];
	wy = -st*c1.w[0] + ct*c1.w[1];
	c1.w[0] = wx;
	c1.w[1] = wy;
	wx =  ct*c2.w[0] + st*c2.w[1];
	wy = -st*c2.w[0] + ct*c2.w[1];
	c2.w[0] = wx;
	c2.w[1] = wy;
#endif

	/* Determine collision outcome */
#ifdef RUBBLE_ZML	

	/* should this be merged with stuff below? (i.e. after "CP->iOutcomes & FRAG" line) */
	if (CP->iRubForcedOutcome == RUB_FORCED_NONE) {

  		extern void rubRubbleCollide(const COLLIDER *col_a,const COLLIDER *col_b,
									 const COLLISION_PARAMS *cp,COLLIDER **col_c,
									 double *dMassInDust,int *nOut);
		printf("c1.id.iOrder = %i c2.id.iOrder = %i\n",
			   c1.id.iOrder, c2.id.iOrder);
		rubRubbleCollide(&c1,&c2,CP,&c,&dMassInDust,&n); /*dMinMass to be part CP*/
      
		if (n == 1)
			iOutcome = MERGE; /* could have produced dust */
		else if (n == 2) {
			assert(0); /*DEBUG THIS SHOULD NEVER HAPPEN*/
			iOutcome = BOUNCE; /* ditto */
			/* need to call pkdBounce() here */
			}
		else if (n > 2)
			iOutcome = FRAG;
		else
			assert(0); /* not possible */
		
/*		printf("iOutcome = %s\n",iOutcome); *//*DEBUG 03.14.04*/
		/* return diagnostic info if necessary */
		if (bDiagInfo) {
			*piOutcome = iOutcome;
			*dT = 0.0; /* kinetic energy loss -- too nasty to worry about */
			for (i=0;i<(n < MAX_NUM_FRAG ? n : MAX_NUM_FRAG);i++) /* only first MAX_NUM_FRAG reported */
			  cOut[i] = c[i];
			*pnOut = n;
			}
		goto finish; /*need to call add or delete particle routine */
		}
#endif /* RUBBLE_ZML */
	
	v2 = 0;
  
	for (i=0;i<3;i++)
		v2 += (c2.v[i] - c1.v[i])*(c2.v[i] - c1.v[i]);

	ve2 = 2*(c1.fMass + c2.fMass)/(c1.fRadius + c2.fRadius);

	dImpactEnergy = 0.5*c1.fMass*c2.fMass/(c1.fMass + c2.fMass)*v2;

	iOutcome = MISS;

	
#ifdef RUBBLE_ZML  
	if (CP->iRubForcedOutcome == RUB_FORCED_MERGE) {
#else
	if (CP->iOutcomes == MERGE ||
		((CP->iOutcomes & MERGE) &&
		 v2 <= CP->dBounceLimit*CP->dBounceLimit*ve2)) {
#endif
		iOutcome = MERGE;
		pkdMerge(pkd,&c1,&c2,CP->dDensity,&c,&n,iTime0,dBaseStep,dTimeNow);
		assert(c != NULL);
		assert(n == 1);
		if (CP->iOutcomes & (BOUNCE | FRAG)) { /* check if spinning too fast */
			double w2max,w2=0;
			w2max = c->fMass/(c->fRadius*c->fRadius*c->fRadius);
			for (i=0;i<3;i++)
				w2 += c->w[i]*c->w[i];
			if (w2 > w2max) {
				int rv;
				free((void *) c);
				iOutcome = BOUNCE;
				rv = pkdBounce(pkd,&c1,&c2,CP->dEpsN,CP->dEpsT,CP->bFixCollapse,&c,&n);
				assert(rv == BOUNCE_OK);
				assert(c != NULL);
				assert(n == 2);
				}
			}
		}	
#ifdef RUBBLE_ZML
	else if (CP->iRubForcedOutcome == RUB_FORCED_BOUNCE) {
#else
    else if (CP->iOutcomes & BOUNCE) {
#endif
		double dEpsN=1,dEpsT=1;
		iOutcome = BOUNCE;
		if (c1.bTinyStep || c2.bTinyStep) {
			dEpsN = CP->dCollapseEpsN;
			dEpsT = CP->dCollapseEpsT;
			}
		else if ((CP->iSlideOption == EscVel && v2 < CP->dSlideLimit2*ve2) ||
				 (CP->iSlideOption == MaxTrv && v2 < CP->dSlideLimit2)) {
			/*#define EXPERIMENTAL*//*DEBUG*/
			dEpsN = CP->dSlideEpsN;
			dEpsT = CP->dSlideEpsT;
			}
		else if (v2 > CP->dCrushLimit) { /* add frag check */
			dEpsN = CP->dCrushEpsN;
			dEpsT = CP->dCrushEpsT;
			}
		else if (CP->iBounceOption == ConstEps) {
			dEpsN = CP->dEpsN;
			dEpsT = CP->dEpsT;
			}
		else {
			const double ls = 1.49597892e13;
			const double ts = 5.0226355648e6;
			double rn[3],vn=0,rr;
			for (i=0;i<3;i++) {
				rn[i] = c2.r[i] - c1.r[i];
				vn += (c2.v[i] - c1.v[i])*rn[i];
				}
			/* note vn >= 0 ==> near miss -- pkdBounce() will deal with it */
			if (vn < 0) {
				vn = - vn/sqrt(rn[0]*rn[0] + rn[1]*rn[1] + rn[2]*rn[2]);
				vn *= ls/ts; /* conversion to cm/s */
				switch (CP->iBounceOption) {
				case Frosty200: /* Hatzes et al 1984 */
					dEpsN = 0.32*pow(vn,-0.234);
					break;
				case Frosty120: /* Hatzes et al 1988 */
					dEpsN = 0.48*pow(vn,-0.20);
					break;
				case Compacted: /* ditto */
					rr = c1.fRadius*c2.fRadius/(c1.fRadius + c2.fRadius)*ls;
					dEpsN = 0.90*exp(-(-0.01*rr + 0.41)*vn);
					break;
				case Glancing: /* Supulver et al 1995 */
					dEpsN = 0.52*pow(vn,-0.14);
					break;
				default:
					assert(0);
					}
				if (dEpsN < 0.25) dEpsN = 0.25; /* could be a parameter */
				else if (dEpsN > 1) dEpsN = 1;
				dEpsT = CP->dEpsT;
				}
			}
		if (pkdBounce(pkd,&c1,&c2,dEpsN,dEpsT,CP->bFixCollapse,&c,&n) == NEAR_MISS)
			iOutcome = MISS;
		assert(c != NULL);
		assert(n == 2);
		}
	else if (CP->iOutcomes & FRAG) {
		iOutcome = FRAG;
		pkdFrag(pkd,&c1,&c2,&c,&n);
		assert(c != NULL);
		assert(n <= MAX_NUM_FRAG);
		}

	if (!CP->bFixCollapse)
		assert(iOutcome != MISS); /* SOMETHING has to happen... */

	if (bDiagInfo) *piOutcome = iOutcome;

	assert(n > 0);

	if (bDiagInfo) {
		if (dT) {
			double moi;
			*dT = 0; /* redundant */
			for (j=0;j<n;j++) {
				moi = 0.4*c[j].fMass*c[j].fRadius*c[j].fRadius;
				for (i=0;i<3;i++)
					*dT += c[j].fMass*(c[j].v[i]*c[j].v[i]) +
						moi*(c[j].w[i]*c[j].w[i]);
				}
			moi = 0.4*c1.fMass*c1.fRadius*c1.fRadius;
			for (i=0;i<3;i++)
				*dT -= c1.fMass*(c1.v[i]*c1.v[i]) + moi*(c1.w[i]*c1.w[i]);
			moi = 0.4*c2.fMass*c2.fRadius*c2.fRadius;
			for (i=0;i<3;i++)
				*dT -= c2.fMass*(c2.v[i]*c2.v[i]) + moi*(c2.w[i]*c2.w[i]);
			*dT *= 0.5;
			}
		/*DEBUG (void) printf("Compare: dT = %e\n",*dT); */
		for (i=0;i<n;i++)
			cOut[i] = c[i];
		*pnOut = n;
		}

#ifdef RUBBLE_ZML
 finish:
#endif
  
	/* Trace particles back to start of step */

	for (i=0;i<n;i++)
		for (j=0;j<3;j++)
			c[i].r[j] -= c[i].v[j]*dt;

#ifdef SLIDING_PATCH
	/* And transform spin(s) back to space frame */

	for (i=0;i<n;i++) {
	  wx = ct*c[i].w[0] - st*c[i].w[1];
	  wy = st*c[i].w[0] + ct*c[i].w[1];
	  c[i].w[0] = wx;
	  c[i].w[1] = wy;
	}
#endif

	/* Handle output cases */

	if (n == 1) /* merge */
	  MergerReorder(pkd,&c1,&c2,c,cOut,dt,CP,bDiagInfo,fOffset,
#ifdef RUBBLE_ZML
					dMassInDust,
#endif
#ifdef SLIDING_PATCH
			fShear,
#endif
			bPeriodic);
	else if (n == 2) { /* bounce or mass transfer */
		if (c1.id.iPid == pkd->idSelf)
			PutColliderInfo(pkd, &c[0],c2.id.iOrder,&pkd->pStore[c1.id.iIndex],dt);
		if (c2.id.iPid == pkd->idSelf) {
			if (bPeriodic) {
				for (i=0;i<3;i++)
					c[1].r[i] -= fOffset[i];
#ifdef SLIDING_PATCH
				c[1].v[1] -= fShear;
#endif
				}
			PutColliderInfo(pkd, &c[1],c1.id.iOrder,&pkd->pStore[c2.id.iIndex],dt);
			}
		}
	else { /* fragmentation */
#ifdef RUBBLE_ZML
		if (CP->iRubForcedOutcome == RUB_FORCED_NONE) {
			PARTICLE p;
			if (c1.id.iPid == pkd->idSelf) {
				PutColliderInfo(pkd, &c[0],c2.id.iOrder,&pkd->pStore[c1.id.iIndex],dt);
				pkd->pStore[c1.id.iIndex].iColor = CP->iRubColor; /* override color */
				}
			if (c2.id.iPid == pkd->idSelf) {
				PutColliderInfo(pkd, &c[1],c1.id.iOrder,&pkd->pStore[c2.id.iIndex],dt);
				pkd->pStore[c2.id.iIndex].iColor = CP->iRubColor;
				}
			/*
			 ** For now, create new particles on 1st processor. Hopefully
			 ** domain decomposition will rebalance quickly. Note that
			 ** rubble pile particles created by 2nd processor are ignored.
			 ** Why do we bother making them at all in that case? Good question.
			 */
			if (c1.id.iPid == pkd->idSelf) {
				p = pkd->pStore[c1.id.iIndex]; /* quick and dirty initialization */
				p.iColor = CP->iRubColor;
				assert(p.iColor != PLANETESIMAL);
				p.dtPrevCol = 0;
				p.iPrevCol = INT_MAX;
				for (i=2;i<n;i++) {
					PutColliderInfo(pkd, &c[i],-1,&p,dt);
					/* fix up some things PutColliderInfo() doesn't take care of */
					p.iOrgIdx = c[i].id.iOrgIdx;
					pkdNewParticle(pkd,p);
					}
				}
			}
#else
		assert(0); /* not implemented yet */
		/* note in this case new iOrders start at pkd->nMaxOrderDark */
#endif
		}

	/* Free resources */

	free((void *) c); /* (if c is NULL, op is ignored) */
	}

void pkdResetColliders(PKD pkd,int iOrder1,int iOrder2)
{
	/*
	 ** Set SMOOTHACTIVE flag for those particles that were involved
	 ** in the collision or that were about to collide with one of the
	 ** particles involved in the collision. Reset otherwise. This
	 ** means new collision circumstances will only be determined for
	 ** these particles, while the remaining (active) particles will
	 ** rely on previously detected potential collisions.
	 */

	PARTICLE *p;
	int i,nLocal = pkdLocal(pkd);

	for (i=0;i<nLocal;i++) {
		p = &pkd->pStore[i];
		if (!TYPEQueryACTIVE(p)) continue;
		if (p->iOrder == iOrder1 || p->iOrder == iOrder2 ||
			p->iOrderCol == iOrder1 || p->iOrderCol == iOrder2)
			TYPESet(p,TYPE_SMOOTHACTIVE);
		else
			TYPEReset(p,TYPE_SMOOTHACTIVE);
		}
	}


#ifdef SLIDING_PATCH
void pkdFindLargeMasses(PKD pkd,double dMass,double dCentralMass,double dHelio,double hill,PARTICLE *p,double *r,int *n)
{
    int i,nLocal = pkdLocal(pkd);

    *n=0;
    for (i=0;i<nLocal;i++) {
	if (pkd->pStore[i].fMass > dMass) {
	    p[*n]=pkd->pStore[i];
	    r[*n]=hill*(pow((pkd->pStore[i].fMass/(3*dCentralMass)),(1.0/3))*(dHelio+pkd->pStore[i].r[0])); /* This assumes z height is negligible */ 
	    *n+= 1;
	    assert(*n < MAXLARGEMASS);
	    }
	}	
    }

void pkdGetNeighborParticles(PKD pkd,double *r,double dRadius2,int id,double dTime,PARTICLE *p,double *sep2,int *n) 
{
    double dDist2,dist2;
    double dx,dx1,dy,dy1,dz,dz1,sx,sy,sz;
    int i,j,nLocal = pkdLocal(pkd);

    /* This subroutine assumes that it is part of the sliding
       patch. Although there may be future uses for this that would be
       more general. Because of the unusual geometry of the sliding
       patch (namely that the ghost cells themselves are shearing),
       this algorithm is surprisingly tricky, and is modelled after
       Derek's INTERSECT macro in smooth.h. --Rory 07/12/04 
    */

    *n=0;
    i=0;
    while (i < nLocal) {
	/* First make sure we're not dealing with the large mass
	   particle */
	if (pkd->pStore[i].iOrder == id) {
	    i++;
	    continue;
	    }
	
	dDist2=0;	
	dx = pkd->pStore[i].r[0] - r[0];
	dx1 = r[0] - pkd->pStore[i].r[0];
	sy = r[1];
	if (dx > 0.0 ) {
	    dx1 += pkd->fPeriod[0];
	    if (dx1 < dx) {
		dDist2 = dx1*dx1;
		sx = r[0] + pkd->fPeriod[0];
		sy += SHEAR(1,dTime,pkd->PP);
		if (sy >= 0.5*pkd->fPeriod[1]) sy -= pkd->fPeriod[1];
		else if (sy < -0.5*pkd->fPeriod[1]) sy += pkd->fPeriod[1];
		}
	    else {
		dDist2 = dx*dx;
		sx = r[0];
		}
	    if (dDist2 > dRadius2) {
		i++;
		continue;
		}	    
	    }
	
	else if (dx1 > 0.0) {
	    dx += pkd->fPeriod[0];
	    if (dx < dx1) {
		dDist2 = dx*dx;
		sx =  r[0] - pkd->fPeriod[0];
		sy += SHEAR(-1,dTime,pkd->PP);
		if (sy >= 0.5*pkd->fPeriod[1]) sy -= pkd->fPeriod[1];
		else if (sy < -0.5*pkd->fPeriod[1]) sy += pkd->fPeriod[1];
		}
	    else {
		dDist2 = dx1*dx1;
		sx = r[0];
		}
	    if (dDist2 > dRadius2) {
		i++;
		continue;
		}	
	    }
	else {
	    dDist2 = 0.0;
	    sx = r[0];
	    }
	
	dy = pkd->pStore[i].r[1] - sy;
	dy1 = sy - pkd->pStore[i].r[1];
	if (dy > 0.0) {
	    dy1 += pkd->fPeriod[1];
	    if (dy1 < dy) {
		dDist2 += dy1*dy1;
		sy += pkd->fPeriod[1];
		}
	    else {
		dDist2 += dy*dy;
		}
	    if (dDist2 > dRadius2) {
		i++;
		continue;
		}	    
	    }
	else if (dy1 > 0) {
	    dy += pkd->fPeriod[1];
	    if (dy < dy1) {
		dDist2 += dy*dy;
		sy -= pkd->fPeriod[1];
		}
	    else {
		dDist2 += dy1*dy1;
		}
	    if (dDist2 > dRadius2) {
		i++;
		continue;
		}
	    }
	else {
	    }
	
	/* Is this necessary in sliding_patch? */
	dz = pkd->pStore[i].r[2] - r[2];
	dz1 = r[2] - pkd->pStore[i].r[2];
	if (dz > 0.0) {
	    dz1 += pkd->fPeriod[2];
	    if (dz1 < dz) {
		dDist2 += dz1*dz1;
		sz = r[2]+pkd->fPeriod[2];
		}
	    else {
		dDist2 += dz*dz;
		sz = r[2];
		}
	    if (dDist2 > dRadius2) {
		i++;
		continue;
		}
	    }
	else if (dz1 > 0.0) {
	    dz += pkd->fPeriod[2];
	    if (dz < dz1) {
		dDist2 += dz*dz;
		sz = r[2]-pkd->fPeriod[2];
		}
	    else {
		dDist2 += dz1*dz1;
		sz = r[2];
		}
	    if (dDist2 > dRadius2) {
		i++;
		continue;
		}	    
	    }
	else {
	    sz = r[2];
	    }
    
	/* The particle is within the sphere */
	assert(dDist2 < dRadius2);
    
	/* Just for safe keeping here is the version in the
	   non-sliding patch, non periodic case */
	dist2=0;	
	for (j=0;j<3;j++) 
	    dist2+= (pkd->pStore[i].r[j]-r[j])*(pkd->pStore[i].r[j]-r[j]);
	   /* if ((dist2 < radius*radius) && (dist2 != 0)) { */
	
	p[*n] = pkd->pStore[i];
	sep2[*n] = dDist2;	
	*n+= 1;
	assert(*n < MAXNEIGHBORS);
	i++;
	}    
    }


#endif /* SLIDING_PATCH */


#endif /* COLLISIONS */
