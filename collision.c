#ifdef COLLISIONS

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "pkd.h"
#include "collision.h"

#define BOUNCE_OK 0
#define NEAR_MISS 1

void
pkdNextCollision(PKD pkd,double *dt,int *iOrder1,int *iOrder2)
{
	/*
	 ** Returns time and iOrder of particles for earliest predicted
	 ** collision on this processor. Argument initialization must
	 ** occur in calling routine.
	 */

	PARTICLE *p;
	int i;

	for (i=0;i<pkdLocal(pkd);i++) {
		p = &pkd->pStore[i];
		if (p->iOrder < 0) continue; /* skip over deleted particles */
		if (p->dtCol < *dt) {
			*dt = p->dtCol;
			*iOrder1 = p->iOrder;
			*iOrder2 = p->iOrderCol;
			}
		}
	}

void
pkdGetColliderInfo(PKD pkd,int iOrder,COLLIDER *c)
{
	/*
	 ** Returns collider info for particle with matching iOrder.
	 */

	PARTICLE *p;
	int i,j;

	for (i=0;i<pkdLocal(pkd);i++) {
		p = &pkd->pStore[i];
		if (p->iOrder == iOrder) {
			c->id.iPid = pkd->idSelf;
			c->id.iOrder = iOrder;
			c->id.iIndex = i;
			c->fMass = p->fMass;
			c->fRadius = 2*p->fSoft;
			for (j=0;j<3;j++) {
				c->r[j] = p->r[j];
				c->v[j] = p->v[j];
				c->w[j] = p->w[j];
				}
			c->dt = p->dt;
			c->iRung = p->iRung;
			c->iColor = p->iColor;
			c->bTinyStep = p->bTinyStep;
			return;
			}
		}
	}

void
PutColliderInfo(const COLLIDER *c,int iOrder2,PARTICLE *p,double dt)
{
	/*
	 ** Stores collider info in particle structure (except id & color).
	 ** Also, dt is stored in dtPrevCol for inelastic collapse checks.
	 **
	 ** NOTE: Because colliding particles have their positions traced back to
	 ** the start of the step using their NEW velocities, it is possible for
	 ** a neighbour to no longer lie inside a collider's search ball. To fix
	 ** this, the ball radius is expanded by the square root (assuming random
	 ** walk) of a conservative estimate of the displacement amount (using
	 ** the "Manhattan metric"). This is probably only needed when the number
	 ** of particles is comparable to nSmooth, but we're going to do it anyway
	 ** just in case...
	 */

	int i;
	double r;

	p->fMass = c->fMass;
	p->fSoft = 0.5*c->fRadius;
	r = fabs(c->r[0] - p->r[0]) +
		fabs(c->r[1] - p->r[1]) +
		fabs(c->r[2] - p->r[2]);
	for (i=0;i<3;i++) {
		p->r[i] = c->r[i];
		p->v[i] = c->v[i];
		p->w[i] = c->w[i];
#ifdef NEED_VPRED
		p->vPred[i] = c->v[i] - dt*p->a[i];
#endif
		}
	p->iRung = c->iRung;
	p->fBall2 += 2*sqrt(p->fBall2)*r + r*r;
	p->dtPrevCol = dt;
	p->iPrevCol = iOrder2; /* stored to avoid false collisions */
	}

void
pkdMerge(PKD pkd,const COLLIDER *c1,const COLLIDER *c2,
		 COLLIDER **cOut,int *pnOut)
{
	/*
	 ** Merges colliders into new centre-of-mass particle (cOut[0]).
	 ** The radius of the merged particle is determined by the bulk
	 ** density of the largest of the colliders.
	 */

	COLLIDER *c;
	FLOAT com_pos[3],com_vel[3],rc1[3],rc2[3],vc1[3],vc2[3],ang_mom[3];
	FLOAT m1,m2,m,r1,r2,r,i1,i2,i;
	int k;

	m1 = c1->fMass;
	m2 = c2->fMass;
	m = m1 + m2;
	r1 = c1->fRadius;
	r2 = c2->fRadius;
	r = (m2 > m1 ? r2*pow(m/m2,1.0/3) : r1*pow(m/m1,1.0/3));
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
	}

	ang_mom[0] = m1*(rc1[1]*vc1[2] - rc1[2]*vc1[1]) + i1*c1->w[0] +
		         m2*(rc2[1]*vc2[2] - rc2[2]*vc2[1]) + i2*c2->w[0];
	ang_mom[1] = m1*(rc1[2]*vc1[0] - rc1[0]*vc1[2]) + i1*c1->w[1] +
				 m2*(rc2[2]*vc2[0] - rc2[0]*vc2[2]) + i2*c2->w[1];
	ang_mom[2] = m1*(rc1[0]*vc1[1] - rc1[1]*vc1[0]) + i1*c1->w[2] +
				 m2*(rc2[0]*vc2[1] - rc2[1]*vc2[0]) + i2*c2->w[2];

	*pnOut = 1;
	*cOut = (COLLIDER *) malloc(*pnOut*sizeof(COLLIDER));
	assert(*cOut);

	c = &(*cOut)[0];

	/* Note: id info set in pkdDoCollision(), used only for log */

	c->fMass = m;
	c->fRadius = r;

	for (k=0;k<3;k++) {
		c->r[k] = com_pos[k];
		c->v[k] = com_vel[k];
		c->w[k] = ang_mom[k]/i;
		}

	/* Set merger's timestep to iRung of largest mass. */
	/* XXX there is a bug in changing timesteps during a collision 
	   but this makes it less bothersome. */
	c->iRung = (c2->fMass > c1->fMass ? c2->iRung : c1->iRung);
	}

int
pkdBounce(PKD pkd,const COLLIDER *c1,const COLLIDER *c2,
		  double dEpsN,double dEpsT,COLLIDER **cOut,int *pnOut)
{
	/* Bounces colliders, preserving particle order */

	COLLIDER *co1,*co2;
	FLOAT n[3],s1[3],s2[3],v[3],s[3],u[3],un[3],ut[3],p[3],q[3];
	FLOAT m1,m2,m,r1,r2,i1,i2,mu,alpha,beta;
	FLOAT a,b,c,d;
	int i;

/*DEBUG
	(void) printf("e_n = %g e_t = %g\n",dEpsN,dEpsT);
*/

	*pnOut = 2;
	*cOut = (COLLIDER *) malloc(*pnOut*sizeof(COLLIDER));
	assert(*cOut);

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
		char ach[256];
		(void) sprintf(ach,"WARNING: %i & %i -- near miss?\n",
					   c1->id.iOrder,c2->id.iOrder);
#ifdef MDL_DIAG
		mdlDiag(pkd->mdl,ach);
#else
		(void) printf("%i: %s",pkd->idSelf,ach);
#endif
#ifdef FIX_COLLAPSE
		return NEAR_MISS; /* particles remain unchanged */
#else
		assert(0);
#endif
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

void
pkdFrag(PKD pkd,const COLLIDER *c1,const COLLIDER *c2,
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
	}

void
pkdDoCollision(PKD pkd,double dt,const COLLIDER *pc1,const COLLIDER *pc2,
			   int bPeriodic,const COLLISION_PARAMS *CP,int *piOutcome,
			   double *dT,COLLIDER *cOut,int *pnOut)
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

	COLLIDER c1,c2,*c;
	FLOAT fOffset[3] = {0,0,0};
	double v2,ve2,dImpactEnergy;
	int bDiagInfo,iOutcome,i,j,n;
#ifdef SLIDING_PATCH
	FLOAT fShear = 0;
#endif

/*DEBUG
	(void) printf("COLLISION %i & %i (dt = %.16e)\n",
				  pc1->id.iOrder,pc2->id.iOrder,dt);
*/

	/* Get local copies of collider info */

	c1 = *pc1;
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
			assert(udotn < 0);
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
		PutColliderInfo(&c1,c2.id.iOrder,&pkd->pStore[c1.id.iIndex],dt);
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
			double m1,r1,i1,nx,nz,a,s1x,s1z,ux,uz,unx,unz,utx,utz;
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
			assert(a < 0);
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
		PutColliderInfo(&c1,c2.id.iOrder,&pkd->pStore[c1.id.iIndex],dt);
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
		for (i=0;i<3;i++) {
			if (c2.r[i] - c1.r[i] >= 0.5*pkd->fPeriod[i]) {
				fOffset[i] -= pkd->fPeriod[i];
				c2.r[i] -= pkd->fPeriod[i];
				}
			else if (c2.r[i] - c1.r[i] < - 0.5*pkd->fPeriod[i]) {
				fOffset[i] += pkd->fPeriod[i];
				c2.r[i] += pkd->fPeriod[i];
				}
#ifdef SLIDING_PATCH
			if (i == 0) {
				if (fOffset[0] < 0) {
					fShear = 1.5*pkd->dOrbFreq*pkd->fPeriod[0];
					fOffset[1] = SHEAR(-1,pkd->fPeriod[0],pkd->fPeriod[1],
									   pkd->dOrbFreq,pkd->dTime);
					}
				else if (fOffset[0] > 0) {
					fShear = - 1.5*pkd->dOrbFreq*pkd->fPeriod[0];
					fOffset[1] = SHEAR(1,pkd->fPeriod[0],pkd->fPeriod[1],
									   pkd->dOrbFreq,pkd->dTime);
					}
				c2.r[1] += fOffset[1];
				c2.v[1] += fShear;
				}
#endif
			}
		}

	/* Advance coordinates to impact time */

	for (i=0;i<3;i++) {
		c1.r[i] += c1.v[i]*dt;
		c2.r[i] += c2.v[i]*dt;
	}

	/* Determine collision outcome */

	v2 = 0;

	for (i=0;i<3;i++)
		v2 += (c2.v[i] - c1.v[i])*(c2.v[i] - c1.v[i]);

	ve2 = 2*(c1.fMass + c2.fMass)/(c1.fRadius + c2.fRadius);

	dImpactEnergy = 0.5*c1.fMass*c2.fMass/(c1.fMass + c2.fMass)*v2;

	/* Determine collision outcome */

	iOutcome = MISS;

	if (CP->iOutcomes == MERGE || ((CP->iOutcomes & MERGE) && v2 <= ve2)) {
		iOutcome = MERGE;
		pkdMerge(pkd,&c1,&c2,&c,&n);
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
				rv = pkdBounce(pkd,&c1,&c2,CP->dEpsN,CP->dEpsT,&c,&n);
				assert(rv == BOUNCE_OK);
				assert(n == 2);
				}
			}
		}
	else if (CP->iOutcomes & BOUNCE) {
		double dEpsN=1,dEpsT=1;
		iOutcome = BOUNCE;
		if (c1.bTinyStep || c2.bTinyStep) {
			dEpsN = CP->dCollapseEpsN;
			dEpsT = CP->dCollapseEpsT;
			}
		else if ((CP->iSlideOption == EscVel && v2 < CP->dSlideLimit2*ve2) ||
				 (CP->iSlideOption == MaxTrv && v2 < CP->dSlideLimit2)) {
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
		if (pkdBounce(pkd,&c1,&c2,dEpsN,dEpsT,&c,&n) == NEAR_MISS)
			iOutcome = MISS;
		assert(n == 2);
		}
	else if (CP->iOutcomes & FRAG) {
		iOutcome = FRAG;
		pkdFrag(pkd,&c1,&c2,&c,&n);
		assert(n <= MAX_NUM_FRAG);
		}

#ifndef FIX_COLLAPSE
	assert(iOutcome != MISS); /* SOMETHING has to happen... */
#endif

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

	/* Trace particles back to start of step */

	for (i=0;i<n;i++)
		for (j=0;j<3;j++)
			c[i].r[j] -= c[i].v[j]*dt;

	/* Handle output cases */

	if (n == 1) { /* merge */
		int iMrg=-1,iDel=-1;
		if (c1.id.iPid == pkd->idSelf) { /* local particle */
			/*
			 ** Keep this particle if it has larger mass, or, if both
			 ** particles are the same mass, keep this particle if it
			 ** has a lower processor number, or, if both particles are
			 ** on same processor, keep particle with lower iOrder.
			 */
			if (c1.fMass > c2.fMass ||
				(c1.fMass == c2.fMass &&
				 (c1.id.iPid < c2.id.iPid ||
				  (c1.id.iPid == c2.id.iPid &&
				   c1.id.iOrder < c2.id.iOrder))))
				iMrg = c1.id.iIndex; /* new centre-of-mass particle */
			else
				iDel = c1.id.iIndex; /* this particle is removed */
			}
		if (c2.id.iPid == pkd->idSelf) { /* same thing for other particle */
			if (c2.fMass > c1.fMass ||
				(c2.fMass == c1.fMass &&
				 (c2.id.iPid < c1.id.iPid ||
				  (c2.id.iPid == c1.id.iPid &&
				   c2.id.iOrder < c1.id.iOrder))))
				iMrg = c2.id.iIndex;
			else
				iDel = c2.id.iIndex;
			}
		if (iMrg > -1) {
			PutColliderInfo(&c[0],INT_MAX,&pkd->pStore[iMrg],dt);
			if (bDiagInfo) {
				cOut[0].id.iPid = pkd->idSelf;
				cOut[0].id.iIndex = iMrg;
				cOut[0].id.iOrder = pkd->pStore[iMrg].iOrder;
				}
			}
		if (iDel > -1) pkdDeleteParticle(pkd,&pkd->pStore[iDel]);
		}
	else if (n == 2) { /* bounce or mass transfer */
		if (c1.id.iPid == pkd->idSelf)
			PutColliderInfo(&c[0],c2.id.iOrder,&pkd->pStore[c1.id.iIndex],dt);
		if (c2.id.iPid == pkd->idSelf) {
			if (bPeriodic) {
				for (i=0;i<3;i++)
					c[1].r[i] -= fOffset[i];
#ifdef SLIDING_PATCH
				c[1].v[1] -= fShear;
#endif
				}
			PutColliderInfo(&c[1],c1.id.iOrder,&pkd->pStore[c2.id.iIndex],dt);
			}
		}
	else { /* fragmentation */
		assert(0); /* not implemented yet */
		/* note in this case new iOrders start at pkd->nMaxOrderDark */
		}

	/* Free resources */

	free((void *) c);
	}

void
pkdResetColliders(PKD pkd,int iOrder1,int iOrder2)
{
	/*
	 ** Set SMOOTHACTIVE flag for those particles that were involved
	 ** in the collision or that were about to collide with one of the
	 ** particles involved in the collision. Reset otherwise.
	 */

	PARTICLE *p;
	int i;

	for (i=0;i<pkdLocal(pkd);i++) {
		p = &pkd->pStore[i];
		if (!TYPEQueryACTIVE(p)) continue;
		if (p->iOrder == iOrder1 || p->iOrder == iOrder2 ||
			p->iOrderCol == iOrder1 || p->iOrderCol == iOrder2)
			TYPESet(p,TYPE_SMOOTHACTIVE);
		else
			TYPEReset(p,TYPE_SMOOTHACTIVE);
		}
	}

#ifdef OLD_KEPLER/*DEBUG broken code*/
/*
 **	dInteract solves a thorny problem in our planetesimal code, namely
 ** when two elliptical orbits come close enough to be considered 
 ** interacting with each other in a non-perturbative sense. That is 
 ** roughly when the two particles come within h of each other, which we
 ** set as the Hill radius of the more massive particle. We name this 
 ** particle 'i' and put particle 'j' (the other one) into a coordinate
 ** frame defined by 'i'. The code will consider only "close approaches"
 ** occuring within the time step 'dDelta'. 'dTime' is the current time
 ** for both particles. If a "close approach" occurs this code will return
 ** a time slightly before the close approach takes place. Some very short
 ** impulsive encounters may be missed by this code, but they are not a 
 ** concern to the integration of the orbits. This means that it is not 
 ** intended to detect the precise moment of hard collisions without some
 ** more refined tests following the time returned by dInteract (and 
 ** suitable reduction of the reduced "Hill" radius).
 **
 ** Written by: Joachim Stadel, 1998 (University of Washington)
 **
 ** Currently assumes G := 1.
 */

double
dInteract(double dTime,double dDelta,double dCentMass,
		  PARTICLE *pi,PARTICLE *pj)
{
	double dEccAnom(double,double);

	/*
	 ** There are a lot of variable definitions for this code.
	 ** "Sorry, but it is just a teeny bit complicated you know."-JS
	 */
	double x,y,z,vx,vy,vz;
	double r2,v2,rv,r,ir,ia,ec,es;
	double hx,hy,hz,ih,ex,ey,ez,ie,yx,yy,yz;
	double imui,ei,ai,ini,Ti,Li;
	double imuj,ej,aj,inj,Tj,Lj;
	double h,lan,aop,i,E,Mend,Eend;
	double amax,amin,a,b,dE,Cx,Cy,x1,eta,y1;
	double l1,m1,n1,l2,m2,n2;
	double t,t1,t2,E0,dE0,E1,E2;
	double A,B,phi0,phi1,a0; /* for phase check */

	/*
	 ** Pick the most massive particle to define the hill sphere and the 
	 ** coordinate system in general.
	 */
	if (pj->fMass > pi->fMass) {
		PARTICLE *pt = pi;
		pi = pj;
		pj = pt;
		}
	/*
	 ** Calculate the constant in the keplerian potentials, that is 
	 ** what mu is. Currently I assume G := 1.
	 */
	imui = 1/(dCentMass + pi->fMass);
	imuj = 1/(dCentMass + pj->fMass);
	/*
	 ** Need the unit vectors defining the x,y,z coordinate direction
	 ** in particle i's orbit, where x is from focus to perihelion
	 ** and z is in the angular momentum direction, y defined then by
	 ** these two unit vectors.
	 */
	x = pi->r[0];
	y = pi->r[1];
	z = pi->r[2];
	vx = pi->v[0];
	vy = pi->v[1];
	vz = pi->v[2];
	r2 = x*x + y*y + z*z;
	v2 = vx*vx + vy*vy + vz*vz;
	rv = x*vx + y*vy + z*vz;
	r = sqrt(r2);
	ir = 1/r;
	hx = y*vz - z*vy;
	hy = z*vx - x*vz;
	hz = x*vy - y*vx;
	assert(hx || hy || hz); /* rectilinear orbits not supported */
	/*
	** Calculate the Runga-Lenz Vector.
	** The Runga-Lenz vector was actually first discovered by Laplace 1792!
	*/
	ex = imui*(vy*hz - vz*hy) - ir*x;
	ey = imui*(vz*hx - vx*hz) - ir*y;
	ez = imui*(vx*hy - vy*hx) - ir*z;
	ei = sqrt(ex*ex + ey*ey + ez*ez);
	/*
	 ** Create some unit vectors, defining the coordinate frame of
	 ** particle i.
	 */
	ih = 1/sqrt(hx*hx + hy*hy + hz*hz);
	hx *= ih;
	hy *= ih;
	hz *= ih;
	if (ei > 0) {
		ie = 1/ei;
		ex *= ie;
		ey *= ie;
		ez *= ie;
		}
	else {
		double t;
		if (hx != 0) {
			ey = 1;
			ex = -hy/hx;
			ez = 0;
			}
		else if (hy != 0) {
			ex = 1;
			ey = -hx/hy;
			ez = 0;
			}
		else if (hz != 0) {
			ey = 1;
			ez = -hy/hz;
			ex = 0;
			}
		else assert(0);
		t = 1/sqrt(ex*ex + ey*ey + ez*ez);
		ex *= t;
		ey *= t;
		ez *= t;
		}
	yx = hy*ez - hz*ey;
	yy = hz*ex - hx*ez;
	yz = hx*ey - hy*ex;
	/*
	 ** now get the other orbital quantities for particle i.
	 */
	ia = 2*ir - v2*imui;
	assert(ia > 0);
	ai = 1/ia;
	ini = sqrt(ai*ai*ai*imui);
	ec = 1 - r*ia;
	es = rv*ini*ia*ia;
	Ti = dTime - ini*(atan2(es,ec) - es);
	Li = ai*ei;
	/*
	 ** Now rotate particle j to i's coord system.
	 ** This is done most easily with the unit vectors we have already
	 ** calculated for particle i, in particular using the Runga-Lenz
	 ** vector, which we had to compute anyway. Also, we avoid trig in
	 ** this way.
	 */
	x = ex*pj->r[0] + ey*pj->r[1] + ez*pj->r[2];
	y = yx*pj->r[0] + yy*pj->r[1] + yz*pj->r[2];
	z = hx*pj->r[0] + hy*pj->r[1] + hz*pj->r[2];
	vx = ex*pj->v[0] + ey*pj->v[1] + ez*pj->v[2];
	vy = yx*pj->v[0] + yy*pj->v[1] + yz*pj->v[2];
	vz = hx*pj->v[0] + hy*pj->v[1] + hz*pj->v[2];
	/*
	 ** Need lan,aop and i for particle j now.
	 */
	r2 = x*x + y*y + z*z;
	v2 = vx*vx + vy*vy + vz*vz;
	rv = x*vx + y*vy + z*vz;
	r = sqrt(r2);
	ir = 1/r;
	hx = y*vz - z*vy;
	hy = z*vx - x*vz;
	hz = x*vy - y*vx;
	assert(hx || hy || hz);
	lan = atan2(hx,-hy);
	i = atan2(sqrt(hx*hx + hy*hy),hz);
	/*
	** Calculate the Runga-Lenz Vector.
	** The Runga-Lenz vector was actually first discovered by Laplace 1792!
	*/
	ex = imuj*(vy*hz - vz*hy) - ir*x;
	ey = imuj*(vz*hx - vx*hz) - ir*y;
	ez = imuj*(vx*hy - vy*hx) - ir*z;
	aop = atan2(ez,ex*sin(i)*cos(lan) + ey*sin(i)*sin(lan));
	/*
	 ** now get the other orbital quantities for particle j.
	 */
	ej = sqrt(ex*ex + ey*ey + ez*ez);
	ia = 2*ir - v2*imuj;
	assert(ia > 0);
	aj = 1/ia;
	inj = sqrt(aj*aj*aj*imuj);
	ec = 1 - r*ia;
	es = rv*inj*ia*ia;
	E = atan2(es,ec);
	if (E < 0) E += 2*M_PI;
	Tj = dTime - inj*(E - es);
	/*
	 ** Reject encounter if safely out of phase.
	 */
	A = (1/inj - 1/ini);
	assert(fabs(A*dDelta) < M_PI);
	B = (Ti/ini - Tj/inj);
	phi0 = A*dTime + B;
	phi1 = A*(dTime + dDelta) + B;
/*	phi0 -= 2*PI*floor(0.5*(1 + phi0*M_1_PI));
	phi1 -= 2*PI*floor(0.5*(1 + phi1*M_1_PI));*/
	phi0 -= 2*M_PI*floor(phi0/(2*M_PI));
	phi1 -= 2*M_PI*floor(phi1/(2*M_PI));
	if (phi0 >= M_PI) phi0 -= 2*M_PI;
	if (phi1 >= M_PI) phi1 -= 2*M_PI;
	assert(phi0 >= -M_PI && phi0 < M_PI);
	assert(phi1 >= -M_PI && phi1 < M_PI);
	a0 = ei + pi->fHill;
	if (A < 0) {
		double phi = phi0;
		phi0 = phi1;
		phi1 = phi;
		}
	phi0 -= ej;
	phi1 += ej;
	if (!((phi1 < a0 && phi1 > -a0) ||
		  (phi0 < a0 && phi0 > -a0) ||
		  (phi1 > a0 && phi0 < -a0))) return 0;
	/*
	 ** Calculate final M & E for particle j.
	 */
	Mend = (dTime + dDelta - Tj)/inj;
	Eend = dEccAnom(Mend,ej);			/* expensive? */
	/*
	 ** Make sure we get the wrapping correct.
	 */
	while (Eend < E) Eend += 2*M_PI;
	Lj = aj*ej;
	/*
	 ** Set h to be the hill radius of particle i.
	 ** Use the precomputed "reduced hill radius" stored with the particle.
	 */
	h = pi->fHill*ai;
	/*
	 ** This might be the right place to make a cut by considering the 
	 ** orbits of i and j in a uniformly rotating frame of the same 
	 ** period of particle i.
	 ** If no overlap of the maximum expected "synodic zones" then 
	 ** return with no interaction.
	 */

	/*
	 ** Set a safe angular window for the orbit of particle i, this is 
	 ** used later in determining the time window.
	 ** dE0 = h/bi;
	 */
	dE0 = h/sqrt(ai*ai - Li*Li);
	/*
	 ** Set up amin and amax which define the elliptical 
	 ** interaction annulus.
	 */
	amax = ai + h;
	amin = ai - h;
	assert(amin > 0);
	/*
	** Compute direction cosines for particle j's orbit in 
	** the coordinate system of particle i.
	*/
	l1 = cos(lan)*cos(aop) - sin(lan)*sin(aop)*cos(i);
	m1 = sin(lan)*cos(aop) + cos(lan)*sin(aop)*cos(i);
	n1 = sin(aop)*sin(i);
	l2 = -cos(lan)*sin(aop) - sin(lan)*cos(aop)*cos(i);
	m2 = -sin(lan)*sin(aop) + cos(lan)*cos(aop)*cos(i);
	n2 = cos(aop)*sin(i);
	/*
	 ** Set up coefficients for the differential ellipse generator for 
	 ** particle j.
	 ** Set dE to be a fraction of particle i's hill radius as an angle.
	 */
	dE = h/aj;
	Cy = sqrt(1 - ej*ej);
    x1 = aj*cos(E);
	y1 = aj*Cy*sin(E);
	Cx = 0.5*dE/Cy;
	Cy *= dE;
	while (E<Eend) {
		/*
		 ** 2nd order (like drift-kick-drift) ellipse generator.
		 */
		x1 -= Cx*y1;
		y1 += Cy*x1;
		x1 -= Cx*y1;
		E += dE;
		/*
		 ** Convert x1,y1 to a point in the coord system of 
		 ** particle i's orbit. First shift calculate eta 
		 ** which is wrt the focus of the ellispe instead of the
		 ** center, i.e., x1 shifted by Lj := a*e.
		 */
		eta = x1 - Lj;
		z = n1*eta + n2*y1;
		/*
		 ** Test if z is greater than h above or below the plane.
		 */
		if (z > h) continue;
		if (z < -h) continue;
		/*
		 ** Calculate the x and y coordinates now.
		 */
	    x = l1*eta + l2*y1;
		y = m1*eta + m2*y1;
		/*
		 ** Now calculate the u coordinate for this point!
		 ** We can actually do this test just as well with a itself.
		 */
		r2 = x*x + y*y;
		a = 0.5*(sqrt(r2) + sqrt(r2 + 4*Li*(x + Li)));
		/*
		 ** Now see if a lies within the limits, or if the point is in the 
		 ** annulus.
		 */
		if (a > amax) continue;
		if (a < amin) continue;
		/*
		 ** Okay, this point GEOMETRICALLY comes within h of particle i's
		 ** orbit. Now we still have to see if particle i is found within
		 ** h of this part of its orbit at this time.
		 */
		/*
		 ** max estimate of the anglular window, 
		 ** use bi intead of ai, see definition of dE0 above.
		 */
		t = inj*(E - ej*sin(E)) + Tj;
		b = sqrt(a*a - Li*Li);
		E0 = atan2(y/b,x/a);
		E1 = E0 - dE0;
		t1 = ini*(E1 - ei*sin(E1)) + Ti;
		E2 = E0 + dE0;
		t2 = ini*(E2 - ei*sin(E2)) + Ti;
		while (t1 < dTime+dDelta) {
			if (t>t1 && t<t2) {
				/*
				 ** we have an interaction at this time.
				 ** Go back to the last step to be conservative.
				 */
				E -= dE;
				t = inj*(E - ej*sin(E)) + Tj;
				return(t);
				}
			/*
			 ** Add period of i to t1 and t2 and test again.
			 */
			t1 += 2*M_PI*ini;
			t2 += 2*M_PI*ini;
			}
		}
	/*
	 ** No Interaction will take place.
	 */
	return(0);
	}
#endif /* OLD_KEPLER */

#endif /* COLLISIONS */
