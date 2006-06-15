/* aggs.c
 * 
 * PKDGRAV Source Code for Aggregate Handling
 * 
 * Author: Kenneth W. Flynn
 *         flynnk@astro.umd.edu
 * Mods:   Derek C. Richardson
 *         dcr@astro.umd.edu
 *
 * Modified: 01/28/01; DCR: 07/10/02, 5/29/03, 7/14/05
 */
 
#ifdef AGGS

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "aggs.h"
#include "collision.h" /* for COLLIDER struct */

void pkdAggsFind(PKD pkd,int *iMaxIdx)
{
	/*
	 ** Returns largest aggregate index found on local processor.
	 */

	PARTICLE *p;
	int i,n,iAggIdx;

	n = pkdLocal(pkd);

	for (i=0;i<n;i++) {
		p = &pkd->pStore[i];
		if (IS_AGG(p)) {
			iAggIdx = AGG_IDX(p);
			if (iAggIdx > *iMaxIdx)
				*iMaxIdx = iAggIdx;
			}
		}
	}

void pkdAggsConfirm(PKD pkd,int iAggIdx,int *bAssigned)
{
	/*
	 ** Sets flag if particle belonging to aggregate found on
	 ** local processor.
	 */

	PARTICLE *p;
	int i,n;

	n = pkdLocal(pkd);

	for (i=0;i<n;i++) {
		p = &pkd->pStore[i];
		if (AGG_IDX(p) == iAggIdx) {
			*bAssigned = 1;
			return;
			}
		}
	}

void pkdAggsMerge(PKD pkd,int iOldIdx,int iNewIdx)
{
	/*
	 ** "Merges" two aggregates by assigning aggregate index of each
	 ** particle in old aggregate to index of new aggregate.  A call
	 ** to msrAggsUpdate() is needed to compute new dynamical
	 ** quantities (COM, spin, etc.; this is done in msrAggsMerge()).
	 */

	/*DEBUG for now, the color of the old aggregate is preserved */

	PARTICLE *p;
	int i,n;

	n = pkdLocal(pkd);

	for (i=0;i<n;i++) {
		p = &pkd->pStore[i];
		if (AGG_IDX(p) == iOldIdx) {
			AGG_SET_IDX(p,iNewIdx);
			}
		}
	}

void pkdAggsBackDrift(PKD pkd,int iAggIdx,double dt)
{
	/*
	 ** Drifts aggregate particle space positions back interval dt
	 ** (should be to start of step) so that collision prediction
	 ** (which assumes particle positions are at start of step) will
	 ** work properly.  Also sets SMOOTHACTIVE for all particles in
	 ** aggregate to force recomputation of collision circumstances.
	 ** Particle "accelerations" are taken to be second-order terms
	 ** in velocity expression, as computed in pkdAggsSetSpaceVel().
	 */

	PARTICLE *p;
	int i,k,n;

	n = pkdLocal(pkd);

	for (i=0;i<n;i++) {
		p = &pkd->pStore[i];
		if (AGG_IDX(p) == iAggIdx) {
			for (k=0;k<3;k++)
				p->r[k] -= (p->v[k] + 0.5*p->a[k]*dt)*dt;
			TYPESet(p,TYPE_SMOOTHACTIVE);
			}
		}
	}

void pkdAggsGetCOM(PKD pkd,int iAggIdx,Scalar *m,Vector mr,Vector mv)
{
	/*
	 ** Computes contribution (moments) of local particles to center-
	 ** of-mass position and velocity of specified aggregate.
	 */

	PARTICLE *p;
	int i,k,n;

	n = pkdLocal(pkd);

	for (i=0;i<n;i++) {
		p = &pkd->pStore[i];
		if (AGG_IDX(p) == iAggIdx) {
			*m += p->fMass;
			for (k=0;k<3;k++) {
				mr[k] += p->fMass*p->r[k];
				mv[k] += p->fMass*p->v[k];
				}
			}
		}
	}

void pkdAggsGetAxesAndSpin(PKD pkd,int iAggIdx,const Vector r_com,
						   const Vector v_com,Matrix I,Vector L)
{
	/*
	 ** Computes contribution of local particles to inertia tensor
	 ** and angular momentum vector relative to center of mass of
	 ** specified aggregate.  Particles belonging to the aggregate
	 ** have a copy of their positions relative to the aggregate COM
	 ** stored in p->r_agg (in space coordinates, not body coordinates,
	 ** since the transformation matrix isn't availble yet---we're
	 ** building it now!; to transform later, see pkdAggsToBodyAxes()).
	 */

	PARTICLE *p;
	Vector r,v;
	Scalar m,R;
	double q;
	int i,n;

	n = pkdLocal(pkd);
 
	for (i=0;i<n;i++) {
		p = &pkd->pStore[i];
		if (AGG_IDX(p) == iAggIdx) {
			/* get pos & vel wrt COM */
			vectorSub(p->r,r_com,r);
			vectorSub(p->v,v_com,v);
			m = p->fMass;
			R = RADIUS(p); /* twice softening length */
			q = AGGS_PARTICLE_INERTIA_PREFACTOR*R*R; /* for convenience */
			/* add inertia tensor contributions */
			/* note caller must fill symmetric elements; Cf. msrAggsGetAxesAndSpin() */
			I[0][0] += m*(q + r[1]*r[1] + r[2]*r[2]); 
			I[0][1] -= m*r[0]*r[1];
			I[0][2] -= m*r[0]*r[2];
			I[1][1] += m*(q + r[0]*r[0] + r[2]*r[2]);
			I[1][2] -= m*r[1]*r[2];
			I[2][2] += m*(q + r[0]*r[0] + r[1]*r[1]);
			/* add angular momentum contributions, L = m (r x v) + I w */
			/* note w for aggregate particles should equal w of aggregate as a whole */
			L[0] += m*(r[1]*v[2] - r[2]*v[1] + p->w[0]*q);
			L[1] += m*(r[2]*v[0] - r[0]*v[2] + p->w[1]*q);
			L[2] += m*(r[0]*v[1] - r[1]*v[0] + p->w[2]*q);
			/* store pos wrt COM */
			vectorCopy(r,p->r_agg);
			}
		}
	}

void pkdAggsSetBodyPos(PKD pkd,int iAggIdx,Matrix spaceToBody)
{
	/*
	 ** Transforms positions of local particles (belonging
	 ** to specified aggregate) from space to body coordinates
	 ** relative to center of mass.  Note pkdAggsGetAxes() must
	 **	be called first.
	 */

	PARTICLE *p;
	Vector tmp;
	int i,n;

	n = pkdLocal(pkd);

	for (i=0;i<n;i++) {
		p = &pkd->pStore[i];
		if (AGG_IDX(p) == iAggIdx) {
			matrixTransform(spaceToBody,p->r_agg,tmp);
			vectorCopy(tmp,p->r_agg);
			}
		}
	}

void pkdAggsSetSpacePos(PKD pkd,int iAggIdx,const Vector r_com,Matrix lambda)
{
	/*
	 ** Transforms positions of local particles (belonging to
	 ** specified aggregate) from body to space coordinates.
	 ** Called by msrAggsAdvance() during drift step.
	 */

	PARTICLE *p;
	int i,n;

	n = pkdLocal(pkd);

	for (i=0;i<n;i++) {
		p = &pkd->pStore[i];
		if (AGG_IDX(p) == iAggIdx) {
			/* transform positions from body to space coords wrt COM */
			matrixTransform(lambda,p->r_agg,p->r);
			/* add center of mass component */
			vectorAdd(p->r,r_com,p->r);
			}
		}
	}

void pkdAggsSetSpaceVel(PKD pkd,int iAggIdx,const Vector v_com,
						const Vector omega,Matrix lambda)
{
	/*
	 ** Computes space velocities of local particles (belonging
	 ** to specified aggregate) to 2nd order.  Called after a
	 ** kick by msrAggsKick() and before back drifting by
	 ** msrAggsBackDrift().
	 */

	PARTICLE *p;
	Vector v,a;
	int i,n;

	n = pkdLocal(pkd);

	for (i=0;i<n;i++) {
		p = &pkd->pStore[i];
		if (AGG_IDX(p) == iAggIdx) {
			/* compute aggregate spin component of particle velocities */
			vectorCross(omega,p->r_agg,v);
			/* transform to space frame */
			matrixTransform(lambda,v,p->v);
			/* add center of mass component */
			vectorAdd(p->v,v_com,p->v);
			/* compute 2nd-order (centripetal) term */
			vectorCross(omega,v,a);
			/* convenient to store this in particle's "a" vector */
			matrixTransform(lambda,v,p->a);
			}
		}
	}

void pkdAggsSetSpaceSpins(PKD pkd,int iAggIdx,const Vector omega)
{
	/* sets particle spins to aggregate spin */

	PARTICLE *p;
	int i,n;

	n = pkdLocal(pkd);

	for (i=0;i<n;i++) {
		p = &pkd->pStore[i];
		if (AGG_IDX(p) == iAggIdx)
			vectorCopy(omega,p->w);
		}
	}

void pkdAggsDelete(PKD pkd,int iAggIdx,int *bFound)
{
	PARTICLE *p;
	int i,n;

	n = pkdLocal(pkd);

	for (i=0;i<n;i++) {
		p = &pkd->pStore[i];
		if (AGG_IDX(p) == iAggIdx) {
			assert(*bFound == 0);
			p->iOrgIdx = INT_MAX; /*DEBUG for now--same as above*/
			p->iColor = 3; /*ditto*/
			*bFound = 1; /* could return, but keep looping as sanity check */
			}
		}
	}

void pkdAggsGetAccel(PKD pkd,int iAggIdx,Scalar *m,Vector ma)
{
	/*
	 ** Computes contribution (moments) of local particles to center
	 ** of mass acceleration of specified aggregate.  Must be called
	 ** after computing interparticle gravitational accelerations.
	 */

	PARTICLE *p;
	int i,k,n;

	n = pkdLocal(pkd);

	for (i=0;i<n;i++) {
		p = &pkd->pStore[i];
		if (AGG_IDX(p) == iAggIdx) {
			*m += p->fMass; /* used as a check in msrAggsGetAccel() */
			for (k=0;k<3;k++)
				ma[k] += p->fMass*p->a[k];
			}
		}
	}

void pkdAggsCheckStress(PKD pkd,int iAggIdx,const Vector r_com,const Vector a_com,
						const Vector omega,FLOAT fTensileStrength,FLOAT fShearStrength,
						int *nLost,int *nLeft)
{
	/*DEBUG COMMENT
	 */

	PARTICLE *p;
	Vector r,r_hat,a,a_rad,a_tan,v,a_cen;
	double adotr;
	FLOAT fTensileStress,fShearStress;
	int i,n;

	n = pkdLocal(pkd);

	for (i=0;i<n;i++) {
		p = &pkd->pStore[i];
		if (AGG_IDX(p) == iAggIdx) {
			/* particle position relative to com */
			vectorSub(p->r,r_com,r);
			/* unit radial vector from com to particle */
			vectorCopy(r,r_hat);
			vectorNorm(r_hat);
			/* differential acceleration wrt com */
			vectorSub(p->a,a_com,a);
			/* construct radial and tangential components */
			adotr = vectorDot(a,r_hat);
			vectorScale(r_hat,adotr,a_rad);
			vectorSub(a,a_rad,a_tan); /* same as r-hat x a */
			/* now add centrifugal term to radial component */
			vectorCross(omega,r,v);
			vectorCross(omega,v,a_cen);
			vectorSub(a_rad,a_cen,v); /* minus: centripetal --> centrifugal */
			/* compute stress */
			fTensileStress = vectorDot(v,r_hat); /* net radial component */
			fShearStress = vectorMag(a_tan); /* tangential component */
			if (fTensileStress > fTensileStrength ||
				fShearStress > fShearStrength) {
				if (fTensileStress > fTensileStrength) {
					/*
					 ** For this instant, apply the excess radial
					 ** acceleration to the particle.  This ensures it
					 ** will separate without immediately recolliding.
					 */
					vectorSub(a,a_rad,a);
					vectorScale(r_hat,fTensileStress - fTensileStrength,v);
					vectorAdd(a,v,a);
					vectorAdd(a_com,a,p->a);
					}
				if (fShearStress > fShearStrength) {
					/*DEBUG do nothing for now!*/
					}
				p->iOrgIdx = INT_MAX; /*DEBUG for now*/
				p->iColor = 3; /*ditto*/
				++(*nLost);
/*DEBUG! release all particles*/
printf("STRESS = %g i=%i\n",fTensileStress,p->iOrder);
for (i=0;i<n;i++) {
	p = &pkd->pStore[i];
	p->iOrgIdx = 9999;
	p->iColor = 3;
	}
*nLost = n; *nLeft = 0; return;
				}
			else
				++(*nLeft);
			}
		}
	}

void pkdAggsGetTorque(PKD pkd,int iAggIdx,const Vector r_com,const Vector a_com,
					  Vector torque)
{
	/*
	 ** Computes contribution of local particles to torque of
	 ** specified aggregate in space coordinates relative to center of
	 ** mass.  Note that the COM position is passed, rather than the
	 ** rotation matrix, since it's simpler this way...
	 */

	PARTICLE *p;
	Vector da,dr,cross,dN;
	int i,n;

	n = pkdLocal(pkd);

	for (i=0;i<n;i++) {
		p = &pkd->pStore[i];
		if (AGG_IDX(p) == iAggIdx) {
			vectorSub(p->r,r_com,dr); /* or could xform body to space */
			vectorSub(p->a,a_com,da);
			vectorCross(dr,da,cross);
			vectorScale(cross,p->fMass,dN);
			/* N += m (a x r) */
			vectorAdd(torque,dN,torque);
			}
		}
	}

static void aggsSingleParticleToAgg(const COLLIDER *c,Aggregate *a)
{
	double inertia;

	a->bAssigned = 1;
	a->mass = c->fMass;
	vectorCopy(c->r,a->r_com);
	vectorCopy(c->v,a->v_com);
	vectorZero(a->a_com); /* unused */
	vectorZero(a->torque); /* unused */
	vectorCopy(c->w,a->omega);
	inertia = AGGS_PARTICLE_INERTIA_PREFACTOR*c->fMass*c->fRadius*c->fRadius;
	vectorSet(a->moments,inertia,inertia,inertia);
	matrixIdentity(a->lambda);
	a->dLastUpdate = 0.0;
	}

static void aggsAggToSingleParticle(const Aggregate *a,COLLIDER *c)
{
	/* id, fRadius, iColor, dt, iRung, bTinyStep, agg unchanged */

	c->fMass = a->mass; /* unchanged in aggsBounce() */
	vectorCopy(a->r_com,c->r); /* unchanged in aggsBounce() */
	vectorCopy(a->v_com,c->v);
	vectorCopy(a->omega,c->w);
	}

static void aggsBounce(const COLLIDER *pc1,const COLLIDER *pc2,
					   const COLLISION_PARAMS *CP,COLLIDER *pcOut[],
					   int *pnOut)
{
	/*
	 ** Collision between aggregates of spheres, from Richardson 1995:
	 ** 	dv1 = gamma (1 + en) (m2/M) un n;
	 ** 	dv2 = - (m1/m2) dv1;
	 ** 	dw1 = m1 I1inv (c1 cross dv1);
	 **		dw2 = - m1 I2inv (c2 cross dv1);
	 ** where M = m1 + m2, un = normal of total relative velocity at
	 ** impact point, and gamma depends on the reduced mass and elements
	 ** of I1inv and I2inv in the ntp basis.
	 **
	 ** NOTE: method assumes NO sliding friction (et = 1).
	 */

	Aggregate *agg1,*agg2;
	Matrix ntpT,ntp,mtmp1,mtmp2,Ibody,spaceToBody1,spaceToBody2,Ispace1,Ispace2,a,b;
	Vector v,n,a1,a2,b1,b2,c1,c2,w1,w2,s1,s2,s,u,t,p,dv1,dv2,dw1,dw2,vtmp;
	double en,m1,m2,M,mu,R1,R2,R,un,c1t,c1p,c2t,c2p,gamma;

	*pnOut = 2;
	*pcOut = (COLLIDER *) malloc((*pnOut)*sizeof(COLLIDER));
	assert(*pcOut != NULL);

	(*pcOut)[0] = *pc1; /* struct copy */
	(*pcOut)[1] = *pc2;

	/*
	 ** Handy pointers.
	 ** NOTE: point to *output* here because we want to change values.
	 */

	agg1 = &((*pcOut)[0].agg);
	agg2 = &((*pcOut)[1].agg);

	/*
	 ** Rather than handle single particles as special cases, we turn
	 ** them into single-particle aggregates for the purpose of solving
	 ** the restitution equations.  They're turned back into single
	 ** particles at the end.
	 */

	if (!COLLIDER_IS_AGG(pc1))
		aggsSingleParticleToAgg(pc1,agg1);
	if (!COLLIDER_IS_AGG(pc2))
		aggsSingleParticleToAgg(pc2,agg2);

	/* sanity check */

	assert(agg1->bAssigned);
	assert(agg2->bAssigned);

	/* convenient shorthand */

	en = CP->dEpsN; /*DEBUG no vel-dep coef/slide/collapse adjust for now*/
	m1 = agg1->mass;
	m2 = agg2->mass;
	M = m1 + m2;
	mu = m1*m2/M;
	R1 = pc1->fRadius;
	R2 = pc2->fRadius;
	R = R1 + R2;

	/* v: relative linear velocity of agg centers of mass */

	vectorSub(agg2->v_com,agg1->v_com,v);

	/*
	 ** n: vector perpendicular to tangent plane at impact site, pointing
	 ** from contact sphere of aggregate 1 to contact sphere of aggregate 2.
	 */

	vectorSub(pc2->r,pc1->r,n);

	/*
	 ** a1,a2: positions of colliding spheres relative to agg centers of mass.
	 ** Note agg rotation during drift will cause some error here...
	 */

	vectorSub(pc1->r,agg1->r_com,a1);
	vectorSub(pc2->r,agg2->r_com,a2);

	/* b1,b2: position of impact site relative to colliding sphere centers */

	vectorScale(n, R1/R,b1);
	vectorScale(n,-R2/R,b2);

	/* c1,c2: position vectors of impact site relative to agg centers of mass */

	vectorAdd(a1,b1,c1);
	vectorAdd(a2,b2,c2);

	/* w1,w2: angular velocities of aggs in space frame */

	matrixTransform(agg1->lambda,agg1->omega,w1);
	matrixTransform(agg2->lambda,agg2->omega,w2);

	/* s1,s2: spin velocity of agg at impact site */

	vectorCross(w1,c1,s1);
	vectorCross(w2,c2,s2);

	/* s: relative spin velocity at impact site */

	vectorSub(s2,s1,s);

	/* u: total relative velocity at impact site */

	vectorAdd(v,s,u);

	/* construct ntp basis */

	vectorGetBasis(n,t,p);

	/* un: normal component of u */

	un = vectorDot(u,n);

	/* c1t,c1p,c2t,c2p: transverse components of c1,c2 */

	c1t = vectorDot(c1,t);
	c1p = vectorDot(c1,p);
	c2t = vectorDot(c2,t);
	c2p = vectorDot(c2,p);

	/* get inverse inertia tensors wrt ntp basis */

	vectorCopy(n,ntpT[0]); /* ntpT: matrix whose rows are n, t, and p */
	vectorCopy(t,ntpT[1]);
	vectorCopy(p,ntpT[2]);

	matrixTranspose(ntpT,ntp); /* ntp: matrix whose columns are n, t, and p */

	matrixDiagonal(agg1->moments,Ibody); /* inertia tensor in body frame */
	matrixMultiply(agg1->lambda,Ibody,mtmp1);
	matrixTranspose(agg1->lambda,spaceToBody1);
	matrixMultiply(mtmp1,spaceToBody1,Ispace1); /* Ispace1: inertia tensor of agg 1 in space frame */
	matrixMultiply(ntpT,Ispace1,mtmp1);
	matrixMultiply(mtmp1,ntp,mtmp2); /* inertia tensor in ntp basis */
	matrixInverse(mtmp2,a); /* a: inverse of inertia tensor of agg 1 in ntp basis */

	matrixDiagonal(agg2->moments,Ibody); /* inertia tensor in body frame */
	matrixMultiply(agg2->lambda,Ibody,mtmp1);
	matrixTranspose(agg2->lambda,spaceToBody2);
	matrixMultiply(mtmp1,spaceToBody2,Ispace2); /* Ispace2: inertia tensor of agg 2 in space frame */
	matrixMultiply(ntpT,Ispace2,mtmp1);
	matrixMultiply(mtmp1,ntp,mtmp2); /* inertia tensor in ntp basis */
	matrixInverse(mtmp2,b); /* b: inverse of inertia tensor of agg 2 in ntp basis */

	/* useful factor */

	gamma = 1.0/(1.0 + mu*(a[1][1]*c1p*c1p - 2.0*a[1][2]*c1p*c1t + a[2][2]*c1t*c1t +
						   b[1][1]*c2p*c2p - 2.0*b[1][2]*c2p*c2t + b[2][2]*c2t*c2t));

	/*DEBUG verbose agg conservation check
	{
	Vector P1,P2,P,L1t,L1rb,L1r,L1,L2t,L2rb,L2r,L2,L;
	Matrix inertia;

	(void) printf("r1 = (%g,%g,%g) v1 = (%g,%g,%g) w1 = (%g,%g,%g)\n",
				  agg1->r_com[0],agg1->r_com[1],agg1->r_com[2],
				  agg1->v_com[0],agg1->v_com[1],agg1->v_com[2],
				  agg1->omega[0],agg1->omega[1],agg1->omega[2]);
	(void) printf("r2 = (%g,%g,%g) v2 = (%g,%g,%g) w2 = (%g,%g,%g)\n",
				  agg2->r_com[0],agg2->r_com[1],agg2->r_com[2],
				  agg2->v_com[0],agg2->v_com[1],agg2->v_com[2],
				  agg2->omega[0],agg2->omega[1],agg2->omega[2]);
	vectorScale(agg1->v_com,agg1->mass,P1);
	vectorScale(agg2->v_com,agg2->mass,P2);
	vectorAdd(P1,P2,P);
	(void) printf("lin mom before = %g %g %g\n",P[0],P[1],P[2]);
	vectorCross(agg1->r_com,agg1->v_com,L1t);
	vectorScale(L1t,agg1->mass,L1t);
	matrixDiagonal(agg1->moments,inertia);
	matrixTransform(inertia,agg1->omega,L1rb);
	matrixTransform(agg1->lambda,L1rb,L1r);
	(void) printf("ang mom before L1t = %g %g %g\n",L1t[0],L1t[1],L1t[2]);
	(void) printf("ang mom before L1r = %g %g %g\n",L1r[0],L1r[1],L1r[2]);
	vectorAdd(L1t,L1r,L1);
	vectorCross(agg2->r_com,agg2->v_com,L2t);
	vectorScale(L2t,agg2->mass,L2t);
	matrixDiagonal(agg2->moments,inertia);
	matrixTransform(inertia,agg2->omega,L2rb);
	matrixTransform(agg2->lambda,L2rb,L2r);
	(void) printf("ang mom before L2t = %g %g %g\n",L2t[0],L2t[1],L2t[2]);
	(void) printf("ang mom before L2r = %g %g %g\n",L2r[0],L2r[1],L2r[2]);
	vectorAdd(L2t,L2r,L2);
	vectorAdd(L1,L2,L);
	(void) printf("ang mom before = %g %g %g\n",L[0],L[1],L[2]);
	}
	*/

	/* compute final velocities and spins */

	vectorScale(n,gamma*(1 + en)*(m2/M)*un,dv1);
	vectorAdd(agg1->v_com,dv1,agg1->v_com);

	vectorScale(dv1,-m1/m2,dv2);
	vectorAdd(agg2->v_com,dv2,agg2->v_com);

	vectorCross(c1,dv1,vtmp);
	matrixInverse(Ispace1,a); /* no longer in ntp basis */
	matrixTransform(a,vtmp,dw1);
	vectorScale(dw1,m1,dw1);
	vectorAdd(w1,dw1,w1);
	matrixTransform(spaceToBody1,w1,agg1->omega); /* new spin in body frame */

	vectorCross(c2,dv1,vtmp);
	matrixInverse(Ispace2,b); /* no longer in ntp basis */
	matrixTransform(b,vtmp,dw2);
	vectorScale(dw2,-m1,dw2);
	vectorAdd(w2,dw2,w2);
	matrixTransform(spaceToBody2,w2,agg2->omega); /* new spin in body frame */

	/*DEBUG verbose agg conservation check
	{
	Vector P1,P2,P,L1t,L1rb,L1r,L1,L2t,L2rb,L2r,L2,L;
	Matrix inertia;

	(void) printf("r1 = (%g,%g,%g) v1 = (%g,%g,%g) w1 = (%g,%g,%g)\n",
				  agg1->r_com[0],agg1->r_com[1],agg1->r_com[2],
				  agg1->v_com[0],agg1->v_com[1],agg1->v_com[2],
				  agg1->omega[0],agg1->omega[1],agg1->omega[2]);
	(void) printf("r2 = (%g,%g,%g) v2 = (%g,%g,%g) w2 = (%g,%g,%g)\n",
				  agg2->r_com[0],agg2->r_com[1],agg2->r_com[2],
				  agg2->v_com[0],agg2->v_com[1],agg2->v_com[2],
				  agg2->omega[0],agg2->omega[1],agg2->omega[2]);
	vectorScale(agg1->v_com,agg1->mass,P1);
	vectorScale(agg2->v_com,agg2->mass,P2);
	vectorAdd(P1,P2,P);
	(void) printf("lin mom  after = %g %g %g\n",P[0],P[1],P[2]);
	vectorCross(agg1->r_com,agg1->v_com,L1t);
	vectorScale(L1t,agg1->mass,L1t);
	matrixDiagonal(agg1->moments,inertia);
	matrixTransform(inertia,agg1->omega,L1rb);
	matrixTransform(agg1->lambda,L1rb,L1r);
	(void) printf("ang mom  after L1t = %g %g %g\n",L1t[0],L1t[1],L1t[2]);
	(void) printf("ang mom  after L1r = %g %g %g\n",L1r[0],L1r[1],L1r[2]);
	vectorAdd(L1t,L1r,L1);
	vectorCross(agg2->r_com,agg2->v_com,L2t);
	vectorScale(L2t,agg2->mass,L2t);
	matrixDiagonal(agg2->moments,inertia);
	matrixTransform(inertia,agg2->omega,L2rb);
	matrixTransform(agg2->lambda,L2rb,L2r);
	(void) printf("ang mom  after L2t = %g %g %g\n",L2t[0],L2t[1],L2t[2]);
	(void) printf("ang mom  after L2r = %g %g %g\n",L2r[0],L2r[1],L2r[2]);
	vectorAdd(L2t,L2r,L2);
	vectorAdd(L1,L2,L);
	(void) printf("ang mom  after = %g %g %g\n",L[0],L[1],L[2]);
	}
	*/

	/* revert back to single particles as needed */

	if (!COLLIDER_IS_AGG(pc1))
		aggsAggToSingleParticle(agg1,&((*pcOut)[0]));
	if (!COLLIDER_IS_AGG(pc2))
		aggsAggToSingleParticle(agg2,&((*pcOut)[1]));
	}

static void aggsPutColliderInfo(const COLLIDER *c,PARTICLE *p,double dt,
								int iAggIdx)
{
	/* used for merging particles with aggs in pkdAggsDoCollision() */

	int i;

	for (i=0;i<3;i++)
		p->r[i] = c->r[i]; /* position at contact */
	AGG_SET_IDX(p,iAggIdx); /* particle now belongs to this agg */
	p->iColor = 4 + iAggIdx%10; /*DEBUG a quick & dirty way to color aggs*/
	}

void pkdAggsDoCollision(PKD pkd,double dt,const COLLIDER *pc1,
						const COLLIDER *pc2,int bPeriodic,
						const COLLISION_PARAMS *CP,int iAggNewIdx,
						int *piOutcome,double *dT,
						COLLIDER *cOut,int *pnOut)
{
	COLLIDER c1,c2;
	double v2,ve2;
	int bReturnOutput,k;

	assert(bPeriodic == 0); /* for now */

	/*DEBUG verbose collision output*/

/*
	(void) printf("COLLISION %i (%i) & %i (%i) (dt = %.16e)\n",
				  pc1->id.iOrder,pc1->id.iOrgIdx,pc2->id.iOrder,pc2->id.iOrgIdx,dt);
*/

	/* get local copies of collider data for manipulation */

	c1 = *pc1; /* struct copy */
	c2 = *pc2;

	/*
	 ** To prevent overwriting data in parallel, only store results
	 ** in output variables if collider 1 is local to this processor.
	 */

	bReturnOutput = (c1.id.iPid == pkd->idSelf);

	if (bReturnOutput && dT != NULL) *dT = 0.0; /*DEBUG not used (change in energy not computed for aggs)*/

	/*
	 ** Advance coordinates of non-aggregate particles to impact time.
	 ** (Aggregate particles already advanced in msrAggsAdvance().)
	 ** Note that the aggregate advance step actually integrates the
	 ** Euler equations of motions to the impact time, taking into
	 ** account gravitational torques on each aggregate, whereas
	 ** collision prediction in CheckForCollision() uses a simpler
	 ** expression good to second order assuming the aggregate spin
	 ** vector(s) remain unchanged over the interval.  This means the
	 ** collision circumstances may be slightly off here (particles
	 ** either overlapping or not touching).  To minimize these
	 ** problems, the timestep should be SHORT.  Even so, it may be
	 ** necessary to turn on the "fix collapse" feature, particularly
	 ** if bouncing is allowed, to circumvent overlap errors.
	 */

	if (!COLLIDER_IS_AGG(&c1)) {
		for (k=0;k<3;k++)
			c1.r[k] += c1.v[k]*dt;
		}
	if (!COLLIDER_IS_AGG(&c2)) {
		for (k=0;k<3;k++)
			c2.r[k] += c2.v[k]*dt;
		}

	/* determine collision outcome */

	v2 = ve2 = 0.0;

	if ((CP->iOutcomes & MERGE) && (CP->iOutcomes & BOUNCE)) {

		/*
		 ** If both merging and bouncing are allowed, determine
		 ** outcome based on rough estimate of mutual escape
		 ** speed of colliders.  Aggregates are treated as giant
		 ** spherical particles for this purpose.
		 */

		Vector r1,r2;
		double m1,m2,d2=0.0;

		if (COLLIDER_IS_AGG(&c1)) {
			m1 = c1.agg.mass;
			vectorCopy(c1.agg.r_com,r1);
			}
		else {
			m1 = c1.fMass;
			vectorCopy(c1.r,r1);
			}

		if (COLLIDER_IS_AGG(&c2)) {
			m2 = c2.agg.mass;
			vectorCopy(c2.agg.r_com,r2);
			}
		else {
			m2 = c2.fMass;
			vectorCopy(c2.r,r2);
			}

		for (k=0;k<3;k++) {
			d2 += (r2[k] - r1[k])*(r2[k] - r1[k]);
			/*DEBUG following ignores 2nd-order terms in v for aggs*/
			v2 += (c2.v[k] - c1.v[k])*(c2.v[k] - c1.v[k]);
			}

		/*
		 ** Since aggregates can have arbitrarily bizarre shapes,
		 ** it's possible (though unlikely) for particles and/or
		 ** aggregates to have exactly overlapping mass centers.
		 ** So we impose a minimum separation for the escape speed
		 ** calculation of the sum of the touching particle radii.
		 */

		if (d2 < c1.fRadius*c1.fRadius + c2.fRadius*c2.fRadius)
			d2 = c1.fRadius*c1.fRadius + c2.fRadius*c2.fRadius;

		assert(d2 > 0.0);

		ve2 = 2*(m1 + m2)/sqrt(d2);
		}

	if (CP->iOutcomes == MERGE ||
		((CP->iOutcomes & MERGE) &&
		 v2 <= CP->dBounceLimit*CP->dBounceLimit*ve2)) {

		/*
		 ** Most of the work of merging aggregates is actually done
		 ** at the master level.  Here we're just concerned with
		 ** updating any unaggregated particles to reflect their new
		 ** aggregate member status.
		 */

		if (bReturnOutput) {
			*piOutcome = MERGE;
			*pnOut = 1;
			/* nothing stored in cOut -- master will take care of this */
			}

		if (COLLIDER_IS_AGG(&c1) && COLLIDER_IS_AGG(&c2)) {
			/* do nothing -- handled in msrAggsMerge() */
			assert(COLLIDER_AGG_IDX(&c1) != COLLIDER_AGG_IDX(&c2));
			}
		else if (COLLIDER_IS_AGG(&c1) && !COLLIDER_IS_AGG(&c2)) {
			/* add single particle at current position to aggregate */
			if (c2.id.iPid == pkd->idSelf) /* only if particle is local */
				aggsPutColliderInfo(&c2,&pkd->pStore[c2.id.iIndex],dt,
									COLLIDER_AGG_IDX(&c1));
			}
		else if (COLLIDER_IS_AGG(&c2) && !COLLIDER_IS_AGG(&c1)) {
			/* ditto */
			if (c1.id.iPid == pkd->idSelf)
				aggsPutColliderInfo(&c1,&pkd->pStore[c1.id.iIndex],dt,
									COLLIDER_AGG_IDX(&c2));
			}
		else { /* i.e., !COLLIDER_IS_AGG(&c1) && !COLLIDER_IS_AGG(&c2) */
			/* make new aggregate from single particles at current positions */
			if (c1.id.iPid == pkd->idSelf)
				aggsPutColliderInfo(&c1,&pkd->pStore[c1.id.iIndex],dt,
									iAggNewIdx);
			if (c2.id.iPid == pkd->idSelf)
				aggsPutColliderInfo(&c2,&pkd->pStore[c2.id.iIndex],dt,
									iAggNewIdx);
			}
		}
	else if (CP->iOutcomes & BOUNCE) {

		/* bounce */

		COLLIDER *c;
		int i,n;

		aggsBounce(&c1,&c2,CP,&c,&n);
		assert(n == 2);

		if (bReturnOutput) {
			*piOutcome = BOUNCE;
			for (i=0;i<n;i++)
				cOut[i] = c[i]; /* struct copy */
			*pnOut = n;
			}

		/*
		 ** Trace unaggregated particles back to start of step
		 ** (particles in aggregates updated in msrAggsBounce()).
		 ** The crazy (i+1)%2 business below is simply shorthand
		 ** for storing the iOrder of the *other* collider, since
		 ** for bouncing there can be only 2 particles involved.
		 */

		for (i=0;i<n;i++)
			if (!COLLIDER_IS_AGG(&c[i]) && c[i].id.iPid == pkd->idSelf) {
				for (k=0;k<3;k++)
					c[i].r[k] -= c[i].v[k]*dt;
				PutColliderInfo(&c[i],c[(i+1)%2].id.iOrder,
								&pkd->pStore[c[i].id.iIndex],dt);
				}

		/* free resources */

		free((void *) c);
		}
	else {
		assert(0);/*DEBUG no other outcomes allowed yet*/
		}

	/*
	 ** For aggs, need to set dtPrevCol and reset iPrevCol.  (For
	 ** non-aggregate particles, this is done in PutColliderInfo().)
	 ** Note we have to set iPrevCol to INT_MAX (infinity) here
	 ** because it's possible for particles inside different aggs to
	 ** collide with one another more than once during the interval.
	 */

	if (COLLIDER_IS_AGG(&c1) && c1.id.iPid == pkd->idSelf) {
		pkd->pStore[c1.id.iIndex].dtPrevCol = dt;
		pkd->pStore[c1.id.iIndex].iPrevCol = INT_MAX;
		}

	if (COLLIDER_IS_AGG(&c2) && c2.id.iPid == pkd->idSelf) {
		pkd->pStore[c2.id.iIndex].dtPrevCol = dt;
		pkd->pStore[c2.id.iIndex].iPrevCol = INT_MAX;
		}
	}

/*** Following routines called (or passed) directly from master ***/

void aggsEulerDerivs(FLOAT t,FLOAT vars[],void *agg_as_void,FLOAT derivs[])
{
	FLOAT *torque = ((Aggregate *) agg_as_void)->torque;
	FLOAT *moments = ((Aggregate *) agg_as_void)->moments;

	/* omega[0, 1, 2] */
	derivs[0] = (torque[0] + vars[1]*vars[2]*(moments[1] - moments[2]))/moments[0];
	derivs[1] = (torque[1] + vars[2]*vars[0]*(moments[2] - moments[0]))/moments[1];
	derivs[2] = (torque[2] + vars[0]*vars[1]*(moments[0] - moments[1]))/moments[2];
 
	/* q1[0, 1, 2] */
	derivs[3] = vars[2]*vars[6] - vars[1]*vars[9];
	derivs[4] = vars[2]*vars[7] - vars[1]*vars[10];
	derivs[5] = vars[2]*vars[8] - vars[1]*vars[11];
 
	/* q2 */
	derivs[6] = vars[0]*vars[9]  - vars[2]*vars[3];
	derivs[7] = vars[0]*vars[10] - vars[2]*vars[4];
	derivs[8] = vars[0]*vars[11] - vars[2]*vars[5];
 
	/* q3 */
	derivs[9]  = vars[1]*vars[3] - vars[0]*vars[6];
	derivs[10] = vars[1]*vars[4] - vars[0]*vars[7];
	derivs[11] = vars[1]*vars[5] - vars[0]*vars[8];
	}

void aggsRungeStep(FLOAT step_size,FLOAT x,FLOAT* y_vals,int n,
				   void* user_data,aggsRungeDerivs func,FLOAT* new_x,
				   FLOAT* new_y_vals)
{
	/*DEBUG could be reconciled with RungeStep() in runge.c*/

	const double one_sixth = 1.0/6.0,one_third = 1.0/3.0;
	FLOAT *k1,*k2,*k3,*k4,*tmp;
	int i;

	k1 = (FLOAT *) malloc(n*sizeof(FLOAT));
	assert(k1 != NULL);
	k2 = (FLOAT *) malloc(n*sizeof(FLOAT));
	assert(k2 != NULL);
	k3 = (FLOAT *) malloc(n*sizeof(FLOAT));
	assert(k3 != NULL);
	k4 = (FLOAT *) malloc(n*sizeof(FLOAT));
	assert(k4 != NULL);

	tmp = (FLOAT *) malloc(n*sizeof(FLOAT));
	assert(tmp != NULL);

	(*func)(x,y_vals,user_data,k1);
	for (i=0;i<n;i++) {
		k1[i] *= step_size;
		tmp[i] = y_vals[i] + 0.5*k1[i];
		}
	(*func)(x + 0.5*step_size,tmp,user_data,k2);
	for (i=0;i<n;i++) {
		k2[i] *= step_size;
		tmp[i] = y_vals[i] + 0.5*k2[i];
		}
	(*func)(x + 0.5*step_size,tmp,user_data,k3);
	for (i=0;i<n;i++) {
		k3[i] *= step_size;
		tmp[i] = y_vals[i] + k3[i];
		}
	(*func)(x + step_size,tmp,user_data,k4);
	for (i=0;i<n;i++) {
		k4[i] *= step_size;
		/* y(n+1) = y(n) + k1/6 + k2/3 + k3/3 + k4/6 */
		new_y_vals[i] = y_vals[i] + one_sixth*(k1[i] + k4[i]) +
			one_third*(k2[i] + k3[i]);
		}

	*new_x = x + step_size;

	free((void *) tmp);

	free((void *) k4);
	free((void *) k3);
	free((void *) k2);
	free((void *) k1);
	}

#endif /* AGGS */
