/* aggs.c
 * 
 * PKDGRAV Source Code
 * 
 * Author: Kenneth W. Flynn
 *         flynnk@astro.umd.edu
 * Mods:   Derek C. Richardson
 *
 * Modified: 01/28/01, DCR: 07/10/02
 */
 
#ifdef AGGS

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "linalg.h"  
#include "collision.h"
#include "aggs.h"
#include "pkd.h"

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
			iAggIdx = AGG_ID(p);
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
		if (AGG_ID(p) == iAggIdx) {
			*bAssigned = 1;
			return;
			}
		}
	}

void pkdAggsMerge(PKD pkd,int iOldID,int iNewID)
{
	/*
	 ** "Merges" two aggregates by assigning aggregate ID of each
	 ** particle in old aggregate to ID of new aggregate. Also sets
	 ** SMOOTHACTIVE for all particles in new aggregate to force
	 ** recomputation of collision circumstances (only needed for the
	 ** case when both colliders are aggregates because
	 ** pkdResetColliders() is sufficient when at least one of the
	 ** colliders is not an aggregate). A call to msrAggsUpdate() is
	 ** needed to compute new dynamical quantities (COM, spin, etc.).
	 */

	PARTICLE *p;
	int i,n;

	n = pkdLocal(pkd);

	for (i=0;i<n;i++) {
		p = &pkd->pStore[i];
		if (AGG_ID(p) == iOldID)
			AGG_SET_ID(p,iNewID);
		if (AGG_ID(p) == iNewID)
			TYPESet(p,TYPE_SMOOTHACTIVE);
		}
	}

void pkdAggsGetCOM(PKD pkd,int iAggIdx,FLOAT *m,FLOAT mr[3],FLOAT mv[3])
{
	/*
	 ** Computes contribution (moments) of local particles to center
	 ** of mass position and velocity of specified aggregate.
	 */

	PARTICLE *p;
	int i,k,n;

	n = pkdLocal(pkd);

	for (i=0;i<n;i++) {
		p = &pkd->pStore[i];
		if (AGG_ID(p) == iAggIdx) {
			*m += p->fMass;
			for (k=0;k<3;k++) {
				mr[k] += p->fMass*p->r[k];
				mv[k] += p->fMass*p->v[k];
				}
			}
		}
	}

void pkdAggsGetAxes(PKD pkd,int iAggIdx,FLOAT r_com[3],FLOAT v_com[3],
					FLOAT I[3][3],FLOAT L[3])
{
	/*
	 ** Computes contribution of local particles to inertia tensor and
	 ** angular momentum vector relative to center of mass of
	 ** specified aggregate.  Particles belonging to the aggregate
	 ** will have a copy of their positions and velocities relative to
	 ** the aggregate COM stored for later use by pkdAggsToBodyAxes().
	 */

	PARTICLE *p;
	Vector r,v;
	Scalar m,R;
	int i,n;

	n = pkdLocal(pkd);
 
	for (i=0;i<n;i++) {
		p = &pkd->pStore[i];
		if (AGG_ID(p) == iAggIdx) {
			/* get pos & vel wrt COM */
			vectorSub(p->r,r_com,r);
			vectorSub(p->v,v_com,v);
			m = p->fMass;
			R = 2*p->fSoft; /* particle radius is twice softening length */
			/* add inertia tensor contributions */
			I[0][0] += m*(GAMMA*R*R + r[1]*r[1] + r[2]*r[2]); 
			I[0][1] -= m*r[0]*r[1];
			I[0][2] -= m*r[0]*r[2];
			I[1][1] += m*(GAMMA*R*R + r[0]*r[0] + r[2]*r[2]);
			I[1][2] -= m*r[1]*r[2];
			I[2][2] += m*(GAMMA*R*R + r[0]*r[0] + r[1]*r[1]);
			/* add angular momentum contributions, L = m (r x v) */
			L[0] += m*(r[1]*v[2] - r[2]*v[1]);
			L[1] += m*(r[2]*v[0] - r[0]*v[2]);
			L[2] += m*(r[0]*v[1] - r[1]*v[0]);
			/* store pos wrt COM */
			vectorCopy(r,p->r_agg);
			}
		}
	}

void pkdAggsToBodyAxes(PKD pkd,int iAggIdx,Matrix spaceToBody)
{/*DEBUG make this more like pkdAggsToSpaceAxes()?*/
	/*
	 ** Transforms positions of local particles
	 ** (belonging to specified aggregate) from space to body
	 ** coordinates relative to center of mass.
	 **DEBUG to be called after GetAxes()
	 */

	PARTICLE *p;
	Vector tmp;
	int i,n;

	n = pkdLocal(pkd);

	for (i=0;i<n;i++) {
		p = &pkd->pStore[i];
		if (AGG_ID(p) == iAggIdx) {
			matrixTransform(spaceToBody,p->r_agg,tmp);
			vectorCopy(tmp,p->r_agg);
			}
		}
	}

void pkdAggsSetSpacePos(PKD pkd,int iAggIdx,Vector r_com,Matrix lambda)
{
	/*
	 ** Transforms positions of local particles (belonging to
	 ** specified aggregate) from body to space coordinates.
	 ** Called after a drift.
	 */

	PARTICLE *p;
	int i,n;

	n = pkdLocal(pkd);

	for (i=0;i<n;i++) {
		p = &pkd->pStore[i];
		if (AGG_ID(p) == iAggIdx) {
			/* transform positions from body to space coords wrt COM */
			matrixTransform(lambda,p->r_agg,p->r);
			/* add center of mass component */
			vectorAdd(p->r,r_com,p->r);
			}
		}
	}

void pkdAggsSetSpaceVel(PKD pkd,int iAggIdx,Vector v_com,
						Vector omega,Matrix lambda)
{
	/*
	 ** Transforms velocities of local particles (belonging to
	 ** specified aggregate) from body to space coordinates.
	 ** Called after a kick. Note depends on omega & lambda!
	 */

	PARTICLE *p;
	Vector v;
	int i,n;

	n = pkdLocal(pkd);

	for (i=0;i<n;i++) {
		p = &pkd->pStore[i];
		if (AGG_ID(p) == iAggIdx) {
			/* compute aggregate spin component of particle velocities */
			vectorCross(omega,p->r_agg,v);
			/* transform to space frame */
			matrixTransform(lambda,v,p->v);
			/* add center of mass component */
			vectorAdd(p->v,v_com,p->v);
			}
		}
	}

void pkdAggsAccel(PKD pkd,int iAggIdx,FLOAT *m,FLOAT ma[3])
{
	/*
	 ** Computes contribution (moments) of local particles to center
	 ** of mass acceleration of specified aggregate.
	 */

	PARTICLE *p;
	int i,k,n;

	n = pkdLocal(pkd);

	for (i=0;i<n;i++) {
		p = &pkd->pStore[i];
		if (AGG_ID(p) == iAggIdx) {
			*m += p->fMass; /*DEBUG store total mass in aggregate?*/
			for (k=0;k<3;k++) ma[k] += p->fMass*p->a[k];			
			}
		}
	}

void pkdAggsTorque(PKD pkd,int iAggIdx,FLOAT r_com[3],FLOAT a_com[3],
				   FLOAT torque[3])
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
		if (AGG_ID(p) == iAggIdx) {
			/* N += m (a x r) */
			vectorSub(p->r,r_com,dr); /* or could xform body to space */
			vectorSub(p->a,a_com,da);
			vectorCross(dr,da,cross);
			vectorScale(cross,p->fMass,dN);
			vectorAdd(torque,dN,torque);
			}
		}
	}

void pkdAggsActivate(PKD pkd)
{
	PARTICLE *p;
	int i,n;

	n = pkdLocal(pkd);

	for (i=0;i<n;i++) {
		p = &pkd->pStore[i];
		if (IS_AGG(p))
			TYPESet(p,TYPE_ALLACTIVE); /* turn on active types */
		}
	}

void pkdAggsDeactivate(PKD pkd)
{
	PARTICLE *p;
	int i,n;

	n = pkdLocal(pkd);

	for (i=0;i<n;i++) {
		p = &pkd->pStore[i];
		if (IS_AGG(p))
			TYPEClearACTIVE(p); /* turn off active types */
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
	assert(k1);
	k2 = (FLOAT *) malloc(n*sizeof(FLOAT));
	assert(k2);
	k3 = (FLOAT *) malloc(n*sizeof(FLOAT));
	assert(k3);
	k4 = (FLOAT *) malloc(n*sizeof(FLOAT));
	assert(k4);

	tmp = (FLOAT *) malloc(n*sizeof(FLOAT));
	assert(tmp);

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

void pkdAggsDoCollision(PKD pkd,double dt,const COLLIDER *pc1,
						const COLLIDER *pc2,int bPeriodic,
						const COLLISION_PARAMS *CP,
						int *piOutcome,double *dT,
						COLLIDER *cOut,int *pnOut)
{
	COLLIDER c1,c2;
	int k;

	assert(bPeriodic == 0); /* for now */

	/* currently only "merging" is allowed */

	c1 = *pc1;
	c2 = *pc2;

	if (!COLLIDER_IS_AGG(&c1))
		for (k=0;k<3;k++)
			c1.r[k] += c1.v[k]*dt;
	if (!COLLIDER_IS_AGG(&c2))
		for (k=0;k<3;k++)
			c2.r[k] += c2.v[k]*dt;

	*piOutcome = MERGE;
	*pnOut = 1;

	if (COLLIDER_IS_AGG(&c1) && COLLIDER_IS_AGG(&c2)) {
		/* do nothing -- handled in msrAggsCollisionUpdate() */
		assert(COLLIDER_AGG_ID(&c1) != COLLIDER_AGG_ID(&c2));
		}
	else if (COLLIDER_IS_AGG(&c1) && !COLLIDER_IS_AGG(&c2)) {
		AGG_SET_ID(&pkd->pStore[c2.id.iIndex],COLLIDER_AGG_ID(&c1));
		}
	else if (COLLIDER_IS_AGG(&c2) && !COLLIDER_IS_AGG(&c1)) {
		AGG_SET_ID(&pkd->pStore[c1.id.iIndex],COLLIDER_AGG_ID(&c2));
		}
	else { /* (!COLLIDER_IS_AGG(&c1) && !COLLIDER_IS_AGG(&c2)) */
		extern void PutColliderInfo(const COLLIDER *c,int
									iOrder2,PARTICLE *p,double dt);/*DEBUG*/
		AGG_SET_ID(&pkd->pStore[c1.id.iIndex],CP->iAggNewID);
		AGG_SET_ID(&pkd->pStore[c2.id.iIndex],CP->iAggNewID);
		/*DEBUG following calls awkward*/
		/*DEBUG!!! following won't work in parallel!!!*/
		PutColliderInfo(&c1,INT_MAX,&pkd->pStore[c1.id.iIndex],0.0);
		PutColliderInfo(&c2,INT_MAX,&pkd->pStore[c2.id.iIndex],0.0);
		}
	/* currently nothing stored in output */
	}

#endif /* AGGS */
