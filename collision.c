#ifdef PLANETS

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <assert.h>

#include "pkd.h"
#include "collision.h"

void partToCollider(PARTICLE *p,int iPid,int iIndex,COLLIDER *c)
{
	int i;

	c->id.iPid = iPid;
	c->id.iIndex = iIndex;
	c->id.iOrder = p->iOrder;
	c->fMass = p->fMass;
	c->fRadius = 2*p->fSoft;
	for (i=0;i<3;++i) {
		c->r[i] = p->r[i];
		c->v[i] = p->v[i];
		c->w[i] = p->w[i];
		}
	c->dt = p->dt;
	}

void colliderToPart(COLLIDER *c,PARTICLE *p)
{
	/* id & dt info ignored */

	int i;

	p->fMass = c->fMass;
	p->fSoft = 0.5*c->fRadius;
	for (i=0;i<3;++i) {
		p->r[i] = c->r[i];
		p->v[i] = c->v[i];
		p->w[i] = c->w[i];
		}
	}

void pkdDoCollision(PKD pkd,int iOutcomes,double dEpsN,double dEpsT,
					COLLIDER *Collider1,COLLIDER *Collider2,
					double *pdImpactEnergy,int *iOutcome,
					COLLIDER *pOut,int *pnOut)
{
	/* Collider2 may or may not be local */

	COLLIDER *p1=&pkd->Collider1,*p2=&pkd->Collider2,*p;
	double dt=pkd->dImpactTime,v2,dImpactEnergy;
	int bReturn = (pkd->idSelf == p1->id.iPid);
	int i,j,n;

#ifdef RUBBLE_TEST
	if (p2->id.iPid == -1) { /* wall collision */
		assert(pkd->pStore[p1->id.iIndex].bStuck == 0);
		for (i=0;i<3;++i) {
			p1->r[i] += p1->v[i]*dt;
			p1->v[i] = 0;
			p1->w[i] = 0;
			}
		colliderToPart(p1,&pkd->pStore[p1->id.iIndex]);
		pkd->pStore[p1->id.iIndex].bStuck = 1;
		return;
		}
#endif /* RUBBLE_TEST */

	assert(p1->id.iPid == pkd->idSelf || p2->id.iPid == pkd->idSelf);

	/* Advance coordinates to impact time */

	for (i=0;i<3;++i) {
		p1->r[i] += p1->v[i]*dt;
		p2->r[i] += p2->v[i]*dt;
	}

	/* Store advanced info for diagnostic output */

	if (bReturn) {
		*Collider1 = *p1;
		*Collider2 = *p2;
		}

	/* Determine collision outcome */

	/*DEBUG impact energy not currently used...*/

	v2 = 0;

	for (i=0;i<3;i++)
		v2 += (p2->v[i] - p1->v[i])*(p2->v[i] - p1->v[i]);

#ifdef RUBBLE_TEST
/*	dEpsN = 0.2 + 0.2*random()/RAND_MAX;/*DEBUG!!!*/
/*	printf("v2 = %e\n",v2);*/
	if (v2 < 0.1) {
/*		printf("   low v2\n");*/
		dEpsN = 1.0; /*DEBUG*/
		dEpsT = 1.0; /*DEBUG*/
		}
#endif /* RUBBLE_TEST */

	dImpactEnergy = 0.5*p1->fMass*p2->fMass/(p1->fMass + p2->fMass)*v2;

	if (bReturn) *pdImpactEnergy = dImpactEnergy;

	/*DEBUG add proper outcome criteria here...following is temporary... */

	if (iOutcomes & MERGE)
		pkdMerge(pkd,&p,&n);
	else if (iOutcomes & BOUNCE)
		pkdBounce(pkd,dEpsN,dEpsT,&p,&n);
	else if (iOutcomes & FRAG)
		pkdFrag(pkd,&p,&n);

	if (bReturn) *iOutcome = iOutcomes;

	assert(n > 0);

	if (bReturn) {
		assert(n <= MAX_NUM_FRAG);
		for (i=0;i<n;++i)
			pOut[i] = p[i];
		*pnOut = n;
		}

	/* Trace output particles back to start of step */

	for (i=0;i<n;++i)
		for (j=0;j<3;++j)
			p[i].r[j] -= p[i].v[j]*dt;

	/* Handle output cases */

	if (n == 1) { /* merge */
		int iMrg=-1,iDel=-1;
		if (p1->id.iPid == pkd->idSelf) { /* local particle */
			/*
			 * Keep this particle if it has larger mass, or, if both particles
			 * are the same mass, keep this particle if it has a lower
			 * processor number, or, if both particles are on same processor,
			 * keep particle with lower iOrder.
			 */
			if (p1->fMass > p2->fMass ||
				(p1->fMass == p2->fMass &&
				 (p1->id.iPid < p2->id.iPid ||
				  (p1->id.iPid == p2->id.iPid &&
				   p1->id.iOrder < p2->id.iOrder))))
				iMrg = p1->id.iIndex; /* new centre of mass particle */
			else
				iDel = p1->id.iIndex; /* this particle is removed */
			}
		if (p2->id.iPid == pkd->idSelf) { /* same thing for other particle */
			if (p2->fMass > p1->fMass ||
				(p2->fMass == p1->fMass &&
				 (p2->id.iPid < p1->id.iPid ||
				  (p2->id.iPid == p1->id.iPid &&
				   p2->id.iOrder < p1->id.iOrder))))
				iMrg = p2->id.iIndex;
			else
				iDel = p2->id.iIndex;
			}
		if (iMrg > -1) {
			colliderToPart(&p[0],&pkd->pStore[iMrg]);
			if (bReturn) {
				pOut[0].id.iPid = pkd->idSelf;
				pOut[0].id.iIndex = iMrg;
				pOut[0].id.iOrder = pkd->pStore[iMrg].iOrder;
				}
			}
		if (iDel > -1) pkdDeleteParticle(pkd,iDel);
		}
	else if (n == 2) { /* bounce or mass transfer */
		if (p1->id.iPid == pkd->idSelf)
			colliderToPart(&p[0],&pkd->pStore[p1->id.iIndex]);
		if (p2->id.iPid == pkd->idSelf)
			colliderToPart(&p[1],&pkd->pStore[p2->id.iIndex]);
		}
	else { /* fragmentation */
		assert(0); /* not implemented yet */
		/* note in this case new iOrders start at pkd->nMaxOrderDark */
		}

	/* Free resources */

	free((void *)p);
	}

void pkdMerge(PKD pkd,COLLIDER **pOut,int *pnOut)
{
	/* note: output particle (pOut[0]) is new centre-of-mass particle */

	COLLIDER *p1=&pkd->Collider1,*p2=&pkd->Collider2,*p;
	FLOAT com_pos[3],com_vel[3],rc1[3],rc2[3],vc1[3],vc2[3],ang_mom[3];
	FLOAT m1,m2,m,r1,r2,r,i1,i2,i;
	int k;

	m1 = p1->fMass;
	m2 = p2->fMass;
	m = m1 + m2;
	r1 = p1->fRadius;
	r2 = p2->fRadius;
	r = r1*pow(m/m1,1.0/3); /*DEBUG assumes equal density!*/
	i1 = 0.4*m1*r1*r1;
	i2 = 0.4*m2*r2*r2;
	i = 0.4*m*r*r;

	for (k=0;k<3;++k) {
		com_pos[k] = (m1*p1->r[k] + m2*p2->r[k])/m;
		rc1[k] = p1->r[k] - com_pos[k];
		rc2[k] = p2->r[k] - com_pos[k];
		com_vel[k] = (m1*p1->v[k] + m2*p2->v[k])/m;
		vc1[k] = p1->v[k] - com_vel[k];
		vc2[k] = p2->v[k] - com_vel[k];
	}

	ang_mom[0] = m1*(rc1[1]*vc1[2] - rc1[2]*vc1[1]) + i1*p1->w[0] +
		         m2*(rc2[1]*vc2[2] - rc2[2]*vc2[1]) + i2*p2->w[0];
	ang_mom[1] = m1*(rc1[2]*vc1[0] - rc1[0]*vc1[2]) + i1*p1->w[1] +
				 m2*(rc2[2]*vc2[0] - rc2[0]*vc2[2]) + i2*p2->w[1];
	ang_mom[2] = m1*(rc1[0]*vc1[1] - rc1[1]*vc1[0]) + i1*p1->w[2] +
				 m2*(rc2[0]*vc2[1] - rc2[1]*vc2[0]) + i2*p2->w[2];

	*pnOut = 1;
	*pOut = (COLLIDER *)malloc(*pnOut*sizeof(COLLIDER));
	assert(*pOut);

	p = &(*pOut)[0];

	/* note: id info set in pkdDoCollision(), used only for log */

	p->fMass = m;
	p->fRadius = r;

	for (k=0;k<3;++k) {
		p->r[k] = com_pos[k];
		p->v[k] = com_vel[k];
		p->w[k] = ang_mom[k]/i;
		}
	}

void pkdBounce(PKD pkd,double dEpsN,double dEpsT,COLLIDER **pOut,int *pnOut)
{
	/* note: particle order preserved */

	COLLIDER *p1=&pkd->Collider1,*p2=&pkd->Collider2;
	FLOAT n[3],s1[3],s2[3],u[3],un[3],ut[3],p[3],q[3];
	FLOAT m1,m2,m,r1,r2,i1,i2,mu,alpha,beta;
	FLOAT a,b,c,d;
	int i;

	m1 = p1->fMass;
	m2 = p2->fMass;
	m = m1 + m2;
	r1 = p1->fRadius;
	r2 = p2->fRadius;
	i1 = 0.4*m1*r1*r1;
	i2 = 0.4*m2*r2*r2;
	mu = m1*m2/m;
	alpha = 2.5*(1/m1 + 1/m2);
	beta = 1/(1 + alpha*mu);

	a = 0;
	for (i=0;i<3;++i) {
		n[i] = (p2->r[i] - p1->r[i]);
		a += n[i]*n[i];
		}
	a = 1/sqrt(a);
	for (i=0;i<3;++i)
		n[i] *= a;

	s1[0] = r1*(p1->w[1]*n[2] - p1->w[2]*n[1]);
	s1[1] = r1*(p1->w[2]*n[0] - p1->w[0]*n[2]);
	s1[2] = r1*(p1->w[0]*n[1] - p1->w[1]*n[0]);

	s2[0] = r2*(p2->w[1]*n[2] - p2->w[2]*n[1]);
	s2[1] = r2*(p2->w[2]*n[0] - p2->w[0]*n[2]);
	s2[2] = r2*(p2->w[0]*n[1] - p2->w[1]*n[0]);

	for (i=0;i<3;++i)
		u[i] = (p2->v[i] - p1->v[i]) + (s2[i] - s1[i]);

	a = u[0]*n[0] + u[1]*n[1] + u[2]*n[2];
	for (i=0;i<3;++i) {
		un[i] = a*n[i];
		ut[i] = u[i] - un[i];
		}

#ifdef RUBBLE_TEST
	if (pkd->pStore[p1->id.iIndex].bStuck) {
		mu = m2;
		beta = 2.0/7;
		}
	else if (pkd->pStore[p2->id.iIndex].bStuck) {
		mu = m1;
		beta = 2.0/7;
		}
#endif /* RUBBLE_TEST */

	a = (1 + dEpsN);
	b = beta*(1 - dEpsT);
	for (i=0;i<3;++i)
		p[i] = a*un[i] + b*ut[i];

	a = mu*b;
	q[0] = a*(n[1]*u[2] - n[2]*u[1]);
	q[1] = a*(n[2]*u[0] - n[0]*u[2]);
	q[2] = a*(n[0]*u[1] - n[1]*u[0]);

#ifdef RUBBLE_TEST
	if (pkd->pStore[p1->id.iIndex].bStuck) {
		a = 0;
		b = -1;
		c = 0;
		d = r2/i2;
		}
	else if (pkd->pStore[p2->id.iIndex].bStuck) {
		a = 1;
		b = 0;
		c = r1/i1;
		d = 0;
		}
	else {
#endif /* RUBBLE_TEST */
	a =   m2/m;
	b = - m1/m;
	c = r1/i1;
	d = r2/i2;
#ifdef RUBBLE_TEST
	}
#endif /* RUBBLE_TEST */

	*pnOut = 2;
	*pOut = (COLLIDER *)malloc(*pnOut*sizeof(COLLIDER));
	assert(*pOut);

	(*pOut)[0] = *p1;
	(*pOut)[1] = *p2;

	p1 = &(*pOut)[0]; /* p1,p2 now point to output structures */
	p2 = &(*pOut)[1];

	for (i=0;i<3;++i) {
		p1->v[i] += a*p[i];
		p2->v[i] += b*p[i];
		p1->w[i] += c*q[i];
		p2->w[i] += d*q[i];
		}
	}

void pkdFrag(PKD pkd,COLLIDER **pOut,int *pnOut)
{
	/* note: be sure new particles have p->iActive = 1 */
	/* may need to assert(*pnOut <= MAX_NUM_FRAG) */
	/* remember to set id info for logging purposes */
	}

#endif /* PLANETS */
