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
	c->fMass = p->fMass;
	c->fRadius = p->fSoft;
#ifdef SOFT_HACK
	c->fRadius *= 2;
#endif
	for (i=0;i<3;++i) {
		c->r[i] = p->r[i];
		c->v[i] = p->v[i];
		c->w[i] = p->w[i];
		}
	}

void colliderToPart(COLLIDER *c,PARTICLE *p)
{
	/* id info ignored */

	int i;

	p->fMass = c->fMass;
	p->fSoft = c->fRadius;
#ifdef SOFT_HACK
	p->fSoft *= 0.5;
#endif
	for (i=0;i<3;++i) {
		p->r[i] = c->r[i];
		p->v[i] = c->v[i];
		p->w[i] = c->w[i];
		}
	}

void pkdDoCollision(PKD pkd,COLLIDER *Collider1,COLLIDER *Collider2,
					double *pdImpactEnergy,COLLIDER *pOut,int *pnOut)
{
	/* Collider2 may or may not be local */

	COLLIDER *p1=&pkd->Collider1,*p2=&pkd->Collider2,*p;
	double dt=pkd->dImpactTime,v2,dImpactEnergy;
	int bReturn = (pkd->idSelf == p1->id.iPid);
	int i,j,n;

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

	v2 = 0;

	for (i=0;i<3;i++)
		v2 += (p2->v[i] - p1->v[i])*(p2->v[i] - p1->v[i]);

	dImpactEnergy = 0.5*(p1->fMass + p2->fMass)*v2;

	if (bReturn) *pdImpactEnergy = dImpactEnergy;

	if (dImpactEnergy < E_BOUNCE) 
		pkdMerge(pkd,&p,&n);
	else if (dImpactEnergy < E_FRAG)
		pkdBounce(pkd,&p,&n);
	else
		pkdFrag(pkd,&p,&n);

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
			 * Keep this particle if it has a lower processor number, or, if
			 * both particles are on same processor, keep particle with lower
			 * index number.
			 */
			if (p1->id.iPid < p2->id.iPid ||
				(p1->id.iPid == p2->id.iPid && p1->id.iIndex < p2->id.iIndex))
				iMrg = p1->id.iIndex; /* new centre of mass particle */
			else
				iDel = p1->id.iIndex; /* this particle is removed */
			}
		if (p2->id.iPid == pkd->idSelf) { /* same thing for other particle */
			if (p2->id.iPid < p1->id.iPid ||
				(p2->id.iPid == p1->id.iPid && p2->id.iIndex < p1->id.iIndex))
				iMrg = p2->id.iIndex;
			else
				iDel = p2->id.iIndex;
			}
		if (iMrg > -1) colliderToPart(&p[0],&pkd->pStore[iMrg]);
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
		}

	/* Free resources */

	free((void *) p);
	}

void pkdMerge(PKD pkd,COLLIDER **pOut,int *pnOut)
{
	/* note: output particle (pOut[0]) is new centre-of-mass particle */

	COLLIDER *p1=&pkd->Collider1,*p2=&pkd->Collider2,*p;
	FLOAT com_pos[3],com_vel[3],rc1[3],rc2[3],vc1[3],vc2[3],ang_mom[3];
	FLOAT m,r,i1,i2,i;
	int k;

	m = p1->fMass + p2->fMass;
	r = p1->fRadius*pow(m/p1->fMass,1.0/3); /*DEBUG assumes equal density!*/
	i1 = 0.4*p1->fMass*p1->fRadius*p1->fRadius;
	i2 = 0.4*p2->fMass*p2->fRadius*p2->fRadius;
	i = 0.4*m*r*r;

	for (k=0;k<3;++k) {
		com_pos[k] = (p1->fMass*p1->r[k] + p2->fMass*p2->r[k])/m;
		rc1[k] = p1->r[k] - com_pos[k];
		rc2[k] = p2->r[k] - com_pos[k];
		com_vel[k] = (p1->fMass*p1->v[k] + p2->fMass*p2->v[k])/m;
		vc1[k] = p1->v[k] - com_vel[k];
		vc2[k] = p2->v[k] - com_vel[k];
	}

	ang_mom[0] = p1->fMass * (rc1[1]*vc1[2] - rc1[2]*vc1[1]) + i1*p1->w[0] +
		         p2->fMass * (rc2[1]*vc2[2] - rc2[2]*vc2[1]) + i2*p2->w[0];
	ang_mom[1] = p1->fMass * (rc1[2]*vc1[0] - rc1[0]*vc1[2]) + i1*p1->w[1] +
				 p2->fMass * (rc2[2]*vc2[0] - rc2[0]*vc2[2]) + i2*p2->w[1];
	ang_mom[2] = p1->fMass * (rc1[0]*vc1[1] - rc1[1]*vc1[0]) + i1*p1->w[2] +
				 p2->fMass * (rc2[0]*vc2[1] - rc2[1]*vc2[0]) + i2*p2->w[2];

	*pnOut = 1;
	*pOut = (COLLIDER *) malloc(*pnOut * sizeof(COLLIDER));
	assert(*pOut);

	p = &(*pOut)[0];

	p->id = p1->id; /* currently redundant -- id info is ignored */

	p->fMass = m;
	p->fRadius = r;

	for (k=0;k<3;++k) {
		p->r[k] = com_pos[k];
		p->v[k] = com_vel[k];
		p->w[k] = ang_mom[k]/i;
		}
	}

void pkdBounce(PKD pkd,COLLIDER **pOut,int *pnOut)
{
	/* note: particle order preserved */

	COLLIDER *p1=&pkd->Collider1,*p2=&pkd->Collider2;
	FLOAT n[3],s1[3],s2[3],u[3],un[3],ut[3],p[3],q[3];
	FLOAT m,i1,i2,mu,alpha,beta;
	FLOAT a,b,c,d;
	int i;

	m = p1->fMass + p2->fMass;
	i1 = 0.4*p1->fMass*p1->fRadius*p1->fRadius;
	i2 = 0.4*p2->fMass*p2->fRadius*p2->fRadius;
	mu = p1->fMass*p2->fMass/m;
	alpha = 2.5*(1/p1->fMass + 1/p2->fMass);
	beta = 1/(1 + alpha*mu);

	a = 0;
	for (i=0;i<3;++i) {
		n[i] = (p2->r[i] - p1->r[i]);
		a += n[i]*n[i];
		}
	a = 1/sqrt(a);
	for (i=0;i<3;++i)
		n[i] *= a;

	s1[0] = p1->fRadius*(p1->w[1]*n[2] - p1->w[2]*n[1]);
	s1[1] = p1->fRadius*(p1->w[2]*n[0] - p1->w[0]*n[2]);
	s1[2] = p1->fRadius*(p1->w[0]*n[1] - p1->w[1]*n[0]);

	s2[0] = p2->fRadius*(p2->w[1]*n[2] - p2->w[2]*n[1]);
	s2[1] = p2->fRadius*(p2->w[2]*n[0] - p2->w[0]*n[2]);
	s2[2] = p2->fRadius*(p2->w[0]*n[1] - p2->w[1]*n[0]);

	for (i=0;i<3;++i)
		u[i] = (p2->v[i] - p1->v[i]) + (s2[i] - s1[i]);

	a = u[0]*n[0] + u[1]*n[1] + u[2]*n[2];
	for (i=0;i<3;++i) {
		un[i] = a*n[i];
		ut[i] = u[i] - un[i];
		}

	a = (1 + EN);
	b = beta*(1 - ET);
	for (i=0;i<3;++i)
		p[i] = a*un[i] + b*ut[i];

	a = mu*b;
	q[0] = a*(n[1]*u[2] - n[2]*u[1]);
	q[1] = a*(n[2]*u[0] - n[0]*u[2]);
	q[2] = a*(n[0]*u[1] - n[1]*u[0]);

	a =   p2->fMass/m;
	b = - p1->fMass/m;
	c = p1->fRadius/i1;
	d = p2->fRadius/i2;

	*pnOut = 2;
	*pOut = (COLLIDER *) malloc(*pnOut * sizeof(COLLIDER));
	assert(*pOut);

	p1 = &(*pOut)[0];
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
	}

#endif /* PLANETS */
