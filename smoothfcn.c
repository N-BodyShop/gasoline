#include <math.h>
#include <assert.h>
#include "smoothfcn.h"

#ifdef COLLISIONS
#include "ssdefs.h"
#include "collision.h"
#endif /* COLLISIONS */

void initDensity(void *p)
{
	((PARTICLE *)p)->fDensity = 0.0;
	}

void combDensity(void *p1,void *p2)
{
	((PARTICLE *)p1)->fDensity += ((PARTICLE *)p2)->fDensity;
	}

void Density(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	FLOAT ih2,r2,rs,fDensity;
	int i;

	ih2 = 4.0/p->fBall2;
	fDensity = 0.0;
	for (i=0;i<nSmooth;++i) {
		r2 = nnList[i].fDist2*ih2;
		rs = 2.0 - sqrt(r2);
		if (r2 < 1.0) rs = (1.0 - 0.75*rs*r2);
		else rs = 0.25*rs*rs*rs;
		fDensity += rs*nnList[i].pPart->fMass;
		}
	p->fDensity = M_1_PI*sqrt(ih2)*ih2*fDensity; 
	}

void DensitySym(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	PARTICLE *q;
	FLOAT fNorm,ih2,r2,rs;
	int i;

	ih2 = 4.0/p->fBall2;
	fNorm = 0.5*M_1_PI*sqrt(ih2)*ih2;
	for (i=0;i<nSmooth;++i) {
		r2 = nnList[i].fDist2*ih2;
		rs = 2.0 - sqrt(r2);
		if (r2 < 1.0) rs = (1.0 - 0.75*rs*r2);
		else rs = 0.25*rs*rs*rs;
		rs *= fNorm;
		q = nnList[i].pPart;
		p->fDensity += rs*q->fMass;
		q->fDensity += rs*p->fMass;
		}
	}

#ifdef SUPERCOOL

void initMeanVel(void *p)
{
	int j;

	for (j=0;j<3;++j) {
		((PARTICLE *)p)->vMean[j] = 0.0;
		}
	}

void combMeanVel(void *p1,void *p2)
{
	int j;

	for (j=0;j<3;++j) {
		((PARTICLE *)p1)->vMean[j] += ((PARTICLE *)p2)->vMean[j];
		}
	}

void MeanVel(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	PARTICLE *q;
	FLOAT fNorm,ih2,r2,rs;
	int i,j;

	ih2 = 4.0/p->fBall2;
	fNorm = M_1_PI*sqrt(ih2)*ih2;
	for (i=0;i<nSmooth;++i) {
		r2 = nnList[i].fDist2*ih2;
		rs = 2.0 - sqrt(r2);
		if (r2 < 1.0) rs = (1.0 - 0.75*rs*r2);
		else rs = 0.25*rs*rs*rs;
		rs *= fNorm;
		q = nnList[i].pPart;
		for (j=0;j<3;++j) {
			p->vMean[j] += rs*q->fMass/q->fDensity*q->v[j];
			}
		}
	}

void MeanVelSym(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	PARTICLE *q;
	FLOAT fNorm,ih2,r2,rs;
	int i,j;

	ih2 = 4.0/p->fBall2;
	fNorm = 0.5*M_1_PI*sqrt(ih2)*ih2;
	for (i=0;i<nSmooth;++i) {
		r2 = nnList[i].fDist2*ih2;
		rs = 2.0 - sqrt(r2);
		if (r2 < 1.0) rs = (1.0 - 0.75*rs*r2);
		else rs = 0.25*rs*rs*rs;
		rs *= fNorm;
		q = nnList[i].pPart;
		for (j=0;j<3;++j) {
			p->vMean[j] += rs*q->fMass/q->fDensity*q->v[j];
			q->vMean[j] += rs*p->fMass/p->fDensity*p->v[j];
			}
		}
	}

#endif


#ifdef GASOLINE

void initHsmDivv(void *p)
{
	((PARTICLE *)p)->fDensity = 0.0;
	((PARTICLE *)p)->fHsmDivv = 0.0;
	}

void combHsmDivv(void *p1,void *p2)
{
	((PARTICLE *)p1)->fDensity += ((PARTICLE *)p2)->fDensity;
	((PARTICLE *)p1)->fHsmDivv += ((PARTICLE *)p2)->fHsmDivv;
	}

void postHsmDivv(PARTICLE *p,SMF *smf)
{
	p->fHsmDivv *= 0.5*sqrt(p->fBall2)/p->fDensity;
	}

void HsmDivv(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	PARTICLE *q;
	FLOAT ih2,r2,r,rs,rs1;
	FLOAT dx,dy,dz,dvx,dvy,dvz,dvdotdr,fHsmDivv,fDensity;
	FLOAT fNorm,fNorm1,fTmp;
	int i;

	ih2 = 4.0/p->fBall2;
	fNorm = M_1_PI*sqrt(ih2)*ih2;
	fNorm1 = fNorm*ih2*smf->a;	/* converts to physical velocities */
	fDensity = 0.0;
	fHsmDivv = 0.0;
	for (i=0;i<nSmooth;++i) {
		q = nnList[i].pPart;
		r2 = nnList[i].fDist2*ih2;
		dx = nnList[i].dx;
		dy = nnList[i].dy;
		dz = nnList[i].dz;
		r = sqrt(r2);
		rs = 2.0 - r;
		if (r2 < 1.0) {
			rs = (1.0 - 0.75*rs*r2);
			rs1 = -3 + 2.25*r;
			}
		else {
			rs1 = -0.75*rs*rs/r;
			rs = 0.25*rs*rs*rs;
			}
		dvx = p->vPred[0] - q->vPred[0];
		dvy = p->vPred[1] - q->vPred[1];
		dvz = p->vPred[2] - q->vPred[2];
		dvdotdr = dvx*dx + dvy*dy + dvz*dz + nnList[i].fDist2*smf->H;
		rs1 *= dvdotdr;
		fDensity += rs*q->fMass;
		fHsmDivv -= rs1*q->fMass;
 		}
	p->fDensity = fNorm*fDensity;
	p->fHsmDivv = fNorm1*fHsmDivv;
	}

void HsmDivvSym(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	PARTICLE *q;
	FLOAT ih2,r2,r,rs,rs1;
	FLOAT dx,dy,dz,dvx,dvy,dvz,dvdotdr;
	FLOAT fNorm,fNorm1;
	int i;

	ih2 = 4.0/p->fBall2;
	fNorm = 0.5*M_1_PI*sqrt(ih2)*ih2;
	fNorm1 = fNorm*ih2*smf->a;		/* converts to physical velocities */
	for (i=0;i<nSmooth;++i) {
		q = nnList[i].pPart;
		r2 = nnList[i].fDist2*ih2;
		dx = nnList[i].dx;
		dy = nnList[i].dy;
		dz = nnList[i].dz;
		r = sqrt(r2);
		rs = 2.0 - r;
		if (r2 < 1.0) {
			rs = (1.0 - 0.75*rs*r2);
			rs1 = -3 + 2.25*r;
			}
		else {
			rs1 = -0.75*rs*rs/r;
			rs = 0.25*rs*rs*rs;
			}
		rs *= fNorm;
		rs1 *= fNorm1;
		dvx = p->vPred[0] - q->vPred[0];
		dvy = p->vPred[1] - q->vPred[1];
		dvz = p->vPred[2] - q->vPred[2];
		dvdotdr = dvx*dx + dvy*dy + dvz*dz + nnList[i].fDist2*smf->H;
		rs1 *= dvdotdr;
		p->fDensity += rs*q->fMass;
		q->fDensity += rs*p->fMass;
		p->fHsmDivv -= rs1*q->fMass;
		q->fHsmDivv -= rs1*p->fMass;
 		}
	}


void initEthdotBV(void *p)
{
	((PARTICLE *)p)->du = 0.0;
	}

void combEthdotBV(void *p1,void *p2)
{
	((PARTICLE *)p1)->du += ((PARTICLE *)p2)->du;
	}

void EthdotBVSym(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	PARTICLE *q;
	FLOAT ih2,r2,r,rs;
	FLOAT dx,dy,dz,dvx,dvy,dvz,dvdotdr;
	FLOAT fNorm;
	int i;
	FLOAT qi,qj,Viscij,Eijp,Eijq;

	ih2 = 4.0/p->fBall2;
	fNorm = 0.5*M_1_PI*sqrt(ih2)*ih2*ih2;
	for (i=0;i<nSmooth;++i) {
		q = nnList[i].pPart;
		if (p == q) continue;
		r2 = nnList[i].fDist2*ih2;
		dx = nnList[i].dx;
		dy = nnList[i].dy;
		dz = nnList[i].dz;
		r = sqrt(r2);
		if (r2 < 1.0) {
			rs = -3 + 2.25*r;
			}
		else {
			rs = 2.0 - r;
			rs *= -0.75*rs/r;
			}
		rs *= fNorm;
		dvx = p->vPred[0] - q->vPred[0];
		dvy = p->vPred[1] - q->vPred[1];
		dvz = p->vPred[2] - q->vPred[2];
		dvdotdr = dvx*dx + dvy*dy + dvz*dz + nnList[i].fDist2*smf->H;
		rs *= dvdotdr;
		if (dvdotdr > 0) {
			Viscij = 0;
			}
		else {
			if (p->fHsmDivv < 0) {
				qi = (smf->beta*p->fHsmDivv - smf->algam*sqrt(p->u))*
					p->fHsmDivv;
				}
			else {
				qi = 0.0;
				}
			if (q->fHsmDivv < 0) {
				qj = (smf->beta*q->fHsmDivv - smf->algam*sqrt(q->u))*
					q->fHsmDivv;
				}
			else {
				qj = 0.0;
				}
			Viscij = qi/p->fDensity + qj/q->fDensity;
			}
		if (smf->bGeometric) {
			Eijp = 0.5*Viscij + (smf->gamma - 1)*
				sqrt((p->u*q->u)/(p->fDensity*q->fDensity));
			Eijq = Eijp;
			}
		else {
			Eijp = 0.5*Viscij + (smf->gamma - 1)*p->u/p->fDensity;
			Eijq = 0.5*Viscij + (smf->gamma - 1)*q->u/q->fDensity;
			}
		p->du += rs*q->fMass*Eijp;
		q->du += rs*p->fMass*Eijq;
 		}
	}


void initAccsph(void *p)
{
	((PARTICLE *)p)->a[0] = 0.0;	
	((PARTICLE *)p)->a[1] = 0.0;	
	((PARTICLE *)p)->a[2] = 0.0;	
	}

void combAccsph(void *p1,void *p2)
{
	((PARTICLE *)p1)->a[0] += ((PARTICLE *)p2)->a[0];
	((PARTICLE *)p1)->a[1] += ((PARTICLE *)p2)->a[1];
	((PARTICLE *)p1)->a[2] += ((PARTICLE *)p2)->a[2];	
	}

void AccsphBVSym(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	PARTICLE *q;
	FLOAT ih2,r2,r,rs;
	FLOAT dx,dy,dz,dvx,dvy,dvz,dvdotdr;
	FLOAT fNorm;
	int i;
	FLOAT qi,qj,Viscij,aij;

	ih2 = 4.0/p->fBall2;
	fNorm = 0.5*M_1_PI*sqrt(ih2)*ih2*ih2;

	/*
	 ** We probably don't want this. This needs to be changed.
	 */
	fNorm *= smf->a;	/* converts to physical accelerations */

	for (i=0;i<nSmooth;++i) {
		q = nnList[i].pPart;
		r2 = nnList[i].fDist2*ih2;
		dx = nnList[i].dx;
		dy = nnList[i].dy;
		dz = nnList[i].dz;
		r = sqrt(r2);
		if (r2 < 1.0) {
			rs = -3 + 2.25*r;
			}
		else {
			rs = 2.0 - r;
			rs *= -0.75*rs/r;
			}
		rs *= fNorm;
		dvx = p->vPred[0] - q->vPred[0];
		dvy = p->vPred[1] - q->vPred[1];
		dvz = p->vPred[2] - q->vPred[2];
		dvdotdr = dvx*dx + dvy*dy + dvz*dz + nnList[i].fDist2*smf->H;
		if (dvdotdr > 0) {
			Viscij = 0;
			}
		else {
			if (p->fHsmDivv < 0) {
				qi = (smf->beta*p->fHsmDivv - smf->algam*sqrt(p->u))*
					p->fHsmDivv;
				}
			else {
				qi = 0.0;
				}
			if (q->fHsmDivv < 0) {
				qj = (smf->beta*q->fHsmDivv - smf->algam*sqrt(q->u))*
					q->fHsmDivv;
				}
			else {
				qj = 0.0;
				}
			Viscij = qi/p->fDensity + qj/q->fDensity;
			}
		if (smf->bGeometric) {
			aij = 2.0*(smf->gamma - 1)*sqrt(p->u*q->u)/
				sqrt(p->fDensity*q->fDensity) + Viscij;
			}
		else {
			aij = (smf->gamma - 1)*p->u/p->fDensity +
				(smf->gamma - 1)*q->u/q->fDensity + Viscij;
			}
		p->a[0] -= rs*q->fMass*aij*dx;
		p->a[1] -= rs*q->fMass*aij*dy;
		p->a[2] -= rs*q->fMass*aij*dz;
		q->a[0] += rs*p->fMass*aij*dx;
		q->a[1] += rs*p->fMass*aij*dy;
		q->a[2] += rs*p->fMass*aij*dz;
 		}
	}

#endif

#ifdef COLLISIONS

void
FindRejects(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	/*
	 ** Checks "nSmooth" neighbours of "p" for Hill radius (or physical
	 ** radius) overlap. Note that only planetesimals with "iOrder"s larger
	 ** than "p->iOrder" are considered. This means at most N/2 particles
	 ** will be rejected.
	 */

	PARTICLE *pn;
	double a,r2,sh,sr;
	int i;

	if (p->iColor != PLANETESIMAL) return;

	a = sqrt(p->r[0]*p->r[0] + p->r[1]*p->r[1]); /* approx. orbital distance */
	for (i=0;i<nSmooth;++i) {
		pn = nnList[i].pPart;
		if (pn->iColor != PLANETESIMAL || pn->iOrder <= p->iOrder) continue;
		r2 = nnList[i].fDist2;
		sh = a*(p->fRedHill + pn->fRedHill);
		sr = 2*(p->fSoft + pn->fSoft); /* radius = 2 * softening */
		if (r2 < sh*sh || r2 < sr*sr) {
			p->dTEnc = -1.0; /* rejects have dTEnc < 0 (see REJECT() macro) */
			return; /* one reject is enough */
			}
		}
	}

void
SetTimeStep(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	/*DEBUG out of date
	 * Determines time-step for particle "p" based on distance and relative
	 * radial velocity of the "nSmooth" nearby particles given in "nnList"
	 * [that are currently approaching the target particle]. The time step is
	 * returned in p->dt and has a value no greater than the initial value of
	 * p->dt on entry.
	 *
	 */

	PARTICLE *pn;
	double r2,dt;
	int i;

	for (i=0;i<nSmooth;++i) {
		pn = nnList[i].pPart;
		if (pn == p)
			continue;
		r2 = nnList[i].fDist2;
		dt = 0.2 * sqrt(r2 * sqrt(r2) / (p->fMass + pn->fMass));
		if (dt < p->dt) p->dt = dt;
		}
	}

void
CheckForCollision(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	/*
	 ** Checks whether particle "p" will collide with any of its "nSmooth"
	 ** nearest neighbours in "nnList" during drift interval smf->dStart to
	 ** smf->dEnd, relative to current time. If any collisions are found,
	 ** the relative time to the one that will occur first is noted in
	 ** pkd->dImpactTime and relevant info of the particle and its collider
	 ** are stored in the pkd->Collider1 & pkd->Collider2 structures.
	 ** NOTE: the iIndex of particle "p" is calculated as (p - pkd->pStore).
	 */

	PKD pkd = smf->pkd;
	PARTICLE *pn;
	FLOAT vx,vy,vz,rdotv,v2,sr,D,dt;
	int i;

	assert(p->iActive);

	for (i=0;i<nSmooth;++i) {
		pn = nnList[i].pPart;
		if (pn == p || !(pn->iActive))
			continue;
		if (COLLISION(pkd->dImpactTime) &&
			(pkd->Collider1.id.iPid == nnList[i].iPid &&
			 pkd->Collider1.id.iIndex == nnList[i].iIndex &&
			 pkd->Collider2.id.iPid == pkd->idSelf &&
			 pkd->Collider2.id.iIndex == p - pkd->pStore))
			continue; /* skip if same colliders but order reversed */
		vx = p->v[0] - pn->v[0];
		vy = p->v[1] - pn->v[1];
		vz = p->v[2] - pn->v[2];
		rdotv = nnList[i].dx*vx + nnList[i].dy*vy + nnList[i].dz*vz;
		if (rdotv >= 0)
			continue; /* skip if particle not approaching */
		v2 = vx*vx + vy*vy + vz*vz;
		sr = 2*(p->fSoft + pn->fSoft); /* softening = 0.5 * particle radius */
		D = 1 - v2*(nnList[i].fDist2 - sr*sr)/(rdotv*rdotv);
		if (D < 0)
			continue; /* no real solutions ==> no collision */
		D = sqrt(D);
		if (smf->dStart == 0) {
#ifdef FIX_COLLAPSE
			if (D >= 1) {
				dt = rdotv*(D - 1)/v2;
				if (dt < pkd->dImpactTime) {/* take most negative */
					int k;
					if (sqrt(nnList[i].fDist2)/sr - 1 > 0.01)
						(void) fprintf(stderr,"***LARGE OVERLAP (%f%%)***\n",
									   100*(sqrt(nnList[i].fDist2)/sr - 1));
					fprintf(stderr,"POSITION FIX %i & %i D=%e dt=%e\n",
							p->iOrder,pn->iOrder,D,dt);
					pkd->dImpactTime = dt;
					partToCollider(p,pkd->idSelf,p - pkd->pStore,&pkd->Collider1);
					partToCollider(pn,nnList[i].iPid,nnList[i].iIndex,&pkd->Collider2);
					continue;
					}
				}
#else /* FIX_COLLAPSE */
			if (D >= 1)
				fprintf(stderr,"OVERLAP! %i (r=%e,%e,%e,iRung=%i) &"
						" %i (r=%e,%e,%e,iRung=%i)"
						" v=%e,%e,%e rv=%e sr=%e d=%e D=%e\n",
						p->iOrder,p->r[0],p->r[1],p->r[2],p->iRung,
						pn->iOrder,pn->r[0],pn->r[1],pn->r[2],pn->iRung,
						vx,vy,vz,rdotv,sr,sqrt(nnList[i].fDist2) - sr,D);
			assert(D < 1); /* particles must not touch or overlap initially */
#endif /* !FIX_COLLAPSE */
			}
		dt = rdotv*(D - 1)/v2; /* minimum time to surface contact */
		if (dt > smf->dStart && dt <= smf->dEnd) {
			if (COLLISION(pkd->dImpactTime)) { /* if not first collision */
				if (dt > pkd->dImpactTime)
					continue; /* skip if this one will happen later */
				/*
				 ** ASSERT: Can't handle multiple simultaneous collisions...
				 */
				assert(dt < pkd->dImpactTime);
				}
			/* currently only support collisions between planetesimals... */
			assert(p->iColor == PLANETESIMAL && pn->iColor == PLANETESIMAL);
			pkd->dImpactTime = dt;
			partToCollider(p,pkd->idSelf,p - pkd->pStore,&pkd->Collider1);
			partToCollider(pn,nnList[i].iPid,nnList[i].iIndex,&pkd->Collider2);
			}
		}

#ifdef SAND_PILE
	if (p->v[2] < 0) { /* check for floor collision */
		dt = (2*p->fSoft - p->r[2]) / p->v[2];
		if (dt > smf->dStart && dt <= smf->dEnd && dt < pkd->dImpactTime) {
			pkd->dImpactTime = dt;
			partToCollider(p,pkd->idSelf,p - pkd->pStore,&pkd->Collider1);
			pkd->Collider2.id.iPid = -1;
			}
		}
#endif /* SAND_PILE */

#ifdef IN_A_BOX
	{
	int j;
	for (i=0;i<3;i++) {
		for (j=-1;j<=1;j+=2) {
			if (j*p->v[i] > 0) { /* possible wall collision */
				dt = (BOX_HALF_SIZE - j*p->r[i] - 2*p->fSoft) / (j*p->v[i]);
				if (dt > smf->dStart && dt <= smf->dEnd && dt < pkd->dImpactTime) {
					pkd->dImpactTime = dt;
					partToCollider(p,pkd->idSelf,p - pkd->pStore,&pkd->Collider1);
					pkd->Collider2.id.iPid = -1 - i;
					}
				}
			}
		}
	}
#endif /* IN_A_BOX */
	}

#endif /* COLLISIONS */
