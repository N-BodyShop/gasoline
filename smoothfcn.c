#include <math.h>
#include <assert.h>
#include "smoothfcn.h"

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
	float ih2,r2,rs,fDensity;
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
	float fNorm,ih2,r2,rs;
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
	float fNorm,ih2,r2,rs;
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
	float fNorm,ih2,r2,rs;
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
	((PARTICLE *)p)->fHsmDivv = 0.0;
	((PARTICLE *)p)->fRhoDivv = 0.0;
	((PARTICLE *)p)->fCutVisc = 0.0;
	}

void combHsmDivv(void *p1,void *p2)
{
	((PARTICLE *)p1)->fHsmDivv += ((PARTICLE *)p2)->fHsmDivv;
	((PARTICLE *)p1)->fRhoDivv += ((PARTICLE *)p2)->fRhoDivv;
	((PARTICLE *)p1)->fCutVisc += ((PARTICLE *)p2)->fCutVisc;
	}

void postHsmDivv(PARTICLE *p,SMF *smf)
{
	p->fHsmDivv *= 0.5*sqrt(p->fBall2);
	}

void HsmDivv(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	PARTICLE *q;
	float ih2,r2,r,rs;
	float dx,dy,dz,dvx,dvy,dvz,dvdotdr,fHsmDivv,fRhoDivv,fCutVisc;
	float fNorm,fTmp;
	int i;

	ih2 = 4.0/p->fBall2;
	fNorm = M_1_PI*sqrt(ih2)*ih2*ih2;
	fNorm *= smf->a;	/* converts to physical velocities */
	fHsmDivv = 0.0;
	fRhoDivv = 0.0;
	fCutVisc = 0.0;
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
		dvx = p->vPred[0] - q->vPred[0];
		dvy = p->vPred[1] - q->vPred[1];
		dvz = p->vPred[2] - q->vPred[2];
		dvdotdr = dvx*dx + dvy*dy + dvz*dz + nnList[i].fDist2*smf->H;
		rs *= dvdotdr;
		fRhoDivv += rs*q->fMass;
		if (dvdotdr < 0) fCutVisc += rs*q->fMass;
		fHsmDivv -= rs*q->fMass/q->fDensity;
 		}
	p->fRhoDivv = fNorm*fRhoDivv;
	p->fCutVisc = fNorm*fCutVisc;
	p->fHsmDivv = fNorm*fHsmDivv;
	}

void HsmDivvSym(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	PARTICLE *q;
	float ih2,r2,r,rs;
	float dx,dy,dz,dvx,dvy,dvz,dvdotdr;
	float fNorm;
	int i;

	ih2 = 4.0/p->fBall2;
	fNorm = 0.5*M_1_PI*sqrt(ih2)*ih2*ih2;
	fNorm *= smf->a;		/* converts to physical velocities */
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
		rs *= dvdotdr;
		p->fRhoDivv += rs*q->fMass;
		p->fHsmDivv -= rs*q->fMass/q->fDensity;
		q->fRhoDivv += rs*p->fMass;
		q->fHsmDivv -= rs*p->fMass/p->fDensity;
		if (dvdotdr < 0) {
			p->fCutVisc += rs*q->fMass;
			q->fCutVisc += rs*p->fMass;
			}
 		}
	}


void initGeomBV(void *p)
{
	((PARTICLE *)p)->A = 0.0;
	((PARTICLE *)p)->B = 0.0;
	}

void combGeomBV(void *p1,void *p2)
{
	((PARTICLE *)p1)->A += ((PARTICLE *)p2)->A;
	((PARTICLE *)p1)->B += ((PARTICLE *)p2)->B;
	}

void postGeomBV(PARTICLE *p,SMF *smf)
{
	p->A *= (smf->gamma-1)/sqrt(p->fDensity);
	if (p->fHsmDivv < 0) {
		p->A -= 0.5*smf->algam*p->fHsmDivv/p->fDensity*p->fCutVisc;
		p->B += 0.5*smf->beta*p->fHsmDivv*p->fHsmDivv/p->fDensity*p->fCutVisc;
		}
	}

void GeomBVSym(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	PARTICLE *q;
	float ih2,r2,r,rs;
	float dx,dy,dz,dvx,dvy,dvz,dvdotdr;
	float fNorm;
	int i;

	ih2 = 4.0/p->fBall2;
	fNorm = 0.5*M_1_PI*sqrt(ih2)*ih2*ih2;
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
		rs *= dvdotdr;
		p->A += rs*q->fMass*sqrt(q->u/q->fDensity);
		q->A += rs*p->fMass*sqrt(p->u/p->fDensity);
		if (dvdotdr < 0) {
			if (q->fHsmDivv < 0) {
				p->B += 0.5*rs*q->fMass/q->fDensity*q->fHsmDivv*
					(smf->beta*q->fHsmDivv - smf->algam*sqrt(q->u));
				}
			if (p->fHsmDivv < 0) {
				q->B += 0.5*rs*p->fMass/p->fDensity*p->fHsmDivv*
					(smf->beta*p->fHsmDivv - smf->algam*sqrt(p->u));
				}
			}
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
	float ih2,r2,r,rs;
	float dx,dy,dz,dvx,dvy,dvz,dvdotdr;
	float fNorm;
	int i;
	float qi,qj,Viscij,Eijp,Eijq;

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
	float ih2,r2,r,rs;
	float dx,dy,dz,dvx,dvy,dvz,dvdotdr;
	float fNorm;
	int i;
	float qi,qj,Viscij,aij;

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

#ifdef PLANETS

#define X 0
#define Y 1
#define Z 2

#define DT_TSC 1.0 /*DEBUG*/

#define MIN(x,y) ((x) < (y) ? (x) : (y)) /*DEBUG*/

void SetTimeStep(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	/*
	 * Determines time-step for particle "p" based on distance and relative
	 * radial velocity of the "nSmooth" nearby particles given in "nnList"
	 * that are currently approaching the target particle. The time step is
	 * returned in p->dt and has a value greater than zero and no more than
	 * the initial value of p->dt on entry.
	 *
	 */

	/*DEBUG currently nSmooth = 63*/

	PARTICLE *pn;
	float vx,vy,vz,rdotv,dt; /* these ought to be doubles! */
	int i;

	p->dt = 10.0; /*DEBUG!!!*/

	for (i=0;i<nSmooth;++i) {
		pn = nnList[i].pPart;
		if (pn == p)
			continue;
		vx = p->v[X] - pn->v[X];
		vy = p->v[Y] - pn->v[Y];
		vz = p->v[Z] - pn->v[Z];
		rdotv = nnList[i].dx * vx + nnList[i].dy * vy + nnList[i].dz * vz;
		if (rdotv >= 0) /* skip if particle not approaching */
			continue;
		dt = - DT_TSC * nnList[i].fDist2 / rdotv;
		p->dt = MIN(p->dt, dt);
		}
	}

#undef MIN /*DEBUG*/

#undef DT_TSC /*DEBUG*/

#define RADIUS(p) ((p)->fSoft) /* softening length <=> particle radius */

#define NO_COLLIDER -1 /*DEBUG*/

#define SQ(x) ((x)*(x)) /*DEBUG*/

void CheckForCollision(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	/*
	 * Checks whether particle "p" will collide with any of its "nSmooth"
	 * nearest neighbours given in "nnList" over interval smf->dStart to
	 * smf->dEnd, relative to current time. If any collisions are found,
	 * the time increment to the one that will occur first is returned in
	 * p->dt and the value of p->iCollider is set to the iOrder of the
	 * colliding particle.
	 *
	 */

	PARTICLE *pn;
	float vx,vy,vz,rdotv,v2,D,dt; /* these ought to be doubles! */
	int i;

	p->dCollTime = smf->dEnd; /* time interval cannot exceed smf->dEnd */
	p->iCollider = NO_COLLIDER;
	for (i=0;i<nSmooth;++i) {
		pn = nnList[i].pPart;
		if (pn == p)
			continue;
		vx = p->v[X] - pn->v[X];
		vy = p->v[Y] - pn->v[Y];
		vz = p->v[Z] - pn->v[Z];
		rdotv = nnList[i].dx * vx + nnList[i].dy * vy + nnList[i].dz * vz;
		if (rdotv >= 0) /* skip if particle not approaching */
			continue;
		v2 = SQ(vx) + SQ(vy) + SQ(vz);
		D = SQ(rdotv) - v2 * (nnList[i].fDist2 - SQ(RADIUS(p) + RADIUS(pn)));
		if (D < 0) /* no real solutions ==> no collision */
			continue;
		D = sqrt(D);
		assert(D <= 1); /* colliding particles must not already overlap */
		dt = rdotv * (D - 1) / v2; /* minimum time to contact */
		if (dt > smf->dStart && dt <= p->dCollTime) {
			assert(p->iCollider == NO_COLLIDER || dt < p->dCollTime); /* can't handle multiple simultaneous collisions */
			p->dCollTime = dt;
			p->iCollider = pn->iOrder;
			}
		}
	}

#undef SQ /*DEBUG*/

#undef NO_COLLIDER /*DEBUG*/

#undef RADIUS

#undef Z
#undef Y
#undef X

#endif /* PLANETS */
