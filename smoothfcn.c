#include <math.h>
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
	fNorm = M_1_PI*sqrt(ih2)*ih2/p->fDensity;
	for (i=0;i<nSmooth;++i) {
		r2 = nnList[i].fDist2*ih2;
		rs = 2.0 - sqrt(r2);
		if (r2 < 1.0) rs = (1.0 - 0.75*rs*r2);
		else rs = 0.25*rs*rs*rs;
		rs *= fNorm;
		q = nnList[i].pPart;
		for (j=0;j<3;++j) {
			p->vMean[j] += rs*q->fMass*q->v[j];
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
			p->vMean[j] += rs/p->fDensity*q->fMass*q->v[j];
			q->vMean[j] += rs/q->fDensity*p->fMass*p->v[j];
			}
		}
	}

#endif


#ifdef GASOLINE

void initDendivv(void *p)
{
	((PARTICLE *)p)->fDensity = 0.0;
	((PARTICLE *)p)->fDivv = 0.0;
	((PARTICLE *)p)->fTmp[0] = 0.0;
	((PARTICLE *)p)->fTmp[1] = 0.0;
	((PARTICLE *)p)->fTmp[2] = 0.0;
	}

void combDendivv(void *p1,void *p2)
{
	((PARTICLE *)p1)->fDensity += ((PARTICLE *)p2)->fDensity;
	((PARTICLE *)p1)->fDivv += ((PARTICLE *)p2)->fDivv;
	((PARTICLE *)p1)->fTmp[0] += ((PARTICLE *)p2)->fTmp[0];
	((PARTICLE *)p1)->fTmp[1] += ((PARTICLE *)p2)->fTmp[1];
	((PARTICLE *)p1)->fTmp[2] += ((PARTICLE *)p2)->fTmp[2];
	}

void postDendivv(PARTICLE *p)
{
	p->fDivv = p->fDivv/p->fDensity;
	p->fCurlv = sqrt(p->fTmp[0]*p->fTmp[0] + p->fTmp[1]*p->fTmp[1] +
					 p->fTmp[2]*p->fTmp[2])/p->fDensity;
	}

void Dendivv(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	PARTICLE *q;
	float ih2,r2,r,rs,rs1,fDensity;
	float dx,dy,dz,dvx,dvy,dvz,dvdotdr,fDivv,fCurlvx,fCurlvy,fCurlvz;
	float fNorm,fNorm1;
	int i;

	ih2 = 4.0/p->fBall2;
	fNorm = M_1_PI*sqrt(ih2)*ih2;
	fNorm1 = fNorm*ih2*smf->a;	/* converts to physical velocities */
	fDensity = 0.0;
	fDivv = 0.0;
	fCurlvx = 0.0;
	fCurlvy = 0.0;
	fCurlvz = 0.0;
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
			rs1 = rs*rs;
			rs = 0.25*rs1*rs;
			rs1 *= -0.75/r;
			}
		dvx = p->v[0] - q->v[0];
		dvy = p->v[1] - q->v[1];
		dvz = p->v[2] - q->v[2];
		dvdotdr = dvx*dx + dvy*dy + dvz*dz + nnList[i].fDist2*smf->H;
		fDensity += rs*q->fMass;
		fDivv += rs1*q->fMass*dvdotdr;
		fCurlvx += rs1*q->fMass*(dvy*dz - dvz*dy);
		fCurlvy += rs1*q->fMass*(dvz*dx - dvx*dz);
		fCurlvx += rs1*q->fMass*(dvx*dy - dvy*dx);
 		}
	p->fDensity = fNorm*fDensity;
	p->fDivv = fNorm1*fDivv;
	p->fTmp[0] = fNorm1*fCurlvx;
	p->fTmp[1] = fNorm1*fCurlvy;
	p->fTmp[2] = fNorm1*fCurlvz;
	}

void DendivvSym(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	PARTICLE *q;
	float ih2,r2,r,rs,rs1;
	float dx,dy,dz,dvx,dvy,dvz,dvdotdr;
	float fNorm,fNorm1;
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
			rs1 = rs*rs;
			rs = 0.25*rs1*rs;
			rs1 *= -0.75/r;
			}
		rs *= fNorm;
		rs1 *= fNorm1;
		dvx = p->v[0] - q->v[0];
		dvy = p->v[1] - q->v[1];
		dvz = p->v[2] - q->v[2];
		dvdotdr = dvx*dx + dvy*dy + dvz*dz + nnList[i].fDist2*smf->H;
		p->fDensity += rs*q->fMass;
		q->fDensity += rs*p->fMass;
		p->fDivv -= rs1*q->fMass*dvdotdr;
		q->fDivv -= rs1*p->fMass*dvdotdr;
		p->fTmp[0] += rs1*q->fMass*(dvy*dz - dvz*dy);
		q->fTmp[0] += rs1*p->fMass*(dvy*dz - dvz*dy);
		p->fTmp[1] += rs1*q->fMass*(dvz*dx - dvx*dz);
		q->fTmp[1] += rs1*p->fMass*(dvz*dx - dvx*dz);
		p->fTmp[2] += rs1*q->fMass*(dvx*dy - dvy*dx);
		q->fTmp[2] += rs1*p->fMass*(dvx*dy - dvy*dx);
 		}
	}


void initEthdotBV(void *p)
{
	((PARTICLE *)p)->fEthdot = 0.0;
	}

void combEthdotBV(void *p1,void *p2)
{
	((PARTICLE *)p1)->fEthdot += ((PARTICLE *)p2)->fEthdot;
	}

void postEthdotBV(PARTICLE *p)
{
	}

void EthdotBV(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	PARTICLE *q;
	float ih2,r2,r,rs,rs1,fDensity;
	float dx,dy,dz,dvx,dvy,dvz,dvdotdr,fDivv,fCurlvx,fCurlvy,fCurlvz;
	float fNorm,fNorm1;
	int i;

	ih2 = 4.0/p->fBall2;
	fNorm = M_1_PI*sqrt(ih2)*ih2;
	fNorm1 = fNorm*ih2*smf->a;	/* converts to physical velocities */
	fDensity = 0.0;
	fDivv = 0.0;
	fCurlvx = 0.0;
	fCurlvy = 0.0;
	fCurlvz = 0.0;
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
			rs1 = rs*rs;
			rs = 0.25*rs1*rs;
			rs1 *= -0.75/r;
			}
		dvx = p->v[0] - q->v[0];
		dvy = p->v[1] - q->v[1];
		dvz = p->v[2] - q->v[2];
		dvdotdr = dvx*dx + dvy*dy + dvz*dz + nnList[i].fDist2*smf->H;
		qi = (alpha*p->csound + beta*hsmDivvp)*hsmDivvp;
		qj = (alpha*q->csound + beta*hsmDivvnbi)*hsmDivvnbi;
		if (bGeometric) {
			eijp = dvdotdr*(p->csound*q->csound/
							(gamma*sqrt(p->fDensity*q->fDensity)) +
							0.5*(qi/p->fDensity + qj/q->fDensity));
			eijnbi = eijp;
			}
		else {
			eijp = dvdotdr*(p->csound*p->csound/(gamma*p->fDensity) +
							0.5*(qi/p->fDensity + qj/q->fDensity));
			eijnbi = dvdotdr*(q->csound*q->csound/(gamma*q->fDensity) +
							0.5*(qi/p->fDensity + qj/q->fDensity));
			}
		fDensity += rs*q->fMass*eijp;
 		}
	p->fEthdot = fNorm*fEthDot;
	}

void EthdotBVSym(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	PARTICLE *q;
	float ih2,r2,r,rs,rs1;
	float dx,dy,dz,dvx,dvy,dvz,dvdotdr;
	float fNorm,fNorm1;
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
			rs1 = rs*rs;
			rs = 0.25*rs1*rs;
			rs1 *= -0.75/r;
			}
		rs *= fNorm;
		rs1 *= fNorm1;
		dvx = p->v[0] - q->v[0];
		dvy = p->v[1] - q->v[1];
		dvz = p->v[2] - q->v[2];
		dvdotdr = dvx*dx + dvy*dy + dvz*dz + nnList[i].fDist2*smf->H;
		p->fDensity += rs*q->fMass;
		q->fDensity += rs*p->fMass;
		p->fDivv -= rs1*q->fMass*dvdotdr;
		q->fDivv -= rs1*p->fMass*dvdotdr;
		p->fTmp[0] += rs1*q->fMass*(dvy*dz - dvz*dy);
		q->fTmp[0] += rs1*p->fMass*(dvy*dz - dvz*dy);
		p->fTmp[1] += rs1*q->fMass*(dvz*dx - dvx*dz);
		q->fTmp[1] += rs1*p->fMass*(dvz*dx - dvx*dz);
		p->fTmp[2] += rs1*q->fMass*(dvx*dy - dvy*dx);
		q->fTmp[2] += rs1*p->fMass*(dvx*dy - dvy*dx);
 		}
	}




#endif









