#ifdef COLLISIONS

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "pkd.h"
#include "collision.h"

void
partToCollider(PARTICLE *p,int iPid,int iIndex,COLLIDER *c)
{
	int i;

	c->id.iPid = iPid;
	c->id.iIndex = iIndex;
	c->id.iOrder = p->iOrder;
	c->fMass = p->fMass;
	c->fRadius = 2*p->fSoft;
	for (i=0;i<3;i++) {
		c->r[i] = p->r[i];
		c->v[i] = p->v[i];
		c->w[i] = p->w[i];
		}
	c->dt = p->dt;
	}

void
colliderToPart(COLLIDER *c,PARTICLE *p)
{
	/* id & dt info ignored */

	int i;

	p->fMass = c->fMass;
	p->fSoft = 0.5*c->fRadius;
	for (i=0;i<3;i++) {
		p->r[i] = c->r[i];
		p->v[i] = c->v[i];
		p->w[i] = c->w[i];
		}
	}

void
pkdDoCollision(PKD pkd,int iOutcomes,double dEpsN,double dEpsT,
			   COLLIDER *Collider1,COLLIDER *Collider2,double *pdImpactEnergy,
			   int *piOutcome,COLLIDER *pOut,int *pnOut)
{
	/* Collider2 may or may not be local */

	COLLIDER *p1=&pkd->Collider1,*p2=&pkd->Collider2,*p;
	double dt=pkd->dImpactTime,v2,ve2,dImpactEnergy;
	int bReturn = (pkd->idSelf == p1->id.iPid),iOutcome;
	int i,j,n;

#ifdef SAND_PILE
	if (p2->id.iPid == -1) { /* floor collision */
		assert(p1->id.iPid == pkd->idSelf);
		assert(pkd->pStore[p1->id.iIndex].bStuck == 0);
		for (i=0;i<3;i++) {
			p1->r[i] += p1->v[i]*dt;
			p1->v[i] = 0;
			p1->w[i] = 0;
			}
		colliderToPart(p1,&pkd->pStore[p1->id.iIndex]);
		pkd->pStore[p1->id.iIndex].bStuck = 1;
		*pdImpactEnergy = 0;
		*piOutcome = MERGE;
		*pnOut = 0;
		return;
		}
#endif /* SAND_PILE */

#ifdef IN_A_BOX
	if (p2->id.iPid < 0) { /* wall collision */
		assert(p1->id.iPid == pkd->idSelf);
		for (i=0;i<3;i++) p1->r[i] += p1->v[i]*dt;
		p1->v[-p2->id.iPid - 1] *= -1; /* reflection */
		for (i=0;i<3;i++) p1->r[i] -= p1->v[i]*dt;
		colliderToPart(p1,&pkd->pStore[p1->id.iIndex]);
		*pdImpactEnergy = 0;
		*piOutcome = BOUNCE;
		*pnOut = 0;
		return;
		}
#endif /* IN_A_BOX */

	/* At least one particle must be local */

	assert(p1->id.iPid == pkd->idSelf || p2->id.iPid == pkd->idSelf);

	/* Advance coordinates to impact time */

	for (i=0;i<3;i++) {
		p1->r[i] += p1->v[i]*dt;
		p2->r[i] += p2->v[i]*dt;
	}

	/* Store info with time-advanced coordinates for diagnostic output */

	if (bReturn) {
		*Collider1 = *p1;
		*Collider2 = *p2;
		}

	/* Determine collision outcome */

	v2 = 0;

	for (i=0;i<3;i++)
		v2 += (p2->v[i] - p1->v[i])*(p2->v[i] - p1->v[i]);

	ve2 = 2*(p1->fMass + p2->fMass)/(p1->fRadius + p2->fRadius);

	/*DEBUG impact energy not currently used...*/

	dImpactEnergy = 0.5*p1->fMass*p2->fMass/(p1->fMass + p2->fMass)*v2;

	if (bReturn) *pdImpactEnergy = dImpactEnergy;

	/* Inhibit sliding phase if only bouncing allowed */

	if (iOutcomes == BOUNCE && v2 < 0.33*ve2) dEpsN = 1.0;

	/* Determine collision outcome */

	iOutcome = NONE;

	if (iOutcomes == MERGE || ((iOutcomes & MERGE) && v2 <= ve2)) {
		iOutcome = MERGE;
		pkdMerge(pkd,&p,&n);
		assert(n == 1);
		if (iOutcomes & (BOUNCE | FRAG)) { /* check if spinning too fast */
			double w2max,w2=0;
			w2max = p->fMass/(p->fRadius*p->fRadius*p->fRadius);
			for (i=0;i<3;i++)
				w2 += p->w[i]*p->w[i];
			if (w2 > w2max) {
				free((void *)p);
				iOutcome = BOUNCE;
				pkdBounce(pkd,dEpsN,dEpsT,&p,&n);
				assert(n == 2);
				}
			}
		}
	else if (iOutcomes & BOUNCE) { /* add strength criterion here */
		iOutcome = BOUNCE;
		pkdBounce(pkd,dEpsN,dEpsT,&p,&n);
		assert(n == 2);
		}
	else if (iOutcomes & FRAG) {
		iOutcome = FRAG;
		pkdFrag(pkd,&p,&n);
		assert(n > 0 && n <= MAX_NUM_FRAG);
		}

	assert(iOutcome != NONE); /* SOMETHING has to happen... */

	if (bReturn) *piOutcome = iOutcome;

	assert(n > 0);

	if (bReturn) {
		for (i=0;i<n;i++)
			pOut[i] = p[i];
		*pnOut = n;
		}

	/* Trace output particles back to start of step */

	for (i=0;i<n;i++)
		for (j=0;j<3;j++)
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

void
pkdMerge(PKD pkd,COLLIDER **pOut,int *pnOut)
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

	for (k=0;k<3;k++) {
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

	for (k=0;k<3;k++) {
		p->r[k] = com_pos[k];
		p->v[k] = com_vel[k];
		p->w[k] = ang_mom[k]/i;
		}
	}

void
pkdBounce(PKD pkd,double dEpsN,double dEpsT,COLLIDER **pOut,int *pnOut)
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
	for (i=0;i<3;i++) {
		n[i] = (p2->r[i] - p1->r[i]);
		a += n[i]*n[i];
		}
	a = 1/sqrt(a);
	for (i=0;i<3;i++)
		n[i] *= a;

	s1[0] = r1*(p1->w[1]*n[2] - p1->w[2]*n[1]);
	s1[1] = r1*(p1->w[2]*n[0] - p1->w[0]*n[2]);
	s1[2] = r1*(p1->w[0]*n[1] - p1->w[1]*n[0]);

	s2[0] = r2*(p2->w[2]*n[1] - p2->w[1]*n[2]);
	s2[1] = r2*(p2->w[0]*n[2] - p2->w[2]*n[0]);
	s2[2] = r2*(p2->w[1]*n[0] - p2->w[0]*n[1]);

	for (i=0;i<3;i++)
		u[i] = (p2->v[i] - p1->v[i]) + (s2[i] - s1[i]);

	a = u[0]*n[0] + u[1]*n[1] + u[2]*n[2];
	for (i=0;i<3;i++) {
		un[i] = a*n[i];
		ut[i] = u[i] - un[i];
		}

#ifdef SAND_PILE
	if (pkd->pStore[p1->id.iIndex].bStuck) {
		mu = m2;
		beta = 2.0/7;
		}
	else if (pkd->pStore[p2->id.iIndex].bStuck) {
		mu = m1;
		beta = 2.0/7;
		}
#endif /* SAND_PILE */

	a = (1 + dEpsN);
	b = beta*(1 - dEpsT);
	for (i=0;i<3;i++)
		p[i] = a*un[i] + b*ut[i];

	a = mu*b;
	q[0] = a*(n[1]*u[2] - n[2]*u[1]);
	q[1] = a*(n[2]*u[0] - n[0]*u[2]);
	q[2] = a*(n[0]*u[1] - n[1]*u[0]);

#ifdef SAND_PILE
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
#endif /* SAND_PILE */
	a =   m2/m;
	b = - m1/m;
	c = r1/i1;
	d = r2/i2;
#ifdef SAND_PILE
	}
#endif /* SAND_PILE */

	*pnOut = 2;
	*pOut = (COLLIDER *)malloc(*pnOut*sizeof(COLLIDER));
	assert(*pOut);

	(*pOut)[0] = *p1;
	(*pOut)[1] = *p2;

	p1 = &(*pOut)[0]; /* p1,p2 now point to output structures */
	p2 = &(*pOut)[1];

	for (i=0;i<3;i++) {
		p1->v[i] += a*p[i];
		p2->v[i] += b*p[i];
		p1->w[i] += c*q[i];
		p2->w[i] += d*q[i];
		}
	}

void
pkdFrag(PKD pkd,COLLIDER **pOut,int *pnOut)
{
	/* note: be sure new particles have p->iActive = 1 */
	/* may need to assert(*pnOut <= MAX_NUM_FRAG) */
	/* remember to set id info for logging purposes */
	}

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
	a0 = ei + pi->fRedHill;
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
	h = pi->fRedHill*ai;
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

#endif /* COLLISIONS */
