#include <stdio.h>
#include <math.h>
#include "pkd.h"
#include "grav.h"


void pkdBucketInteract(PKD pkd,int iBucket)
{
	PARTICLE *p;
	KDN *pkdn;
	ILP *ilp;
	ILC *ilc;
	int n,i,j;
	double fPot,ax,ay,az;
	double x,y,z,dx,dy,dz,d2,h,twoh,dir,dir3,dir4,dir5;
	double qirx,qiry,qirz,qir,tr,qir3;
#if 0
        double fIMass;
#endif

	/*
	 ** Now process the two interaction lists for each particle.
	 */
	pkdn = &pkd->kdNodes[iBucket];
	p = &pkd->pStore[pkdn->pLower];
	n = pkdn->pUpper - pkdn->pLower + 1;
	ilp = pkd->ilp;
	ilc = pkd->ilc;
	for (i=0;i<n;++i) {
		ax = 0.0;
		ay = 0.0;
		az = 0.0;
		fPot = 0.0;
#if 0
                fIMass = p[i].fMass;
#endif
		x = p[i].r[0];
		y = p[i].r[1];
		z = p[i].r[2];
		h = p[i].fSoft;
		for (j=0;j<pkd->nPart;++j) {
			dx = x - ilp[j].x;
			dy = y - ilp[j].y;
			dz = z - ilp[j].z;
			d2 = dx*dx + dy*dy + dz*dz;
			twoh = h + ilp[j].h;
			SPLINE(d2,twoh,dir,dir3);
			dir *= ilp[j].m;
			dir3 *= ilp[j].m;
			fPot -= dir;
			ax -= dx*dir3;
			ay -= dy*dir3;
			az -= dz*dir3;
#if 0
                        fIMass += ilp[j].m;
#endif
			}
		for (j=0;j<pkd->nCell;++j) {
			dx = x - ilc[j].x;
			dy = y - ilc[j].y;
			dz = z - ilc[j].z;
			d2 = dx*dx + dy*dy + dz*dz;
			twoh = h + ilc[j].h;
			SPLINE(d2,twoh,dir,dir3);
			/*
			 ** Monopole interaction.
			 */
			fPot -= dir*ilc[j].m;
			ax -= dx*dir3*ilc[j].m;
			ay -= dy*dir3*ilc[j].m;
			az -= dz*dir3*ilc[j].m;
#if 0
                        fIMass += ilc[j].m;
#endif
			/*
			 ** Normalize dx,dy and dz for the higher order moments.
			 */
			dx *= dir;
			dy *= dir;
			dz *= dir;
			/*
			 ** Quadrupole interaction.
			 */
			qirx = ilc[j].xx*dx + ilc[j].xy*dy + ilc[j].xz*dz;
			qiry = ilc[j].xy*dx + ilc[j].yy*dy + ilc[j].yz*dz;
			qirz = ilc[j].xz*dx + ilc[j].yz*dy + ilc[j].zz*dz;
			qir = qirx*dx + qiry*dy + qirz*dz;
			tr = ilc[j].xx + ilc[j].yy + ilc[j].zz;
			qir3 = 15*qir - 3*tr;
			dir4 = dir3*dir;
			fPot -= 0.5*dir*dir*dir*(3*qir - tr);
			ax -= 0.5*dir4*(qir3*dx - 6*qirx);
			ay -= 0.5*dir4*(qir3*dy - 6*qiry);
			az -= 0.5*dir4*(qir3*dz - 6*qirz);
			}
		p[i].fPot = fPot;
		p[i].a[0] = ax;
		p[i].a[1] = ay;
		p[i].a[2] = az;
#if 0
		p[i].fIMass = fIMass;
#endif
		}
	/*
	 ** Do the inter-bucket interactions.
	 */
	for (i=0;i<n-1;++i) {
		for (j=i+1;j<n;++j) {
			dx = p[j].r[0] - p[i].r[0];
			dy = p[j].r[1] - p[i].r[1];
			dz = p[j].r[2] - p[i].r[2];
			d2 = dx*dx + dy*dy + dz*dz;
			twoh = p[i].fSoft + p[j].fSoft;
			SPLINE(d2,twoh,dir,dir3);
			p[j].fPot -= dir*p[i].fMass;
			p[j].a[0] -= dx*dir3*p[i].fMass;
			p[j].a[1] -= dy*dir3*p[i].fMass;
			p[j].a[2] -= dz*dir3*p[i].fMass;
#if 0
                        p[j].fIMass += p[i].fMass;
#endif
			p[i].fPot -= dir*p[j].fMass;
			p[i].a[0] += dx*dir3*p[j].fMass;
			p[i].a[1] += dy*dir3*p[j].fMass;
			p[i].a[2] += dz*dir3*p[j].fMass;
#if 0
                        p[i].fIMass += p[j].fMass;
#endif
			}
		}
	}






