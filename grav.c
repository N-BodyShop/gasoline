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
	double x,y,z,dx,dy,dz,d2,h,twoh,a,b,c,d;
	double qirx,qiry,qirz,qir,tr,qir3;

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
			SPLINE(d2,twoh,a,b,c,d);
			a *= ilp[j].m;
			b *= ilp[j].m;
			fPot -= a;
			ax -= dx*b;
			ay -= dy*b;
			az -= dz*b;
			}
		if (pkd->iOrder == 2) {	
			/*
			 ** Monopole and Quadrupole interactions.
			 */
			for (j=0;j<pkd->nCell;++j) {
				dx = x - ilc[j].x;
				dy = y - ilc[j].y;
				dz = z - ilc[j].z;
				d2 = dx*dx + dy*dy + dz*dz;
				twoh = h + ilc[j].h;
				SPLINE(d2,twoh,a,b,c,d);
				qirx = ilc[j].xx*dx + ilc[j].xy*dy + ilc[j].xz*dz;
				qiry = ilc[j].xy*dx + ilc[j].yy*dy + ilc[j].yz*dz;
				qirz = ilc[j].xz*dx + ilc[j].yz*dy + ilc[j].zz*dz;
				qir = 0.5*(qirx*dx + qiry*dy + qirz*dz);
				tr = 0.5*(ilc[j].xx + ilc[j].yy + ilc[j].zz);
				qir3 = b*ilc[j].m + d*qir - c*tr;
				fPot -= a*ilc[j].m + c*qir - b*tr;
				ax -= qir3*dx - c*qirx;
				ay -= qir3*dy - c*qiry;
				az -= qir3*dz - c*qirz;
				}
			}
		else {
			/*
			 ** Monopole interactions.
			 */
			for (j=0;j<pkd->nCell;++j) {
				dx = x - ilc[j].x;
				dy = y - ilc[j].y;
				dz = z - ilc[j].z;
				d2 = dx*dx + dy*dy + dz*dz;
				twoh = h + ilc[j].h;
				SPLINE(d2,twoh,a,b,c,d);
				a *= ilc[j].m;
				b *= ilc[j].m;
				fPot -= a;
				ax -= dx*b;
				ay -= dy*b;
				az -= dz*b;
				}
			}
		p[i].fPot = fPot;
		p[i].a[0] = ax;
		p[i].a[1] = ay;
		p[i].a[2] = az;
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
			SPLINE(d2,twoh,a,b,c,d);
			p[j].fPot -= a*p[i].fMass;
			p[j].a[0] -= dx*b*p[i].fMass;
			p[j].a[1] -= dy*b*p[i].fMass;
			p[j].a[2] -= dz*b*p[i].fMass;
			p[i].fPot -= a*p[j].fMass;
			p[i].a[0] += dx*b*p[j].fMass;
			p[i].a[1] += dy*b*p[j].fMass;
			p[i].a[2] += dz*b*p[j].fMass;
			}
		}
	}











