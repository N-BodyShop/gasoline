#include <stdio.h>
#include <math.h>
#include "pkd.h"
#include "grav.h"


void pkdBucketInteract(PKD pkd,int iBucket)
{
	PARTICLE *p;
	KDN *pkdn;
	ILP *ilp;
	ILCS *ilcs;
	ILCN *ilcn;
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
	ilcs = pkd->ilcs;
	ilcn = pkd->ilcn;
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
			for (j=0;j<pkd->nCellSoft;++j) {
				dx = x - ilcs[j].x;
				dy = y - ilcs[j].y;
				dz = z - ilcs[j].z;
				d2 = dx*dx + dy*dy + dz*dz;
				twoh = h + ilcs[j].h;
				SPLINE(d2,twoh,a,b,c,d);
				qirx = ilcs[j].xx*dx + ilcs[j].xy*dy + ilcs[j].xz*dz;
				qiry = ilcs[j].xy*dx + ilcs[j].yy*dy + ilcs[j].yz*dz;
				qirz = ilcs[j].xz*dx + ilcs[j].yz*dy + ilcs[j].zz*dz;
				qir = 0.5*(qirx*dx + qiry*dy + qirz*dz);
				tr = 0.5*(ilcs[j].xx + ilcs[j].yy + ilcs[j].zz);
				qir3 = b*ilcs[j].m + d*qir - c*tr;
				fPot -= a*ilcs[j].m + c*qir - b*tr;
				ax -= qir3*dx - c*qirx;
				ay -= qir3*dy - c*qiry;
				az -= qir3*dz - c*qirz;
				}
			for (j=0;j<pkd->nCellNewt;++j) {
				dx = x - ilcn[j].x;
				dy = y - ilcn[j].y;
				dz = z - ilcn[j].z;
				d2 = dx*dx + dy*dy + dz*dz;
				a = 1.0/sqrt(d2);
				b = a*a*a;
				c = 3.0*b*a*a;
				d = 5.0*c*a*a;
				qirx = ilcn[j].xx*dx + ilcn[j].xy*dy + ilcn[j].xz*dz;
				qiry = ilcn[j].xy*dx + ilcn[j].yy*dy + ilcn[j].yz*dz;
				qirz = ilcn[j].xz*dx + ilcn[j].yz*dy + ilcn[j].zz*dz;
				qir = 0.5*(qirx*dx + qiry*dy + qirz*dz);
				tr = 0.5*(ilcn[j].xx + ilcn[j].yy + ilcn[j].zz);
				qir3 = b*ilcn[j].m + d*qir - c*tr;
				fPot -= a*ilcn[j].m + c*qir - b*tr;
				ax -= qir3*dx - c*qirx;
				ay -= qir3*dy - c*qiry;
				az -= qir3*dz - c*qirz;
				}
			}
		else {
			/*
			 ** Monopole-only interactions.
			 */
			for (j=0;j<pkd->nCellSoft;++j) {
				dx = x - ilcs[j].x;
				dy = y - ilcs[j].y;
				dz = z - ilcs[j].z;
				d2 = dx*dx + dy*dy + dz*dz;
				twoh = h + ilcs[j].h;
				SPLINE(d2,twoh,a,b,c,d);
				b *= ilcs[j].m;
				fPot -= ilcs[j].m*a;
				ax -= dx*b;
				ay -= dy*b;
				az -= dz*b;
				}
			for (j=0;j<pkd->nCellNewt;++j) {
				dx = x - ilcn[j].x;
				dy = y - ilcn[j].y;
				dz = z - ilcn[j].z;
				a = 1.0/sqrt(dx*dx + dy*dy + dz*dz);
				b = ilcn[j].m*a*a*a;
				fPot -= ilcn[j].m*a;
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











