#include <stdio.h>
#include <math.h>
#include "mdl.h"
#include "pkd.h"
#include "grav.h"


void pkdBucketInteract(PKD pkd,int iBucket,int iOrder)
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
	double hirx,hiry,hirz,hir,xxyy,xxzz,yyzz,xyzz,xzzz,yzzz;
	double oirx,oiry,oirz,oir,xxy,xyy,xyz,xzz,yzz;
	double dir,dir2,dir3,dir4,dir5,dir6;
	double *sqrttmp,*d2a;

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
			pkd->d2a[j] = dx*dx + dy*dy + dz*dz;
			}
		if (pkd->nPart>0) v_sqrt1(pkd->nPart,pkd->d2a,pkd->sqrttmp);
		for (j=0;j<pkd->nPart;++j) {
			dx = x - ilp[j].x;
			dy = y - ilp[j].y;
			dz = z - ilp[j].z;
			twoh = h + ilp[j].h;
			SPLINE(pkd->sqrttmp[j],pkd->d2a[j],twoh,a,b,c,d);
			a *= ilp[j].m;
			b *= ilp[j].m;
			fPot -= a;
			ax -= dx*b;
			ay -= dy*b;
			az -= dz*b;
			}
		for (j=0;j<pkd->nCellSoft;++j) {
			dx = x - ilcs[j].x;
			dy = y - ilcs[j].y;
			dz = z - ilcs[j].z;
			pkd->d2a[j] = dx*dx + dy*dy + dz*dz;
			}
		if (pkd->nCellSoft>0) v_sqrt1(pkd->nCellSoft,pkd->d2a,pkd->sqrttmp);
		for (j=0;j<pkd->nCellSoft;++j) {
			dx = x - ilcs[j].x;
			dy = y - ilcs[j].y;
			dz = z - ilcs[j].z;
			twoh = h + ilcs[j].h;
			SPLINE(pkd->sqrttmp[j],pkd->d2a[j],twoh,a,b,c,d);
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
			pkd->d2a[j] = dx*dx + dy*dy + dz*dz;
			}
		if (pkd->nCellNewt>0) v_sqrt1(pkd->nCellNewt,pkd->d2a,pkd->sqrttmp);
		for (j=0;j<pkd->nCellNewt;++j) {
			dir = pkd->sqrttmp[j];
			dir2 = dir*dir;
			dir3 = dir2*dir;
			dir4 = dir3*dir;
			dir5 = dir4*dir;
			dir6 = dir5*dir;
			/*
			 ** Normalize dx,dy and dz for the higher order moments.
			 */
			dx = dir*(x - ilcn[j].x);
			dy = dir*(y - ilcn[j].y);
			dz = dir*(z - ilcn[j].z);
			switch (iOrder) {
			case 4:
				/*
				 ** Hexadecapole interaction.
				 */
				xxyy = ilcn[j].xxyy;
				xxzz = -ilcn[j].xxxx - xxyy;
				yyzz = xxyy - ilcn[j].yyyy;
				xyzz = -ilcn[j].xxxy - ilcn[j].xyyy;
				xzzz = -ilcn[j].xxxz - ilcn[j].xyyz;
				yzzz = -ilcn[j].xxyz - ilcn[j].yyyz;
				hirx = dx*(dx*(ilcn[j].xxxx*dx + 3*ilcn[j].xxxy*dy + 3*ilcn[j].xxxz*dz) +
						   3*xxyy*dy*dy + 3*xxzz*dz*dz + 6*ilcn[j].xxyz*dy*dz) +
							   dy*(3*xyzz*dz*dz + 3*ilcn[j].xyyz*dy*dz + ilcn[j].xyyy*dy*dy) +
								   xzzz*dz*dz*dz;
				hiry = dy*(dy*(3*ilcn[j].xyyy*dx + ilcn[j].yyyy*dy + 3*ilcn[j].yyyz*dz) +
						   3*xxyy*dx*dx + 3*yyzz*dz*dz + 6*ilcn[j].xyyz*dx*dz) +
							   dx*(3*xyzz*dz*dz + 3*ilcn[j].xxyz*dx*dz + ilcn[j].xxxy*dx*dx) +
								   yzzz*dz*dz*dz;
				hirz = dz*(dz*(3*xzzz*dx + 3*yzzz*dy - (xxzz+yyzz)*dz) +
						   3*xxzz*dx*dx + 3*yyzz*dy*dy + 6*xyzz*dx*dy) +
							   dx*(3*ilcn[j].xyyz*dy*dy + 3*ilcn[j].xxyz*dx*dy + ilcn[j].xxxy*dx*dx) +
								   ilcn[j].yyyz*dy*dy*dy;
				hir = hirx*dx + hiry*dy + hirz*dz;
				fPot -= (1.0/24.0)*dir5*(105*hir);
				ax -= (1.0/24.0)*dir6*(945*hir*dx - 420*hirx);
				ay -= (1.0/24.0)*dir6*(945*hir*dy - 420*hiry);
				az -= (1.0/24.0)*dir6*(945*hir*dz - 420*hirz);
			case 3:
				/*
				 ** Octopole interaction.
				 */
				xxy = ilcn[j].xxy;
				xyy = ilcn[j].xyy;
				xyz = ilcn[j].xyz;
				xzz = -ilcn[j].xxx - xyy;
				yzz = -xxy - ilcn[j].yyy;
				oirx = dx*(ilcn[j].xxx*dx + 2*xxy*dy + 2*ilcn[j].xxz*dz) +
					xyy*dy*dy + xzz*dz*dz + 2*xyz*dy*dz;
				oiry = dy*(2*xyy*dx + ilcn[j].yyy*dy + 2*ilcn[j].yyz*dz) +
					xxy*dx*dx + yzz*dz*dz + 2*xyz*dx*dz;
				oirz = dz*(2*xzz*dx + 2*yzz*dy - (ilcn[j].xxz + ilcn[j].yyz)*dz) +
					ilcn[j].xxz*dx*dx + ilcn[j].yyz*dy*dy + 2*xyz*dx*dy;
				oir = oirx*dx + oiry*dy + oirz*dz;
				fPot -= (1.0/6.0)*dir4*(15*oir);
				ax -= (1.0/6.0)*dir5*(105*oir*dx - 45*oirx);
				ay -= (1.0/6.0)*dir5*(105*oir*dy - 45*oiry);
				az -= (1.0/6.0)*dir5*(105*oir*dz - 45*oirz);
			case 2:
				/*
				 ** Quadrupole interaction.
				 */
				qirx = ilcn[j].xx*dx + ilcn[j].xy*dy + ilcn[j].xz*dz;
				qiry = ilcn[j].xy*dx + ilcn[j].yy*dy + ilcn[j].yz*dz;
				qirz = ilcn[j].xz*dx + ilcn[j].yz*dy - (ilcn[j].xx + ilcn[j].yy)*dz;
				qir = qirx*dx + qiry*dy + qirz*dz;
				fPot -= 0.5*dir3*(3*qir);
				ax -= 0.5*dir4*(15*qir*dx - 6*qirx);
				ay -= 0.5*dir4*(15*qir*dy - 6*qiry);
				az -= 0.5*dir4*(15*qir*dz - 6*qirz);
			case 1:
				/*
				 ** Monopole interaction.
				 */
				fPot -= dir*ilcn[j].m;
				ax -= dx*dir2*ilcn[j].m;
				ay -= dy*dir2*ilcn[j].m;
				az -= dz*dir2*ilcn[j].m;
				}
			}
		p[i].fPot = fPot;
		p[i].a[0] = ax;
		p[i].a[1] = ay;
		p[i].a[2] = az;
		/*
		 ** Try a cache check to improve responsiveness.
		 */
		mdlCacheCheck(pkd->mdl);
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
			SPLINE(1.0/sqrt(d2),d2,twoh,a,b,c,d);
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











