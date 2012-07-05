#ifndef QEVAL_HINCLUDED
#define QEVAL_HINCLUDED
/*
 ** QEVAL() Macro, reduced multipole supporting up to l=4
 **             Joachim Stadel (May 1995)
 ** Parameters to QEVAL MACRO are:
 **   iOrder:   (int) order of expansion.
 **   mom:      (struct ilCellNewt) contains the reduced multipoles.
 **   gam:      (double *) array of gamma numbers.
 **   dx,dy,dz: (double) displacements FROM the expansion center of mass.
 **   ax,ay,az: (double) acceleration (accumulated to)
 **   fPot:             (double) potential (accumulated to)
 ** Scoring: (+,*)
 **		Hexa=(34,123) Octu=(22,57) Quad=(13,18) Mono=(5,5) Total=(74,203) = 277 Flops/Hexa-Qeval
 ** Case 4: (74,203)	total = 277
 ** Case 3:	(40,80)		total = 120
 ** Case 2: (18,23)		total = 41
 ** Case 1: (5,5)		total = 10
 */
#define QEVAL(iOrder,mom,gam,dx,dy,dz,ax,ay,az,fPot)\
{\
double Qmirx,Qmiry,Qmirz,Qmir,Qta;\
Qta = 0.0;\
switch (iOrder) {\
case 4:\
Qmirx = (1.0/6.0)*(mom.xzzz*dz*dz*dz + 3*mom.xyzz*dy*dz*dz + 3*mom.xyyz*dy*dy*dz + mom.xyyy*dy*dy*dy + 3*mom.xxzz*dx*dz*dz + 6*mom.xxyz*dx*dy*dz + 3*mom.xxyy*dx*dy*dy + 3*mom.xxxz*dx*dx*dz + 3*mom.xxxy*dx*dx*dy + mom.xxxx*dx*dx*dx);\
Qmiry = (1.0/6.0)*(mom.yzzz*dz*dz*dz + 3*mom.xyzz*dx*dz*dz + 3*mom.xxyz*dx*dx*dz + mom.xxxy*dx*dx*dx + 3*mom.yyzz*dy*dz*dz + 6*mom.xyyz*dx*dy*dz + 3*mom.xxyy*dx*dx*dy + 3*mom.yyyz*dy*dy*dz + 3*mom.xyyy*dx*dy*dy + mom.yyyy*dy*dy*dy);\
Qmirz = (1.0/6.0)*(mom.yyyz*dy*dy*dy + 3*mom.xyyz*dx*dy*dy + 3*mom.xxyz*dx*dx*dy + mom.xxxz*dx*dx*dx + 3*mom.yyzz*dy*dy*dz + 6*mom.xyzz*dx*dy*dz + 3*mom.xxzz*dx*dx*dz + 3*mom.yzzz*dy*dz*dz + 3*mom.xzzz*dx*dz*dz + mom.zzzz*dz*dz*dz);\
Qmir = (1.0/4.0)*(Qmirx*dx + Qmiry*dy + Qmirz*dz);\
fPot -= gam[4]*Qmir;\
Qta += gam[5]*Qmir;\
ax += gam[4]*Qmirx;\
ay += gam[4]*Qmiry;\
az += gam[4]*Qmirz;\
case 3:\
Qmirx = (1.0/2.0)*(mom.xzz*dz*dz + 2*mom.xyz*dy*dz + mom.xyy*dy*dy + 2*mom.xxz*dx*dz + 2*mom.xxy*dx*dy + mom.xxx*dx*dx);\
Qmiry = (1.0/2.0)*(mom.yzz*dz*dz + 2*mom.xyz*dx*dz + mom.xxy*dx*dx + 2*mom.yyz*dy*dz + 2*mom.xyy*dx*dy + mom.yyy*dy*dy);\
Qmirz = (1.0/2.0)*(mom.yyz*dy*dy + 2*mom.xyz*dx*dy + mom.xxz*dx*dx + 2*mom.yzz*dy*dz + 2*mom.xzz*dx*dz + mom.zzz*dz*dz);\
Qmir = (1.0/3.0)*(Qmirx*dx + Qmiry*dy + Qmirz*dz);\
fPot -= gam[3]*Qmir;\
Qta += gam[4]*Qmir;\
ax += gam[3]*Qmirx;\
ay += gam[3]*Qmiry;\
az += gam[3]*Qmirz;\
case 2:\
Qmirx = (1.0/1.0)*(mom.xz*dz + mom.xy*dy + mom.xx*dx);\
Qmiry = (1.0/1.0)*(mom.yz*dz + mom.xy*dx + mom.yy*dy);\
Qmirz = (1.0/1.0)*(mom.yz*dy + mom.xz*dx + mom.zz*dz);\
Qmir = (1.0/2.0)*(Qmirx*dx + Qmiry*dy + Qmirz*dz);\
fPot -= gam[2]*Qmir;\
Qta += gam[3]*Qmir;\
ax += gam[2]*Qmirx;\
ay += gam[2]*Qmiry;\
az += gam[2]*Qmirz;\
case 1:\
default:\
fPot -= gam[0]*mom.m;\
Qta += gam[1]*mom.m;\
ax -= dx*Qta;\
ay -= dy*Qta;\
az -= dz*Qta;\
}\
}

#define QEVAL_FLOP	{10,10,41,120,277}

#endif
