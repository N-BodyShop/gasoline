#ifndef MEVAL_HINCLUDED
#define MEVAL_HINCLUDED
/*
 ** MEVAL() Macro, complete multipole supporting up to l=4
 **             Joachim Stadel (May 1995)
 ** Parameters to MEVAL MACRO are:
 **   iOrder:   (int) order of expansion.
 **   mom:      (struct ilCellNewt) contains the COMPLETE multipoles.
 **   gam:      (double *) array of gamma numbers.
 **   dx,dy,dz: (double) displacements FROM the expansion center of mass.
 **   ax,ay,az: (double) acceleration (accumulated to)
 **   fPot:             (double) potential (accumulated to)
 **
 ** Quick scoring guide: (+,*)
 ** Case 4:	(102,241) 	total = 343
 ** Case 3:	(57,94)  	total = 151
 ** Case 2:	(22,26)		total = 48
 ** Case 1:	(5,5)		total = 10
 */
#define MEVAL(iOrder,mom,gam,dx,dy,dz,ax,ay,az,fPot)\
{\
double Qmirx,Qmiry,Qmirz,Qmir,Qta;\
double Qxx, Qyy, Qzz, Qxy, Qxz, Qyz, Qx, Qy, Qz, Qtr, Qhx, Qhy, Qhz, Qh;\
Qta = 0.0;\
switch (iOrder) {\
case 4:\
Qmirx = (1.0/6.0)*(mom.xzzz*dz*dz*dz + 3*mom.xyzz*dy*dz*dz + 3*mom.xyyz*dy*dy*dz + mom.xyyy*dy*dy*dy + 3*mom.xxzz*dx*dz*dz + 6*mom.xxyz*dx*dy*dz + 3*mom.xxyy*dx*dy*dy + 3*mom.xxxz*dx*dx*dz + 3*mom.xxxy*dx*dx*dy + mom.xxxx*dx*dx*dx);\
Qmiry = (1.0/6.0)*(mom.yzzz*dz*dz*dz + 3*mom.xyzz*dx*dz*dz + 3*mom.xxyz*dx*dx*dz + mom.xxxy*dx*dx*dx + 3*mom.yyzz*dy*dz*dz + 6*mom.xyyz*dx*dy*dz + 3*mom.xxyy*dx*dx*dy + 3*mom.yyyz*dy*dy*dz + 3*mom.xyyy*dx*dy*dy + mom.yyyy*dy*dy*dy);\
Qmirz = (1.0/6.0)*(mom.yyyz*dy*dy*dy + 3*mom.xyyz*dx*dy*dy + 3*mom.xxyz*dx*dx*dy + mom.xxxz*dx*dx*dx + 3*mom.yyzz*dy*dy*dz + 6*mom.xyzz*dx*dy*dz + 3*mom.xxzz*dx*dx*dz + 3*mom.yzzz*dy*dz*dz + 3*mom.xzzz*dx*dz*dz + mom.zzzz*dz*dz*dz);\
Qmir = (1.0/4.0)*(Qmirx*dx + Qmiry*dy + Qmirz*dz);\
Qxx = mom.xxxx + mom.xxyy + mom.xxzz;\
Qxy = mom.xxxy + mom.xyyy + mom.xyzz;\
Qxz = mom.xxxz + mom.xyyz + mom.xzzz;\
Qyy = mom.xxyy + mom.yyyy + mom.yyzz;\
Qyz = mom.xxyz + mom.yyyz + mom.yzzz;\
Qzz = mom.xxzz + mom.yyzz + mom.zzzz;\
Qtr = (1.0/8.0)*(Qxx + Qyy + Qzz);\
Qhx = 0.5*(Qxx*dx + Qxy*dy + Qxz*dz);\
Qhy = 0.5*(Qxy*dx + Qyy*dy + Qyz*dz);\
Qhz = 0.5*(Qxz*dx + Qyz*dy + Qzz*dz);\
Qh = 0.5*(Qhx*dx + Qhy*dy + Qhz*dz);\
fPot -= gam[4]*Qmir - gam[3]*Qh + gam[2]*Qtr;\
Qta += gam[5]*Qmir - gam[4]*Qh + gam[3]*Qtr;\
ax += gam[4]*Qmirx - gam[3]*Qhx;\
ay += gam[4]*Qmiry - gam[3]*Qhy;\
az += gam[4]*Qmirz - gam[3]*Qhz;\
case 3:\
Qmirx = (1.0/2.0)*(mom.xzz*dz*dz + 2*mom.xyz*dy*dz + mom.xyy*dy*dy + 2*mom.xxz*dx*dz + 2*mom.xxy*dx*dy + mom.xxx*dx*dx);\
Qmiry = (1.0/2.0)*(mom.yzz*dz*dz + 2*mom.xyz*dx*dz + mom.xxy*dx*dx + 2*mom.yyz*dy*dz + 2*mom.xyy*dx*dy + mom.yyy*dy*dy);\
Qmirz = (1.0/2.0)*(mom.yyz*dy*dy + 2*mom.xyz*dx*dy + mom.xxz*dx*dx + 2*mom.yzz*dy*dz + 2*mom.xzz*dx*dz + mom.zzz*dz*dz);\
Qmir = (1.0/3.0)*(Qmirx*dx + Qmiry*dy + Qmirz*dz);\
Qx = 0.5*(mom.xxx + mom.xyy + mom.xzz);\
Qy = 0.5*(mom.xxy + mom.yyy + mom.yzz);\
Qz = 0.5*(mom.xxz + mom.yyz + mom.zzz);\
Qtr = Qx*dx + Qy*dy + Qz*dz;\
fPot -= gam[3]*Qmir - gam[2]*Qtr;\
Qta += gam[4]*Qmir - gam[3]*Qtr;\
ax += gam[3]*Qmirx - gam[2]*Qx;\
ay += gam[3]*Qmiry - gam[2]*Qy;\
az += gam[3]*Qmirz - gam[2]*Qz;\
case 2:\
Qmirx = (1.0/1.0)*(mom.xz*dz + mom.xy*dy + mom.xx*dx);\
Qmiry = (1.0/1.0)*(mom.yz*dz + mom.xy*dx + mom.yy*dy);\
Qmirz = (1.0/1.0)*(mom.yz*dy + mom.xz*dx + mom.zz*dz);\
Qmir = (1.0/2.0)*(Qmirx*dx + Qmiry*dy + Qmirz*dz);\
Qtr = 0.5*(mom.xx + mom.yy + mom.zz);\
fPot -= gam[2]*Qmir - gam[1]*Qtr;\
Qta += gam[3]*Qmir - gam[2]*Qtr;\
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

#define MEVAL_FLOP	{10,10,48,151,343}

#endif
