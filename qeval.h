#ifndef QEVAL_HINCLUDED
#define QEVAL_HINCLUDED

/*
 ** QEVAL() Macro, reduced quadrupole evaluator currently supports up to Hexadecapole order.
 **		Joachim Stadel (May 1995)
 ** Parameters to QEVAL MACRO are:
 **   iOrder:   (int) order of expansion.
 **   mom:	(struct ilCellNewt) contains the reduced multipoles.
 **   gam:	(double *) array of gamma numbers.
 **   dx,dy,dz:	(double) displacements FROM the expansion center of mass.	
 **   ax,ay,az: (double) acceleration (accumulated to) 
 **   fPot:		(double) potential (accumulated to)
 */
#define QEVAL(iOrder,mom,gam,dx,dy,dz,ax,ay,az,fPot)\
{\
	double Qmirx,Qmiry,Qmirz,Qmir,Qta;\
	double Qxxyy,Qxxzz,Qyyzz,Qxyzz,Qxzzz,Qyzzz,Qxzz,Qyzz;\
	Qta = 0.0;\
	switch (iOrder) {\
	case 4:\
		/*\
		 ** Hexadecapole interaction.\
		 */\
		Qxxyy = mom.xxyy;\
		Qxxzz = -mom.xxxx - Qxxyy;\
		Qyyzz = -Qxxyy - mom.yyyy;\
		Qxyzz = -mom.xxxy - mom.xyyy;\
		Qxzzz = -mom.xxxz - mom.xyyz;\
		Qyzzz = -mom.xxyz - mom.yyyz;\
		Qmirx = (1.0/6.0)*(dx*(dx*(mom.xxxx*dx + 3*mom.xxxy*dy + 3*mom.xxxz*dz) +\
							   3*Qxxyy*dy*dy + 3*Qxxzz*dz*dz + 6*mom.xxyz*dy*dz) +\
						   dy*(3*Qxyzz*dz*dz + 3*mom.xyyz*dy*dz + mom.xyyy*dy*dy) +\
						   Qxzzz*dz*dz*dz);\
		Qmiry = (1.0/6.0)*(dy*(dy*(3*mom.xyyy*dx + mom.yyyy*dy + 3*mom.yyyz*dz) +\
							   3*Qxxyy*dx*dx + 3*Qyyzz*dz*dz + 6*mom.xyyz*dx*dz) +\
						   dx*(3*Qxyzz*dz*dz + 3*mom.xxyz*dx*dz + mom.xxxy*dx*dx) +\
						   Qyzzz*dz*dz*dz);\
		Qmirz = (1.0/6.0)*(dz*(dz*(3*Qxzzz*dx + 3*Qyzzz*dy - (Qxxzz + Qyyzz)*dz) +\
							   3*Qxxzz*dx*dx + 3*Qyyzz*dy*dy + 6*Qxyzz*dx*dy) +\
						   dx*(3*mom.xyyz*dy*dy + 3*mom.xxyz*dx*dy + mom.xxxz*dx*dx) +\
						   mom.yyyz*dy*dy*dy);\
		Qmir = 0.25*(Qmirx*dx + Qmiry*dy + Qmirz*dz);\
		fPot -= gam[4]*Qmir;\
		Qta += gam[5]*Qmir;\
		ax += gam[4]*Qmirx;\
		ay += gam[4]*Qmiry;\
		az += gam[4]*Qmirz;\
	case 3:\
		/*\
		 ** Octopole interaction.\
		 */\
		Qxzz = -mom.xxx - mom.xyy;\
		Qyzz = -mom.xxy - mom.yyy;\
		Qmirx = 0.5*(dx*(mom.xxx*dx + 2*mom.xxy*dy + 2*mom.xxz*dz) +\
					 mom.xyy*dy*dy + Qxzz*dz*dz + 2*mom.xyz*dy*dz);\
		Qmiry = 0.5*(dy*(2*mom.xyy*dx + mom.yyy*dy + 2*mom.yyz*dz) +\
					 mom.xxy*dx*dx + Qyzz*dz*dz + 2*mom.xyz*dx*dz);\
		Qmirz = 0.5*(dz*(2*Qxzz*dx + 2*Qyzz*dy - (mom.xxz + mom.yyz)*dz) +\
					 mom.xxz*dx*dx + mom.yyz*dy*dy + 2*mom.xyz*dx*dy);\
		Qmir = (1.0/3.0)*(Qmirx*dx + Qmiry*dy + Qmirz*dz);\
		fPot -= gam[3]*Qmir;\
		Qta += gam[4]*Qmir;\
		ax += gam[3]*Qmirx;\
		ay += gam[3]*Qmiry;\
		az += gam[3]*Qmirz;\
	case 2:\
		/*\
		 ** Quadrupole interaction.\
		 */\
		Qmirx = mom.xx*dx + mom.xy*dy + mom.xz*dz;\
		Qmiry = mom.xy*dx + mom.yy*dy + mom.yz*dz;\
		Qmirz = mom.xz*dx + mom.yz*dy - (mom.xx + mom.yy)*dz;\
		Qmir = 0.5*(Qmirx*dx + Qmiry*dy + Qmirz*dz);\
		fPot -= gam[2]*Qmir;\
		Qta += gam[3]*Qmir;\
		ax += gam[2]*Qmirx;\
		ay += gam[2]*Qmiry;\
		az += gam[2]*Qmirz;\
	case 1:\
	default:\
		/*\
		 ** Monopole interaction.\
		 */\
		fPot -= gam[0]*mom.m;\
		Qta += gam[1]*mom.m;\
		ax -= dx*Qta;\
		ay -= dy*Qta;\
		az -= dz*Qta;\
		}\
	}

#define QEVAL_H(iOrder,mom,gam,dx,dy,dz,fPotc,fPots)\
{\
	double Qmirx,Qmiry,Qmirz,Qmir;\
	double Qxxyy,Qxxzz,Qyyzz,Qxyzz,Qxzzz,Qyzzz,Qxzz,Qyzz;\
	switch (iOrder) {\
	case 4:\
		/*\
		 ** Hexadecapole interaction.\
		 */\
		Qxxyy = mom.xxyy;\
		Qxxzz = -mom.xxxx - Qxxyy;\
		Qyyzz = -Qxxyy - mom.yyyy;\
		Qxyzz = -mom.xxxy - mom.xyyy;\
		Qxzzz = -mom.xxxz - mom.xyyz;\
		Qyzzz = -mom.xxyz - mom.yyyz;\
		Qmirx = (1.0/6.0)*(dx*(dx*(mom.xxxx*dx + 3*mom.xxxy*dy + 3*mom.xxxz*dz) +\
							   3*Qxxyy*dy*dy + 3*Qxxzz*dz*dz + 6*mom.xxyz*dy*dz) +\
						   dy*(3*Qxyzz*dz*dz + 3*mom.xyyz*dy*dz + mom.xyyy*dy*dy) +\
						   Qxzzz*dz*dz*dz);\
		Qmiry = (1.0/6.0)*(dy*(dy*(3*mom.xyyy*dx + mom.yyyy*dy + 3*mom.yyyz*dz) +\
							   3*Qxxyy*dx*dx + 3*Qyyzz*dz*dz + 6*mom.xyyz*dx*dz) +\
						   dx*(3*Qxyzz*dz*dz + 3*mom.xxyz*dx*dz + mom.xxxy*dx*dx) +\
						   Qyzzz*dz*dz*dz);\
		Qmirz = (1.0/6.0)*(dz*(dz*(3*Qxzzz*dx + 3*Qyzzz*dy - (Qxxzz + Qyyzz)*dz) +\
							   3*Qxxzz*dx*dx + 3*Qyyzz*dy*dy + 6*Qxyzz*dx*dy) +\
						   dx*(3*mom.xyyz*dy*dy + 3*mom.xxyz*dx*dy + mom.xxxz*dx*dx) +\
						   mom.yyyz*dy*dy*dy);\
		Qmir = 0.25*(Qmirx*dx + Qmiry*dy + Qmirz*dz);\
		fPotc -= gam[4]*Qmir;\
	case 3:\
		/*\
		 ** Octopole interaction.\
		 */\
		Qxzz = -mom.xxx - mom.xyy;\
		Qyzz = -mom.xxy - mom.yyy;\
		Qmirx = 0.5*(dx*(mom.xxx*dx + 2*mom.xxy*dy + 2*mom.xxz*dz) +\
					 mom.xyy*dy*dy + Qxzz*dz*dz + 2*mom.xyz*dy*dz);\
		Qmiry = 0.5*(dy*(2*mom.xyy*dx + mom.yyy*dy + 2*mom.yyz*dz) +\
					 mom.xxy*dx*dx + Qyzz*dz*dz + 2*mom.xyz*dx*dz);\
		Qmirz = 0.5*(dz*(2*Qxzz*dx + 2*Qyzz*dy - (mom.xxz + mom.yyz)*dz) +\
					 mom.xxz*dx*dx + mom.yyz*dy*dy + 2*mom.xyz*dx*dy);\
		Qmir = (1.0/3.0)*(Qmirx*dx + Qmiry*dy + Qmirz*dz);\
		fPots += gam[3]*Qmir;\
	case 2:\
		/*\
		 ** Quadrupole interaction.\
		 */\
		Qmirx = mom.xx*dx + mom.xy*dy + mom.xz*dz;\
		Qmiry = mom.xy*dx + mom.yy*dy + mom.yz*dz;\
		Qmirz = mom.xz*dx + mom.yz*dy - (mom.xx + mom.yy)*dz;\
		Qmir = 0.5*(Qmirx*dx + Qmiry*dy + Qmirz*dz);\
		fPotc -= gam[2]*Qmir;\
	case 1:\
	default:\
		/*\
		 ** Monopole interaction.\
		 */\
		fPotc -= gam[0]*mom.m;\
		}\
	}
#endif





