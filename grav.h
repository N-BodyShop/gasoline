#ifndef GRAV_HINCLUDED
#define GRAV_HINCLUDED


/*
 ** see (A1) and (A2) of TREESPH: A UNIFICATION OF SPH WITH THE 
 ** HIERARCHICAL TREE METHOD by Lars Hernquist and Neal Katz.
 ** APJ Supplemant Series 70:416-446, 1989
 ** 
 ** Higher derivative terms c and d for use with quadrupole spline
 ** softening (Joachim Stadel, Dec. 94).
 ** Had to append SPLINE_ to the local variables to avoid the case that one
 ** of the passed variables to the macro has the same name. It looks pretty
 ** unreadable now, oh well. (Joachim Stadel, Jan. 95)
 */
#define SPLINEQ(invr,r2,twoh,a,b,c,d)\
{\
	double SPLINE_u,SPLINE_dih,SPLINE_dir=(invr);\
	if ((r2) < (twoh)*(twoh)) {\
		SPLINE_dih = 2.0/(twoh);\
		SPLINE_u = SPLINE_dih/SPLINE_dir;\
		if (SPLINE_u < 1.0) {\
			a = SPLINE_dih*(7.0/5.0 - 2.0/3.0*SPLINE_u*SPLINE_u + 3.0/10.0*SPLINE_u*SPLINE_u*SPLINE_u*SPLINE_u\
					 - 1.0/10.0*SPLINE_u*SPLINE_u*SPLINE_u*SPLINE_u*SPLINE_u);\
			b = SPLINE_dih*SPLINE_dih*SPLINE_dih*(4.0/3.0 - 6.0/5.0*SPLINE_u*SPLINE_u + 1.0/2.0*SPLINE_u*SPLINE_u*SPLINE_u);\
		    c = SPLINE_dih*SPLINE_dih*SPLINE_dih*SPLINE_dih*SPLINE_dih*(12.0/5.0 - 3.0/2.0*SPLINE_u);\
			d = 3.0/2.0*SPLINE_dih*SPLINE_dih*SPLINE_dih*SPLINE_dih*SPLINE_dih*SPLINE_dih*SPLINE_dir;\
			}\
		else {\
			a = -1.0/15.0*SPLINE_dir + SPLINE_dih*(8.0/5.0 - 4.0/3.0*SPLINE_u*SPLINE_u + SPLINE_u*SPLINE_u*SPLINE_u\
			              - 3.0/10.0*SPLINE_u*SPLINE_u*SPLINE_u*SPLINE_u + 1.0/30.0*SPLINE_u*SPLINE_u*SPLINE_u*SPLINE_u*SPLINE_u);\
			b = -1.0/15.0*SPLINE_dir*SPLINE_dir*SPLINE_dir + SPLINE_dih*SPLINE_dih*SPLINE_dih*(8.0/3.0 - 3.0*SPLINE_u + 6.0/5.0*SPLINE_u*SPLINE_u - 1.0/6.0*SPLINE_u*SPLINE_u*SPLINE_u);\
			c = -1.0/5.0*SPLINE_dir*SPLINE_dir*SPLINE_dir*SPLINE_dir*SPLINE_dir + 3.0*SPLINE_dih*SPLINE_dih*SPLINE_dih*SPLINE_dih*SPLINE_dir\
				+ SPLINE_dih*SPLINE_dih*SPLINE_dih*SPLINE_dih*SPLINE_dih*(-12.0/5.0 + 1.0/2.0*SPLINE_u);\
			d = -SPLINE_dir*SPLINE_dir*SPLINE_dir*SPLINE_dir*SPLINE_dir*SPLINE_dir*SPLINE_dir\
				+ 3.0*SPLINE_dih*SPLINE_dih*SPLINE_dih*SPLINE_dih*SPLINE_dir*SPLINE_dir*SPLINE_dir\
					- 1.0/2.0*SPLINE_dih*SPLINE_dih*SPLINE_dih*SPLINE_dih*SPLINE_dih*SPLINE_dih*SPLINE_dir;\
			}\
		}\
	else {\
		a = SPLINE_dir;\
		b = a*a*a;\
		c = 3.0*b*a*a;\
		d = 5.0*c*a*a;\
		}\
	}


#define SPLINEM(invr,r2,twoh,a,b)\
{\
	double SPLINE_u,SPLINE_dih,SPLINE_dir=(invr);\
	if ((r2) < (twoh)*(twoh)) {\
		SPLINE_dih = 2.0/(twoh);\
		SPLINE_u = SPLINE_dih/SPLINE_dir;\
		if (SPLINE_u < 1.0) {\
			a = SPLINE_dih*(7.0/5.0 - 2.0/3.0*SPLINE_u*SPLINE_u + 3.0/10.0*SPLINE_u*SPLINE_u*SPLINE_u*SPLINE_u\
					 - 1.0/10.0*SPLINE_u*SPLINE_u*SPLINE_u*SPLINE_u*SPLINE_u);\
			b = SPLINE_dih*SPLINE_dih*SPLINE_dih*(4.0/3.0 - 6.0/5.0*SPLINE_u*SPLINE_u + 1.0/2.0*SPLINE_u*SPLINE_u*SPLINE_u);\
			}\
		else {\
			a = -1.0/15.0*SPLINE_dir + SPLINE_dih*(8.0/5.0 - 4.0/3.0*SPLINE_u*SPLINE_u + SPLINE_u*SPLINE_u*SPLINE_u\
			              - 3.0/10.0*SPLINE_u*SPLINE_u*SPLINE_u*SPLINE_u + 1.0/30.0*SPLINE_u*SPLINE_u*SPLINE_u*SPLINE_u*SPLINE_u);\
			b = -1.0/15.0*SPLINE_dir*SPLINE_dir*SPLINE_dir + SPLINE_dih*SPLINE_dih*SPLINE_dih*(8.0/3.0 - 3.0*SPLINE_u + 6.0/5.0*SPLINE_u*SPLINE_u - 1.0/6.0*SPLINE_u*SPLINE_u*SPLINE_u);\
			}\
		}\
	else {\
		a = SPLINE_dir;\
		b = a*a*a;\
		}\
	}


#define SPLINE(r2,twoh,a,b)\
{\
	double SPLINE_r,SPLINE_u,SPLINE_dih,SPLINE_dir;\
	SPLINE_r = sqrt(r2);\
	if (SPLINE_r < (twoh)) {\
		SPLINE_dih = 2.0/(twoh);\
		SPLINE_u = SPLINE_r*SPLINE_dih;\
		if (SPLINE_u < 1.0) {\
			a = SPLINE_dih*(7.0/5.0 - 2.0/3.0*SPLINE_u*SPLINE_u + 3.0/10.0*SPLINE_u*SPLINE_u*SPLINE_u*SPLINE_u\
					 - 1.0/10.0*SPLINE_u*SPLINE_u*SPLINE_u*SPLINE_u*SPLINE_u);\
			b = SPLINE_dih*SPLINE_dih*SPLINE_dih*(4.0/3.0 - 6.0/5.0*SPLINE_u*SPLINE_u + 1.0/2.0*SPLINE_u*SPLINE_u*SPLINE_u);\
			}\
		else {\
			SPLINE_dir = 1.0/SPLINE_r;\
			a = -1.0/15.0*SPLINE_dir + SPLINE_dih*(8.0/5.0 - 4.0/3.0*SPLINE_u*SPLINE_u + SPLINE_u*SPLINE_u*SPLINE_u\
			              - 3.0/10.0*SPLINE_u*SPLINE_u*SPLINE_u*SPLINE_u + 1.0/30.0*SPLINE_u*SPLINE_u*SPLINE_u*SPLINE_u*SPLINE_u);\
			b = -1.0/15.0*SPLINE_dir*SPLINE_dir*SPLINE_dir + SPLINE_dih*SPLINE_dih*SPLINE_dih*(8.0/3.0 - 3.0*SPLINE_u + 6.0/5.0*SPLINE_u*SPLINE_u - 1.0/6.0*SPLINE_u*SPLINE_u*SPLINE_u);\
			}\
		}\
	else {\
		a = 1.0/SPLINE_r;\
		b = a*a*a;\
		}\
	}


int pkdBucketInteract(PKD,int,int);

#endif
