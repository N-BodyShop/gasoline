#ifndef GRAV_HINCLUDED
#define GRAV_HINCLUDED


/*
 ** see (A1) and (A2) of TREESPH: A UNIFICATION OF SPH WITH THE 
 ** HIERARCHICAL TREE METHOD by Lars Hernquist and Neal Katz.
 ** APJ Supplemant Series 70:416-446, 1989
 */
#define SPLINE(r2,twoh,dir,dir3)\
{\
	double r,u,dih;\
	r = sqrt(r2);\
	if (r < twoh) {\
		dih = 2.0/twoh;\
		u = r*dih;\
		if (u < 1.0) {\
			dir = -2.0*dih*u*u*(1.0/3.0 - u*u*(3.0/20.0 - 1.0/20.0*u)) + 7.0/5.0*dih;\
			dir3 = dih*dih*dih*(4.0/3.0 - u*u*(6.0/5.0 - 1.0/2.0*u));\
			}\
		else {\
			dir = 1.0/r;\
			dir3 = dir*dir*dir;\
			dir = -1.0/15.0*dir - dih*u*u*(4.0/3.0 - u*(1.0 - u*(3.0/10.0 - 1.0/30.0*u))) + 8.0/5.0*dih;\
			dir3 *= (-1.0/15.0 + u*u*u*(8.0/3.0 - u*(3.0 - u*(6.0/5.0 - 1.0/6.0*u))));\
			}\
		}\
	else {\
		dir = 1.0/r;\
		dir3 = dir*dir*dir;\
		}\
	}


#if (0)
#define SPLINE(d2,twoh,dir,dir3)\
{\
	dir = 1.0/sqrt(d2);\
	dir3 = dir*dir*dir;\
	}
#endif


void pkdBucketInteract(PKD,int);

#endif
