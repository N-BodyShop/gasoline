#include <math.h>

void v_sqrt1(double *pd2,double *pdir,int n)
{
	int i;
	
	for (i=0;i<n;++i) pdir[i] = 1.0/sqrt(pd2[i]);
	}

