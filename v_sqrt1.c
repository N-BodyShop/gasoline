#include <math.h>

void v_sqrt1(int n,double *pd2,double *pdir)
{
	int i;
	
	for (i=0;i<n;++i) pdir[i] = 1.0/sqrt(pd2[i]);
	}

