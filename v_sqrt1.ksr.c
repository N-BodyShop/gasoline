/*
 ** Vector 1/sqrt for the KSR computer, using the approximation instruction
 ** instead of the chebychev approximation as the first step.
 */
void v_sqrt1(int n,double *r2,double *a)
{
    int i;
    double s;
    
    i = 0;
    switch (n&7) {
    loop:
	s = _fsqrta(r2[i]);
	a[i++] = s*(1.5 - 0.5*r2[i]*s*s);
    case 7:
	s = _fsqrta(r2[i]);
	a[i++] = s*(1.5 - 0.5*r2[i]*s*s);
    case 6:
	s = _fsqrta(r2[i]);
	a[i++] = s*(1.5 - 0.5*r2[i]*s*s);
    case 5:
	s = _fsqrta(r2[i]);
	a[i++] = s*(1.5 - 0.5*r2[i]*s*s);
    case 4:
	s = _fsqrta(r2[i]);
	a[i++] = s*(1.5 - 0.5*r2[i]*s*s);
    case 3:
	s = _fsqrta(r2[i]);
	a[i++] = s*(1.5 - 0.5*r2[i]*s*s);
    case 2:
	s = _fsqrta(r2[i]);
	a[i++] = s*(1.5 - 0.5*r2[i]*s*s);
    case 1:
	s = _fsqrta(r2[i]);
	a[i++] = s*(1.5 - 0.5*r2[i]*s*s);
    case 0:
	if (i<n) goto loop;
	}
    }





