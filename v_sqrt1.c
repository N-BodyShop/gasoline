#include <math.h>

#define NCHEB		200

/*
 ** Function shifted from 1 to 2 to -1 to 1 (for 1/sqrt).
 */
#define FUNC(r) 	pow(1.5 + 0.5*r,-0.5)

void v_sqrt1(int n,double *r2,double *a)
{
	static int bInit = 0;
	static double u[2048],t[2048];
	double g0,g1,g2,g3,g4,g5;
	static double sg0,sg1,sg2,sg3,sg4,sg5;
	static int sbLittleEndian;
	int bLittleEndian;
	int i,it;
	double x,s;

	if (!bInit) {
		double xi,pi,zero,factor,sum;		
		double c[NCHEB],f[NCHEB+1],z[NCHEB+1];
		double d0,d1,d2,d3,d4,d5;
		int j,k;
		union
			{
			short s;
			char c[sizeof (short)];
			} ue;
		/*
		 ** Determine the Endian'ness of the architecture.
		 */
		ue.s = 1;
		sbLittleEndian = (ue.c[sizeof(short) - 1] == 1)?0:1;
		/*
		 ** Initialize exponent tables.
		 */
		xi = 1.0;
		t[1023] = 1.0;
		u[1023] = 1.0;
		for (i=1;i<1024;++i) {
			xi *= 0.5;
			t[1023+i] = sqrt(xi);
			if(xi != 0.0)
			t[1023-i] = 1.0/t[1023+i];
			u[1023+i] = xi;
			if(xi != 0.0)
			u[1023-i] = 1.0/xi;
			}
		/*
		 ** Precompute zeros of Chebychev polynomials and function.
		 */
		pi = M_PI;
		for (k=1;k<=NCHEB;++k) {
			z[k] = (pi*(k - 0.5))/NCHEB;
			zero = cos(z[k]);
			f[k] = FUNC(zero);
			}
		/*
		 ** Get coefficients of Chebychev fit
		 */
		factor = 2.0/NCHEB;
		for (j=0;j<NCHEB;++j) {
			sum = 0.0;
			for (k=1;k<=NCHEB;++k) sum += f[k]*cos(z[k]*j);
			c[j] = factor*sum;
			}
		/*
		 ** Get coefficients of powers of x on -1 < x < 1
		 */
		d0 = 0.5*c[0] - c[2] + c[4];
		d1 = c[1] - 3*c[3] + 5*c[5];
		d2 = 2*c[2] - 8*c[4];
		d3 = 4*c[3] - 20*c[5];
		d4 = 8*c[4];
		d5 = 16*c[5];
		/*
		 ** Get coefficients of powers of x on 1 < x < 2.
		 */
		sg0 = d0 - 3*d1 +  9*d2 - 27*d3 +  81*d4 -  243*d5;
		sg1 =      2*d1 - 12*d2 + 54*d3 - 216*d4 +  810*d5;
		sg2 =              4*d2 - 36*d3 + 216*d4 - 1080*d5;
		sg3 =                      8*d3 -  96*d4 +  720*d5;
		sg4 =                              16*d4 -  240*d5;
		sg5 =                                        32*d5;
		/*
		 ** Flag that the static variables are initialized.
		 */
		bInit = 1;
		}
	/*
	 * Load coefficients.
	 */
	g0 = sg0;
	g1 = sg1;
	g2 = sg2;
	g3 = sg3;
	g4 = sg4;
	g5 = sg5;
	bLittleEndian = sbLittleEndian;
#ifndef UNROLLED_VSQRT
	/*
	 ** Actual Calculation. Rolled loop.
	 */
	for(i = 0; i < n; i++) {
		it = (((int *)r2)[(i<<1)|bLittleEndian])>>20;
		x = r2[i]*u[it];
		s = g0 + x*(g1 + x*(g2 + x*(g3 + x*(g4 + x*g5))));
		s *= 1.5 - 0.5*x*s*s;
		a[i] = s*t[it];
		}
#else
	/*
	 ** Actual Calculation. Loop unrooled 8 times.
	 */
	i = 0;
	switch (n&7) {
	loop:
		it = (((int *)r2)[(i<<1)|bLittleEndian])>>20;
		x = r2[i]*u[it];
		s = g0 + x*(g1 + x*(g2 + x*(g3 + x*(g4 + x*g5))));
		s *= 1.5 - 0.5*x*s*s;
		a[i++] = s*t[it];
	case 7:
		it = (((int *)r2)[(i<<1)|bLittleEndian])>>20;
		x = r2[i]*u[it];
		s = g0 + x*(g1 + x*(g2 + x*(g3 + x*(g4 + x*g5))));
		s *= 1.5 - 0.5*x*s*s;
		a[i++] = s*t[it];
	case 6:
		it = (((int *)r2)[(i<<1)|bLittleEndian])>>20;
		x = r2[i]*u[it];
		s = g0 + x*(g1 + x*(g2 + x*(g3 + x*(g4 + x*g5))));
		s *= 1.5 - 0.5*x*s*s;
		a[i++] = s*t[it];
	case 5:
		it = (((int *)r2)[(i<<1)|bLittleEndian])>>20;
		x = r2[i]*u[it];
		s = g0 + x*(g1 + x*(g2 + x*(g3 + x*(g4 + x*g5))));
		s *= 1.5 - 0.5*x*s*s;
		a[i++] = s*t[it];
	case 4:
		it = (((int *)r2)[(i<<1)|bLittleEndian])>>20;
		x = r2[i]*u[it];
		s = g0 + x*(g1 + x*(g2 + x*(g3 + x*(g4 + x*g5))));
		s *= 1.5 - 0.5*x*s*s;
		a[i++] = s*t[it];
	case 3:
		it = (((int *)r2)[(i<<1)|bLittleEndian])>>20;
		x = r2[i]*u[it];
		s = g0 + x*(g1 + x*(g2 + x*(g3 + x*(g4 + x*g5))));
		s *= 1.5 - 0.5*x*s*s;
		a[i++] = s*t[it];
	case 2:
		it = (((int *)r2)[(i<<1)|bLittleEndian])>>20;
		x = r2[i]*u[it];
		s = g0 + x*(g1 + x*(g2 + x*(g3 + x*(g4 + x*g5))));
		s *= 1.5 - 0.5*x*s*s;
		a[i++] = s*t[it];
	case 1:
		it = (((int *)r2)[(i<<1)|bLittleEndian])>>20;
		x = r2[i]*u[it];
		s = g0 + x*(g1 + x*(g2 + x*(g3 + x*(g4 + x*g5))));
		s *= 1.5 - 0.5*x*s*s;
		a[i++] = s*t[it];
	case 0:
		if (i<n) goto loop;
		}
#endif
	}





