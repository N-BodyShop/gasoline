/** @brief Context for stiff integration
 */
typedef struct StiffContextStructure {
    double epsmin;		/* relative accuracy criterion */
    double sqreps;		/* parameter to calculate initial timestep */
    double epscl;		/* 1/epsmin */
    double epsmax;		/* repeat timestep if correction is
				   greater than epsmax*epsmin*y(i) */   
    double dtmin;		/* minimum timestep allowed */
    int itermax;		/* number of corrector iterations */
    int nv; 			/* number of dependent variables */
    double *ymin;		/* minimum value for y */
    double *y0;			/* initial y for global timestep */
    double *y1;			/* predicted value */
    double *q;			/* scratch for creation rate */
    double *d;			/* scratch for destruction rate */
    double *rtau;		/* ratio of timestep to timescale */
    double *ys;			/* initial y for individual timesteps */
    double *qs;			/* initial production rate */
    double *rtaus;		/* initial rtau */
    double *scrarray;
    void *Data; /* Pointer to fixed data used by derivs */
    void (*derivs)(void *Data, double, const double [], double [], double []);
} STIFF;


/*
 * Integrator Step Headers
 */

STIFF *StiffInit(double eps, 	/* relative accuracy parameter */
		 int nv,	/* number of dependent variables */
		 void *Data, 	/* pointer to extra data */
		 void (*derivs)(void *Data,
				double t, const double yin[],  /* input */
				double yheat[],	/* heating or creation
						   rate */
			       double ycool[]	/* cooling or
						   destruction rate */
				)
		 );
		   
void StiffFinalize( STIFF *s );
void StiffStep(STIFF *s, double y[], double tstart, double htry) ;

/* 
 * Root Finder Header
 */

double RootFind(double (*func)(void *Data, double), void *Data, double x1, double x2, double tol);

/*
 * Utils
 */

static double SQR(double a) 
{
    return a*a;
}

static double FMAX(double maxarg1, double maxarg2)
{
    return (maxarg1 > maxarg2 ? maxarg1 : maxarg2);
}

static double FMIN(double minarg1, double minarg2)
{
    return (minarg1 < minarg2 ? minarg1 : minarg2);
}

static double MAX(double maxarg1, double maxarg2)
{
    return (maxarg1 > maxarg2 ? maxarg1 : maxarg2);
}

static double MIN(double minarg1, double minarg2)
{
    return (minarg1 < minarg2 ? minarg1 : minarg2);
}
#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))

void nrerror(char error_text[]);
double *vector(long nl, long nh);
int *ivector(long nl, long nh);
unsigned char *cvector(long nl, long nh);
unsigned long *lvector(long nl, long nh);
double *dvector(long nl, long nh);
double **matrix(long nrl, long nrh, long ncl, long nch);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
int **imatrix(long nrl, long nrh, long ncl, long nch);
double **submatrix(double **a, long oldrl, long oldrh, long oldcl, long oldch,
  long newrl, long newcl);
double **convert_matrix(double *a, long nrl, long nrh, long ncl, long nch);
double ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_vector(double *v, long nl, long nh);
void free_ivector(int *v, long nl, long nh);
void free_cvector(unsigned char *v, long nl, long nh);
void free_lvector(unsigned long *v, long nl, long nh);
void free_dvector(double *v, long nl, long nh);
void free_matrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
void free_submatrix(double **b, long nrl, long nrh, long ncl, long nch);
void free_convert_matrix(double **b, long nrl, long nrh, long ncl, long nch);
void free_f3tensor(double ***t, long nrl, long nrh, long ncl, long nch,
  long ndl, long ndh);


