#define KMAXX 7 
#define IMAXX (KMAXX+1)

typedef struct StiffContextStructure {
  double **d, *x;
  double eps,eps1;
  double     hMin; /* Trying to keep the integrator sane */

  /* stifbs */
  double *dfdx,**dfdy,*err,*yerr,*ysav,*yseq;

  int first,kmax,kopt;
  double xnew;
  double a[IMAXX+1];
  double alf[KMAXX+1][KMAXX+1];
  int nseq[IMAXX+1];

  /* simpr */
  int *indx;
  double **simpr_a, *del, *ytemp;

  /* pzextr */
  double *c;
  int nv;

  void *Data; /* Pointer to fixed data used by derivs/jacobn */
  void (*derivs)(void *Data, double, double [], double []);
  void (*jacobn)(void *Data, double x, double y[], double dfdx[], double **dfdy);
} STIFF;

/*
 * Integerator Step Headers
 */

STIFF *StiffInit( double eps, int nv, void *Data, 
		   void (*derivs)(void *Data, double, double [], double []),
		   void (*jacobn)(void *Data, double x, double y[], double dfdx[], double **dfdy) );
		   
void StiffFinalize( STIFF *s );
void StiffStep(STIFF *s, double y[], double dydx[], double *xx, double htry, 
		double yscal[], double *hdid, double *hnext );

void simpr(STIFF *s, double y[], double dydx[], double dfdx[], double **dfdy,
	   double xs, double htot, int nstep, double yout[] );
void pzextr(STIFF *s, int iest, double xest, double yest[], double yz[], double dy[] );

void lubksb(double **a, int n, int *indx, double b[]);
void ludcmp(double **a, int n, int *indx, double *d);

/* 
 * Root Finder Header
 */

double RootFind(double (*func)(void *Data, double), void *Data, double x1, double x2, double tol);

/*
 * Utils
 */

static inline double SQR(double a) 
{
    return a*a;
}

static inline double FMAX(double maxarg1, double maxarg2)
{
    return (maxarg1 > maxarg2 ? maxarg1 : maxarg2);
}

static inline double FMIN(double minarg1, double minarg2)
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


