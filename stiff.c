#ifdef GASOLINE
#ifndef NOCOOLING

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "stiff.h"

#define NR_END 1
#define FREE_ARG char*



/* 
 * Stiff Bulirsch-Stoer
 */ 

#define SAFE1 0.25
#define SAFE2 0.7
#define REDMAX 1.0e-5
#define REDMIN 0.7
#define TINY 1.0e-30
#define SCALMX 0.1

STIFF *StiffInit( double eps, int nv, void *Data,
		   void (*derivs)(void *Data, double, double [], double []),
		   void (*jacobn)(void *Data, double x, double y[], double dfdx[], double **dfdy) 
		   ) 
{
  int k,iq;

  STIFF *s;

  s = (STIFF *) malloc(sizeof(STIFF));
  assert(s!=NULL);

  s->eps = eps;
  s->eps1=SAFE1*eps;

  s->nv = nv;  

  /* stifbs */
  s->d=matrix(1,KMAXX,1,KMAXX);
  s->dfdx=vector(1,nv);
  s->dfdy=matrix(1,nv,1,nv);
  s->err=vector(1,KMAXX);
  s->x=vector(1,KMAXX);
  s->yerr=vector(1,nv);
  s->ysav=vector(1,nv);
  s->yseq=vector(1,nv);

  /* stifbs setup */
  s->first = 1;
  s->kmax;
  s->kopt;
  s->xnew = -1.0e29;
  assert(IMAXX >= 8);
  s->nseq[0] = 0;
  s->nseq[1] = 2;
  s->nseq[2] = 6;
  s->nseq[3] = 10;
  s->nseq[4] = 14;
  s->nseq[5] = 22;
  s->nseq[6] = 34;
  s->nseq[7] = 50;
  s->nseq[8] = 70;
  s->a[1]=s->nseq[1]+1;
  for (k=1;k<=KMAXX;k++) s->a[k+1]=s->a[k]+s->nseq[k+1];
  for (iq=2;iq<=KMAXX;iq++) {
    for (k=1;k<iq;k++)
      s->alf[k][iq]=pow(s->eps1,((s->a[k+1]-s->a[iq+1])/
			      ((s->a[iq+1]-s->a[1]+1.0)*(2*k+1))));
  }
  s->a[1] += nv;
  for (k=1;k<=KMAXX;k++) s->a[k+1]=s->a[k]+s->nseq[k+1];
  for (s->kopt=2;s->kopt<KMAXX;s->kopt++)
    if (s->a[s->kopt+1] > s->a[s->kopt]*s->alf[s->kopt-1][s->kopt]) break;
  s->kmax=s->kopt;

  /* simpr */
  s->indx=ivector(1,nv);
  s->simpr_a=matrix(1,nv,1,nv);
  s->del=vector(1,nv);
  s->ytemp=vector(1,nv);

  /* pzextr */
  s->c=vector(1,nv);

  s->Data = Data;
  s->derivs = derivs;
  s->jacobn = jacobn;

  return s;
}

void StiffFinalize( STIFF *s ) 
{
  free_vector(s->ytemp,1,s->nv);
  free_vector(s->del,1,s->nv);
  free_matrix(s->simpr_a,1,s->nv,1,s->nv);
  free_ivector(s->indx,1,s->nv);

  free_vector(s->c,1,s->nv);

  free_vector(s->yseq,1,s->nv);
  free_vector(s->ysav,1,s->nv);
  free_vector(s->yerr,1,s->nv);
  free_vector(s->x,1,KMAXX);
  free_vector(s->err,1,KMAXX);
  free_matrix(s->dfdy,1,s->nv,1,s->nv);
  free_vector(s->dfdx,1,s->nv);
  free_matrix(s->d,1,KMAXX,1,KMAXX);

  free(s);
}

void StiffStep(STIFF *s, double y[], double dydx[], double *xx, double htry, 
		double yscal[], double *hdid, double *hnext )
{
  int i,iq,k,kk,km, nv = s->nv;
  double errmax,fact,h,red,scale,work,wrkmin,xest;
  double *dfdx = s->dfdx,**dfdy = s->dfdy,*err = s->err;
  double *yerr = s->yerr,*ysav = s->ysav,*yseq = s->yseq;
  int reduct,exitflag=0;

  h=htry;
  for (i=1;i<=nv;i++) ysav[i]=y[i];
  (*(s->jacobn))(s->Data,*xx,y,s->dfdx,s->dfdy);
  if (*xx != s->xnew || h != (*hnext)) {
    s->first=1;
    s->kopt=s->kmax;
  }
  reduct=0;
  for (;;) {
    for (k=1;k<=s->kmax;k++) {
      s->xnew=(*xx)+h;
      if (s->xnew == (*xx)) assert(0); /* nrerror("step size underflow in stifbs"); */
      simpr(s,ysav,dydx,s->dfdx,s->dfdy,*xx,h,s->nseq[k],yseq);
      xest=SQR(h/s->nseq[k]);
      pzextr(s,k,xest,yseq,y,yerr);
      if (k != 1) {
        errmax=TINY;
        for (i=1;i<=nv;i++) errmax=FMAX(errmax,fabs(yerr[i]/yscal[i]));
        errmax /= s->eps;
        km=k-1;
        err[km]=pow(errmax/SAFE1,1.0/(2*km+1));
      }
      if (k != 1 && (k >= s->kopt-1 || s->first)) {
#ifdef USEHMIN
        if (errmax < 1.0 || h < s->hMin) {
	  if (h < s->hMin) h=s->hMin;
#else 
        if (errmax < 1.0) {
#endif
          exitflag=1;
          break;
        }
        if (k == s->kmax || k == s->kopt+1) {
          red=SAFE2/err[km];
          break;
        }
        else if (k == s->kopt && s->alf[s->kopt-1][s->kopt] < err[km]) {
            red=1.0/err[km];
            break;
          }
        else if (s->kopt == s->kmax && s->alf[km][s->kmax-1] < err[km]) {
            red=s->alf[km][s->kmax-1]*SAFE2/err[km];
            break;
          }
        else if (s->alf[km][s->kopt] < err[km]) {
          red=s->alf[km][s->kopt-1]/err[km];
          break;
        }
      }
    }
    if (exitflag) break;
    red=FMIN(red,REDMIN);
    red=FMAX(red,REDMAX);
    h *= red;
    reduct=1;
  }
  *xx=s->xnew;
  *hdid=h;
  s->first=0;
  wrkmin=1.0e35;
  for (kk=1;kk<=km;kk++) {
    fact=FMAX(err[kk],SCALMX);
    work=fact*s->a[kk+1];
    if (work < wrkmin) {
      scale=fact;
      wrkmin=work;
      s->kopt=kk+1;
    }
  }
  *hnext=h/scale;
  if (s->kopt >= k && s->kopt != s->kmax && !reduct) {
    fact=FMAX(scale/s->alf[s->kopt-1][s->kopt],SCALMX);
    if (s->a[s->kopt+1]*fact <= wrkmin) {
      *hnext=h/fact;
      s->kopt++;
    }
  }
}
#undef KMAXX
#undef IMAXX
#undef SAFE1
#undef SAFE2
#undef REDMAX
#undef REDMIN
#undef TINY
#undef SCALMX

/*
 * Simpr
 */

void simpr(STIFF *s, double y[], double dydx[], double dfdx[], double **dfdy, 
	   double xs, double htot, int nstep, double yout[] )
{
  int i,j,nn,*indx = s->indx;
  double d,h,x,**a = s->simpr_a,*del = s->del,*ytemp = s->ytemp;
  int n = s->nv;

  h=htot/nstep;
  for (i=1;i<=n;i++) {
    for (j=1;j<=n;j++) a[i][j] = -h*dfdy[i][j];
    ++a[i][i];
  }
  ludcmp(a,n,indx,&d);
  for (i=1;i<=n;i++)
    yout[i]=h*(dydx[i]+h*dfdx[i]);
  lubksb(a,n,indx,yout);
  for (i=1;i<=n;i++)
    ytemp[i]=y[i]+(del[i]=yout[i]);
  x=xs+h;
  (*(s->derivs))(s->Data,x,ytemp,yout);
  for (nn=2;nn<=nstep;nn++) {
    for (i=1;i<=n;i++)
      yout[i]=h*yout[i]-del[i];
    lubksb(a,n,indx,yout);
    for (i=1;i<=n;i++)
      ytemp[i] += (del[i] += 2.0*yout[i]);
    x += h;
    (*(s->derivs))(s->Data,x,ytemp,yout);
  }
  for (i=1;i<=n;i++)
    yout[i]=h*yout[i]-del[i];
  lubksb(a,n,indx,yout);
  for (i=1;i<=n;i++)
    yout[i] += ytemp[i];
}

/*
 * Ludcmp
 */ 

#define TINY 1.0e-20;

void ludcmp(double **a, int n, int *indx, double *d)
{
  int i,imax,j,k;
  double big,dum,sum,temp;
  double *vv;

  vv=vector(1,n);
  *d=1.0;
  for (i=1;i<=n;i++) {
    big=0.0;
    for (j=1;j<=n;j++)
      if ((temp=fabs(a[i][j])) > big) big=temp;
    if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
    vv[i]=1.0/big;
  }
  for (j=1;j<=n;j++) {
    for (i=1;i<j;i++) {
      sum=a[i][j];
      for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    big=0.0;
    for (i=j;i<=n;i++) {
      sum=a[i][j];
      for (k=1;k<j;k++)
        sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
      if ( (dum=vv[i]*fabs(sum)) >= big) {
        big=dum;
        imax=i;
      }
    }
    if (j != imax) {
      for (k=1;k<=n;k++) {
        dum=a[imax][k];
        a[imax][k]=a[j][k];
        a[j][k]=dum;
      }
      *d = -(*d);
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if (a[j][j] == 0.0) a[j][j]=TINY;
    if (j != n) {
      dum=1.0/(a[j][j]);
      for (i=j+1;i<=n;i++) a[i][j] *= dum;
    }
  }
  free_vector(vv,1,n);
}
#undef TINY

/*
 * Lubksb
 */

void lubksb(double **a, int n, int *indx, double b[])
{
  int i,ii=0,ip,j;
  double sum;

  for (i=1;i<=n;i++) {
    ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    if (ii)
      for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
    else if (sum) ii=i;
    b[i]=sum;
  }
  for (i=n;i>=1;i--) {
    sum=b[i];
    for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
    b[i]=sum/a[i][i];
  }
}

/* 
 * Pzextr
 */

void pzextr(STIFF *s, int iest, double xest, double yest[], double yz[], double dy[])
{
  int k1,j;
  double q,f2,f1,delta;
  double **d = s->d, *x = s->x, *c = s->c;
  int nv = s->nv;

  x[iest]=xest;
  for (j=1;j<=nv;j++) dy[j]=yz[j]=yest[j];
  if (iest == 1) {
    for (j=1;j<=nv;j++) d[j][1]=yest[j];
  } else {
    for (j=1;j<=nv;j++) c[j]=yest[j];
    for (k1=1;k1<iest;k1++) {
      delta=1.0/(x[iest-k1]-xest);
      f1=xest*delta;
      f2=x[iest-k1]*delta;
      for (j=1;j<=nv;j++) {
        q=d[j][k1];
        d[j][k1]=dy[j];
        delta=c[j]-q;
        dy[j]=f1*delta;
        c[j]=f2*delta;
        yz[j] += dy[j];
      }
    }
    for (j=1;j<=nv;j++) d[j][iest]=dy[j];
  }
}

/*
 * Jacobn test code
 */ 

void jacobn(double x, double y[], double dfdx[], double **dfdy, int n)
{
  int i;

  for (i=1;i<=n;i++) dfdx[i]=0.0;
  dfdy[1][1] = -0.013-1000.0*y[3];
  dfdy[1][2]=0.0;
  dfdy[1][3] = -1000.0*y[1];
  dfdy[2][1]=0.0;
  dfdy[2][2] = -2500.0*y[3];
  dfdy[2][3] = -2500.0*y[2];
  dfdy[3][1] = -0.013-1000.0*y[3];
  dfdy[3][2] = -2500.0*y[3];
  dfdy[3][3] = -1000.0*y[1]-2500.0*y[2];
}

void derivs(double x, double y[], double dydx[])
{
  dydx[1] = -0.013*y[1]-1000.0*y[1]*y[3];
  dydx[2] = -2500.0*y[2]*y[3];
  dydx[3] = -0.013*y[1]-1000.0*y[1]*y[3]-2500.0*y[2]*y[3];
}

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
  fprintf(stderr,"Numerical Recipes run-time error...\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  assert(0);
}

double *vector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
  double *v;

  v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
  if (!v) nrerror("allocation failure in vector()");
  return v-nl+NR_END;
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
  int *v;

  v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
  if (!v) nrerror("allocation failure in ivector()");
  return v-nl+NR_END;
}

unsigned char *cvector(long nl, long nh)
/* allocate an unsigned char vector with subscript range v[nl..nh] */
{
  unsigned char *v;

  v=(unsigned char *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(unsigned char)));
  if (!v) nrerror("allocation failure in cvector()");
  return v-nl+NR_END;
}

unsigned long *lvector(long nl, long nh)
/* allocate an unsigned long vector with subscript range v[nl..nh] */
{
  unsigned long *v;

  v=(unsigned long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long)));
  if (!v) nrerror("allocation failure in lvector()");
  return v-nl+NR_END;
}

double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
  double *v;

  v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
  if (!v) nrerror("allocation failure in dvector()");
  return v-nl+NR_END;
}

double **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  double **m;

  /* allocate pointers to rows */
  m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
  if (!m) nrerror("allocation failure 1 in matrix()");
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
  if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  double **m;

  /* allocate pointers to rows */
  m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
  if (!m) nrerror("allocation failure 1 in matrix()");
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
  if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  int **m;

  /* allocate pointers to rows */
  m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
  if (!m) nrerror("allocation failure 1 in matrix()");
  m += NR_END;
  m -= nrl;


  /* allocate rows and set pointers to them */
  m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
  if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

double **submatrix(double **a, long oldrl, long oldrh, long oldcl, long oldch,
  long newrl, long newcl)
/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
{
  long i,j,nrow=oldrh-oldrl+1,ncol=oldcl-newcl;
  double **m;

  /* allocate array of pointers to rows */
  m=(double **) malloc((size_t) ((nrow+NR_END)*sizeof(double*)));
  if (!m) nrerror("allocation failure in submatrix()");
  m += NR_END;
  m -= newrl;

  /* set pointers to rows */
  for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

double **convert_matrix(double *a, long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix m[nrl..nrh][ncl..nch] that points to the matrix
declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
and ncol=nch-ncl+1. The routine should be called with the address
&a[0][0] as the first argument. */
{
  long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1;
  double **m;

  /* allocate pointers to rows */
  m=(double **) malloc((size_t) ((nrow+NR_END)*sizeof(double*)));
  if (!m)  nrerror("allocation failure in convert_matrix()");
  m += NR_END;
  m -= nrl;

  /* set pointers to rows */
  m[nrl]=a-ncl;
  for(i=1,j=nrl+1;i<nrow;i++,j++) m[j]=m[j-1]+ncol;
  /* return pointer to array of pointers to rows */
  return m;
}

double ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
  long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
  double ***t;

  /* allocate pointers to pointers to rows */
  t=(double ***) malloc((size_t)((nrow+NR_END)*sizeof(double**)));
  if (!t) nrerror("allocation failure 1 in f3tensor()");
  t += NR_END;
  t -= nrl;

  /* allocate pointers to rows and set pointers to them */
  t[nrl]=(double **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double*)));
  if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
  t[nrl] += NR_END;
  t[nrl] -= ncl;

  /* allocate rows and set pointers to them */
  t[nrl][ncl]=(double *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(double)));
  if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
  t[nrl][ncl] += NR_END;
  t[nrl][ncl] -= ndl;

  for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
  for(i=nrl+1;i<=nrh;i++) {
    t[i]=t[i-1]+ncol;
    t[i][ncl]=t[i-1][ncl]+ncol*ndep;
    for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
  }

  /* return pointer to array of pointers to rows */
  return t;
}

void free_vector(double *v, long nl, long nh)
/* free a double vector allocated with vector() */
{
  free((FREE_ARG) (v+nl-NR_END));
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
  free((FREE_ARG) (v+nl-NR_END));
}

void free_cvector(unsigned char *v, long nl, long nh)
/* free an unsigned char vector allocated with cvector() */
{
  free((FREE_ARG) (v+nl-NR_END));
}

void free_lvector(unsigned long *v, long nl, long nh)
/* free an unsigned long vector allocated with lvector() */
{
  free((FREE_ARG) (v+nl-NR_END));
}

void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
  free((FREE_ARG) (v+nl-NR_END));
}

void free_matrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by matrix() */
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by imatrix() */
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

void free_submatrix(double **b, long nrl, long nrh, long ncl, long nch)
/* free a submatrix allocated by submatrix() */
{
  free((FREE_ARG) (b+nrl-NR_END));
}

void free_convert_matrix(double **b, long nrl, long nrh, long ncl, long nch)
/* free a matrix allocated by convert_matrix() */
{
  free((FREE_ARG) (b+nrl-NR_END));
}

void free_f3tensor(double ***t, long nrl, long nrh, long ncl, long nch,
  long ndl, long ndh)
/* free a double f3tensor allocated by f3tensor() */
{
  free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
  free((FREE_ARG) (t[nrl]+ncl-NR_END));
  free((FREE_ARG) (t+nrl-NR_END));
}

#define ITMAX 100
#define EPS 3.0e-8

double RootFind(double (*func)(void *Data, double), void *Data, double x1, double x2, double tol)
{
  void nrerror(char error_text[]);
  int iter;
  double a=x1,b=x2,c=x2,d,e,min1,min2;
  double fa=(*func)(Data, a),fb=(*func)(Data, b),fc,p,q,r,s,tol1,xm;

  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
    nrerror("Root must be bracketed in zbrent");
  fc=fb;
  for (iter=1;iter<=ITMAX;iter++) {
    if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
      c=a;
      fc=fa;
      e=d=b-a;
    }
    if (fabs(fc) < fabs(fb)) {
      a=b;
      b=c;
      c=a;
      fa=fb;
      fb=fc;
      fc=fa;
    }
    tol1=2.0*EPS*fabs(b)+0.5*tol;
    xm=0.5*(c-b);
    if (fabs(xm) <= tol1 || fb == 0.0) return b;
    if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
      s=fb/fa;
      if (a == c) {
        p=2.0*xm*s;
        q=1.0-s;
      } else {
        q=fa/fc;
        r=fb/fc;
        p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
        q=(q-1.0)*(r-1.0)*(s-1.0);
      }
      if (p > 0.0) q = -q;
      p=fabs(p);
      min1=3.0*xm*q-fabs(tol1*q);
      min2=fabs(e*q);
      if (2.0*p < (min1 < min2 ? min1 : min2)) {
        e=d;
        d=p/q;
      } else {
        d=xm;
        e=d;
      }
    } else {
      d=xm;
      e=d;
    }
    a=b;
    fa=fb;
    if (fabs(d) > tol1)
      b += d;
    else
      b += (xm > 0.0 ? fabs(tol1) : -fabs(tol1));
    fb=(*func)(Data, b);
  }
  nrerror("Maximum number of iterations exceeded in zbrent");
  return 0.0;
}
#undef ITMAX
#undef EPS
#endif
#endif

