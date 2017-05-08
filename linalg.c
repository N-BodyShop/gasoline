
/** linalg.c
 * 
 *  PKDGRAV Source Code
 * 
 *  Author: Kenneth W. Flynn
 *          flynnk@astro.umd.edu
 *  Mods:   Derek C. Richardson
 *          dcr@astro.umd.edu
 *          [3-D hardwired for speed]
 */


#include <math.h>
#include <assert.h>
#include "linalg.h"

#ifdef AGGS /*DEBUG for now*/

void vectorCopy(const Vector u,Vector v)
{
	v[0] = u[0];
	v[1] = u[1];
	v[2] = u[2];
	}

void vectorScale(const Vector u,FLOAT scalar,Vector v)
{
	v[0] = u[0]*scalar;
	v[1] = u[1]*scalar;
	v[2] = u[2]*scalar;
	}

void vectorAdd(const Vector v1,const Vector v2,Vector v)
{
	v[0] = v1[0] + v2[0];
	v[1] = v1[1] + v2[1];
	v[2] = v1[2] + v2[2];
	}

void vectorSub(const Vector v1,const Vector v2,Vector v)
{
	v[0] = v1[0] - v2[0];
	v[1] = v1[1] - v2[1];
	v[2] = v1[2] - v2[2];
	}

FLOAT vectorDot(const Vector v1,const Vector v2)
{
	return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
	}

FLOAT vectorMagSq(const Vector v)
{
	return vectorDot(v,v);
	}

FLOAT vectorMag(const Vector v)
{
	return sqrt(vectorMagSq(v));
	}

void vectorNorm(Vector v)
{
	double mag = vectorMag(v);
	assert(mag > 0.0);
	vectorScale(v,1.0/mag,v);
	}

void vectorCross(const Vector v1,const Vector v2,Vector v)
{
	v[0] = v1[1]*v2[2] - v1[2]*v2[1];
	v[1] = v1[2]*v2[0] - v1[0]*v2[2];
	v[2] = v1[0]*v2[1] - v1[1]*v2[0];
	}

void vectorSet(Vector v,double x,double y,double z)
{
	v[0] = x;
	v[1] = y;
	v[2] = z;
	}

void vectorZero(Vector v)
{
	vectorSet(v,0.0,0.0,0.0);
	}

void vectorGetBasis(Vector a,Vector b,Vector c)
{
	/* Given vector "a", this routine returns orthonormal basis (a,b,c) */

	Vector v,ctmp;
	Matrix I;
	double proj;

	matrixIdentity(I); /* unit matrix */

	/* Get spanning set...first guess: choose I[1] (y-hat) and I[2] (z_hat) as 2nd & 3rd vecs */

	vectorCopy(I[1],b);
	vectorCopy(I[2],c);

	/* If "a" is actually null, set it to I[0] (x-hat) and return */

	if (a[0] == 0.0 && a[1] == 0.0 && a[2] == 0.0) {
		vectorCopy(I[0],a);
		return;
	}

	/*
	** If "a" does not have an x component, make 2nd vector I[0] (x-hat). If in
	** addition "a" does not have a y component, make 3rd vector I[1] (y-hat).
	** Now a, b, and c span 3-space.
	*/

	if (a[0] == 0.0) {
		vectorCopy(I[0],b);
		if (a[1] == 0.0)
			vectorCopy(I[1],c);
	}

	/* Construct orthonormal basis using the Gram-Schmidt orthonormalization process */

	vectorNorm(a); /* first basis vector */

	/* Construct second basis vector */

	proj = vectorDot(a,b);
	vectorScale(a,proj,v);
	vectorSub(b,v,b);
	vectorNorm(b);

	/* Construct third basis vector */

	proj = vectorDot(a,c);
	vectorScale(a,proj,v);
	vectorSub(c,v,ctmp);
	proj = vectorDot(b,c);
	vectorScale(b,proj,v);
	vectorSub(ctmp,v,c);
	vectorNorm(c);
	}

void matrixCopy(Matrix a,Matrix b)
{
	vectorCopy(a[0],b[0]);
	vectorCopy(a[1],b[1]);
	vectorCopy(a[2],b[2]);
	}

void matrixTransform(Matrix m,const Vector u,Vector v)
{
	v[0] = vectorDot(m[0],u);
	v[1] = vectorDot(m[1],u);
	v[2] = vectorDot(m[2],u);
	}

void matrixTranspose(Matrix a,Matrix b)
{
	b[0][0] = a[0][0];
	b[0][1] = a[1][0];
	b[0][2] = a[2][0];
	b[1][0] = a[0][1];
	b[1][1] = a[1][1];
	b[1][2] = a[2][1];
	b[2][0] = a[0][2];
	b[2][1] = a[1][2];
	b[2][2] = a[2][2];
	}

void matrixSwapRows(Matrix m,int row1,int row2)
{
	Vector tmp;

	vectorCopy(m[row2],tmp);
	vectorCopy(m[row1],m[row2]);
	vectorCopy(tmp,m[row1]);
	}

void matrixIdentity(Matrix m)
{
	m[0][0] = 1.0;
	m[0][1] = 0.0;
	m[0][2] = 0.0;
	m[1][0] = 0.0;
	m[1][1] = 1.0;
	m[1][2] = 0.0;
	m[2][0] = 0.0;
	m[2][1] = 0.0;
	m[2][2] = 1.0;
	}

void matrixDiagonal(const Vector v,Matrix m)
{
	m[0][0] = v[0];
	m[0][1] = 0.0;
	m[0][2] = 0.0;
	m[1][0] = 0.0;
	m[1][1] = v[1];
	m[1][2] = 0.0;
	m[2][0] = 0.0;
	m[2][1] = 0.0;
	m[2][2] = v[2];
	}

FLOAT matrixSumOffDiagElem(Matrix m)
{
	return m[0][1] + m[0][2] + m[1][0] +
		m[1][2] + m[2][0] + m[2][1];
	}

FLOAT matrixSumAbsOffDiagElem(Matrix m)
{
	return fabs(m[0][1]) + fabs(m[0][2]) + fabs(m[1][0]) +
		fabs(m[1][2]) + fabs(m[2][0]) + fabs(m[2][1]);
	}

void matrixScale(Matrix a,Scalar s,Matrix b)
{
	vectorScale(a[0],s,b[0]);
	vectorScale(a[1],s,b[1]);
	vectorScale(a[2],s,b[2]);
	}

void matrixMultiply(Matrix a,Matrix b,Matrix c)
{
	c[0][0] = a[0][0]*b[0][0] + a[0][1]*b[1][0] + a[0][2]*b[2][0];
	c[0][1] = a[0][0]*b[0][1] + a[0][1]*b[1][1] + a[0][2]*b[2][1];
	c[0][2] = a[0][0]*b[0][2] + a[0][1]*b[1][2] + a[0][2]*b[2][2];
	c[1][0] = a[1][0]*b[0][0] + a[1][1]*b[1][0] + a[1][2]*b[2][0];
	c[1][1] = a[1][0]*b[0][1] + a[1][1]*b[1][1] + a[1][2]*b[2][1];
	c[1][2] = a[1][0]*b[0][2] + a[1][1]*b[1][2] + a[1][2]*b[2][2];
	c[2][0] = a[2][0]*b[0][0] + a[2][1]*b[1][0] + a[2][2]*b[2][0];
	c[2][1] = a[2][0]*b[0][1] + a[2][1]*b[1][1] + a[2][2]*b[2][1];
	c[2][2] = a[2][0]*b[0][2] + a[2][1]*b[1][2] + a[2][2]*b[2][2];
	}

void matrixInverse(Matrix mat_in,Matrix mat_out)
{
	int row_to_pivot; /* Actual row that will be pivoted upon. */
	int pivot_row;    /* Current row of the matrix.  We will make row_to_pivot
						 equal to this by moving that row to this position. */
	FLOAT max;        /* The maximum value in the correct place.  The row we
						 pivot is based on the row with the maximum value. */
	int x,y;
	FLOAT scale,custom_scale;
	Matrix m;

	matrixCopy(mat_in,m);
	matrixIdentity(mat_out);

	/* Which row of the matrix are we on? */
	for (pivot_row = 0; pivot_row < 3; ++pivot_row) {
		/* First, we identify the largest element in the column in question and
		 *  then move it to be in the current position. */

		/* Start looking at the current row. */
		row_to_pivot = pivot_row;
		max = fabs(m[pivot_row][pivot_row]);
		for (y = pivot_row + 1; y < 3; ++y) {
			/* If element is large, mark it as maximum. */
			if (fabs(m[y][pivot_row]) > max) {
				row_to_pivot = y;
				max = fabs(m[y][pivot_row]);
				}
			}

		/* Okay, we now know what row to pivot.  Move it to the right place. */
		if (row_to_pivot != pivot_row) {
			matrixSwapRows(m,row_to_pivot,pivot_row);
			matrixSwapRows(mat_out,row_to_pivot,pivot_row);
			}

		/* Row is in place, now we move on to the actual pivot. */

		/* First we compute how much we need to scale the pivot row by, which is
		 *  1 / pivot_element. */
		scale = 1.0 / m[pivot_row][pivot_row];

		/* Next, we compute the base pivot_row.  We will just leave this in the
		 *  matrix at the appropriate spot. */
		for (x = 0; x < 3; ++x) {
			m[pivot_row][x] *= scale;
			mat_out[pivot_row][x] *= scale;
			}

		/* Finally, for all rows except the pivot row, we subtract the pivot row 
		 *  of the matrix.  This is the actual pivot.  We will need a custom scale
		 *  factor for each row to eliminate the element. */
		for (y = 0; y < 3; ++y) {
			if (y != pivot_row) {
				/* Get the custom scale for this row. */
				custom_scale = m[y][pivot_row];

				/* Now pivot. */
				for (x = 0; x < 3; ++x) {
					m[y][x] -= custom_scale * m[pivot_row][x];
					mat_out[y][x] -= custom_scale * mat_out[pivot_row][x];
					}
				}
			}
		}
	}

#define JACOBI_N 3 /*DEBUG hard-wired for 3D*/

void jacobi(Matrix m,Vector eig_vals,Matrix eig_vecs)
{
	int y,x;                  /* Row, column of matrix for element to be 
							   *  eliminated. */
	int j;                    /* Column of matrix. */
	int sweep_count;          /* Current number of sweeps made */
	FLOAT off_diag_sum;       /* Sum of off diagonal elements. */
	FLOAT threshold;          /* Threshold for performing the rotation.  For the
							   *  first three sweeps, this is set equal to 1/5
							   *  the off diagonal sum divided by n^2.  After
							   *  four sweeps, this is 0. */
	Matrix a;                 /* Copy of input matrix. */
	FLOAT s;                  /* Rotation angle sine */
	FLOAT c;                  /* Rotation angle cosine */
	FLOAT t;
	FLOAT tau;
	FLOAT theta; 
	FLOAT temp;

	/* Make copy of input matrix. */
	matrixCopy(m,a);

	/* Initialize result. */
	matrixIdentity(eig_vecs);
 
	/* Do Jacobi rotations */
	for (sweep_count = 0; sweep_count < MAX_JACOBI_SWEEPS; ++sweep_count) {
		/* First, check to see if we are done. */
		off_diag_sum = matrixSumAbsOffDiagElem(a);
		if (off_diag_sum == 0.0) { /* Make that underflow... */
			/* Store results */
			for (x = 0; x < JACOBI_N; ++x)
				eig_vals[x] = a[x][x];
			/* Eig vectors is ready */
			return;
			}

		/* Not done yet, so determine the threshold value for performing the
		 *  rotation. */
		if (sweep_count < 3)
			threshold = 0.2 * off_diag_sum / (JACOBI_N * JACOBI_N);
		else
			threshold = 0.0;
  
		/* Do a sweep; this means do a rotation for each off diag element in the 
		 *  matrix.  Here, (x,y) are the coordinates for the element being 
		 *  "eliminated."  We only need eliminate the top triangle of the symmetric
		 *  matrix. */
		for (y = 0; y < JACOBI_N - 1; ++y)
			for (x = y + 1; x < JACOBI_N; ++x) {
    
				/* Now, we check to see if we should bother with the rotation.  If the
				 *  element being "eliminated" is too small, then we don't bother.  Too
				 *  small in this case is when the element * 100 is small enough that
				 *  adding it to either the diagonal element in its row or the diagonal
				 *  element in its column is not enough to change the diagonal element
				 *  within machine precision.  We only do this check if we are on the
				 *  fourth or higher sweep.  If we "don't bother" we set the element to
				 *  0 (eliminating it) and then move on. */
				if (sweep_count > 3)
					if (((100.0 * a[y][x] + fabs(a[x][x])) == fabs(a[x][x])) &&
						((100.0 * a[y][x] + fabs(a[y][y])) == fabs(a[y][y])))
						{
						a[y][x] = a[x][y] = 0.0;
						continue;
						}

				/* Now we check the threshold value. */
				if (fabs(a[y][x]) <= threshold)
					continue;
     
				/* Okay, we must perform the rotation.
				 * theta -> t -> c, s, tau -> rotation */

				/* theta, t */

				/* Start theta calculation.  Check size; if too big, get t without
				 *  the need to square theta, otherwise, get t the normal way. */
				theta = a[x][x] - a[y][y];
				/* If theta too big... */
				if ((fabs(theta) + 100.0 * a[y][x]) == fabs(theta)) {
					/* Do t = 1/2/theta (note that two formulas are combined here) */
					t = a[y][x] / theta;
					}
				else { /* theta just right... */
					theta = 0.5 * theta / a[y][x];
					/* t = sgn(theta)/( |theta| + sqrt(theta^2 + 1) ). */
					t = 1.0/(fabs(theta) + sqrt(1.0+theta*theta));
					if (theta < 0.0)
						t = -t;
					}
				/* c, s, tau */
				c = 1.0 / sqrt(1 + t*t);
				s = t * c;
				tau = s / (1.0 + c);

				/* Really do rotation, finally! */
				a[y][y] -= t * a[y][x];
				a[x][x] += t * a[y][x];
				a[y][x] = a[x][y] = 0.0;

				for (j = 0; j < JACOBI_N; ++j)
					if (j != x && j != y) {
						temp = a[j][x];

						/* Do column w/ index x */
						a[j][x] += s*(a[j][y] - tau * a[j][x]);
						a[x][j] = a[j][x];

						/* Do column w/ index y */
						a[j][y] -= s*(temp + tau * a[j][y]);
						a[y][j] = a[j][y];
						}

				/* Update result. */
				for (j = 0; j < JACOBI_N; ++j) {
					temp = eig_vecs[j][x];

					/* Do column w/ index x */
					eig_vecs[j][x] += s*(eig_vecs[j][y] - tau * eig_vecs[j][x]);

					/* Do column w/ index y */
					eig_vecs[j][y] -= s*(temp + tau * eig_vecs[j][y]);  
					}
				}
		}

	/* Should never get here.  Means we had too many sweeps. */
	assert(0);
	}

#endif /* AGGS */
