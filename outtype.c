#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "pkd.h"
#include "outtype.h"
#include "floattype.h"

#ifdef COLLISIONS
#include "collision.h"
#endif /* COLLISIONS */

FLOAT ArrType(PARTICLE *p,int iType)
{
	switch (iType) {
	case OUT_DENSITY_ARRAY:
		return(p->fDensity);
	case OUT_COLOR_ARRAY:
#ifdef COLORCODE
		return(p->fColor);
#endif
	case OUT_POT_ARRAY:
		return(p->fPot);
	case OUT_AMAG_ARRAY:
		return(sqrt(p->a[0]*p->a[0] + p->a[1]*p->a[1] + p->a[2]*p->a[2]));
	case OUT_RUNG_ARRAY:
		return(p->iRung);
	case OUT_DT_ARRAY:
		return(p->dt);
#ifdef GASOLINE
	case OUT_U_ARRAY:
		return(p->u);
#ifndef NOCOOLING
	case OUT_UDOT_ARRAY:
		return(p->uDot);
	case OUT_HI_ARRAY:
		return(p->Y_HI);
	case OUT_HeI_ARRAY:
		return(p->Y_HeI);
	case OUT_HeII_ARRAY:
		return(p->Y_HeII);
#endif
#endif
	case OUT_H_ARRAY:
		return(sqrt(p->fBall2*0.25));
#ifdef COLLISIONS
	case OUT_REJECTS_ARRAY:
		/* Rejected particles indicated by their iOrder, otherwise -1 */
		return(REJECT(p) ? p->iOrder : -1);
#endif /* COLLISIONS */
	default:
		return(0.0);
		}
	}


FLOAT VecType(PARTICLE *p,int iDim,int iType)
{
	switch (iType) {
	case OUT_POS_VECTOR:
		return(p->r[iDim]);
	case OUT_VEL_VECTOR:
		return(p->v[iDim]);
	case OUT_ACCEL_VECTOR:
		return(p->a[iDim]);
#if defined(GASOLINE) || defined(SLIDING_PATCH)
	case OUT_VPRED_VECTOR:
		return(p->vPred[iDim]);
#endif
	default:
		return(0.0);
		}
	}


void pkdOutArray(PKD pkd,char *pszFileName,int iArrType)
{
	FILE *fp;
	FLOAT fOut;
	int i;

	fp = fopen(pszFileName,"a");
	assert(fp != NULL);
	/*
	 ** Write Array Elements!
	 */
	for (i=0;i<pkd->nLocal;++i) {
		fOut = ArrType(&pkd->pStore[i],iArrType);
		fprintf(fp,"%.14g\n",fOut);
		}
	i = fclose(fp);
	if (i != 0) {
		perror("pkdOutArray: could not close file");
		exit(1);
		}
	}


void pkdOutVector(PKD pkd,char *pszFileName,int iDim,int iVecType)
{
	FILE *fp;
	FLOAT fOut;
	int i;

	fp = fopen(pszFileName,"a");
	assert(fp != NULL);
	/*
	 ** Write Vector Elements!
	 */
	for (i=0;i<pkd->nLocal;++i) {
		fOut = VecType(&pkd->pStore[i],iDim,iVecType);
		fprintf(fp,"%.14g\n",fOut);
		}
	i = fclose(fp);
	if (i != 0) {
		perror("pkdOutVector: could not close file");
		exit(1);
		}
	}
