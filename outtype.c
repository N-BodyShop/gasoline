#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "pkd.h"
#include "outtype.h"

float ArrType(PARTICLE *p,int iType)
{
	switch (iType) {
#ifdef SMOOTH_CODE
	case OUT_DENSITY_ARRAY:
		return(p->fDensity);
#endif
	case OUT_POT_ARRAY:
		return(p->fPot);
	case OUT_AMAG_ARRAY:
		return(sqrt(p->a[0]*p->a[0] + p->a[1]*p->a[1] + p->a[2]*p->a[2]));
	default:
		return(0.0);
		}
	}


float VecType(PARTICLE *p,int iDim,int iType)
{
	switch (iType) {
	case OUT_POS_VECTOR:
		return(p->r[iDim]);
	case OUT_VEL_VECTOR:
		return(p->v[iDim]);
	case OUT_ACCEL_VECTOR:
		return(p->a[iDim]);
	default:
		return(0.0);
		}
	}


void pkdOutArray(PKD pkd,char *pszFileName,int iArrType)
{
	FILE *fp;
	float fOut;
	int i;

	fp = fopen(pszFileName,"a");
	assert(fp != NULL);
	/*
	 ** Write Array Elements!
	 */
	for (i=0;i<pkd->nLocal;++i) {
		fOut = ArrType(&pkd->pStore[i],iArrType);
		fprintf(fp,"%.8g\n",fOut);
		}
	fclose(fp);
	}


void pkdOutVector(PKD pkd,char *pszFileName,int iDim,int iVecType)
{
	FILE *fp;
	float fOut;
	int i;

	fp = fopen(pszFileName,"a");
	assert(fp != NULL);
	/*
	 ** Write Vector Elements!
	 */
	for (i=0;i<pkd->nLocal;++i) {
		fOut = VecType(&pkd->pStore[i],iDim,iVecType);
		fprintf(fp,"%.8g\n",fOut);
		}
	fclose(fp);
	}









