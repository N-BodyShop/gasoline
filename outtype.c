#include <stdio.h>
#include <stdlib.h>
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
	case OUT_IORDER_ARRAY:
	    return((FLOAT) p->iOrder);
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
	case OUT_SOFT_ARRAY:
	        return(p->fSoft);
#ifdef GASOLINE
	case OUT_U_ARRAY:
		return(p->u);
#ifndef NOCOOLING
	case OUT_UDOT_ARRAY:
		return(p->uDot);
	case OUT_COOL_ARRAY0:
		return(COOL_ARRAY0(&p->CoolParticle));
	case OUT_COOL_ARRAY1:
		return(COOL_ARRAY1(&p->CoolParticle));
	case OUT_COOL_ARRAY2:
		return(COOL_ARRAY2(&p->CoolParticle));
#endif
	case OUT_BALSARASWITCH_ARRAY:
		return(p->BalsaraSwitch);
	case OUT_DIVV_ARRAY:
		return(p->divv);
	case OUT_MUMAX_ARRAY:
		return(p->mumax);
	case OUT_DIVONCONH_ARRAY:
	        return(p->divv/(p->c/sqrt(p->fBall2*0.25)));
	case OUT_DIVONCONX_ARRAY:
	        return(p->divv/(p->c/pow(p->fMass/p->fDensity,1./3.)));
        case OUT_PDV_ARRAY:
  	        return(p->PdV);
#ifdef PDVDEBUG
	case OUT_PDVPRES_ARRAY:
	        return(p->PdVpres);
	case OUT_PDVVISC_ARRAY:
	        return(p->PdVvisc);
#endif
#ifdef SHOCKTRACK
	case OUT_SHOCKTRACKER_ARRAY:
	        return(p->ShockTracker);
	case OUT_DIVRHOV_ARRAY:
	        return(p->divrhov);
#endif
#ifdef STARFORM
	case OUT_IGASORDER_ARRAY:
	    return((FLOAT) p->iGasOrder);
	case OUT_MASSFORM_ARRAY:
	    return((FLOAT) p->fMassForm);
#endif
#ifdef SIMPLESF
	case OUT_TCOOLAGAIN_ARRAY:
		return(p->fTimeForm);
	case OUT_MSTAR_ARRAY:
		return(p->fMassStar);
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
#ifdef NEED_VPRED
	case OUT_VPRED_VECTOR:
		return(p->vPred[iDim]);
#endif
#if defined(GASOLINE)
#if defined(SHOCKTRACK)
	case OUT_GRADRHO_VECTOR:
		return(p->gradrho[iDim]);
	case OUT_ACCELPRES_VECTOR:
		return(p->aPres[iDim]);
#endif
#if defined(STARFORM)
	case OUT_RFORM_VECTOR:
	    return(p->rForm[iDim]);
	case OUT_VFORM_VECTOR:
	    return(p->vForm[iDim]);
#endif
#endif
	default:
		return(0.0);
		}
	}


void pkdOutArray(PKD pkd,char *pszFileName,int iArrType, int iBinaryOutput)
{
	FILE *fp;
	FLOAT fOut;
	float FloatOut;
	double DoubleOut;
	int IntOut;
	long long LongOut;
	int i;

	fp = fopen(pszFileName,"a");
	assert(fp != NULL);
	/*
	 ** Write Array Elements!
	 */
	switch (iBinaryOutput) {
	case 0:
		for (i=0;i<pkd->nLocal;++i) {
			fOut = ArrType(&pkd->pStore[i],iArrType);
			fprintf(fp,"%.14g\n",fOut);
			}
		break;
	case 1:
		for (i=0;i<pkd->nLocal;++i) {
			switch (iArrType) {
			case OUT_IGASORDER_ARRAY:
			case OUT_IORDER_ARRAY:
				IntOut = pkd->pStore[i].iOrder;
				fwrite(&IntOut, sizeof(IntOut), 1, fp );
				break;
			default:
				FloatOut = ArrType(&pkd->pStore[i],iArrType);
				fwrite(&FloatOut, sizeof(float), 1, fp );
				}
			}
		break;
	case 2:
		for (i=0;i<pkd->nLocal;++i) {
			switch (iArrType) {
			case OUT_IGASORDER_ARRAY:
			case OUT_IORDER_ARRAY:
				LongOut = pkd->pStore[i].iOrder;
				fwrite(&LongOut, sizeof(LongOut), 1, fp );
				break;
			default:
				DoubleOut = ArrType(&pkd->pStore[i],iArrType);
				fwrite(&DoubleOut, sizeof(double), 1, fp );
				}
			}

		break;
	case 3:
		for (i=0;i<pkd->nLocal;++i) {
			switch (iArrType) {
			case OUT_IGASORDER_ARRAY:
			case OUT_IORDER_ARRAY:
				fwrite(&(pkd->pStore[i].iOrder), sizeof(pkd->pStore[i].iOrder), 1, fp );
				break;
			default:
				fOut = ArrType(&pkd->pStore[i],iArrType);
				fwrite(&fOut, sizeof(FLOAT), 1, fp );
				}
			}
		break;
		}

	i = fclose(fp);
	if (i != 0) {
		perror("pkdOutArray: could not close file");
		exit(1);
		}
	}


void pkdOutVector(PKD pkd,char *pszFileName,int iDim,int iVecType,int iBinaryOutput)
{
	FILE *fp;
	FLOAT fOut;
	float FloatOut;
	double DoubleOut;
	int i;

	fp = fopen(pszFileName,"a");
	assert(fp != NULL);
	/*
	 ** Write Vector Elements!
	 */
	switch (iBinaryOutput) {
	case 0:
		if (iDim < 0) {
			for (i=0;i<pkd->nLocal;++i) {
				fOut = VecType(&pkd->pStore[i],0,iVecType);
				fprintf(fp,"%.14g\n",fOut);
				fOut = VecType(&pkd->pStore[i],1,iVecType);
				fprintf(fp,"%.14g\n",fOut);
				fOut = VecType(&pkd->pStore[i],2,iVecType);
				fprintf(fp,"%.14g\n",fOut);
				}
			
			}
		else {
			for (i=0;i<pkd->nLocal;++i) {
				fOut = VecType(&pkd->pStore[i],iDim,iVecType);
				fprintf(fp,"%.14g\n",fOut);
				}
			}
		break;
	case 1:
		if (iDim < 0) {
			for (i=0;i<pkd->nLocal;++i) {
				FloatOut = VecType(&pkd->pStore[i],0,iVecType);
				fwrite(&FloatOut, sizeof(float), 1, fp );
				FloatOut = VecType(&pkd->pStore[i],1,iVecType);
				fwrite(&FloatOut, sizeof(float), 1, fp );
				FloatOut = VecType(&pkd->pStore[i],2,iVecType);
				fwrite(&FloatOut, sizeof(float), 1, fp );
				}
			
			}
		else {
			for (i=0;i<pkd->nLocal;++i) {
				FloatOut = VecType(&pkd->pStore[i],iDim,iVecType);
				fwrite(&FloatOut, sizeof(float), 1, fp );
				}
			}
		break;
	case 2:
		if (iDim < 0) {
			for (i=0;i<pkd->nLocal;++i) {
				DoubleOut = VecType(&pkd->pStore[i],0,iVecType);
				fwrite(&DoubleOut, sizeof(double), 1, fp );
				DoubleOut = VecType(&pkd->pStore[i],1,iVecType);
				fwrite(&DoubleOut, sizeof(double), 1, fp );
				DoubleOut = VecType(&pkd->pStore[i],2,iVecType);
				fwrite(&DoubleOut, sizeof(double), 1, fp );
				}
			
			}
		else {
			for (i=0;i<pkd->nLocal;++i) {
				DoubleOut = VecType(&pkd->pStore[i],iDim,iVecType);
				fwrite(&DoubleOut, sizeof(double), 1, fp );
				}
			}
		break;
	case 3:
		if (iDim < 0) {
			for (i=0;i<pkd->nLocal;++i) {
				fOut = VecType(&pkd->pStore[i],0,iVecType);
				fwrite(&fOut, sizeof(FLOAT), 1, fp );
				fOut = VecType(&pkd->pStore[i],1,iVecType);
				fwrite(&fOut, sizeof(FLOAT), 1, fp );
				fOut = VecType(&pkd->pStore[i],2,iVecType);
				fwrite(&fOut, sizeof(FLOAT), 1, fp );
				}
			
			}
		else {
			for (i=0;i<pkd->nLocal;++i) {
				fOut = VecType(&pkd->pStore[i],iDim,iVecType);
				fwrite(&fOut, sizeof(FLOAT), 1, fp );
				}
			}
		break;
		}

	i = fclose(fp);
	if (i != 0) {
		perror("pkdOutVector: could not close file");
		exit(1);
		}
	}
