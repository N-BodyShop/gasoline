#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "pkd.h"
#include "outtype.h"
#include "floattype.h"
#ifdef CRAY_XT3
#include "../xdr/types.h"
#include "../xdr/xdr.h"
#else
#include <rpc/types.h>
#include <rpc/xdr.h>
#endif

#ifdef COLLISIONS
#include "collision.h"
#endif /* COLLISIONS */

FLOAT VecType(PKD pkd, PARTICLE *p,int iDim,int iType)
{
#ifdef GASOLINE
  FLOAT vTemp;
#endif

	switch (iType) {
	case OUT_IORDER_ARRAY:
	    return((FLOAT) p->iOrder);
#ifdef STARFORM
	case OUT_DENSITYFORM_ARRAY:
#endif
#ifdef GASOLINE
        case OUT_GASDENSITY_ARRAY:
#endif
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
	case OUT_SPHDT_ARRAY:
	case OUT_DT_ARRAY:
		return(p->dt);
	case OUT_SOFT_ARRAY:
	        return(p->fSoft);
#ifdef GASOLINE
        case OUT_TEMPFORM_ARRAY:
	case OUT_U_ARRAY:
		return(p->u);
	case OUT_TEMP_ARRAY:
#ifndef NOCOOLING
                vTemp = CoolCodeEnergyToTemperature( pkd->Cool, &p->CoolParticle, p->u );
#else
                vTemp = pkd->duTFac*p->u;
#endif
		return(vTemp);
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
	case OUT_METALS_ARRAY:
  	        return(p->fMetals);
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
	case OUT_TIMEFORM_ARRAY:
	    return((FLOAT) p->fTimeForm);
	case OUT_MASSFORM_ARRAY:
	    return((FLOAT) p->fMassForm);
	case OUT_COOLTURNONTIME_ARRAY:
	    return((FLOAT) p->fTimeCoolIsOffUntil);
	case OUT_OXYGENMASSFRAC_ARRAY:
	    return((FLOAT) p->fMFracOxygen);
	case OUT_IRONMASSFRAC_ARRAY:
	    return((FLOAT) p->fMFracIron);
#endif
#ifdef SIMPLESF
	case OUT_TCOOLAGAIN_ARRAY:
		return(p->fTimeForm);
	case OUT_MSTAR_ARRAY:
		return(p->fMassStar);
#endif
	case OUT_SPHH_ARRAY:
#endif
	case OUT_H_ARRAY:
		return(sqrt(p->fBall2*0.25));
#ifdef COLLISIONS
	case OUT_REJECTS_ARRAY:
		/* Rejected particles indicated by their iOrder, otherwise -1 */
		return(REJECT(p) ? p->iOrder : -1);
#endif /* COLLISIONS */
	case OUT_POS_VECTOR:
		return(p->r[iDim]);
	case OUT_VEL_VECTOR:
		return(p->v[iDim]);
	case OUT_MASS_ARRAY:
		return(p->fMass);
	case OUT_ACCELG_VECTOR:
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

void VecFilename(char *achFile, int iType)
{
	switch (iType) {
	case OUT_IORDER_ARRAY:
	    strncat(achFile,"iord",256);
            break;
        case OUT_DENSITY_ARRAY:
		strncat(achFile,"den",256);
            break;
	case OUT_COLOR_ARRAY:
#ifdef COLORCODE
		return(p->fColor);
#endif
	case OUT_POT_ARRAY:
            strncat(achFile,"pot",256);
            break;
	case OUT_AMAG_ARRAY:
            strncat(achFile,"amag",256);
            break;
	case OUT_RUNG_ARRAY:
            strncat(achFile,"rung",256);
            break;

	case OUT_MASS_ARRAY:
            strncat(achFile,"mass",256);
            break;
	case OUT_DT_ARRAY:
		strncat(achFile,"dt",256);
            break;
	case OUT_SPHDT_ARRAY:
		strncat(achFile,"SPHdt",256);
            break;
	case OUT_SOFT_ARRAY:
	        strncat(achFile,"soft",256);
            break;
#ifdef GASOLINE
	case OUT_TEMPFORM_ARRAY:
		strncat(achFile,"Tform",256);
            break;
	case OUT_TEMP_ARRAY:
		strncat(achFile,"temperature",256);
            break;
        case OUT_GASDENSITY_ARRAY:
		strncat(achFile,"GasDensity",256);
            break;
#ifndef NOCOOLING
/*	case OUT_UDOT_ARRAY:
		return(p->uDot);
*/
	case OUT_COOL_ARRAY0:
		strncat(achFile,"HI",256);
            break;
	case OUT_COOL_ARRAY1:
		strncat(achFile,"HeI",256);
            break;
	case OUT_COOL_ARRAY2:
		strncat(achFile,"HeII",256);
            break;
#endif
	case OUT_BALSARASWITCH_ARRAY:
		strncat(achFile,"BSw",256);
            break;
	case OUT_DIVV_ARRAY:
            strncat(achFile,"divv",256);
            break;
	case OUT_MUMAX_ARRAY:
            strncat(achFile,"mumax",256);
            break;
	case OUT_DIVONCONH_ARRAY:
            strncat(achFile,"dch",256);
            break;
	case OUT_DIVONCONX_ARRAY:
            strncat(achFile,"dcx",256);
            break;
        case OUT_PDV_ARRAY:
            strncat(achFile,"PdV",256);
            break;
	case OUT_PDVPRES_ARRAY:
	        strncat(achFile,"PdVpres",256);
            break;
	case OUT_PDVVISC_ARRAY:
	        strncat(achFile,"PdVvisc",256);
            break;
	case OUT_METALS_ARRAY:
	    strncat(achFile,"Metals",256);
            break;
#ifdef SHOCKTRACK
	case OUT_SHOCKTRACKER_ARRAY:
	        strncat(achFile,"ST",256);
            break;
	case OUT_DIVRHOV_ARRAY:
	        strncat(achFile,"divrhov",256);
            break;
#endif
#ifdef STARFORM
	case OUT_IGASORDER_ARRAY:
	    strncat(achFile,"igasorder",256);
            break;
	case OUT_TIMEFORM_ARRAY:
	    strncat(achFile,"timeform",256);
            break;
	case OUT_MASSFORM_ARRAY:
	    strncat(achFile,"massform",256);
            break;
        case OUT_DENSITYFORM_ARRAY:
		strncat(achFile,"rhoform",256);
            break;
	case OUT_COOLTURNONTIME_ARRAY:
	    strncat(achFile,"coolontime",256);
            break;
	case OUT_OXYGENMASSFRAC_ARRAY:
	    strncat(achFile,"OxMassFrac",256);
            break;
	case OUT_IRONMASSFRAC_ARRAY:
	    strncat(achFile,"FeMassFrac",256);
            break;
#endif
#ifdef SIMPLESF
	case OUT_TCOOLAGAIN_ARRAY:
		strncat(achFile,"tCoolAgain",256);
            break;
	case OUT_MSTAR_ARRAY:
		strncat(achFile,"mStar",256);
            break;
#endif
	case OUT_SPHH_ARRAY:
		strncat(achFile,"SPHH",256);
            break;
#endif
	case OUT_H_ARRAY:
		strncat(achFile,"smoothlength",256);
            break;
	case OUT_POS_VECTOR:
	    strncat(achFile,"pos",256);
            break;
	case OUT_VEL_VECTOR:
	    strncat(achFile,"vel",256);
            break;
	case OUT_ACCEL_VECTOR:
	    strncat(achFile,"acc",256);
            break;
	case OUT_ACCELG_VECTOR:
	    strncat(achFile,"accg",256);
            break;
#ifdef NEED_VPRED
	case OUT_VPRED_VECTOR:
	    strncat(achFile,"vpred",256);
            break;
#endif
#if defined(GASOLINE)
#if defined(SHOCKTRACK)
	case OUT_GRADRHO_VECTOR:
		strncat(achFile,"gradrho",256);
            break;
	case OUT_ACCELPRES_VECTOR:
	    strncat(achFile,"accp",256);
            break;
#endif
#if defined(STARFORM)
	case OUT_RFORM_VECTOR:
	    strncat(achFile,"rform",256);
            break;
	case OUT_VFORM_VECTOR:
	    strncat(achFile,"vform",256);
            break;
#endif
#endif
	default:
            assert(1);
		}
	}

void pkdOutNChilada(PKD pkd,char *pszFileName,int nGasStart, int nDarkStart, int nStarStart, int iVecType, float minValue[3][3], float maxValue[3][3], double duTFac)
{
    FILE *gasFp = NULL, *darkFp = NULL, *starFp = NULL;
    float fOut, min[3][3], max[3][3];
    int i, iDim, nDim, headerlength;
    int nGas, nDark, nStar;
    char darkFileName[256],gasFileName[256],starFileName[256];
    XDR gasXdrs, darkXdrs, starXdrs;
    for(i=0;i<3;i++){
        for(iDim=0;iDim<3;iDim++){
            min[i][iDim] = 1e30;
            max[i][iDim] = -1e30;
            }
        }
    
    nGas = pkd->nGas; nDark = pkd->nDark; nStar = pkd->nStar;
    pkd->duTFac = duTFac;
    switch (iVecType){
        /* Gas only floats */
        case OUT_COOLTURNONTIME_ARRAY:
        case OUT_DIVV_ARRAY:
        case OUT_TCOOLAGAIN_ARRAY:
        case OUT_MSTAR_ARRAY:
        case OUT_COOL_ARRAY0:
        case OUT_COOL_ARRAY1:
        case OUT_COOL_ARRAY2:
        case OUT_SPHH_ARRAY:
        case OUT_TEMP_ARRAY:
        case OUT_GASDENSITY_ARRAY:
        case OUT_PDVPRES_ARRAY:
        case OUT_PDVVISC_ARRAY:
            nDark=nStar=0;
            break;
        case OUT_IGASORDER_ARRAY:
        case OUT_TIMEFORM_ARRAY:
        case OUT_MASSFORM_ARRAY:
        case OUT_RFORM_VECTOR:
        case OUT_VFORM_VECTOR:
        case OUT_DENSITYFORM_ARRAY:
        case OUT_TEMPFORM_ARRAY:
            nGas=nDark=0;
            break;
        case OUT_METALS_ARRAY:
        case OUT_OXYGENMASSFRAC_ARRAY:
        case OUT_IRONMASSFRAC_ARRAY:
            nDark=0;
            break;
        }
        
    /*
     * N-Chilada has a 28 byte header (see FieldHeader in
     * structures/tree_xdr.h): 4 byte magic + 8 byte time + 8 byte
     * numParticles + 4 byte dimension + 4 byte data type code.  The
     * header is followed by a min and max field 6 numbers for vectors
     * and 2 for scalars.
     */
    /*
     * XXX WARNING:  Seeks below are assuming sizeof(float) = 4, since
     * RFCs for XDR assume 4 byte floats.
     */
    
    if(iVecType > OUT_1D3DSPLIT) {
        nDim=3;
        headerlength = 28 + 6*4;
        }
    else {
        nDim=1;
        headerlength = 28 + 2*4;
        }
    
    if(nGas){
        sprintf(gasFileName, "%s/gas/",pszFileName);
        VecFilename(gasFileName,iVecType);
        gasFp = fopen(gasFileName,"r+");
        mdlassert(pkd->mdl,gasFp != NULL);
        xdrstdio_create(&gasXdrs,gasFp,XDR_ENCODE);
        }
    
    if(nDark){
        sprintf(darkFileName, "%s/dark/",pszFileName);
        VecFilename(darkFileName,iVecType);
        darkFp = fopen(darkFileName,"r+");
        mdlassert(pkd->mdl,darkFp != NULL);
        xdrstdio_create(&darkXdrs,darkFp,XDR_ENCODE);
        }
        
    if(nStar){
        sprintf(starFileName, "%s/star/",pszFileName);
        VecFilename(starFileName,iVecType);
        starFp = fopen(starFileName,"r+");
        mdlassert(pkd->mdl,starFp != NULL);
        xdrstdio_create(&starXdrs,starFp,XDR_ENCODE);
        }
    
    switch (iVecType) {
#ifdef STARFORM
        case OUT_IGASORDER_ARRAY:
	  if(nStar) {
	    pkdGenericSeek(pkd,starFp, nStarStart, headerlength, 4);
            for (i=0;i<pkd->nLocal;++i) {
                if (pkdIsStar(pkd,&pkd->pStore[i])) {
                    xdr_int(&starXdrs,&(pkd->pStore[i].iGasOrder));
                    min[2][0] = (pkd->pStore[i].iOrder > min[2][0]) ? min[2][0]: pkd->pStore[i].iOrder;
                    max[2][0] = (pkd->pStore[i].iOrder < max[2][0]) ? max[2][0]: pkd->pStore[i].iOrder;
                    }
                }
	    }
            break;
#endif
        case OUT_IORDER_ARRAY:
            if(nGas) pkdGenericSeek(pkd,gasFp, nGasStart,headerlength, 4);
            if(nDark) pkdGenericSeek(pkd,darkFp, nDarkStart,headerlength, 4);
            if(nStar) pkdGenericSeek(pkd,starFp, nStarStart,headerlength, 4);
            for (i=0;i<pkd->nLocal;++i) {
                if (nGas && pkdIsGas(pkd,&pkd->pStore[i])) {
                    xdr_int(&gasXdrs,&(pkd->pStore[i].iOrder));
                    min[0][0] = (pkd->pStore[i].iOrder > min[0][0]) ? min[0][0]: pkd->pStore[i].iOrder;
                    max[0][0] = (pkd->pStore[i].iOrder < max[0][0]) ? max[0][0]: pkd->pStore[i].iOrder;
                    }
                if (nDark && pkdIsDark(pkd,&pkd->pStore[i])) { 
                    xdr_int(&darkXdrs,&(pkd->pStore[i].iOrder));
                    min[1][0] = (pkd->pStore[i].iOrder > min[1][0]) ? min[1][0]: pkd->pStore[i].iOrder;
                    max[1][0] = (pkd->pStore[i].iOrder < max[1][0]) ? max[1][0]: pkd->pStore[i].iOrder;
                    }
                if (nStar && pkdIsStar(pkd,&pkd->pStore[i])) {
                    xdr_int(&starXdrs,&(pkd->pStore[i].iOrder));
                    min[2][0] = (pkd->pStore[i].iOrder > min[2][0]) ? min[2][0]: pkd->pStore[i].iOrder;
                    max[2][0] = (pkd->pStore[i].iOrder < max[2][0]) ? max[2][0]: pkd->pStore[i].iOrder;
                    }
                }
            break;
#ifdef GASOLINE				
        case OUT_TEMP_ARRAY:
        case OUT_GASDENSITY_ARRAY:
	  if(nGas) {
            pkdGenericSeek(pkd,gasFp, nGasStart,headerlength,sizeof(pkd->pStore[i].iOrder));
            for (i=0;i<pkd->nLocal;++i) {
                if (pkdIsGas(pkd,&pkd->pStore[i])) {
                    fOut = VecType(pkd, &pkd->pStore[i],0,iVecType);
                    min[0][0] = (fOut > min[0][0]) ? min[0][0]: fOut;
                    max[0][0] = (fOut < max[0][0]) ? max[0][0]: fOut;
                    xdr_float(&gasXdrs,&fOut);
                    }
                }
	    }
            break;
#endif
        default:
            if(nGas) pkdGenericSeek(pkd,gasFp, nGasStart*nDim,
                                headerlength,sizeof(fOut));
            if(nDark) pkdGenericSeek(pkd,darkFp, nDarkStart*nDim,
                                headerlength,sizeof(fOut));
            if(nStar) pkdGenericSeek(pkd,starFp, nStarStart*nDim,
                                headerlength,sizeof(fOut));
            for (i=0;i<pkd->nLocal;++i) {
                for (iDim = 0; iDim <nDim; iDim++) {
                    fOut = VecType(pkd, &pkd->pStore[i],iDim,iVecType);
                    if (nGas && pkdIsGas(pkd,&pkd->pStore[i])) {
                        xdr_float(&gasXdrs,&fOut);
                        min[0][iDim] = (fOut > min[0][iDim]) ? min[0][iDim]: fOut;
                        max[0][iDim] = (fOut < max[0][iDim]) ? max[0][iDim]: fOut;
                        }
                    if (nDark && pkdIsDark(pkd,&pkd->pStore[i])) {
                        xdr_float(&darkXdrs,&fOut);
                        min[1][iDim] = (fOut > min[1][iDim]) ? min[1][iDim]: fOut;
                        max[1][iDim] = (fOut < max[1][iDim]) ? max[1][iDim]: fOut;
                        }
                    if (nStar && pkdIsStar(pkd,&pkd->pStore[i])) {
                        xdr_float(&starXdrs,&fOut);
                        min[2][iDim] = (fOut > min[2][iDim]) ? min[2][iDim]: fOut;
                        max[2][iDim] = (fOut < max[2][iDim]) ? max[2][iDim]: fOut;
                        }
                    }
                }
        }
    if(nDark) xdr_destroy(&darkXdrs);
    if(nGas) xdr_destroy(&gasXdrs);
    if(nStar) xdr_destroy(&starXdrs);
    if(nDark) fclose(darkFp);
    if(nGas) fclose(gasFp);
    if(nStar) fclose(starFp);
    for(i=0;i<3;i++){
        for(iDim=0;iDim<nDim;iDim++){
            minValue[i][iDim] = min[i][iDim];
            maxValue[i][iDim] = max[i][iDim];
            }
        }
    }

void xdr_FLOAT(XDR *xdrs, FLOAT *fIn) 
{
#ifdef SINGLE
    xdr_float(xdrs, fIn);
#else
    xdr_double(xdrs, fIn);
#endif
    }

void pkdOutVector(PKD pkd,char *pszFileName,int nStart, int iDim,int iVecType,int iBinaryOutput, int N, int bStandard)
{
    FILE *fp;
    FLOAT fOut;
    float FloatOut;
    double DoubleOut;
    int IntOut;
    long LongOut;
    int i;
    char vecFileName[256];

    strcpy(vecFileName, pszFileName);
    VecFilename(vecFileName,iVecType);
    if(iBinaryOutput){ 
        fp = fopen(vecFileName,"r+");
        mdlassert(pkd->mdl,fp != NULL);
        } else {
        fp = fopen(vecFileName,"a");
        mdlassert(pkd->mdl,fp != NULL);
        }
    /*
     ** Write Vector Elements!
     */
    if( bStandard && iBinaryOutput ){
        XDR xdrs;
        xdrstdio_create(&xdrs,fp,XDR_ENCODE);
        switch (iBinaryOutput) {
            case 1:
                if (iDim < 0) {
                    pkdGenericSeek(pkd,fp, nStart*3,sizeof(int),sizeof(FloatOut));
                    for (i=0;i<pkd->nLocal;++i) {
                        FloatOut = VecType(pkd, &pkd->pStore[i],0,iVecType);
                        xdr_float(&xdrs,&FloatOut);
                        FloatOut = VecType(pkd, &pkd->pStore[i],1,iVecType);
                        xdr_float(&xdrs,&FloatOut);
                        FloatOut = VecType(pkd, &pkd->pStore[i],2,iVecType);
                        xdr_float(&xdrs,&FloatOut);
                        }
                    } else {
                    switch (iVecType) {
#ifdef STARFORM
                        case OUT_IGASORDER_ARRAY:
                            pkdGenericSeek(pkd,fp, nStart,sizeof(int),sizeof(IntOut));
                            for (i=0;i<pkd->nLocal;++i) {
                                xdr_int(&xdrs,&(pkd->pStore[i].iGasOrder));
                                }
                            break;
#endif
                        case OUT_IORDER_ARRAY:
                            pkdGenericSeek(pkd,fp,nStart, sizeof(int), sizeof(IntOut));
                            for (i=0;i<pkd->nLocal;++i) {
                                xdr_int(&xdrs,&(pkd->pStore[i].iOrder));
                                }
                            break;
                        default:
                            pkdGenericSeek(pkd,fp, N*iDim+nStart,sizeof(int),sizeof(FloatOut));
                            for (i=0;i<pkd->nLocal;++i) {
                                FloatOut = VecType(pkd, &pkd->pStore[i],iDim,iVecType);
                                xdr_float(&xdrs,&FloatOut);
                                }
                        }
                    }
                    break;
            case 2:
                if (iDim < 0) {
                    pkdGenericSeek(pkd,fp, nStart*3,sizeof(int),sizeof(DoubleOut));
                    for (i=0;i<pkd->nLocal;++i) {
                        DoubleOut = VecType(pkd, &pkd->pStore[i],0,iVecType);
                        xdr_double(&xdrs,&DoubleOut);
                        DoubleOut = VecType(pkd, &pkd->pStore[i],1,iVecType);
                        xdr_double(&xdrs,&DoubleOut);
                        DoubleOut = VecType(pkd, &pkd->pStore[i],2,iVecType);
                        xdr_double(&xdrs,&DoubleOut);
                        }
                        
                    } else {
                    switch (iVecType) {
#ifdef STARFORM
                        case OUT_IGASORDER_ARRAY:
                            pkdGenericSeek(pkd,fp, nStart,sizeof(int),sizeof(LongOut));
                            for (i=0;i<pkd->nLocal;++i) {
                                LongOut = pkd->pStore[i].iGasOrder;
                                xdr_long(&xdrs,&LongOut);
                                }
                            break;
#endif
                        case OUT_IORDER_ARRAY:
                            pkdGenericSeek(pkd,fp, nStart,sizeof(int),sizeof(LongOut));
                            for (i=0;i<pkd->nLocal;++i) {
                                LongOut = pkd->pStore[i].iOrder;
                                xdr_long(&xdrs,&LongOut);
                                }
                            break;
                        default:
                            pkdGenericSeek(pkd,fp, N*iDim+nStart,sizeof(int),sizeof(DoubleOut));
                            for (i=0;i<pkd->nLocal;++i) {
                                DoubleOut = VecType(pkd, &pkd->pStore[i],iDim,iVecType);
                                xdr_double(&xdrs,&DoubleOut);
                                }
                        }
                    }
                break;
            case 3:
                if (iDim < 0) {
                    pkdGenericSeek(pkd,fp, nStart*3,sizeof(int),sizeof(fOut));
                    for (i=0;i<pkd->nLocal;++i) {
                        fOut = VecType(pkd, &pkd->pStore[i],0,iVecType);
                        xdr_FLOAT(&xdrs,&fOut);
                        fOut = VecType(pkd, &pkd->pStore[i],1,iVecType);
                        xdr_FLOAT(&xdrs,&fOut);
                        fOut = VecType(pkd, &pkd->pStore[i],2,iVecType);
                        xdr_FLOAT(&xdrs,&fOut);
                        }
                    } else {
                    switch (iVecType) {
#ifdef STARFORM
                        case OUT_IGASORDER_ARRAY:
                            pkdGenericSeek(pkd,fp, nStart, sizeof(int), sizeof(pkd->pStore[i].iGasOrder));
                            for (i=0;i<pkd->nLocal;++i) {
                                xdr_int(&xdrs,&(pkd->pStore[i].iGasOrder));
                                }
                            break;
#endif
                        case OUT_IORDER_ARRAY:
                            pkdGenericSeek(pkd,fp, nStart,sizeof(int),sizeof(pkd->pStore[i].iOrder));
                            for (i=0;i<pkd->nLocal;++i) {
                                xdr_int(&xdrs,&(pkd->pStore[i].iOrder));
                                }
                            break;
                        default:
                                pkdGenericSeek(pkd,fp, N*iDim+nStart,sizeof(int),sizeof(fOut));
                                for (i=0;i<pkd->nLocal;++i) {
                                        fOut = VecType(pkd, &pkd->pStore[i],iDim,iVecType);
                                        xdr_FLOAT(&xdrs,&fOut);
                                        }
                        }
                    }
                break;
                }
        }
    else {
        switch (iBinaryOutput) {
            case 0:
                if (iDim < 0) {
                    for (i=0;i<pkd->nLocal;++i) {
                        fOut = VecType(pkd, &pkd->pStore[i],0,iVecType);
                        fprintf(fp,"%.14g\n",fOut);
                        fOut = VecType(pkd, &pkd->pStore[i],1,iVecType);
                        fprintf(fp,"%.14g\n",fOut);
                        fOut = VecType(pkd, &pkd->pStore[i],2,iVecType);
                        fprintf(fp,"%.14g\n",fOut);
                        }
                    } else {
                    for (i=0;i<pkd->nLocal;++i) {
                        fOut = VecType(pkd, &pkd->pStore[i],iDim,iVecType);
                        fprintf(fp,"%.14g\n",fOut);
                        }
                    }
                    break;
            case 1:
                if (iDim < 0) {
                    pkdGenericSeek(pkd,fp, nStart*3,sizeof(int),sizeof(FloatOut));
                    for (i=0;i<pkd->nLocal;++i) {
                        FloatOut = VecType(pkd, &pkd->pStore[i],0,iVecType);
                        fwrite(&FloatOut, sizeof(float), 1, fp );
                        FloatOut = VecType(pkd, &pkd->pStore[i],1,iVecType);
                        fwrite(&FloatOut, sizeof(float), 1, fp );
                        FloatOut = VecType(pkd, &pkd->pStore[i],2,iVecType);
                        fwrite(&FloatOut, sizeof(float), 1, fp );
                        }
                    } else {
                    switch (iVecType) {
    #ifdef STARFORM
                        case OUT_IGASORDER_ARRAY:
                            pkdGenericSeek(pkd,fp, nStart,sizeof(int),sizeof(IntOut));
                            for (i=0;i<pkd->nLocal;++i) {
                                fwrite(&(pkd->pStore[i].iGasOrder), sizeof(IntOut), 1, fp );
                                }
                            break;
    #endif
                        case OUT_IORDER_ARRAY:
                            pkdGenericSeek(pkd,fp,nStart, sizeof(int), sizeof(IntOut));
                            for (i=0;i<pkd->nLocal;++i) {
                                fwrite(&(pkd->pStore[i].iOrder), sizeof(IntOut), 1, fp );
                                }
                            break;
                        default:
                            pkdGenericSeek(pkd,fp, N*iDim+nStart,sizeof(int),sizeof(FloatOut));
                            for (i=0;i<pkd->nLocal;++i) {
                                FloatOut = VecType(pkd, &pkd->pStore[i],iDim,iVecType);
                                fwrite(&FloatOut, sizeof(float), 1, fp );
                                }
                        }
                    }
                    break;
            case 2:
                if (iDim < 0) {
                    pkdGenericSeek(pkd,fp, nStart*3,sizeof(int),sizeof(DoubleOut));
                    for (i=0;i<pkd->nLocal;++i) {
                        DoubleOut = VecType(pkd, &pkd->pStore[i],0,iVecType);
                        fwrite(&DoubleOut, sizeof(double), 1, fp );
                        DoubleOut = VecType(pkd, &pkd->pStore[i],1,iVecType);
                        fwrite(&DoubleOut, sizeof(double), 1, fp );
                        DoubleOut = VecType(pkd, &pkd->pStore[i],2,iVecType);
                        fwrite(&DoubleOut, sizeof(double), 1, fp );
                        }
                        
                    } else {
                    switch (iVecType) {
    #ifdef STARFORM
                        case OUT_IGASORDER_ARRAY:
                            pkdGenericSeek(pkd,fp, nStart,sizeof(int),sizeof(LongOut));
                            for (i=0;i<pkd->nLocal;++i) {
                                LongOut = pkd->pStore[i].iGasOrder;
                                fwrite(&LongOut, sizeof(LongOut), 1, fp );
                                }
                            break;
    #endif
                        case OUT_IORDER_ARRAY:
                            pkdGenericSeek(pkd,fp, nStart,sizeof(int),sizeof(LongOut));
                            for (i=0;i<pkd->nLocal;++i) {
                                LongOut = pkd->pStore[i].iOrder;
                                fwrite(&LongOut, sizeof(LongOut), 1, fp );
                                }
                            break;
                        default:
                            pkdGenericSeek(pkd,fp, N*iDim+nStart,sizeof(int),sizeof(DoubleOut));
                            for (i=0;i<pkd->nLocal;++i) {
                                DoubleOut = VecType(pkd, &pkd->pStore[i],iDim,iVecType);
                                fwrite(&DoubleOut, sizeof(double), 1, fp );
                                }
                        }
                    }
                break;
            case 3:
                if (iDim < 0) {
                    pkdGenericSeek(pkd,fp, nStart*3,sizeof(int),sizeof(fOut));
                    for (i=0;i<pkd->nLocal;++i) {
                        fOut = VecType(pkd, &pkd->pStore[i],0,iVecType);
                        fwrite(&fOut, sizeof(FLOAT), 1, fp );
                        fOut = VecType(pkd, &pkd->pStore[i],1,iVecType);
                        fwrite(&fOut, sizeof(FLOAT), 1, fp );
                        fOut = VecType(pkd, &pkd->pStore[i],2,iVecType);
                        fwrite(&fOut, sizeof(FLOAT), 1, fp );
                        }
                    } else {
                    switch (iVecType) {
    #ifdef STARFORM
                        case OUT_IGASORDER_ARRAY:
                            pkdGenericSeek(pkd,fp, nStart, sizeof(int), sizeof(pkd->pStore[i].iGasOrder));
                            for (i=0;i<pkd->nLocal;++i) {
                                fwrite(&(pkd->pStore[i].iGasOrder), sizeof(pkd->pStore[i].iGasOrder), 1, fp );
                                }
                            break;
    #endif
                        case OUT_IORDER_ARRAY:
                            pkdGenericSeek(pkd,fp, nStart,sizeof(int),sizeof(pkd->pStore[i].iOrder));
                            for (i=0;i<pkd->nLocal;++i) {
                                fwrite(&(pkd->pStore[i].iOrder), sizeof(pkd->pStore[i].iOrder), 1, fp );
                                }
                            break;
                        default:
                                pkdGenericSeek(pkd,fp, N*iDim+nStart,sizeof(int),sizeof(fOut));
                                for (i=0;i<pkd->nLocal;++i) {
                                        fOut = VecType(pkd, &pkd->pStore[i],iDim,iVecType);
                                        fwrite(&fOut, sizeof(FLOAT), 1, fp );
                                        }
                        }
                    }
                break;
                }
            }

	i = fclose(fp);
	if (i != 0) {
		perror("pkdOutVector: could not close file");
		exit(1);
		}
	}
