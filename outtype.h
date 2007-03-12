#ifndef OUTTYPE_HINCLUDED
#define OUTTYPE_HINCLUDED

#include "pkd.h"

#define NUMOUTPUTS  128

enum DataTypeCode {
	INT8 = 1,
	UINT8,
	INT16,
	UINT16,
	INT32,
	UINT32,
	INT64,
	UINT64,
	FLOAT32,
	FLOAT64
};

enum outtype_arraytype {
	OUT_NULL,
        BIG_FILE,
	OUT_COLOR_ARRAY,
	OUT_DENSITY_ARRAY,
	OUT_DENSITYU_ARRAY,
	OUT_PRES_ARRAY,
	OUT_POT_ARRAY,
	OUT_AMAG_ARRAY,
	OUT_IMASS_ARRAY,
        OUT_MASS_ARRAY,
	OUT_RUNG_ARRAY,
	OUT_SPHH_ARRAY,
	OUT_U_ARRAY,
        OUT_TEMP_ARRAY,
        OUT_GASDENSITY_ARRAY,
	OUT_HSMDIVV_ARRAY,
	OUT_UDOT_ARRAY,
	OUT_COOL_ARRAY0,
	OUT_COOL_ARRAY1,
	OUT_COOL_ARRAY2,
	OUT_BALSARASWITCH_ARRAY,
	OUT_DIVV_ARRAY,
	OUT_MUMAX_ARRAY,
	OUT_SHOCKTRACKER_ARRAY,
	OUT_DIVONCONH_ARRAY,
	OUT_DIVONCONX_ARRAY,
	OUT_DIVRHOV_ARRAY,
	OUT_PDV_ARRAY,
	OUT_PDVPRES_ARRAY,
	OUT_PDVVISC_ARRAY,
	OUT_PDVSN_ARRAY,
	OUT_USN_ARRAY,
	OUT_METALS_ARRAY,
	OUT_IGASORDER_ARRAY,
        OUT_TIMEFORM_ARRAY,
        OUT_TEMPFORM_ARRAY,
	OUT_MASSFORM_ARRAY,
	OUT_DENSITYFORM_ARRAY,
	OUT_COOLTURNONTIME_ARRAY,
        OUT_OXYGENMASSFRAC_ARRAY,
        OUT_IRONMASSFRAC_ARRAY,
	OUT_TCOOLAGAIN_ARRAY,
	OUT_MSTAR_ARRAY,
	OUT_SOFT_ARRAY,
	OUT_IORDER_ARRAY,
	OUT_H_ARRAY,
        OUT_SPHDT_ARRAY,
	OUT_DT_ARRAY,
	OUT_REJECTS_ARRAY,
	OUT_TOFF_YR_ARRAY,
	OUT_TCOOL_YR_ARRAY,
	OUT_TDYN_YR_ARRAY,
	OUT_RATIOSOUNDDYN_ARRAY,
	OUT_L_JEANS_ARRAY,
	OUT_ISMALL_JEANS_ARRAY,
        OUT_1D3DSPLIT,  /* NOTICE!!
                       * Everything above here is 1D 
                       * Everything below here is 3D
                       */
	OUT_POS_VECTOR,
	OUT_VEL_VECTOR,
	OUT_ACCEL_VECTOR,
        OUT_ACCELG_VECTOR,
	OUT_VPRED_VECTOR,
	OUT_GRADRHO_VECTOR,
	OUT_ACCELPRES_VECTOR,
	OUT_RFORM_VECTOR,
	OUT_VFORM_VECTOR
	};

/*void pkdOutArray(PKD pkd,char *pszFileName,int nStart, int iArrType, int iBinaryOutput);*/
void VecFilename(char *achFile, int iType);
void pkdOutVector(PKD pkd,char *pszFileName,int nStart, int iDim,int iVecType,int iBinaryOutput, int N, int bStandard);
void pkdGenericSeek(PKD pkd,FILE *fp,int nStart,int iHeader, int iElement);
void pkdOutNChilada(PKD pkd,char *pszFileName,int nGasStart, int nDarkStart, int nStarStart, int iVecType, float minValue[3][3], float maxValue[3][3], double duTFac, double dvFac);
#endif
