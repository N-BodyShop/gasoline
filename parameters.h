#ifndef PARAMETERS_HINCLUDED
#define PARAMETERS_HINCLUDED

#include "cosmo.h"

struct parameters {
	/*
	 ** Parameters for PKDGRAV.
	 */
	int nThreads;
	int bDiag;
	int bOverwrite;
	int bVWarnings;
	int bVStart;
	int bVStep;
	int bVRungStat;
	int bVDetails;
	int bPeriodic;
	int bRestart;
	int bParaRead;
	int bParaWrite;
	int bCannonical;
	int bStandard;
	int bKDK;
	int bBinary;
	int bGravStep;
	int bEpsAccStep;
	int bSqrtPhiStep;
	int bAccelStep; /* true if bEpsAccStep or bSqrtPhiStep */
	int bDensityStep;
	int nTruncateRung;
	int bNonSymp;
	int bDoDensity;
 	int bDohOutput;
 	int bDoSphhOutput;
	int bDodtOutput;
	int bDoIonOutput;
	int bSymCool;
	int bDoGravity;
	int bFandG;
	int bHeliocentric;
	int bLogHalo;
	int bHernquistSpheroid;
	int bMiyamotoDisk;
	int bRotFrame;
	double dOmega;
	double dOmegaDot;
	int bPatch;
	double dOrbFreq;
	int nBucket;
	int iOutInterval;
	int iLogInterval;
	int iCheckInterval;
	int iOrder;
	int bEwald;
	int iEwOrder;
	int nReplicas;
	int iStartStep;
	int nSteps;
	int nSmooth;
	int iMaxRung;
	int nSuperCool;
	int nGrowMass;
	int iWallRunTime;
        int bPhysicalSoft;  
        int bSoftMaxMul;
        int bVariableSoft;
        int nSoftNbr;
        int bSoftByType;
        int bDoSoftOutput;
	double dEta;
	double dExtraStore;
	double dSoft;
        double dSoftMax;
	double dDelta;
	double dEwCut;
	double dEwhCut;
	double dTheta;
	double dTheta2;
	double daSwitchTheta;
	double dAbsPartial;
	double dRelPartial;
	double dAbsTotal;
	double dRelTotal;
	double dPeriod;
	double dxPeriod;
	double dyPeriod;
	double dzPeriod;
	CSM csm;
	double dRedTo;
	double dCentMass;
	char achDigitMask[256];
	char achInFile[256];
	char achOutName[256];
	char achDataSubPath[256];
	double dCoolFac;
	double dCoolDens;
	double dCoolMaxDens;
	double dGrowDeltaM;
	double dGrowStartT;
	double dGrowEndT;
	double dFracNoDomainDecomp;
	double dFracNoDomainDimChoice;
	int    bRungDD;
	double dRungDDWeight;
	/*
	 ** Additional parameters for GASOLINE.
	 */
	int bGeometric;
	int iGasModel;
	double dEtaCourant;
	double dEtauDot;
	double duDotLimit;
        double dShockTrackerA;
        double dShockTrackerB;
	double dConstAlpha;
	double dConstBeta;
	double dConstGamma;
	double dMeanMolWeight;
	double dGasConst;
	double dMsolUnit;
	double dKpcUnit;
	double ddHonHLimit;
	double dGmPerCcUnit;
	double dComovingGmPerCcUnit;
	double dErgPerGmUnit;
	double dSecUnit;
	double dMassFracHelium;
	double dCoolingTmin;
	double dCoolingTmax;
	int    nCoolingTable;
	int    bViscosityLimiter;
	int    bViscosityLimitdt;
        int    bShockTracker;
	int    bBulkViscosity;
	int    bGasDomainDecomp;
        int    bLowerSoundSpeed;
	int    bFastGas;
	double dFracFastGas;
	double dhMinOverSoft;
	int    bDoGas;
	int    bSphStep;
	int    bUV;
	int    bSN;
	double dSNRhoCut;
 	double dSNTMin;
        double dSNTMax;
	double dSNMetalCut;
	double dSNHeatFraction;
	int    bStarForm;
	int    bFeedBack;
#ifdef STARFORM
	STFM   stfm;
	FB     fb;
#endif
#ifdef GLASS
	/*
	 ** Additional parameters for GLASS.
	 */
	double dGlassDamper;
	/* Hack for shock tube glasses */
	double dGlassPoverRhoL;
	double dGlassPoverRhoR;
	double dGlassxL;
	double dGlassxR;
	double dGlassVL;
	double dGlassVR;
#endif
#ifdef COLLISIONS
	/*
	 ** Additional parameters for collision code...
	 */
	int bFindRejects;
	int bDoCollLog;
	char achCollLog[256];
	double dSmallStep;
	double dxUnifGrav;
	double dyUnifGrav;
	double dzUnifGrav;
	COLLISION_PARAMS CP;
#endif /* COLLISIONS */
	};

#endif
