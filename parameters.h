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
	int bEpsVel;
	int bSqrtPhi;
	int bISqrtRho;
	int bNonSymp;
	int bDoDensity;
	int bDohOutput;
	int bDoIonOutput;
	int bSymCool;
	int bDoGravity;
	int bFandG;
	int bHeliocentric;
	int bLogHalo;
	int bHernquistSpheroid;
	int bMiyamotoDisk;
	int nBucket;
	int iOutInterval;
	int iLogInterval;
	int iCheckInterval;
	int iOrder;
	int iEwOrder;
	int nReplicas;
        int iStartStep;
	int nSteps;
	int nSmooth;
	int iMaxRung;
	int nSuperCool;
	int nGrowMass;
        int iWallRunTime;
	double dEta;
	double dExtraStore;
	double dSoft;
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
	/*
	 ** Additional parameters for GASOLINE.
	 */
	int bGeometric;
        int iGasModel;
	double dEtaCourant;
	double dEtauDot;
        double duDotLimit;
	double dConstAlpha;
	double dConstBeta;
	double dConstGamma;
	double dMeanMolWeight;
	double dGasConst;
	double dMsolUnit;
	double dKpcUnit;
	double dGmPerCcUnit;
	double dComovingGmPerCcUnit;
	double dErgPerGmUnit;
	double dSecUnit;
        double dMassFracHelium;
        double dCoolingTmin;
        double dCoolingTmax;
        int    nCoolingTable;
        int    bViscosityLimiter;
        int    bBulkViscosity;
        int    bDoGas;
        int    bUV;
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
	double dSmallStep;
	int iOutcomes;
	double dEpsN;
	double dEpsT;
	int bDoCollLog;
	char achCollLog[256];
#endif /* COLLISIONS */
	};


#endif




