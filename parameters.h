#ifndef PARAMETERS_HINCLUDED
#define PARAMETERS_HINCLUDED


struct parameters {
	int nThreads;
	int bDiag;
	int bVerbose;
	int bPeriodic;
	int bGatherScatter;
	int bRestart;
	int bComove;
	int bParaRead;
	int bParaWrite;
	int bCannonical;
	int bStandard;
	int bKDK;
	int bBinary;
	int nBucket;
	int iOutInterval;
	int iLogInterval;
	int iCheckInterval;
	int iOrder;
	int iEwOrder;
	int nReplicas;
	int nSteps;
	int nSmooth;
	double dExtraStore;
	double dSoft;
	double dDelta;
	double dEwCut;
	double dEwhCut;
	double dTheta;
	double dAbsPartial;
	double dRelPartial;
	double dAbsTotal;
	double dRelTotal;
	double dPeriod;
	double dHubble0;
	double dOmega0;
	double dRedTo;
	char achInFile[256];
	char achOutName[256];
	char achDataSubPath[256];
	};


#endif



