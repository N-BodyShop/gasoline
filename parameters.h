#ifndef PARAMETERS_HINCLUDED
#define PARAMETERS_HINCLUDED


struct parameters {
	int nThreads;
	int bDiag;
	int bVerbose;
	int bPeriodic;
	int bGatherScatter;
	int bRestart;
	int nBucket;
	int iOutInterval;
	int iLogInterval;
	int iCheckInterval;
	int iOrder;
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
	char achInFile[256];
	char achOutName[256];
	char achDataSubPath[256];
	double dHubble0;
	double dOmega0;
	int bComove;
	};


#endif

