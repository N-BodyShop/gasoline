#ifndef COSMO_HINCLUDED
#define COSMO_HINCLUDED

typedef struct csmContext {
    int bComove;	   
    double dHubble0;
    double dOmega0;
    double dLambda;
    double dOmegaRad;
    } * CSM;

void csmInitialize(CSM *pcsm);
double csmExp2Hub(CSM csm, double dExp);
double csmTime2Hub(CSM csm,double dTime);
double csmExp2Time(CSM csm,double dExp);
double csmTime2Exp(CSM csm,double dTime);
double csmComoveDriftInt(CSM csm, double dIExp);
double csmComoveKickInt(CSM csm, double dIExp);
double csmComoveDriftFac(CSM csm,double dTime,double dDelta);
double csmComoveKickFac(CSM csm,double dTime,double dDelta);

     
#endif
