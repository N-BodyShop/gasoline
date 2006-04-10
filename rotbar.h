#ifndef ROTBAR_HINCLUDED
#define ROTBAR_HINCLUDED

#include "fdl.h"
#include "param.h"

typedef struct rotbarContext {
    double bratio;
    double cratio;
  double amplitude; /* dAmp times dAmpFac (gets passed to pkdRotatingBar() */
  double dAmplitude; /* relative amplitude */
    double Fcorot;
    int soft;
    double dPos[3];
    double dVel[3];
    double dAcc[3];
    double dTorque[3];
    double (*getMass)(struct rotbarContext *, double);
    double (*getPot)(struct rotbarContext *, double);
    double dMass;	/* parameters for rotating bar */
    double dLength;	
    double dCorotFac;
    double dTurnOff;
    double dTurnOn;
    double dAmpFac;
    double dDuration;	/* time scale for Turn Off/On */
    double dPosAng;	/* current angle */
    double dTime0;	/* time at beginning of big timestep */
    double dOmega;
    double dB5;
    double dIz;
    double dLz;
    double dLz0;		/* Initial Angular moment of bar */
    double dLzPart;		/* Angular momentum of particles */
    double *pdMmono;
    double *pdRmono;
    double *pdPmono;
    int nBinMono;
    double A;
    double B;
    double C;
    int bFixedBar;
    int bMonopole;
    } *ROTBAR;

void rotbarAddParams(ROTBAR rotbar, PRM prm);
void rotbarLogParams(ROTBAR rotbar, FILE *fp );
void rotbarCheckWrite(ROTBAR rotbar, FDL_CTX *fdl);
void rotbarCheckRead(ROTBAR rotbar, FDL_CTX *fdl);
void rotbarInitialize(ROTBAR *protbar);
void rotbarInitValues(ROTBAR rotbar);
void rotbarDrift(ROTBAR rotbar, double dTime, double dDelta);
void rotbarKick(ROTBAR rotbar, double dvFacOne, double dvFacTwo);
#endif
