#ifndef OUTURB_HINCLUDED
#define OUTURB_HINCLUDED

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

typedef struct OUturbparam {
//Input Parameters
    double StDecay;      // decay time
    double StEnergy;     // target vel disp?
    double StDtFreq;     // frequency of updates to turb stirring field
    double StKmin;       // k range
    double StKmax;
    double StSolWeight;  // solenoidal weighting 0-1 (for projection of driver)
    double StAmplFac;    // target vel disp?
    double StScaleHeight; // Gaussian fall off for forcing
    double StStartTime;  // starting time for OU sequence 
    int StSpectForm;
    int StSeed;
    } OUTURBPARAM;

typedef struct OUturb {
//Input Parameters
    double StDecay;      // decay time
    double StEnergy;     // target vel disp?
    double StDtFreq;     // frequency of updates to turb stirring field
    double StKmin;       // k range
    double StKmax;
    double StSolWeight;  // solenoidal weighting 0-1 (for projection of driver)
    double StAmplFac;    // target vel disp?
    double StScaleHeight; // Gaussian fall off for forcing
    int StSpectForm;
    int StSeed;

//Proc Id
    int idSelf;
    
//Ornstein-Uhlenbeck variables
    double StOUVar;
    double *StOUPhases;
    gsl_rng *StRng;

//forcing field in fourier space
    double *StAmpl;
    double *StAka; //phases (real part)
    double *StAkb; //phases (imag part)
    double *StMode;
    int StNModes;

//    int StTPrev; // Something like step number
    double dTimePrev; // Gasoline uses time directly
    double StSolWeightNorm;

/*double StEnergyAcc;
double StEnergyDeacc;

double StLastStatTime;

FILE *FdTurb;*/

    } * OUTURB;

struct inInitouturb {
    OUTURBPARAM outurbparam;
    double BoxSize;
    double dTime;
    int bDetails;
    int bRestart;
    };

struct inAccelouturb {
    double dTime;
    };

void outurb_AddParams( OUTURBPARAM *outurbparam, PRM prm );
void outurbLogParams( OUTURBPARAM *outurbparam, FILE *fp );

void outurb_init(OUTURB *pouturb, OUTURBPARAM outurbparam, int idSelf, int bDetails, int bRestart, double BoxSize, double dTime);
void outurb_st_init_ouseq(OUTURB outurb);
void outurb_st_update_ouseq(OUTURB outurb);
double outurb_st_grn(OUTURB outurb);
void outurb_st_calc_phases(OUTURB outurb);
void outurb_set_turb_ampl(OUTURB outurb, double dTime);
void outurb_add_turb_accel(OUTURB outurb, double dTime, PARTICLE *p, int nParticle);
void outurb_driving_step_first_half(OUTURB outurb, double dTime, double dDelta, PARTICLE *p, int nParticle);
void outurb_driving_step_second_half(OUTURB outurb, double dTime, double dDelta, PARTICLE *p, int nParticle);

#endif
