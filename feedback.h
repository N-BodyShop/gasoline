#ifndef FEEDBACK_HINCLUDED
#define FEEDBACK_HINCLUDED
#include "pkd.h"

struct snContext;

/*
 * Structure to Characterize Star Formation event.  This is input
 * needed for all feedback effects.
 */
typedef struct sfevent {
    double dMass;            /* mass in star formation event in solar masses */
    double dTimeForm;      /* time of star formation event in years */
    double dMetals;           /*  metallicity of stars in event */
    double dMFracOxygen;           /*  metallicity of stars in event */
    double dMFracIron;           /*  metallicity of stars in event */
    } SFEvent;

/*
 * Structure to return feedback effects.
 */
typedef struct fbeffects {
    double dEnergy;		/* Energy produced (ergs); note that
				   this might need to be divided up
				   into different forms of energy. */
    double dMassLoss;		/* Mass lost in Solar Masses */
    double dMetals;		/* Fraction of the mass lost in
				   elements heavier than Helium */
    double dMIron;              /* Solar masses of iron ejected */
    double dMOxygen;            /* Solar masses of oxygen ejected */
    } FBEffects;

#define FB_SNII 0
#define FB_SNIA 1
#define FB_WIND 2
#define FB_UV 	3
#define FB_NFEEDBACKS 4

typedef struct fbContext 
{
    double dGmUnit;		/* system mass in grams */
    double dSecUnit;		/* system time in seconds */
    double dErgPerGmUnit;	/* system specific energy in ergs/gm */
	double dInitStarMass; 
    }  * FB;

void fbInitialize(FB *pfb);

void pkdFeedback(PKD pkd, FB fb, struct snContext * sn, double dTime, double dDelta,
		 FBEffects *fbTotals);

	/* solar mass in grams */
#define MSOLG 1.99e33
        /* seconds per year = 3600*24*365.25 */
#define SEC_YR 3.15576e07;

#endif

    
