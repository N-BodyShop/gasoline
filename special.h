#ifndef SPECIAL_HINCLUDED
#define SPECIAL_HINCLUDED

#ifdef SPECIAL_PARTICLES

#define MAX_NUM_SPECIAL_PARTICLES 3

/* bit masks */

#define SPECIAL_OBLATE	(1<<0)
#define SPECIAL_GR		(1<<1) /* not implemented */
#define SPECIAL_FORCE	(1<<2)
#define SPECIAL_NOGHOSTPERT	(1<<3)	/* No seeds in ghost cells */

struct oblate_s {
	double dRadEq;	/* equatorial radius (system units) */
	double J2,J4;	/* J2, J4 oblateness coefficients */
	double p[3];	/* pole of the oblate object */
	};

struct force_s {
	double dMag;	/* magnitude of force in system units */
	};

typedef struct {
	int iType;
	struct oblate_s oblate;
	struct force_s force;
	} SPECIAL_PARTICLE_DATA;

typedef struct {
        int bNonInertial;
        double dCentMass;
        double dxPeriod;
        double dyPeriod;
        double dzPeriod;
        double dOmega;
        double dTime;
        int nReplicas;
        } SPECIAL_MASTER_INFO;

typedef struct {
	int iOrder;
	FLOAT fMass,fRadius;
	FLOAT r[3];
	} SPECIAL_PARTICLE_INFO;

void pkdGetSpecialParticles(PKD pkd, int nSpecial, int iId[],
							SPECIAL_MASTER_INFO *mInfo, SPECIAL_PARTICLE_INFO sInfo[]);
void pkdDoSpecialParticles(PKD pkd, int nSpecial, SPECIAL_MASTER_INFO *mInfo,
						   SPECIAL_PARTICLE_DATA sData[],
						   SPECIAL_PARTICLE_INFO sInfo[], FLOAT aFrame[]);

#endif /* SPECIAL_PARTICLES */

#endif
