#ifndef SPECIAL_HINCLUDED
#define SPECIAL_HINCLUDED

#ifdef SPECIAL_PARTICLES

#define MAX_NUM_SPECIAL_PARTICLES 3

/* bit masks */

#define SPECIAL_OBLATE	(1<<0)
#define SPECIAL_GR		(1<<1) /* not implemented */

struct oblate_s {
	double dRadEq;	/* equatorial radius (system units) */
	double J2,J4;	/* J2, J4 oblateness coefficients */
	};

typedef struct {
	int iType;
	struct oblate_s oblate;
	} SPECIAL_PARTICLE_DATA;

typedef struct {
	int iOrder;
	FLOAT fMass;
	FLOAT r[3];
	} SPECIAL_PARTICLE_INFO;

void pkdGetSpecialParticles(PKD pkd, int nSpecial, int iId[],
							FLOAT fCentMass, SPECIAL_PARTICLE_INFO sInfo[]);
void pkdDoSpecialParticles(PKD pkd, int nSpecial, int bNonInertial,
						   SPECIAL_PARTICLE_DATA sData[],
						   SPECIAL_PARTICLE_INFO sInfo[], FLOAT aFrame[]);

#endif /* SPECIAL_PARTICLES */

#endif
