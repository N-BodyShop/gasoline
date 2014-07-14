
#ifndef TREEZIP_HINCLUDED
#define TREEZIP_HINCLUDED

#include "treeziptypes.h"

/* very long integer emulation macros -- uses definitions above */
#include "treezipkey.h"

/* 
   Gives a safety valve if more than nPerBucket particles
   have identical keys 

   Probably should leave ON
*/
#define DEPTHCHECK

/* 
   Remove Label info in particle data

   Probably should leave OFF
*/
/* #define NOLABEL */
typedef unsigned int LABELTYPE;

/* debug info */
/* #define DEBUG  */

/* very verbose debug info */
/* #define DEBUG2 */

typedef struct {
	tzkey k;
	LABELTYPE label;
#ifdef DEBUG
	float pos[3];
#endif
	} tzparticle;

/* should be multiple of nPerBucket */
#define TZNPLIST 32768

typedef struct tzParticleBlock {
	tzparticle particledata[TZNPLIST];
	void *nextblock;
	} tzparticleblock;

typedef struct tzNodeStruct {
	struct tzNodeStruct *child;
	int n;
	tzparticle *p;
	} tznode;

/* should be multiple of 2 */
#define TZNCNODE 32768

typedef struct tzNodeBlock {
	tznode nodedata[TZNCNODE];
	void *nextblock;
	} tznodeblock;

typedef struct TreeZipContext {
	/* Tree properties */
	int nBits_nPerBucket;
	int nPerBucket ;
	int iMaskBucket;
	int nBits_Position ;
	int iMaskPosition ;
	int nBits_PerDirection;
	TZ_UINT64 iBigInt;
	int nBitsMaxPrecision;
	int nBitsMinPrecision;
	int nBitsMaxPrecision_PerDirection;
	int nBitsMinPrecision_PerDirection;
	TZ_UINT16 nBitsExtra_PerParticle;

	int nParticle;
	int nNode;
	int nBucket;
#ifdef DEPTHCHECK
	int maxdepth;
#endif

	/* Tree */
	tznode root;

	/* Dataset info */
	double dmin[3], dmax[3], dIsize[3];
	
	/* Node data */
	tznodeblock *cnodeblock;
	tznode *cnode;
	int icnode;
	
	/* Particle DATA */
	tzparticleblock *pblock;
	tzparticle *plist;
	int iplist;

	/* Output variables */
	FILE *fpout;
	TZ_UINT64 bitstream;
	int nbits;
	int nTotalbits;
	int nWritebits;
	int nParticlebits;
	int nLabelbits;

  } TZX;


/* Function prototypes */

TZX *tzInit( double dmin[3], double dmax[3], int nBits_nPerBucket, int nBits_Position, int nMaxBitsPrecision, int nMinBitsPrecision );

void tzFinalize( TZX *tz );

void tzEmptyTree( TZX *tz );

tznode *tzNodeAllocation( TZX *tz, int n );

tzparticle *tzParticleAllocation( TZX *tz, int n );

void tzWriteBits( TZX *tz, int bits, int nbitwrite );

void tzWriteNode( TZX *tz, tznode *c, int l );

void tzWriteTreeZip( TZX *tz );

void tzAddPos( TZX *tz, double *r, LABELTYPE label );

void tzOutputFile( TZX *tz, FILE *fpout );

void tzWriteHeader( TZX *tz );

#endif
