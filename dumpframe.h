#ifndef DUMPFRAME_HINCLUDED
#define DUMPFRAME_HINCLUDED

#define DF_NXPIXMAX 1024
#define DF_NYPIXMAX 1024
/* PST */

#ifdef USE_PNG
#include "png.h"        /* libpng header; includes zlib.h and setjmp.h */
#include "writepng.h"   /* typedefs, common macros, public prototypes */
#endif

typedef struct dfImage {
	float r,g,b;
	} DFIMAGE;

#define DF_NBYTEDUMPFRAME (3*sizeof(DFIMAGE)*DF_NXPIXMAX*DF_NYPIXMAX)

/* 
   Projection can be 2D or voxels

   In principle you can render voxels
   and encode them in different ways.
   I will just build a treezip every time.
   */

enum df_dimension {
	DF_2D, /* Generic 2D */
	DF_3D /* Voxel */
	};

enum df_projectstyle {
	DF_PROJECT_NULL,
	DF_PROJECT_ORTHO,
	DF_PROJECT_PERSPECTIVE
	};

/* in principle voxels can have encoding options 
   For now it will be treezip 
   */
enum df_encodestyle {
	DF_ENCODE_NULL,
	DF_ENCODE_PPM,
	DF_ENCODE_PNG,
	DF_ENCODE_RLE,
	DF_ENCODE_TREEZIP
	};

/* in principle voxels can have rendering options 
   For now this is ignored 
   */
enum df_renderstyle {
	DF_RENDER_NULL,
	DF_RENDER_POINT,
	DF_RENDER_TSC,
	DF_RENDER_SOLID,
	DF_RENDER_SHINE
	};

enum df_numbering {
	DF_NUMBERING_FRAME,
	DF_NUMBERING_STEP,
	DF_NUMBERING_TIME
	};

/* PST */
/* There is a common subset with framesetup that needs it's own structure -- parent class even */
struct inDumpFrame {
	double dTime;
	double dStep;
	/* 2D Projection info */
	double r[3]; /* Centre */
	double x[3]; /* Projection vectors */
	double y[3];
	double z[3];
	double zClipNear,zClipFar,zEye;
	int bPeriodic; /* Is it periodic? */
	int nxRepNeg,nxRepPos; /* Replicas for Periodic */
	int nyRepNeg,nyRepPos;
	int nzRepNeg,nzRepPos;
	int nxPix,nyPix;    /* Pixmap dimensions */
	int iProject;      /* Projection */

	/* Rendering */
	double pScale1,pScale2;
	double dGasSoftMul,dDarkSoftMul,dStarSoftMul;
	double xlim,ylim,hmul;
	DFIMAGE ColStar,ColGas,ColDark;
	int bColMassWeight;
	int bLogScale;
	int bGasSph;
	int iRender;       
	/* Render Particle interface */
	int offsetp_r,offsetp_fMass,offsetp_fSoft,offsetp_fBall2,offsetp_iActive;
	int sizeofp;
	int iTypeGas,iTypeDark,iTypeStar;
	
	double dMassGasMin, dMassGasMax;
	double dMassDarkMin,dMassDarkMax;
	double dMassStarMin,dMassStarMax;

	int bNonLocal; /* Is this data going non-local? */
	int bVDetails;
    };

/* This is a table entry for later interpolation, 
   changing properties of frames */
struct dfFrameSetup {
	double dTime;
	/* Projection info */
	double target[3]; /* Centre */
	double eye[3]; /* Eye Position */
	double up[3]; /* up vector */
	double FOV;
	double zEye1,zEye2,zEye;
	double zClipNear,zClipFar; /* clipping */
    double eye2[3]; /* a second vector to add to the eye vector */
	int bEye2;
	int bzClipFrac; /* Is z clipping a fraction of eye to target distance? */
	int bzEye,bzEye1,bzEye2; /* Is zEye a fixed distance? */
	int bEyeAbsolute; /* Eye position is in absolute coordinates (default relative to target point) */
	int bAnchor;
	int nxPix,nyPix;    /* Pixmap dimensions */
	int bPeriodic; /* Is it periodic? */
	int iProject;      /* Projection */

	/* Rendering controllers */
	double pScale1,pScale2;
	DFIMAGE ColStar,ColGas,ColDark;
	int bColMassWeight;
	int bLogScale;
	int bGasSph;
	double dGasSoftMul,dDarkSoftMul,dStarSoftMul;
	int iRender;       /* Rendering */
        


	int bNonLocal; /* Is this data going non-local? */
    };

/* MSR */

struct DumpFrameContext {
	int bAllocated; /* Was this malloc'ed? */

	int nFrame;
    int iMaxRung;	
	double dStep;
	double dTime;
	double dDumpFrameStep;
	double dDumpFrameTime;

	/* Particle Filters */
	double dMassGasMin, dMassGasMax;
	double dMassDarkMin,dMassDarkMax;
	double dMassStarMin,dMassStarMax;

	/* Time dependent Frame setup data */
	int iFrameSetup;
	int nFrameSetup;
	struct dfFrameSetup *fs;
	struct dfFrameSetup a;
	struct dfFrameSetup b;
	struct dfFrameSetup c;
	struct dfFrameSetup d;
	double rdt;
	double dTimeMod;
	double dTimeLoop,dPeriodLoop;
	int bLoop;

	int iDimension; /* 2D Pixel or 3D Voxel */
	int iEncode;
	int iNumbering;
	int bVDetails;
	char FileName[161];
	};

void dfInitialize( struct DumpFrameContext **pdf, double dTime, 
				  double dDumpFrameTime, double dStep, double dDumpFrameStep,
				  double dDelta, int iMaxRung, int bVDetails, char* );
void dfFinalize( struct DumpFrameContext *df );

void *dfAllocateImage( int nxPix, int nyPix );
void dfFreeImage( void *Image );

void dfParseOptions( struct DumpFrameContext *df, char * filename );

void dfParseCameraDirections( struct DumpFrameContext *df, char * filename );

void dfSetupFrame( struct DumpFrameContext *df, double dTime, double dStep, struct inDumpFrame *in );

void dfMergeImage( struct inDumpFrame *in, void *Image1, int *nImage1, void *Image2, int *nImage2 );
void dfClearImage( struct inDumpFrame *in, void *Image, int *nImage );

void dfRenderParticlePoint( struct inDumpFrame *in, void *vImage, 
						   double *r, double fMass, double fSoft, double fBall2, int iActive );
void dfRenderParticleTSC( struct inDumpFrame *in, void *vImage, 
						 double *r, double fMass, double fSoft, double fBall2, int iActive );
void dfRenderParticleSolid( struct inDumpFrame *in, void *vImage, 
						 double *r, double fMass, double fSoft, double fBall2, int iActive );
void dfRenderParticleInit( struct inDumpFrame *in, int iTypeGas, int iTypeDark, int iTypeStar );
void dfRenderParticle( struct inDumpFrame *in, void *vImage, 
					  double *r, double fMass, double fSoft, double fBall2, int iActive );

void dfRenderParticlesInit( struct inDumpFrame *in, int iTypeGas, int iTypeDark, int iTypeStar,
						   double *pr, double *pfMass, double *pfSoft, double *pfBall2, unsigned int *piActive, 
						   void *p, int sizeofp );
void dfRenderParticles( struct inDumpFrame *in, void *vImage, void *pStore, int n );

/* Relies on PKD definition -- not portable */
#ifdef PKD_HINCLUDED
void dfRenderImage( PKD pkd, struct inDumpFrame *in, void *vImage );
void dfRenderImageOld( PKD pkd, struct inDumpFrame *in, void *Image );
#endif

void dfFinishFrame( struct DumpFrameContext *df, double dTime, double dStep, struct inDumpFrame *in, void *Image );

/* Interpolation Functions */

void dfGetCoeff4( struct DumpFrameContext *df, int ifs );
void dfGetCoeff3L( struct DumpFrameContext *df, int ifs );
void dfGetCoeff3R( struct DumpFrameContext *df, int ifs );
void dfGetCoeff2( struct DumpFrameContext *df, int ifs );

void dfGetCoeff( struct DumpFrameContext *df, int ifs );

void dfInterp( struct DumpFrameContext *df, struct dfFrameSetup *pfs, double x );

/* #include "dumpvoxel.h" */

#endif
