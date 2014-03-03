#include "define.h"
/*
 * This alternative dumpframe was created to:
 * 1) allow coloring by all particle properties
 * 2) change brightness scaling to be be a number that 
 *   a| increases with brightness
 *   b| is in the range 0-1
 * 
 * To allow access to more properties, I removed a lot of the
 * portability features, so that now I pass PKDs and PARTICLEs 
 * into the dfRender_____ functions.
 * 
 * The brightness scaling now relies on the dMinGasMass.
 */

#ifdef GSS_DUMPFRAME
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>

/* #define DEBUGTRACK 10  */

#ifdef DEBUGTRACK
int trackp = 0;
int trackcnt = 0;
int ntrack = 0, ntrack2 = 0;
#endif

/* Using pkd.h deprecated for greater portability */
#include "pkd.h"
#include "outtype.h"

#ifndef DBL_MAX
#define DBL_MAX 1.7976931348623157E+308
#endif

#include "dumpframe.h"

void dfInitialize( struct DumpFrameContext **pdf, double dYearUnit, 
		   double dTime,double dDumpFrameTime,double dStep,
		   double dDumpFrameStep,double dDelta,int iMaxRung, 
		   int bVDetails, char* filename) {
    double tock=0.0;/*initialized to suppress warning: DCR 12/19/02*/
    struct DumpFrameContext *df;

    if (*pdf!=NULL) (*pdf)->fs = NULL;
    if (dDumpFrameTime <= 0 && dDumpFrameStep <=0) return;

    if (*pdf == NULL) {
	*pdf = (struct DumpFrameContext *) malloc(sizeof(struct DumpFrameContext ));
	(*pdf)->bAllocated = 1;
	}
    else {
	(*pdf)->bAllocated = 0;
	}

    df = *pdf;
    df->iMaxRung = 0;
    df->dYearUnit = dYearUnit;

    df->dTime = 0;
    df->dDumpFrameTime = 0;
    if (dDumpFrameTime > 0) { 
	df->dDumpFrameTime = dDumpFrameTime;
	df->dTime = dTime-dDumpFrameTime*1e-10;
	/* Round to a nearby power of two */
	tock = -log((dDumpFrameTime-floor(dDumpFrameTime))/dDelta*0.50001+1e-60)
	    /log(2.0);
	if (tock <= iMaxRung) {
	    df->iMaxRung = tock;
	    }
	}

    df->dStep = 0;
    df->dDumpFrameStep = 0;
    if (dDumpFrameStep > 0) { 
	df->dDumpFrameStep = dDumpFrameStep;
	df->dStep = dStep-dDumpFrameStep*1e-10;
	/* Round to a nearby power of two */
	tock = -log((dDumpFrameStep-floor(dDumpFrameStep))*0.50001+1e-60)
	    /log(2.0);
	if (tock <= iMaxRung && tock >df->iMaxRung) {
	    df->iMaxRung = floor(tock);
	    }
	}
    /* General options */
    df->bVDetails = bVDetails;
    dfParseOptions( *pdf, filename );

    /* Parse time dependent camera options (mostly 2D) */
    dfParseCameraDirections( *pdf, filename );

    if (df->bVDetails) printf("DF Initialized Frame dumping: Time Interval %g [%g] (Step %i) Step Interval %g -> MaxRung %i\n",dDumpFrameTime,floor(dDumpFrameTime)+dDelta*pow(2.0,-floor(tock)),(int) floor(dDumpFrameTime/dDelta),dDumpFrameStep,df->iMaxRung);
    }


void dfFinalize( struct DumpFrameContext *df ) {
    if (df != NULL) {
	if (df->fs != NULL) free( df->fs );
	if (df->bAllocated) free( df );
	}
    }

void *dfAllocateImage( int nxPix, int nyPix ) {
    DFIMAGE *Image;

    assert ( nxPix*nyPix*sizeof(DFIMAGE) <= DF_NBYTEDUMPFRAME );

    Image = (DFIMAGE *) malloc( nxPix*nyPix*sizeof(DFIMAGE) );
    assert (Image != NULL );
	

    return Image;
    }

void dfFreeImage( void *Image ) {
    free( Image );
    }


#define SET( a, b ) {				\
	a[0]=b[0];				\
	a[1]=b[1];				\
	a[2]=b[2];				\
	} 
#define ADD( a, b ) {				\
	a[0]+=b[0];				\
	a[1]+=b[1];				\
	a[2]+=b[2];				\
	} 
#define DIFF( a, b, c ) {			\
	c[0] = b[0]-a[0];			\
	c[1] = b[1]-a[1];			\
	c[2] = b[2]-a[2];			\
	}
#define DOT( a, b ) ( a[0]*b[0] + a[1]*b[1] + a[2]*b[2] )
#define DIST( a, b ) sqrt( (a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])+(a[2]-b[2])*(a[2]-b[2]) )
#define LEN( a ) sqrt( (a[0])*(a[0])+(a[1])*(a[1])+(a[2])*(a[2]) )
#define NORM( a ) {					\
	double ilen;					\
	ilen = 1./LEN ( a );				\
	a[0] *= ilen; a[1] *= ilen; a[2] *= ilen;	\
	}
#define SIZEVEC( a, f ) {				\
	double ilen;					\
	ilen = f/LEN ( a );				\
	a[0] *= ilen; a[1] *= ilen; a[2] *= ilen;	\
	}
#define CURL( a, b, c ) {			\
	c[0] = a[1]*b[2] - a[2]*b[1];		\
	c[1] = a[2]*b[0] - a[0]*b[2];		\
	c[2] = a[0]*b[1] - a[1]*b[0];		\
	}

void dfProjection( struct inDumpFrame *in, struct dfFrameSetup *fs ) {
    double width,height;
    double vec[3];
    /*double norm[3]; -- not used: DCR 12/19/02*/

    in->r[0] = fs->target[0];
    in->r[1] = fs->target[1];
    in->r[2] = fs->target[2];
    in->nxPix = fs->nxPix;
    in->nyPix = fs->nyPix;
    in->bExpansion = fs->bExpansion;
    in->bPeriodic = fs->bPeriodic;
    in->iProject = fs->iProject;
    in->pScale1 = fs->pScale1;
    in->pScale2 = fs->pScale2;
    {
	int j;
	for (j=0;j<10;j++) {
	    in->ColStar[j] = fs->ColStar[j];
	    in->ColGas[j] = fs->ColGas[j];
	    in->ColDark[j] = fs->ColDark[j];
	    in->dColStar[j] = fs->dColStar[j];
	    in->dColGas[j] = fs->dColGas[j];
	    in->dColDark[j] = fs->dColDark[j];
	    }
	}
    in->nColStar = fs->nColStar;
    in->nColGas = fs->nColGas;
    in->nColDark = fs->nColDark;
	    
    in->bColMassWeight = fs->bColMassWeight;
    in->ColGasProperty = fs->ColGasProperty;
    in->ColGasPropMin = fs->ColGasPropMin;
    in->ColGasPropMax = fs->ColGasPropMax;
    in->ColDarkProperty = fs->ColDarkProperty;
    in->ColStarProperty = fs->ColStarProperty;
    in->bGasSph = fs->bGasSph;
    in->iColStarAge = fs->iColStarAge;
    in->bColLogInterp = fs->bColLogInterp;
    in->iLogScale = fs->iLogScale;
    in->iTarget = fs->iTarget;
    in->dGasSoftMul = fs->dGasSoftMul;
    in->dDarkSoftMul = fs->dDarkSoftMul;
    in->dStarSoftMul = fs->dStarSoftMul;
    in->iRender = fs->iRender;
	
    if (fs->bEyeAbsolute) {
	DIFF( fs->eye, in->r, in->z );
	}
    else {
	double zero[3]={ 0, 0, 0 };
	DIFF( fs->eye, zero, in->z );
	}

    /*	fprintf(stderr,"eye2: %i  %lf %lf %lf,  %lf %lf %lf",fs->bEye2,fs->eye2[0],fs->eye2[1],fs->eye2[2], in->z[0],in->z[1],in->z[2] ); */

    if (fs->bzEye1) SIZEVEC( in->z, fs->zEye1 );

    /* 	fprintf(stderr,"eye2: %i  %lf %lf %lf,  %lf %lf %lf",fs->bEye2,fs->eye2[0],fs->eye2[1],fs->eye2[2], in->z[0],in->z[1],in->z[2] ); */


    if (fs->bEye2) {
	SET( vec, fs->eye2 );
	if (fs->bzEye2) SIZEVEC( vec, fs->zEye2 );
	ADD ( in->z, vec );
	}

    /*	fprintf(stderr,"eye2: %i  %lf %lf %lf,  %lf %lf %lf",fs->bEye2,fs->eye2[0],fs->eye2[1],fs->eye2[2], in->z[0],in->z[1],in->z[2] ); */

    if (fs->bzEye) {
	in->zEye = fs->zEye;
	}
    else {
	in->zEye = LEN( in->z );
	}

    if (fs->bExpansion) in->zEye/= in->dExp; /* Use physical units for sizing? */
    /* zEye is used to scale viewport */

    assert( in->zEye > 0.0 );
    if (fs->bzClipFrac) {
	in->zClipNear = in->zEye*fs->zClipNear;
	in->zClipFar = in->zEye*fs->zClipFar;
	}
    else {
	in->zClipNear = fs->zClipNear;
	in->zClipFar = fs->zClipFar;
	}

    NORM( in->z );

    width = 2*tan( fs->FOV*M_PI/180.*0.5 )*in->zEye;
    height = width*in->nyPix/in->nxPix;

    CURL( fs->up, in->z, in->x );
    if (fs->iProject == DF_PROJECT_PERSPECTIVE) {
	SIZEVEC( in->x, (in->nxPix*0.5*in->zEye/width) );
	}
    else {
	SIZEVEC( in->x, (in->nxPix*0.5/width) );
	}

    CURL( in->z, in->x, in->y );
    if (fs->iProject == DF_PROJECT_PERSPECTIVE) {
	SIZEVEC( in->y, (in->nyPix*0.5*in->zEye/height) );
	}
    else {
	SIZEVEC( in->y, (in->nyPix*0.5/height) );
	}

    if (fs->bPeriodic) {
	in->nxRepNeg = 0;  /* Replicas for Periodic: Will setup sensibly ultimately */
	in->nxRepPos = 0;
	in->nyRepNeg = 0;  
	in->nyRepPos = 0;
	in->nzRepNeg = 0;  
	in->nzRepPos = 0;
	}
    else {
	in->nxRepNeg = 0;  /* Not Periodic */
	in->nxRepPos = 0;
	in->nyRepNeg = 0;  
	in->nyRepPos = 0;
	in->nzRepNeg = 0;  
	in->nzRepPos = 0;
	}
	
    /* Render helper variables */
    in->xlim = (in->nxPix-1)*.5;
    in->ylim = (in->nyPix-1)*.5;
    in->hmul = 4*sqrt(in->x[0]*in->x[0] + in->x[1]*in->x[1] + in->x[2]*in->x[2]);

    /* in->bNonLocal not set (an internal use only variable) */

    if (in->bVDetails) {
	printf("DF Projection: %i x %i FOV %f  width %f height %f\n",fs->nxPix,fs->nyPix,fs->FOV,width,height);
	printf("DF eye %f %f %f, target %f %f %f, (separation %f)\n",fs->eye[0],fs->eye[1],fs->eye[2],fs->target[0],fs->target[1],fs->target[2],in->zEye );
	printf("DF up %f %f %f  Z-Clipping: Near %f Far %f\n",fs->up[0],fs->up[1],fs->up[2],in->zClipNear,in->zClipFar);
	printf("DF Vectors: x %f %f %f, y %f %f %f z %f %f %f\n",in->x[0],in->x[1],in->x[2], in->y[0],in->y[1],in->y[2], in->z[0],in->z[1],in->z[2] );
	printf("DF Colours: star %f %f %f gas %f %f %f dark %f %f %f\n",in->ColStar[0].r,in->ColStar[0].g,in->ColStar[0].b,in->ColGas[0].r,in->ColGas[0].g,in->ColGas[0].b,in->ColDark[0].r,in->ColDark[0].g,in->ColDark[0].b);
	if (in->nColGas > 1 || in->nColStar > 1 || in->nColDark > 1) {
	    int iCol,nCol;
	    nCol = in->nColGas;
	    if (in->nColDark > nCol) nCol=in->nColDark;
	    if (in->nColStar > nCol) nCol=in->nColStar;
	    for (iCol=1;iCol<nCol;iCol++) {
		printf("DF Colours:");
		if (iCol < in->nColStar) printf(" star %f %f %f",in->ColStar[iCol].r,in->ColStar[iCol].g,in->ColStar[iCol].b);
		else printf(" star  -  -  -");
		if (iCol < in->nColGas) printf(" gas %f %f %f",in->ColGas[iCol].r,in->ColGas[iCol].g,in->ColGas[iCol].b);
		else printf(" gas  -  -  -");
		if (iCol < in->nColDark) printf(" dark %f %f %f\n",in->ColDark[iCol].r,in->ColDark[iCol].g,in->ColDark[iCol].b);
		else printf(" dark  -  -  -\n");
		}
	    }		
	}
    }

void dfParseOptions( struct DumpFrameContext *df, char * filename ) {
    FILE *fp;
    int n,nitem;
    char line[81],command[40],word[40];
    char FileBaseName[81]="Frame";
    char NumberingFormat[41]="";

    df->iDimension = DF_2D;
    df->iEncode = DF_ENCODE_PPM;
    df->iNumbering = DF_NUMBERING_FRAME;
    df->dMassGasMin  = 0;
    df->dMassGasMax  = DBL_MAX;
    df->dMassDarkMin = 0;
    df->dMassDarkMax = DBL_MAX;
    df->dMassStarMin = 0;
    df->dMassStarMax = DBL_MAX;
    df->nFrame = 0;
    df->bGetCentreOfMass = 0;
    df->bGetOldestStar = 0;
    df->bGetPhotogenic = 0;

    fp = fopen( filename, "r" );
    if (fp==NULL) return;

    for ( ;; ) {
	if (fgets( line, 81, fp ) == NULL) break;
	sscanf( line, "%s", command );
	if (!strcmp( command, "t" ) || !strcmp( command, "time" )) df->nFrameSetup++; /* count times */
	}

    rewind( fp );

    n=0;
    for ( ;; ) {
	if (fgets( line, 81, fp ) == NULL) break;
	nitem = sscanf( line, "%s", command );
	if (nitem != 1 || command[0]=='#') continue;
	else if (!strcmp( command, "dim") || !strcmp( command, "dimension") ) {
	    nitem = sscanf( line, "%s %s", command, word );
	    assert( nitem == 2 );
	    if (!strcmp( word, "2") ) {
		df->iDimension = DF_2D;
		}
	    else if (!strcmp( word, "3") ) {
		df->iDimension = DF_3D;
		}
	    else {
		fprintf(stderr,"DF Unknown dimensions: %s\n",word);
		assert( 0 );
		} 
	    }
	else if (!strcmp( command, "file") ) {
	    nitem = sscanf( line, "%s %80s", command, FileBaseName );
	    assert( strlen(FileBaseName) > 0 );
	    }
	else if (!strcmp( command, "gas") ) {
	    nitem = sscanf( line, "%s %s", command, word );
	    if (nitem == 2 && (!strcmp( word, "no") ||!strcmp( word, "off" ))) {
		df->dMassGasMin = DBL_MAX;
		}
	    else {
		nitem = sscanf( line, "%s %lf %lf", command,
				&df->dMassGasMin, &df->dMassGasMax );
		assert( nitem == 3 );
		}
	    }
	else if (!strcmp( command, "dark") ) {
	    nitem = sscanf( line, "%s %s", command, word );
	    if (nitem == 2 && (!strcmp( word, "no") ||!strcmp( word, "off" ))) {
		df->dMassDarkMin = DBL_MAX;
		}
	    else {
		nitem = sscanf( line, "%s %lf %lf", command,
				&df->dMassDarkMin, &df->dMassDarkMax );
		assert( nitem == 3 );
		}
	    }
	else if (!strcmp( command, "star") ) {
	    nitem = sscanf( line, "%s %s", command, word );
	    if (nitem == 2 && (!strcmp( word, "no") ||!strcmp( word, "off" ))) {
		df->dMassStarMin = DBL_MAX;
		}
	    else {
		nitem = sscanf( line, "%s %lf %lf", command,
				&df->dMassStarMin, &df->dMassStarMax );
		assert( nitem == 3 );
		}
	    }
	else if (!strcmp( command, "encode") ) {
	    nitem = sscanf( line, "%s %s", command, word );
	    assert( nitem == 2 );
	    if (!strcmp( word, "ppm") ) {
		df->iEncode = DF_ENCODE_PPM;
		}
	    else if (!strcmp( word, "png") ) {
		df->iEncode = DF_ENCODE_PNG;
#ifndef USE_PNG
		fprintf(stderr,"DF PNG encoding support not compiled in\n");
		assert(0);
#endif				
		}
	    else if (!strcmp( word, "rle") ) {
		df->iEncode = DF_ENCODE_RLE;
		fprintf(stderr,"DF RLE encoding not supported yet\n");
		assert(0);
		}
	    else if (!strcmp( word, "treezip") ) {
		df->iEncode = DF_ENCODE_TREEZIP;
		fprintf(stderr,"DF Voxel encoding not supported yet\n");
		assert(0);
		}
	    else {
		fprintf(stderr,"DF Unknown encoding: %s\n",word);
		assert(0);
		} 
	    }
	else if (!strcmp( command, "numbering") ) {
	    nitem = sscanf( line, "%s %s %40s", command, word, NumberingFormat );
	    if (nitem == 2 ) {
		nitem = sscanf( line, "%s %s", command, word );
		assert( nitem == 2 );
		}

	    if (!strcmp( word, "frame") ) {
		df->iNumbering = DF_NUMBERING_FRAME;
		}
	    else if (!strcmp( word, "step") ) {
		df->iNumbering = DF_NUMBERING_STEP;
		}
	    else if (!strcmp( word, "time") ) {
		df->iNumbering = DF_NUMBERING_TIME;
		}
	    else {
		fprintf(stderr,"DF Unknown numbering type: %s\n",word);
		assert( 0 );
		} 
	    }
	else if (!strcmp( command, "frame") ) {
	    nitem = sscanf( line, "%s %i", command, &df->nFrame );
	    assert( nitem == 2 );
	    }
	}
	
    if (df->iDimension == DF_3D) {
	fprintf(stderr,"DF Voxel projection/encoding not supported yet\n");
	assert(0);

	df->iEncode = DF_ENCODE_TREEZIP; /* Encode same thing */
	}
    else assert (df->iEncode != DF_ENCODE_TREEZIP);
	
    if (strlen(NumberingFormat)==0) {
	switch( df->iNumbering ) {
	case DF_NUMBERING_FRAME:
	    strcpy(NumberingFormat,"%09i");
	    break;
	case DF_NUMBERING_STEP:
	case DF_NUMBERING_TIME:
	    strcpy(NumberingFormat,"%012.6f");
	    break;
	    }
	}
	
    switch ( df->iEncode ) {
    case DF_ENCODE_PPM:
	sprintf(df->FileName,"%s.%s.ppm", FileBaseName, NumberingFormat );
	break;
    case DF_ENCODE_PNG:
	sprintf(df->FileName,"%s.%s.png", FileBaseName, NumberingFormat );
	break;
    default:
	sprintf(df->FileName,"%s.%s.null", FileBaseName, NumberingFormat );
	break;
	}

    fclose(fp);
    }

int dfPropertyToOuttype(char *word){
    if(!strcmp(word,"T") || !strcmp(word,"temp") || !strcmp(word,"Temp") ||
       !strcmp(word,"Temperature") || !strcmp(word,"temperature") )
	return OUT_TEMP_ARRAY;
    else if(!strcmp(word,"udot") || !strcmp(word,"uDot"))
	return OUT_UDOT_ARRAY;
    else if(!strcmp(word,"u") || !strcmp(word,"U"))
	return OUT_U_ARRAY;
    else if(!strcmp(word,"metals") || !strcmp(word,"Metals") ||
       !strcmp(word,"metallicity") || !strcmp(word,"Z"))
	return OUT_METALS_ARRAY;
    else return 0;
    }

void dfParseCameraDirections( struct DumpFrameContext *df, char * filename ) {
    FILE *fp;
    /*	struct inDumpFrame in; */
    struct dfFrameSetup fs;
    int n,nitem;
#define CMWOFF 1
#define CMWON  2
    int bCMW = 0;
    char line[81],command[40],word[40];

    df->bLoop = 0;
    df->dTimeMod = 0;
    df->dTimeLoop = 1e20;
    df->dPeriodLoop = 0;
	
    /* Defaults -- most of this is for 2D only */
    fs.dTime = 0;
    fs.nxPix = 800;
    fs.nyPix = 600;
    fs.target[0] = 0; 
    fs.target[1] = 0;
    fs.target[2] = 0;
    fs.eye[0] = 0;
    fs.eye[1] = 0;
    fs.eye[2] = -.5;
    fs.up[0] = 0;
    fs.up[1] = 1;
    fs.up[2] = 0;
    fs.zEye = 0.0;
    fs.zEye1 = 0.0;
    fs.zEye2 = 0.0;
    fs.eye2[0] = 0;
    fs.eye2[1] = 0;
    fs.eye2[2] = 0;
    fs.bEye2 = 0;
    fs.bzEye = 0;
    fs.bzEye1 = 0;
    fs.bzEye2 = 0;
    fs.bEyeAbsolute = 0;
    fs.FOV = 90.;
    fs.zClipNear = 0.01;
    fs.zClipFar = 2.0;
    fs.bzClipFrac = 1;
    fs.bExpansion = 0;  /* to physical? */
    fs.bPeriodic = 0;  /* Periodic? */
    fs.iProject = DF_PROJECT_PERSPECTIVE;
    /* Render */
    {
	int j;
	for (j=0;j<10;j++) {
	    fs.ColDark[j].r = 1.0;
	    fs.ColDark[j].g = 0.7;
	    fs.ColDark[j].b = 0.5;
	    fs.ColGas[j].r = 0.5;
	    fs.ColGas[j].g = 1.0;
	    fs.ColGas[j].b = 0.7;
	    fs.ColStar[j].r = 0.5;
	    fs.ColStar[j].g = 0.7;
	    fs.ColStar[j].b = 1.0;
	    }
	}
    fs.nColDark = 0;
    fs.nColGas = 0;
    fs.nColStar = 0;
    fs.bColMassWeight = 0;
    fs.ColGasProperty = OUT_TEMP_ARRAY;
    fs.ColGasPropMin = 2;
    fs.ColGasPropMax = 7;
    fs.ColDarkProperty = 0;
    fs.ColStarProperty = OUT_TIMEFORM_ARRAY;
    fs.bGasSph = 1;
    fs.iColStarAge = DF_STAR_AGE_BASIC;
    fs.bColLogInterp = 1;
    fs.iLogScale = DF_LOG_NULL;
    fs.iTarget = DF_TARGET_USER;
    fs.dGasSoftMul = 1;
    fs.dDarkSoftMul = 1;
    fs.dStarSoftMul = 1;
	

    fs.iRender = DF_RENDER_POINT;

    fp = fopen( filename, "r" );
    if (fp==NULL) {
	fprintf(stderr,"DF Could not open camera director file: %s\n",filename );

	df->iFrameSetup = 0;
	df->nFrameSetup = 1;
	df->fs = (struct dfFrameSetup *) malloc(sizeof(struct dfFrameSetup)*df->nFrameSetup);
	assert( df->fs != NULL );

	if (df->bVDetails) printf("DF Default Frame Setup\n" );
	/* dfProjection( &in, &fs ); */ /* Redundant now -- done again later */
	df->fs[0] = fs;
	return;
	}

    if (df->bVDetails) printf("DF Reading camera directions from: %s\n",filename );
	
    df->iFrameSetup = -1;
    df->nFrameSetup=1; /* Defaults */
    for ( ;; ) {
	if (fgets( line, 81, fp ) == NULL) break;
	sscanf( line, "%s", command );
	if (!strcmp( command, "t" ) || !strcmp( command, "time" )) df->nFrameSetup++; /* count times */
	}

    rewind( fp );
    df->fs = (struct dfFrameSetup *) malloc(sizeof(struct dfFrameSetup)*df->nFrameSetup);
    assert( df->fs != NULL );

    n=0;
    if (df->bVDetails) printf("DF Default Frame Setup\n" );
    for ( ;; ) {
	if (fgets( line, 81, fp ) == NULL) break;
	nitem = sscanf( line, "%s", command );
	if (nitem != 1 || command[0]=='#') continue;
	else if (!strcmp( command, "t" ) || !strcmp( command, "time" )) {
	    /* dfProjection( &in, &fs );  */ /* Redundant now -- done again later */
	    df->fs[n] = fs;
	    n++; /* next setup */

	    nitem = sscanf( line, "%s %lf", command, &fs.dTime );
	    assert( nitem == 2 );
			
	    if (df->bVDetails) printf("DF Frame Setup from time: %f\n",fs.dTime );
	    }
	else if (!strcmp( command, "target") ) {
	    nitem = sscanf( line, "%s %lf %lf %lf", command, &fs.target[0], &fs.target[1], &fs.target[2] );
	    if (nitem != 4) {
		nitem = sscanf( line, "%s %s", command, word );
		assert( nitem == 2);
		if (!strcmp( word, "gas") ) {
		    fs.iTarget = DF_TARGET_COM_GAS;
		    df->bGetCentreOfMass = 1;
		    }
		else if (!strcmp( word, "dark") ) {
		    fs.iTarget = DF_TARGET_COM_DARK;
		    df->bGetCentreOfMass = 1;
		    }
		else if (!strcmp( word, "star") ) {
		    fs.iTarget = DF_TARGET_COM_STAR;
		    df->bGetCentreOfMass = 1;
		    }
		else if (!strcmp( word, "all") ) {
		    fs.iTarget = DF_TARGET_COM_ALL;
		    df->bGetCentreOfMass = 1;
		    }
		else if (!strcmp( word, "oldeststar" ) ) {
		    fs.iTarget = DF_TARGET_OLDEST_STAR;
		    df->bGetOldestStar = 1;
		    }
		else if (!strcmp( word, "photogenic" ) ) {
		    fs.iTarget = DF_TARGET_PHOTOGENIC;
		    df->bGetPhotogenic = 1;
		    }
		else {
		    assert(0);
		    }
		}
	    }
	else if (!strcmp( command, "eye") || !strcmp( command, "eye1") ) {
	    nitem = sscanf( line, "%s %lf %lf %lf", command, &fs.eye[0], &fs.eye[1], &fs.eye[2] );
	    if ( nitem != 4 ) {
		nitem = sscanf( line, "%s %s", command, word );
		assert( nitem == 2 );
		if (!strcmp( word, "rel") ) {
		    fs.bEyeAbsolute = 0;
		    }
		else if (!strcmp( word, "abs") ) {
		    fs.bEyeAbsolute = 1;
		    }
		else {
		    fprintf(stderr,"DF Unknown eye offset type: %s\n",word);
		    assert( 0 );
		    }
		}
	    }
	else if (!strcmp( command, "eye2") ) {
	    nitem = sscanf( line, "%s %lf %lf %lf", command, &fs.eye2[0], &fs.eye2[1], &fs.eye2[2] );
	    fs.bEye2 = 1;
	    assert( nitem == 4 );
	    }
	else if (!strcmp( command, "up") ) {
	    nitem = sscanf( line, "%s %lf %lf %lf", command, &fs.up[0], &fs.up[1], &fs.up[2] );
	    assert( nitem == 4 );
	    }
	else if (!strcmp( command, "size") ) {
	    nitem = sscanf( line, "%s %i %i", command, &fs.nxPix, &fs.nyPix );
	    assert( nitem == 3 );
	    assert( fs.nxPix*fs.nyPix <= DF_NXPIXMAX*DF_NYPIXMAX );
	    }
	else if (!strcmp( command, "zeye") ) {
	    fs.bzEye = 1;
	    nitem = sscanf( line, "%s %lf", command, &fs.zEye );
	    assert( nitem == 2 );
	    }
	else if (!strcmp( command, "zeye1") ) {
	    fs.bzEye1 = 1;
	    nitem = sscanf( line, "%s %lf", command, &fs.zEye1 );
	    assert( nitem == 2 );
	    }
	else if (!strcmp( command, "zeye2") ) {
	    fs.bzEye2 = 1;
	    nitem = sscanf( line, "%s %lf", command, &fs.zEye2 );
	    assert( nitem == 2 );
	    }
	else if (!strcmp( command, "exp") || !strcmp( command, "physical") ) {
	    fs.bExpansion = 1;
	    }
	else if (!strcmp( command, "fov") || !strcmp( command, "FOV")  ) {
	    nitem = sscanf( line, "%s %lf", command, &fs.FOV );
	    assert( nitem == 2 );
	    }
	else if (!strcmp( command, "loop") ) {
	    nitem = sscanf( line, "%s %lf %lf", command, &df->dTimeLoop, &df->dPeriodLoop );
	    assert( nitem == 3 );
	    }
	else if (!strcmp( command, "clip") ) {
	    df->bLoop = 1;
	    nitem = sscanf( line, "%s %lf %lf", command, &fs.zClipNear, &fs.zClipFar );
	    assert( nitem == 3 );
	    }
	else if (!strcmp( command, "clipabs") ) {
	    fs.bzClipFrac = 0;
	    nitem = sscanf( line, "%s %lf %lf", command, &fs.zClipNear, &fs.zClipFar );
	    assert( nitem == 3 );
	    }
	else if (!strcmp( command, "clipfrac") ) {
	    fs.bzClipFrac = 1;
	    nitem = sscanf( line, "%s %lf %lf", command, &fs.zClipNear, &fs.zClipFar );
	    assert( nitem == 3 );
	    }
	else if (!strcmp( command, "project") ) {
	    nitem = sscanf( line, "%s %s", command, word );
	    assert( nitem == 2 );
	    if (!strcmp( word, "ortho") ) {
		fs.iProject = DF_PROJECT_ORTHO;
		}
	    else if (!strcmp( word, "perspective") ) {
		fs.iProject = DF_PROJECT_PERSPECTIVE;
		}
	    else {
		fprintf(stderr,"DF Unknown projection: %s\n",word);
		assert( 0 );
		} 
	    }
	else if (!strcmp( command, "colgasbyproperty" )) {
	    nitem = sscanf( line, "%s %s %f %f", command, word,
			    &fs.ColGasPropMin, &fs.ColGasPropMax);
	    if ((nitem==2)||(nitem==4)) 
		fs.ColGasProperty=dfPropertyToOuttype(word);
	    else {
		fprintf(stderr,"DF wrong argument number to colgasbyproperty property [min max]: %s\n",word);
		assert( 0 );
		} 
	    }
	else if (!strcmp( command, "gascolormap" )) {
	    float scaler=1.0;
	    nitem = sscanf( line, "%s %s %f", command, word, &scaler);
	    if (nitem==3 && !strcmp(word,"jet")) {
		fs.ColGas[fs.nColGas].r*=scaler*1.0;
		fs.ColGas[fs.nColGas].g*=0.0;
		fs.ColGas[fs.nColGas].b*=scaler*0.5;
		fs.dColGas[fs.nColGas]=0.0;
		fs.nColGas++;
		fs.ColGas[fs.nColGas].r*=scaler;
		fs.ColGas[fs.nColGas].g=0;
		fs.ColGas[fs.nColGas].b*=scaler;
		fs.dColGas[fs.nColGas]=0.15;
		fs.nColGas++;
		fs.ColGas[fs.nColGas].r*=0;
		fs.ColGas[fs.nColGas].g*=scaler;
		fs.ColGas[fs.nColGas].b*=scaler;
		fs.dColGas[fs.nColGas]=0.4;
		fs.nColGas++;
		fs.ColGas[fs.nColGas].r*=scaler;
		fs.ColGas[fs.nColGas].g*=scaler;
		fs.ColGas[fs.nColGas].b*=0;
		fs.dColGas[fs.nColGas]=0.7;
		fs.nColGas++;
		fs.ColGas[fs.nColGas].r*=scaler;
		fs.ColGas[fs.nColGas].g*=0;
		fs.ColGas[fs.nColGas].b*=0;
		fs.dColGas[fs.nColGas]=0.9;
		fs.nColGas++;
		fs.ColGas[fs.nColGas].r*=0.5*scaler;
		fs.ColGas[fs.nColGas].g*=0;
		fs.ColGas[fs.nColGas].b*=0;
		fs.dColGas[fs.nColGas]=1.0;
		fs.nColGas++;
		}
	    }
	else if (!strcmp( command, "coldarkbyproperty" )) {
	    nitem = sscanf( line, "%s %s", command, word);
	    if (nitem==2) fs.ColDarkProperty=dfPropertyToOuttype(word);
	    }
	else if (!strcmp( command, "colstarbyproperty" )) {
	    nitem = sscanf( line, "%s %s", command, word);
	    if (nitem==2) fs.ColStarProperty=dfPropertyToOuttype(word);
	    }
	else if (!strcmp( command, "colgasint" )) {
	    float scaler;
	    assert(fs.nColGas < 10);
	    nitem = sscanf( line, "%s %f %f %f %f %f", command, &fs.ColGas[fs.nColGas].r, &fs.ColGas[fs.nColGas].g, &fs.ColGas[fs.nColGas].b, &fs.dColGas[fs.nColGas], &scaler );
	    if (nitem == 6) {
		assert( !(bCMW & CMWOFF) );
		bCMW |= CMWON;
		fs.bColMassWeight = 1;
		fs.ColGas[fs.nColGas].r*=scaler;
		fs.ColGas[fs.nColGas].g*=scaler;
		fs.ColGas[fs.nColGas].b*=scaler;
		}
	    else {
		nitem = sscanf( line, "%s %f %f %f %f", command, &fs.ColGas[fs.nColGas].r, &fs.ColGas[fs.nColGas].g, &fs.ColGas[fs.nColGas].b, &fs.dColGas[fs.nColGas] );
		assert( nitem == 5 );
		assert( !(bCMW & CMWON) );
		bCMW |= CMWOFF;
		}
	    if (fs.nColGas > 0) assert(fs.dColGas[fs.nColGas] > fs.dColGas[fs.nColGas-1]);
	    fs.nColGas++;
	    }
	else if (!strcmp( command, "colgas" )) {
	    float scaler;
	    nitem = sscanf( line, "%s %f %f %f %f", command, &fs.ColGas[0].r, &fs.ColGas[0].g, &fs.ColGas[0].b, &scaler );
	    if (nitem == 5) {
		assert( !(bCMW & CMWOFF) );
		bCMW |= CMWON;
		fs.bColMassWeight = 1;
		fs.ColGas[0].r*=scaler;
		fs.ColGas[0].g*=scaler;
		fs.ColGas[0].b*=scaler;
		}
	    else {
		nitem = sscanf( line, "%s %f %f %f", command, &fs.ColGas[0].r, &fs.ColGas[0].g, &fs.ColGas[0].b );
		assert( nitem == 4 );
		assert( !(bCMW & CMWON) );
		bCMW |= CMWOFF;
		}
	    fs.nColGas = 0;
	    }
	else if (!strcmp( command, "coldarkint" )) {
	    float scaler;
	    assert(fs.nColDark < 10);
	    nitem = sscanf( line, "%s %f %f %f %f %f", command, &fs.ColDark[fs.nColDark].r, &fs.ColDark[fs.nColDark].g, &fs.ColDark[fs.nColDark].b, &fs.dColDark[fs.nColDark], &scaler );
	    if (nitem == 6) {
		assert( !(bCMW & CMWOFF) );
		bCMW |= CMWON;
		fs.bColMassWeight = 1;
		fs.ColDark[fs.nColDark].r*=scaler;
		fs.ColDark[fs.nColDark].g*=scaler;
		fs.ColDark[fs.nColDark].b*=scaler;
		}
	    else {
		nitem = sscanf( line, "%s %f %f %f %f", command, &fs.ColDark[fs.nColDark].r, &fs.ColDark[fs.nColDark].g, &fs.ColDark[fs.nColDark].b, &fs.dColDark[fs.nColDark] );
		assert( nitem == 5 );
		assert( !(bCMW & CMWON) );
		bCMW |= CMWOFF;
		}
	    if (fs.nColDark > 0) assert(fs.dColDark[fs.nColDark] > fs.dColDark[fs.nColDark-1]);
	    fs.nColDark++;
	    }
	else if (!strcmp( command, "coldark" )) {
	    float scaler;
	    nitem = sscanf( line, "%s %f %f %f %f", command, &fs.ColDark[0].r, &fs.ColDark[0].g, &fs.ColDark[0].b, &scaler );
	    if (nitem == 5) {
		assert( !(bCMW & CMWOFF) );
		bCMW |= CMWON;
		fs.bColMassWeight = 1;
		fs.ColDark[0].r*=scaler;
		fs.ColDark[0].g*=scaler;
		fs.ColDark[0].b*=scaler;
		}
	    else {
		nitem = sscanf( line, "%s %f %f %f", command, &fs.ColDark[0].r, &fs.ColDark[0].g, &fs.ColDark[0].b );
		assert( nitem == 4 );
		assert( !(bCMW & CMWON) );
		bCMW |= CMWOFF;
		}
	    fs.nColDark = 0;
	    }
	else if (!strcmp( command, "colstarint" )) {
	    float scaler;
	    assert(fs.nColStar < 10);
	    nitem = sscanf( line, "%s %f %f %f %f %f", command, &fs.ColStar[fs.nColStar].r, &fs.ColStar[fs.nColStar].g, &fs.ColStar[fs.nColStar].b, &fs.dColStar[fs.nColStar], &scaler );
	    if (nitem == 6) {
		assert( !(bCMW & CMWOFF) );
		bCMW |= CMWON;
		fs.bColMassWeight = 1;
		fs.ColStar[fs.nColStar].r*=scaler;
		fs.ColStar[fs.nColStar].g*=scaler;
		fs.ColStar[fs.nColStar].b*=scaler;
		}
	    else {
		nitem = sscanf( line, "%s %f %f %f %f", command, &fs.ColStar[fs.nColStar].r, &fs.ColStar[fs.nColStar].g, &fs.ColStar[fs.nColStar].b, &fs.dColStar[fs.nColStar] );
		assert( nitem == 5 );
		assert( !(bCMW & CMWON) );
		bCMW |= CMWOFF;
		}
	    if (fs.nColStar > 0) assert(fs.dColStar[fs.nColStar] > fs.dColStar[fs.nColStar-1]);
	    fs.nColStar++;
	    }
	else if (!strcmp( command, "colstar" ) 
		 || !strcmp( command, "colstarbri" ) 
		 || !strcmp( command, "colstarcol" )
		 || !strcmp( command, "colstarbricol" )) {

	    float scaler;
	    if (!strcmp( command, "colstar" ))
		fs.iColStarAge = DF_STAR_AGE_BASIC;
	    else if(!strcmp( command, "colstarbri" ))
		fs.iColStarAge = DF_STAR_AGE_BRIGHT;
	    else if(!strcmp( command, "colstarcol" ))
		fs.iColStarAge = DF_STAR_AGE_COLOUR;
	    else if(!strcmp( command, "colstarbricol" ))
		fs.iColStarAge = DF_STAR_AGE_BRIGHT_COLOUR;

	    fs.nColStar = 0;
	    nitem = sscanf( line, "%s %f %f %f %f", command, &fs.ColStar[0].r, &fs.ColStar[0].g, &fs.ColStar[0].b, &scaler );
	    if (nitem == 5) {
		assert( !(bCMW & CMWOFF) );
		bCMW |= CMWON;
		fs.bColMassWeight = 1;
		fs.ColStar[0].r*=scaler;
		fs.ColStar[0].g*=scaler;
		fs.ColStar[0].b*=scaler;
		}
	    else {
		nitem = sscanf( line, "%s %f %f %f", command, &fs.ColStar[0].r, &fs.ColStar[0].g, &fs.ColStar[0].b );
		assert( nitem == 4 );
		assert( !(bCMW & CMWON) );
		bCMW |= CMWOFF;
		}
	    }
	else if (!strcmp( command, "colmass" )) {
	    assert( !(bCMW & CMWOFF) );
	    bCMW |= CMWON;
	    fs.bColMassWeight = 1;
	    }
	else if (!strcmp( command, "colloginterp" )) {
	    nitem = sscanf( line, "%s %s", command, word );
	    if (nitem == 2 && (!strcmp( word, "no") ||!strcmp( word, "off" ))) {
		fs.bColLogInterp = 0;
		}
	    fs.bColLogInterp = 1;
	    }
	else if (!strcmp( command, "logscale" )) {
	    nitem = sscanf( line, "%s %lf %lf", command, &fs.pScale1, &fs.pScale2 );
	    assert( nitem == 3 );
	    fs.iLogScale = DF_LOG_SATURATE;
	    }
	else if (!strcmp( command, "logscalecoloursafe" )) {
	    nitem = sscanf( line, "%s %lf %lf", command, &fs.pScale1, &fs.pScale2 );
	    assert( nitem == 3 );
	    fs.iLogScale = DF_LOG_COLOURSAFE;
	    }
	else if (!strcmp( command, "softgassph" )) {
	    fs.bGasSph = 1;
	    }
	else if (!strcmp( command, "softgas" )) {
	    fs.bGasSph = 0;
	    nitem = sscanf( line, "%s %lf", command, &fs.dGasSoftMul );
	    assert( nitem == 2 );
	    }
	else if (!strcmp( command, "softdark" )) {
	    nitem = sscanf( line, "%s %lf", command, &fs.dDarkSoftMul );
	    assert( nitem == 2 );
	    }
	else if (!strcmp( command, "softstar" )) {
	    nitem = sscanf( line, "%s %lf", command, &fs.dStarSoftMul );
	    assert( nitem == 2 );
	    }
	else if (!strcmp( command, "render") ) {
	    nitem = sscanf( line, "%s %s", command, word );
	    assert( nitem == 2 );
	    if (!strcmp( word, "point") ) {
		fs.iRender = DF_RENDER_POINT;
		}
	    else if (!strcmp( word, "tsc") ) {
		fs.iRender = DF_RENDER_TSC;
		}
	    else if (!strcmp( word, "solid") ) {
		fs.iRender = DF_RENDER_SOLID;
		}
	    else if (!strcmp( word, "shine") ) {
		fs.iRender = DF_RENDER_SHINE;
		}
	    else {
		fprintf(stderr,"DF Unknown rendering: %s\n",word);
		assert( 0 );
		} 
	    }
	}
    if (fs.nColGas == 1) fs.nColGas = 0;
    if (fs.nColDark == 1) fs.nColDark = 0;
    if (fs.nColStar == 1) fs.nColStar = 0;
	
    /* Redundant now -- done again later */
    /*	
	in.bVDetails = df->bVDetails;
	dfProjection( &in, &fs );
    */ 

    df->fs[n] = fs;
    n++; /* next setup */
    if (n==1) df->iFrameSetup = 0;
    assert ( n == df->nFrameSetup );

    fclose(fp);
    }

void dfSetupFrame( struct DumpFrameContext *df, double dTime, double dStep, 
		   double dExp, double *com, double dMinGasMass,
		   struct inDumpFrame *vin ) {
    struct dfFrameSetup fs;

    int ifs = df->iFrameSetup;

    vin->dTime = dTime;
    vin->dYearUnit = df->dYearUnit;
    vin->dStep = dStep;
    vin->dExp = dExp;
    vin->bVDetails = df->bVDetails;
    vin->dMassStarMin = df->dMassStarMin;
    vin->dMassStarMax = df->dMassStarMax;
    vin->dMassGasMin  = df->dMassGasMin;
    vin->dMassGasMax  = df->dMassGasMax;
    vin->dMassDarkMin = df->dMassDarkMin;
    vin->dMassDarkMax = df->dMassDarkMax;
    /* Note slight name difference from above properties. This is from
       min gas mass star formation and is being used to normalize the
       brightness scale better.  The above masses would be set
       based on configuration file.  Why are they set every time a
       frame dumps? -GSS March 2013 */
    vin->dMinGasMass = dMinGasMass;

    if (df->bLoop) {
	dTime -= df->dTimeMod;
	if (dTime > df->dPeriodLoop+df->dTimeLoop) {
	    dTime -= df->dPeriodLoop;
	    df->dTimeMod += df->dPeriodLoop;
	    ifs = -1;
	    }
	}
	
    if (ifs == -1) {
	assert( df->nFrameSetup > 1 );
	if (dTime < df->fs[1].dTime) { /* Outside range */
	    fprintf(stderr,"DF WARNING Initial Time outside camera direction table: %g < %g\n",dTime,df->fs[1].dTime);
	    df->iFrameSetup = ifs = 0;
	    }
	else {
	    ifs = 1;
	    while (ifs < df->nFrameSetup && dTime >= df->fs[ifs+1].dTime ) ifs++;
	    if (ifs >= df->nFrameSetup-1) { /* Outside Range */
		if (dTime == df->fs[df->nFrameSetup-1].dTime && ifs > 1) ifs--;
		else {
		    fprintf(stderr,"DF WARNING Initial Time outside camera direction table: %g > %g\n",dTime,df->fs[df->nFrameSetup-1].dTime);
		    df->iFrameSetup = ifs = 0;
		    }
		}
	    df->iFrameSetup = ifs;
	    dfGetCoeff( df, ifs );
	    }
	}
    else if (df->nFrameSetup > 1) {
	while (dTime > df->fs[ifs+1].dTime || (!ifs && dTime >= df->fs[ifs+1].dTime)) { 
	    ifs++;
	    if (ifs >= df->nFrameSetup-1) {
		fprintf(stderr,"DF WARNING Time outside camera direction table: %g > %g\n",dTime,df->fs[df->nFrameSetup-1].dTime);
		df->iFrameSetup = ifs = 0;
		break;
		}
	    df->iFrameSetup = ifs;
	    dfGetCoeff( df, ifs );
	    }
	}

	
    if (df->bVDetails) printf("DF Interpolating at t=%g Setups: %i (t=%g) %i (t=%g)\n",dTime,ifs,df->fs[ifs].dTime,ifs+1,df->fs[ifs+1].dTime);

    fs = df->fs[ifs];

    /* Nothing to interpolate? */
    if (ifs) {
	/* 
	   Interpolate Eye position, FOV etc... 
	   from df->fs[ifs].dTime <= dTime <= df->fs[ifs+1].dTime
	*/
	dfInterp( df, &fs, (dTime-fs.dTime)*df->rdt );
	}

    switch (fs.iTarget) {
    case DF_TARGET_COM_GAS:
	if (com[3]>0) {
	    fs.target[0] = com[0]/com[3];
	    fs.target[1] = com[1]/com[3];
	    fs.target[2] = com[2]/com[3];
	    }
	break;
    case DF_TARGET_COM_DARK:
	if (com[7]>0) {
	    fs.target[0] = com[4]/com[7];
	    fs.target[1] = com[5]/com[7];
	    fs.target[2] = com[6]/com[7];
	    }
	break;
    case DF_TARGET_COM_STAR:
	if (com[11]>0) {
	    fs.target[0] = com[8]/com[11];
	    fs.target[1] = com[9]/com[11];
	    fs.target[2] = com[10]/com[11];
	    }
	break;
    case DF_TARGET_COM_ALL:
	if (com[3]+com[7]+com[11]>0) {
	    fs.target[0] = (com[0]+com[4]+com[8])/(com[3]+com[7]+com[11]);
	    fs.target[1] = (com[1]+com[5]+com[9])/(com[3]+com[7]+com[11]);
	    fs.target[2] = (com[2]+com[6]+com[10])/(com[3]+com[7]+com[11]);
	    }
	break;
    case DF_TARGET_OLDEST_STAR:
	if (com[3] < FLT_MAX) {
	    fs.target[0] = com[0];
	    fs.target[1] = com[1];
	    fs.target[2] = com[2];
	    }
	break;
    case DF_TARGET_PHOTOGENIC:
	if (com[3] < FLT_MAX) {
	    fs.target[0] = com[0]/com[3];
	    fs.target[1] = com[1]/com[3];
	    fs.target[2] = com[2]/com[3];
	    }
	break;
	}

    dfProjection( vin, &fs ); 
    }


void dfClearImage( struct inDumpFrame *in, void *vImage, int *nImage ) {
    DFIMAGE *Image = vImage;
    int i;
    DFIMAGE blank;

    blank.r = 0;
    blank.g = 0;
    blank.b = 0;

    *nImage = in->nxPix * in->nyPix * sizeof(DFIMAGE);

    for (i=in->nxPix*in->nyPix-1;i>=0;i--) Image[i] = blank;
    }

DFIMAGE dfColInterp( int nCol, float *dCol, DFIMAGE *Col, FLOAT fVar, int bColLogInterp ) 
{
    int iCol;
    float xCol;
    DFIMAGE col;
    
    if (!nCol) return Col[0];
    if (bColLogInterp) fVar = log10(fVar);
    if (fVar <= dCol[0]) {
	return Col[0];
	}
    else if (fVar >= dCol[nCol-1]) {
	return Col[nCol-1];
	}
    iCol = 1; /* iCol specifies which range in colormap we're in */
    for (;;) { /*find where we are in the colormap */
	if (fVar <= dCol[iCol]) break;
	iCol++;
        if (iCol >= nCol) return Col[nCol-1];
	}
    /* xCol is the fraction distance through colormap range */
    xCol = (fVar - dCol[iCol-1])/(dCol[iCol]-dCol[iCol-1]);
    col.r = xCol*Col[iCol].r + (1-xCol)*Col[iCol-1].r;
    col.g = xCol*Col[iCol].g + (1-xCol)*Col[iCol-1].g;
    col.b = xCol*Col[iCol].b + (1-xCol)*Col[iCol-1].b;
    return col;
    }

DFIMAGE dfInterpOnColorMap(int nCol,float *dCol,DFIMAGE *Col,FLOAT fVar)
{
    /* fVar should come in between 0 and 1 or else the values will pin
       to the min/max colors -GSS March 2013 */
    int iCol;
    float xCol;
    DFIMAGE col;
    
    if (!nCol) return Col[0];
    if (fVar <= dCol[0]) {
	return Col[0];
	}
    else if (fVar >= dCol[nCol-1]) {
	return Col[nCol-1];
	}
    iCol = 1; /* iCol specifies which range in colormap we're in */
    for (;;) { /*find where we are in the colormap */
	if (fVar <= dCol[iCol]) break;
	iCol++;
        if (iCol >= nCol) return Col[nCol-1];
	}
    /* xCol is the fraction distance through colormap range */
    xCol = (fVar - dCol[iCol-1])/(dCol[iCol]-dCol[iCol-1]);
    col.r = xCol*Col[iCol].r + (1-xCol)*Col[iCol-1].r;
    col.g = xCol*Col[iCol].g + (1-xCol)*Col[iCol-1].g;
    col.b = xCol*Col[iCol].b + (1-xCol)*Col[iCol-1].b;
    return col;
    }

void dfRenderParticlePoint( struct inDumpFrame *in, void *vImage, PKD pkd, 
			    PARTICLE *p)
{
    DFIMAGE *Image = vImage;
    DFIMAGE col; /* Colour */
    double x,y,z,dr[3];
    int j;
    int xp,yp;
    int iActive=p->iActive;
    double *r=p->r;
    double fMass=p->fMass;

    if (pkdIsGas(pkd,p)) {
	if (fMass < in->dMassGasMin || fMass > in->dMassGasMax) return;
	col = in->ColGas[0];
	}
    else if (pkdIsDark(pkd,p)) {
	if (fMass < in->dMassDarkMin || fMass > in->dMassDarkMax) return;
	col = in->ColDark[0];
	}
    else if (pkdIsStar(pkd,p)) {
	if (fMass < in->dMassStarMin || fMass > in->dMassStarMax) return;
	col = in->ColStar[0];
	}

    for (j=0;j<3;j++) {
	dr[j] = r[j]-in->r[j];
	}
    z = dr[0]*in->z[0] + dr[1]*in->z[1] + dr[2]*in->z[2] + in->zEye;
    if (z >= in->zClipNear && z <= in->zClipFar) {
	x = dr[0]*in->x[0] + dr[1]*in->x[1] + dr[2]*in->x[2];
	if (in->iProject == DF_PROJECT_PERSPECTIVE) x/=z;
	if (fabs(x)<in->xlim) {
	    y = dr[0]*in->y[0] + dr[1]*in->y[1] + dr[2]*in->y[2];
	    if (in->iProject == DF_PROJECT_PERSPECTIVE) y/=z;
	    if (fabs(y)<in->ylim) {
		xp = x+in->xlim;
		yp = in->ylim-y; /* standard screen convention */
		Image[ xp + yp*in->nxPix ].r += col.r;
		Image[ xp + yp*in->nxPix ].g += col.g;
		Image[ xp + yp*in->nxPix ].b += col.b;
		}
	    }
	}
    }


void dfRenderParticleTSC( struct inDumpFrame *in, void *vImage, PKD pkd, 
			  PARTICLE *p)
{
    DFIMAGE *Image = vImage;
    DFIMAGE col; /* Colour */
    double h;
    double x,y,z,dr[3],br0=1;
    int hint;
    int j;
    int xp,yp;
    int iActive=p->iActive;
    double *r=p->r;
    double fMass=p->fMass;
    double fSoft=p->fSoft;
    double fAge=(in->dTime-p->fTimeForm)*in->dYearUnit;
    FLOAT fVar;

#ifdef DEBUGTRACK
    if (fAge > 3) {
	ntrack++;
	trackp++;
	if ((trackp<DEBUGTRACK)) {
	    printf("p%d %f %f %f  %f %f\n",trackp,r[0],r[1],r[2],fMass,fSoft);
	    }
	}
#endif

    if (pkdIsGas(pkd,p)) {
	if (fMass < in->dMassGasMin || fMass > in->dMassGasMax) return;
	if (in->bGasSph) h = sqrt(p->fBall2)*0.5*in->dGasSoftMul;
	else h = fSoft*in->dGasSoftMul;
	fVar = VecType(pkd,p,1,in->ColGasProperty);
	fVar = in->bColLogInterp ? log10(fVar):fVar;
	fVar = (fVar-in->ColGasPropMin)/(in->ColGasPropMax-in->ColGasPropMin);
	/* For temperature, try scaling brightness by the temperature
	   since there is so much less hot gas. */
	if (in->ColGasProperty == OUT_TEMP_ARRAY) {
	    br0 *= exp(3.0*fVar);
	    }
	col = dfInterpOnColorMap(in->nColGas,&in->dColGas[0],&in->ColGas[0],
				 fVar);
	}
    else if (pkdIsDark(pkd,p)) {
	if (fMass < in->dMassDarkMin || fMass > in->dMassDarkMax) return;
	h = fSoft*in->dDarkSoftMul;
	col = dfColInterp( in->nColDark, &in->dColDark[0], &in->ColDark[0], 
		       VecType(pkd,p,1,in->ColDarkProperty),in->bColLogInterp);
	}
    else if (pkdIsStar(pkd,p)) {
	float al;
	if (fMass < in->dMassStarMin || fMass > in->dMassStarMax) return;
	h = fSoft*in->dStarSoftMul;
	switch (in->iColStarAge) {
	case DF_STAR_AGE_BRIGHT_COLOUR:
	case DF_STAR_AGE_COLOUR:
	    al = (log(fabs(fAge+1e6))-13.815)*(1/9.6);
	    col.b = in->ColStar[0].b*fabs(1-0.7*al);
	    col.g = in->ColStar[0].g*0.4;
	    col.r = in->ColStar[0].r*(0.4+0.32*al);
	    break;
	default:
	    col = dfColInterp( in->nColStar, &in->dColStar[0], &in->ColStar[0],
			       VecType(pkd,p,1,in->ColDarkProperty), 
			       in->bColLogInterp);
	    }
	switch (in->iColStarAge) {
	case DF_STAR_AGE_BASIC:
	    break;
	case DF_STAR_AGE_BRIGHT:
	case DF_STAR_AGE_BRIGHT_COLOUR:
	    br0 = 1./(fabs(fAge)/1e6+1.);
	    break;
	case DF_STAR_AGE_COLOUR:
	    break;
	    }
	}

    if (in->bColMassWeight) {
	if (in->dMinGasMass>0) br0*=fMass/in->dMinGasMass;
	else br0*=fMass;
	}
	
    for (j=0;j<3;j++) {
	dr[j] = r[j]-in->r[j];
	}

    z = dr[0]*in->z[0] + dr[1]*in->z[1] + dr[2]*in->z[2] + in->zEye;
    if (z >= in->zClipNear && z <= in->zClipFar) {

	if (in->iProject == DF_PROJECT_PERSPECTIVE) h = h*in->hmul/z;
	else h = h*in->hmul;
	hint = h;
	x = dr[0]*in->x[0] + dr[1]*in->x[1] + dr[2]*in->x[2];
	if (in->iProject == DF_PROJECT_PERSPECTIVE) x/=z;
	if (fabs(x)<in->xlim+hint) {
	    y = dr[0]*in->y[0] + dr[1]*in->y[1] + dr[2]*in->y[2];
	    if (in->iProject == DF_PROJECT_PERSPECTIVE) y/=z;
	    if (fabs(y)<in->ylim+hint) {
		x = x+in->xlim;
		xp = x;
		y = in->ylim-y;  /* standard screen convention */
		yp = y;
		if (hint < 1) {
		    Image[ xp + yp*in->nxPix ].r += br0*col.r;
		    Image[ xp + yp*in->nxPix ].g += br0*col.g;
		    Image[ xp + yp*in->nxPix ].b += br0*col.b;
#ifdef DEBUGTRACK
		    if (fAge > 3) {
			trackcnt++;
			if ((trackcnt<DEBUGTRACK)) {
			    printf("%i %i %f %f %f %f %f*\n",xp,yp,col.r,col.g,col.b,br0, fAge );
			    }
			}
#endif
		    }
		else {
		    int xpmin,xpmax,ypmin,ypmax,ix,iy;
		    DFIMAGE *Imagey;
		    double br,br1,r2,ih2;
		    ih2 = 1./(h*h);
		    br1 = br0*(6/(2.0*3.1412))*ih2;
		    xpmin = xp - hint; if (xpmin<0) xpmin=0;
		    xpmax = xp + hint; if (xpmax>=in->nxPix) xpmax=in->nxPix-1;
		    ypmin = yp - hint; if (ypmin<0) ypmin=0;
		    ypmax = yp + hint; if (ypmax>=in->nyPix) ypmax=in->nyPix-1;
		    for (iy=ypmin,Imagey = Image + iy*in->nxPix;iy<=ypmax;iy++,Imagey += in->nxPix) {
			for (ix=xpmin;ix<=xpmax;ix++) {
			    r2 = ((ix-x)*(ix-x)+(iy-y)*(iy-y))*ih2;
			    if (r2 > 1) continue;
			    br = br1*(1.0-sqrt(r2));
			    Imagey[ ix ].r += br*col.r;
			    Imagey[ ix ].g += br*col.g;
			    Imagey[ ix ].b += br*col.b;
#ifdef DEBUGTRACK
			    if (fAge > 3) {
				trackcnt++;
				if ((trackcnt<DEBUGTRACK)) {
				    printf("%i %i %f %f %f\n",ix,iy,br*col.r,br*col.g,br*col.b );
				    }
				}
#endif
			    }
			}
		    }
		}
	    }
	}
	    

    }

void dfRenderParticleSolid( struct inDumpFrame *in, void *vImage, PKD pkd,
			    PARTICLE *p)
{
    DFIMAGE *Image = vImage;
    DFIMAGE col; /* Colour */
    double h;
    double x,y,z,dr[3],br0;
    int hint;
    int j;
    int xp,yp;
    int iActive=p->iActive;
    double *r=p->r;
    double fMass=p->fMass;
    double fSoft=p->fSoft;
    double fBall2=p->fBall2;
    double fDensity=p->fDensity;
    double fAge=(in->dTime-p->fTimeForm)*in->dYearUnit;
		
    br0=1;

    if (pkdIsGas(pkd,p)) {
	if (fMass < in->dMassGasMin || fMass > in->dMassGasMax) return;
	if (in->bGasSph) h = sqrt(fBall2)*0.5*in->dGasSoftMul;
	else h = fSoft*in->dGasSoftMul;
	col = in->ColGas[0];
	}
    else if (pkdIsDark(pkd,p)) {
	if (fMass < in->dMassDarkMin || fMass > in->dMassDarkMax) return;
	h = fSoft*in->dDarkSoftMul;
	col = in->ColDark[0];
	}
    else if (pkdIsStar(pkd,p)) {
	if (fMass < in->dMassStarMin || fMass > in->dMassStarMax) return;
	h = fSoft*in->dStarSoftMul;
	switch (in->iColStarAge) {
	case DF_STAR_AGE_BRIGHT_COLOUR:
	case DF_STAR_AGE_COLOUR:
	    {
		float al;
		al = (log(fabs(fAge+1e6))-13.815)*(1/9.6);
		col.b = in->ColStar[0].b*fabs(1-0.7*al);
		col.g = in->ColStar[0].g*0.4;
		col.r = in->ColStar[0].r*(0.4+0.32*al);
		/*		  printf("star: %g %g %g %g %g  %g %g %g\n",fAge,al,col.r,col.g,col.b,in->ColStar.r,in->ColStar.g,in->ColStar.b); */
		}
	    break;
	default:
	    col = in->ColStar[0];
	    }
	switch (in->iColStarAge) {
	case DF_STAR_AGE_BASIC:
	    break;
	case DF_STAR_AGE_BRIGHT:
	case DF_STAR_AGE_BRIGHT_COLOUR:
	    br0 = 1./(fabs(fAge)/1e6+1.);
	    break;
	    }
	}
	
    if (in->bColMassWeight) br0*=fMass;
	
    for (j=0;j<3;j++) {
	dr[j] = r[j]-in->r[j];
	}
    z = dr[0]*in->z[0] + dr[1]*in->z[1] + dr[2]*in->z[2] + in->zEye;
    if (z >= in->zClipNear && z <= in->zClipFar) {
	if (in->iProject == DF_PROJECT_PERSPECTIVE) h = h*in->hmul/z;
	else h = h*in->hmul;
	hint = h;
	x = dr[0]*in->x[0] + dr[1]*in->x[1] + dr[2]*in->x[2];
	if (in->iProject == DF_PROJECT_PERSPECTIVE) x/=z;
	if (fabs(x)<in->xlim+hint) {
	    y = dr[0]*in->y[0] + dr[1]*in->y[1] + dr[2]*in->y[2];
	    if (in->iProject == DF_PROJECT_PERSPECTIVE) y/=z;
	    if (fabs(y)<in->ylim+hint) {
		xp = x+in->xlim;
		yp = in->ylim-y; /* standard screen convention */
		if (hint < 1) {
		    double br = 0.523599*br0*h*h; /* integral of TSC to h */
		    Image[ xp + yp*in->nxPix ].r += br*col.r;
		    Image[ xp + yp*in->nxPix ].g += br*col.g;
		    Image[ xp + yp*in->nxPix ].b += br*col.b;
		    }
		else {
		    int xpmin,xpmax,ypmin,ypmax,ix,iy;
		    DFIMAGE *Imagey;
		    double br,r2,ih2;
		    ih2 = 1./(h*h);
		    xpmin = xp - hint; if (xpmin<0) xpmin=0;
		    xpmax = xp + hint; if (xpmax>=in->nxPix) xpmax=in->nxPix-1;
		    ypmin = yp - hint; if (ypmin<0) ypmin=0;
		    ypmax = yp + hint; if (ypmax>=in->nyPix) ypmax=in->nyPix-1;
		    for (iy=ypmin,Imagey = Image + iy*in->nxPix;iy<=ypmax;iy++,Imagey += in->nxPix) {
			for (ix=xpmin;ix<=xpmax;ix++) {
			    if (ix==xp && iy==yp) {
				br = br0*(1.57080-1.04720/h); /* Integral of tsc disc to r=1 */
				}
			    else {
				r2 = ((ix-xp)*(ix-xp)+(iy-yp)*(iy-yp))*ih2;
				if (r2 > 1) continue;
				br = br0*(1.0-sqrt(r2));
				}
			    Imagey[ ix ].r += br*col.r;
			    Imagey[ ix ].g += br*col.g;
			    Imagey[ ix ].b += br*col.b;
			    }
			}
		    }
		}
	    }
	}
    }


void pkdDumpFrame( PKD pkd, struct inDumpFrame *in, void *vImage)
{
    int i;
    PARTICLE *p;
    int n = pkdLocal(pkd);

    for(i = 0; i < n; ++i) {
        p = &pkd->pStore[i];
	switch (in->iRender) {
	case DF_RENDER_POINT: 
	    dfRenderParticlePoint(in, vImage, pkd, p);
	    break;
	case DF_RENDER_TSC:
	    dfRenderParticleTSC( in, vImage, pkd, p);
				 break;
	case DF_RENDER_SOLID:
	    dfRenderParticleSolid( in, vImage, pkd, p);
	    break;
	case DF_RENDER_SHINE: /* Not implemented -- just does point */
	    dfRenderParticlePoint( in, vImage, pkd, p);
	    break;
	    }
	}
    }


void dfMergeImage( struct inDumpFrame *in, void *vImage1, int *nImage1, void *vImage2, int *nImage2 ) {
    int i;
    DFIMAGE *Image1 = vImage1, *Image2 = vImage2;

    switch (in->iRender) {
    case DF_RENDER_POINT:
    case DF_RENDER_TSC:
    case DF_RENDER_SOLID:
	assert( *nImage1 == in->nxPix*in->nyPix*sizeof(DFIMAGE) );
	assert( *nImage1 == *nImage2 );

	for (i=in->nxPix*in->nyPix-1;i>=0;i--) {
	    Image1[i].r += Image2[i].r;
	    Image1[i].g += Image2[i].g;
	    Image1[i].b += Image2[i].b;
	    }
	break;
    case DF_RENDER_SHINE:
	assert( *nImage1 == in->nxPix*in->nyPix*sizeof(DFIMAGE) );
	assert( *nImage1 == *nImage2 );

	for (i=in->nxPix*in->nyPix-1;i>=0;i--) {
	    if (Image2[i].r > Image1[i].r) Image1[i].r = Image2[i].r;
	    if (Image2[i].g > Image1[i].g) Image1[i].g = Image2[i].g;
	    if (Image2[i].b > Image1[i].b) Image1[i].b = Image2[i].b;
	    }
	break;
    default:
	break;
	}
    }

void dfFinishFrame( struct DumpFrameContext *df, double dTime, double dStep, struct inDumpFrame *in, void  *vImage ) {
    DFIMAGE *Image = vImage;
    char fileout[160];
    FILE *fp;
    int i;
    int iMax;
    unsigned char *gray,*g;
    /*char number[40]; -- not used: DCR 12/19/02*/

#ifdef DEBUGTRACK
    printf("ntracked %d\n",ntrack);
#endif

    switch( df->iNumbering ) {
    case DF_NUMBERING_FRAME:
	sprintf(fileout, df->FileName, df->nFrame );
	break;
    case DF_NUMBERING_STEP:
	sprintf(fileout, df->FileName, dStep );
	break;
    case DF_NUMBERING_TIME:
	sprintf(fileout, df->FileName, dTime );
	break;
	}

    df->nFrame++; /* NB: need to sort out something for restarts */
	
    iMax = in->nxPix*in->nyPix;
    gray = (unsigned char *) malloc(sizeof(unsigned char)*3*iMax);
    assert( gray != NULL );

    if (in->iRender == DF_RENDER_POINT) {
	for (i=0,g=gray;i<iMax;i++) {
	    if ( Image[i].r > 0.0 ) {
		*g=255; g++; 
		*g=255; g++; 
		*g=255; g++; 
		}
	    else {
		*g=0; g++; 
		*g=0; g++; 
		*g=0; g++; 
		}
	    }
	}
    else if (in->iRender == DF_RENDER_TSC || in->iRender == DF_RENDER_SOLID || in->iRender == DF_RENDER_SHINE) {
	if (in->iLogScale == DF_LOG_COLOURSAFE) {
	    double lmin,factor;
	    lmin = log(in->pScale1);
	    factor = 255.999/(log(in->pScale2)-lmin);
	    for (i=0,g=gray;i<iMax;i++) {
		int bing;
		float tot = Image[i].r;
		if (Image[i].g > tot) tot = Image[i].g;
		if (Image[i].b > tot) tot = Image[i].b;
		if (tot <= 0) {
		    *g = 0; g++;
		    *g = 0; g++;
		    *g = 0; g++;
		    }
		else {
		    tot = factor/tot*(log(tot)-lmin);
		    bing = Image[i].r*tot;
		    *g = (bing < 255 ? (bing < 0 ? 0 : bing) : 255 );
		    g++;
		    bing = Image[i].g*tot;
		    *g = (bing < 255 ? (bing < 0 ? 0 : bing) : 255 );
		    g++;
		    bing = Image[i].b*tot;
		    *g = (bing < 255 ? (bing < 0 ? 0 : bing) : 255 );
		    g++;
		    }
		}
	    }
	else if (in->iLogScale == DF_LOG_SATURATE) {
	    double lmin,factor;
	    lmin = log(in->pScale1);
	    factor = 255.999/(log(in->pScale2)-lmin);
	    for (i=0,g=gray;i<iMax;i++) {
		int bing;
		if (Image[i].r <= 0) *g = 0;
		else { 
		    bing = factor*(log(Image[i].r)-lmin);
		    *g = (bing < 255 ? (bing < 0 ? 0 : bing) : 255 );
		    }
		g++;
		if (Image[i].g <= 0) *g = 0;
		else { 
		    bing = factor*(log(Image[i].g)-lmin);
		    *g = (bing < 255 ? (bing < 0 ? 0 : bing) : 255 );
		    }
		g++;
		if (Image[i].b <= 0) *g = 0;
		else { 
		    bing = factor*(log(Image[i].b)-lmin);
		    *g = (bing < 255 ? (bing < 0 ? 0 : bing) : 255 );
		    }
		g++;
		}
	    }
	else {
	    for (i=0,g=gray;i<iMax;i++) {
		int bing;
		bing = 260*(1.-1./(0.1*Image[i].r+1));
		*g = (bing < 255 ? bing : 255 );
		g++;
		bing = 260*(1.-1./(0.1*Image[i].g+1));
		*g = (bing < 255 ? bing : 255 );
		g++;
		bing = 260*(1.-1./(0.1*Image[i].b+1));
		*g = (bing < 255 ? bing : 255 );
		g++;
		}
	    }
	}
    fp = fopen(fileout,"w");
    assert(fp!=NULL);

    if (df->iEncode == DF_ENCODE_PPM) {
	fprintf(fp,"P6\n#T=%20.10f\n%5i %5i\n255\n",dTime,in->nxPix,in->nyPix);
	fwrite( gray, 3*iMax, sizeof(char), fp);
	}

    else if (df->iEncode == DF_ENCODE_PNG) {
#ifdef USE_PNG
	static mainprog_info wpng_info;
	int rowbytes;
	int iErr;
		
	wpng_info.outfile = fp;
	wpng_info.infile = NULL; 
	wpng_info.image_data = gray; 
	wpng_info.pnmtype = 6; /* RGB */

	wpng_info.sample_depth = 8;  /* <==> maxval 255 */

	wpng_info.width = in->nxPix;
	wpng_info.height = in->nyPix;
	rowbytes = 3*wpng_info.width; /* 8 bit RGB */
	wpng_info.row_pointers = (uch **)malloc(wpng_info.height*sizeof(uch *)); 

	wpng_info.filter = FALSE; 
	wpng_info.interlaced = FALSE; 
	wpng_info.have_bg = FALSE; 
	wpng_info.have_time = FALSE; 
	wpng_info.have_text = 0; 
	wpng_info.gamma = 0.0; 

	iErr =  writepng_init(&wpng_info);
	assert (!iErr);

	for (i = 0;  i < wpng_info.height;  ++i) 
            wpng_info.row_pointers[i] = wpng_info.image_data + i*rowbytes; 
	
        iErr = writepng_encode_image(&wpng_info);
	assert (!iErr);

	writepng_cleanup(&wpng_info);
#endif
	}

    free( gray );
    fclose(fp);

    if (df->dDumpFrameTime > 0 && dTime >= df->dTime)
        df->dTime = df->dTime + df->dDumpFrameTime;

    if (df->dDumpFrameStep > 0 && dStep >= df->dStep) 
        df->dStep = df->dStep + df->dDumpFrameStep;

    }
  
#endif /*GSS_DUMPFRAME*/
