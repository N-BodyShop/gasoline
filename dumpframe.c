#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "pkd.h"
#include "dumpframe.h"

/* #define IORDER_BLUETOORANGE (32768+32768) */
#define IORDER_BLUETOORANGE 2000000000


void dfInitialize( struct DumpFrameContext **pdf, double dTime, 
				  double dDumpFrameInterval, double dDelta, int iMaxRung, char *filename ) {
	double tock;

	(*pdf)->iMaxRung = 0;
	if (dDumpFrameInterval > 0) { 
		(*pdf)->dDumpFrameInterval = dDumpFrameInterval;
		(*pdf)->dTime = dTime-dDumpFrameInterval*1e-5;
		/* Round to a nearby power of two */
		tock = -log((dDumpFrameInterval-floor(dDumpFrameInterval))/dDelta*0.50001+1e-60)
			/log(2.0);
		if (tock <= iMaxRung) {
			(*pdf)->iMaxRung = tock;
			}

		dfParseCameraDirections( *pdf, filename );

		printf("DF Initialized Frame dumping: Time Interval %g [%g] (Step %i Rung %i)\n",dDumpFrameInterval,floor(dDumpFrameInterval)+dDelta*pow(2.0,-floor(tock)),(int) floor(dDumpFrameInterval/dDelta),(*pdf)->iMaxRung);
		}

	}

void dfFinalize( struct DumpFrameContext *df ) {
	free( df->fs );
	free( df );
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


#define DIFF( a, b, c ) { \
							   c[0] = b[0]-a[0]; \
							   c[1] = b[1]-a[1]; \
							   c[2] = b[2]-a[2]; \
						   }
#define DOT( a, b ) ( a[0]*b[0] + a[1]*b[1] + a[2]*b[2] )
#define DIST( a, b ) sqrt( (a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])+(a[2]-b[2])*(a[2]-b[2]) )
#define LEN( a ) sqrt( (a[0])*(a[0])+(a[1])*(a[1])+(a[2])*(a[2]) )
#define NORM( a ) { \
						double ilen; \
						ilen = 1./LEN ( a ); \
						a[0] *= ilen; a[1] *= ilen; a[2] *= ilen; \
					}
#define SIZE( a, f ) { \
						double ilen; \
						ilen = f/LEN ( a ); \
						a[0] *= ilen; a[1] *= ilen; a[2] *= ilen; \
					}
#define CURL( a, b, c ) { \
							  c[0] = a[1]*b[2] - a[2]*b[1]; \
							  c[1] = a[2]*b[0] - a[0]*b[2]; \
							  c[2] = a[0]*b[1] - a[1]*b[0]; \
						  }

void dfProjection( struct inDumpFrame *in, struct dfFrameSetup *fs ) {
	double width,height;
	double vec[3];
	double norm[3];

	in->dTime = fs->dTime;
    in->r[0] = fs->target[0];
    in->r[1] = fs->target[1];
    in->r[2] = fs->target[2];
	in->nxPix = fs->nxPix;
	in->nyPix = fs->nyPix;
	in->bPeriodic = fs->bPeriodic;
	in->iProject = fs->iProject;
	in->iEncode = fs->iEncode;
	in->iRender = fs->iRender;
	
	DIFF( fs->eye, in->r, in->z );
	if (fs->bzEye) {
		in->zEye = fs->zEye;
		}
	else {
		in->zEye = LEN( in->z );
		}

	assert( in->zEye > 0 );
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
		SIZE( in->x, (in->nxPix*0.5*in->zEye/width) );
		}
	else {
		SIZE( in->x, (in->nxPix*0.5/width) );
		}

	CURL( in->z, in->x, in->y );
    if (fs->iProject == DF_PROJECT_PERSPECTIVE) {
		SIZE( in->y, (in->nyPix*0.5*in->zEye/height) );
		}
	else {
		SIZE( in->y, (in->nyPix*0.5/height) );
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

	/* in->bNonLocal not set (an internal use only variable) */

	printf("DF Projection: %i x %i FOV %f  width %f height %f\n",fs->nxPix,fs->nyPix,fs->FOV,width,height);
	printf("DF eye %f %f %f, target %f %f %f, (separation %f)\n",fs->eye[0],fs->eye[1],fs->eye[2],fs->target[0],fs->target[1],fs->target[2],in->zEye );
    printf("DF up %f %f %f  Z-Clipping: Near %f Far %f\n",fs->up[0],fs->up[1],fs->up[2],in->zClipNear,in->zClipFar);
	printf("DF Vectors: x %f %f %f, y %f %f %f z %f %f %f\n",in->x[0],in->x[1],in->x[2], in->y[0],in->y[1],in->y[2], in->z[0],in->z[1],in->z[2] );

	}

void dfParseCameraDirections( struct DumpFrameContext *df, char * filename ) {
	FILE *fp;
	struct inDumpFrame in;
	struct dfFrameSetup fs;
	int n,nitem;
	char line[81],command[40],word[40];

	df->bLoop = 0;
	df->dTimeMod = 0;
	df->dTimeLoop = 1e20;
	df->dPeriodLoop = 0;
	
/* Defaults */
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
	fs.bzEye = 0;
	fs.FOV = 90.;
	fs.zClipNear = 0.01;
	fs.zClipFar = 2.0;
	fs.bzClipFrac = 1;
    fs.bPeriodic = 0;  /* Periodic? */
	fs.iProject = DF_PROJECT_PERSPECTIVE;
	fs.iEncode = DF_ENCODE_PPM;
	fs.iRender = DF_RENDER_POINT;

	fp = fopen( filename, "r" );
	if (fp==NULL) {
		printf("DF Could not open camera director file: %s\n",filename );

		df->iFrameSetup = 0;
		df->nFrameSetup = 1;
		df->fs = (struct dfFrameSetup *) malloc(sizeof(struct dfFrameSetup)*df->nFrameSetup);
		assert( df->fs != NULL );

		printf("DF Default Frame Setup\n" );
		dfProjection( &in, &fs ); /* Redundant now -- done again later */
		df->fs[0] = fs;
		return;
		}

	printf("DF Reading camera directions from: %s\n",filename );
	
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
	printf("DF Default Frame Setup\n" );
	for ( ;; ) {
		if (fgets( line, 81, fp ) == NULL) break;
		nitem = sscanf( line, "%s", command );
		if (nitem != 1 || command[0]=='#') continue;
		else if (!strcmp( command, "t" ) || !strcmp( command, "time" )) {
			dfProjection( &in, &fs ); /* Redundant now -- done again later */
			df->fs[n] = fs;
			n++; /* next setup */

			nitem = sscanf( line, "%s %lf", command, &fs.dTime );
			assert( nitem == 2 );
			
			printf("DF Frame Setup from time: %f\n",fs.dTime );
			}
		else if (!strcmp( command, "target") ) {
			nitem = sscanf( line, "%s %lf %lf %lf", command, &fs.target[0], &fs.target[1], &fs.target[2] );
			assert( nitem == 4 );
			}
		else if (!strcmp( command, "eye") ) {
			nitem = sscanf( line, "%s %lf %lf %lf", command, &fs.eye[0], &fs.eye[1], &fs.eye[2] );
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
		else if (!strcmp( command, "encode") ) {
			nitem = sscanf( line, "%s %s", command, word );
			assert( nitem == 2 );
		    if (!strcmp( word, "ppm") ) {
				fs.iEncode = DF_ENCODE_PPM;
			}
		    else if (!strcmp( word, "png") ) {
				fs.iEncode = DF_ENCODE_PNG;
#ifndef USE_PNG
				fprintf(stderr,"DF PNG encoding support not compiled in\n",word);
				assert(0);
#endif				
			}
		    else if (!strcmp( word, "rle") ) {
				fs.iEncode = DF_ENCODE_RLE;
				fprintf(stderr,"DF RLE encoding not supported yet\n",word);
				assert(0);
			}
			else {
				fprintf(stderr,"DF Unknown encoding: %s\n",word);
				assert( 0 );
				} 
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
			else {
				fprintf(stderr,"DF Unknown rendering: %s\n",word);
				assert( 0 );
				} 
			}
		else if (!strcmp( command, "zEye") ) {
			nitem = sscanf( line, "%s %s", command, word );
			assert( nitem == 2 );
		    if (!strcmp( word, "yes") || !strcmp( word, "on") ) {
				fs.bzEye = 1;
			}
		    else if (!strcmp( word, "no") || !strcmp( word, "off") ) {
				fs.bzEye = 0;
			}
			else {
				fprintf(stderr,"DF Unknown zEye setting: %s\n",word);
				assert( 0 );
				} 
			}
		else {
			fprintf(stderr,"DF Unknown command: %s\n",command);
			assert( 0 );
			}
		}
		
	dfProjection( &in, &fs ); /* Redundant now -- done again later */
	df->fs[n] = fs;
	n++; /* next setup */
	if (n==1) df->iFrameSetup = 0;
	assert ( n == df->nFrameSetup );
	}

#define DFINTERP_SCALAR( x, xx, xxx, a, b, c, d, f ) { \
f = xxx*a + xx*b + x*c + d; \
}

#define DFINTERP_VEC( x, xx, xxx, a, b, c, d, f ) { \
DFINTERP_SCALAR( x, xx, xxx, a[0], b[0], c[0], d[0], f[0] ); \
DFINTERP_SCALAR( x, xx, xxx, a[1], b[1], c[1], d[1], f[1] ); \
DFINTERP_SCALAR( x, xx, xxx, a[2], b[2], c[2], d[2], f[2] ); \
}

void dfSetupFrame( struct DumpFrameContext *df, double dTime, struct inDumpFrame *vin ) {
	struct dfFrameSetup fs;
	double x, xx, xxx;

	int ifs = df->iFrameSetup;

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
			fprintf(stderr,"DF Initial Time outside camera direction table: %g < %g\n",dTime,df->fs[1].dTime);
			df->iFrameSetup = ifs = 0;
			}
		else {
			ifs = 1;
			while (ifs < df->nFrameSetup && dTime < df->fs[ifs].dTime ) ifs++;
			if (ifs >= df->nFrameSetup-1) { /* Outside Range */
				if (dTime == df->fs[df->nFrameSetup-1].dTime && ifs > 1) ifs--;
				else {
					fprintf(stderr,"DF Initial Time outside camera direction table: %g > %g\n",dTime,df->fs[df->nFrameSetup-1].dTime);
					df->iFrameSetup = ifs = 0;
					}
				}
			df->iFrameSetup = ifs;
			dfGetCoeff( df, ifs );
			}
		}
	else if (ifs) {
		while (dTime > df->fs[ifs+1].dTime) { 
			ifs++;
			if (ifs >= df->nFrameSetup-1) {
				fprintf(stderr,"DF Time outside camera direction table: %g > %g\n",dTime,df->fs[df->nFrameSetup-1].dTime);
				df->iFrameSetup = ifs = 0;
				break;
				}
			df->iFrameSetup = ifs;
			dfGetCoeff( df, ifs );
			}
		}

	printf("DF Interpolating at t=%g Setups: %i (t=%g) %i (t=%g)\n",dTime,ifs,df->fs[ifs].dTime,ifs+1,df->fs[ifs+1].dTime);

	/* Nothing to interpolate? */
	if (!ifs) {
		dfProjection( vin, &df->fs[0] ); 
		return;
		}

	/* 
	   Interpolate Eye position, FOV etc... 
	   from df->fs[ifs].dTime <= dTime <= df->fs[ifs+1].dTime
	   */
	fs = df->fs[ifs];

	x = (dTime-fs.dTime)*df->rdt;
	xx = x*x;
	xxx = x*xx;

	DFINTERP_VEC( x, xx, xxx, df->a.eye, df->b.eye, df->c.eye, df->d.eye, fs.eye );
	DFINTERP_VEC( x, xx, xxx, df->a.target, df->b.target, df->c.target, df->d.target, fs.target );
	DFINTERP_VEC( x, xx, xxx, df->a.up, df->b.up, df->c.up, df->d.up, fs.up );
	DFINTERP_SCALAR( x, xx, xxx, df->a.FOV, df->b.FOV, df->c.FOV, df->d.FOV, fs.FOV );
	DFINTERP_SCALAR( x, xx, xxx, df->a.zClipNear, df->b.zClipNear, df->c.zClipNear, df->d.zClipNear, fs.zClipNear );
	DFINTERP_SCALAR( x, xx, xxx, df->a.zClipFar, df->b.zClipFar, df->c.zClipFar, df->d.zClipFar, fs.zClipFar );

    dfProjection( vin, &fs ); 
    }


void dfRenderImage( PKD pkd, struct inDumpFrame *in, void *vImage, int *nImage ) {
	PARTICLE *p;
	DFIMAGE *Image = vImage;
	int i,j;
	double x,y,z,dr[3];
	double xlim = (in->nxPix-1)*.5;
	double ylim = (in->nyPix-1)*.5;
	int xp,yp;
	DFIMAGE blank;
	blank.r = 0;
	blank.g = 0;
    blank.b = 0;

	*nImage = in->nxPix * in->nyPix * sizeof(DFIMAGE);

	for (i=in->nxPix*in->nyPix-1;i>=0;i--) Image[i] = blank;

	p = pkd->pStore;
	if (in->iRender == DF_RENDER_POINT) {
		for (i=0;i<pkd->nLocal;i++) {
			if (!TYPETest( &p[i], TYPE_STAR )) continue;
			for (j=0;j<3;j++) {
				dr[j] = p[i].r[j]-in->r[j];
				}
			z = dr[0]*in->z[0] + dr[1]*in->z[1] + dr[2]*in->z[2] + in->zEye;
			if (z >= in->zClipNear && z <= in->zClipFar) {
				x = dr[0]*in->x[0] + dr[1]*in->x[1] + dr[2]*in->x[2];
				if (in->iProject == DF_PROJECT_PERSPECTIVE) x/=z;
				if (fabs(x)<xlim) {
					y = dr[0]*in->y[0] + dr[1]*in->y[1] + dr[2]*in->y[2];
					if (in->iProject == DF_PROJECT_PERSPECTIVE) y/=z;
					if (fabs(y)<ylim) {
						xp = x+xlim;
						yp = ylim-y; /* standard screen convention */
						Image[ xp + yp*in->nxPix ].r += 1.0;
						Image[ xp + yp*in->nxPix ].g += 1.0;
						Image[ xp + yp*in->nxPix ].b += 1.0;
						}
					}
				}
			}
		}
	else if (in->iRender == DF_RENDER_TSC) {
		double hmul = 4*sqrt(in->x[0]*in->x[0] + in->x[1]*in->x[1] + in->x[2]*in->x[2]),h;
		int hint;

		for (i=0;i<pkd->nLocal;i++) {
			if (!TYPETest( &p[i], TYPE_STAR )) continue;
			for (j=0;j<3;j++) {
				dr[j] = p[i].r[j]-in->r[j];
				}
			z = dr[0]*in->z[0] + dr[1]*in->z[1] + dr[2]*in->z[2] + in->zEye;
			if (z >= in->zClipNear && z <= in->zClipFar) {
				if (in->iProject == DF_PROJECT_PERSPECTIVE) h = (p[i].fSoft*hmul/z);
				else h = (p[i].fSoft*hmul);
				hint = h;
				x = dr[0]*in->x[0] + dr[1]*in->x[1] + dr[2]*in->x[2];
				if (in->iProject == DF_PROJECT_PERSPECTIVE) x/=z;
				if (fabs(x)<xlim+hint) {
					y = dr[0]*in->y[0] + dr[1]*in->y[1] + dr[2]*in->y[2];
					if (in->iProject == DF_PROJECT_PERSPECTIVE) y/=z;
					if (fabs(y)<ylim+hint) {
						xp = x+xlim;
						yp = ylim-y; /* standard screen convention */
						if (hint < 1) {
							double br = 0.523599*h*h; /* integral of TSC to h */
							if (p[i].iOrder < IORDER_BLUETOORANGE) {
								Image[ xp + yp*in->nxPix ].r += br*0.5;
								Image[ xp + yp*in->nxPix ].g += br*0.7;
								Image[ xp + yp*in->nxPix ].b += br*1.0;
								}
							else {
								Image[ xp + yp*in->nxPix ].r += br*1.0;
								Image[ xp + yp*in->nxPix ].g += br*0.7;
								Image[ xp + yp*in->nxPix ].b += br*0.5;
								}
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
							if (p[i].iOrder < IORDER_BLUETOORANGE) {
								for (iy=ypmin,Imagey = Image + iy*in->nxPix;iy<=ypmax;iy++,Imagey += in->nxPix) {
									for (ix=xpmin;ix<=xpmax;ix++) {
										if (ix==xp && iy==yp) {
											br = 1.57080-1.04720/h; /* Integral of TSC to r=1 */
											}
										else {
											r2 = ((ix-xp)*(ix-xp)+(iy-yp)*(iy-yp))*ih2;
											if (r2 > 1) continue;
											br = 1.0-sqrt(r2);
											}
										Imagey[ ix ].r += 0.5*br;
										Imagey[ ix ].g += br*.7;
										Imagey[ ix ].b += br;
										}
									}
								}
							else {
								for (iy=ypmin,Imagey = Image + iy*in->nxPix;iy<=ypmax;iy++,Imagey += in->nxPix) {
									for (ix=xpmin;ix<=xpmax;ix++) {
										if (ix==xp && iy==yp) {
											br = 1.57080-1.04720/h; /* Integral of TSC to r=1 */
											}
										else {
											r2 = ((ix-xp)*(ix-xp)+(iy-yp)*(iy-yp))*ih2;
											if (r2 > 1) continue;
											br = 1.0-sqrt(r2);
											}
										Imagey[ ix ].r += br;
										Imagey[ ix ].g += br*.7;
										Imagey[ ix ].b += br*.5;
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}

void dfMergeImage( struct inDumpFrame *in, void *vImage1, int *nImage1, void *vImage2, int *nImage2 ) {
	int i;
	DFIMAGE *Image1 = vImage1, *Image2 = vImage2;

	assert( *nImage1 == in->nxPix*in->nyPix*sizeof(DFIMAGE) );
	assert( *nImage1 == *nImage2 );

	for (i=in->nxPix*in->nyPix-1;i>=0;i--) {
		Image1[i].r += Image2[i].r;
		Image1[i].g += Image2[i].g;
		Image1[i].b += Image2[i].b;
		}
	}

void dfFinishFrame( struct DumpFrameContext *df, double dTime, struct inDumpFrame *in, void  *vImage ) {
	DFIMAGE *Image = vImage;
	char fileout[80];
	FILE *fp;
	int i;
	int iMax;

	unsigned char *gray,*g;

	if (in->iEncode == DF_ENCODE_PPM) {
		sprintf(fileout,"%s%09i.ppm", "Frame", df->nFrame );
		}
	else if (in->iEncode == DF_ENCODE_PNG) {
		sprintf(fileout,"%s%09i.png", "Frame", df->nFrame );
		}
	else {
		sprintf(fileout,"%s%09i.null", "Frame", df->nFrame );
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
	else if (in->iRender == DF_RENDER_TSC) {

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

	fp = fopen(fileout,"w");
	assert(fp!=NULL);

	if (in->iEncode == DF_ENCODE_PPM) {
		fprintf(fp,"P6\n#T=%20.10f\n%5i %5i\n255\n",dTime,in->nxPix,in->nyPix);
		fwrite( gray, 3*iMax, sizeof(char), fp);
		}

	else if (in->iEncode == DF_ENCODE_PNG) {
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
	
	df->dTime = df->dTime + df->dDumpFrameInterval;

	}
