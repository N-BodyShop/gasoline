
#ifndef SSIO_HINCLUDED
#define SSIO_HINCLUDED

/*
 ** ssio.h -- DCR 98-09-16
 ** ------
 ** Header file specific to Solar System data I/O routines.
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/param.h>	/* for MAXPATHLEN */
#include <rpc/rpc.h>	/* for XDR routines */
#include "mode.h"

#ifndef MAXPATHLEN
#define MAXPATHLEN 256
#endif

#ifndef N_DIM
#define N_DIM 3
#endif

/*
 ** Following structures intended for use with xdr I/O routines. Note that
 ** internal representations may differ from format written to disk. As a
 ** result, for now use hardwired sizes for seeking. Is there an xdr routine
 ** to figure this out automatically?...
 */

typedef struct ss_head {
	double	time;
	int		n_data;
	int		pad;	/* unused: pad to 8-byte boundary */
	} SSHEAD;

#define SSHEAD_SIZE 16	/* XDR assumes 8-byte doubles and 4-byte ints */

typedef struct ss_data {
	double	mass;
	double	radius;
	double	pos[N_DIM];
	double	vel[N_DIM];
	double	spin[N_DIM];
	int		color;
	int		org_idx;
	} SSDATA;

#define SSDATA_SIZE 96	/* XDR assumes 8-byte doubles and 4-byte ints */

typedef struct ssio {
	FILE *fp;
	XDR xdrs;
	} SSIO;

#define SSIO_READ	0
#define SSIO_WRITE	1
#define SSIO_UPDATE	2

#define SS_EXT ".ss"

int ssioNewExt(const char *infile, const char *inext,
			   char *outfile, const char *outext);
int ssioOpen(const char *filename, SSIO *ssio, u_int mode);
int ssioHead(SSIO *ssio, SSHEAD *head);
int ssioData(SSIO *ssio, SSDATA *data);
int ssioClose(SSIO *ssio);
int ssioSetPos(SSIO *ssio, u_int pos);
void ssioRewind(SSIO *ssio);

#endif
