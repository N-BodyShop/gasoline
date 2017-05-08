#include "ssio.h" /* also defines MAXPATHLEN */
#ifdef COLLISIONS /* pkdgrav version should otherwise match with ss version */

/*
 ** ssio.c -- DCR 98-09-16
 ** ------
 ** Solar System data I/O routines.
 */

#include <string.h>
#include <assert.h>

static const char *
ssioBasename(const char *path)
{
	char *p;

	assert(path != NULL);
	p = strrchr(path,'/');
	if (p) return p + 1;
	else return path;
	}

int
ssioNewExt(const char *infile,const char *inext,
		   char *outfile,const char *outext)
{
	const char *basename;
	char *c;
	size_t n;

	assert(infile != NULL && inext != NULL && outfile != NULL && outext != NULL);
	basename = ssioBasename(infile);
	if ((c = strrchr(basename,'.')) && strstr(c,inext))
		n = c - basename;
	else
		n = strlen(basename);
	if (n + strlen(outext) >= (size_t) MAXPATHLEN)
		return 1;
	(void) strncpy(outfile,basename,n); /* not null terminated */
	(void) strcpy(outfile + n,outext);
	return 0;
	}

int
ssioOpen(const char *filename,SSIO *ssio,const u_int mode)
{
	const char type[][3] = {"r","w","r+"};
	const enum xdr_op op[] = {XDR_DECODE,XDR_ENCODE,XDR_ENCODE};

	assert(filename != NULL && ssio != NULL);
	assert(mode == SSIO_READ || mode == SSIO_WRITE || mode == SSIO_UPDATE);
	if (!(ssio->fp = fopen(filename,type[mode])))
		return 1;
	xdrstdio_create(&ssio->xdrs,ssio->fp,op[mode]);
	return 0;
	}

int
ssioHead(SSIO *ssio,SSHEAD *head)
{
	assert(ssio != NULL && head != NULL);
	if (!xdr_double(&ssio->xdrs,&head->time)) return 1;
	if (!xdr_int(&ssio->xdrs,&head->n_data)) return 1;
	if (!xdr_int(&ssio->xdrs,&head->pad)) return 1;
	return 0;
	}

int
ssioData(SSIO *ssio,SSDATA *data)
{
	int i;

	assert(ssio != NULL && data != NULL);
	if (!xdr_double(&ssio->xdrs,&data->mass)) return 1;
	if (!xdr_double(&ssio->xdrs,&data->radius)) return 1;
	for (i=0;i<N_DIM;i++)
		if (!xdr_double(&ssio->xdrs,&data->pos[i])) return 1;
	for (i=0;i<N_DIM;i++)
		if (!xdr_double(&ssio->xdrs,&data->vel[i])) return 1;
	for (i=0;i<N_DIM;i++)
		if (!xdr_double(&ssio->xdrs,&data->spin[i])) return 1;
	if (!xdr_int(&ssio->xdrs,&data->color)) return 1;
	if (!xdr_int(&ssio->xdrs,&data->org_idx)) return 1;
	return 0;
	}

int
ssioClose(SSIO *ssio)
{
	assert(ssio != NULL);
	xdr_destroy(&ssio->xdrs);
	if (fclose(ssio->fp)) return 1;
	return 0;
	}

int
ssioSetPos(SSIO *ssio,const u_int pos)
{
	assert(ssio != NULL);
	if (!xdr_setpos(&ssio->xdrs,pos)) return 1;
	return 0;
	}

void
ssioRewind(SSIO *ssio)
{
	assert(ssio != NULL);
	rewind(ssio->fp);
	}

/* ssio.c */

#endif
