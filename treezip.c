#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <rpc/types.h>
#include <rpc/xdr.h>
#include <rpc/rpc.h>
#include <assert.h>
#include "tipsydefs.h"

/* four bytes version string */
#define VERSION "1.3 \0"

/*
   The file treezipkey.h handles the keys
   You can do "faster" things with 64 bit key or go 
   deeper with 128 bit keys
*/
//#define TZKEY64
#include "treezip.h"

/*

   There are pathological cases of many particles in the same
   location to float precision (eg. from a tipsy file).
   A logical solution is sort equally into left/right nodes beyond some depth: nBitsMaxPrecision 
   This can be the full tree depth -- there is no info left in the keys there anyway
   This lets you to go arbitrarily far beyond eg. 63 bits for as many buckets as you need.

   Definition of nBitsMaxPrecision changed to be per direction*3
*/

TZX *tzInit( double dmin[3], double dmax[3], int nBits_nPerBucket, int nBits_Position, int nBitsMaxPrecision_PerDirection, int nBitsMinPrecision_PerDirection ) {
    TZX *tz;

	/* test that the very long int emulation is working correctly */
	TZKEY_ASSERT_CORRECT_SHIFTING();

	tz = (TZX *) malloc(sizeof(TZX));
	assert(tz!=NULL);

	tz->nParticle = 0;
	tz->nNode = 1;
	tz->nBucket = 1;
#ifdef DEPTHCHECK
	tz->maxdepth = 0;
#endif

	tz->nBits_PerDirection = TZNBITS_PER_DIRECTION;
	tz->iBigInt = (((TZ_UINT64) 1)<<tz->nBits_PerDirection);
	/* Allows particles to have identical keys without breaking the code */
	tz->nBitsMaxPrecision_PerDirection = nBitsMaxPrecision_PerDirection;
	tz->nBitsMaxPrecision = tz->nBitsMaxPrecision_PerDirection*3;
    assert(nBitsMaxPrecision_PerDirection*3 <= TZNBITS_IN_KEY);

	tz->nBitsMinPrecision_PerDirection = nBitsMinPrecision_PerDirection;
	tz->nBitsMinPrecision = tz->nBitsMinPrecision_PerDirection*3;
    assert(nBitsMinPrecision_PerDirection <= nBitsMaxPrecision_PerDirection);

	tz->nBits_nPerBucket = nBits_nPerBucket;
	tz->nPerBucket = 1<<nBits_nPerBucket;
	tz->iMaskBucket = tz->nPerBucket-1;

	tz->nBits_Position = nBits_Position;
	tz->iMaskPosition = (1 << nBits_Position)-1;
	tz->nBitsExtra_PerParticle = 0;
	
	tz->dmin[0] = dmin[0];
	tz->dmin[1] = dmin[1];
	tz->dmin[2] = dmin[2];
	tz->dmax[0] = dmax[0];
	tz->dmax[1] = dmax[1];
	tz->dmax[2] = dmax[2];
	tz->dIsize[0] = (tz->iBigInt-1)/(dmax[0]-dmin[0]);
	tz->dIsize[1] = (tz->iBigInt-1)/(dmax[1]-dmin[1]);
	tz->dIsize[2] = (tz->iBigInt-1)/(dmax[2]-dmin[2]);

	tz->cnodeblock = NULL;
	tz->cnode = NULL;
	tz->icnode = TZNCNODE; /* forces new allocation */
	
	tz->pblock = NULL;
	tz->plist = NULL;
	tz->iplist = TZNPLIST; /* forces new allocation */
	
	tz->nTotalbits = 0;
	tz->nWritebits = 0;
	tz->nParticlebits = 0;
	tz->nLabelbits = 0;

	tz->root.child = NULL;
	tz->root.n = 0;
	tz->root.p = tzParticleAllocation( tz, tz->nPerBucket);
	assert( tz->root.p != NULL );

	return tz;
	}

void tzFinalize( TZX *tz ) {
	tznodeblock *nb,*nbnew;
	tzparticleblock *pb, *pbnew;
	
	assert(tz!=NULL);

	nb = tz->cnodeblock;

	while(nb!=NULL) {
		nbnew = nb->nextblock;
		free(nb);
		nb = nbnew;
		}

	pb = tz->pblock;
	
	while(pb!=NULL) {
		pbnew = pb->nextblock;
		free(pb);
		pb = pbnew;
		}

	free(tz->root.p);
	free(tz);
	}


void tzEmptyTree( TZX *tz ) {
	tznodeblock *nb,*nbnew;
	tzparticleblock *pb, *pbnew;
	
	assert(tz!=NULL);

	nb = tz->cnodeblock;

	while(nb!=NULL) {
		nbnew = nb->nextblock;
		free(nb);
		nb = nbnew;
		}

	pb = tz->pblock;
	
	while(pb!=NULL) {
		pbnew = pb->nextblock;
		free(pb);
		pb = pbnew;
		}

	free(tz->root.p);

	/* Reset all the storage related counters */

	tz->nParticle = 0;
	tz->nNode = 1;
	tz->nBucket = 1;
#ifdef DEPTHCHECK
	tz->maxdepth = 0;
#endif

	tz->cnodeblock = NULL;
	tz->cnode = NULL;
	tz->icnode = TZNCNODE; /* forces new allocation */
	
	tz->pblock = NULL;
	tz->plist = NULL;
	tz->iplist = TZNPLIST; /* forces new allocation */
	
	tz->root.child = NULL;
	tz->root.n = 0;
	tz->root.p = tzParticleAllocation( tz, tz->nPerBucket);
	assert( tz->root.p != NULL );

	}


/* Note: storage freed using link list */
tznode *tzNodeAllocation( TZX *tz, int n) {
	tznode *nodetmp;
	tznodeblock *nb;

	if (tz->icnode+n > TZNCNODE) {
		nb = (tznodeblock *) malloc(sizeof(tznodeblock));
		assert(nb!=NULL);
		nb->nextblock = tz->cnodeblock;
        tz->cnodeblock = nb;
		tz->cnode = &(nb->nodedata[0]);
		tz->icnode = 0;
		}
	nodetmp = &tz->cnode[tz->icnode];
	tz->icnode += n;
	return nodetmp;
	}

/* Note: storage freed using link list */
tzparticle *tzParticleAllocation( TZX *tz, int n) {
	tzparticle *ptmp;
	tzparticleblock *pb;

	if (tz->iplist+n > TZNPLIST) {
		pb = (tzparticleblock *) malloc(sizeof(tzparticleblock));
		assert(pb!=NULL);
		pb->nextblock = tz->pblock;
        tz->pblock = pb;
		tz->plist = &(pb->particledata[0]);
		tz->iplist = 0;
		}
	ptmp = &tz->plist[tz->iplist];
	tz->iplist += n;
	return ptmp;
	}

void tzWriteBits( TZX *tz, int bits, int nbitwrite )
{
#ifdef DEBUG2
	printf("Data to write: %i  %i\n",bits,nbitwrite);
#endif
	tz->bitstream <<= nbitwrite;
	tz->bitstream |= bits;
	tz->nbits += nbitwrite;
	tz->nTotalbits += nbitwrite;
	while (tz->nbits>=8) {
		TZ_UINT8 out;
		tz->nbits -= 8;
		tz->nWritebits += 8;
		out = (tz->bitstream>>tz->nbits)&255;
#ifdef DEBUG2
		{
		int j;
		for (j=0;j<8;j++) {
			printf("%1i",(out>>j)&1);
			}
		printf(" %8i\n",nWritebits-8);
		}
#endif
		fwrite( &out, sizeof(TZ_UINT8), 1, tz->fpout );
		}
    }

void tzWriteNode( TZX *tz, tznode *c, int l ) {
	int bits, i, nBits_Position,iMaskPosition;

	tz->nNode--;
	if (c->p == NULL) {
		bits = 0;
		if (c->child[0].n) bits |=1;
		if (c->child[1].n) bits |=2;
#ifdef DEBUG2
		printf("bits %i\n",bits );
#endif
		tzWriteBits( tz, bits, 2 );
		if (bits&1) tzWriteNode( tz, &c->child[0], l+1 );
		if (bits&2) tzWriteNode( tz, &c->child[1], l+1 );
		return;
		}
	else {
/* testing label compression */
#if (0)
		LABELTYPE labelmin,labelmax;
		unsigned int j,k;

		labelmax = labelmin = c->p[0].label;
		for (i=1;i<c->n;i++) {
			if (labelmin > c->p[i].label) labelmin = c->p[i].label;
			else if (labelmax < c->p[i].label) labelmax = c->p[i].label;
			}
		j = labelmax-labelmin;
		k = 0;
		while (j) { k++; j>>=1; }
		if (k<17) tz->nLabelbits += 17+5+c->n*k; 
		else tz->nLabelbits += 5+c->n*k;
/*
		fprintf(stderr,"label: %i-%i, range %i, %i\n",labelmin,labelmax,labelmax-labelmin,k);
*/
#endif
		bits = 0;
#ifdef DEBUG2
		printf("bits %i\n",bits );
#endif
		tzWriteBits( tz, bits, 2 );

		assert(c->n <= tz->nPerBucket);
		tz->nBucket--;
#ifdef DEBUG2
		printf("nBucket %i iMaskBucket %i\n",c->n,iMaskBucket );
#endif
	    tzWriteBits( tz, c->n&tz->iMaskBucket, tz->nBits_nPerBucket );
		/* Min precision plays a role here 
		   Need to know depth to determine if min precision not yet reached */
		nBits_Position = tz->nBitsMinPrecision-l;
		if (nBits_Position < tz->nBits_Position) nBits_Position = tz->nBits_Position;
		iMaskPosition = (1 << nBits_Position)-1;
	
		for (i=0;i<c->n;i++) {
			tz->nParticle--;
			tz->nParticlebits+=nBits_Position;
//			printf("nParticlebits %i\n",nParticlebits);
			tzWriteBits( tz, TZKEY_RET_ANDINT(c->p[i].k,iMaskPosition), nBits_Position );
			}
		}
	}

void tzOutputFile( TZX *tz, FILE *fpout ) {
	assert(fpout != NULL);
	tz->fpout = fpout;
	}

void tzWriteHeader( TZX *tz ) {
	TZ_UINT8  tmp;
	TZ_UINT16 tmps;
	char version[5]=VERSION;
	float endian=1;
	
	fwrite(version,sizeof(char),4,tz->fpout); /* Version */
	fwrite(&endian,sizeof(float),1,tz->fpout); /* Endian checker */
	
	tmp = tz->nBits_Position;
	fwrite(&tmp,sizeof(TZ_UINT8),1,tz->fpout);
	tmp = tz->nBits_nPerBucket;
	fwrite(&tmp,sizeof(TZ_UINT8),1,tz->fpout);
	tmp = tz->nBitsMaxPrecision;
	fwrite(&tmp,sizeof(TZ_UINT8),1,tz->fpout);
	tmp = tz->nBitsMinPrecision;
	fwrite(&tmp,sizeof(TZ_UINT8),1,tz->fpout);

	fwrite(&tz->dmin[0],sizeof(double),1,tz->fpout);
	fwrite(&tz->dmin[1],sizeof(double),1,tz->fpout);
	fwrite(&tz->dmin[2],sizeof(double),1,tz->fpout);

	fwrite(&tz->dmax[0],sizeof(double),1,tz->fpout);
	fwrite(&tz->dmax[1],sizeof(double),1,tz->fpout);
	fwrite(&tz->dmax[2],sizeof(double),1,tz->fpout);
	
	tmps = tz->nBitsExtra_PerParticle;
	fwrite(&tmps,sizeof(TZ_UINT16),1,tz->fpout);

	fwrite(&tz->nParticle,sizeof(int),1,tz->fpout);

	

#ifdef DEBUG
	fprintf(stderr,"HEADER: Version \"%4s\", Endian test float: %f\n"
			"n %i dmin %f %f %f dmax %f %f %f\n"
			"nBitsPosition %i nBitsPerBucket %i\n",
			&version[0],endian,tz->nParticle,
			tz->dmin[0],tz->dmin[1],tz->dmin[2],
			tz->dmax[0],tz->dmax[1],tz->dmax[2],
			tz->nBits_Position,
			tz->nBits_nPerBucket );
#endif

	}

void tzWriteTreeZip( TZX *tz ) {
	tznode *c ;

	/* Initialize bitstream */
	tz->nbits = 0;
	tz->bitstream = 0;

	tzWriteNode( tz, &tz->root, 0 );

	/* flush bitstream */
	tzWriteBits( tz, 0, 7 );
	}

void tzAddPos( TZX *tz, double *r, LABELTYPE label ) {
    tzkey k;
	TZ_UINT64 x,y,z;  /* 64 bit unsigned integers */
	int i;
	tzparticle p;
	tznode *c ;
#ifdef DEPTHCHECK
	int depth;
#endif

	assert( r[0] >= tz->dmin[0] );
	assert( r[1] >= tz->dmin[1] );
	assert( r[2] >= tz->dmin[2] );
	x = (((double) r[0]) - tz->dmin[0])*tz->dIsize[0];
	y = (((double) r[1]) - tz->dmin[1])*tz->dIsize[1];
	z = (((double) r[2]) - tz->dmin[2])*tz->dIsize[2];

	assert(x < tz->iBigInt);
	assert(y < tz->iBigInt);
	assert(z < tz->iBigInt);

	/* Make a TZNBITS_IN_KEY (see treezipkey.h for exact number) bit interleaved key */    

        TZKEY_ZERO(k);
	for (i=0;i<TZNBITS_PER_DIRECTION;i++) {
		TZKEY_LSHIFT1(k); 
		TZKEY_ORINT(k,(x&1));  x>>=1;
		TZKEY_LSHIFT1(k); 
		TZKEY_ORINT(k,(y&1));  y>>=1;
		TZKEY_LSHIFT1(k); 
		TZKEY_ORINT(k,(z&1));  z>>=1;
		}

	tz->nParticle++;

#ifdef DEPTHCHECK
	depth = 0;
#endif
	c = &tz->root;
	for (;;) {
		if (c->p != NULL) {
			/* We are at a leaf (bucket) */
			if (c->n < tz->nPerBucket) {
				/* Bucket isn't full -- add this particle to it and exit */
				p.k = k;
				p.label = label;
#ifdef DEBUG
				p.pos[0] = r[0];
				p.pos[1] = r[1];
				p.pos[2] = r[2];
#endif
				c->p[c->n] = p;
				c->n++;
#ifdef DEPTHCHECK
				if (depth > tz->maxdepth) tz->maxdepth = depth;
#endif
				break;
				}
			else {
				/* Bucket is full -- turn this leaf node into interior node
				   and create two child leaf nodes 
				   then continue as for interior node below...
				   */
				tzparticle *p0, *p1, *pend;
#ifdef DEPTHCHECK
/* This is not a problem as long as nBitsMaxPrecision is less than TZBITS_IN_KEY -- added an assert earlier */
/*
				if (depth > TZNBITS_IN_KEY) {
					int j;
					fprintf(stderr,"Hit max depth %i %i\n",depth,TZNBITS_IN_KEY);
					fprintf(stderr,"New particle: %15.12f %15.12f %15.12f\n",r[0],r[1],r[2]);
#ifdef DEBUG
					for (j=0;j<c->n;j++) {
						fprintf(stderr,"Bucket %3i: %15.12f %15.12f %15.12f\n",j,c->p[j].pos[0],c->p[j].pos[1],c->p[j].pos[2]);
						}
#endif
					assert(0);
					}
*/
#endif

				assert(c->child == NULL);
				c->child = tzNodeAllocation( tz, 2);
				assert(c->child != NULL);
				tz->nNode+= 2;
				c->child[0].child = NULL;
				c->child[0].n = 0;
				c->child[0].p = c->p;
				c->child[1].child = NULL;
				c->child[1].n = 0;
				c->child[1].p = tzParticleAllocation( tz, tz->nPerBucket);
				assert(c->child[1].p != NULL);
				tz->nBucket++;
				
				pend = c->p + c->n-1;
				p0 = c->child[0].p;
				p1 = c->child[1].p;

#ifdef DEPTHCHECK
				if (depth > tz->nBitsMaxPrecision) {
					/* Probably a cluster of near identical positions 
					   -- can't separate so sort equally into sub-bins 
					   */
					while (p0 <= pend) {
#ifdef DEBUG
						fprintf(stderr,"%i: depth %i, sort: %i\n",tz->nParticle,depth,(pend-p0)&1);
#endif
						if ((pend-p0)&1) {
							TZKEY_RSHIFT1(p0->k);
							*p1 = *p0;
							p1++;
							c->child[1].n++;
							*p0 = *pend;
							*pend--;
							}
						else {
							TZKEY_RSHIFT1(p0->k);
							p0++;
							c->child[0].n++;
							}
						}
					}
				else 
#endif
					{
					while (p0 <= pend) {
						if (TZKEY_RET_AND1(p0->k)) {
							TZKEY_RSHIFT1(p0->k);
							*p1 = *p0;
							p1++;
							c->child[1].n++;
							*p0 = *pend;
						*pend--;
							}
						else {
							TZKEY_RSHIFT1(p0->k);
							p0++;
							c->child[0].n++;
							}
						}
					}
				assert(c->child[0].n + c->child[1].n == c->n);
				c->p = NULL;
				}
			}

		assert(c->n == tz->nPerBucket);
		assert(c->child != NULL);

		/* Go down one level of the tree */
		c = &c->child[ TZKEY_RET_AND1(k) ];
		TZKEY_RSHIFT1(k);
#ifdef DEPTHCHECK
		depth++;
#endif

		}
		

//	printf("x y z %g %g %g %i %i %i  k %lli %g\n",r[0],r[1],r[2],x,y,z,k,log(k*1.0)/log(2.0), k&1);
	}


