
#ifndef TREEZIPKEY_HINCLUDED
#define TREEZIPKEY_HINCLUDED

#include "treeziptypes.h"

/* set this outside using -DTZKEY64 or in some header file for 64 bit keys:
   #define TZKEY64
*/

/* 64 bit key (default is 128 bit -- slower but more precision possible) */
#ifdef TZKEY64

/* bits precision for each direction <= 63 
   Has to be 1/3 bits in key
   21 is equivalent to 2.09e6 cells across 
*/
#define TZNBITS_PER_DIRECTION  21

/* key is 64 bit unsigned int 
   We use only 63 = 21*3 of them 
*/
#define TZNBITS_IN_KEY     (TZNBITS_PER_DIRECTION*3)

typedef TZ_UINT64 tzkey;

#define TZKEY_AND(a,b) { a&=b; }
#define TZKEY_RET_LOW64(a) (a)
#define TZKEY_RET_AND1(a) (a&((TZ_UINT64) 1))
#define TZKEY_RET_ANDINT(a,i) (a&i)

#define TZKEY_ZERO(a) { a=0; }
#define TZKEY_SETINT(a,i) { a=i; }

#define TZKEY_SETBIT(a,n) { a=(((TZ_UINT64) 1)<<n); }
#define TZKEY_TESTBIT(a,n) { a&=(((TZ_UINT64) 1)<<n); }
#define TZKEY_ORBIT(a,n) { a|=(((TZ_UINT64) 1)<<n); } 

#define TZKEY_ORINT(a,i) { a|=i; }
#define TZKEY_OR(a,b) { a|=b; }
#define TZKEY_INTOR(i,a) { i|=a; }

#define TZKEY_RSHIFT1(a) { a>>=1; }
#define TZKEY_LSHIFT1(a) { a<<=1; }

/* n must be < 64 */
#define TZKEY_RSHIFT_LT64(a,n) { a>>=n; }
#define TZKEY_LSHIFT_LT64(a,n) { a<<=n; }
#define TZKEY_RSHIFT(a,n) { a>>=n; }
#define TZKEY_LSHIFT(a,n) { a<<=n; }

#define TZKEY_PRINTRULER( out ) { int j; \
				   for (j=63;j>=0;j--) fprintf(out,"%1i",(j/10)%10 ); \
                   fprintf(out,"\n"); \
				   for (j=63;j>=0;j--) fprintf(out,"%1i",j%10 ); \
                   fprintf(out,"\n"); \
                }
#define TZKEY_PRINTKEY( out, a ) { int j; \
				   for (j=63;j>=0;j--) fprintf(out,"%1i",(a>>j)&1); \
                }

#define TZKEY_ASSERT_CORRECT_SHIFTING() { \
	assert( sizeof(TZ_UINT16) == 2 ); \
	assert( sizeof(TZ_UINT32) == 4 ); \
	assert( sizeof(TZ_UINT64) == 8 ); \
										  }

#else 

/* bits precision for each direction <= 63 
   Has to be 1/3 bits in key
   42 is equivalent to 4.39e12 cells across 
*/
#define TZNBITS_PER_DIRECTION  42

/* key is a 128 bit unsigned int (emulated) 
   We use only 126 = 42*3 of them 
*/
#define TZNBITS_IN_KEY     (TZNBITS_PER_DIRECTION*3)

typedef struct tzKeyStruct {
	TZ_UINT64 i1,i0;
	} tzkey;
/* typedef unsigned long long key; */

#define TZKEY_AND(a,b) { a.i0&=b.i0; a.i1&=b.i1; }
#define TZKEY_RET_LOW64(a) (a.i0)
#define TZKEY_RET_AND1(a) (a.i0&1)
#define TZKEY_RET_ANDINT(a,i) (a.i0&(i))

#define TZKEY_ZERO(a) { a.i0=0; a.i1=0; }
#define TZKEY_SETINT(a,i) { a.i0=(i); a.i1=0; }

#define TZKEY_SETBIT(a,n) { if (n<64) { a.i0=(((TZ_UINT64) 1)<<n); a.i1=0; } else { a.i0=0; a.i1=(((TZ_UINT64) 1)<<(n-64)); } }
#define TZKEY_TESTBIT(a,n) { if (n<64) { a.i0&=(((TZ_UINT64) 1)<<n); a.i1=0; } else { a.i0=0; a.i1&=(((TZ_UINT64) 1)<<(n-64)); } }
#define TZKEY_ORBIT(a,n) { if (n<64) { a.i0|=(((TZ_UINT64) 1)<<n); } else { a.i1|=(((TZ_UINT64) 1)<<(n-64)); } }

#define TZKEY_ORINT(a,i) { a.i0|=(i); }
#define TZKEY_OR(a,b) { a.i0|=b.i0; a.i1|=b.i1; }
#define TZKEY_INTOR(i,a) { i|=(a).i0; }
#define TZKEY_RSHIFT1(a) { a.i0>>=1; a.i0|=(a.i1<<63); a.i1>>=1; }
#define TZKEY_LSHIFT1(a) { a.i1<<=1; a.i1|=(a.i0>>63); a.i0<<=1; }

/* n must be < 64 */
#define TZKEY_RSHIFT_LT64(a,n) { a.i0>>=(n); a.i0|= (a.i1<<(64-(n)) ); a.i1>>=(n); }
#define TZKEY_LSHIFT_LT64(a,n) { a.i1<<=(n); a.i1|= (a.i0>>(64-(n)) ); a.i0<<=(n); }
/* n must be < 128 */
#define TZKEY_RSHIFT(a,n) { TZ_UINT32 ntz=n; if (ntz>=64) { a.i0 = a.i1; a.i1 = 0; ntz-=64; TZKEY_RSHIFT_LT64(a, ntz); } \
                                  else { TZKEY_RSHIFT_LT64(a,ntz) } }
#define TZKEY_LSHIFT(a,n) { TZ_UINT32 ntz=n; if (ntz>=64) { a.i1 = a.i0; a.i0 = 0; ntz-=64; TZKEY_LSHIFT_LT64(a, ntz); } \
								  else { TZKEY_LSHIFT_LT64(a,ntz) } }

#define TZKEY_PRINTRULER( out ) { int j; \
				   for (j=127;j>=0;j--) fprintf(out,"%1i",(j/10)%10 ); \
                   fprintf(out,"\n"); \
				   for (j=127;j>=0;j--) fprintf(out,"%1i",j%10 ); \
                   fprintf(out,"\n"); \
                }
#define TZKEY_PRINTKEY( out, a ) { int j; \
				   for (j=63;j>=0;j--) fprintf(out,"%1i",(a.i1>>j)&1); \
				   for (j=63;j>=0;j--) fprintf(out,"%1i",(a.i0>>j)&1); \
                }


#define TZKEY_ASSERT_CORRECT_SHIFTING() { \
						   tzkey key_shifting_test;  \
											  assert( sizeof(TZ_UINT16) == 2 ); \
											  assert( sizeof(TZ_UINT32) == 4 ); \
											  assert( sizeof(TZ_UINT64) == 8 ); \
						   assert( sizeof(key_shifting_test.i0) == 8 );  \
						   key_shifting_test.i0 = 63;  \
						   key_shifting_test.i1 = 0;  \
					       TZKEY_LSHIFT( key_shifting_test, 60 );  \
					       assert( key_shifting_test.i0 == ( ((TZ_UINT64) 63)<<60 ) );  \
					       assert( key_shifting_test.i1 == ( ((TZ_UINT64) 63)>>(64-60) ) );  \
					       TZKEY_RSHIFT( key_shifting_test, 59 );  \
					       assert( key_shifting_test.i0 == 63<<1 );  \
					       assert( key_shifting_test.i1 == 0 );  \
								  }


#endif

#endif
