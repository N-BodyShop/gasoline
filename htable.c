#include <stdlib.h>
#include "htable.h"

#define N	10000
#define M	100

/*
** This function returns the first prime greater than min.
** If min is greater than the Nth prime the function returns 0.
*/
int prime(int min)
{
    int i,k,x,inc,lim,square,prime;
    int p[M+1],v[M+1];

    x = 1;
    inc = 4;
    lim = 1;
    square = 9;
    for (i=3;i<=N;++i) {
        /* find the next prime number p[i] */
        prime = 0;
        while (!prime) {
            x += inc;
            inc = 6-inc;
            if (square <= x) {
                v[++lim] = square;
                square = p[lim+1] * p[lim+1];
                }
            k = 2;
            prime = 1;
            while (prime && (k < lim)) {
                ++k;
                if (v[k] < x) v[k] += 2*p[k];
                prime = (x != v[k]);
                }
            }
        if (x>min) return(x);
        if (i <= M) p[i] = x;
        }
    return(0);
    }


/*
** This function initializes a hash table.
** 'tbl_size' is the desired size of the table. The table
** sizes which are actually chosen by init are primes. So
** if you have passed a prime already pass 1 for 'isprime',
** otherwise pass 0 for 'isprime' and init will choose the
** prime just larger than the value you give for 'tbl_size'.
** 'key_compare' is a user supplied function which returns 0
** if the keys within two data elements match.
** 'hash_func' is a user supplied function takes a data element
** and produces a hash value from the key within the data.
*/ 
HTBL *HTBL_init(int tbl_size,int isprime,int (*key_compare)(void *,void *),
				unsigned int (*hash_func)(void *))
{
    HTBL *T;

    T = (HTBL *)malloc(sizeof(HTBL));
    T->key_compare = key_compare;
    T->hash_func = hash_func;
    if (isprime) T->tbl_size = tbl_size;
    else T->tbl_size = prime(tbl_size);
    if (!T->tbl_size) {
        exit(1);
        }
    T->tbl = (void **)calloc(T->tbl_size,sizeof(void *));
    T->keys = 0;
    return(T);
    }


void HTBL_finish(HTBL *T,void (*free_entry)(void *))
{
	int i;

	if (!T) return;
	if (free_entry != NULL) {
		for (i=0;i<T->tbl_size;++i) {
			if (T->tbl[i] != NULL) free_entry(T->tbl[i]);
			}
		}
	if (T->tbl) free(T->tbl);
	free(T);
	}


int HTBL_insert(HTBL *T,void *data)
{
    unsigned int ts = T->tbl_size,o,i,hv;

    if (T->keys == ts) return(0);

    hv = (*(T->hash_func))(data);
    o = hv%(ts-1)+1;
    i = hv%ts;
    while (T->tbl[i]) {
        /*
        ** Make sure the key is not already in the table.
        */
        if (!(*(T->key_compare))(data,T->tbl[i])) return(0);
        i = (i+o)%ts;
        }
    T->tbl[i] = data;
    ++(T->keys);
    return(1);
    }


void *HTBL_lookup(HTBL *T,void *data)
{
    unsigned int ts = T->tbl_size,o,i,hv;
    int c=0;

    hv = (*(T->hash_func))(data);
    o = hv%(ts-1)+1;
    i = hv%ts;
    while (T->tbl[i]) {
        if (!(*(T->key_compare))(data,T->tbl[i])) return(T->tbl[i]);
        if (++c == T->keys) return(0);
        i = (i+o)%ts;
        }
    return(0);
    }
