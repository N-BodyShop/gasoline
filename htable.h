#ifndef HTABLE_HINCLUDED

typedef struct hash_table {
    int tbl_size;
    int keys;
    int (*key_compare)(void *, void *);
    unsigned int (*hash_func)(void *);
    void **tbl;
    } HTBL;

HTBL *HTBL_init(int, int, int (*)(void *, void *),unsigned int (*)(void *));
void HTBL_finish(HTBL *,void (*)(void *));
int HTBL_insert(HTBL *, void *);
void *HTBL_lookup(HTBL *, void *);

#define HTABLE_HINCLUDED
#endif
