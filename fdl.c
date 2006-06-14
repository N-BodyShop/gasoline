#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#include <string.h>

#include "htable.h"
#include "fdl.h"

/* 
 ** The following two following constants were arbitrarily picked to
 ** put some randomness into the hash function.  The only real requirement
 ** is that SCATTER should be bigger than the hash table size.
 */
#define SCATTER	31856531
#define MAGIC	17


typedef struct tbl_entry {
    char *ident;
    FDL_NODE *point;
    int value;
    int curr_value;
    } ENTRY;


int FDL_error(FDL_CTX *CTX,int err,char *file,int line)
{
    static const char *message[] = FDL_MESSAGE_INITIALIZE;
    char c;
    
    if (!err) return(0);
    else if (err < FDL_MESSAGES) {
        c = message[err][0];
        printf("%c%03d %s(%d):%s\n",c,err,file,line,message[err]+1);
        if (c == 'E') {
            FDL_finish(CTX);
            exit(1);
            }
        if (c == 'W') return(err);
        if (c == 'I') return(0);
        FDL_finish(CTX);
        exit(1);
        }
    else {
        c = message[0][0];
        printf("%c%03d %s(%d):%s\n",c,0,file,line,message[0]+1);
        FDL_finish(CTX);
        exit(1);
        }
	assert(0); /* can't be here */
	return 0;
    }


FDL_NODE *FDL_newnode(FDL_CTX *CTX,int type)
{
    FDL_NODE *new;
    
    new = (FDL_NODE *)malloc(sizeof(FDL_NODE));
    if (!new)
        FDL_error(CTX,FDL_MSG_MEMNODE,__FILE__,__LINE__);
    new->parent = NULL;
    new->next = NULL;
    new->prev = NULL;
    new->child = NULL;
    new->type = type;
    new->ident = NULL;
    if (type == FLT_SIZE || type == FIX_SIZE || type == INDEX)
        new->size = 0;
    else new->size = -1;
    return(new);
    }


void FDL_killtree(FDL_NODE *p)
{
    FDL_NODE *q,*t;
    
    if (!p) return;
    q = p->child;
    while (q) {
        t = q->next;
        FDL_killtree(q);
        q = t;
        }
    if (p->type != FIX_SIZE) free(p->ident);
    free(p);
    }


void FDL_addcom(FDL_CTX *CTX,char *com,int comnum)
{
    int l;
    
    l = strlen(com);
    CTX->command[comnum] = (char *)malloc(l+1);
    if (!CTX->command[comnum])
        FDL_error(CTX,FDL_MSG_MEMCTX,__FILE__,__LINE__);
    strcpy(CTX->command[comnum],com);
    }


int key_compare(e1,e2)
void *e1,*e2;
{
    char *key1,*key2;
    
    key1 = ((ENTRY *)e1)->ident;
    key2 = ((ENTRY *)e2)->ident;
    return(strcmp(key1,key2));
    }


unsigned int hash_func(e)
void *e;
{
    unsigned char *p;
    unsigned int i;
    
    p = (unsigned char *) ((ENTRY *)e)->ident;
    i = SCATTER+MAGIC;
    while (*p) i += (*(p++))*SCATTER;
    return(i);
    }

void free_entry(e)
void *e;
{
	if (e != NULL) free(e);
	}


FDL_CTX *FDL_initialize(FILE *fp,FILE *fdlp,int bCopyFdl)
{
    FDL_CTX *CTX;
    
    CTX = (FDL_CTX *)malloc(sizeof(FDL_CTX));
    assert(CTX != NULL);
    /*
     ** Set up the FDL command strings.
     */
    FDL_addcom(CTX,"ELEMENT",ELEMENT);
    FDL_addcom(CTX,"ARRAY",ARRAY);
    FDL_addcom(CTX,"FDL",FDL);
    FDL_addcom(CTX,"TYPE",TYPE);
    FDL_addcom(CTX,"SIZE",FDL_SIZE);
    FDL_addcom(CTX,"INDEX",INDEX);
    
    CTX->fileptr = fp;
    CTX->fdlp = fdlp;
    CTX->identifiers = 0;
	CTX->bCopyFdl = bCopyFdl;
    return(CTX);
    }


void FDL_finish(FDL_CTX *CTX)
{
    int i;
    
    if (CTX->fdlp) fclose(CTX->fdlp);
    fclose(CTX->fileptr);
    for (i=0;i<NUM_COMMANDS;++i) free(CTX->command[i]);
    FDL_killtree(CTX->ast);
    FDL_killtree(CTX->trailer);
	/*
	 ** Note that HTBL_finish *only* removes the hash table and
	 ** None of the entries that were inserted into it. This is 
	 ** really poor design but we can make sure by forcing it giving 
	 ** a free_entry function.
	 */
    HTBL_finish(CTX->T,free_entry);
    free(CTX);
    }


int FDL_token(FDL_CTX *CTX)
{
    int ic;
    
    if (CTX->fdlp) {
        ic = fgetc(CTX->fdlp);
        if (ferror(CTX->fdlp))
            FDL_error(CTX,FDL_MSG_FILEREAD,__FILE__,__LINE__);
		if (CTX->bCopyFdl) {
			fputc(ic,CTX->fileptr);
			if (ferror(CTX->fileptr))
				FDL_error(CTX,FDL_MSG_FILEWRITE,__FILE__,__LINE__);
			}
        }
    else {
        ic = fgetc(CTX->fileptr);
        if (ferror(CTX->fileptr))
            FDL_error(CTX,FDL_MSG_FILEREAD,__FILE__,__LINE__);
        }
    if (isspace(ic)) ic = ' ';
    else if (ic == EOF) ic = 0;
    CTX->curr_token = ic;
    return(ic);
    }


int FDL_command(FDL_CTX *CTX)
{
    int tok,i;
    char *p;
    
    tok = CTX->curr_token;
    while (tok == ' ') tok = FDL_token(CTX);
    for (i=0;i<NUM_COMMANDS;++i) if (tok == CTX->command[i][0]) break;
    if (i == NUM_COMMANDS) {
        return(-1);
        }
    else {
        p = CTX->command[i];
        ++p;
        while (*p) {
            tok = FDL_token(CTX);
            if (tok != *p) return(-1);
            ++p;
            }
        tok = FDL_token(CTX);
        if (tok != ' ') return(-1);
        else return(i);
        }
    }


int isoneof(s,c)
char *s,c;
{
    while (*s) {
        if (*s == '@' && c == 0) return(1);
        if (*s == c) return(1);
        ++s;
        }
    return(0);
    }


int FDL_ident(FDL_CTX *CTX)
{
    int tok,count;
    char *p;
    
    count = 0;
    tok = CTX->curr_token;
    while (tok == ' ') tok = FDL_token(CTX);
    if (isalpha(tok)) {
        p = CTX->curr_ident;
        if (++count == MAX_IDENT)
            return(FDL_error(CTX,FDL_MSG_LONGIDENT,__FILE__,__LINE__));
        *p = tok;
        ++p;
        tok = FDL_token(CTX);
        while (isalnum(tok)||isoneof("_",tok)) {
            if (++count == MAX_IDENT)
                return(FDL_error(CTX,FDL_MSG_LONGIDENT,__FILE__,__LINE__));
            *p = tok;
            ++p;
            tok = FDL_token(CTX);
            }
        *p = 0;
        return(0);
        }
    else return(FDL_MSG_IDENT);
    }


int FDL_numint(FDL_CTX *CTX)
{
    int tok;
    
    tok = CTX->curr_token;
    while (tok == ' ') tok = FDL_token(CTX);
    if (isdigit(tok)) {
        CTX->curr_numint = tok - '0';
        tok = FDL_token(CTX);
        while (isdigit(tok)) {
            CTX->curr_numint *= 10;
            CTX->curr_numint += tok - '0';
            tok = FDL_token(CTX);
            }
        return(0);
        }
    else return(FDL_MSG_NUMINT);
    }


void FDL_copyident(FDL_CTX *CTX,char **string)
{
    int l;
    
    l = strlen(CTX->curr_ident);
    *string = (char *)malloc(l+1);
    if (*string == NULL)
        FDL_error(CTX,FDL_MSG_MEMIDENT,__FILE__,__LINE__);
    strcpy(*string,CTX->curr_ident);
    }


int FDL_sizeof(char *p)
{
    if (!strcmp(p,"float")) return(sizeof(float));
    if (!strcmp(p,"double")) return(sizeof(double));
    if (!strcmp(p,"int")) return(sizeof(int));
    if (!strcmp(p,"char")) return(sizeof(char));
    return(0);
    }


void FDL_trailer(FDL_CTX *CTX)
{
    FDL_NODE *new;
    
    new = FDL_newnode(CTX,TYPE);
    new->size = 0;
    CTX->ast->next = new;
    new->prev = CTX->ast;
    CTX->trailer = new;
    }


int FDL_expr(FDL_CTX *CTX,FDL_NODE **E,char *retsyms)
{
    int command,tok,err;
    FDL_NODE *p,*index,*size,*q;
    
    command = FDL_command(CTX);	
    if (command < 0)
        return(FDL_error(CTX,FDL_MSG_COMMAND,__FILE__,__LINE__));
    switch (command) {
        
    case ELEMENT:
    case ARRAY:
    case FDL:
        
        err = FDL_ident(CTX);
        if (err)
            return(FDL_error(CTX,err,__FILE__,__LINE__));
        *E = FDL_newnode(CTX,command);
        FDL_copyident(CTX,&(*E)->ident);
        ++CTX->identifiers;
        tok = CTX->curr_token;
        while (tok == ' ') tok = FDL_token(CTX);
        if (isoneof(retsyms,tok)) return(0);
        else if (tok == '(') {
            FDL_token(CTX);
            err = FDL_expr(CTX,&(*E)->child,")");
            if (err) return(FDL_error(CTX,err,__FILE__,__LINE__));
			/*
			 ** If it is an FDL node then we can stop processing here.
			 */
			if ((*E)->type == FDL) return(0);
            /*
             ** If it is an array node make sure there is an INDEX and
             ** a SIZE specified.
             */
            if ((*E)->type == ARRAY) {
                q = (*E)->child;
                index = NULL;
                size = NULL;
                while (q) {
                    if (q->type == FLT_SIZE || q->type == FIX_SIZE) size = q;
                    else if (q->type == INDEX) index = q;
                    q = q->next;
                    }
                if (!index)
                    return(FDL_error(CTX,FDL_MSG_NOINDEX,__FILE__,__LINE__));
                if (!size)
                    return(FDL_error(CTX,FDL_MSG_NOSIZE,__FILE__,__LINE__));
                }
            tok = FDL_token(CTX);
            while (tok == ' ') tok = FDL_token(CTX);
            if (isoneof(retsyms,tok)) return(0);
            }
        err = FDL_expr(CTX,&(*E)->next,retsyms);
        if (err) return(FDL_error(CTX,err,__FILE__,__LINE__));
        else return(0);
        
    case TYPE:
    case INDEX:
        
        err = FDL_ident(CTX);
        if (err) return(FDL_error(CTX,err,__FILE__,__LINE__));
        p = FDL_newnode(CTX,command);
        FDL_copyident(CTX,&p->ident);
        p->size = FDL_sizeof(p->ident);
        *E = p;
        tok = CTX->curr_token;
        while (tok == ' ') tok = FDL_token(CTX);
        if (isoneof(retsyms,tok)) return(0);
        err = FDL_expr(CTX,&(*E)->next,retsyms);
        if (err) return(FDL_error(CTX,err,__FILE__,__LINE__));
        else return(0);
        
    case FDL_SIZE:
        
        tok = CTX->curr_token;
        while (tok == ' ') tok = FDL_token(CTX);
        if (isdigit(tok)) {
            err = FDL_numint(CTX);
            if (err) return(FDL_error(CTX,err,__FILE__,__LINE__));
            *E = FDL_newnode(CTX,FIX_SIZE);
            (*E)->numint = CTX->curr_numint;
            }
        else if (isalpha(tok)) {
            err = FDL_ident(CTX);
            if (err) return(FDL_error(CTX,err,__FILE__,__LINE__));
            *E = FDL_newnode(CTX,FLT_SIZE);
            FDL_copyident(CTX,&(*E)->ident);
            }
        tok = CTX->curr_token;
        while (tok == ' ') tok = FDL_token(CTX);
        if (isoneof(retsyms,tok)) return(0);
        err = FDL_expr(CTX,&(*E)->next,retsyms);
        if (err) return(FDL_error(CTX,err,__FILE__,__LINE__));
        else return(0);
        }
	return -1; /* should really be an error check here */
    }


void FDL_showtree(FDL_CTX *CTX,FDL_NODE *p,int lvl)
{
    int i;
    
    while (p) {
        for (i=0;i<lvl;++i) printf("  ");
        if (p->type >= NUM_COMMANDS) {
            if (p->type == FIX_SIZE) printf("FIX_SIZE:%d\n",p->numint);
            else if (p->type == FLT_SIZE) {
                printf("FLT_SIZE:");
                if (p->ident) printf("%s %d\n",p->ident,p->size);
                else printf("%d\n",p->size);
                }
            else printf("UNKNOWN:\n");
            }
        else {
            printf("%s:",CTX->command[p->type]);
            if (p->ident) printf("%s %d\n",p->ident,p->size);
            else printf("%d\n",p->size);
            }
        FDL_showtree(CTX,p->child,lvl+1);
        p = p->next;
        }
    }


void FDL_polish(FDL_CTX *CTX,FDL_NODE *p,FDL_NODE *parent,int lvl)
{
    FDL_NODE *q,*prev;
    ENTRY *e,*el;
    
    if (!p) return;
    q = p->child;
    prev = NULL;
    while (q) {
        FDL_polish(CTX,q,p,lvl+1);
        q->prev = prev;
        prev = q;
        q = q->next;
        }
    p->parent = parent;
    p->lvl = lvl;
    if (p->type == ELEMENT || p->type == ARRAY || p->type == FDL) {
        /*
         ** Add an entry to the hash table.
         */
        e = (ENTRY *)malloc(sizeof(ENTRY));
        if (!e)
            FDL_error(CTX,FDL_MSG_MEMSYM,__FILE__,__LINE__);
        e->ident = p->ident;
        e->point = p;
        if (!HTBL_insert(CTX->T,(void *)e)) {
            /*
             ** Identifier already exists.
             ** Must have unique identifiers for nodes.
             */
            printf("The identifier '%s' was previously defined.\n",p->ident);
            printf("Please choose unique identifiers for the commands\n");
            printf("ELEMENT, ARRAY and FDL.\n");
            exit(1);
            }
        }
    else if (p->type == INDEX) {
        /*
         ** Add an entry to the hash table.
         */
        e = (ENTRY *)malloc(sizeof(ENTRY));
        if (!e) FDL_error(CTX,FDL_MSG_MEMBUILDSYM,__FILE__,__LINE__);
        e->ident = p->ident;
        e->point = p;
        /*
         ** Set initial index values to 0.
         */
        e->value = 0;
        e->curr_value = 0;
        HTBL_insert(CTX->T,(void *)e);
        }
    else if (p->type == FLT_SIZE) {
        /*
         ** Make sure that the size for a floating array is okay.
         */
        e = (ENTRY *)malloc(sizeof(ENTRY));
        if (!e) FDL_error(CTX,FDL_MSG_MEMBUILDSYM,__FILE__,__LINE__);
        e->ident = p->ident;
        el = (ENTRY *)HTBL_lookup(CTX->T,(void *)e);
        free(e);
        if (!el) { 
            /*
             ** Identifier does not exist.
             */
            printf("The identifier '%s' has not been defined.\n",
                   p->ident);
            printf("When using floating arrays, SIZE must be defined\n");
            printf("by an integer element in a fixed portion of the\n");
            printf("file.\n");
            exit(1);
            }
        q = el->point;
        if (q->type != ELEMENT) {
            printf("When using floating arrays, SIZE must be defined\n");
            printf("by an integer element in a fixed portion of the\n");
            printf("file.\n");
            exit(1);
            }
        else if (q->child->type != TYPE) {
            printf("When using floating arrays, SIZE must be defined\n");
            printf("by an integer element in a fixed portion of the\n");
            printf("file.\n");
            exit(1);
            }
        else if (strcmp(q->child->ident,"int")) {
            printf("When using floating arrays, SIZE must be defined\n");
            printf("by an integer element in a fixed portion of the\n");
            printf("file.\n");
            exit(1);
            }
        /*
         ** Invalidate the size by setting to negative value.
         */
        el->value = -1;
        }		    
    }


void FDL_init_htbl(FDL_CTX *CTX)
{
    int size;
    
    /* Impose at least 75% loading factor for the hash table */
    size = (int)(CTX->identifiers/0.75);
    CTX->T = HTBL_init(size,0,key_compare,hash_func);
    FDL_polish(CTX,CTX->ast,NULL,0);
    }


FDL_CTX *FDL_openwith(char *filename,char *fdlname)
{
    FILE *fp,*fdlp;
    FDL_CTX *CTX;
    int err;
    
    fdlp = fopen(fdlname,"r");
    if (!fdlp) {
	perror("FDL_openwith(fdlname)");
        printf("Could not open file %s.\n",fdlname);
        exit(1);
        }
    fp = fopen(filename,"r");
    if (!fp) {
	perror("FDL_openwith(filename)");
        printf("Could not open file %s.\n",filename);
        exit(1);
        }
    CTX = FDL_initialize(fp,fdlp,0);
    FDL_token(CTX);
    err = FDL_expr(CTX,&CTX->ast,"@");
    if (err) FDL_error(CTX,err,__FILE__,__LINE__);
    FDL_trailer(CTX);
    FDL_init_htbl(CTX);
    CTX->curr_point = CTX->ast;
	CTX->iOffset = 0;
    return(CTX);
    }


FDL_CTX *FDL_open(char *filename)
{
    FILE *fp;
    FDL_CTX *CTX;
    int err;
    
    fp = fopen(filename,"r");
    if (!fp) {
        printf("Could not open file %s.\n",filename);
        exit(1);
        }
    CTX = FDL_initialize(fp,NULL,0);
    FDL_token(CTX);
    err = FDL_expr(CTX,&CTX->ast,"@");
    if (err) FDL_error(CTX,err,__FILE__,__LINE__);
    FDL_trailer(CTX);
    FDL_init_htbl(CTX);
    CTX->curr_point = CTX->ast;
	CTX->iOffset = ftell(CTX->fileptr);
    return(CTX);
    }


FDL_CTX *FDL_createwith(char *filename,char *fdlname)
{
    FILE *fp,*fdlp;
    FDL_CTX *CTX;
    int err;
    
    fdlp = fopen(fdlname,"r");
    if (!fdlp) {
        printf("Could not open file %s.\n",fdlname);
        exit(1);
        }
    fp = fopen(filename,"w+");
    if (!fp) {
        printf("Could not open file %s.\n",filename);
        exit(1);
        }
    CTX = FDL_initialize(fp,fdlp,0);
    FDL_token(CTX);
    err = FDL_expr(CTX,&CTX->ast,"@");
    if (err) FDL_error(CTX,err,__FILE__,__LINE__);
    FDL_trailer(CTX);
    FDL_init_htbl(CTX);
    CTX->curr_point = CTX->ast;
	CTX->iOffset = 0;
    return(CTX);
    }


FDL_CTX *FDL_create(char *filename,char *fdlname)
{
    FILE *fp,*fdlp;
    FDL_CTX *CTX;
    int err;
    
    fdlp = fopen(fdlname,"r");
    if (!fdlp) {
        printf("Could not open file %s.\n",fdlname);
        exit(1);
        }
    fp = fopen(filename,"w+");
    if (!fp) {
        printf("Could not open file %s.\n",filename);
        exit(1);
        }
    CTX = FDL_initialize(fp,fdlp,1);
    FDL_token(CTX);
    err = FDL_expr(CTX,&CTX->ast,"@");
    if (err) FDL_error(CTX,err,__FILE__,__LINE__);
    FDL_trailer(CTX);
    FDL_init_htbl(CTX);
    CTX->curr_point = CTX->ast;
	CTX->iOffset = ftell(CTX->fileptr);
    return(CTX);
    }


FDL_CTX *FDL_modifywith(char *filename,char *fdlname)
{
    FILE *fp,*fdlp;
    FDL_CTX *CTX;
    int err;
    
    fdlp = fopen(fdlname,"r");
    if (!fdlp) {
        printf("Could not open file %s.\n",fdlname);
        exit(1);
        }
    fp = fopen(filename,"r+");
    if (!fp) {
        printf("Could not open file %s.\n",filename);
        exit(1);
        }
    CTX = FDL_initialize(fp,fdlp,0);
    FDL_token(CTX);
    err = FDL_expr(CTX,&CTX->ast,"@");
    if (err) FDL_error(CTX,err,__FILE__,__LINE__);
    FDL_trailer(CTX);
    FDL_init_htbl(CTX);
    CTX->curr_point = CTX->ast;
	CTX->iOffset = 0;
    return(CTX);
    }


FDL_CTX *FDL_modify(char *filename)
{
    FILE *fp;
    FDL_CTX *CTX;
    int err;
    
    fp = fopen(filename,"r+");
    if (!fp) {
        printf("Could not open file %s.\n",filename);
        exit(1);
        }
    CTX = FDL_initialize(fp,NULL,0);
    FDL_token(CTX);
    err = FDL_expr(CTX,&CTX->ast,"@");
    if (err) FDL_error(CTX,err,__FILE__,__LINE__);
    FDL_trailer(CTX);
    FDL_init_htbl(CTX);
    CTX->curr_point = CTX->ast;
	CTX->iOffset = ftell(CTX->fileptr);
    return(CTX);
    }


void FDL_fread(FDL_CTX *,void *);

void FDL_seek(FDL_CTX *,FDL_NODE *);

int FDL_elemsize(FDL_CTX *CTX,FDL_NODE *p)
{
    int size,mult;
    FDL_NODE *q,*sp,*old;
    ENTRY e,*el;
    
    if (p->size == -1) {
        /*
         ** Element size is invalidated.
         */
        if (p->type != ARRAY) {
            size = 0;
            q = p->child;
            while (q) {
                size += FDL_elemsize(CTX,q);
                q = q->next;
                }
            }
        else {
            size = 0;
            q = p->child;
            sp = NULL;
            while (q) {
                if (q->type == FLT_SIZE || q->type == FIX_SIZE) {
                    sp = q;
                    break;
                    }
                q = q->next;
                }
            assert(sp != NULL);
            q = sp;
            if (q->type == FIX_SIZE) {
                mult = q->numint;
                size = 0;
                q = p->child;
                while (q) {
                    size += FDL_elemsize(CTX,q);
                    q = q->next;
                    }
                size = size*mult;
                }
            else {
                /*
                 ** This is the tricky bit, floating array.
                 */
                e.ident = q->ident;
                el = (ENTRY *)HTBL_lookup(CTX->T,(void *)(&e));
                assert(el != NULL);
                mult = el->value;
                if (mult == -1) {
                    /*
                     ** Worst possible case.
                     ** But we must be able to find this value from the
                     ** CURRENT file position pointer.
                     ** Note: Mutual recursion!
                     */
					old = CTX->curr_point;
                    FDL_seek(CTX,el->point);
					FDL_fread(CTX,&mult);
					FDL_seek(CTX,old);
                    el->value = mult;
                    }
                size = 0;
                q = p->child;
                while (q) {
                    size += FDL_elemsize(CTX,q);
                    q = q->next;
                    }
                size = size*mult;
                }
            }
        p->size = size;
        }
    return(p->size);
    }


void FDL_seek(FDL_CTX *CTX,FDL_NODE *p)
{
    FDL_NODE *pp,*cc,*q;
    ENTRY e,*el=NULL;
    int clvl,plvl;
    int pdist,cdist;
    int dist;
    
    pdist = 0;
    cdist = 0;
    pp = p;
    cc = CTX->curr_point;
    plvl = pp->lvl;
    clvl = cc->lvl;
    while (plvl > clvl) {
        while (pp->prev) {
            pp = pp->prev;
            pdist = pdist + FDL_elemsize(CTX,pp);
            }
        pp = pp->parent;
        --plvl;
        if (pp->type == ARRAY) {
            dist = 0;
            q = pp->child;
            while (q) {
                if (q->type == INDEX) {
                    e.ident = q->ident;
                    el = (ENTRY *)HTBL_lookup(CTX->T,(void *)(&e));
                    assert(el != NULL);
                    }
                else dist += FDL_elemsize(CTX,q);
                q = q->next;
                }
            pdist += dist*el->value;
            }
        }
    if (cc->type == ARRAY) {
        dist = 0;
        q = cc->child;
        while (q) {
            if (q->type == INDEX) {
                e.ident = q->ident;
                el = (ENTRY *)HTBL_lookup(CTX->T,(void *)(&e));
                assert(el != NULL);
                }
            else dist += FDL_elemsize(CTX,q);
            q = q->next;
            }
        cdist += dist*el->curr_value;
        }
    while (clvl > plvl) {
        while (cc->prev) {
            cc = cc->prev;
            cdist = cdist + FDL_elemsize(CTX,cc);
            }
        cc = cc->parent;
        --clvl;
        if (cc->type == ARRAY) {
            dist = 0;
            q = cc->child;
            while (q) {
                if (q->type == INDEX) {
                    e.ident = q->ident;
                    el = (ENTRY *)HTBL_lookup(CTX->T,(void *)(&e));
                    assert(el != NULL);
                    }
                else dist += FDL_elemsize(CTX,q);
                q = q->next;
                }
            cdist += dist*el->curr_value;
            }
        }
    while (cc != pp) {
        while (pp->prev) {
            pp = pp->prev;
            pdist = pdist + FDL_elemsize(CTX,pp);
            }
        pp = pp->parent;
        --plvl;
        while (cc->prev) {
            cc = cc->prev;
            cdist = cdist + FDL_elemsize(CTX,cc);
            }
        cc = cc->parent;
        --clvl;
        if (pp->type == ARRAY) {
            dist = 0;
            q = pp->child;
            while (q) {
                if (q->type == INDEX) {
                    e.ident = q->ident;
                    el = (ENTRY *)HTBL_lookup(CTX->T,(void *)(&e));
                    assert(el != NULL);
                    }
                else dist += FDL_elemsize(CTX,q);
                q = q->next;
                }
            pdist += dist*el->value;
            }
        if (cc->type == ARRAY) {
            dist = 0;
            q = cc->child;
            while (q) {
                if (q->type == INDEX) {
                    e.ident = q->ident;
                    el = (ENTRY *)HTBL_lookup(CTX->T,(void *)(&e));
                    assert(el != NULL);
                    }
                else dist += FDL_elemsize(CTX,q);
                q = q->next;
                }
            cdist += dist*el->curr_value;
            }
        }
    fseek(CTX->fileptr,pdist-cdist,SEEK_CUR);
    /*
     ** Now update the current values.
     */
    CTX->curr_point = p;
    while (p->parent) {
		p = p->parent;
        if (p->type == ARRAY) {
            q = p->child;
            while (q) {
                if (q->type == INDEX) {
                    e.ident = q->ident;
                    el = (ENTRY *)HTBL_lookup(CTX->T,(void *)(&e));
                    assert(el != NULL);
                    }
                q = q->next;
                }
            el->curr_value = el->value;
            }
        }
#ifdef DEBUG
    printf("seek:%d\n",pdist-cdist);
#endif
    }


int FDL_offsetp(FDL_CTX *CTX,FDL_NODE *p)
{
    FDL_NODE *pp,*cc,*q;
    ENTRY e,*el=NULL;
    int clvl,plvl;
    int pdist,cdist;
    int dist;
    
    pdist = 0;
    cdist = 0;
    pp = p;
    cc = CTX->ast;
    plvl = pp->lvl;
    clvl = cc->lvl;
    while (plvl > clvl) {
        while (pp->prev) {
            pp = pp->prev;
            pdist = pdist + FDL_elemsize(CTX,pp);
            }
        pp = pp->parent;
        --plvl;
        if (pp->type == ARRAY) {
            dist = 0;
            q = pp->child;
            while (q) {
                if (q->type == INDEX) {
                    e.ident = q->ident;
                    el = (ENTRY *)HTBL_lookup(CTX->T,(void *)(&e));
                    assert(el != NULL);
                    }
                else dist += FDL_elemsize(CTX,q);
                q = q->next;
                }
            pdist += dist*el->value;
            }
        }
    if (cc->type == ARRAY) {
        dist = 0;
        q = cc->child;
        while (q) {
            if (q->type == INDEX) {
                e.ident = q->ident;
                el = (ENTRY *)HTBL_lookup(CTX->T,(void *)(&e));
                assert(el != NULL);
                }
            else dist += FDL_elemsize(CTX,q);
            q = q->next;
            }
        cdist += dist*el->curr_value;
        }
    while (clvl > plvl) {
        while (cc->prev) {
            cc = cc->prev;
            cdist = cdist + FDL_elemsize(CTX,cc);
            }
        cc = cc->parent;
        --clvl;
        if (cc->type == ARRAY) {
            dist = 0;
            q = cc->child;
            while (q) {
                if (q->type == INDEX) {
                    e.ident = q->ident;
                    el = (ENTRY *)HTBL_lookup(CTX->T,(void *)(&e));
                    assert(el != NULL);
                    }
                else dist += FDL_elemsize(CTX,q);
                q = q->next;
                }
            cdist += dist*el->curr_value;
            }
        }
    while (cc != pp) {
        while (pp->prev) {
            pp = pp->prev;
            pdist = pdist + FDL_elemsize(CTX,pp);
            }
        pp = pp->parent;
        --plvl;
        while (cc->prev) {
            cc = cc->prev;
            cdist = cdist + FDL_elemsize(CTX,cc);
            }
        cc = cc->parent;
        --clvl;
        if (pp->type == ARRAY) {
            dist = 0;
            q = pp->child;
            while (q) {
                if (q->type == INDEX) {
                    e.ident = q->ident;
                    el = (ENTRY *)HTBL_lookup(CTX->T,(void *)(&e));
                    assert(el != NULL);
                    }
                else dist += FDL_elemsize(CTX,q);
                q = q->next;
                }
            pdist += dist*el->value;
            }
        if (cc->type == ARRAY) {
            dist = 0;
            q = cc->child;
            while (q) {
                if (q->type == INDEX) {
                    e.ident = q->ident;
                    el = (ENTRY *)HTBL_lookup(CTX->T,(void *)(&e));
                    assert(el != NULL);
                    }
                else dist += FDL_elemsize(CTX,q);
                q = q->next;
                }
            cdist += dist*el->curr_value;
            }
        }
	return(pdist-cdist);
    }


void FDL_next(FDL_CTX *CTX)
{
    FDL_NODE *size,*index,*q,*old,*p;
    ENTRY e,*el;
    int max;
        
    p = CTX->curr_point;
    if (p->type == ARRAY) {
        q = p->child;
        index = NULL;
        size = NULL;
        while (q) {
            if (q->type == FLT_SIZE || q->type == FIX_SIZE) size = q;
            else if (q->type == INDEX) index = q;
            if (index && size) break;
            q = q->next;
            }
        assert(index != NULL);
        assert(size != NULL);
        if (size->type == FIX_SIZE) max = q->numint;
        else {
            e.ident = size->ident;
            el = (ENTRY *)HTBL_lookup(CTX->T,(void *)(&e));
            assert(el != NULL);
            max = el->value;
            if (max == -1) {
                /*
                 ** We must be able to find this value from the
                 ** CURRENT file position pointer.
                 ** Note: Mutual recursion!
                 */
				old = CTX->curr_point;
                FDL_seek(CTX,el->point);
				FDL_fread(CTX,&max);
				FDL_seek(CTX,old);
                el->value = max;
                }
            }
        e.ident = index->ident;
        el = (ENTRY *)HTBL_lookup(CTX->T,(void *)(&e));
        assert(el != NULL);
        if (el->value < (max-1)) {
            el->curr_value = el->value + 1;
            CTX->curr_point = p;
			return;
            }
        }
    if (p->next) {
		p = p->next;
		if (p->type == ARRAY) {
			q = p->child;
			index = NULL;
			while (q) {
				if (q->type == INDEX) {
					index = q;
					break;
					}
				q = q->next;
				}
			assert(index != NULL);
			e.ident = index->ident;
			el = (ENTRY *)HTBL_lookup(CTX->T,(void *)(&e));
			assert(el != NULL);
			el->curr_value = 0;
			}
		CTX->curr_point = p;
		return;
		}
    else {
        p = p->parent;
		CTX->curr_point = p;
        FDL_next(CTX);
		return;
        }
    }


void FDL_fread(FDL_CTX *CTX,void *value)
{
    fread(value,FDL_elemsize(CTX,CTX->curr_point),1,CTX->fileptr);
    if (ferror(CTX->fileptr))
		FDL_error(CTX,FDL_MSG_FILEREAD,__FILE__,__LINE__);
    FDL_next(CTX);
    }


void FDL_fwrite(FDL_CTX *CTX,void *value)
{
    fwrite(value,FDL_elemsize(CTX,CTX->curr_point),1,CTX->fileptr);
    if (ferror(CTX->fileptr))
		FDL_error(CTX,FDL_MSG_FILEWRITE,__FILE__,__LINE__);
    FDL_next(CTX);
    }


int FDL_read(FDL_CTX *CTX,char *ident,void *value)
{
    ENTRY e,*el;
    FDL_NODE *p;
    
    e.ident = ident;
    el = (ENTRY *)HTBL_lookup(CTX->T,(void *)(&e));
    if (!el) {
        /*
         ** This is not an identifier in the table, i.e.
         ** this file does not contain this symbol.
         */
		printf("Identifier:%s not found in symbol table.\n",ident);
        return(FDL_error(CTX,FDL_MSG_NOTFOUND,__FILE__,__LINE__));
        }
    p = el->point;
    /*
     ** Make sure it is a base type
     */
    if (p->child->type != TYPE)
        return(FDL_error(CTX,FDL_MSG_NOTTYPE,__FILE__,__LINE__));
    FDL_seek(CTX,p);
    FDL_fread(CTX,value);
    return(0);
    }


/*
 ** Note: All sizes of FLT_ARRAYs must be written to the file prior
 ** to writing any other data.
 ** Eg, for the rv.fld structure the variable "number_of_particles" must
 ** be written to first as it is used to determine the sizes of the
 ** mass_array and the rv_array.
 */
int FDL_write(FDL_CTX *CTX,char *ident,void *value)
{
    ENTRY e,*el;
    FDL_NODE *p;
    
    e.ident = ident;
    el = (ENTRY *)HTBL_lookup(CTX->T,(void *)(&e));
    if (!el) {
        /*
         ** This is not an identifier in the table, i.e.
         ** this file does not contain this symbol.
         */
		printf("Identifier:%s not found in symbol table.\n",ident);
        return(FDL_error(CTX,FDL_MSG_NOTFOUND,__FILE__,__LINE__));
        }
    p = el->point;
    /*
     ** Make sure it is a base type
     */
    if (p->child->type != TYPE)
        return(FDL_error(CTX,FDL_MSG_NOTTYPE,__FILE__,__LINE__));
    FDL_seek(CTX,p);
    FDL_fwrite(CTX,value);
    return(0);
    }


int FDL_index(FDL_CTX *CTX,char *ident,int index)
{
    FDL_NODE *p,*old;
    ENTRY e,*el,e1,*el1;
    int max;
    
    e.ident = ident;
    el = (ENTRY *)HTBL_lookup(CTX->T,(void *)(&e));
    if (!el) {
        /*
         ** This is not an identifier in the table, i.e.
         ** this file does not contain this symbol.
         */
		printf("Identifier:%s not found in symbol table.\n",ident);
        return(FDL_error(CTX,FDL_MSG_NOTFOUND,__FILE__,__LINE__));
        }
    /*
     ** Now make sure the index is not greater than the size.
     */
    p = el->point->parent->child;
    while (p) {
        if (p->type == FIX_SIZE || p->type == FLT_SIZE) break;
        p = p->next;
        }
    assert(p != NULL);
    if (p->type == FIX_SIZE) {
        if (index >= p->numint)
            FDL_error(CTX,FDL_MSG_BIGINDEX,__FILE__,__LINE__);
        }
    else {
        /*
         ** This is the tricky bit, floating array.
         */
        e1.ident = p->ident;
        el1 = (ENTRY *)HTBL_lookup(CTX->T,(void *)(&e1));
        assert(el1 != NULL);
        max = el1->value;
        if (max == -1) {
            /*
             ** Worst possible case.
             ** But we must be able to find this value from the
             ** CURRENT file position pointer.
             ** Note: Mutual recursion!
             */
			old = CTX->curr_point;
            FDL_seek(CTX,el1->point);
			FDL_fread(CTX,&max);
			FDL_seek(CTX,old);
            el1->value = max;
            }
        if (index >= max)
            FDL_error(CTX,FDL_MSG_BIGINDEX,__FILE__,__LINE__);
        }
    el->value = index;
    return(0);
    }


int FDL_offset(FDL_CTX *CTX,char *ident)
{
    ENTRY e,*el;
    FDL_NODE *p;
    
    e.ident = ident;
    el = (ENTRY *)HTBL_lookup(CTX->T,(void *)(&e));
    if (!el) {
        /*
         ** This is not an identifier in the table, i.e.
         ** this file does not contain this symbol.
         */
		printf("Identifier:%s not found in symbol table.\n",ident);
        return(FDL_error(CTX,FDL_MSG_NOTFOUND,__FILE__,__LINE__));
        }
    p = el->point;
	return(FDL_offsetp(CTX,p) + CTX->iOffset);
	}
