#ifndef FDL_HINCLUDED

#include <stdio.h>

#include "htable.h"

#define ELEMENT			0
#define ARRAY			1
#define FDL				2
#define TYPE			3
#define FDL_SIZE		4
#define INDEX			5
#define NUM_COMMANDS	6
#define FIX_SIZE		7
#define FLT_SIZE		8

#define MAX_IDENT		60
#define MAX_IDS			1000


#define FDL_MSG_SUCCESS		0
#define FDL_MSG_COMMAND		1
#define FDL_MSG_LONGIDENT	2
#define FDL_MSG_IDENT		3
#define FDL_MSG_UNEXCH		4
#define FDL_MSG_NUMINT		5
#define FDL_MSG_NOTTYPE		6
#define FDL_MSG_NOTFOUND	7
#define FDL_MSG_NOINDEX		8
#define FDL_MSG_NOSIZE		9
#define FDL_MSG_MEMCTX		10
#define FDL_MSG_MEMNODE		11
#define FDL_MSG_FILEREAD	12
#define FDL_MSG_MEMIDENT	13
#define FDL_MSG_MEMSYM		14
#define FDL_MSG_MEMBUILDSYM	15
#define FDL_MSG_BIGINDEX	16
#define FDL_MSG_FILEWRITE	17

#define FDL_MESSAGES		18

#define FDL_MESSAGE_INITIALIZE	{\
"E Undefined Error.",\
"E Unrecognized command.",\
"E Identifier is too long.",\
"E Error scanning identifier.",\
"E Unexpected character encountered while scanning.",\
"E Error scanning number.",\
"E The value requested is not a base type.",\
"E The identifier requested was not defined.",\
"E No INDEX field specified in ARRAY construct.",\
"E No SIZE field specified in ARRAY construct.",\
"E Not enough memory to allocate FDL context.",\
"E Not enough memory to allocate FDL node.",\
"E Error reading from file.",\
"E Not enough memory to allocate identifiers.",\
"E Not enough memory to allocate symbol table.",\
"E Not enough memory to build symbol table.",\
"E The index specified exceeds the bounds of the array.",\
"E Error writing to file."\
}

typedef struct fdl_node {
    struct fdl_node *child;
    struct fdl_node *parent;
    struct fdl_node *next;
    struct fdl_node *prev;
    int type;
    int lvl;
    int size;
    char *ident;
    int numint;
    } FDL_NODE;

typedef struct fdl_context {
    FILE *fileptr;
    FILE *fdlp;
    FDL_NODE *ast;
    FDL_NODE *trailer;
    HTBL *T;
    FDL_NODE *curr_point;
    int identifiers;
    int curr_token;
    char curr_ident[MAX_IDENT+1];
    int curr_numint;
    char *command[NUM_COMMANDS];
	int bCopyFdl;
	int iOffset;
    } FDL_CTX;

void FDL_finish(FDL_CTX *);
void FDL_showtree(FDL_CTX *,FDL_NODE *,int);
FDL_CTX *FDL_openwith(char *,char *);
FDL_CTX *FDL_open(char *);
FDL_CTX *FDL_createwith(char *,char *);
FDL_CTX *FDL_create(char *,char *);
FDL_CTX *FDL_modifywith(char *,char *);
FDL_CTX *FDL_modify(char *);
int FDL_read(FDL_CTX *,char *,void *);
int FDL_write(FDL_CTX *,char *,void *);
int FDL_index(FDL_CTX *,char *,int);
int FDL_offset(FDL_CTX *,char *);

#define FDL_HINCLUDED
#endif
