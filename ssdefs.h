#ifndef SSDEFS_HINCLUDED
#define SSDEFS_HINCLUDED

/*
 ** ssdefs.h -- DCR 98-09-16
 ** ========
 ** Header file for external interfacing with Solar System data and routines.
 */

#include "ssio.h"

/* Object color identifiers */

#define SUN				5	/* Yellow */
#define JUPITER			2	/* Red */
#define SATURN			11	/* Khaki */
#define URANUS			6	/* Magenta */
#define NEPTUNE			7	/* Cyan */
#define PLANETESIMAL	3	/* Green */
#define TEST			4	/* Blue (test particle) */

#define RESERVED_COLOR(c) ((c) == SUN || (c) == JUPITER)

/* Filename for list of particles with rejected initial positions */

#define REJECTS_FILE "rejects.out"

/* Filenames for collision logs */

#define COLL_LOG_TXT "ss.coll.txt"
#define COLL_LOG_BIN "ss.coll.bin"

#endif
