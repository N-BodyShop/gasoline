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

/* Filename for list of particles with rejected initial positions */

#define REJECTS_FILE "rejects.out"

/* Filename for collision log */

#define COLLISION_LOG "ss.collisions"

#endif /* !SSDEFS_HINCLUDED */

/* ssdefs.h */
