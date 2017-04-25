
#ifndef COOLING_HINCLUDED
#define COOLING_HINCLUDED

#ifdef GASOLINE
#ifndef NOCOOLING

#ifdef COOLING_DISK
#include "cooling_disk.h"
#else

#ifdef COOLING_PLANET
#include "cooling_planet.h"
#else

#ifdef COOLING_COSMO
#include "cooling_cosmo.h"
#else

#ifdef COOLING_METAL
#include "cooling_metal.h"
#else

#ifdef COOLING_MOLECULARH
#include "cooling_metal_H2.h"
#else

#ifdef COOLING_POLY
#include "cooling_poly.h"
#else

#ifdef COOLING_GRACKLE
#include "cooling_grackle.h"
#else

#error "No valid cooling function specified"

#endif
#endif
#endif
#endif
#endif
#endif

#endif
#endif

#endif
#endif
#endif

