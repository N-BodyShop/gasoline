#ifndef COOLING_HINCLUDED
#define COOLING_HINCLUDED

#ifdef GASOLINE
#ifndef NOCOOLING

#ifdef COOLING_PLANET
#include "cooling_planet.h"
#else

#ifdef COOLING_COSMO
#include "cooling_cosmo.h"
#else

#error "No valid cooling function specified"

#endif
#endif

#endif
#endif

#endif
