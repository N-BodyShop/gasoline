#ifndef MODE_HINCLUDED
#define MODE_HINCLUDED

#ifdef GALAXY
#define GASOLINE              //Do hydro 
#define JEANSSOFT             //Calculate presure floor to keep jeans length larger than softening or sph h (default is just sph h)
#define PARTICLESPLIT         //Particles will now divide if they get too heavy
#define SUPERBUBBLE	          //Use the Keller+ 2014 Superbubble feedback 
#define STARFORM              //Make new stars according to the starformation recipe
#define CHABRIER              //Use the Chabrier 2003 IMF for starformation (See DOI:10.1086/376392)  (SHOULD BE DEFAULT)
#endif 

#ifdef COLLISIONS
#define AGGS                  //Include support for aggregates (like asteroids and rubble piles) See aggs.c
#define COLLISIONS            //Use solid-body collisions (not compatible with -DGASOLINE)
#endif 

#ifdef GLASS_IC
#define GASOLINE              //Do hydro 
#define GLASS                 //Use this to damp for glass initial conditions
#endif 

#ifdef TURBULENT
#define GASOLINE              //Do hydro 
#define OUTURBDRIVER          //Drive turbulence as in Price & Federrath 2010 with OU variables
#endif 

#ifdef NBODY
#undef GASOLINE               //Disable SPH
#endif


#endif
