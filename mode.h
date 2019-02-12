#ifndef MODE_HINCLUDED
#define MODE_HINCLUDED

#ifdef GALAXY
#define GASOLINE              //Do hydro 
#define NSMOOTHINNER          //Ensure at least 18 neighbours inside 1.41 h
#define DIVVOFF               //Remove the convergent flow requirement for starformation
#define JEANSSOFT             //Calculate presure floor to keep jeans length larger than softening or sph h (default is just sph h)
#define PARTICLESPLIT         //Particles will now divide if they get too heavy
#define SUPERBUBBLE	          //Use the Keller+ 2014 Superbubble feedback 
#define STARFORM              //Make new stars according to the starformation recipe
#define CHABRIER              //Use the Chabrier 2003 IMF for starformation (See DOI:10.1086/376392)  (SHOULD BE DEFAULT)
#define DIFFUSION             //Enable Thermal Diffusion
#define GDFORCE               //Use the geometric density force expression
#endif 

#ifdef JUSTHYDRO
#define GASOLINE              //Do hydro 
#define DIFFUSION             //Enable Thermal Diffusion
#define NSMOOTHINNER          //Ensure at least 18 neighbours inside 1.41 h
#define GDFORCE               //Use the geometric density force expression
#endif

#ifdef COLLISIONAL
#define COLLISIONS            //Use solid-body collisions (not compatible with -DGASOLINE)
#endif 

#ifdef GLASS_IC
#define GASOLINE              //Do hydro 
#define GLASS                 //Use this to damp for glass initial conditions
#endif 

#ifdef NBODY
#undef GASOLINE               //Disable SPH
#define NOCOOLING             //Disable Cooling
#undef COOLING_METAL          
#undef COOLING_COSMO
#undef COOLING_POLY
#undef COOLING_GRACKLE
#undef COOLING_MOLECULARH
#undef COOLING_PLANET
#endif

#ifdef BLASTWAVE
#undef SUPERBUBBLE
#endif

#ifdef KROUPA
#undef CHABRIER
#endif

#ifdef KROUPA01
#undef CHABRIER
#endif

#ifdef MILLERSCALO
#undef CHABRIER
#endif

#ifdef NONTHERMAL
#undef SUPERBUBBLE
#define UNONCOOL
#endif

#ifdef JEANSSMOOTH
#undef JEANSSOFT
#undef JEANSSOFTONLY
#endif

#endif
