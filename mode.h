#ifndef MODE_HINCLUDED
#define MODE_HINCLUDED

#ifdef GALAXY
GASOLINE              //Do hydro 
JEANSSOFT             //Calculate presure floor to keep jeans length larger than softening or sph h (default is just sph h)
PARTICLESPLIT         //Particles will now divide if they get too heavy
SUPERBUBBLE	          //Use the Keller+ 2014 Superbubble feedback 
STARFORM              //Make new stars according to the starformation recipe
CHABRIER              //Use the Chabrier 2003 IMF for starformation (See DOI:10.1086/376392)  (SHOULD BE DEFAULT)
#endif 

#ifdef COLLISIONS
#endif

#ifdef NBODY
#endif

#endif
