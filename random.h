
#ifndef RANDOM_HINCLUDED
#define RANDOM_HINCLUDED

/* define your favorite random number generator here */
/* synopsis: rangen() returns random value between 0 and 1 */

#include <stdlib.h> /* prototypes and RAND_MAX for rand() & srand() */

double rangen(void);
void srangen(long int seed);

#define rangen() (rand()/((double) RAND_MAX))
#define srangen(s) srand(s)

#endif /* RANDOM_H */
