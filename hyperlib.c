#include <math.h>

double asinh(double a)
{
    return log(a+sqrt(a*a+1));
}

double acosh(double a)
{
    return log(a+sqrt(a*a-1));
}
