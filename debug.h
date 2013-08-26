#ifndef DEBUG_HINCLUDED
#define DEBUG_HINCLUDED
#include <stdio.h>
#include <signal.h>
#ifdef DBGPRINT
#define dbgprint(MSG, ...) fprintf(stderr, "DEBUG $s %s()%d: " MSG "\n", __FILE__, __FUNCTION__, __LINE__, ##__VA_ARGS__)
#else
#define dbgprint(MSG, ...)
#endif
void catch_FPE(int sig);
#endif
