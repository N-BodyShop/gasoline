#ifndef LOG_HINCLUDED
#define LOG_HINCLUDED
#endif

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

// How many characters wide should our logs be?
#define LOGCOL 160

typedef struct LogFormatter {
    int labelCnt;
    int *lineMem;
    char **label;
    char **line;
} LOGGER;

void LogParams(LOGGER *lgr, char *label, char *name, ...);
void LogFlush(LOGGER *lgr, FILE *fp);
LOGGER *initLog();
