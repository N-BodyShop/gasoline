#ifndef LOG_HINCLUDED
#define LOG_HINCLUDED
#endif

#include <stdio.h>
#include <string.h>

// How many characters wide should our logs be?
#define LOGCOL 160

typedef struct LogFormatter {
    int labelCnt;
    char **label;
    char **line;
} * LOGGER;

void LogParams(LOGGER lgr, char *label, char *param, FILE *fp);
void LogFlush(LOGGER lgr, FILE *fp);

