#include "log.h"
#include "assert.h"
LOGGER *initLog()
{
    LOGGER *lgr = malloc(sizeof(LOGGER));
    lgr->labelCnt = 0;
    return lgr;
}
void LogParams(LOGGER * lgr, char *label, char *name,  ...)
{
    assert(strlen(label) > 0);
    assert(strlen(name) > 0);
    assert(strlen(name) < LOGCOL);
    assert(strlen(name) < LOGCOL);
    va_list argp;
    va_start(argp, name);
    int i,found=0,curlabel=0;
    void *c;
    char param[160];
    vsprintf(param,name,argp);
    //Find where the label is in our label array.
    for(i=0;i<lgr->labelCnt;i++)
    {
        if(strcmp(lgr->label[i], label) == 0)
        {
            found=1;
            curlabel=i;
            break;
        }
    }
    //New label, add it to the label list
    if(!found)
    {
        curlabel=lgr->labelCnt;
        lgr->labelCnt++;
        if(lgr->labelCnt == 1)
        {
            lgr->label = malloc(sizeof(char*));
            lgr->line = malloc(sizeof(char*));
            lgr->lineMem = malloc(sizeof(int));
        }
        else
        {
            lgr->label = realloc(lgr->label, lgr->labelCnt*sizeof(char*));
            lgr->line = realloc(lgr->line, lgr->labelCnt*sizeof(char*));
            lgr->lineMem = realloc(lgr->lineMem, lgr->labelCnt*sizeof(int));
        }
        lgr->label[curlabel] = malloc(LOGCOL+1);
        lgr->line[curlabel] = malloc(LOGCOL+1);
        lgr->lineMem[curlabel] = LOGCOL;
        strcpy(lgr->label[curlabel], label);
        lgr->line[curlabel][0] = '\0';
    }
    //Check to see if line is blank, and if so, append the label to it.
    if(strlen(lgr->line[curlabel]) == 0)
    {
        strcat(lgr->line[curlabel], "\n# ");
        strcat(lgr->line[curlabel], lgr->label[curlabel]);
        strcat(lgr->line[curlabel], ": ");
    }
    //Check to make sure the line will be under the column limit
    if((strlen(lgr->line[curlabel]) + strlen(param)) < lgr->lineMem[curlabel])
    {
        strcat(lgr->line[curlabel], " ");
        strcat(lgr->line[curlabel], param);
    }
    //If the line is over the limit, write the contents out and start a new line
    //with the requested label.
    else
    {
        lgr->line[curlabel] = realloc(lgr->line[curlabel], lgr->lineMem[curlabel]+LOGCOL);
        lgr->lineMem[curlabel] += LOGCOL;
        strcat(lgr->line[curlabel], "\n# ");
        strcat(lgr->line[curlabel], lgr->label[curlabel]);
        strcat(lgr->line[curlabel], ": ");
        strcat(lgr->line[curlabel], param);
    }
}

void LogFlush(LOGGER * lgr, FILE *fp)
{
    int i;
    void *c;
    for(i=0;i<lgr->labelCnt;i++)
    {
        fprintf(fp, lgr->line[i]);
        c = lgr->line[i];
        free(c);
        c = lgr->label[i];
        free(c);
    }
    free(lgr->label);
    free(lgr->line);
    free(lgr->lineMem);
    free(lgr);
}
