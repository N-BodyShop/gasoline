#include "log.h"
#include "assert.h"
LOGGER *initLog()
{
    LOGGER *lgr = malloc(sizeof(LOGGER));
    lgr->labelCnt = 0;
    return lgr;
}
void LogParams(LOGGER * lgr, char *label, char *param)
{
    assert(strlen(label) > 0);
    assert(strlen(param) > 0);
    assert(strlen(param) < LOGCOL);
    assert(strlen(param) < LOGCOL);
    int i;
    int found=0;
    int curlabel=0;
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
        lgr->label = memcpy(malloc(lgr->labelCnt*sizeof(char*)), lgr->label, curlabel*sizeof(char*));
        lgr->line = memcpy(malloc(lgr->labelCnt*sizeof(char*)), lgr->line, curlabel*sizeof(char*));
        lgr->label[curlabel] = malloc(LOGCOL*sizeof(char));
        lgr->line[curlabel] = malloc(LOGCOL*sizeof(char));
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
    if(((strlen(lgr->line[curlabel]) % LOGCOL) + strlen(param)) < LOGCOL-1)
    {
        strcat(lgr->line[curlabel], " ");
        strcat(lgr->line[curlabel], param);
    }
    //If the line is over the limit, write the contents out and start a new line
    //with the requested label.
    else
    {
        lgr->line[curlabel] = memcpy(malloc((strlen(lgr->line[curlabel])+LOGCOL)*sizeof(char)), lgr->line[curlabel], 1+strlen(lgr->line[curlabel])*sizeof(char));
        strcat(lgr->line[curlabel], "\n# ");
        strcat(lgr->line[curlabel], lgr->label[curlabel]);
        strcat(lgr->line[curlabel], ": ");
        strcat(lgr->line[curlabel], param);
    }
}

void LogFlush(LOGGER * lgr, FILE *fp)
{
    int i;
    for(i=0;i<lgr->labelCnt;i++)
    {
        fprintf(fp, lgr->line[i]);
        /*free(lgr->line[i]);*/
        /*free(lgr->label[i]);*/
    }
    /*free(lgr);*/
}
