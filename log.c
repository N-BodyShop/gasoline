#include "log.h"

void LogParams(LOGGER lgr, char *label, char *param, FILE *fp)
{
    int i;
    int found=0;
    int curlabel=0;
    //Find where the label is in our label array.
    for(i=0;i<lgr->labelCnt;i++)
    {
        if(strcmp(lgr->label[i], label))
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
        lgr->label[lgr->labelCnt] = malloc(LOGCOL*sizeof(char));
        lgr->line[lgr->labelCnt] = malloc(LOGCOL*sizeof(char));
        strcpy(lgr->label[lgr->labelCnt], label);
        lgr->labelCnt++;
    }
    //Check to see if line is blank, and if so, append the label to it.
    if(strlen(lgr->line[curlabel]) == 0)
    {
        strcat(lgr->line[curlabel], "\n# ");
        strcat(lgr->line[curlabel], lgr->label[curlabel]);
        strcat(lgr->line[curlabel], ": ");
    }
    //Check to make sure the line will be under the column limit
    if(strlen(lgr->line[curlabel]) + strlen(label) < LOGCOL-1)
    {
        strcat(lgr->line[curlabel], " ");
        strcat(lgr->line[curlabel], label);
    }
    //If the line is over the limit, write the contents out and start a new line
    //with the requested label.
    else
    {
        fprintf(fp, lgr->line[curlabel]);
        lgr->line[curlabel][0] = '\0';
        strcat(lgr->line[curlabel], "\n# ");
        strcat(lgr->line[curlabel], lgr->label[curlabel]);
        strcat(lgr->line[curlabel], ": ");
        strcat(lgr->line[curlabel], label);
    }
}

void LogFlush(LOGGER lgr, FILE *fp)
{
    int i;
    for(i=0;i<lgr->labelCnt;i++)
    {
    if(strlen(lgr->line[i]) > strlen(lgr->label[i])+3)
    {
        fprintf(fp, lgr->line[i]);
        lgr->line[i][0] = '\0';
        lgr->label[i][0] = '\0';
    }
    }
}
