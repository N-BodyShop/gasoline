typedef struct inMSIMFSec
{
    struct MillerScaloContext ms;
    struct snContext sn;
} *MSSN;

double dMSIMFSec(MSSN p, double mass);
double dMSIMFSecM(MSSN p, double mass);
double dNSNIa(MSSN p, double mass1, double mass2);
double dMSNIa(MSSN p, double mass1, double mass2);

