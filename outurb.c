//#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "pkd.h"
#include "outurb.h"

// Dimensions -- hardwired to 3
#define NUMDIMS 3

//#include "allvars.h"
//#include "proto.h"

/*
//Ornstein-Uhlenbeck variables
static double StOUVar;
static double *StOUPhases;
static gsl_rng *StRng;

//forcing field in fourie space
static double *StAmpl;
static double *StAka; //phases (real part)
static double *StAkb; //phases (imag part)
static double *StMode;
static int StNModes;

static int StTPrev;
static double StSolWeightNorm;
*/

/*double StEnergyAcc;
  double StEnergyDeacc;

  double StLastStatTime;

  FILE *FdTurb;*/

void outurb_AddParams( OUTURBPARAM *outurbparam, PRM prm ) {
    outurbparam->StDecay = 1.;
    prmAddParam(prm,"StDecay",2,&outurbparam->StDecay,
                sizeof(double),"StDecay",
                "<OU Turb driving decay time> = 1.");
    outurbparam->StEnergy = 0.0002;
    prmAddParam(prm,"StEnergy",2,&outurbparam->StEnergy,
                sizeof(double),"StEnergy",
                "<OU Turb driving Energy> = 0.0002");
    outurbparam->StDtFreq = 0.005;
    prmAddParam(prm,"StDtFreq",2,&outurbparam->StDtFreq,
                sizeof(double),"StDtFreq",
                "<OU Turb driving mode update timeStep> = 0.005");
    outurbparam->StKmin = 6.27;
    prmAddParam(prm,"StKmin",2,&outurbparam->StKmin,
                sizeof(double),"StKmin",
                "<OU Turb driving k min> = 6.27");
    outurbparam->StKmax = 12.57;
    prmAddParam(prm,"StKmax",2,&outurbparam->StKmax,
                sizeof(double),"StKmax",
                "<OU Turb driving k max> = 12.57");
    outurbparam->StSolWeight = 1.;
    prmAddParam(prm,"StSolWeight",2,&outurbparam->StSolWeight,
                sizeof(double),"StSolWeight",
                "<OU Turb driving mode Solenoidal fraction> = 1");
    outurbparam->StAmplFac = 1.;
    prmAddParam(prm,"StAmplFac",2,&outurbparam->StAmplFac,
                sizeof(double),"StAmplFac",
                "<OU Turb driving amplitude factor> = 1.");
    outurbparam->StScaleHeight = 1e30;
    prmAddParam(prm,"StScaleHeight",2,&outurbparam->StScaleHeight,
                sizeof(double),"StScaleHeight",
                "<OU Turb driving Scale Height (Gaussian)> = 1e30");
    outurbparam->StStartTime = 2e30;
    prmAddParam(prm,"StStartTime",2,&outurbparam->StStartTime,
        sizeof(double),"StStartTime",
                "<OU Turb driving start time -- for setting correct sequence> = ic time");

    outurbparam->StSeed = 42;
    prmAddParam(prm,"StSeed",1,&outurbparam->StSeed,
                sizeof(int),"StSeed",
                "<OU Turb driving seed> = 42");
    outurbparam->StSpectForm = 0;
    prmAddParam(prm,"StSpectForm",1,&outurbparam->StSpectForm,
                sizeof(int),"StSpectForm",
                "<OU Turb driving Spectral Form> = 2");
    }

void outurbLogParams( OUTURBPARAM *outurbparam, LOGGER *lgr ) {
  LogParams(lgr,"TURBULENT DRIVER","\n# OUturb: StDecay: %g",outurbparam->StDecay);
  LogParams(lgr,"TURBULENT DRIVER"," StEnergy: %g",outurbparam->StEnergy);
  LogParams(lgr,"TURBULENT DRIVER"," StDtFreq: %g",outurbparam->StDtFreq);
  LogParams(lgr,"TURBULENT DRIVER"," StKmin: %g",outurbparam->StKmin);
  LogParams(lgr,"TURBULENT DRIVER"," StKmax: %g",outurbparam->StKmax);
  LogParams(lgr,"TURBULENT DRIVER"," StSolWeight: %g",outurbparam->StSolWeight);
  LogParams(lgr,"TURBULENT DRIVER"," StAmplFac: %g",outurbparam->StAmplFac);
  LogParams(lgr,"TURBULENT DRIVER"," StScaleHeight: %g",outurbparam->StScaleHeight);
  LogParams(lgr,"TURBULENT DRIVER"," StSeed: %d",outurbparam->StSeed);
  LogParams(lgr,"TURBULENT DRIVER"," StSpectForm: %d",outurbparam->StSpectForm);
}


void outurb_init(OUTURB *pouturb, OUTURBPARAM outurbparam, int idSelf, int bDetails, int bRestart,
    double BoxSize, double dTime)
    {
    int ikx, iky, ikz;
    double kx,ky,kz,k;
    double ampl;
    OUTURB outurb;
    int ikxmax,ikymax,ikzmax;
    double kc, amin;

    if (outurbparam.StStartTime > 1e30) {
        if (bRestart && idSelf == 0)
            fprintf(stderr,"OUturb: ERROR: Restarting without StStartTime specified:  Random sequence inconsistent\n");
        assert(!bRestart);
        if (idSelf == 0)
            fprintf(stderr,"OUturb: WARNING: StStartTime not specified:  Cannot restart this run.\n");
        }

    *pouturb = outurb = (OUTURB) malloc(sizeof(struct OUturb));
    assert(outurb!=NULL);

    if (bDetails==2) 
        outurb->idSelf = idSelf;
    else
        outurb->idSelf = idSelf+10000; // suppress prints

    outurb->StDecay = outurbparam.StDecay;
    outurb->StEnergy = outurbparam.StEnergy;
    outurb->StDtFreq = outurbparam.StDtFreq;
    outurb->StKmin = outurbparam.StKmin;
    outurb->StKmax = outurbparam.StKmax;
    outurb->StSolWeight = outurbparam.StSolWeight;
    outurb->StAmplFac = outurbparam.StAmplFac;
    outurb->StScaleHeight = outurbparam.StScaleHeight;
    outurb->StSeed = outurbparam.StSeed;
    outurb->StSpectForm = outurbparam.StSpectForm;

    ikxmax = BoxSize*outurb->StKmax/2./M_PI;
#if NUMDIMS > 1
    ikymax = BoxSize*outurb->StKmax/2./M_PI;
#if NUMDIMS > 2
    ikzmax = BoxSize*outurb->StKmax/2./M_PI;
#else
    ikzmax = 0;
#endif

#else
    ikymax = 0;
    ikzmax = 0;
#endif

    outurb->StNModes = 0;
    for(ikx = 0;ikx <= ikxmax; ikx++)
        {
        kx = 2.*M_PI*ikx/BoxSize;
        for(iky = 0;iky <= ikymax; iky++)
            {
            ky = 2.*M_PI*iky/BoxSize;
            for(ikz = 0;ikz <= ikzmax; ikz++)
                {
                kz = 2.*M_PI*ikz/BoxSize;
                k = sqrt(kx*kx+ky*ky+kz*kz);
                if(k>=outurb->StKmin && k<=outurb->StKmax)
                    {
#if NUMDIMS ==1
                    outurb->StNModes+=1;
#endif
#if NUMDIMS == 2

                    outurb->StNModes+=2;
#endif
#if NUMDIMS == 3

                    outurb->StNModes+=4;

#endif
                    }
                }
            }
        }

    if(outurb->idSelf == 0)
        {
        printf("OUturb: Using %d modes, %d %d %d\n",outurb->StNModes,ikxmax,ikymax,ikzmax);
        }

    outurb->StMode = (double*) malloc(outurb->StNModes * 3 * sizeof(double));
    outurb->StAka = (double*) malloc(outurb->StNModes * 3 * sizeof(double));
    outurb->StAkb = (double*) malloc(outurb->StNModes * 3 * sizeof(double));
    outurb->StAmpl = (double*) malloc(outurb->StNModes * sizeof(double));
    outurb->StOUPhases = (double*) malloc(outurb->StNModes * 6 * sizeof(double));

    outurb->StOUVar = sqrt(outurb->StEnergy/outurb->StDecay);
    kc = 0.5*(outurb->StKmin+outurb->StKmax);
    amin = 0.;

#if NUMDIMS == 3
    outurb->StSolWeightNorm = sqrt(3.0/3.0)*sqrt(3.0)*1.0/sqrt(1.0-2.0*outurb->StSolWeight+3.0*outurb->StSolWeight*outurb->StSolWeight);
#endif
#if NUMDIMS == 2
    outurb->StSolWeightNorm = sqrt(3.0/2.0)*sqrt(3.0)*1.0/sqrt(1.0-2.0*outurb->StSolWeight+2.0*outurb->StSolWeight*outurb->StSolWeight);
#endif
#if NUMDIMS == 1
    outurb->StSolWeightNorm = sqrt(3.0/1.0)*sqrt(3.0)*1.0/sqrt(1.0-2.0*outurb->StSolWeight+1.0*outurb->StSolWeight*outurb->StSolWeight);
#endif

    outurb->StNModes = 0;

    for(ikx = 0;ikx <= ikxmax; ikx++)
        {
        kx = 2.*M_PI*ikx/BoxSize;

        for(iky = 0;iky <= ikymax; iky++)
            {
            ky = 2.*M_PI*iky/BoxSize;
         
            for(ikz = 0;ikz <= ikzmax; ikz++)
                {
                kz = 2.*M_PI*ikz/BoxSize;

                k = sqrt(kx*kx+ky*ky+kz*kz);
                if(k>=outurb->StKmin && k<=outurb->StKmax)
                    {
                    ampl = 0.;
                    if(outurb->StSpectForm == 0)
                        {
                        ampl = 1.; 
                        }
                    else if(outurb->StSpectForm ==  1)
                        {
                        ampl = 4.0*(amin-1.0)/((outurb->StKmax-outurb->StKmin)*(outurb->StKmax-outurb->StKmin))*((k-kc)*(k-kc))+1.0;
                        }
                    else if(outurb->StSpectForm == 2)
                        {
                        ampl = pow(outurb->StKmin,5./3)/pow(k,5./3);
                        }
                    else if(outurb->StSpectForm == 3)
                        {
                        ampl = pow(outurb->StKmin,2.)/pow(k,2.);
                        }
                    else
                        {
//                        endrun(66611666);    
                        assert(0);
                        }


                    outurb->StAmpl[outurb->StNModes+0] = ampl;
                    outurb->StMode[3*outurb->StNModes+0] = kx;
                    outurb->StMode[3*outurb->StNModes+1] = ky;
                    outurb->StMode[3*outurb->StNModes+2] = kz;
                    if(outurb->idSelf==0)
                        {
                        printf("OUturb: Mode: %d, ikx=%d, iky=%d, ikz=%d, kx=%f, ky=%f, kz=%f, ampl=%f\n",outurb->StNModes,ikx,iky,ikz,kx,ky,kz,ampl);

                        }
                    outurb->StNModes++;


#if NUMDIMS > 1
                    outurb->StAmpl[outurb->StNModes+0] = ampl;
                    outurb->StMode[3*outurb->StNModes+0] = kx;
                    outurb->StMode[3*outurb->StNModes+1] = -ky;
                    outurb->StMode[3*outurb->StNModes+2] = kz;
                    if(outurb->idSelf==0)
                        {
                        printf("OUturb: Mode: %d, ikx=%d, iky=%d, ikz=%d, kx=%f, ky=%f, kz=%f, ampl=%f\n",outurb->StNModes,ikx,-iky,ikz,kx,-ky,kz,ampl);

                        }
                    outurb->StNModes++;

#if NUMDIMS > 2
                    outurb->StAmpl[outurb->StNModes+0] = ampl;
                    outurb->StMode[3*outurb->StNModes+0] = kx;
                    outurb->StMode[3*outurb->StNModes+1] = ky;
                    outurb->StMode[3*outurb->StNModes+2] = -kz;
                    if(outurb->idSelf==0)
                        {
                        printf("OUturb: Mode: %d, ikx=%d, iky=%d, ikz=%d, kx=%f, ky=%f, kz=%f, ampl=%f\n",outurb->StNModes,ikx,iky,-ikz,kx,ky,-kz,ampl);

                        }
                    outurb->StNModes++;

                    outurb->StAmpl[outurb->StNModes+0] = ampl;
                    outurb->StMode[3*outurb->StNModes+0] = kx;
                    outurb->StMode[3*outurb->StNModes+1] = -ky;
                    outurb->StMode[3*outurb->StNModes+2] = -kz;
                    if(outurb->idSelf==0)
                        {
                        printf("OUturb: Mode: %d, ikx=%d, iky=%d, ikz=%d, kx=%f, ky=%f, kz=%f, ampl=%f\n",outurb->StNModes,ikx,-iky,-ikz,kx,-ky,-kz,ampl);

                        }
                    outurb->StNModes++;
#endif
#endif
                    }
                }
            }
        }

//    outurb->StTPrev = -1;

    /*  outurb->StEnergyDeacc = 0;
        outurb->StEnergyAcc = 0;
        outurb->StLastouturb->StatTime = outurb->TimeBegin - outurb->TimeBetouturb->Statistics;*/

    outurb->StRng = gsl_rng_alloc(gsl_rng_ranlxd1);
    gsl_rng_set(outurb->StRng, outurb->StSeed);

    outurb_st_init_ouseq(outurb);
    outurb_st_calc_phases(outurb);

//    printf("OUturb %d: calling set_turb_ampl in init_turb\n",(outurb->idSelf%10000));
//    outurb_set_turb_ampl(outurb, dTime);
    //gsl_rng_free(outurb->StRng);

    outurb->dTimePrev = (outurbparam.StStartTime < 1e30 ? outurbparam.StStartTime : dTime);
    outurb_set_turb_ampl(outurb, dTime);
    }

void outurb_st_init_ouseq(OUTURB outurb)
    {
    int i;

    for(i = 0;i<6*outurb->StNModes;i++)
        {
        outurb->StOUPhases[i] = outurb_st_grn(outurb)*outurb->StOUVar;
        }
    }

void outurb_st_update_ouseq(OUTURB outurb)
    {
    int i;
    double damping = exp( -outurb->StDtFreq/outurb->StDecay);

    for(i = 0;i<6*outurb->StNModes;i++)
        {
        outurb->StOUPhases[i] = outurb->StOUPhases[i] * damping + outurb->StOUVar * sqrt(1.-damping*damping)*outurb_st_grn(outurb);
        }
    }

double outurb_st_grn(OUTURB outurb)
    {
    double r0 = gsl_rng_uniform(outurb->StRng);
    double r1 = gsl_rng_uniform(outurb->StRng);

    return sqrt(2. * log(1. / r0) ) * cos(2. * M_PI * r1);
    }



void outurb_st_calc_phases(OUTURB outurb)
    {
    int i,j;
    for(i = 0; i < outurb->StNModes;i++)
        {
        double ka = 0.;
        double kb = 0.;
        double kk = 0.;

        int dim = NUMDIMS;

        for(j = 0; j<dim;j++)
            {
            kk += outurb->StMode[3*i+j]*outurb->StMode[3*i+j];
            ka += outurb->StMode[3*i+j]*outurb->StOUPhases[6*i+2*j+1];
            kb += outurb->StMode[3*i+j]*outurb->StOUPhases[6*i+2*j+0];
            }
        for(j = 0; j<dim;j++)
            {
            double diva = outurb->StMode[3*i+j]*ka/kk;
            double divb = outurb->StMode[3*i+j]*kb/kk;
            double curla = outurb->StOUPhases[6*i+2*j+0] - divb;
            double curlb = outurb->StOUPhases[6*i+2*j+1] - diva;

            outurb->StAka[3*i+j] = outurb->StSolWeight*curla+(1.-outurb->StSolWeight)*divb;
            outurb->StAkb[3*i+j] = outurb->StSolWeight*curlb+(1.-outurb->StSolWeight)*diva;
            }
        }
    }

void outurb_set_turb_ampl(OUTURB outurb, double dTime)
    {
    int i;

//  double delta = (Ti_Current - outurb->StTPrev) * outurb->Timebase_interval;
    double delta = (dTime - outurb->dTimePrev);

    if(outurb->idSelf == 0) {
        printf("OUturb: Update? %f %f %f %f\n",delta,outurb->StDtFreq,dTime,outurb->dTimePrev);
        }

    if(delta >= outurb->StDtFreq*0.9999)
        {
// outurb->StTPrev = outurb->StTPrev + outurb->StDtFreq/outurb->Timebase_interval;
// JW: Potential bug from original this is not called often enough compared to StDtFreq, added this while
// This should allow continue from restarts by regenerating full sequence from the start time
        while (delta >= outurb->StDtFreq*0.9999) {
            outurb_st_update_ouseq(outurb);
            delta -= outurb->StDtFreq;
            outurb->dTimePrev += outurb->StDtFreq;
            if ((outurb->idSelf%10000)==0)
                printf("OUturb: %d Time %f, New phases:  %f %f %f %f %f %f\n",outurb->idSelf,outurb->dTimePrev,
                outurb->StOUPhases[0],outurb->StOUPhases[1],outurb->StOUPhases[2],
                outurb->StOUPhases[3],outurb->StOUPhases[4],outurb->StOUPhases[5]);
            }

        outurb_st_calc_phases(outurb);


        if(outurb->idSelf == 0) {
//    printf("Updated turbulent stirring field at time %f.\n", outurb->StTPrev * outurb->Timebase_interval);
            printf("OUturb: Updated turbulent stirring field at time %f (%f).\n", outurb->dTimePrev, dTime );
            }
        }
    }

void outurb_add_turb_accel(OUTURB outurb, double dTime, PARTICLE *p, int nParticle)
    {
    int i, j, m;
    double acc[3];// dt_kick;
    double volume = 0;
    double a2x=0,a2y=0,a2z=0;
    int na=0;

    if(outurb->idSelf == 0) {
        printf("OUturb: add_turb_accel %f %d.\n", dTime,outurb->StNModes );
        }

    outurb_set_turb_ampl(outurb, dTime);

    for(i = 0; i < nParticle; i++) { 
        if (TYPEQueryACTIVE(&(p[i]))) {
                {
                double fx = 0;
                double fy = 0;
                double fz = 0;

                for(m = 0;m<outurb->StNModes;m++) //calc force
                    {
                    double kxx = outurb->StMode[3*m+0]*p[i].r[0];
                    double kyy = outurb->StMode[3*m+1]*p[i].r[1];
                    double kzz = outurb->StMode[3*m+2]*p[i].r[2];
                    double kdotx = kxx+kyy+kzz;
                    double ampl = outurb->StAmpl[m];

                    double realt = cos(kdotx);
                    double imagt = sin(kdotx);

                    fx += ampl*(outurb->StAka[3*m+0]*realt - outurb->StAkb[3*m+0]*imagt);
                    fy += ampl*(outurb->StAka[3*m+1]*realt - outurb->StAkb[3*m+1]*imagt);
                    fz += ampl*(outurb->StAka[3*m+2]*realt - outurb->StAkb[3*m+2]*imagt);
//                    if (i==0) printf("%d %d ax+= %g %g %g %g %g\n",outurb->idSelf,m,ampl,outurb->StAka[3*m+0],realt,outurb->StAkb[3*m+0],imagt);
                    }


                fx *= 2.*outurb->StAmplFac*outurb->StSolWeightNorm;//*volume;
                fy *= 2.*outurb->StAmplFac*outurb->StSolWeightNorm;//*volume;
                fz *= 2.*outurb->StAmplFac*outurb->StSolWeightNorm;//*volume;
//                if (i==0) printf("%d %d ax*= %g %g\n",outurb->idSelf,m,outurb->StAmplFac*outurb->StSolWeightNorm);
//            if(P[i].Mass > 0.) // Gadgetism
                    {
                    double u = p[i].r[2]/outurb->StScaleHeight;
                    double vfac = exp(-0.5*u*u);
                    fx *= vfac;
                    fy *= vfac;
                    fz *= vfac;
                    //FIXME correct for mean force
                    na++;
                    acc[0] = fx;///P[i].Mass;
                    a2x += fx*fx;
                    acc[1] = fy;///P[i].Mass;
                    a2y += fy*fy;
#if NUMDIMS > 2
                    acc[2] = fz;///P[i].Mass;
                    a2z += fz*fz;
#else
                    acc[2] = 0;
#endif

                    for(j = 0; j < 3; j++)
                        {
                        p[i].a[j] += acc[j]; // Gasoline Particle
                        }
                    }
                }
            }
        
        }
    if(outurb->idSelf == 0)
        {
            printf("OUturb: Finished turbulent accel computation n=%d a_rms %f %f %f.\n",na,sqrt(a2x/(na+1e-10)),sqrt(a2y/(na+1e-10)),sqrt(a2z/(na+1e-10)));
        }
    
    }



/* These are not useful for Gasoline */
void outurb_driving_step_first_half(OUTURB outurb, double dTime, double dDelta, PARTICLE *p, int nParticle)
        {
//    CPU_outurb->Step[CPU_MISC] += measure_time();
        int i, j;
        int ti_step, tstart, tend;
        double dvel[3], dt_gravkick, atime;

        outurb_add_turb_accel(outurb, dTime, p, nParticle);

        assert(0);
#if (0)
        for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i]) // Gadget uses RLE , we don't
            {
            ti_step = P[i].TimeBin ? (((int) 1) << P[i].TimeBin) : 0;

            tstart = P[i].Ti_begstep;   /* beginning of step */
            tend = P[i].Ti_begstep + ti_step / 2;   /* midpoint of step */

            if(outurb->ComovingIntegrationOn)
                dt_gravkick = get_gravkick_factor(tstart, tend);
            else
                dt_gravkick = (tend - tstart) * outurb->Timebase_interval;

            if(P[i].Type == 0)
                {


                for(j = 0; j < 3; j++)
                    {
                    dvel[j] = SphP[i].TurbAccel[j] * dt_gravkick;
                    P[i].Vel[j] += dvel[j];
                    //P[i].dp[j] += P[i].Mass * dvel[j];
                    }


                }
            }
#endif

//    CPU_outurb->Step[CPU_DRIFT] += measure_time();
        }

/* These are not useful for Gasoline */
void outurb_driving_step_second_half(OUTURB outurb, double dTime, double dDelta, PARTICLE *p, int nParticle)
        {
//    CPU_outurb->Step[CPU_MISC] += measure_time();
        int i, j;
        int ti_step, tstart, tend;
        double dvel[3], dt_gravkick, atime;

        assert(0);
#if (0)
        if(outurb->ComovingIntegrationOn)
            atime = outurb->Time;
        else
            atime = 1.0;


        for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
            {
            ti_step = P[i].TimeBin ? (((int) 1) << P[i].TimeBin) : 0;

            tstart = P[i].Ti_begstep + ti_step / 2; /* midpoint of step */
            tend = P[i].Ti_begstep + ti_step;   /* end of step */

            if(outurb->ComovingIntegrationOn)
                dt_gravkick = get_gravkick_factor(tstart, tend);
            else
                dt_gravkick = (tend - tstart) * outurb->Timebase_interval;

            if(P[i].Type == 0)
                {

                for(j = 0; j < 3; j++)
                    {
                    dvel[j] = SphP[i].TurbAccel[j] * dt_gravkick;
                    P[i].Vel[j] += dvel[j];
                    //P[i].dp[j] += P[i].Mass * dvel[j];

                    SphP[i].VelPred[j] = P[i].Vel[j];
                    }


                }
            }
#endif

//  CPU_outurb->Step[CPU_DRIFT] += measure_time();
        }



