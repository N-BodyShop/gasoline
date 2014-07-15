/*
  With molecular H
  icc -g -I../mdl/null -I../pkdgravGIT -DNOCOMPTON -DCOOLING_MOLECULARH -DGASOLINE -DVERBOSE ../pkdgravGIT/cooling_metal_H2.c ../pkdgravGIT/stiff.c ../pkdgravGIT/param.c testcoolGIT.c -o testcoolGIT -lm 

  With old cooling
icc -g -I../mdl/null -I../pkdgravGIT -DNOCOMPTON -DCOOLING_METAL -DGASOLINE -DVERBOSE ../pkdgravGIT/cooling_metal.c ../pkdgravGIT/stiff.c ../pkdgravGIT/param.c testcoolGIT.c -o testcoolGIT -lm

  icc -g -I../mdl/null -I../pkdgravGIT -DCOOLING_METAL -DGASOLINE -DVERBOSE ../pkdgravGIT/cooling_metal.c ../pkdgravGIT/stiff.c ../pkdgravGIT/param.c testcoolGIT.c -o testcoolGIT -lm 

  Charlotte added a broken extension for the table (nH > 1e3)
  It doesn't work -- turn it off with -DNOCOOLTABLEFIT
  e */

#include <fenv.h>
#include <fpu_control.h>
#include <signal.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include "stiff.h"
#ifdef COOLING_COSMO
#include "cooling_cosmo.h"
#endif
#ifdef COOLING_METAL
#include "cooling_metal.h"
#endif
#ifdef COOLING_BATE
#include "cooling_bate.h"
#endif
#ifdef COOLING_MOLECULARH
#include "cooling_metal_H2.h"
#endif
#ifdef COOLING_METAL_NOH2
#include "cooling_metal_noH2.h"
#endif

//Temperature range
#define MINTEMP 20
#define MAXTEMP 1e9
#define NTEMPSTEPS 400

//Density range
#define MINRHO 2e-5
#define MAXRHO 1e5
#define NRHOSTEPS 450


//Metallicity
#define METAL ZSOLAR

//External PdV Work
#define EXTHEAT 0

//Redshift
#define ZRED 0
#define TIME 3.541544e-01

//Constants
#define MH_G 1.66e-24


//Integration Time (seconds)
#define DT 3e16 //~1 Gyr

int main() {
    int i;
    FILE *fpout;
    double E,E_copy,newT,tstep,rhostep, dens, dt;
    clock_t t;//For timing
    //Cooling-specific variables
    COOL *cl, *cl_copy;
    PERBARYON Y, Y_copy;
    COOLPARTICLE cp;
    COOLPARAM clParam;

    // Turn on floating point exceptions!
    feenableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);

    //Initialize cooling context objects.
    //Use more-or-less default values here.
    clParam.bIonNonEqm = 1;
    clParam.bUV = 1;
    clParam.bUVTableUsesTime = 0;
    clParam.bDoIonOutput = 1;
    clParam.dMassFracHelium = 0.236;
    clParam.dCoolingTmin= 10;
    clParam.dCoolingTmax = 1e9;
    clParam.nCoolingTable = 15001;
    clParam.bLowTCool = 0;
    clParam.bSelfShield = 0;
#ifndef COOLING_COSMO
    clParam.bMetal = 1;
    strcpy(clParam.CoolInFile, "../../cooltable_xdr");
#endif
#ifdef COOLING_METAL
    clParam.dPhotoelectricHeating=0;
    clParam.dPhotoelectricnMin=0.01;
    clParam.dzTimeClampUV =  -1;
#endif

    cl = CoolInit();
    CoolInitRatesTable( cl, clParam );  /* reading UV and metal table */ 
    clInitConstants(cl, 1,1,1,1,1,clParam); //Won't need these, but they need initializing
    CoolSetTime( cl, TIME, ZRED); 
    //Let the stiff integrator take smaller steps than usual
    STIFF *st = ((clDerivsData *) (cl->DerivsData))->IntegratorContext;
    st->dtmin = 1e-17;
    
    //Print header
    fpout=fopen("testcool.out","w");
    fprintf(fpout, "Temperature (K)\tDensity (H/cc)\tEquilibrium T\tWalltime (s)\n");

    for(tstep=log10(MINTEMP);tstep<log10(MAXTEMP);tstep+=(log10(MAXTEMP)-log10(MINTEMP))/NTEMPSTEPS)
    {
        for(rhostep=log10(MINRHO);rhostep<log10(MAXRHO);rhostep+=(log10(MAXRHO)-log10(MINRHO))/NRHOSTEPS)
        {
            printf("Calculating Equlibrium T for %5.4eK %5.4eH/cc\n", pow(10.0,tstep), pow(10.0,rhostep));
            //Calculate density in code units
            dens = pow(10.0,rhostep)*MH_G;
            
            //Calculate ionization states
            // Explicit values: Fairly ionized initial guess
            Y.e = 1.0000000000000001e-17; 
            Y.Total = 0.8299999999999973; 
            Y.HI = 0.1;
            Y.HII = 0.6;
            Y.HeI = 0.058999999999999997; 
            Y.HeII = 9.9999999999999998e-17; 
            Y.HeIII = -9.9999999999999998e-17; 

            //Calculate new energy, temp
            t=clock();
            E = clThermalEnergy(Y.Total, pow(10.0,tstep));
            clIntegrateEnergy(cl, &Y, &E, EXTHEAT, dens, METAL,0, DT);
            CoolPERBARYONtoPARTICLE(cl, &Y, &cp, METAL);
            newT = CoolEnergyToTemperature(cl, &cp, E,dens, METAL);
            t = clock()-t;
            
            //Print results
            fprintf(fpout, "%7.6e\t%7.6e\t%7.6e\t%4.3e\n", pow(10.0,tstep), pow(10.0,rhostep), newT, ((float)t)/CLOCKS_PER_SEC);
        }
    }
    return 0;
}
