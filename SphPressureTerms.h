
/*
   Macros for function:
       void SphPressureTermsSym(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
   in smoothfcn.c

   This code is used 3 times.  Constrained to use Macros (not functions) because
   Macro definitions change between uses as follows:
   1) with p and q active
     hash define PACTIVE(xxx) xxx
     hash define QACTIVE(xxx) xxx
   2) with just p active
     hash define PACTIVE(xxx) xxx
     hash define QACTIVE(xxx)
   3) with just q active
     hash define PACTIVE(xxx)
     hash define QACTIVE(xxx) xxx

   All Macros and Variables not defined here are defined in smoothfcn.c
 */

#ifdef DRHODT
    #define DRHODTACTIVE(xxx) xxx
    #ifdef RTFORCE
        #define RHO_DIVV(a,b) (b)
    #else
        #define RHO_DIVV(a,b) (a)
    #endif
#else
    #define DRHODTACTIVE(xxx)
#endif

#ifdef DIFFUSION

    #ifdef FEEDBACKDIFFLIMIT
        #define DIFFUSIONLimitTest() (diffSum == 0 || smf->dTime < p->fTimeCoolIsOffUntil || smf->dTime < q->fTimeCoolIsOffUntil)
    #else
        #define DIFFUSIONLimitTest() (diffSum == 0)
    #endif


        #define DIFFUSIONBase() double diffSum = (p->diff+q->diff); \
                                double diffBase = (DIFFUSIONLimitTest() ? 0 : diffSum);
        #define MASSDIFFFAC(pother_) 1
        #define DIFFUSIONMetalsBase() double diffMetalsBase = 2*smf->dMetalDiffusionCoeff*diffBase \
             /(p->fDensity+q->fDensity);
        #define DIFFUSIONMass()
        #define DIFFUSIONVelocity()
#else
    #define DIFFUSIONBase()
    #define DIFFUSIONMetalsBase()
    #define DIFFUSIONMass()
    #define DIFFUSIONVelocity()
#endif


#ifdef DIFFUSION
    #if defined(UNONCOOL) && !defined(TWOPHASE)
        #define DIFFUSIONThermaluHot() \
                { double diffuNc = diffTh*(p->uHotPred-q->uHotPred); \
                PACTIVE( p->uHotDotDiff += diffuNc*rq );        \
                QACTIVE( q->uHotDotDiff -= diffuNc*rp );        \
                }
    #else
        #define DIFFUSIONThermaluHot()
    #endif
    #ifdef DIFFUSIONPRICE
        #define DIFFUSIONThermal(dt_) \
            { double irhobar = 2/(p->fDensity+q->fDensity);     \
             double vsig = sqrt(fabs(qPoverRho2*q->fDensity*q->fDensity - pPoverRho2*p->fDensity*p->fDensity)*irhobar); \
             double diffTh = smf->dThermalDiffusionCoeff*0.5*(ph+sqrt(0.25*BALL2(q)))*irhobar*vsig; \
             double diffu = diffTh*(p->uPred-q->uPred);             \
             PACTIVE( p->uDotDiff += diffu*rq );                     \
             QACTIVE( q->uDotDiff-= diffu*rp );                     \
             DIFFUSIONThermaluHot(); }
    #else
        #ifndef NODIFFUSIONTHERMAL
            /* Default -- thermal diffusion */
            #ifdef THERMALCOND
                    /* Harmonic average coeff */
                    #define DIFFUSIONThermalCondBase(dt_) double dThermalCondSum = p->fThermalCond + q->fThermalCond; \
                        double dThermalCond = ( dThermalCondSum <= 0 ? 0 : 4*p->fThermalCond*q->fThermalCond/(dThermalCondSum*p->fDensity*q->fDensity) ); \
                        if (dThermalCond > 0 && (dt_diff = smf->dtFacDiffusion*ph*ph/(dThermalCond*p->fDensity)) < dt_) dt_ = dt_diff;
            #else
                #define DIFFUSIONThermalCondBase(dt_) double dThermalCond=0;
            #endif
                #define DIFFUSIONShockCondBase() double dShockCond = 0;
                #define DIFFUSIONShockCond() double dShockCond = 0;

                #define DIFFUSIONThermal(dt_) \
                    { double diffTh = (2*smf->dThermalDiffusionCoeff*diffBase/(p->fDensity+q->fDensity)); \
                      double dt_diff, diffu;                                                  \
                      DIFFUSIONThermalCondBase(dt_);                                    \
                      if (diffTh > 0 && (dt_diff= smf->dtFacDiffusion*ph*ph/(diffTh*p->fDensity)) < dt_) dt_ = dt_diff; \
                      diffu = (diffTh+dThermalCond+dShockCond)*(p->uPred-q->uPred);              \
                      PACTIVE( p->uDotDiff += diffu*rq*MASSDIFFFAC(q) );                \
                      QACTIVE( q->uDotDiff -= diffu*rp*MASSDIFFFAC(p) );                \
                      DIFFUSIONThermaluHot(); }
        #else
            #define DIFFUSIONThermal(dt_)
        #endif
    #endif

    #define DIFFUSIONMetals() \
        { double diff = diffMetalsBase*(p->fMetals - q->fMetals); \
          PACTIVE( p->fMetalsDot += diff*rq*MASSDIFFFAC(q) ); \
          QACTIVE( q->fMetalsDot -= diff*rp*MASSDIFFFAC(p) ); }
    #ifdef STARFORM
        #define DIFFUSIONMetalsOxygen() \
            { double diff = diffMetalsBase*(p->fMFracOxygen - q->fMFracOxygen); \
              PACTIVE( p->fMFracOxygenDot += diff*rq*MASSDIFFFAC(q) ); \
              QACTIVE( q->fMFracOxygenDot -= diff*rp*MASSDIFFFAC(p) ); }
        #define DIFFUSIONMetalsIron() \
            { double diff = diffMetalsBase*(p->fMFracIron - q->fMFracIron); \
              PACTIVE( p->fMFracIronDot += diff*rq*MASSDIFFFAC(q) ); \
              QACTIVE( q->fMFracIronDot -= diff*rp*MASSDIFFFAC(p) ); }
    #else
        #define DIFFUSIONMetalsOxygen()
        #define DIFFUSIONMetalsIron()
    #endif /* STARFORM */
#else /* No diffusion */
    #define DIFFUSIONShockCondBase()
    #define DIFFUSIONShockCond()
    #define DIFFUSIONThermal(dt_)
    #define DIFFUSIONMetals()
    #define DIFFUSIONMetalsOxygen()
    #define DIFFUSIONMetalsIron()
#endif

#if defined(VARALPHA) || defined(CULLENDEHNEN)
    #define ALPHA (smf->alpha*0.5*(p->alpha+q->alpha))
    #define BETA  (smf->beta*0.5*(p->alpha+q->alpha))
#else
    #define ALPHA (smf->alpha)
    #define BETA  (smf->beta)
#endif
#define SETDTNEW_PQ(dt_)  { if (dt_ < p->dtNew) p->dtNew=dt_; \
                            if (dt_ < q->dtNew) q->dtNew=dt_; \
                            if (4*q->dt < p->dtNew) p->dtNew = 4*q->dt; \
                            if (4*p->dt < q->dtNew) q->dtNew = 4*p->dt; }

#ifdef VSIGVISC
    #define ARTIFICIALVISCOSITY(visc_,dt_) { absmu = -dvdotdr*smf->a            \
                /sqrt(nnList[i].fDist2); /* mu multiply by a to be consistent with physical c */ \
            if (absmu>p->mumax) p->mumax=absmu; /* mu terms for gas time step */ \
    		if (absmu>q->mumax) q->mumax=absmu; \
    		visc_ = (ALPHA*(pc + q->c) + BETA*1.5*absmu); \
    		dt_ = smf->dtFacCourant*ph/(0.625*(pc + q->c)+0.375*visc_);     \
    		visc_ = SWITCHCOMBINE(p,q)*visc_ \
    		    *absmu/(pDensity + q->fDensity); }
#else
    #define ARTIFICIALVISCOSITY(visc_,dt_) { double hav=0.5*(ph+sqrt(0.25*BALL2(q)));  /* h mean */ \
    		absmu = -hav*dvdotdr*smf->a  \
    		    /(nnList[i].fDist2+0.01*hav*hav); /* mu multiply by a to be consistent with physical c */ \
    		if (absmu>p->mumax) p->mumax=absmu; /* mu terms for gas time step */ \
    		if (absmu>q->mumax) q->mumax=absmu; \
    		visc_ = (ALPHA*(pc + q->c) + BETA*2*absmu);	\
    		dt_ = smf->dtFacCourant*hav/(0.625*(pc + q->c)+0.375*visc_); \
    		visc_ = SWITCHCOMBINE(p,q)*visc_ \
    		    *absmu/(pDensity + q->fDensity); }
#endif

    /* Force Calculation between particles p and q */
DRHODTACTIVE( PACTIVE( p->fDivv_PdV -= rq/p->fDivv_Corrector/RHO_DIVV(pDensity,q->fDensity)*dvdotdr; ));
DRHODTACTIVE( QACTIVE( q->fDivv_PdV -= rp/p->fDivv_Corrector/RHO_DIVV(q->fDensity,pDensity)*dvdotdr; ));
DRHODTACTIVE( PACTIVE( p->fDivv_PdVcorr -= rq/RHO_DIVV(pDensity,q->fDensity)*dvdotdr; ));
DRHODTACTIVE( QACTIVE( q->fDivv_PdVcorr -= rp/RHO_DIVV(q->fDensity,pDensity)*dvdotdr; ));
PACTIVE( p->uDotPdV += rq*PRES_PDV(pPoverRho2,qPoverRho2)*dvdotdr; );
QACTIVE( q->uDotPdV += rp*PRES_PDV(qPoverRho2,pPoverRho2)*dvdotdr; );
//if (p->iOrder == 865177) printf("sphp %d: %g %g %g %g\n",p->iOrder,p->u,pPoverRho2,dvdotdr,rq*PRES_PDV(pPoverRho2,qPoverRho2)*dvdotdr );
//if (q->iOrder == 865177) printf("sphq %q: %g %g %g %g\n",q->iOrder,q->u,qPoverRho2,dvdotdr,rp*PRES_PDV(qPoverRho2,pPoverRho2)*dvdotdr );
PACTIVE( Accp = (PRES_ACC(pPoverRho2f,qPoverRho2f)); );
QACTIVE( Accq = (PRES_ACC(qPoverRho2f,pPoverRho2f)); );
//DIFFUSIONShockCondBase();
if (dvdotdr>=0.0) {
    dt = smf->dtFacCourant*ph/(2*(pc > q->c ? pc : q->c));
    }
else {
    ARTIFICIALVISCOSITY(visc,dt); /* Calculate Artificial viscosity terms and associated dt */
        PACTIVE( p->uDotAV += rq*(0.5*visc)*dvdotdr; );
        QACTIVE( q->uDotAV += rp*(0.5*visc)*dvdotdr; );
        PACTIVE( Accp += visc; );
        QACTIVE( Accq += visc; );
	}
PACTIVE( Accp *= rq*aFac; );/* aFac - convert to comoving acceleration */
QACTIVE( Accq *= rp*aFac; );
PACTIVE( ACCEL(p,0) -= Accp * dx; );
PACTIVE( ACCEL(p,1) -= Accp * dy; );
PACTIVE( ACCEL(p,2) -= Accp * dz; );
QACTIVE( ACCEL(q,0) += Accq * dx; );
QACTIVE( ACCEL(q,1) += Accq * dy; );
QACTIVE( ACCEL(q,2) += Accq * dz; );

DIFFUSIONBase();
DIFFUSIONShockCondBase();
if (dvdotdr < 0.0) {DIFFUSIONShockCond();}
DIFFUSIONThermal(dt);
DIFFUSIONMetalsBase();
DIFFUSIONMetals();
DIFFUSIONMetalsOxygen();
DIFFUSIONMetalsIron();
DIFFUSIONMass();
DIFFUSIONVelocity();
SETDTNEW_PQ(dt);
