#include "pkd.h"
#include "cooling.h"

void print_coolparticle(COOLPARTICLE cp) {
	printf("* %9s: %13e %10s: %13e %10s: %13e *\n", "f_HI", cp.f_HI, "f_HeI", cp.f_HeI, "F_HeII", cp.f_HeII);
}
void print_particle(PKD pkd,PARTICLE p) {
	printf("********************************************************************************\n");
	printf("*                                                                              *\n");
	if (pkdIsGasByOrder(pkd, &p)) printf("*                             GAS PARTICLE DEBUG                               *\n");
	if (pkdIsStarByOrder(pkd, &p)) printf("*                            STAR PARTICLE DEBUG                               *\n");
	if (pkdIsDarkByOrder(pkd, &p)) printf("*                            DARK PARTICLE DEBUG                               *\n");
	printf("*                                                                              *\n");
	printf("*         iOrder: %9d iActive: %9d iRung: %2d cpStart: %3d          *\n", p.iOrder, p.iActive, p.iRung, p.cpStart);
	printf("*                                                                              *\n");
	printf("* %9s: %13e %10s: %13e %10s: %13e *\n", "fWeight", p.fWeight, "fMass", p.fMass, "fSoft", p.fSoft);
	printf("*                                                                              *\n");
	printf("* %9s: %13e %10s: %13e %10s: %13e *\n", "r_x", p.r[0], "r_y", p.r[1], "r_z", p.r[2]);
	printf("* %9s: %13e %10s: %13e %10s: %13e *\n", "v_x", p.v[0], "v_y", p.v[1], "v_z", p.v[2]);
	printf("* %9s: %13e %10s: %13e %10s: %13e *\n", "a_x", p.a[0], "a_y", p.a[1], "a_z", p.a[2]);
	printf("*                                                                              *\n");
	printf("* %9s: %13e %10s: %13e %10s: %13e *\n", "dt", p.dt, "dtNew", p.dtNew, "dtOld", p.dtOld);
	printf("* %9s: %13e %10s: %13e %10s: %13e *\n", "dtGrav", p.dtGrav, "fPot", p.fPot, "fDensity", p.fDensity);
	printf("* %9s: %13e %10s: %13e                           *\n", "fBall2", p.fBall2, "fBallMax", p.fBallMax);
	printf("*                                                                              *\n");
#ifdef GASOLINE
	printf("* %9s: %13e %10s: %13e %10s: %13e *\n", "u", p.u, "uPred", p.uPred, "uDotAV", p.uDotAV);
	printf("* %9s: %13e %10s: %13e %10s: %13e *\n", "PoverRho2", p.PoverRho2, "mumax", p.mumax, "uDotPdV", p.uDotPdV);
	printf("* %9s: %13e %10s: %13e %10s: %13e *\n", "curlv_x", p.curlv[0], "curlv_y", p.curlv[1], "curlv_z", p.curlv[2]);
	printf("* %9s: %13e %10s: %13e %10s: %13e *\n", "divv", p.divv, "fDivv_Corr", p.fDivv_Corrector, "fDivv_t", p.fDivv_t);
	printf("* %9s: %13e %10s: %13e %10s: %13e *\n", "fMetals", p.fMetals, "fTimeForm", p.fTimeForm, "c", p.c);
	printf("* %9s: %13e                                                     *\n", "uDotDiff", p.uDotDiff); 
	print_coolparticle(p.CoolParticle);
#ifdef DIFFUSION
	printf("* %9s: %13e %10s: %13e %10s: %13e *\n", "diff", p.diff, "MetalsDot", p.fMetalsDot, "MetalsPred", p.fMetalsPred);
#ifdef MASSDIFF
	printf("* %9s: %13e %10s: %13e                           *\n", "fMassDot", p.fMassDot, "fMass0", p.fMass0);
#endif
#ifdef DRHODT
	printf("* %9s: %13e %10s: %13e %10s: %13e *\n", "fDens_t", p.fDensity_t, "fDens_PdV", p.fDensity_PdV, "fDens_PdVc", p.fDensity_PdVcorr);
	printf("* %9s: %13e %10s: %13e                           *\n", "fDivv_PdV", p.fDivv_PdV, "fDivv_PdVc", p.fDivv_PdVcorr);
#endif
#endif
	printf("* %9s: %13e                                                     *\n", "BalsaraSW", p.BalsaraSwitch); 
	printf("*                                                                              *\n");
#endif
#ifdef UNONCOOL
	printf("* %9s: %13e %10s: %13e %10s: %13e *\n", "uNoncool", p.uNoncool, "uNCPred", p.uNoncoolPred, "uNCDot", p.uNoncoolDot);
	printf("* %9s: %13e                                                     *\n", "uNCDotDiff", p.uNoncoolDotDiff); 
	printf("*                                                                              *\n");
#endif
#ifdef SHOCKTRACK
	printf("* %9s: %13e %10s: %13e %10s: %13e *\n", "aPres_x", p.aPres[0], "aPres_y", p.aPres[1], "aPres_z", p.aPres[2]);
	printf("* %9s: %13e %10s: %13e %10s: %13e *\n", "gradrho_x", p.gradrho[0], "gradrho_y", p.gradrho[1], "gradrho_z", p.gradrho[2]);
	printf("* %9s: %13e %10s: %13e                           *\n", "divrhov", p.divrhov, "ShockTrack", p.ShockTracker);
#endif
#ifdef STARFORM
	printf("* %9s: %13e %10s: %13e %10s: %13e *\n", "uDotFB", p.uDotFB, "fMSN", p.fMSN, "fNSN", p.fNSN);
	printf("* %9s: %13e %10s: %13e %10s: %13e *\n", "fMOOut", p.fMOxygenOut, "fMFeOut", p.fMIronOut, "fMFracO", p.fMFracOxygen);
	printf("* %9s: %13e %10s: %13e %10s: %13e *\n", "fMFracFe", p.fMFracIron, "fSNMetals", p.fSNMetals, "fNSNtot", p.fNSNtot);
	printf("* %9s: %13e %10s: %13e %10s: %13d *\n", "fCoolOff", p.fTimeCoolIsOffUntil, "fMassForm", p.fMassForm, "iGasOrder", p.iGasOrder);
#ifdef DIFFUSION
	printf("* %9s: %13e %10s: %13e %10s: %13e *\n", "MFracODot", p.fMFracOxygenDot, "MFracFeDot", p.fMFracIronDot, "MFracOPred", p.fMFracOxygenPred);
	printf("* %9s: %13e                                                     *\n", "MFrFePred", p.fMFracIronPred); 
#endif
	printf("*                                                                              *\n");
#endif
#ifdef SINKING
	printf("* %9s: %13e %10s: %13e %10s: %13e *\n", "rS0Unit_x", p.rSinking0Unit[0], "rS0Unit_y", p.rSinking0Unit[0], "rS0Unit_z", p.rSinking0Unit[0]);
	printf("* %9s: %13e %10s: %13e %10s: %13e *\n", "vSTang0_x", p.vSinkingTang0Unit[0], "vSTang0_y", p.vSinkingTang0Unit[0], "vSTang0_z", p.vSinkingTang0Unit[0]);
	printf("* %9s: %13e %10s: %13e %10s: %13e *\n", "rS0Mag", p.rSinking0Mag, "vSTang0Mag", p.rSinkingTang0Mag, "vSinkr0", p.vSinkingr0);
	printf("* %9s: %13e %10s: %13e %10s: %13d *\n", "SinkTime", p.fSinkingTIme, "fTrueMass", p.fTrueMass, "SinkOnto", p.iSinginkOnto);
	printf("*                                                                              *\n");
#endif
#ifdef SIMPLESF
	printf("* %9s: %13e %10s: %13e %10s: %13d *\n", "fESN", p.fESN, "fMassStar", p.fMassStar, "iGasOrder", p.iGasOrder);
#endif
#if defined(SIMPLESF) || defined(EXTRASINKDATA)
	printf("* %9s: %13e %10s: %13e %10s: %13e *\n", "rForm_x", p.rForm[0], "rForm_y", p.rForm[1], "rForm_z", p.rForm[2]);
	printf("* %9s: %13e %10s: %13e %10s: %13e *\n", "vForm_x", p.vForm[0], "vForm_y", p.vForm[1], "vForm_z", p.vForm[2]);
	printf("*                                                                              *\n");
#endif
#ifdef CHECKSF
	printf("* %9s: %13e %10s: %13e %10s: %13e *\n", "tOff", p.tOff, "tcool", p.tcool, "tdyn", p.tdyn);
	printf("* %9s: %13e %10s: %13e %10s: %13d *\n", "l_jeans", p.tOff, "sound/dyn", p.ratiosounddyn, "small_jeans", p.small_jeans);
#endif
#ifdef COLLISIONS
	printf("* %9s: %13d %10s: %13d %10s: %13d *\n", "iOrgIdx", p.iOrgIdx, "iColor", p.iColor, "iDriftType", p.iDriftType);
	printf("* %9s: %13d %10s: %13d %10s: %13d *\n", "iOrderCol", p.iOrderCol, "iPrevCol", p.iPrevCol, "bTinyStep", p.bTinyStep);
	printf("* %9s: %13e %10s: %13e %10s: %13e *\n", "w_x", p.w[0], "w_y", p.w[1], "w_z", p.w[2]);
	printf("* %9s: %13e                                                     *\n", "mindist2", p.mindist2); 
	printf("*                                                                              *\n");
#endif
#ifdef RUBBLE_ZML
	printf("* %9s: %13e %10s: %13d %10s: %13d *\n", "dDustMass", p.dDustMass, "iBin", p.iBin, "MayCollide", p.bMayCollide);
#endif
#ifdef SPECIAL_PARTICLES
	printf("*     SPECIAL_PARTICLES:                    bGhostExclue: %13d        *\n", p.bGhostExclude);
#endif
#ifdef CHANGESOFT
	printf("*     CHANGESOFT:                                 fSoft0: %13e        *\n", p.fSoft0);
#endif
#ifdef DENSITYU
	printf("*       DENSITYU:                              fDensityU: %13e        *\n", p.fDensityU);
#endif
#ifdef DODVDS
	printf("*     DODVDS:                                       dvds: %13e        *\n", p.dvds);
#endif
#ifdef SLIDING_PATCH
	printf("*     SLIDING_PATCH:           bAzWrap %7d       dPy: %13e        *\n", p.dPy, p.bAzWrap);
#endif
#ifdef SURFACE_AREA
	printf("*     SURFACE_AREA:                                fArea: %13e        *\n", p.fArea);
#endif
#ifdef NEED_VPRED
	printf("*    NEED_VPRED: vPred: %15e %15e %15e        * \n", p.vPred[0], p.vPred[1], p.vPred[2]);
#endif
#ifdef SUPERCOOL
	printf("*     SUPERCOOL: vMean: %15e %15e %15e        * \n", p.vMean[0], p.vMean[1], p.vMean[2]);
#endif
#ifdef NORMAL
	printf("*       NORMAL: normal: %15e %15e %15e        * \n", p.normal[0], p.normal[1], p.normal[2]);
#endif
#ifdef AGGS
	printf("*       AGGS:    r_agg: %15e %15e %15e        * \n", p.r_agg[0], p.r_agg[1], p.r_agg[2]);
#endif
#ifdef SAND_PILE
	printf("*     SAND_PILE:                                  bStuck: %13d        *\n", p.bStuck);
#endif
#ifdef COLORCODE
	printf("*     COLORCODE:                                  fColor: %13e        *\n", p.fColor);
#endif
#ifdef PDVDEBUG
	printf("*     PDVDEBUG:     PdVvisc: %19e PdVpres: %18e   *\n", p.PdVvisc, p.PdVpres);
#endif
#ifdef VARALPHA
	printf("*     VARALPHA:       alpha: %17e alphaPred: %18e   *\n", p.fSoft, p.alphaPred);
#endif
	printf("*                                                                              *\n");
	printf("********************************************************************************\n");
	}
