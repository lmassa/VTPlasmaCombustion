#ifndef _PLASMAPHYSICSF_F_H_
#define _PLASMAPHYSICSF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDPSOURCEFAB 
#define GUARDPSOURCEFAB 
// Prototype for Fortran procedure PSOURCEFAB ...
//
void FORTRAN_NAME( PSOURCEFAB ,psourcefab )(
      CHFp_FRA(dst)
      ,CHFp_CONST_FRA(x)
      ,CHFp_CONST_FRA(y)
      ,CHFp_CONST_VR(coef)
      ,CHFp_CONST_INT(destcomp)
      ,CHFp_BOX(region) );

#define FORT_PSOURCEFAB FORTRAN_NAME( inlinePSOURCEFAB, inlinePSOURCEFAB)
#define FORTNT_PSOURCEFAB FORTRAN_NAME( PSOURCEFAB, psourcefab)

inline void FORTRAN_NAME(inlinePSOURCEFAB, inlinePSOURCEFAB)(
      CHFp_FRA(dst)
      ,CHFp_CONST_FRA(x)
      ,CHFp_CONST_FRA(y)
      ,CHFp_CONST_VR(coef)
      ,CHFp_CONST_INT(destcomp)
      ,CHFp_BOX(region) )
{
 CH_TIMELEAF("FORT_PSOURCEFAB");
 FORTRAN_NAME( PSOURCEFAB ,psourcefab )(
      CHFt_FRA(dst)
      ,CHFt_CONST_FRA(x)
      ,CHFt_CONST_FRA(y)
      ,CHFt_CONST_VR(coef)
      ,CHFt_CONST_INT(destcomp)
      ,CHFt_BOX(region) );
}
#endif  // GUARDPSOURCEFAB 

#ifndef GUARDPESOURCEFAB 
#define GUARDPESOURCEFAB 
// Prototype for Fortran procedure PESOURCEFAB ...
//
void FORTRAN_NAME( PESOURCEFAB ,pesourcefab )(
      CHFp_FRA(dst)
      ,CHFp_CONST_FRA(rr)
      ,CHFp_CONST_FRA(nd)
      ,CHFp_CONST_FRA(te)
      ,CHFp_CONST_FRA(jh)
      ,CHFp_CONST_VR(ecoef)
      ,CHFp_CONST_VR(mass)
      ,CHFp_CONST_VR(sigebg)
      ,CHFp_CONST_VR(sigeo2)
      ,CHFp_CONST_VI(spec)
      ,CHFp_CONST_REAL(nbg)
      ,CHFp_CONST_REAL(mte)
      ,CHFp_CONST_INT(tgcomp)
      ,CHFp_CONST_INT(eind)
      ,CHFp_CONST_INT(bgind)
      ,CHFp_CONST_INT(o2ind)
      ,CHFp_CONST_INT(cstte)
      ,CHFp_CONST_INT(cstbg)
      ,CHFp_BOX(region) );

#define FORT_PESOURCEFAB FORTRAN_NAME( inlinePESOURCEFAB, inlinePESOURCEFAB)
#define FORTNT_PESOURCEFAB FORTRAN_NAME( PESOURCEFAB, pesourcefab)

inline void FORTRAN_NAME(inlinePESOURCEFAB, inlinePESOURCEFAB)(
      CHFp_FRA(dst)
      ,CHFp_CONST_FRA(rr)
      ,CHFp_CONST_FRA(nd)
      ,CHFp_CONST_FRA(te)
      ,CHFp_CONST_FRA(jh)
      ,CHFp_CONST_VR(ecoef)
      ,CHFp_CONST_VR(mass)
      ,CHFp_CONST_VR(sigebg)
      ,CHFp_CONST_VR(sigeo2)
      ,CHFp_CONST_VI(spec)
      ,CHFp_CONST_REAL(nbg)
      ,CHFp_CONST_REAL(mte)
      ,CHFp_CONST_INT(tgcomp)
      ,CHFp_CONST_INT(eind)
      ,CHFp_CONST_INT(bgind)
      ,CHFp_CONST_INT(o2ind)
      ,CHFp_CONST_INT(cstte)
      ,CHFp_CONST_INT(cstbg)
      ,CHFp_BOX(region) )
{
 CH_TIMELEAF("FORT_PESOURCEFAB");
 FORTRAN_NAME( PESOURCEFAB ,pesourcefab )(
      CHFt_FRA(dst)
      ,CHFt_CONST_FRA(rr)
      ,CHFt_CONST_FRA(nd)
      ,CHFt_CONST_FRA(te)
      ,CHFt_CONST_FRA(jh)
      ,CHFt_CONST_VR(ecoef)
      ,CHFt_CONST_VR(mass)
      ,CHFt_CONST_VR(sigebg)
      ,CHFt_CONST_VR(sigeo2)
      ,CHFt_CONST_VI(spec)
      ,CHFt_CONST_REAL(nbg)
      ,CHFt_CONST_REAL(mte)
      ,CHFt_CONST_INT(tgcomp)
      ,CHFt_CONST_INT(eind)
      ,CHFt_CONST_INT(bgind)
      ,CHFt_CONST_INT(o2ind)
      ,CHFt_CONST_INT(cstte)
      ,CHFt_CONST_INT(cstbg)
      ,CHFt_BOX(region) );
}
#endif  // GUARDPESOURCEFAB 

#ifndef GUARDFSOURCEFAB 
#define GUARDFSOURCEFAB 
// Prototype for Fortran procedure FSOURCEFAB ...
//
void FORTRAN_NAME( FSOURCEFAB ,fsourcefab )(
      CHFp_FRA(dst)
      ,CHFp_CONST_FRA(rr)
      ,CHFp_CONST_FRA(nd)
      ,CHFp_CONST_FRA(te)
      ,CHFp_CONST_FRA(gp)
      ,CHFp_CONST_FRA(velo)
      ,CHFp_CONST_VR(gcoef)
      ,CHFp_CONST_VR(zcoef)
      ,CHFp_CONST_VR(mass)
      ,CHFp_CONST_VR(sigebg)
      ,CHFp_CONST_VR(sigeo2)
      ,CHFp_CONST_VI(spec)
      ,CHFp_CONST_REAL(nbg)
      ,CHFp_CONST_REAL(mte)
      ,CHFp_CONST_REAL(eta)
      ,CHFp_CONST_INT(tgcomp)
      ,CHFp_CONST_INT(eind)
      ,CHFp_CONST_INT(bgind)
      ,CHFp_CONST_INT(o2ind)
      ,CHFp_CONST_INT(cstbg)
      ,CHFp_CONST_INT(cstte)
      ,CHFp_CONST_INT(mombegin)
      ,CHFp_CONST_INT(enrgIndx)
      ,CHFp_BOX(region) );

#define FORT_FSOURCEFAB FORTRAN_NAME( inlineFSOURCEFAB, inlineFSOURCEFAB)
#define FORTNT_FSOURCEFAB FORTRAN_NAME( FSOURCEFAB, fsourcefab)

inline void FORTRAN_NAME(inlineFSOURCEFAB, inlineFSOURCEFAB)(
      CHFp_FRA(dst)
      ,CHFp_CONST_FRA(rr)
      ,CHFp_CONST_FRA(nd)
      ,CHFp_CONST_FRA(te)
      ,CHFp_CONST_FRA(gp)
      ,CHFp_CONST_FRA(velo)
      ,CHFp_CONST_VR(gcoef)
      ,CHFp_CONST_VR(zcoef)
      ,CHFp_CONST_VR(mass)
      ,CHFp_CONST_VR(sigebg)
      ,CHFp_CONST_VR(sigeo2)
      ,CHFp_CONST_VI(spec)
      ,CHFp_CONST_REAL(nbg)
      ,CHFp_CONST_REAL(mte)
      ,CHFp_CONST_REAL(eta)
      ,CHFp_CONST_INT(tgcomp)
      ,CHFp_CONST_INT(eind)
      ,CHFp_CONST_INT(bgind)
      ,CHFp_CONST_INT(o2ind)
      ,CHFp_CONST_INT(cstbg)
      ,CHFp_CONST_INT(cstte)
      ,CHFp_CONST_INT(mombegin)
      ,CHFp_CONST_INT(enrgIndx)
      ,CHFp_BOX(region) )
{
 CH_TIMELEAF("FORT_FSOURCEFAB");
 FORTRAN_NAME( FSOURCEFAB ,fsourcefab )(
      CHFt_FRA(dst)
      ,CHFt_CONST_FRA(rr)
      ,CHFt_CONST_FRA(nd)
      ,CHFt_CONST_FRA(te)
      ,CHFt_CONST_FRA(gp)
      ,CHFt_CONST_FRA(velo)
      ,CHFt_CONST_VR(gcoef)
      ,CHFt_CONST_VR(zcoef)
      ,CHFt_CONST_VR(mass)
      ,CHFt_CONST_VR(sigebg)
      ,CHFt_CONST_VR(sigeo2)
      ,CHFt_CONST_VI(spec)
      ,CHFt_CONST_REAL(nbg)
      ,CHFt_CONST_REAL(mte)
      ,CHFt_CONST_REAL(eta)
      ,CHFt_CONST_INT(tgcomp)
      ,CHFt_CONST_INT(eind)
      ,CHFt_CONST_INT(bgind)
      ,CHFt_CONST_INT(o2ind)
      ,CHFt_CONST_INT(cstbg)
      ,CHFt_CONST_INT(cstte)
      ,CHFt_CONST_INT(mombegin)
      ,CHFt_CONST_INT(enrgIndx)
      ,CHFt_BOX(region) );
}
#endif  // GUARDFSOURCEFAB 

#ifndef GUARDREACRATEFAB 
#define GUARDREACRATEFAB 
// Prototype for Fortran procedure REACRATEFAB ...
//
void FORTRAN_NAME( REACRATEFAB ,reacratefab )(
      CHFp_FRA(dst)
      ,CHFp_CONST_FRA(x)
      ,CHFp_CONST_FRA(y)
      ,CHFp_CONST_FRA(gp)
      ,CHFp_CONST_VR(coef)
      ,CHFp_CONST_VR(stoe)
      ,CHFp_CONST_VR(actreac)
      ,CHFp_CONST_VR(transpe)
      ,CHFp_CONST_VI(spec)
      ,CHFp_CONST_REAL(press)
      ,CHFp_CONST_REAL(nbg)
      ,CHFp_CONST_REAL(mte)
      ,CHFp_CONST_INT(destcomp)
      ,CHFp_CONST_INT(varcomp)
      ,CHFp_CONST_INT(ecomp)
      ,CHFp_CONST_INT(bgind)
      ,CHFp_CONST_INT(o2ind)
      ,CHFp_CONST_INT(cstbg)
      ,CHFp_CONST_INT(cstte)
      ,CHFp_BOX(region) );

#define FORT_REACRATEFAB FORTRAN_NAME( inlineREACRATEFAB, inlineREACRATEFAB)
#define FORTNT_REACRATEFAB FORTRAN_NAME( REACRATEFAB, reacratefab)

inline void FORTRAN_NAME(inlineREACRATEFAB, inlineREACRATEFAB)(
      CHFp_FRA(dst)
      ,CHFp_CONST_FRA(x)
      ,CHFp_CONST_FRA(y)
      ,CHFp_CONST_FRA(gp)
      ,CHFp_CONST_VR(coef)
      ,CHFp_CONST_VR(stoe)
      ,CHFp_CONST_VR(actreac)
      ,CHFp_CONST_VR(transpe)
      ,CHFp_CONST_VI(spec)
      ,CHFp_CONST_REAL(press)
      ,CHFp_CONST_REAL(nbg)
      ,CHFp_CONST_REAL(mte)
      ,CHFp_CONST_INT(destcomp)
      ,CHFp_CONST_INT(varcomp)
      ,CHFp_CONST_INT(ecomp)
      ,CHFp_CONST_INT(bgind)
      ,CHFp_CONST_INT(o2ind)
      ,CHFp_CONST_INT(cstbg)
      ,CHFp_CONST_INT(cstte)
      ,CHFp_BOX(region) )
{
 CH_TIMELEAF("FORT_REACRATEFAB");
 FORTRAN_NAME( REACRATEFAB ,reacratefab )(
      CHFt_FRA(dst)
      ,CHFt_CONST_FRA(x)
      ,CHFt_CONST_FRA(y)
      ,CHFt_CONST_FRA(gp)
      ,CHFt_CONST_VR(coef)
      ,CHFt_CONST_VR(stoe)
      ,CHFt_CONST_VR(actreac)
      ,CHFt_CONST_VR(transpe)
      ,CHFt_CONST_VI(spec)
      ,CHFt_CONST_REAL(press)
      ,CHFt_CONST_REAL(nbg)
      ,CHFt_CONST_REAL(mte)
      ,CHFt_CONST_INT(destcomp)
      ,CHFt_CONST_INT(varcomp)
      ,CHFt_CONST_INT(ecomp)
      ,CHFt_CONST_INT(bgind)
      ,CHFt_CONST_INT(o2ind)
      ,CHFt_CONST_INT(cstbg)
      ,CHFt_CONST_INT(cstte)
      ,CHFt_BOX(region) );
}
#endif  // GUARDREACRATEFAB 

#ifndef GUARDRHSFAB 
#define GUARDRHSFAB 
// Prototype for Fortran procedure RHSFAB ...
//
void FORTRAN_NAME( RHSFAB ,rhsfab )(
      CHFp_FRA1(dst)
      ,CHFp_CONST_FRA(x)
      ,CHFp_CONST_VR(coef)
      ,CHFp_CONST_VI(spec)
      ,CHFp_CONST_REAL(epsg)
      ,CHFp_CONST_REAL(epsd)
      ,CHFp_CONST_INT(ebj)
      ,CHFp_CONST_INT(varcomp)
      ,CHFp_BOX(region) );

#define FORT_RHSFAB FORTRAN_NAME( inlineRHSFAB, inlineRHSFAB)
#define FORTNT_RHSFAB FORTRAN_NAME( RHSFAB, rhsfab)

inline void FORTRAN_NAME(inlineRHSFAB, inlineRHSFAB)(
      CHFp_FRA1(dst)
      ,CHFp_CONST_FRA(x)
      ,CHFp_CONST_VR(coef)
      ,CHFp_CONST_VI(spec)
      ,CHFp_CONST_REAL(epsg)
      ,CHFp_CONST_REAL(epsd)
      ,CHFp_CONST_INT(ebj)
      ,CHFp_CONST_INT(varcomp)
      ,CHFp_BOX(region) )
{
 CH_TIMELEAF("FORT_RHSFAB");
 FORTRAN_NAME( RHSFAB ,rhsfab )(
      CHFt_FRA1(dst)
      ,CHFt_CONST_FRA(x)
      ,CHFt_CONST_VR(coef)
      ,CHFt_CONST_VI(spec)
      ,CHFt_CONST_REAL(epsg)
      ,CHFt_CONST_REAL(epsd)
      ,CHFt_CONST_INT(ebj)
      ,CHFt_CONST_INT(varcomp)
      ,CHFt_BOX(region) );
}
#endif  // GUARDRHSFAB 

#ifndef GUARDPRINT 
#define GUARDPRINT 
// Prototype for Fortran procedure print ...
//
void FORTRAN_NAME( PRINT ,print )(
      CHFp_FRA1(facemldata)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(facebox)
      ,CHFp_CONST_INT(flag) );

#define FORT_PRINT FORTRAN_NAME( inlinePRINT, inlinePRINT)
#define FORTNT_PRINT FORTRAN_NAME( PRINT, print)

inline void FORTRAN_NAME(inlinePRINT, inlinePRINT)(
      CHFp_FRA1(facemldata)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(facebox)
      ,CHFp_CONST_INT(flag) )
{
 CH_TIMELEAF("FORT_PRINT");
 FORTRAN_NAME( PRINT ,print )(
      CHFt_FRA1(facemldata)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_BOX(facebox)
      ,CHFt_CONST_INT(flag) );
}
#endif  // GUARDPRINT 

#ifndef GUARDTRANSPCOEF 
#define GUARDTRANSPCOEF 
// Prototype for Fortran procedure transpcoef ...
//
void FORTRAN_NAME( TRANSPCOEF ,transpcoef )(
      CHFp_FRA1(facemdata)
      ,CHFp_FRA1(facemldata)
      ,CHFp_FRA1(faceddata)
      ,CHFp_FRA1(facedldata)
      ,CHFp_CONST_FRA(cellpdata)
      ,CHFp_CONST_FRA1(facegpdata)
      ,CHFp_CONST_FRA(celltdata)
      ,CHFp_CONST_FRA(cellndata)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(facebox)
      ,CHFp_CONST_VR(coef)
      ,CHFp_CONST_VI(spec)
      ,CHFp_CONST_REAL(zsign)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(press)
      ,CHFp_CONST_REAL(nbg)
      ,CHFp_CONST_REAL(mte)
      ,CHFp_CONST_INT(tcomp)
      ,CHFp_CONST_INT(ecomp)
      ,CHFp_CONST_INT(bgcomp)
      ,CHFp_CONST_INT(o2comp)
      ,CHFp_CONST_INT(cstbg)
      ,CHFp_CONST_INT(cstte)
      ,CHFp_CONST_INT(s) );

#define FORT_TRANSPCOEF FORTRAN_NAME( inlineTRANSPCOEF, inlineTRANSPCOEF)
#define FORTNT_TRANSPCOEF FORTRAN_NAME( TRANSPCOEF, transpcoef)

inline void FORTRAN_NAME(inlineTRANSPCOEF, inlineTRANSPCOEF)(
      CHFp_FRA1(facemdata)
      ,CHFp_FRA1(facemldata)
      ,CHFp_FRA1(faceddata)
      ,CHFp_FRA1(facedldata)
      ,CHFp_CONST_FRA(cellpdata)
      ,CHFp_CONST_FRA1(facegpdata)
      ,CHFp_CONST_FRA(celltdata)
      ,CHFp_CONST_FRA(cellndata)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(facebox)
      ,CHFp_CONST_VR(coef)
      ,CHFp_CONST_VI(spec)
      ,CHFp_CONST_REAL(zsign)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(press)
      ,CHFp_CONST_REAL(nbg)
      ,CHFp_CONST_REAL(mte)
      ,CHFp_CONST_INT(tcomp)
      ,CHFp_CONST_INT(ecomp)
      ,CHFp_CONST_INT(bgcomp)
      ,CHFp_CONST_INT(o2comp)
      ,CHFp_CONST_INT(cstbg)
      ,CHFp_CONST_INT(cstte)
      ,CHFp_CONST_INT(s) )
{
 CH_TIMELEAF("FORT_TRANSPCOEF");
 FORTRAN_NAME( TRANSPCOEF ,transpcoef )(
      CHFt_FRA1(facemdata)
      ,CHFt_FRA1(facemldata)
      ,CHFt_FRA1(faceddata)
      ,CHFt_FRA1(facedldata)
      ,CHFt_CONST_FRA(cellpdata)
      ,CHFt_CONST_FRA1(facegpdata)
      ,CHFt_CONST_FRA(celltdata)
      ,CHFt_CONST_FRA(cellndata)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_BOX(facebox)
      ,CHFt_CONST_VR(coef)
      ,CHFt_CONST_VI(spec)
      ,CHFt_CONST_REAL(zsign)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_REAL(press)
      ,CHFt_CONST_REAL(nbg)
      ,CHFt_CONST_REAL(mte)
      ,CHFt_CONST_INT(tcomp)
      ,CHFt_CONST_INT(ecomp)
      ,CHFt_CONST_INT(bgcomp)
      ,CHFt_CONST_INT(o2comp)
      ,CHFt_CONST_INT(cstbg)
      ,CHFt_CONST_INT(cstte)
      ,CHFt_CONST_INT(s) );
}
#endif  // GUARDTRANSPCOEF 

#ifndef GUARDMAXDT 
#define GUARDMAXDT 
// Prototype for Fortran procedure maxDt ...
//
void FORTRAN_NAME( MAXDT ,maxdt )(
      CHFp_FRA1(source)
      ,CHFp_CONST_FRA1(Tc)
      ,CHFp_CONST_FRA(Nc)
      ,CHFp_CONST_REAL(epsg)
      ,CHFp_BOX(cellBox)
      ,CHFp_CONST_VR(coef)
      ,CHFp_CONST_REAL(nbg)
      ,CHFp_CONST_REAL(mte)
      ,CHFp_CONST_INT(eind)
      ,CHFp_CONST_INT(bgcomp)
      ,CHFp_CONST_INT(o2comp)
      ,CHFp_CONST_INT(cnbg)
      ,CHFp_CONST_INT(cte) );

#define FORT_MAXDT FORTRAN_NAME( inlineMAXDT, inlineMAXDT)
#define FORTNT_MAXDT FORTRAN_NAME( MAXDT, maxdt)

inline void FORTRAN_NAME(inlineMAXDT, inlineMAXDT)(
      CHFp_FRA1(source)
      ,CHFp_CONST_FRA1(Tc)
      ,CHFp_CONST_FRA(Nc)
      ,CHFp_CONST_REAL(epsg)
      ,CHFp_BOX(cellBox)
      ,CHFp_CONST_VR(coef)
      ,CHFp_CONST_REAL(nbg)
      ,CHFp_CONST_REAL(mte)
      ,CHFp_CONST_INT(eind)
      ,CHFp_CONST_INT(bgcomp)
      ,CHFp_CONST_INT(o2comp)
      ,CHFp_CONST_INT(cnbg)
      ,CHFp_CONST_INT(cte) )
{
 CH_TIMELEAF("FORT_MAXDT");
 FORTRAN_NAME( MAXDT ,maxdt )(
      CHFt_FRA1(source)
      ,CHFt_CONST_FRA1(Tc)
      ,CHFt_CONST_FRA(Nc)
      ,CHFt_CONST_REAL(epsg)
      ,CHFt_BOX(cellBox)
      ,CHFt_CONST_VR(coef)
      ,CHFt_CONST_REAL(nbg)
      ,CHFt_CONST_REAL(mte)
      ,CHFt_CONST_INT(eind)
      ,CHFt_CONST_INT(bgcomp)
      ,CHFt_CONST_INT(o2comp)
      ,CHFt_CONST_INT(cnbg)
      ,CHFt_CONST_INT(cte) );
}
#endif  // GUARDMAXDT 

#ifndef GUARDJOULEHEAT 
#define GUARDJOULEHEAT 
// Prototype for Fortran procedure JOULEHEAT ...
//
void FORTRAN_NAME( JOULEHEAT ,jouleheat )(
      CHFp_FRA1(source)
      ,CHFp_CONST_FRA(Ec)
      ,CHFp_CONST_FRA1(Tc)
      ,CHFp_CONST_FRA(Nc)
      ,CHFp_CONST_FRA(GNc)
      ,CHFp_CONST_FRA(Vel)
      ,CHFp_BOX(cellBox)
      ,CHFp_CONST_VR(coef)
      ,CHFp_CONST_REAL(nbg)
      ,CHFp_CONST_REAL(mte)
      ,CHFp_CONST_INT(eind)
      ,CHFp_CONST_INT(bgcomp)
      ,CHFp_CONST_INT(o2comp)
      ,CHFp_CONST_INT(cnbg)
      ,CHFp_CONST_INT(cte) );

#define FORT_JOULEHEAT FORTRAN_NAME( inlineJOULEHEAT, inlineJOULEHEAT)
#define FORTNT_JOULEHEAT FORTRAN_NAME( JOULEHEAT, jouleheat)

inline void FORTRAN_NAME(inlineJOULEHEAT, inlineJOULEHEAT)(
      CHFp_FRA1(source)
      ,CHFp_CONST_FRA(Ec)
      ,CHFp_CONST_FRA1(Tc)
      ,CHFp_CONST_FRA(Nc)
      ,CHFp_CONST_FRA(GNc)
      ,CHFp_CONST_FRA(Vel)
      ,CHFp_BOX(cellBox)
      ,CHFp_CONST_VR(coef)
      ,CHFp_CONST_REAL(nbg)
      ,CHFp_CONST_REAL(mte)
      ,CHFp_CONST_INT(eind)
      ,CHFp_CONST_INT(bgcomp)
      ,CHFp_CONST_INT(o2comp)
      ,CHFp_CONST_INT(cnbg)
      ,CHFp_CONST_INT(cte) )
{
 CH_TIMELEAF("FORT_JOULEHEAT");
 FORTRAN_NAME( JOULEHEAT ,jouleheat )(
      CHFt_FRA1(source)
      ,CHFt_CONST_FRA(Ec)
      ,CHFt_CONST_FRA1(Tc)
      ,CHFt_CONST_FRA(Nc)
      ,CHFt_CONST_FRA(GNc)
      ,CHFt_CONST_FRA(Vel)
      ,CHFt_BOX(cellBox)
      ,CHFt_CONST_VR(coef)
      ,CHFt_CONST_REAL(nbg)
      ,CHFt_CONST_REAL(mte)
      ,CHFt_CONST_INT(eind)
      ,CHFt_CONST_INT(bgcomp)
      ,CHFt_CONST_INT(o2comp)
      ,CHFt_CONST_INT(cnbg)
      ,CHFt_CONST_INT(cte) );
}
#endif  // GUARDJOULEHEAT 

#ifndef GUARDJOULEHEATFLUID 
#define GUARDJOULEHEATFLUID 
// Prototype for Fortran procedure JOULEHEATFLUID ...
//
void FORTRAN_NAME( JOULEHEATFLUID ,jouleheatfluid )(
      CHFp_FRA1(source)
      ,CHFp_CONST_FRA(Ec)
      ,CHFp_CONST_FRA1(Tc)
      ,CHFp_CONST_FRA(Nc)
      ,CHFp_CONST_FRA(GNc)
      ,CHFp_BOX(cellBox)
      ,CHFp_CONST_VR(coef)
      ,CHFp_CONST_INT(ecomp)
      ,CHFp_CONST_INT(bgcomp)
      ,CHFp_CONST_INT(o2comp)
      ,CHFp_CONST_INT(cnbg)
      ,CHFp_CONST_INT(cte)
      ,CHFp_CONST_INT(s)
      ,CHFp_CONST_REAL(zsign)
      ,CHFp_CONST_REAL(nbg)
      ,CHFp_CONST_REAL(mte)
      ,CHFp_CONST_REAL(etaT) );

#define FORT_JOULEHEATFLUID FORTRAN_NAME( inlineJOULEHEATFLUID, inlineJOULEHEATFLUID)
#define FORTNT_JOULEHEATFLUID FORTRAN_NAME( JOULEHEATFLUID, jouleheatfluid)

inline void FORTRAN_NAME(inlineJOULEHEATFLUID, inlineJOULEHEATFLUID)(
      CHFp_FRA1(source)
      ,CHFp_CONST_FRA(Ec)
      ,CHFp_CONST_FRA1(Tc)
      ,CHFp_CONST_FRA(Nc)
      ,CHFp_CONST_FRA(GNc)
      ,CHFp_BOX(cellBox)
      ,CHFp_CONST_VR(coef)
      ,CHFp_CONST_INT(ecomp)
      ,CHFp_CONST_INT(bgcomp)
      ,CHFp_CONST_INT(o2comp)
      ,CHFp_CONST_INT(cnbg)
      ,CHFp_CONST_INT(cte)
      ,CHFp_CONST_INT(s)
      ,CHFp_CONST_REAL(zsign)
      ,CHFp_CONST_REAL(nbg)
      ,CHFp_CONST_REAL(mte)
      ,CHFp_CONST_REAL(etaT) )
{
 CH_TIMELEAF("FORT_JOULEHEATFLUID");
 FORTRAN_NAME( JOULEHEATFLUID ,jouleheatfluid )(
      CHFt_FRA1(source)
      ,CHFt_CONST_FRA(Ec)
      ,CHFt_CONST_FRA1(Tc)
      ,CHFt_CONST_FRA(Nc)
      ,CHFt_CONST_FRA(GNc)
      ,CHFt_BOX(cellBox)
      ,CHFt_CONST_VR(coef)
      ,CHFt_CONST_INT(ecomp)
      ,CHFt_CONST_INT(bgcomp)
      ,CHFt_CONST_INT(o2comp)
      ,CHFt_CONST_INT(cnbg)
      ,CHFt_CONST_INT(cte)
      ,CHFt_CONST_INT(s)
      ,CHFt_CONST_REAL(zsign)
      ,CHFt_CONST_REAL(nbg)
      ,CHFt_CONST_REAL(mte)
      ,CHFt_CONST_REAL(etaT) );
}
#endif  // GUARDJOULEHEATFLUID 

#ifndef GUARDELECTEMP 
#define GUARDELECTEMP 
// Prototype for Fortran procedure ELECTEMP ...
//
void FORTRAN_NAME( ELECTEMP ,electemp )(
      CHFp_FRA1(Tele)
      ,CHFp_CONST_FRA1(eden)
      ,CHFp_CONST_FRA1(enrg)
      ,CHFp_CONST_REAL(mte)
      ,CHFp_CONST_INT(cte)
      ,CHFp_BOX(cellBox) );

#define FORT_ELECTEMP FORTRAN_NAME( inlineELECTEMP, inlineELECTEMP)
#define FORTNT_ELECTEMP FORTRAN_NAME( ELECTEMP, electemp)

inline void FORTRAN_NAME(inlineELECTEMP, inlineELECTEMP)(
      CHFp_FRA1(Tele)
      ,CHFp_CONST_FRA1(eden)
      ,CHFp_CONST_FRA1(enrg)
      ,CHFp_CONST_REAL(mte)
      ,CHFp_CONST_INT(cte)
      ,CHFp_BOX(cellBox) )
{
 CH_TIMELEAF("FORT_ELECTEMP");
 FORTRAN_NAME( ELECTEMP ,electemp )(
      CHFt_FRA1(Tele)
      ,CHFt_CONST_FRA1(eden)
      ,CHFt_CONST_FRA1(enrg)
      ,CHFt_CONST_REAL(mte)
      ,CHFt_CONST_INT(cte)
      ,CHFt_BOX(cellBox) );
}
#endif  // GUARDELECTEMP 

#ifndef GUARDMASK 
#define GUARDMASK 
// Prototype for Fortran procedure MASK ...
//
void FORTRAN_NAME( MASK ,mask )(
      CHFp_FRA1(DMask)
      ,CHFp_FRA(src)
      ,CHFp_CONST_FRA1(eden)
      ,CHFp_CONST_FRA1(enrg)
      ,CHFp_BOX(cellBox) );

#define FORT_MASK FORTRAN_NAME( inlineMASK, inlineMASK)
#define FORTNT_MASK FORTRAN_NAME( MASK, mask)

inline void FORTRAN_NAME(inlineMASK, inlineMASK)(
      CHFp_FRA1(DMask)
      ,CHFp_FRA(src)
      ,CHFp_CONST_FRA1(eden)
      ,CHFp_CONST_FRA1(enrg)
      ,CHFp_BOX(cellBox) )
{
 CH_TIMELEAF("FORT_MASK");
 FORTRAN_NAME( MASK ,mask )(
      CHFt_FRA1(DMask)
      ,CHFt_FRA(src)
      ,CHFt_CONST_FRA1(eden)
      ,CHFt_CONST_FRA1(enrg)
      ,CHFt_BOX(cellBox) );
}
#endif  // GUARDMASK 

}

#endif
