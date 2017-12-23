#ifndef _EBSPECIESFLUXOPF_F_H_
#define _EBSPECIESFLUXOPF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSPECIESFLUXGSRB 
#define GUARDSPECIESFLUXGSRB 
// Prototype for Fortran procedure speciesfluxgsrb ...
//
void FORTRAN_NAME( SPECIESFLUXGSRB ,speciesfluxgsrb )(
      CHFp_FRA1(phi)
      ,CHFp_CONST_FRA1(rhs)
      ,CHFp_CONST_FRA1(relcoef)
      ,CHFp_CONST_FRA1(acoef)
      ,CHFp_CONST_FRA1(b0)
      ,CHFp_CONST_FRA1(b1)
      ,CHFp_CONST_FRA1(b2)
      ,CHFp_CONST_FRA1(W0)
      ,CHFp_CONST_FRA1(W1)
      ,CHFp_CONST_FRA1(W2)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(region)
      ,CHFp_CONST_INT(redBlack) );

#define FORT_SPECIESFLUXGSRB FORTRAN_NAME( inlineSPECIESFLUXGSRB, inlineSPECIESFLUXGSRB)
#define FORTNT_SPECIESFLUXGSRB FORTRAN_NAME( SPECIESFLUXGSRB, speciesfluxgsrb)

inline void FORTRAN_NAME(inlineSPECIESFLUXGSRB, inlineSPECIESFLUXGSRB)(
      CHFp_FRA1(phi)
      ,CHFp_CONST_FRA1(rhs)
      ,CHFp_CONST_FRA1(relcoef)
      ,CHFp_CONST_FRA1(acoef)
      ,CHFp_CONST_FRA1(b0)
      ,CHFp_CONST_FRA1(b1)
      ,CHFp_CONST_FRA1(b2)
      ,CHFp_CONST_FRA1(W0)
      ,CHFp_CONST_FRA1(W1)
      ,CHFp_CONST_FRA1(W2)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(region)
      ,CHFp_CONST_INT(redBlack) )
{
 CH_TIMELEAF("FORT_SPECIESFLUXGSRB");
 FORTRAN_NAME( SPECIESFLUXGSRB ,speciesfluxgsrb )(
      CHFt_FRA1(phi)
      ,CHFt_CONST_FRA1(rhs)
      ,CHFt_CONST_FRA1(relcoef)
      ,CHFt_CONST_FRA1(acoef)
      ,CHFt_CONST_FRA1(b0)
      ,CHFt_CONST_FRA1(b1)
      ,CHFt_CONST_FRA1(b2)
      ,CHFt_CONST_FRA1(W0)
      ,CHFt_CONST_FRA1(W1)
      ,CHFt_CONST_FRA1(W2)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(region)
      ,CHFt_CONST_INT(redBlack) );
}
#endif  // GUARDSPECIESFLUXGSRB 

#ifndef GUARDEBSFOREGAPPLYDOMAINFLUX 
#define GUARDEBSFOREGAPPLYDOMAINFLUX 
// Prototype for Fortran procedure ebsforegapplydomainflux ...
//
void FORTRAN_NAME( EBSFOREGAPPLYDOMAINFLUX ,ebsforegapplydomainflux )(
      CHFp_FRA1(phi)
      ,CHFp_CONST_FRA1(faceflux)
      ,CHFp_CONST_FRA1(bc)
      ,CHFp_CONST_FRA1(Wc)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(side)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_EBSFOREGAPPLYDOMAINFLUX FORTRAN_NAME( inlineEBSFOREGAPPLYDOMAINFLUX, inlineEBSFOREGAPPLYDOMAINFLUX)
#define FORTNT_EBSFOREGAPPLYDOMAINFLUX FORTRAN_NAME( EBSFOREGAPPLYDOMAINFLUX, ebsforegapplydomainflux)

inline void FORTRAN_NAME(inlineEBSFOREGAPPLYDOMAINFLUX, inlineEBSFOREGAPPLYDOMAINFLUX)(
      CHFp_FRA1(phi)
      ,CHFp_CONST_FRA1(faceflux)
      ,CHFp_CONST_FRA1(bc)
      ,CHFp_CONST_FRA1(Wc)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(side)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_EBSFOREGAPPLYDOMAINFLUX");
 FORTRAN_NAME( EBSFOREGAPPLYDOMAINFLUX ,ebsforegapplydomainflux )(
      CHFt_FRA1(phi)
      ,CHFt_CONST_FRA1(faceflux)
      ,CHFt_CONST_FRA1(bc)
      ,CHFt_CONST_FRA1(Wc)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(side)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDEBSFOREGAPPLYDOMAINFLUX 

#ifndef GUARDSPECIESFLUXINPLACE 
#define GUARDSPECIESFLUXINPLACE 
// Prototype for Fortran procedure speciesfluxinplace ...
//
void FORTRAN_NAME( SPECIESFLUXINPLACE ,speciesfluxinplace )(
      CHFp_FRA1(opphidir)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_FRA1(b0)
      ,CHFp_CONST_FRA1(b1)
      ,CHFp_CONST_FRA1(b2)
      ,CHFp_CONST_FRA1(W0)
      ,CHFp_CONST_FRA1(W1)
      ,CHFp_CONST_FRA1(W2)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_SPECIESFLUXINPLACE FORTRAN_NAME( inlineSPECIESFLUXINPLACE, inlineSPECIESFLUXINPLACE)
#define FORTNT_SPECIESFLUXINPLACE FORTRAN_NAME( SPECIESFLUXINPLACE, speciesfluxinplace)

inline void FORTRAN_NAME(inlineSPECIESFLUXINPLACE, inlineSPECIESFLUXINPLACE)(
      CHFp_FRA1(opphidir)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_FRA1(b0)
      ,CHFp_CONST_FRA1(b1)
      ,CHFp_CONST_FRA1(b2)
      ,CHFp_CONST_FRA1(W0)
      ,CHFp_CONST_FRA1(W1)
      ,CHFp_CONST_FRA1(W2)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SPECIESFLUXINPLACE");
 FORTRAN_NAME( SPECIESFLUXINPLACE ,speciesfluxinplace )(
      CHFt_FRA1(opphidir)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_CONST_FRA1(b0)
      ,CHFt_CONST_FRA1(b1)
      ,CHFt_CONST_FRA1(b2)
      ,CHFt_CONST_FRA1(W0)
      ,CHFt_CONST_FRA1(W1)
      ,CHFt_CONST_FRA1(W2)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSPECIESFLUXINPLACE 

#ifndef GUARDINCRAPPLYEBSFO 
#define GUARDINCRAPPLYEBSFO 
// Prototype for Fortran procedure incrapplyebsfo ...
//
void FORTRAN_NAME( INCRAPPLYEBSFO ,incrapplyebsfo )(
      CHFp_FRA1(lhs)
      ,CHFp_CONST_FRA1(interiorflux)
      ,CHFp_CONST_FRA1(domainfluxlo)
      ,CHFp_CONST_FRA1(domainfluxhi)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(lobox)
      ,CHFp_BOX(hibox)
      ,CHFp_BOX(centerbox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_CONST_INT(hashi)
      ,CHFp_CONST_INT(facedir) );

#define FORT_INCRAPPLYEBSFO FORTRAN_NAME( inlineINCRAPPLYEBSFO, inlineINCRAPPLYEBSFO)
#define FORTNT_INCRAPPLYEBSFO FORTRAN_NAME( INCRAPPLYEBSFO, incrapplyebsfo)

inline void FORTRAN_NAME(inlineINCRAPPLYEBSFO, inlineINCRAPPLYEBSFO)(
      CHFp_FRA1(lhs)
      ,CHFp_CONST_FRA1(interiorflux)
      ,CHFp_CONST_FRA1(domainfluxlo)
      ,CHFp_CONST_FRA1(domainfluxhi)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(lobox)
      ,CHFp_BOX(hibox)
      ,CHFp_BOX(centerbox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_CONST_INT(hashi)
      ,CHFp_CONST_INT(facedir) )
{
 CH_TIMELEAF("FORT_INCRAPPLYEBSFO");
 FORTRAN_NAME( INCRAPPLYEBSFO ,incrapplyebsfo )(
      CHFt_FRA1(lhs)
      ,CHFt_CONST_FRA1(interiorflux)
      ,CHFt_CONST_FRA1(domainfluxlo)
      ,CHFt_CONST_FRA1(domainfluxhi)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(lobox)
      ,CHFt_BOX(hibox)
      ,CHFt_BOX(centerbox)
      ,CHFt_CONST_INT(haslo)
      ,CHFt_CONST_INT(hashi)
      ,CHFt_CONST_INT(facedir) );
}
#endif  // GUARDINCRAPPLYEBSFO 

#ifndef GUARDDECRINVRELCOEFEBSFO 
#define GUARDDECRINVRELCOEFEBSFO 
// Prototype for Fortran procedure decrinvrelcoefebsfo ...
//
void FORTRAN_NAME( DECRINVRELCOEFEBSFO ,decrinvrelcoefebsfo )(
      CHFp_FRA1(relcoef)
      ,CHFp_FRA1(Dk)
      ,CHFp_FRA1(VW)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_BOX(box)
      ,CHFp_REAL(dx)
      ,CHFp_INT(idir) );

#define FORT_DECRINVRELCOEFEBSFO FORTRAN_NAME( inlineDECRINVRELCOEFEBSFO, inlineDECRINVRELCOEFEBSFO)
#define FORTNT_DECRINVRELCOEFEBSFO FORTRAN_NAME( DECRINVRELCOEFEBSFO, decrinvrelcoefebsfo)

inline void FORTRAN_NAME(inlineDECRINVRELCOEFEBSFO, inlineDECRINVRELCOEFEBSFO)(
      CHFp_FRA1(relcoef)
      ,CHFp_FRA1(Dk)
      ,CHFp_FRA1(VW)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_BOX(box)
      ,CHFp_REAL(dx)
      ,CHFp_INT(idir) )
{
 CH_TIMELEAF("FORT_DECRINVRELCOEFEBSFO");
 FORTRAN_NAME( DECRINVRELCOEFEBSFO ,decrinvrelcoefebsfo )(
      CHFt_FRA1(relcoef)
      ,CHFt_FRA1(Dk)
      ,CHFt_FRA1(VW)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_BOX(box)
      ,CHFt_REAL(dx)
      ,CHFt_INT(idir) );
}
#endif  // GUARDDECRINVRELCOEFEBSFO 

#ifndef GUARDGETFLUXEBSFO 
#define GUARDGETFLUXEBSFO 
// Prototype for Fortran procedure getfluxebsfo ...
//
void FORTRAN_NAME( GETFLUXEBSFO ,getfluxebsfo )(
      CHFp_FRA1(flux)
      ,CHFp_CONST_FRA1(bcoef)
      ,CHFp_CONST_FRA1(Wcoef)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_BOX(opphibox)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_INT(idir) );

#define FORT_GETFLUXEBSFO FORTRAN_NAME( inlineGETFLUXEBSFO, inlineGETFLUXEBSFO)
#define FORTNT_GETFLUXEBSFO FORTRAN_NAME( GETFLUXEBSFO, getfluxebsfo)

inline void FORTRAN_NAME(inlineGETFLUXEBSFO, inlineGETFLUXEBSFO)(
      CHFp_FRA1(flux)
      ,CHFp_CONST_FRA1(bcoef)
      ,CHFp_CONST_FRA1(Wcoef)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_BOX(opphibox)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_INT(idir) )
{
 CH_TIMELEAF("FORT_GETFLUXEBSFO");
 FORTRAN_NAME( GETFLUXEBSFO ,getfluxebsfo )(
      CHFt_FRA1(flux)
      ,CHFt_CONST_FRA1(bcoef)
      ,CHFt_CONST_FRA1(Wcoef)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_BOX(opphibox)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_INT(idir) );
}
#endif  // GUARDGETFLUXEBSFO 

#ifndef GUARDGSRBEBSFO 
#define GUARDGSRBEBSFO 
// Prototype for Fortran procedure gsrbebsfo ...
//
void FORTRAN_NAME( GSRBEBSFO ,gsrbebsfo )(
      CHFp_FRA1(phi)
      ,CHFp_CONST_FRA1(lph)
      ,CHFp_CONST_FRA1(rhs)
      ,CHFp_CONST_FRA1(lam)
      ,CHFp_BOX(coloredbox) );

#define FORT_GSRBEBSFO FORTRAN_NAME( inlineGSRBEBSFO, inlineGSRBEBSFO)
#define FORTNT_GSRBEBSFO FORTRAN_NAME( GSRBEBSFO, gsrbebsfo)

inline void FORTRAN_NAME(inlineGSRBEBSFO, inlineGSRBEBSFO)(
      CHFp_FRA1(phi)
      ,CHFp_CONST_FRA1(lph)
      ,CHFp_CONST_FRA1(rhs)
      ,CHFp_CONST_FRA1(lam)
      ,CHFp_BOX(coloredbox) )
{
 CH_TIMELEAF("FORT_GSRBEBSFO");
 FORTRAN_NAME( GSRBEBSFO ,gsrbebsfo )(
      CHFt_FRA1(phi)
      ,CHFt_CONST_FRA1(lph)
      ,CHFt_CONST_FRA1(rhs)
      ,CHFt_CONST_FRA1(lam)
      ,CHFt_BOX(coloredbox) );
}
#endif  // GUARDGSRBEBSFO 

#ifndef GUARDINCRAVEFACETOCELL 
#define GUARDINCRAVEFACETOCELL 
// Prototype for Fortran procedure incravefacetocell ...
//
void FORTRAN_NAME( INCRAVEFACETOCELL ,incravefacetocell )(
      CHFp_FRA1(cellgrad)
      ,CHFp_CONST_FRA1(facegrad)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(cellbox) );

#define FORT_INCRAVEFACETOCELL FORTRAN_NAME( inlineINCRAVEFACETOCELL, inlineINCRAVEFACETOCELL)
#define FORTNT_INCRAVEFACETOCELL FORTRAN_NAME( INCRAVEFACETOCELL, incravefacetocell)

inline void FORTRAN_NAME(inlineINCRAVEFACETOCELL, inlineINCRAVEFACETOCELL)(
      CHFp_FRA1(cellgrad)
      ,CHFp_CONST_FRA1(facegrad)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(cellbox) )
{
 CH_TIMELEAF("FORT_INCRAVEFACETOCELL");
 FORTRAN_NAME( INCRAVEFACETOCELL ,incravefacetocell )(
      CHFt_FRA1(cellgrad)
      ,CHFt_CONST_FRA1(facegrad)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_BOX(cellbox) );
}
#endif  // GUARDINCRAVEFACETOCELL 

#ifndef GUARDSG1FD 
#define GUARDSG1FD 
// Prototype for Fortran procedure SG1FD ...
//
void FORTRAN_NAME( SG1FD ,sg1fd )(
      CHFp_FRA1(MobData)
      ,CHFp_FRA1(DiffData)
      ,CHFp_CONST_FRA1(AMobData)
      ,CHFp_CONST_FRA1(ADiffData)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) );

#define FORT_SG1FD FORTRAN_NAME( inlineSG1FD, inlineSG1FD)
#define FORTNT_SG1FD FORTRAN_NAME( SG1FD, sg1fd)

inline void FORTRAN_NAME(inlineSG1FD, inlineSG1FD)(
      CHFp_FRA1(MobData)
      ,CHFp_FRA1(DiffData)
      ,CHFp_CONST_FRA1(AMobData)
      ,CHFp_CONST_FRA1(ADiffData)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_SG1FD");
 FORTRAN_NAME( SG1FD ,sg1fd )(
      CHFt_FRA1(MobData)
      ,CHFt_FRA1(DiffData)
      ,CHFt_CONST_FRA1(AMobData)
      ,CHFt_CONST_FRA1(ADiffData)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDSG1FD 

#ifndef GUARDSG1 
#define GUARDSG1 
// Prototype for Fortran procedure SG1 ...
//
void FORTRAN_NAME( SG1 ,sg1 )(
      CHFp_FRA1(MobData)
      ,CHFp_FRA1(DiffData)
      ,CHFp_CONST_FRA1(AMobData)
      ,CHFp_CONST_FRA1(ADiffData)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) );

#define FORT_SG1 FORTRAN_NAME( inlineSG1, inlineSG1)
#define FORTNT_SG1 FORTRAN_NAME( SG1, sg1)

inline void FORTRAN_NAME(inlineSG1, inlineSG1)(
      CHFp_FRA1(MobData)
      ,CHFp_FRA1(DiffData)
      ,CHFp_CONST_FRA1(AMobData)
      ,CHFp_CONST_FRA1(ADiffData)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_SG1");
 FORTRAN_NAME( SG1 ,sg1 )(
      CHFt_FRA1(MobData)
      ,CHFt_FRA1(DiffData)
      ,CHFt_CONST_FRA1(AMobData)
      ,CHFt_CONST_FRA1(ADiffData)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDSG1 

}

#endif
