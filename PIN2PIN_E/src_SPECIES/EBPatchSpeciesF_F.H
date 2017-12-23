#ifndef _EBPATCHSPECIESF_F_H_
#define _EBPATCHSPECIESF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSPECRIEMANN 
#define GUARDSPECRIEMANN 
// Prototype for Fortran procedure specriemann ...
//
void FORTRAN_NAME( SPECRIEMANN ,specriemann )(
      CHFp_BOX(dcalc)
      ,CHFp_FRA(primgdnv)
      ,CHFp_CONST_FRA(primleft)
      ,CHFp_CONST_FRA(primrigh)
      ,CHFp_CONST_FRA1(advectvel)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(nprim) );

#define FORT_SPECRIEMANN FORTRAN_NAME( inlineSPECRIEMANN, inlineSPECRIEMANN)
#define FORTNT_SPECRIEMANN FORTRAN_NAME( SPECRIEMANN, specriemann)

inline void FORTRAN_NAME(inlineSPECRIEMANN, inlineSPECRIEMANN)(
      CHFp_BOX(dcalc)
      ,CHFp_FRA(primgdnv)
      ,CHFp_CONST_FRA(primleft)
      ,CHFp_CONST_FRA(primrigh)
      ,CHFp_CONST_FRA1(advectvel)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(nprim) )
{
 CH_TIMELEAF("FORT_SPECRIEMANN");
 FORTRAN_NAME( SPECRIEMANN ,specriemann )(
      CHFt_BOX(dcalc)
      ,CHFt_FRA(primgdnv)
      ,CHFt_CONST_FRA(primleft)
      ,CHFt_CONST_FRA(primrigh)
      ,CHFt_CONST_FRA1(advectvel)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_CONST_INT(nprim) );
}
#endif  // GUARDSPECRIEMANN 

#ifndef GUARDGETFLUXSPEC 
#define GUARDGETFLUXSPEC 
// Prototype for Fortran procedure GETFLUXSPEC ...
//
void FORTRAN_NAME( GETFLUXSPEC ,getfluxspec )(
      CHFp_BOX(dcalc)
      ,CHFp_CONST_FRA(primgdnv)
      ,CHFp_CONST_FRA1(advectvel)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_FRA(flux)
      ,CHFp_CONST_INT(nprim) );

#define FORT_GETFLUXSPEC FORTRAN_NAME( inlineGETFLUXSPEC, inlineGETFLUXSPEC)
#define FORTNT_GETFLUXSPEC FORTRAN_NAME( GETFLUXSPEC, getfluxspec)

inline void FORTRAN_NAME(inlineGETFLUXSPEC, inlineGETFLUXSPEC)(
      CHFp_BOX(dcalc)
      ,CHFp_CONST_FRA(primgdnv)
      ,CHFp_CONST_FRA1(advectvel)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_FRA(flux)
      ,CHFp_CONST_INT(nprim) )
{
 CH_TIMELEAF("FORT_GETFLUXSPEC");
 FORTRAN_NAME( GETFLUXSPEC ,getfluxspec )(
      CHFt_BOX(dcalc)
      ,CHFt_CONST_FRA(primgdnv)
      ,CHFt_CONST_FRA1(advectvel)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_FRA(flux)
      ,CHFt_CONST_INT(nprim) );
}
#endif  // GUARDGETFLUXSPEC 

#ifndef GUARDFLOORSPECIES 
#define GUARDFLOORSPECIES 
// Prototype for Fortran procedure FLOORSPECIES ...
//
void FORTRAN_NAME( FLOORSPECIES ,floorspecies )(
      CHFp_BOX(dcalc)
      ,CHFp_FRA(state)
      ,CHFp_CONST_REAL(mx)
      ,CHFp_CONST_REAL(mn) );

#define FORT_FLOORSPECIES FORTRAN_NAME( inlineFLOORSPECIES, inlineFLOORSPECIES)
#define FORTNT_FLOORSPECIES FORTRAN_NAME( FLOORSPECIES, floorspecies)

inline void FORTRAN_NAME(inlineFLOORSPECIES, inlineFLOORSPECIES)(
      CHFp_BOX(dcalc)
      ,CHFp_FRA(state)
      ,CHFp_CONST_REAL(mx)
      ,CHFp_CONST_REAL(mn) )
{
 CH_TIMELEAF("FORT_FLOORSPECIES");
 FORTRAN_NAME( FLOORSPECIES ,floorspecies )(
      CHFt_BOX(dcalc)
      ,CHFt_FRA(state)
      ,CHFt_CONST_REAL(mx)
      ,CHFt_CONST_REAL(mn) );
}
#endif  // GUARDFLOORSPECIES 

#ifndef GUARDSCALARFABMULT 
#define GUARDSCALARFABMULT 
// Prototype for Fortran procedure SCALARFABMULT ...
//
void FORTRAN_NAME( SCALARFABMULT ,scalarfabmult )(
      CHFp_FRA1(flux)
      ,CHFp_CONST_FRA1(advectvel)
      ,CHFp_CONST_REAL(value)
      ,CHFp_BOX(box) );

#define FORT_SCALARFABMULT FORTRAN_NAME( inlineSCALARFABMULT, inlineSCALARFABMULT)
#define FORTNT_SCALARFABMULT FORTRAN_NAME( SCALARFABMULT, scalarfabmult)

inline void FORTRAN_NAME(inlineSCALARFABMULT, inlineSCALARFABMULT)(
      CHFp_FRA1(flux)
      ,CHFp_CONST_FRA1(advectvel)
      ,CHFp_CONST_REAL(value)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SCALARFABMULT");
 FORTRAN_NAME( SCALARFABMULT ,scalarfabmult )(
      CHFt_FRA1(flux)
      ,CHFt_CONST_FRA1(advectvel)
      ,CHFt_CONST_REAL(value)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSCALARFABMULT 

#ifndef GUARDSPECIESEXTRAPBC 
#define GUARDSPECIESEXTRAPBC 
// Prototype for Fortran procedure SPECIESEXTRAPBC ...
//
void FORTRAN_NAME( SPECIESEXTRAPBC ,speciesextrapbc )(
      CHFp_FRA(flux)
      ,CHFp_CONST_FRA(wextrap)
      ,CHFp_CONST_FRA1(advectvel)
      ,CHFp_CONST_INT(nspec)
      ,CHFp_BOX(box) );

#define FORT_SPECIESEXTRAPBC FORTRAN_NAME( inlineSPECIESEXTRAPBC, inlineSPECIESEXTRAPBC)
#define FORTNT_SPECIESEXTRAPBC FORTRAN_NAME( SPECIESEXTRAPBC, speciesextrapbc)

inline void FORTRAN_NAME(inlineSPECIESEXTRAPBC, inlineSPECIESEXTRAPBC)(
      CHFp_FRA(flux)
      ,CHFp_CONST_FRA(wextrap)
      ,CHFp_CONST_FRA1(advectvel)
      ,CHFp_CONST_INT(nspec)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SPECIESEXTRAPBC");
 FORTRAN_NAME( SPECIESEXTRAPBC ,speciesextrapbc )(
      CHFt_FRA(flux)
      ,CHFt_CONST_FRA(wextrap)
      ,CHFt_CONST_FRA1(advectvel)
      ,CHFt_CONST_INT(nspec)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSPECIESEXTRAPBC 

}

#endif
