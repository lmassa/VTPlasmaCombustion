#ifndef _EBSOLIDF_F_H_
#define _EBSOLIDF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSOLIDBC 
#define GUARDSOLIDBC 
// Prototype for Fortran procedure SOLIDBC ...
//
void FORTRAN_NAME( SOLIDBC ,solidbc )(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(Wextrap)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(box) );

#define FORT_SOLIDBC FORTRAN_NAME( inlineSOLIDBC, inlineSOLIDBC)
#define FORTNT_SOLIDBC FORTRAN_NAME( SOLIDBC, solidbc)

inline void FORTRAN_NAME(inlineSOLIDBC, inlineSOLIDBC)(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(Wextrap)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SOLIDBC");
 FORTRAN_NAME( SOLIDBC ,solidbc )(
      CHFt_FRA(F)
      ,CHFt_CONST_FRA(Wextrap)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(dir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSOLIDBC 

#ifndef GUARDSOLIDBCRZ 
#define GUARDSOLIDBCRZ 
// Prototype for Fortran procedure SOLIDBCRZ ...
//
void FORTRAN_NAME( SOLIDBCRZ ,solidbcrz )(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(Wextrap)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(box) );

#define FORT_SOLIDBCRZ FORTRAN_NAME( inlineSOLIDBCRZ, inlineSOLIDBCRZ)
#define FORTNT_SOLIDBCRZ FORTRAN_NAME( SOLIDBCRZ, solidbcrz)

inline void FORTRAN_NAME(inlineSOLIDBCRZ, inlineSOLIDBCRZ)(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(Wextrap)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SOLIDBCRZ");
 FORTRAN_NAME( SOLIDBCRZ ,solidbcrz )(
      CHFt_FRA(F)
      ,CHFt_CONST_FRA(Wextrap)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(dir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSOLIDBCRZ 

#ifndef GUARDSLOPEBCS 
#define GUARDSLOPEBCS 
// Prototype for Fortran procedure SLOPEBCS ...
//
void FORTRAN_NAME( SLOPEBCS ,slopebcs )(
      CHFp_FRA(dW)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi) );

#define FORT_SLOPEBCS FORTRAN_NAME( inlineSLOPEBCS, inlineSLOPEBCS)
#define FORTNT_SLOPEBCS FORTRAN_NAME( SLOPEBCS, slopebcs)

inline void FORTRAN_NAME(inlineSLOPEBCS, inlineSLOPEBCS)(
      CHFp_FRA(dW)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi) )
{
 CH_TIMELEAF("FORT_SLOPEBCS");
 FORTRAN_NAME( SLOPEBCS ,slopebcs )(
      CHFt_FRA(dW)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(dir)
      ,CHFt_BOX(loBox)
      ,CHFt_CONST_INT(hasLo)
      ,CHFt_BOX(hiBox)
      ,CHFt_CONST_INT(hasHi) );
}
#endif  // GUARDSLOPEBCS 

#ifndef GUARDEXTRAPBC 
#define GUARDEXTRAPBC 
// Prototype for Fortran procedure extrapbc ...
//
void FORTRAN_NAME( EXTRAPBC ,extrapbc )(
      CHFp_FRA(flux)
      ,CHFp_CONST_FRA(wextrap)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(box) );

#define FORT_EXTRAPBC FORTRAN_NAME( inlineEXTRAPBC, inlineEXTRAPBC)
#define FORTNT_EXTRAPBC FORTRAN_NAME( EXTRAPBC, extrapbc)

inline void FORTRAN_NAME(inlineEXTRAPBC, inlineEXTRAPBC)(
      CHFp_FRA(flux)
      ,CHFp_CONST_FRA(wextrap)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_EXTRAPBC");
 FORTRAN_NAME( EXTRAPBC ,extrapbc )(
      CHFt_FRA(flux)
      ,CHFt_CONST_FRA(wextrap)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDEXTRAPBC 

#ifndef GUARDOUTFLOW 
#define GUARDOUTFLOW 
// Prototype for Fortran procedure outflow ...
//
void FORTRAN_NAME( OUTFLOW ,outflow )(
      CHFp_FRA(flux)
      ,CHFp_CONST_FRA(wextrap)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(box) );

#define FORT_OUTFLOW FORTRAN_NAME( inlineOUTFLOW, inlineOUTFLOW)
#define FORTNT_OUTFLOW FORTRAN_NAME( OUTFLOW, outflow)

inline void FORTRAN_NAME(inlineOUTFLOW, inlineOUTFLOW)(
      CHFp_FRA(flux)
      ,CHFp_CONST_FRA(wextrap)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_OUTFLOW");
 FORTRAN_NAME( OUTFLOW ,outflow )(
      CHFt_FRA(flux)
      ,CHFt_CONST_FRA(wextrap)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDOUTFLOW 

#ifndef GUARDDIRICHLETBC 
#define GUARDDIRICHLETBC 
// Prototype for Fortran procedure dirichletbc ...
//
void FORTRAN_NAME( DIRICHLETBC ,dirichletbc )(
      CHFp_FRA(flux)
      ,CHFp_CONST_FRA(wextrap)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(box) );

#define FORT_DIRICHLETBC FORTRAN_NAME( inlineDIRICHLETBC, inlineDIRICHLETBC)
#define FORTNT_DIRICHLETBC FORTRAN_NAME( DIRICHLETBC, dirichletbc)

inline void FORTRAN_NAME(inlineDIRICHLETBC, inlineDIRICHLETBC)(
      CHFp_FRA(flux)
      ,CHFp_CONST_FRA(wextrap)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_DIRICHLETBC");
 FORTRAN_NAME( DIRICHLETBC ,dirichletbc )(
      CHFt_FRA(flux)
      ,CHFt_CONST_FRA(wextrap)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDDIRICHLETBC 

#ifndef GUARDFIXPRIMS 
#define GUARDFIXPRIMS 
// Prototype for Fortran procedure fixprims ...
//
void FORTRAN_NAME( FIXPRIMS ,fixprims )(
      CHFp_VR(primitive) );

#define FORT_FIXPRIMS FORTRAN_NAME( inlineFIXPRIMS, inlineFIXPRIMS)
#define FORTNT_FIXPRIMS FORTRAN_NAME( FIXPRIMS, fixprims)

inline void FORTRAN_NAME(inlineFIXPRIMS, inlineFIXPRIMS)(
      CHFp_VR(primitive) )
{
 CH_TIMELEAF("FORT_FIXPRIMS");
 FORTRAN_NAME( FIXPRIMS ,fixprims )(
      CHFt_VR(primitive) );
}
#endif  // GUARDFIXPRIMS 

#ifndef GUARDEXTRAPBCRZ 
#define GUARDEXTRAPBCRZ 
// Prototype for Fortran procedure extrapbcrz ...
//
void FORTRAN_NAME( EXTRAPBCRZ ,extrapbcrz )(
      CHFp_FRA(flux)
      ,CHFp_CONST_FRA(wextrap)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(box) );

#define FORT_EXTRAPBCRZ FORTRAN_NAME( inlineEXTRAPBCRZ, inlineEXTRAPBCRZ)
#define FORTNT_EXTRAPBCRZ FORTRAN_NAME( EXTRAPBCRZ, extrapbcrz)

inline void FORTRAN_NAME(inlineEXTRAPBCRZ, inlineEXTRAPBCRZ)(
      CHFp_FRA(flux)
      ,CHFp_CONST_FRA(wextrap)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_EXTRAPBCRZ");
 FORTRAN_NAME( EXTRAPBCRZ ,extrapbcrz )(
      CHFt_FRA(flux)
      ,CHFt_CONST_FRA(wextrap)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDEXTRAPBCRZ 

}

#endif
