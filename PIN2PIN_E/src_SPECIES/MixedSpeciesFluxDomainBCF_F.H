#ifndef _MIXEDSPECIESFLUXDOMAINBCF_F_H_
#define _MIXEDSPECIESFLUXDOMAINBCF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSETMIXEDFACEFLUX 
#define GUARDSETMIXEDFACEFLUX 
// Prototype for Fortran procedure SETMIXEDFACEFLUX ...
//
void FORTRAN_NAME( SETMIXEDFACEFLUX ,setmixedfaceflux )(
      CHFp_FRA1(faceFlux)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_FRA1(regD)
      ,CHFp_CONST_FRA1(regW)
      ,CHFp_CONST_REAL(value)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(iside)
      ,CHFp_BOX(box) );

#define FORT_SETMIXEDFACEFLUX FORTRAN_NAME( inlineSETMIXEDFACEFLUX, inlineSETMIXEDFACEFLUX)
#define FORTNT_SETMIXEDFACEFLUX FORTRAN_NAME( SETMIXEDFACEFLUX, setmixedfaceflux)

inline void FORTRAN_NAME(inlineSETMIXEDFACEFLUX, inlineSETMIXEDFACEFLUX)(
      CHFp_FRA1(faceFlux)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_FRA1(regD)
      ,CHFp_CONST_FRA1(regW)
      ,CHFp_CONST_REAL(value)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(iside)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SETMIXEDFACEFLUX");
 FORTRAN_NAME( SETMIXEDFACEFLUX ,setmixedfaceflux )(
      CHFt_FRA1(faceFlux)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_CONST_FRA1(regD)
      ,CHFt_CONST_FRA1(regW)
      ,CHFt_CONST_REAL(value)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_INT(iside)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSETMIXEDFACEFLUX 

}

#endif
