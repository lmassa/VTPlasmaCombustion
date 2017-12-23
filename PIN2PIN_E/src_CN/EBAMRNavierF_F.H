#ifndef _EBAMRNAVIERF_F_H_
#define _EBAMRNAVIERF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDLIMITBAR 
#define GUARDLIMITBAR 
// Prototype for Fortran procedure limitbar ...
//
void FORTRAN_NAME( LIMITBAR ,limitbar )(
      CHFp_FRA(qbar)
      ,CHFp_CONST_FRA(q)
      ,CHFp_CONST_REAL(Dlimit)
      ,CHFp_BOX(box) );

#define FORT_LIMITBAR FORTRAN_NAME( inlineLIMITBAR, inlineLIMITBAR)
#define FORTNT_LIMITBAR FORTRAN_NAME( LIMITBAR, limitbar)

inline void FORTRAN_NAME(inlineLIMITBAR, inlineLIMITBAR)(
      CHFp_FRA(qbar)
      ,CHFp_CONST_FRA(q)
      ,CHFp_CONST_REAL(Dlimit)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_LIMITBAR");
 FORTRAN_NAME( LIMITBAR ,limitbar )(
      CHFt_FRA(qbar)
      ,CHFt_CONST_FRA(q)
      ,CHFt_CONST_REAL(Dlimit)
      ,CHFt_BOX(box) );
}
#endif  // GUARDLIMITBAR 

}

#endif
