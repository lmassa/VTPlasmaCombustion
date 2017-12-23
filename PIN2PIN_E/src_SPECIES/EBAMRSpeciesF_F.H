#ifndef _EBAMRSPECIESF_F_H_
#define _EBAMRSPECIESF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDCELLGRADEBSPEC 
#define GUARDCELLGRADEBSPEC 
// Prototype for Fortran procedure cellgradebspec ...
//
void FORTRAN_NAME( CELLGRADEBSPEC ,cellgradebspec )(
      CHFp_FRA1(grad)
      ,CHFp_CONST_FRA1(vel)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(lobox)
      ,CHFp_BOX(hibox)
      ,CHFp_BOX(centerbox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_CONST_INT(hashi)
      ,CHFp_CONST_INT(divdir) );

#define FORT_CELLGRADEBSPEC FORTRAN_NAME( inlineCELLGRADEBSPEC, inlineCELLGRADEBSPEC)
#define FORTNT_CELLGRADEBSPEC FORTRAN_NAME( CELLGRADEBSPEC, cellgradebspec)

inline void FORTRAN_NAME(inlineCELLGRADEBSPEC, inlineCELLGRADEBSPEC)(
      CHFp_FRA1(grad)
      ,CHFp_CONST_FRA1(vel)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(lobox)
      ,CHFp_BOX(hibox)
      ,CHFp_BOX(centerbox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_CONST_INT(hashi)
      ,CHFp_CONST_INT(divdir) )
{
 CH_TIMELEAF("FORT_CELLGRADEBSPEC");
 FORTRAN_NAME( CELLGRADEBSPEC ,cellgradebspec )(
      CHFt_FRA1(grad)
      ,CHFt_CONST_FRA1(vel)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(lobox)
      ,CHFt_BOX(hibox)
      ,CHFt_BOX(centerbox)
      ,CHFt_CONST_INT(haslo)
      ,CHFt_CONST_INT(hashi)
      ,CHFt_CONST_INT(divdir) );
}
#endif  // GUARDCELLGRADEBSPEC 

#ifndef GUARDSHARFETTERGUMMEL 
#define GUARDSHARFETTERGUMMEL 
// Prototype for Fortran procedure SHARFETTERGUMMEL ...
//
void FORTRAN_NAME( SHARFETTERGUMMEL ,sharfettergummel )(
      CHFp_FRA(MobData)
      ,CHFp_FRA(DiffData)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(nspec) );

#define FORT_SHARFETTERGUMMEL FORTRAN_NAME( inlineSHARFETTERGUMMEL, inlineSHARFETTERGUMMEL)
#define FORTNT_SHARFETTERGUMMEL FORTRAN_NAME( SHARFETTERGUMMEL, sharfettergummel)

inline void FORTRAN_NAME(inlineSHARFETTERGUMMEL, inlineSHARFETTERGUMMEL)(
      CHFp_FRA(MobData)
      ,CHFp_FRA(DiffData)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(nspec) )
{
 CH_TIMELEAF("FORT_SHARFETTERGUMMEL");
 FORTRAN_NAME( SHARFETTERGUMMEL ,sharfettergummel )(
      CHFt_FRA(MobData)
      ,CHFt_FRA(DiffData)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(nspec) );
}
#endif  // GUARDSHARFETTERGUMMEL 

#ifndef GUARDSHARFETTERGUMMEL_FD 
#define GUARDSHARFETTERGUMMEL_FD 
// Prototype for Fortran procedure SHARFETTERGUMMEL_FD ...
//
void FORTRAN_NAME( SHARFETTERGUMMEL_FD ,sharfettergummel_fd )(
      CHFp_FRA(MobData)
      ,CHFp_FRA(DiffData)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(nspec) );

#define FORT_SHARFETTERGUMMEL_FD FORTRAN_NAME( inlineSHARFETTERGUMMEL_FD, inlineSHARFETTERGUMMEL_FD)
#define FORTNT_SHARFETTERGUMMEL_FD FORTRAN_NAME( SHARFETTERGUMMEL_FD, sharfettergummel_fd)

inline void FORTRAN_NAME(inlineSHARFETTERGUMMEL_FD, inlineSHARFETTERGUMMEL_FD)(
      CHFp_FRA(MobData)
      ,CHFp_FRA(DiffData)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(nspec) )
{
 CH_TIMELEAF("FORT_SHARFETTERGUMMEL_FD");
 FORTRAN_NAME( SHARFETTERGUMMEL_FD ,sharfettergummel_fd )(
      CHFt_FRA(MobData)
      ,CHFt_FRA(DiffData)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(nspec) );
}
#endif  // GUARDSHARFETTERGUMMEL_FD 

}

#endif
