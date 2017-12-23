#ifndef _EBPIN2PINF_F_H_
#define _EBPIN2PINF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDGETPOSTSHOCKVEL 
#define GUARDGETPOSTSHOCKVEL 
// Prototype for Fortran procedure getpostshockvel ...
//
void FORTRAN_NAME( GETPOSTSHOCKVEL ,getpostshockvel )(
      CHFp_REAL(aupostshockvel) );

#define FORT_GETPOSTSHOCKVEL FORTRAN_NAME( inlineGETPOSTSHOCKVEL, inlineGETPOSTSHOCKVEL)
#define FORTNT_GETPOSTSHOCKVEL FORTRAN_NAME( GETPOSTSHOCKVEL, getpostshockvel)

inline void FORTRAN_NAME(inlineGETPOSTSHOCKVEL, inlineGETPOSTSHOCKVEL)(
      CHFp_REAL(aupostshockvel) )
{
 CH_TIMELEAF("FORT_GETPOSTSHOCKVEL");
 FORTRAN_NAME( GETPOSTSHOCKVEL ,getpostshockvel )(
      CHFt_REAL(aupostshockvel) );
}
#endif  // GUARDGETPOSTSHOCKVEL 

#ifndef GUARDGETSHOCKLOCX 
#define GUARDGETSHOCKLOCX 
// Prototype for Fortran procedure getshocklocx ...
//
void FORTRAN_NAME( GETSHOCKLOCX ,getshocklocx )(
      CHFp_CONST_REAL(aushockloc)
      ,CHFp_CONST_REAL(time) );

#define FORT_GETSHOCKLOCX FORTRAN_NAME( inlineGETSHOCKLOCX, inlineGETSHOCKLOCX)
#define FORTNT_GETSHOCKLOCX FORTRAN_NAME( GETSHOCKLOCX, getshocklocx)

inline void FORTRAN_NAME(inlineGETSHOCKLOCX, inlineGETSHOCKLOCX)(
      CHFp_CONST_REAL(aushockloc)
      ,CHFp_CONST_REAL(time) )
{
 CH_TIMELEAF("FORT_GETSHOCKLOCX");
 FORTRAN_NAME( GETSHOCKLOCX ,getshocklocx )(
      CHFt_CONST_REAL(aushockloc)
      ,CHFt_CONST_REAL(time) );
}
#endif  // GUARDGETSHOCKLOCX 

#ifndef GUARDSETPIN2PIN 
#define GUARDSETPIN2PIN 
// Prototype for Fortran procedure setpin2pin ...
//
void FORTRAN_NAME( SETPIN2PIN ,setpin2pin )(
      CHFp_CONST_REAL(augamma)
      ,CHFp_CONST_REAL(aums)
      ,CHFp_CONST_REAL(aucenter)
      ,CHFp_CONST_REAL(aupreshockpress)
      ,CHFp_CONST_REAL(aupreshockdense)
      ,CHFp_CONST_INT(aushocknorm)
      ,CHFp_CONST_INT(aushockback) );

#define FORT_SETPIN2PIN FORTRAN_NAME( inlineSETPIN2PIN, inlineSETPIN2PIN)
#define FORTNT_SETPIN2PIN FORTRAN_NAME( SETPIN2PIN, setpin2pin)

inline void FORTRAN_NAME(inlineSETPIN2PIN, inlineSETPIN2PIN)(
      CHFp_CONST_REAL(augamma)
      ,CHFp_CONST_REAL(aums)
      ,CHFp_CONST_REAL(aucenter)
      ,CHFp_CONST_REAL(aupreshockpress)
      ,CHFp_CONST_REAL(aupreshockdense)
      ,CHFp_CONST_INT(aushocknorm)
      ,CHFp_CONST_INT(aushockback) )
{
 CH_TIMELEAF("FORT_SETPIN2PIN");
 FORTRAN_NAME( SETPIN2PIN ,setpin2pin )(
      CHFt_CONST_REAL(augamma)
      ,CHFt_CONST_REAL(aums)
      ,CHFt_CONST_REAL(aucenter)
      ,CHFt_CONST_REAL(aupreshockpress)
      ,CHFt_CONST_REAL(aupreshockdense)
      ,CHFt_CONST_INT(aushocknorm)
      ,CHFt_CONST_INT(aushockback) );
}
#endif  // GUARDSETPIN2PIN 

#ifndef GUARDPIN2PININIT 
#define GUARDPIN2PININIT 
// Prototype for Fortran procedure pin2pininit ...
//
void FORTRAN_NAME( PIN2PININIT ,pin2pininit )(
      CHFp_FRA(u)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_PIN2PININIT FORTRAN_NAME( inlinePIN2PININIT, inlinePIN2PININIT)
#define FORTNT_PIN2PININIT FORTRAN_NAME( PIN2PININIT, pin2pininit)

inline void FORTRAN_NAME(inlinePIN2PININIT, inlinePIN2PININIT)(
      CHFp_FRA(u)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_PIN2PININIT");
 FORTRAN_NAME( PIN2PININIT ,pin2pininit )(
      CHFt_FRA(u)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDPIN2PININIT 

#ifndef GUARDPOINTPIN2PININIT 
#define GUARDPOINTPIN2PININIT 
// Prototype for Fortran procedure pointpin2pininit ...
//
void FORTRAN_NAME( POINTPIN2PININIT ,pointpin2pininit )(
      CHFp_REAL(density)
      ,CHFp_REALVECT(momentum)
      ,CHFp_REAL(energy)
      ,CHFp_CONST_INTVECT(iv)
      ,CHFp_CONST_REAL(dx) );

#define FORT_POINTPIN2PININIT FORTRAN_NAME( inlinePOINTPIN2PININIT, inlinePOINTPIN2PININIT)
#define FORTNT_POINTPIN2PININIT FORTRAN_NAME( POINTPIN2PININIT, pointpin2pininit)

inline void FORTRAN_NAME(inlinePOINTPIN2PININIT, inlinePOINTPIN2PININIT)(
      CHFp_REAL(density)
      ,CHFp_REALVECT(momentum)
      ,CHFp_REAL(energy)
      ,CHFp_CONST_INTVECT(iv)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_POINTPIN2PININIT");
 FORTRAN_NAME( POINTPIN2PININIT ,pointpin2pininit )(
      CHFt_REAL(density)
      ,CHFt_REALVECT(momentum)
      ,CHFt_REAL(energy)
      ,CHFt_CONST_INTVECT(iv)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDPOINTPIN2PININIT 

#ifndef GUARDBLINIT 
#define GUARDBLINIT 
// Prototype for Fortran procedure blinit ...
//
void FORTRAN_NAME( BLINIT ,blinit )(
      CHFp_FRA1(u)
      ,CHFp_CONST_INT(ivar)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(x0)
      ,CHFp_CONST_REAL(ReRoot2)
      ,CHFp_CONST_VR(y)
      ,CHFp_CONST_VR(xdata)
      ,CHFp_CONST_VR(ydata)
      ,CHFp_CONST_VR(B)
      ,CHFp_CONST_VR(C)
      ,CHFp_CONST_VR(D)
      ,CHFp_CONST_VR(E)
      ,CHFp_CONST_VR(F)
      ,CHFp_BOX(box) );

#define FORT_BLINIT FORTRAN_NAME( inlineBLINIT, inlineBLINIT)
#define FORTNT_BLINIT FORTRAN_NAME( BLINIT, blinit)

inline void FORTRAN_NAME(inlineBLINIT, inlineBLINIT)(
      CHFp_FRA1(u)
      ,CHFp_CONST_INT(ivar)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(x0)
      ,CHFp_CONST_REAL(ReRoot2)
      ,CHFp_CONST_VR(y)
      ,CHFp_CONST_VR(xdata)
      ,CHFp_CONST_VR(ydata)
      ,CHFp_CONST_VR(B)
      ,CHFp_CONST_VR(C)
      ,CHFp_CONST_VR(D)
      ,CHFp_CONST_VR(E)
      ,CHFp_CONST_VR(F)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_BLINIT");
 FORTRAN_NAME( BLINIT ,blinit )(
      CHFt_FRA1(u)
      ,CHFt_CONST_INT(ivar)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_REAL(x0)
      ,CHFt_CONST_REAL(ReRoot2)
      ,CHFt_CONST_VR(y)
      ,CHFt_CONST_VR(xdata)
      ,CHFt_CONST_VR(ydata)
      ,CHFt_CONST_VR(B)
      ,CHFt_CONST_VR(C)
      ,CHFt_CONST_VR(D)
      ,CHFt_CONST_VR(E)
      ,CHFt_CONST_VR(F)
      ,CHFt_BOX(box) );
}
#endif  // GUARDBLINIT 

#ifndef GUARDSPLINE5 
#define GUARDSPLINE5 
// Prototype for Fortran procedure spline5 ...
//
void FORTRAN_NAME( SPLINE5 ,spline5 )(
      CHFp_REAL(yeval)
      ,CHFp_CONST_REAL(xeval)
      ,CHFp_CONST_VR(xdata)
      ,CHFp_CONST_VR(ydata)
      ,CHFp_CONST_VR(b)
      ,CHFp_CONST_VR(c)
      ,CHFp_CONST_VR(d)
      ,CHFp_CONST_VR(e)
      ,CHFp_CONST_VR(f) );

#define FORT_SPLINE5 FORTRAN_NAME( inlineSPLINE5, inlineSPLINE5)
#define FORTNT_SPLINE5 FORTRAN_NAME( SPLINE5, spline5)

inline void FORTRAN_NAME(inlineSPLINE5, inlineSPLINE5)(
      CHFp_REAL(yeval)
      ,CHFp_CONST_REAL(xeval)
      ,CHFp_CONST_VR(xdata)
      ,CHFp_CONST_VR(ydata)
      ,CHFp_CONST_VR(b)
      ,CHFp_CONST_VR(c)
      ,CHFp_CONST_VR(d)
      ,CHFp_CONST_VR(e)
      ,CHFp_CONST_VR(f) )
{
 CH_TIMELEAF("FORT_SPLINE5");
 FORTRAN_NAME( SPLINE5 ,spline5 )(
      CHFt_REAL(yeval)
      ,CHFt_CONST_REAL(xeval)
      ,CHFt_CONST_VR(xdata)
      ,CHFt_CONST_VR(ydata)
      ,CHFt_CONST_VR(b)
      ,CHFt_CONST_VR(c)
      ,CHFt_CONST_VR(d)
      ,CHFt_CONST_VR(e)
      ,CHFt_CONST_VR(f) );
}
#endif  // GUARDSPLINE5 

#ifndef GUARDADDKETOU 
#define GUARDADDKETOU 
// Prototype for Fortran procedure addketou ...
//
void FORTRAN_NAME( ADDKETOU ,addketou )(
      CHFp_FRA(u)
      ,CHFp_BOX(box) );

#define FORT_ADDKETOU FORTRAN_NAME( inlineADDKETOU, inlineADDKETOU)
#define FORTNT_ADDKETOU FORTRAN_NAME( ADDKETOU, addketou)

inline void FORTRAN_NAME(inlineADDKETOU, inlineADDKETOU)(
      CHFp_FRA(u)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_ADDKETOU");
 FORTRAN_NAME( ADDKETOU ,addketou )(
      CHFt_FRA(u)
      ,CHFt_BOX(box) );
}
#endif  // GUARDADDKETOU 

#ifndef GUARDQUINAT 
#define GUARDQUINAT 
// Prototype for Fortran procedure QUINAT ...
//
void FORTRAN_NAME( QUINAT ,quinat )(
      CHFp_CONST_VR(X)
      ,CHFp_CONST_VR(Y)
      ,CHFp_VR(B)
      ,CHFp_VR(C)
      ,CHFp_VR(D)
      ,CHFp_VR(E)
      ,CHFp_VR(F) );

#define FORT_QUINAT FORTRAN_NAME( inlineQUINAT, inlineQUINAT)
#define FORTNT_QUINAT FORTRAN_NAME( QUINAT, quinat)

inline void FORTRAN_NAME(inlineQUINAT, inlineQUINAT)(
      CHFp_CONST_VR(X)
      ,CHFp_CONST_VR(Y)
      ,CHFp_VR(B)
      ,CHFp_VR(C)
      ,CHFp_VR(D)
      ,CHFp_VR(E)
      ,CHFp_VR(F) )
{
 CH_TIMELEAF("FORT_QUINAT");
 FORTRAN_NAME( QUINAT ,quinat )(
      CHFt_CONST_VR(X)
      ,CHFt_CONST_VR(Y)
      ,CHFt_VR(B)
      ,CHFt_VR(C)
      ,CHFt_VR(D)
      ,CHFt_VR(E)
      ,CHFt_VR(F) );
}
#endif  // GUARDQUINAT 

}

#endif
