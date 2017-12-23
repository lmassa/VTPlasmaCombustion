#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EBPin2PinIBCFactory.H"
#include "EBPin2PinIBC.H"
/******************/
/******************/
EBPin2PinIBCFactory::
EBPin2PinIBCFactory()
  :EBPhysIBCFactory()
{
  m_BLPtr = NULL;
}
/*****************/
EBPin2PinIBCFactory::
EBPin2PinIBCFactory(const Real&     a_gamma,
                        const Real&     a_ms,
                        const Real&     a_center,
                        const int&      a_shocknorm,
                        const bool&     a_shockbackward,
                        const bool&     a_doRZCoords,
			const BoundaryLayer* const a_BLPtr)
  :EBPhysIBCFactory()
{
  
  define(a_gamma,a_ms,a_center,a_shocknorm,a_shockbackward,a_doRZCoords,a_BLPtr);

}
/******************/
void EBPin2PinIBCFactory:: 
define(const Real&     a_gamma,
       const Real&     a_ms,
       const Real&     a_center,
       const int&      a_shocknorm,
       const bool&     a_shockbackward,
       const bool&     a_doRZCoords,
       const BoundaryLayer* const a_BLPtr)
{
    if (a_gamma < 1.0)
    {
      MayDay::Error("bogus gamma in EBPin2Pin");
    }
  m_gamma         = a_gamma;
  m_ms            = a_ms;
  m_center        = a_center;
  m_shocknorm     = a_shocknorm;
  m_shockbackward = a_shockbackward;
  m_doRZCoords    = a_doRZCoords;
  m_BLPtr         = a_BLPtr;
}
/******************/

EBPin2PinIBCFactory::
~EBPin2PinIBCFactory()
{
}
/******************/
/******************/
EBPhysIBC*
EBPin2PinIBCFactory::
create() const
{
  EBPin2PinIBC* retval =
    new EBPin2PinIBC(m_gamma,
                         m_ms,
                         m_center,
                         m_shocknorm,
                         m_shockbackward,
                         m_doRZCoords,
			 m_BLPtr);
  return static_cast<EBPhysIBC*>(retval);
}
/******************/
/******************/

