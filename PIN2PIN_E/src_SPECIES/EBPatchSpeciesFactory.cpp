#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EBPatchSpeciesFactory.H"


/******************/
EBPatchSpeciesFactory::
EBPatchSpeciesFactory() :EBPatchGodunovFactory()
{
}

EBPatchSpeciesFactory::
EBPatchSpeciesFactory(const EBPhysIBCFactory* a_bcFactoryPtr,
		      PlasmaPhysics*    a_PlasmaPhysics,
		      const bool&             a_useLimiting,
		      const int&              a_nSpec,
		      const Real&             a_maxVal,
		      const Real&             a_minVal,
		      const bool&             a_setMaxMin)
  :EBPatchGodunovFactory()
{
  
  define(a_bcFactoryPtr,a_PlasmaPhysics, a_useLimiting,a_nSpec,a_maxVal, a_minVal,a_setMaxMin);
}
//-----------------------------------------------------------------------
void EBPatchSpeciesFactory::
define(const EBPhysIBCFactory*     a_bcFactoryPtr,
       PlasmaPhysics*        a_PlasmaPhysics,
       const bool&                 a_useLimiting,
       const int&                  a_nSpec,
       const Real&                 a_maxVal,
       const Real&                 a_minVal,
       const bool&                 a_setMaxMin)

{
  
  m_bcFactoryPtr = a_bcFactoryPtr;
  m_PlasmaPhysics = a_PlasmaPhysics;
  m_useLimiting  = a_useLimiting;
  m_setMaxMin    = a_setMaxMin;
  m_maxVal       = a_maxVal;
  m_minVal       = a_minVal;
  m_nSpec        = a_nSpec;
}
/******************/
EBPatchSpeciesFactory::
~EBPatchSpeciesFactory()
{
}
/******************/
EBPatchSpecies*
EBPatchSpeciesFactory::
create() const
{
  EBPatchSpecies* retval= new EBPatchSpecies();
  retval->setEBPhysIBC(*m_bcFactoryPtr);
  retval->useLimiting(m_useLimiting);
  retval->setSpeciesNumber(m_nSpec);
  if (m_setMaxMin)
    {
      retval->setMaxMin(m_maxVal, m_minVal);
    }
  return static_cast<EBPatchSpecies*>(retval);
}
/******************/



