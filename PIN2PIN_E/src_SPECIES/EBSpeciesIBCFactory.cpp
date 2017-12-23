#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EBSpeciesIBCFactory.H"
#include "EBSpeciesIBC.H"
/******************/
/******************/
EBSpeciesIBCFactory::
EBSpeciesIBCFactory()
  :EBPhysIBCFactory()
{
  m_PlasmaPhysics = NULL;
}
/******************/
/******************/

EBSpeciesIBCFactory::
~EBSpeciesIBCFactory()
{
}
/******************/

void EBSpeciesIBCFactory::
define( const int&      a_nSpec, PlasmaPhysics*  a_PlasmaPhysics)
{ 
  m_nSpec = a_nSpec;
  m_PlasmaPhysics=a_PlasmaPhysics;
}
/******************/
EBPhysIBC*
EBSpeciesIBCFactory::
create() const
{
  EBSpeciesIBC* retval =
    new EBSpeciesIBC(m_nSpec, m_PlasmaPhysics);
  return static_cast<EBPhysIBC*>(retval);
}
/******************/
/******************/

