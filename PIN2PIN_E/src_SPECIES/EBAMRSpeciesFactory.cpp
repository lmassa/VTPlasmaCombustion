#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "AMRLevel.H"
#include "EBAMRSpecies.H"
#include "EBAMRSpeciesFactory.H"
#include "NamespaceHeader.H"
/****************/
/****************/
EBAMRSpeciesFactory::
EBAMRSpeciesFactory(const Real&                   a_initialDtMultiplier,
		    const Real&                   a_cfl,
		    const int &                   a_redistRad,
		    const RealVect&               a_domainLength,
		    const Real&                   a_refineThresh,
		    const int &                   a_tagBufferSize,
		    const int &                   a_verbosity,
		    const bool&                   a_useLimiting,
		    const bool&                   a_doSmushing,
		    const bool&                   a_doRZCoords,
		    const bool&                   a_hasSourceTerm,
		    const EBPatchSpeciesFactory* const a_patchSpeciesFactory,
		    bool                          a_tagAll)
{
  m_tagAll = a_tagAll;

  m_initialDtMultiplier = a_initialDtMultiplier;
  m_cfl                 = a_cfl;
  m_useLimiting         = a_useLimiting;
  m_redistRad           = a_redistRad;
  m_domainLength        = a_domainLength;
  m_refineThresh        = a_refineThresh;
  m_tagBufferSize       = a_tagBufferSize;
  m_verbosity           = a_verbosity;
  m_doSmushing          = a_doSmushing;
  m_doRZCoords          = a_doRZCoords;
  m_hasSourceTerm       = a_hasSourceTerm;
  m_patchSpeciesFactory = a_patchSpeciesFactory;
}
/****************/
/****************/
AMRLevel*
EBAMRSpeciesFactory::
new_amrlevel() const
{
  EBAMRSpecies* amrg_ptr = new EBAMRSpecies();

  amrg_ptr->define(m_cfl,m_domainLength,m_refineThresh,m_tagBufferSize,
		   m_initialDtMultiplier,m_useLimiting);
  amrg_ptr->patchSpecies(m_patchSpeciesFactory);
  amrg_ptr->redistRadius(m_redistRad);
  amrg_ptr->verbosity(m_verbosity);
  amrg_ptr->setPhysics(m_patchSpeciesFactory->m_PlasmaPhysics);

  return (static_cast <AMRLevel*> (amrg_ptr));
}
/****************/
/****************/
EBAMRSpeciesFactory::
~EBAMRSpeciesFactory()
{
}
/****************/
/****************/
#include "NamespaceFooter.H"
