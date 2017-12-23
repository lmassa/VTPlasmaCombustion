#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBAMRSPECIESFACTORY_H_
#define _EBAMRSPECIESFACTORY_H_

#include "EBLevelSpecies.H"
#include "AMRLevelFactory.H"
#include "EBPatchSpeciesFactory.H"
#include "NamespaceHeader.H"

///
/**
 */
class EBAMRSpeciesFactory : public AMRLevelFactory
{
public:
  ///
  AMRLevel* new_amrlevel() const;

  ///
  EBAMRSpeciesFactory(const Real&                        a_initialDtMultiplier,
                      const Real&                        a_cfl,
                      const int &                        a_redistRad,
                      const RealVect&                    a_domainLength,
                      const Real&                        a_refineThresh,
                      const int &                        a_tagBufferSize,
                      const int &                        a_verbosity,
                      const bool&                        a_useLimiting,
                      const bool&                        a_doSmushing,
                      const bool&                        a_doRZCoords,
                      const bool&                        a_hasSourceTerm,
                      const EBPatchSpeciesFactory* const a_patchSpeciesFactory,
                      bool                               a_tagAll = false);

  ///
  virtual ~EBAMRSpeciesFactory();

protected:
  bool                   m_tagAll;
  Real                   m_initialDtMultiplier;
  Real                   m_cfl;
  int                    m_redistRad;
  RealVect               m_domainLength;
  Real                   m_refineThresh;
  int                    m_tagBufferSize;
  int                    m_verbosity;
  bool                   m_useLimiting;
  bool                   m_doSmushing;
  bool                   m_doRZCoords;
  bool                   m_hasSourceTerm;
  const EBPatchSpeciesFactory* m_patchSpeciesFactory;


private:
  EBAMRSpeciesFactory()
  {
    MayDay::Error("invalid operator");
  }

  //disallowed for all the usual reasons
  void operator=(const EBAMRSpeciesFactory& a_input)
  {
    MayDay::Error("invalid operator");
  }
  EBAMRSpeciesFactory(const EBAMRSpeciesFactory& a_input)
  {
    MayDay::Error("invalid operator");
  }

};

#include "NamespaceFooter.H"
#endif
