#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBAMRNAVIERFACTORY_H_
#define _EBAMRNAVIERFACTORY_H_

#include "EBAMRGodunov.H"
#include "EBAMRNavier.H"
#include "AMRLevelFactory.H"
#include "NamespaceHeader.H"

///
/**
 */
class EBAMRNavierFactory : public AMRLevelFactory
{
public:
  ///
  AMRLevel* new_amrlevel() const;

  ///
  EBAMRNavierFactory(const Real&                   a_initialDtMultiplier,
		     const Real&                   a_cfl,
		     const int &                   a_redistRad,
		     const RealVect&               a_domainLength,
		     const Real&                   a_refineThresh,
		     const int &                   a_tagBufferSize,
		     const int &                   a_verbosity,
		     const bool&                   a_useMassRedist,
		     const bool&                   a_doSmushing,
		     const bool&                   a_doRZCoords,
		     const bool&                   a_hasSourceTerm,
		     const EBPatchGodunovFactory*  const a_patchGodunov,
		     int                           a_symmetric,
		     Vector<RefCountedPtr<LevelData<EBCellFAB> > >*        a_aco,
		     Vector<RefCountedPtr<LevelData<EBFluxFAB> > >*        a_eta,
		     Vector<RefCountedPtr<LevelData<EBFluxFAB> > >*        a_lambda,
		     Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >* a_etaIrreg,
		     Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >* a_lambdaIrreg,
		     Vector<RefCountedPtr<LevelData<EBCellFAB> > >*        a_rcv,
		     bool                                                        a_tagAll=false);

  ///
  virtual ~EBAMRNavierFactory();

protected:
  bool                   m_tagAll;
  Real                   m_initialDtMultiplier;
  Real                   m_cfl;
  int                    m_redistRad;
  RealVect               m_domainLength;
  Real                   m_refineThresh;
  int                    m_tagBufferSize;
  int                    m_verbosity;
  bool                   m_useMassRedist;
  bool                   m_doSmushing;
  bool                   m_doRZCoords;
  bool                   m_hasSourceTerm;
  const EBPatchGodunovFactory* m_patchGodunov;
  int m_symmetric;

  Vector<RefCountedPtr<LevelData<EBCellFAB> > >*         m_aco;
  Vector<RefCountedPtr<LevelData<EBFluxFAB> > >*         m_eta;
  Vector<RefCountedPtr<LevelData<EBFluxFAB> > >*         m_lambda;
  Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >*  m_etaIrreg;
  Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >*  m_lambdaIrreg;
  Vector<RefCountedPtr<LevelData<EBCellFAB> > >*         m_rcv;


private:
  EBAMRNavierFactory()
  {
    MayDay::Error("invalid operator");
  }

  //disallowed for all the usual reasons
  void operator=(const EBAMRNavierFactory& a_input)
  {
    MayDay::Error("invalid operator");
  }
  EBAMRNavierFactory(const EBAMRNavierFactory& a_input)
  {
    MayDay::Error("invalid operator");
  }

};

#include "NamespaceFooter.H"
#endif
