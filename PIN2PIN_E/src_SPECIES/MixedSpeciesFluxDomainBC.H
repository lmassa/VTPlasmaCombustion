#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MIXEDSPECIESFLUXDOMAINBC_H_
#define _MIXEDSPECIESFLUXDOMAINBC_H_

#include "RefCountedPtr.H"
#include "BaseDomainBC.H"
#include "BaseBCFuncEval.H"
#include "PlasmaPhysics.H"
#include "NamespaceHeader.H"


///
/**
 */
class SpeciesFluxBaseDomainBC: public BaseDomainBC
{
public:
  ///
  /**
   */
  SpeciesFluxBaseDomainBC()
  {
    m_coefSet  = false;
  }

  virtual ~SpeciesFluxBaseDomainBC()
  {
  }

  ///
  /**
   */ 
  //from BaseDomainBC
  void setCoef(EBLevelGrid                         &  a_eblg,
               Real                                &  a_beta,
               RefCountedPtr<LevelData<EBFluxFAB> >&  a_Dk, 
	       RefCountedPtr<LevelData<EBFluxFAB> >&  a_vecW)
  {
    m_coefSet = true;
    m_eblg    = a_eblg;
    m_beta    = a_beta;
    m_Dk      = a_Dk;
    m_vecW    = a_vecW;
  }
  // for PlasmaPhysics
  void setPhysics(PlasmaPhysics*  a_PlasmaPhysics, const int&  a_ispec);

protected:
  //from BaseDomainBC
  bool                                  m_coefSet;
  EBLevelGrid                           m_eblg;
  Real                                  m_beta;
  RefCountedPtr<LevelData<EBFluxFAB> >  m_Dk;
  RefCountedPtr<LevelData<EBFluxFAB> >  m_vecW;
  // proper to plasma physics
  int m_ispec;
  PlasmaPhysics*  m_PlasmaPhysics;
  bool m_onlyHomogeneous;
};
/**************************************///////////////////////
/**************************************///////////////////////
/**************************************///////////////////////
/**************************************///////////////////////
/**************************************///////////////////////
/**************************************///////////////////////
/**************************************///////////////////////
/**************************************///////////////////////
/**************************************///////////////////////
class MixedSpeciesFluxDomainBC: public SpeciesFluxBaseDomainBC
{
public:
  MixedSpeciesFluxDomainBC();

  virtual ~MixedSpeciesFluxDomainBC();


  //from Conductivity...

  virtual void getFaceFlux(BaseFab<Real>&        a_faceFlux,
                           const BaseFab<Real>&  a_phi,
                           const RealVect&       a_probLo,
                           const RealVect&       a_dx,
                           const int&            a_idir,
                           const Side::LoHiSide& a_side,
                           const DataIndex&      a_dit,
                           const Real&           a_time,
                           const bool&           a_useHomogeneous);

  virtual void getFaceFlux(Real&                 a_faceFlux,
                           const VolIndex&       a_vof,
                           const int&            a_comp,
                           const EBCellFAB&      a_phi,
                           const RealVect&       a_probLo,
                           const RealVect&       a_dx,
                           const int&            a_idir,
                           const Side::LoHiSide& a_side,
                           const DataIndex&      a_dit,
                           const Real&           a_time,
                           const bool&           a_useHomogeneous);


private:
  virtual void getFaceGradPhi(Real&                 a_faceFlux,
                              const FaceIndex&      a_face,
                              const int&            a_comp,
                              const EBCellFAB&      a_phi,
                              const RealVect&       a_probLo,
                              const RealVect&       a_dx,
                              const int&            a_idir,
                              const Side::LoHiSide& a_side,
                              const DataIndex&      a_dit,
                              const Real&           a_time,
                              const bool&           a_useAreaFrac,
                              const RealVect&       a_centroid,
                              const bool&           a_useHomogeneous);

  virtual void getFaceVel(Real&                 a_faceFlux,
                          const FaceIndex&      a_face,
                          const EBFluxFAB&      a_vel,
                          const RealVect&       a_probLo,
                          const RealVect&       a_dx,
                          const int&            a_idir,
                          const int&            a_icomp,
                          const Real&           a_time,
                          const Side::LoHiSide& a_side,
                          const bool&           a_doDivFreeOutflow);

  void getInsideVal(Real&                 a_faceVal,
		    const FaceIndex&      a_face,
		    const int&            a_comp,
		    const EBCellFAB&      a_phi,
		    const int&            a_idir,
		    const Side::LoHiSide& a_side,
		    const DataIndex&      a_dit,
		    const Real&           a_time,
		    const bool&           a_useAreaFrac);

private:

};

/**************************************///////////////////////
/**************************************///////////////////////
/**************************************///////////////////////
/**************************************///////////////////////
/**************************************///////////////////////
/**************************************///////////////////////
/**************************************///////////////////////
/**************************************///////////////////////
/**************************************///////////////////////
class MixedSpeciesFluxDomainBCFactory: public BaseDomainBCFactory
{
public:
  MixedSpeciesFluxDomainBCFactory();
  
  virtual ~MixedSpeciesFluxDomainBCFactory();

  virtual MixedSpeciesFluxDomainBC* create(const ProblemDomain& a_domain,
					   const EBISLayout&    a_layout,
					   const RealVect&      a_dx);

  void setPhysics(PlasmaPhysics*  a_PlasmaPhysics, const int&  a_ispec);
private:
  bool m_onlyHomogeneous;

  int m_ispec;
  PlasmaPhysics*  m_PlasmaPhysics;
};
#include "NamespaceFooter.H"
#endif
