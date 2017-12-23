#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MIXEDCONDUCTIVITYDOMAINBC_H_
#define _MIXEDCONDUCTIVITYDOMAINBC_H_

#include "RefCountedPtr.H"
#include "BaseDomainBC.H"

#include "DirichletPoissonDomainBC.H"
#include "BaseMixBCValue.H"
#include "NamespaceHeader.H"

class MixedConductivityDomainBC: public ConductivityBaseDomainBC
{
public:
MixedConductivityDomainBC();

   virtual ~MixedConductivityDomainBC();

   void setArguments(Vector<Vector<int> >&  a_domMixBc, Real a_dvalue, Real a_nvalue, RefCountedPtr<BaseMixBCValue> a_func);  

  virtual void setEBOrder(int a_ebOrder);

  virtual void setFunction(RefCountedPtr<BaseMixBCValue> a_flux);

  virtual void getFaceFlux(BaseFab<Real>&        a_faceFlux,
                           const BaseFab<Real>&  a_phi,
                           const RealVect&       a_probLo,
                           const RealVect&       a_dx,
                           const int&            a_idir,
                           const Side::LoHiSide& a_side,
                           const DataIndex&      a_dit,
                           const Real&           a_time,
                           const bool&           a_useHomogeneous);


  virtual bool isDirichletDom(const VolIndex&   a_ivof,
                              const VolIndex&   a_jvof,
                              const EBCellFAB&  a_phi)
  {
    return true;
  }

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

  ///
  /**
     This is used by the projections to get grad(phi) at domain faces (usually outflow).
     Called by InflowOutflowPoissonDomainBC::getFaceFlux.
   */
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
                              const bool&           a_useHomogeneous);/// 
  /**
     From Neumann:
	 This is used by the projections to get extrapolated
     grad(phi) at domain faces for cell-centered gradient.
     Called by getFaceFlux.
   */
  virtual void getFaceFluxGradPhi(Real&                 a_faceFlux,
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

  bool m_isDefined;

  void setTime(Real* a_time)
  {
    if(m_timePtr == NULL) m_isDefined = false;
    m_timePtr = a_time;
  }
private:
  bool m_isFunction;

  Vector<Vector<int> >  m_domMixBc;
  Real m_dvalue;
  Real m_nvalue;
  RefCountedPtr<BaseMixBCValue> m_func;
  bool m_onlyHomogeneous;
  int m_ebOrder;
  Real* m_timePtr;
};

class MixedConductivityDomainBCFactory: public BaseDomainBCFactory
{
public:
  MixedConductivityDomainBCFactory();

  virtual ~MixedConductivityDomainBCFactory();

  void setArguments(Vector<Vector<int> >& a_domMixBc, Real a_dvalue, Real a_nvalue, RefCountedPtr<BaseMixBCValue> a_func);
  virtual void setEBOrder(int a_ebOrder);

  void setTime(Real* a_time);

  virtual void setFunction(RefCountedPtr<BaseMixBCValue> a_flux);

  virtual MixedConductivityDomainBC* create(const ProblemDomain& a_domain,
                                                const EBISLayout&    a_layout,
                                                const RealVect&      a_dx);

private:
  Vector<Vector<int> >   m_domMixBc;
  bool m_onlyHomogeneous;
  Real m_dvalue;
  Real m_nvalue;
  bool m_isFunction;
  RefCountedPtr<BaseMixBCValue> m_func;
  int m_ebOrder;
  Real* m_timePtr;
};
#include "NamespaceFooter.H"
#endif
