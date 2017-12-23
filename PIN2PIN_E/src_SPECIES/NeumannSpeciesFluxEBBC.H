#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _NEUMANNSPECIESFLUXEBBC_H_
#define _NEUMANNSPECIESFLUXEBBC_H_

#include "RefCountedPtr.H"

#include "NeumannPoissonEBBC.H"
#include "SpeciesFluxBaseEBBC.H"
#include "BaseBCValue.H"
#include "NamespaceHeader.H"

///
/**
 */



class NeumannSpeciesFluxEBBC: public SpeciesFluxBaseEBBC
{
public:
  /// constructor
  NeumannSpeciesFluxEBBC(const ProblemDomain& a_domain,
			 const EBISLayout&    a_layout,
			 const RealVect&      a_dx);
  ///


  virtual void applyEBFlux(EBCellFAB&                    a_lphi,
                           const EBCellFAB&              a_phi,
                           VoFIterator&                  a_vofit,
                           const LayoutData<IntVectSet>& a_cfivs,
                           const DataIndex&              a_dit,
                           const RealVect&               a_probLo,
                           const RealVect&               a_dx,
                           const Real&                   a_factor,
                           const bool&                   a_useHomogeneous,
                           const Real&                   a_time);

  ///MD added for compatibility with mixed BC 
  virtual void applyEBFlux(EBCellFAB&                    a_lphi,
                           const EBCellFAB&              a_phi,
			   EBCellFAB&                    a_fphi,
                           VoFIterator&                  a_vofit,
                           const LayoutData<IntVectSet>& a_cfivs,
                           const DataIndex&              a_dit,
                           const RealVect&               a_probLo,
                           const RealVect&               a_dx,
                           const Real&                   a_factor,
                           const bool&                   a_useHomogeneous,
                           const Real&                   a_time)
  {
    applyEBFlux( a_lphi,a_phi, a_vofit,a_cfivs,a_dit, a_probLo,a_dx,a_factor,a_useHomogeneous, a_time);
  }


  ///

  ///
  virtual ~NeumannSpeciesFluxEBBC();

  ///
  virtual void setValue(Real a_value);

  ///
  virtual int whichBC()
  {
    return 0;
  };
  
  ///
  virtual void setFunction(RefCountedPtr<BaseBCValue> a_func);


  ///
  virtual void getEBFlux(Real&                         a_flux,
                         const VolIndex&               a_vof,
                         const LevelData<EBCellFAB>&   a_phi,
                         const LayoutData<IntVectSet>& a_cfivs,
                         const DataIndex&              a_dit,
                         const RealVect&               a_probLo,
                         const RealVect&               a_dx,
                         const bool&                   a_useHomogeneous,
                         const Real&                   a_time,
                         const pair<int,Real>*         a_cacheHint=0 );

  virtual void define(const LayoutData<IntVectSet>& a_cfivs,
		      const Real&                   a_factor);
  virtual LayoutData<BaseIVFAB<Real> >& getFluxWeight()
  {
    return m_fluxWeight; //for the secondary emission coeff
  }
  virtual LayoutData<BaseIVFAB<VoFStencil> >* getFluxStencil(int ivar)
  {
    return &m_fluxStencil;
  }

  // m_fluxWeight is the weight of the Dirichlet BC value used when
  // calculating the EB flux as a linear combination of the stencil
  // cell values and the BC value.
  LayoutData<BaseIVFAB<Real> >       m_fluxWeight;
  LayoutData<BaseIVFAB<VoFStencil> > m_fluxStencil;
  ProblemDomain m_domain;
  EBISLayout    m_layout;
  RealVect m_dx;
  
  bool m_isDefined;

private:
  NeumannPoissonEBBC m_bc;
};



/*******************************************/////////////////
/*******************************************/////////////////
/*******************************************/////////////////
/*******************************************/////////////////
/*******************************************/////////////////
/*******************************************/////////////////
/*******************************************/////////////////
/*******************************************/////////////////
/*******************************************/////////////////
/*******************************************/////////////////

///
class NeumannSpeciesFluxEBBCFactory: public BaseEBBCFactory
{
public:
  ///
  NeumannSpeciesFluxEBBCFactory();

  ///
  virtual ~NeumannSpeciesFluxEBBCFactory();

  ///
  virtual void setValue(Real a_value);


  ///
  virtual void setFunction(RefCountedPtr<BaseBCValue> a_func);

  ///
  virtual void setData(RefCountedPtr<LevelData<BaseIVFAB<Real> > >& a_data, int a_ispec, Vector<int>& a_wallBCInfo)
  {
    m_data = a_data;
    m_dataBased = true;
    m_ispec = a_ispec;
    m_wallBCInfo = a_wallBCInfo;
  }

  ///
  virtual NeumannSpeciesFluxEBBC* create(const ProblemDomain& a_domain,
                                     const EBISLayout&    a_layout,
                                     const RealVect&      a_dx,
                                     const IntVect*       a_ghostCellsPhi=0,
                                     const IntVect*       a_ghostCellsRhs=0);

private:
  bool m_onlyHomogeneous;
  bool m_isFunction;
  int m_ispec;
  RefCountedPtr<LevelData<BaseIVFAB<Real> > >m_data;
  bool m_dataBased;

  Real m_value;
  RefCountedPtr<BaseBCValue> m_func;
  Vector<int> m_wallBCInfo;



};

#include "NamespaceFooter.H"
#endif
