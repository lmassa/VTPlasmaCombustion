#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _BASESPECIESFLUXEBBC_H_
#define _BASESPECIESFLUXEBBC_H_

#include "RefCountedPtr.H"

#include "NeumannPoissonEBBC.H"
#include "BaseBCValue.H"
#include "NamespaceHeader.H"

///
/**
 */

class SpeciesFluxBaseEBBC: public BaseEBBC
{
public:
  ///
  /**
   */
  SpeciesFluxBaseEBBC()
  {
    m_coefSet  = false;
    m_dataBased = false;
  }

  virtual ~SpeciesFluxBaseEBBC()
  {
  }

  virtual void define(const LayoutData<IntVectSet>& a_cfivs,
                      const Real&                   a_factor)
  {
    ;
  }
  virtual LayoutData<BaseIVFAB<VoFStencil> >* getFluxStencil(int ivar)
  {
    return NULL;
  }

  ///
  /**
   */
  void setCoef(EBLevelGrid                                &  a_eblg,
               Real                                       &  a_beta,
               RefCountedPtr<LevelData<BaseIVFAB<Real> > >&  a_Dk,
               RefCountedPtr<LevelData<BaseIVFAB<Real> > >&  a_vecW)
  {
    m_coefSet = true;
    //this is really fakeBeta = 1
    m_beta    = a_beta;
    m_eblg    = a_eblg;
    m_Dk      = a_Dk;
    m_vecW    = a_vecW;
  }

  
  ///
  virtual void setData( RefCountedPtr<LevelData<BaseIVFAB<Real> > >& a_data, int a_ispec, Vector<int>& a_wallBCInfo)
  {
    m_data = a_data;
    m_dataBased = true;
    m_wallBCInfo = a_wallBCInfo;
    m_ispec = a_ispec;
  }

protected:

 
  //from the base conductivity
  EBLevelGrid                                  m_eblg;
  bool                                         m_coefSet;
  Real                                         m_beta;
  RefCountedPtr<LevelData<BaseIVFAB<Real> > >  m_Dk;  //should not be used in the flux BC
  RefCountedPtr<LevelData<BaseIVFAB<Real> > >  m_vecW;
  int m_ispec;
  RefCountedPtr<BaseBCValue> m_func;
  Vector<int> m_wallBCInfo;

};




#include "NamespaceFooter.H"
#endif
