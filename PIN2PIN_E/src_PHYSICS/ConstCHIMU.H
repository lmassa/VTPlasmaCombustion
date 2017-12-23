#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CONSTCHIMU_H_
#define _CONSTCHIMU_H_

#include  <iostream>

#include "gasModel.H"
#include "NamespaceHeader.H"

class ConstCHIMU : public gasModel
{
public:

  

  virtual void setPlasmaParameters(plasmaParameters& a_ppars, Vector<Real>& a_minData, Vector<Real>& a_maxData,  string& a_bgSpec, Vector<string>& a_charSpec);
  
  virtual void initialize(LevelData<EBCellFAB>&       a_conState,
			  const plasmaParameters&     a_ppars,
			  const EBISLayout&           a_ebisl,
			  const Real&                 a_dx,
			  const RealVect&             a_domainLength);
  
  virtual Real setImplicitTimeStep() const;
  virtual void setWallBCinfo(Vector<int>& a_wallBCInfo) const;

  Real m_domainLength;
  int m_nSpecies, m_eIndex, m_iIndex;
  
};

#include "NamespaceFooter.H"
#endif
