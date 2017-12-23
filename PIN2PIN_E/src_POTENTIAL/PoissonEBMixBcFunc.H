#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _POISSONEBMIXBCFUNC_H_
#define _POISSONEBMIXBCFUNC_H_

#include "RealVect.H"
#include "BaseMixBCValue.H"
#include "MayDay.H"
#include "UsingNamespace.H"

// LM changed  BaseBCValue -> BaseMixBCValue
class PoissonEBMixBcFunc: public BaseMixBCValue
{
public:
  PoissonEBMixBcFunc();

  virtual ~PoissonEBMixBcFunc();

  virtual void define(const vector<Real>& a_Eb_Dir_xint, const Real& a_EbDirBcValue, const Real& a_EbNeuBcValue, const Real& a_omega);

  virtual bool isDirichlet(const RealVect& a_point) const;

  virtual Real valueDir(const RealVect& a_point, const RealVect& a_normal, const Real& a_time, const int& a_comp) const;

  virtual Real valueNeu(const RealVect& a_point, const RealVect& a_normal, const Real& a_time, const int& a_comp) const;
  


private:
  bool m_isDefined;
  Real m_EbDirBcValue;
  Real m_EbNeuBcValue;
  Real m_omega;
  vector<Real> m_Eb_Dir_xint;
};

#endif
