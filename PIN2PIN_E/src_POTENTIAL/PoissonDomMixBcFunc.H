#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _POISSONDOMMIXBCFUNC_H_
#define _POISSONDOMMIXBCFUNC_H_

#include "RealVect.H"
#include "BaseMixBCValue.H"
#include "MayDay.H"
#include "UsingNamespace.H"

// LM changed  BaseBCValue -> BaseMixBCValue
class PoissonDomMixBcFunc: public BaseMixBCValue
{
public:
  PoissonDomMixBcFunc();

  virtual ~PoissonDomMixBcFunc();

  virtual void define(const vector<Real>& a_Dirich_xint, const vector<Real>& a_dom_yint, const Real& a_domDirBcValue, const Real& a_domNeuBcValue, const vector<Real>& a_DirEB_xint, const Real& a_EbDirBcValue, const Real& a_omega);

  virtual bool isDirichlet(const RealVect& a_point, const Real& a_dx) const;
  virtual bool isDirichlet(const RealVect& a_point) const;

  virtual Real valueDir(const RealVect& a_point, const RealVect& a_normal, const Real& a_time, const int& a_comp) const;

  virtual Real valueNeu(const RealVect& a_point, const RealVect& a_normal, const Real& a_time, const int& a_comp) const;


private:
  bool m_isDefined;
  Real m_domDirBcValue;
  Real m_domNeuBcValue;
  Real m_EbDirBcValue;
  Real m_omega;
  vector<Real> m_Dirich_xint;
  vector<Real> m_dom_yint;
  vector<Real> m_DirEB_xint;
};

#endif
