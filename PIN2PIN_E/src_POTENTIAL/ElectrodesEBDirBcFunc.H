#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _ELECTRODESEBDIRBCFUNC_H_
#define _ELECTRODESEBDIRBCFUNC_H_

#include "RealVect.H"
#include "BaseBCValue.H"
#include "MayDay.H"
#include "UsingNamespace.H"

class ElectrodesEBDirBcFunc: public BaseBCValue
{
public:
  ElectrodesEBDirBcFunc();

  virtual ~ElectrodesEBDirBcFunc();

  virtual void define(const Real& a_Potential, const Real& a_omega, const Real& a_xloc, Real* a_time);


  virtual Real value(const RealVect& a_point, const RealVect& a_normal, const Real& a_time, const int& a_comp) const;
  


private:
  bool m_isDefined;
  Real m_Potential;
  Real m_xloc;
  Real m_omega;
  Real* m_time;
};

#endif
