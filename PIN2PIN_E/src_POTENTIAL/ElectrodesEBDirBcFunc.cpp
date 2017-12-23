#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "RealVect.H"
#include "EBArith.H"
#include "EBArithF_F.H"
#include "functionsF_F.H"
#include "ElectrodesEBDirBcFunc.H"
// #include "EBAMRPoissonOp.H"


ElectrodesEBDirBcFunc::ElectrodesEBDirBcFunc()
{
  m_isDefined = false;
}

ElectrodesEBDirBcFunc::~ElectrodesEBDirBcFunc()
{
}

void ElectrodesEBDirBcFunc::define(const Real& a_Potential, const Real& a_omega, const Real& a_xloc, Real* a_time)
{
  m_isDefined = true;
  m_Potential = a_Potential;
  m_xloc = a_xloc;
  m_omega = a_omega;
  m_time = a_time;
}


// called for bottom boundary only
Real ElectrodesEBDirBcFunc::value(const RealVect& a_point,
				  const RealVect& a_normal,
				  const Real&     a_time,
				  const int&      a_comp) const
{
  CH_assert(m_isDefined);
  Real time = *m_time;
  Real value;
  if (a_point[0]>=m_xloc)
    // RHS electrode
    if(m_omega <= 0)
      value=m_Potential*cos(m_omega*time);
    else
      value=m_Potential*sin(m_omega*time);
  else 
    // LHS electrode - MD this is only for direct discharge 
    value = 0;
  return value;
}
