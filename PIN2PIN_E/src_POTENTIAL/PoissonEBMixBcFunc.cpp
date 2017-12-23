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
#include "PoissonEBMixBcFunc.H"
#include "EBAMRPoissonOp.H"


PoissonEBMixBcFunc::PoissonEBMixBcFunc()
{
  m_isDefined = false;
  vector<Real> m_Eb_Dir_xint;
}

PoissonEBMixBcFunc::~PoissonEBMixBcFunc()
{
}

void PoissonEBMixBcFunc::define(const vector<Real>& a_Eb_Dir_xint, const Real& a_EbDirBcValue, const Real& a_EbNeuBcValue, const Real& a_omega)
{
  m_isDefined = true;
  m_Eb_Dir_xint = a_Eb_Dir_xint;
  m_EbDirBcValue = a_EbDirBcValue;
  m_EbNeuBcValue = a_EbNeuBcValue;
  m_omega = a_omega;
}

bool PoissonEBMixBcFunc::isDirichlet(const RealVect& a_point) const
{
   //CH_assert(m_isDefined);
   bool isDirichlet;
   if (((a_point[0]>=m_Eb_Dir_xint[0]) && (a_point[0]<=m_Eb_Dir_xint[1])) || ((a_point[0]>=m_Eb_Dir_xint[2]) && (a_point[0]<=m_Eb_Dir_xint[3])))
     {
		  //Dirichlet BC
       isDirichlet=true;
     }
   else
     {
       //Neumann BC
       isDirichlet=false;
     }
   return isDirichlet;
}

// called for bottom boundary only
Real PoissonEBMixBcFunc::valueDir(const RealVect& a_point,
				  const RealVect& a_normal,
				  const Real&     a_time,
				  const int&      a_comp) const
{
  CH_assert(m_isDefined);
  Real value;
  if ((a_point[0]>=m_Eb_Dir_xint[0]) && (a_point[0]<=m_Eb_Dir_xint[1]))
    // LHS electrode
    if(m_omega <= 0)
      value=m_EbDirBcValue*cos(m_omega*a_time);
    else
      value=m_EbDirBcValue*sin(m_omega*a_time);
  else 
    // RHS electrode - MD this is only for direct discharge 
   // for DBD use DIrichlet on bottom domain (plasma.inputs), pass only 2 x-coordinates to this function in ElectricPotential and remove second "if" condition in isDirichlet above
    value = 0;
  return value;
}

// called for bottom boundary only
Real PoissonEBMixBcFunc::valueNeu(
						const RealVect& a_point,
                        const RealVect& a_normal,
                        const Real&     a_time,
                        const int&      a_comp) const
{
	CH_assert(m_isDefined);
	Real value;
	value=m_EbNeuBcValue;	
    return value;
}


