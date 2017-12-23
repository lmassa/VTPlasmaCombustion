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
#include "PoissonDomMixBcFunc.H"
#include "EBAMRPoissonOp.H"


PoissonDomMixBcFunc::PoissonDomMixBcFunc()
{
  m_isDefined = false;
  vector<Real> m_Dirich_xint;
}

PoissonDomMixBcFunc::~PoissonDomMixBcFunc()
{
}

void PoissonDomMixBcFunc::define(const vector<Real>& a_Dirich_xint, const vector<Real>& a_dom_yint, const Real& a_domDirBcValue, const Real& a_domNeuBcValue, const vector<Real>& a_DirEB_xint, const Real& a_EbDirBcValue, const Real& a_omega)
{
  m_isDefined = true;
  m_Dirich_xint = a_Dirich_xint;
  m_dom_yint = a_dom_yint;
  m_domDirBcValue = a_domDirBcValue;
  m_domNeuBcValue = a_domNeuBcValue;
  m_DirEB_xint = a_DirEB_xint;
  m_EbDirBcValue = a_EbDirBcValue;
  m_omega = a_omega;
}

bool PoissonDomMixBcFunc::isDirichlet(const RealVect& a_point, const Real& a_dx) const
{
   //CH_assert(m_isDefined);
  
  return isDirichlet(a_point);
}
bool PoissonDomMixBcFunc::isDirichlet(const RealVect& a_point) const
{
   //CH_assert(m_isDefined);
   bool flag;
	if (((a_point[0]>=m_Dirich_xint[0]) && (a_point[0]<=m_Dirich_xint[1]) && (a_point[1]<= m_dom_yint[0])) || (a_point[1]>=m_dom_yint[2]))
		{
			//Dirichlet BC
			flag=true;
		}
	else
		{
			//Neumann BC
			flag=false;
		}
   return flag;
}

// called for bottom boundary only
Real PoissonDomMixBcFunc::valueDir(
						const RealVect& a_point,
                        const RealVect& a_normal,
                        const Real&     a_time,
                        const int&      a_comp) const
{
  CH_assert(m_isDefined);
  Real value;
  //if(a_point[1]<1e-4)
    value=0;
  //else //MDCancel: for Soloviev top Dirichlet 
  //  value = m_EbDirBcValue*sin(m_omega*a_time)*(0.5-0.25/atan(1.0)*atan((a_point[0]-m_DirEB_xint[1])/(a_point[1]-m_dom_yint[1])));
  
  return value;
}

// called for bottom boundary only
Real PoissonDomMixBcFunc::valueNeu(
						const RealVect& a_point,
                        const RealVect& a_normal,
                        const Real&     a_time,
                        const int&      a_comp) const
{
	CH_assert(m_isDefined);
	Real value;
	value=m_domNeuBcValue;	
    return value;
}


