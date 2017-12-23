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
#include "IonizationWaveEBDirBcFunc.H"
// #include "EBAMRPoissonOp.H"


IonizationWaveEBDirBcFunc::IonizationWaveEBDirBcFunc()
{
  m_isDefined = false;

  ParmParse pp;
  std::vector<Real> dLArray(SpaceDim);
  pp.getarr("domain_length",dLArray,0,SpaceDim);
  m_maxvalXYZ1 = dLArray[0]/2;
  m_momentum_frequency = -1;
  m_E0 = 1.0; 
  m_Theta = 1e0;
  m_Ehat = 0;
  
  
}

IonizationWaveEBDirBcFunc::~IonizationWaveEBDirBcFunc()
{
}

void IonizationWaveEBDirBcFunc::define(const Real& a_Potential, const Real& a_omega, const Real& a_xloc, Real* a_time)
{
  m_isDefined = true;
  m_Potential = a_Potential; //Not used:: Using the hardwired m_E0
  m_xloc = m_maxvalXYZ1;
  m_omega = a_omega;
  m_time = a_time;
}


// called for bottom boundary only
Real IonizationWaveEBDirBcFunc::value(const RealVect& a_point,
				  const RealVect& a_normal,
				  const Real&     a_time,
				  const int&      a_comp) const
{
  CH_assert(m_isDefined);
  Real time = *m_time;
  Real exposedElectrodeVoltage=0;

  
  Real xR = m_maxvalXYZ1 + m_momentum_frequency*time;
  
  if(m_Theta == 1e0) 
    exposedElectrodeVoltage=exp(min(xR,0e0)) + max(xR,0e0);
  else
    if(m_Ehat == 0e0)
      {
	Real eXR = exp(xR);
	Real seXR= sqrt(eXR);
	exposedElectrodeVoltage=(-eXR + seXR*sqrt(4e0 + eXR) + 4e0*asinh( seXR/2e0 ))/2e0;
      }
    else
      MayDay::Error("Ehat >0 not implemented in IonizationWaveEBDirBcFunc");
  exposedElectrodeVoltage=2*m_maxvalXYZ1 * m_E0;


	
  Real value;
  if (a_point[0]>=m_xloc)
      value=exposedElectrodeVoltage;
  else 
    // LHS electrode - MD this is only for direct discharge 
    value = 0;

  
  return value;
}
