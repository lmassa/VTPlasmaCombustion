
#ifndef _BOUNDARYLAYER_H_
#define _BOUNDARYLAYER_H_

#define N_ELEMENTS  5  // number of state vector elements

#include "Vector.H"
#include "REAL.H"
#include "rksuite.H"

class BoundaryLayer
{

public:

  // weak construction is bad
  BoundaryLayer()
  {
     m_isDefined = false;
  }

  ///
  /**
   */
  BoundaryLayer(const Real a_rm, const Real a_ga, const Real a_Tinf, const Real a_pr, const Real a_boundLoc, const Real a_Rex, const Real a_ReCH);
  void define(const Real a_rm, const Real a_ga, const Real a_Tinf, const Real a_pr, const Real a_boundLoc, const Real a_Rex, const Real a_ReCH);
  int nvars() const;
  int numPts() const;
  int densityIndex() const;
  int momentumIndex() const;
  int energyIndex() const;
  Real getX0() const;
  Real getReLRoot2() const;
  void getVars(Vector<Real> & a_eta, Vector<Vector<Real> >& a_BLcons);
  const Vector<Real>& exportEta() const;
  const Vector<Vector<Real> >& exportCons() const;
  const Vector<Vector<Vector<Real> > >& exportCoeffs() const;
  void QUINAT(Vector<Real>& a_eta, Vector<Real>& a_fun, Vector<Vector<Real> >& coeffs);

  // Real getViscosity(Real y[N_ELEMENTS] );
  //void computeDerivatives(Real time,Real y[N_ELEMENTS], Real fr[N_ELEMENTS] );
  Real thetaInt(Real y[N_ELEMENTS] );
  Real thetaFun(int n, Vector<Real> y);
  void y2cons(Real& a_rho, Real& a_momx, Real& a_momy, Real& a_eng,   const Real& a_eta, const Real a_y[N_ELEMENTS]);
  

  Real mp_boundLoc, m_thetaMom;
private:
  Real m_rm, m_ga, m_pr, m_Tinf, m_rmuinf, m_xend, m_Lstar, m_ReRoot2, m_X0;
  Vector<Real> m_y, m_yend;
  Vector<Vector<Vector<Real> > > m_BLspline;
  Vector<Vector<Real> > m_BLcons;
  Vector<Real> m_eta;
  int m_nvars, m_ny;
  int BL_RHO, BL_MOMX, BL_MOMY, BL_ENG;
  Parameters m_BLparams;
  bool m_isDefined;
};

#endif
