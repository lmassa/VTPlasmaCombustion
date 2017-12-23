#include <math.h>

#include <iostream>
#include <stdlib.h>
#include "rksuiteI.h"
#include "gauss_legendre.h"
#include "BoundaryLayer.H"


// helper functions for exponentiation to integer powers
#define SQR(x) ((x)*(x))
#define P2(x) ((x)*(x))
#define P3(x) (P2(x)*(x))
#define P4(x) (P2(x)*P2(x))
#define P5(x) (P3(x)*P2(x))
#define P7(x) (P5(x)*P2(x))
//
static Real getViscosity(Real y[N_ELEMENTS],const Real rm,const Real ga,const Real tinf, const Real rmuinf) {

//For Helium:

  const Real  a = 1.4844e-6;//1.458e0;
  const Real  b = 79.4; //110.4e0;
  const Real  z = 8.1039e-8; //.693873e-6/1.0e-5;

  Real rm2 = rm*rm;
  Real t = 1.0e0 + (ga-1.0e0)*rm2*y[2]/2.00;
  Real tstar = fabs(t*tinf);
  Real rmu;

  if (tstar >= b)
    rmu = (a*pow(tstar,1.5e0)/(tstar+b))/rmuinf;
  else
    rmu = z*tstar/rmuinf;


  return rmu;
}
static void computeDerivatives(Real time,
			       Real y[N_ELEMENTS],
			       Real fr[N_ELEMENTS],
			       Parameters& rkparams) {

  Real rm = rkparams.R1;
  Real ga = rkparams.R2;
  Real pr = rkparams.R3;
  Real tinf = rkparams.R4;
  Real rmuinf = rkparams.R5;

  Real rm2 = rm*rm;
  Real t = 1.0e0 + (ga-1.0e0)*rm2*y[2]/2.00;
  Real rmu = getViscosity(y, rm, ga,  tinf,  rmuinf);

  
  if (! ( fabs(rmu) < 1e5)) exit(1);

  //y[0] = u/u_inf
  //y[1] = mu/mu_inf d(u/u_inf)/d(eta)
  //y[2] = (T-T_inf)/(T0-T_inf)
  //y[4] = g

  fr[0] = y[1]/rmu;
  fr[1] = -y[4]*y[1]/rmu;
  fr[2] = pr*y[3]/rmu;
  fr[3] = -pr*y[4]*y[3]/rmu - 2.0e0*y[1]*y[1]/rmu;
  fr[4] = 0.5e0*y[0]/t;
}


Real BoundaryLayer::thetaInt(Real y[N_ELEMENTS] ) {


  Real rm2 = m_rm*m_rm;
  Real t = 1.0e0 + (m_ga-1.0e0)*rm2*y[2]/2.00;
  
  Real fun = y[0]/t *(1-y[0]);

  return fun;
}


Real BoundaryLayer::thetaFun(int n, Vector<Real> y) {

  Real xstart = 0;
  Real *mp[2] = {NULL,NULL};
  Real* w = NULL;
  Real result=0;
  int m, dtbl;
  

  gauss_points( m, mp[0], mp[1], w, dtbl, n, xstart, m_xend);


  Real thres[N_ELEMENTS] = {0};
  Real yl[N_ELEMENTS] = {0};
  Real yp[N_ELEMENTS] = {0};
  Real ymax[N_ELEMENTS] = {0};
  int   cflag;  // status return 
  bool isodd = n&1;
  bool broken = false;

  RKSUITE    rksuite;
  for (int k = 0; k < N_ELEMENTS; k++){
    thres[k] = 1e-12;yl[k] = y[k];}
  Real TOL = 1e-12; 
  int    method = 3;  // RK(6,7)

  
  Real tgot = xstart;
  computeDerivatives( tgot, yl, yp, m_BLparams );
  for (int itime = 0; itime < 2; itime++)
    {  
      Real* x = mp[itime];
      int ii = (itime>0) ? (m-1) : 0;
      Real twant = x[ii];
      rksuite.setup(N_ELEMENTS, tgot, yl, twant, TOL, thres, method, "UT", false, 0.0f, false, m_BLparams);
      int j;
      for (int i = 0; i < m; i++)
	{
	  if (itime == 0)
	    j = m-1-i;
	  else
	    j=i;

	  cflag = 0;
	  twant = x[j];
	  if ( fabs(twant - tgot) > 1e-12)
	    {
	      do{
		rksuite.ut( computeDerivatives, twant, tgot, yl, yp, ymax, cflag );
	      }while (cflag>= 2 && cflag <=4 );
	    }
	  if (cflag >= 5) {
	    std::cout << cflag << ", "<< TOL << " UT " << tgot <<", " << yl[0] << ", " << yl[1] << ", " <<yl[2] << std::endl;
	    broken = true;
	    break;}
	  if (itime == 0 || (! isodd) || j > 0 )
	    result += w[j]*thetaInt(yl);
	}
    }
	
  for (int k = 0; k < 2; k++)
    free(mp[k]);
  if (dtbl) free(w);
  return result*0.5*(m_xend-xstart);

}

Real BoundaryLayer::getReLRoot2() const
{
  return m_ReRoot2;
}

Real BoundaryLayer::getX0() const
{
  return m_X0;
}

int BoundaryLayer::nvars() const
{
  return m_nvars;
}
int BoundaryLayer::numPts() const
{
  return m_ny;
}
int BoundaryLayer::densityIndex() const
{
  return BL_RHO;
}
int BoundaryLayer::momentumIndex() const
{
  return BL_MOMX;
}
int BoundaryLayer::energyIndex() const
{
  return BL_ENG;
}


const Vector<Real>& BoundaryLayer::exportEta() const
{
  return m_eta;
}
const Vector<Vector<Real> >& BoundaryLayer::exportCons() const
{
  return m_BLcons;
}
const Vector<Vector<Vector<Real> > >& BoundaryLayer::exportCoeffs() const
{
  return m_BLspline;
}

void BoundaryLayer::y2cons(Real& a_rho, Real& a_momx, Real& a_momy, Real& a_eng,  const Real& a_eta, const Real a_y[N_ELEMENTS]) 
{
  Real rm2 = P2(m_rm);
  Real t = 1.0e0 + (m_ga-1.0e0)*rm2*a_y[2]/2.00;
  Real dimFact = m_ga*m_rm;
  a_rho = m_ga/t;
  Real rhoU =a_y[0]/t;
  Real rhoV = (a_eta*rhoU/2.0 - a_y[4]);
  a_momx = dimFact*rhoU;
  a_momy = dimFact*rhoV;
  a_eng = 1.0/(m_ga-1); //+ (P2(a_momx)+P2(a_momy))/(2.0*a_rho);
  
  
  //std::cout <<a_eta << ", " << rhoU << std::endl;
  //Real yp[N_ELEMENTS] = {0};
  //Real yl[N_ELEMENTS] = {0};
  //for(int i=0;i<N_ELEMENTS;i++) yl[i]=a_y[i];
  //computeDerivatives( 0.0, yl, yp, m_BLparams );
  //std::cout <<a_eta << ", " << - yp[4]*a_y[1]*a_eta/(0.280 - m_X0)*m_rm*m_rm*m_ga << std::endl;
}
//
void BoundaryLayer::getVars(Vector<Real> & a_eta, Vector<Vector<Real> >& a_BLcons)
{
  
  // unravel the inp/out var
  Vector<Real>& rho = a_BLcons[BL_RHO];
  Vector<Real>& momx = a_BLcons[BL_MOMX];
  Vector<Real>& momy = a_BLcons[BL_MOMY];
  Vector<Real>& eng = a_BLcons[BL_ENG];

  int ny = a_eta.size();
  Real xstart = 0;
  Real thres[N_ELEMENTS] = {0};
  Real yl[N_ELEMENTS] = {0};
  Real yp[N_ELEMENTS] = {0};
  Real ymax[N_ELEMENTS] = {0};
  int  cflag;  // status return 
  bool broken = false;
  RKSUITE    rksuite;
  for (int k = 0; k < N_ELEMENTS; k++){
    thres[k] = 1e-9;yl[k] = m_y[k];}
  Real TOL = 1e-12; 
  int    method = 3;  // RK(6,7)

  Real tgot = xstart;
  computeDerivatives( tgot, yl, yp, m_BLparams );
  Real twant = a_eta[ny-1];
  rksuite.setup(N_ELEMENTS, tgot, yl, twant, TOL, thres, method, "UT", false, 0.0f, false, m_BLparams );
  for (int i = 0; i < ny; i++)
    {
      cflag = 0;
      twant = a_eta[i];   // eta variable
      if (twant > m_xend)
	for (int k = 0; k < N_ELEMENTS; k++)
	  yl[k] = m_yend[k];
      else
	{
	  if ( fabs(twant - tgot) > 1e-12)
	    {
	      do{
		rksuite.ut( computeDerivatives, twant, tgot, yl, yp, ymax, cflag );
	      }while (cflag>= 2 && cflag <=4 );
	    }
	  if (cflag >= 5) {
	    std::cout << cflag << ", "<< TOL << " UT " << tgot <<", " << yl[0] << ", " << yl[1] << ", " <<yl[2] << std::endl;
	    broken = true;
	    break;}
	}
      y2cons(rho[i],momx[i],momy[i],eng[i],a_eta[i],yl);
    }
	
}


void BoundaryLayer::QUINAT(Vector<Real>& a_eta, Vector<Real>& a_fun, Vector<Vector<Real> >& coeffs)
{

  int NPT = a_eta.size()-1;
  if (NPT<=1) { std::cout << " error NPT <=1"; return;}


  // unravel the inp/out var
  Vector<Real>& X = a_eta;
  Vector<Real>& Y = a_fun;
  Vector<Real>& B = coeffs[0];
  Vector<Real>& C = coeffs[1];
  Vector<Real>& D = coeffs[2];
  Vector<Real>& E = coeffs[3];
  Vector<Real>& F = coeffs[4];


  int M = NPT - 2;
  Real Q = X[1] - X[0];
  Real R = X[2] - X[1];
  Real Q2 = Q*Q;
  Real R2 = R*R;
  Real QR = Q + R;
  Real P, P2, PQ, PR, P3, Q3, S, T, U, V;
  int I;
  D[0] = 0e0;
  E[0] = 0e0;
  D[1] = 0e0;
  if (Q != 0e0) D[1] = 6e0*Q*Q2/(QR*QR);

  if (M>=1) {
    for (I=1;I<=M;I++)
    {
      P = Q;
      Q = R;
      R = X[I+2] - X[I+1];
      P2 = Q2;
      Q2 = R2;
      R2 = R*R;
      PQ = QR;
      QR = Q + R;
      if (Q == 0.0){D[I+1] = 0e0;E[I] = 0e0;F[I-1] = 0e0;continue;}
      Q3 = Q2*Q;
      PR = P*R;
      Real PQQR = PQ*QR;
      D[I+1] = 6e0*Q3/(QR*QR);
      D[I] += (Q+Q)*(15e0*PR*PR+(P+R)*Q*(20e0*PR+7e0*Q2)+Q2*(8.*(P2+R2)+21.*PR+Q2+Q2) )/(PQQR*PQQR);
      D[I-1] += 6e0*Q3/(PQ*PQ);
      E[I] = Q2*(P*QR+3e0*PQ*(QR+R+R))/(PQQR*QR);
      E[I-1] = E[I-1] + Q2*(R*PQ+3e0*QR*(PQ+P+P))/(PQQR*PQ);
      F[I-1] = Q3/PQQR;
      }
  }

  if (R!=0e0) D[M] += 6e0*R*R2/(QR*QR);





  for(I=1;I<=NPT;I++)
  {
    if (X[I]==X[I-1]){
      B[I] = Y[I];
      Y[I] = Y[I-1];
    }
    else
      B[I] = (Y[I]-Y[I-1])/(X[I]-X[I-1]);
  }
  for(I=2;I<=NPT;I++)
    {
    if (X[I] == X[I-2]) 
    {
    C[I] = B[I]*0.5e0;
    B[I] = B[I-1];
    }
    else
      C[I] = (B[I]-B[I-1])/(X[I]-X[I-2]);
      }



      if (M >= 1) 
       {
       P = 0e0;
       C[0] = 0e0;
       E[M] = 0e0;
       F[0] = 0e0;
       F[M-1] = 0e0;
       F[M] = 0e0;
       C[1] = C[3] - C[2];
       D[1] = 1e0/D[1];
       }

       if (M >= 2){
       for(I=2;I<=M;I++)
       {
       Q = D[I-1]*E[I-1];
       D[I] = 1e0/(D[I]-P*F[I-2]-Q*E[I-1]);
       E[I] = E[I] - Q*F[I-1];
       C[I] = C[I+2] - C[I+1] - P*C[I-2] - Q*C[I-1];
       P = D[I-1]*F[I-1];
       }
       }

       I = NPT - 1;
       C[NPT-1] = 0e0;
       C[NPT] = 0e0;
       if (NPT >= 3) {
       for (M=3;M<=NPT;M++){
       I = I - 1;
       C[I] = (C[I]-E[I]*C[I+1]-F[I]*C[I+2])*D[I];
       }
       }



       M = NPT - 1;
      Q = X[1] - X[0];
      R = X[2] - X[1];
      Real B1 = B[1];
      Q3 = Q*Q*Q;
      QR = Q + R;
      if (QR ==0.) {
      V = 0e0;
      T = 0e0;
      }
      else
      { 
	V = C[1]/QR;
	T = V;
      }
      F[0] = 0e0;
      if (Q != 0e0) F[0] = V/Q;
      for(I=1;I<=M;I++){
	P = Q;
	Q = R;
	R = 0e0;
	if (I!=M) R = X[I+2] - X[I+1];
	P3 = Q3;
	Q3 = Q*Q*Q;
	PQ = QR;
	QR = Q + R;
	S = T;
	T = 0e0;
	if (QR!=0e0) T = (C[I+1]-C[I])/QR;
	U = V;
	V = T - S;
         if (PQ==0) {
	   C[I] = C[I-1];
	   D[I] = 0e0;
	   E[I] = 0e0;
	   F[I] = 0e0;
	}
         else
	   {
	     F[I] = F[I-1];
	     if  (Q != 0e0) F[I] = V/Q;
	     E[I] = 5e0*S;
	     D[I] = 10e0*(C[I]-Q*S);
	     C[I] = D[I]*(P-Q) + (B[I+1]-B[I]+(U-E[I])*P3-(V+E[I])*Q3)/PQ;
	     B[I] = (P*(B[I+1]-V*Q3)+Q*(B[I]-U*P3))/PQ - P*Q*(D[I]+E[I]* (Q-P) );
									      }
      }

      P = X[1] - X[0];
      S = F[0]*P*P*P;
      E[0] = 0e0;
      D[0] = 0e0;
      C[0] = C[1] - 10e0*S;
      B[0] = B1 - (C[0]+S)*P;

      Q = X[NPT] - X[NPT-1];
      T = F[NPT-1]*Q*Q*Q;
      E[NPT] = 0e0;
      D[NPT] = 0e0;
      C[NPT] = C[NPT-1] + 10e0*T;
      B[NPT] = B[NPT] + (C[NPT]-T)*Q;
}
     

BoundaryLayer::BoundaryLayer(const Real a_rm, const Real a_ga, const Real a_Tinf, const Real a_pr, const Real a_boundLoc, const Real a_Rex, const Real a_ReCH) {

  define(a_rm, a_ga, a_Tinf, a_pr, a_boundLoc, a_Rex, a_ReCH);
}

void BoundaryLayer::define(const Real a_rm, const Real a_ga, const Real a_Tinf, const Real a_pr, const Real a_boundLoc, const Real a_Rex, const Real a_ReCH) {

  m_isDefined = true;
  m_rm =a_rm;
  m_ga = a_ga;
  m_pr = a_pr;
  m_Tinf = a_Tinf; 
  m_rmuinf = 1.4844e-6*(pow(m_Tinf,1.5e0))/(m_Tinf+79.4);
  m_y.resize(N_ELEMENTS,0.0);
  m_yend.resize(N_ELEMENTS,0.0);
  mp_boundLoc = a_boundLoc;
  m_nvars = 4; // always 2dimensional
  Real ReL = a_ReCH*m_ga*m_rm;
  BL_RHO=0;BL_MOMX=1;BL_MOMY=2;BL_ENG=3;
  m_BLparams = (Parameters){a_rm,a_ga,a_pr,m_Tinf,m_rmuinf};
  RKSUITE rksuite;  // create the integrator object instance

  Real epx = 1.e-3;

  Real vtarget[] = {1e0,0e0};
  Real vin[] = {1e0 , 1e0};
  Real pert[3][2] = { {epx,0}, {0,epx}, {0,0}};
  Real vout[3][2] = {0};
  int jcomp[] = {0,2};

  Real resid = 100;
  int icount = 0;

  Real xstart = 0;
  m_xend = 120;
  Real xstop;
  if (m_xend > 7)  xstop = 7; else xstop = m_xend; 

  while (resid > 1.e-9 || xstop < m_xend)
    {
      icount = icount+1;
      if(icount> 100) break;
      Real yp[N_ELEMENTS] = {0}; 
      int   cflag;  // status return  

      if(resid <  1e-5) 
	{
	  xstop = 1.2*xstop;xstop = min(xstop,m_xend); icount =0;
	}
      //xspan = [ xstart xstop];
      for (int i = 0; i <= 2; i++) {
    
	Real y[N_ELEMENTS] = {0};
	for (int k = 1; k <= 2; k++)
	  y[k] = vin[k-1] + pert[i][k-1];
	Real ymax[N_ELEMENTS] = {0};
 
	Real TOL;
	if(xstop < m_xend)
	  {
	    Real thres[N_ELEMENTS] = {0};
	    for (int k = 0; k < N_ELEMENTS; k++)
	      thres[k] = 1e-7;
	    TOL = 1e-6; 
	    int    method = 2;  // RK(4,5)
	    rksuite.setup(N_ELEMENTS, xstart, y, m_xend, TOL, thres, 
			  method, "UT", false, 0.0f, false, m_BLparams );
	  }
	else
	  {
	    Real thres[N_ELEMENTS] = {0};
	    for (int k = 0; k < N_ELEMENTS; k++)
	      thres[k] = 1e-12;
	    TOL = 1e-12; 
	    int    method = 3;  // RK(6,7)
	    rksuite.setup(N_ELEMENTS, xstart, y, m_xend, TOL, thres, 
			  method, "UT", false, 0.0f, false, m_BLparams );
	  }
	Real tgot = xstart;
	Real twant = xstop;
	computeDerivatives( tgot, y, yp, m_BLparams );
	cflag = 0;
	do {
	  rksuite.ut( computeDerivatives, twant, tgot, y, yp, ymax, cflag );}  while (cflag>= 2 && cflag <=4 );
	
	if (cflag > 5) {
	  std::cout << cflag << ", "<< TOL << " UT " << tgot <<", " << y[0] << ", " << y[1] << ", " <<y[2] << std::endl;	
	  break;}

	for (int k = 0; k <= 1; k++)
	  vout[i][k] = y[jcomp[k]]- vtarget[k];
	if (i==2){ for (int k = 0; k < N_ELEMENTS; k++) m_yend[k] = y[k];}
      }


      Real a = vout[0][0] - vout[2][0];Real c = vout[0][1] - vout[2][1];Real b = vout[1][0] - vout[2][0];Real d = vout[1][1] - vout[2][1];Real r = vout[2][0]*epx;Real s = vout[2][1]*epx;
      Real deltav[] = {(d*r - b*s)/(b*c - a*d),(c*r - a*s)/(-(b*c) + a*d)};
      for (int k = 0; k <= 1; k++)
	vin[k] += deltav[k];


      resid = sqrt(P2(vout[2][0]) + P2(vout[2][1])); 
    }

  //std::cout << " Resid " << resid << ", " << vin[0] << ", " <<vin[1] << std::endl;
  for (int k = 1; k <= 2; k++)
    m_y[k] = vin[k-1];


  m_X0 = -a_Rex/ReL; //negative because it is the origin wrt to the inflow
  m_Lstar = sqrt(a_Rex)/ReL;  //delta in MD worksheet
  m_ReRoot2 = sqrt(ReL); // assuming m_thetaAssign is given at x=0


  m_thetaMom = thetaFun( 100, m_y)*m_Lstar;
  //std::cout <<100<< " thetaOlstar " << m_thetaMom << ", Lstar= " << m_Lstar  << std::endl;


  //std::cout << " X0 " << m_X0 << std::endl;
  //std::cout << " ReL " << ReL << std::endl;

  m_ny = 300;

  m_eta.resize(m_ny,0.0);
  Real deta = (m_xend - 0.0)/Real(m_ny-1.0);
  for (int i = 0; i < m_ny; i++) m_eta[i] = i*deta;

  m_BLspline.resize(m_nvars);
  m_BLcons.resize(m_nvars);
  for (int i = 0; i < m_nvars; i++)
    {
      m_BLcons[i].resize(m_ny,0.0);
    }
  getVars(m_eta, m_BLcons);


  for (int i = 0; i < m_nvars; i++)
    {
      m_BLspline[i].resize(5); // 5==spline order
      for (int k = 0; k < 5; k++) m_BLspline[i][k].resize(m_ny,0.0);
      QUINAT(m_eta,m_BLcons[i],m_BLspline[i]);
    }


}
