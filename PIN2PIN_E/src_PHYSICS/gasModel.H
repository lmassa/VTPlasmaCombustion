#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _GASMODEL_H_
#define _GASMODEL_H_

#include "PolyGeom.H"
#include "LevelData.H"
#include "EBISLayout.H"
#include "RealVect.H"

#include "NamespaceHeader.H"




inline Real getDistanceFromAxis(const RealVect& a_xloc,
                         const RealVect& a_axis,
                         const RealVect& a_orig)
{
  RealVect trueNorm, closestPt;
  PolyGeom::pointToLine(closestPt, trueNorm, a_xloc,
                        a_orig, a_axis);
  RealVect distVec = a_xloc - closestPt;
  Real distance=1.0;
  PolyGeom::unifyVector(distVec, distance);
  return distance;
}
// To test it, in poisson driver you simply type

// PlasmaPhysics obj; // here (in the constructor) the parameters are read

class plasmaParameters
{
public:
  plasmaParameters() : isConstTe(true), isConstBg(true), isConstMobDiff(0), isDimensionless(false), isStiff(false), constMobIon(0), constDiffIon(0), initNe(1e-12), charToSpecThrFac(1.0) {}
  Vector<Vector<Real> > reacRateCoef;
  Vector<Vector<Real> > transpCoef;
  Vector<Vector<Real> > reactant;
  Vector<Vector<Real> > netReac;
  Vector<Real> secEmElec, secEmDiel;
  Vector<Real> actReac;
  Vector<Real> secEmEnergy;
  Vector<Real> eInelEn;
  Vector<Real> gInelEn;
  Vector<Real> chargeNr;  
  Vector<Real> sigebg, sigeo2, stickFac2bg, stickFac2o2;  
  Vector<int> species;
  Vector<int> reactions;
  Vector<Real> mass;
  Real kb, q, NA, etaT, mol2pm3;
  int porder;
  //LM:: New parameters added @VT
  bool isConstTe, isConstBg, isConstMobDiff, isDimensionless, isStiff;
  Real constMobEle, constDiffEle, constAlpha;
  Real constMobIon, constDiffIon;
  vector<bool> isNeutral;
  Real initNe, charToSpecThrFac;
  int eIndex, iIndex;
};




class gasModel
{
public:
  
  virtual ~gasModel() {}
  virtual void setPlasmaParameters(plasmaParameters& a_ppars, Vector<Real>& a_minData, Vector<Real>& a_maxData,  string& a_bgSpec, Vector<string>& a_charSpec)=0;

  virtual void initialize(LevelData<EBCellFAB>&       a_conState,
			  const plasmaParameters&     a_ppars,
			  const EBISLayout&           a_ebisl,
			  const Real&                 a_dx,
			  const RealVect&             a_domainLength) {}
  
  virtual Real setImplicitTimeStep() const
  {return 0.0;}

  virtual void setWallBCinfo(Vector<int>& a_wallBCInfo) const
  {}

};
#include "NamespaceFooter.H"
#endif
