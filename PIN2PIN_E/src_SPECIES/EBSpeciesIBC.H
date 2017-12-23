#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBSPECIESIBC_H_
#define _EBSPECIESIBC_H_

#include  <iostream>

#include "LevelData.H"
#include "FArrayBox.H"
#include "Vector.H"
#include "RealVect.H"
#include "AMRIO.H"

#include "EBPhysIBC.H"
#include "EBSolidF_F.H" // contains SLOPEBCS
#include "EBPatchSpeciesF_F.H" // contains SCALARFABMULT, SPECIESEXTRAPBC
#include "PlasmaPhysics.H"

class EBSpeciesIBC : public EBPhysIBC
{
public:
  virtual ~EBSpeciesIBC();

  EBSpeciesIBC(const Real&     a_nSpec, PlasmaPhysics*  a_PlasmaPhysics);

  ///
  void define(const ProblemDomain&  a_domain,
              const RealVect&       a_dx);

  ///  For every box in this level, this function is called.
  void fluxBC(EBFluxFAB&            a_flux,
              const EBCellFAB&      a_Wcenter,
              const EBCellFAB&      a_Wextrap,
              const Side::LoHiSide& a_sd,
              const Real&           a_time,
              const EBISBox&        a_ebisBox,
              const DataIndex&      a_dit,
              const Box&            a_box,
              const Box&            a_faceBox,
              const int&            a_dir);

  /// Initialize
  void initialize(LevelData<EBCellFAB>& a_conState,
                  const EBISLayout& a_ebisl) const;

  /// Set boundary slopes
  /**
     The boundary slopes in a_dW are already set to one sided difference
     approximations.  If this function doesn't change them they will be
     used for the slopes at the boundaries.
  */
  void setBndrySlopes(EBCellFAB&       a_deltaPrim,
                      const EBCellFAB& a_primState,
                      const EBISBox&   a_ebisBox,
                      const Box&       a_box,
                      const int&       a_dir);
  void setVelocities(const EBCellFAB& a_normalVel,
		     const EBFluxFAB& a_advectionVel); 

protected:
  bool m_isFortranCommonSet;
  bool m_isDefined;
  int m_nSpec;
  PlasmaPhysics* m_PlasmaPhysics;

  bool m_isVelSet;
  const EBCellFAB* m_normalVelPtr;
  const EBFluxFAB* m_advectionVelPtr;
  
  

  Real m_dx;
  ProblemDomain  m_domain;
  RealVect m_domainLength;

private:
  EBSpeciesIBC()
  {
    MayDay::Error("invalid operator");
  }

  //disallowed for all the usual reasons
  void operator=(const EBSpeciesIBC& a_input)
  {
    MayDay::Error("invalid operator");
  }
  EBSpeciesIBC(const EBSpeciesIBC& a_input)
  {
    MayDay::Error("invalid operator");
  }


};
#endif
