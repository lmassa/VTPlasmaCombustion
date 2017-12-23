#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif


#ifndef _EBPIN2PINIBCFACTORY_H_
#define _EBPIN2PINIBCFACTORY_H_

#include "EBPin2PinIBC.H"
#include "EBCellFAB.H"
#include "EBFaceFAB.H"
#include "REAL.H"
#include "LevelData.H"
#include "ProblemDomain.H"
#include "EBPhysIBC.H"
#include "EBPin2PinIBC.H"
#include "EBPhysIBCFactory.H"

///
/**
 */
class EBPin2PinIBCFactory: public EBPhysIBCFactory
{
public:

  EBPin2PinIBCFactory();
  ///
  ~EBPin2PinIBCFactory();

  ///
  EBPin2PinIBCFactory(const Real&     a_gamma,
                          const Real&     a_ms,
                          const Real&     a_center,
                          const int&      a_shocknorm,
                          const bool&     a_shockbackward,
                          const bool&     a_doRZCoords,
			  const BoundaryLayer* const a_BLPtr);

  void define(const Real&     a_gamma,
	      const Real&     a_ms,
	      const Real&     a_center,
	      const int&      a_shocknorm,
	      const bool&     a_shockbackward,
	      const bool&     a_doRZCoords,
	      const BoundaryLayer* const a_BLPtr);
  ///
  EBPhysIBC* create() const;

protected:
  Real m_gamma;
  Real m_ms;
  Real m_center;
  int  m_shocknorm;
  bool m_shockbackward;
  bool m_doRZCoords;
  const BoundaryLayer* m_BLPtr;
private:

  void operator=(const EBPin2PinIBCFactory& a_input)
  {
    MayDay::Error("Invalid operator");
  }

  EBPin2PinIBCFactory(const EBPin2PinIBCFactory& a_input)
  {
    MayDay::Error("Invalid operator");
  }


};

#endif
