#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBPATCHSPECIESFACTORY_H_
#define _EBPATCHSPECIESFACTORY_H_

#include "EBSpeciesMacros.H"
#include "EBPatchSpecies.H"
#include "EBPatchGodunovFactory.H"
#include "PlasmaPhysics.H"

///
/**
 */
///
/**
 */
class EBPatchSpeciesFactory: public EBPatchGodunovFactory
{
public: 
  EBPatchSpeciesFactory();
  ///
  /**
   */
  virtual ~EBPatchSpeciesFactory();
  
  void define(const EBPhysIBCFactory* a_bcFactoryPtr,
	      PlasmaPhysics*    a_PlasmaPhysics,
	      const bool&             a_useLimiting,
	      const int&              a_nSpec,
	      const Real&             a_maxVal,
	      const Real&             a_minVal,
	      const bool&             a_setMaxMin);
  ///
  /**
   */
  EBPatchSpecies* create() const;

  ///
  /**
   */
  EBPatchSpeciesFactory(const EBPhysIBCFactory*  a_bc,
			PlasmaPhysics*     a_PlasmaPhysics,
			const bool&              a_useLimiting,
			const int&               a_nSpec,
			const Real&              a_maxVal =-1.e99,
			const Real&              a_minVal = 1.e99,
			const bool&              a_setMaxMin = false);

  PlasmaPhysics*     m_PlasmaPhysics;// unprotected
protected:
  const EBPhysIBCFactory*  m_bcFactoryPtr;
  bool                     m_useLimiting;
  bool                     m_setMaxMin;
  Real                     m_maxVal;
  Real                     m_minVal;
  Real                     m_nSpec;

private:

  void operator=(const EBPatchSpeciesFactory& a_input)
  {
    MayDay::Error("invalid operator");
  }
  EBPatchSpeciesFactory(const EBPatchSpeciesFactory& a_input)
  {
    MayDay::Error("invalid operator");
  }
};
#endif
