#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif


#ifndef _EBSPECIESIBCFACTORY_H_
#define _EBSPECIESIBCFACTORY_H_

#include "EBCellFAB.H"
#include "EBFaceFAB.H"
#include "REAL.H"
#include "LevelData.H"
#include "ProblemDomain.H"
#include "EBPhysIBC.H"
#include "EBPhysIBCFactory.H"
#include "EBSpeciesIBC.H"
#include "PlasmaPhysics.H"

///
/**
 */
class EBSpeciesIBCFactory: public EBPhysIBCFactory
{
public:
  ///
  EBSpeciesIBCFactory();

  ///
  ~EBSpeciesIBCFactory();


  void  define( const int&      a_nSpec, PlasmaPhysics*  a_PlasmaPhysics);

  ///
  EBPhysIBC* create() const;

protected:
  int m_nSpec;
  PlasmaPhysics*  m_PlasmaPhysics;
private:

  void operator=(const EBSpeciesIBCFactory& a_input)
  {
    MayDay::Error("Invalid operator");
  }

  EBSpeciesIBCFactory(const EBSpeciesIBCFactory& a_input)
  {
    MayDay::Error("Invalid operator");
  }


};
#endif
