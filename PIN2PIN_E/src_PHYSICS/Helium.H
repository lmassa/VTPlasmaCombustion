#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _Helium_H_
#define _Helium_H_

#include  <iostream>

#include "gasModel.H"
#include "NamespaceHeader.H"

class Helium : public gasModel
{
public:

  virtual void setPlasmaParameters(plasmaParameters& a_ppars, Vector<Real>& a_minData, Vector<Real>& a_maxData,  string& a_bgSpec, Vector<string>& a_charSpec);
};

#include "NamespaceFooter.H"
#endif
