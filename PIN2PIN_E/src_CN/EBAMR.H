#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBAMR_H_
#define _EBAMR_H_


#include "AMR.H"
#include "DisjointBoxLayout.H"
#include "NamespaceHeader.H"

class EBAMR: public AMR
{

public:
  void setupForGivenDBLRun(const Vector<DisjointBoxLayout>& a_grids); 
  void setupForGivenICRun(const Vector<DisjointBoxLayout>& a_grids);
  void writePlotFile(const string& a_prefix, const int a_ival) const;
  void writeCheckpointData(HDF5Handle& a_handle) const;
  void setupForRestartS(HDF5Handle& a_handle);
  void setFinestLevel(int a_finest_level);
  void setTime(const Real& a_time);
};

#include "NamespaceFooter.H"
#endif
