#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "ParmParse.H"
#include "CH_HDF5.H"
#include "parstream.H"

//Physics classes
#include "AMRIO.H"
#include "LevelDataOps.H"
#include "Vector.H"
#include "FArrayBox.H"
#include "DisjointBoxLayout.H"
#include "Box.H"
#include "FineInterp.H"
#include "RVect.H"
#include "IVect.H"

//addtional classes to test PETSC
#include <petscmat.h>
#include <cstdarg>
#include "Tuple.H"
#include "BaseEBFaceFAB.H"
#include "BaseEBCellFAB.H"
#include "EBLGIntegrator.H" //for CNUM
#include "EBAMRPlasma.H" //for CNUM


#include <iostream>
#include <fstream>
#include <deque>

#define CHECK_FLOATING_PT 0

#if CHECK_FLOATING_PT==1
#include <fenv.h>
#endif

using std::ifstream;
using std::ios;


/***************/
int
main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  MPI_Init(&a_argc,&a_argv);
#endif
  PetscInitializeNoArguments();
  { //scoping trick
#if CHECK_FLOATING_PT==1
  //  int except =  FE_DIVBYZERO | FE_UNDERFLOW | FE_OVERFLOW |  FE_INVALID ;
  int except =  FE_DIVBYZERO |  FE_OVERFLOW |  FE_INVALID ;
  feenableexcept(except);
#endif

 // Check for an input file
  char* inFile = NULL;

  if (a_argc > 1)
    {
      inFile = a_argv[1];
    }
  else
    {
      pout() << "Usage: <executable name> <inputfile>" << endl;
      pout() << "No input file specified" << endl;
      return -1;
    }
  // Parse the command line and the input file (if any)
  ParmParse pp(a_argc-2,a_argv+2,NULL,inFile);
  int maxSteps = 1;
  pp.get("max_step", maxSteps);
  Real maxTime = 1;
  pp.get("max_time", maxTime);
  
  EBAMRPlasma solver;
  solver.runPlasmaOnly(maxTime,maxSteps,"IW");
  
  }
#ifdef CH_MPI
  CH_TIMER_REPORT();
  dumpmemoryatexit();
  MPI_Finalize();
#endif
}
