#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBAMRPLASMA_H_
#define _EBAMRPLASMA_H_

#include "AMR.H"
#include "EBAMR.H"
#include "AMRMultiGrid.H"
#include "LevelData.H"
#include "EBCellFAB.H"
#include "EBISLayout.H"
#include "EBCoarseAverage.H"
#include "EBSimpleSolver.H"
#include "UsingNamespace.H"
#include "EBViscousTensorOpFactory.H"
#include "EBConductivityOpFactory.H"
#include "EBAMRNavierFactory.H"
#include "BoundaryLayer.H"
#include "EBPatchPolytropicFactory.H"
#include "EBPin2PinIBCFactory.H"
//species stuff
#include "EBSpeciesIBCFactory.H"
#include "EBPatchSpeciesFactory.H"
#include "EBAMRSpecies.H"
#include "ElectricPotential.H"

///
/**
  Agglomeration of the many parameters of an AMR calc
*/
class PlasmaParameters
{
public:
  IntVect       nCells;
  int           maxGridSize;
  int           blockFactor;
  int           bufferSize;
  Real          fillRatio;
  int           ebBcType;
  int           domBcType;
  int           ebBcTypeT;
  Real          ebBCValueT;
  int           orderEB;
  int           maxLevel;    //maxlevel allowed in the grid
  int           finestLevel;  // not necessary identical to max level in case somelevel are not populated
  int           numLevels;
  Vector<int>   refRatio;
  ProblemDomain coarsestDomain;
  RealVect      coarsestDx;
  RealVect      domainLength;
  RealVect      probLo;
  RealVect      probHi;
  IntVect       ghostPhi;
  IntVect       ghostRHS;
  bool          noRefCorners;
  void coarsen(int a_factor);
  void  refine(int a_factor);

};


class Sutherland
{
public:
  Sutherland(Real tinfinity);
  Real t0;
  Real a;
  Real b;
  Real tinf;
  Real muinf;
  Real exp;
  Real factor;
};

/// Class to manage non-subcycled advance for EBAMRINS.
/**
*/
class EBAMRPlasma
{

public:


  EBAMRPlasma();




  void definegrids();


  void run(Real a_max_period, int a_dumps_per_period, const char* a_inFile);

  void regrid(const bool a_isPlasmaOnly=false);
  void printGridSize();
  void regridRhos(Vector<LevelData<EBCellFAB>* >& a_rhosData, Vector<DisjointBoxLayout>& a_oldDBL, const int& flag);
  void getVelTemp(const bool& a_dimensionalScaling =true);
  void getSpecies();
  void getPlasmaParameters(PlasmaParameters&  a_params, bool a_forceSingleLevel = false);
  
  
 
  void updateGhost(Vector<LevelData<EBCellFAB>* >&  a_phi, const Real& a_limit);

  void tagCells(IntVectSet& a_tags, const int & a_lev);

  void writeCheckpointFile(const char* a_inFile) const;

  void outputHDF(const Vector<LevelData<EBCellFAB>* >& a_Ptr, const string& a_name, const int& a_ival);

  Real currentTime() const;
  void stdOut(const Real& a_EPresid) const;

  void conclude() const;

  //LMadd@VT
  void runPlasmaOnly(Real a_max_time, int a_max_step, const char* a_inFile);
  inline void grabPotentials() {
    m_phi = m_ElectricPotential.getPhi();
    m_gradPhi = m_ElectricPotential.getGradPhi();
    m_cellGradPhi = m_ElectricPotential.getCellGradPhi();
    m_IVGradPhi = m_ElectricPotential.getIVGradPhi();
  }
  

private:
  //factories for the NS equations
  EBPin2PinIBCFactory m_bcGammaFactory;
  EBPatchPolytropicFactory m_patchGammaFactory;
  //factories for the species equations
  EBSpeciesIBCFactory m_bcSpeciesFactory;
  EBPatchSpeciesFactory m_patchSpeciesFactory;


  PlasmaPhysics m_PlasmaPhysics;
  
  // two AMR objects to split the NS and species computations
  EBAMR m_amr;
  EBAMR m_amrS; 

  //stuff
  PlasmaParameters m_params;
  Vector<ProblemDomain>  m_domain;
  Vector<Real>           m_dx;
  AMRMultiGrid<LevelData<EBCellFAB> > m_solver, m_solverT;
  int m_verbosity;
  int m_ghost_Dk;
  Real m_alpha, m_beta, m_gamma;
  int m_tagShrinkDomain, m_tagVortComps, m_maxRefLevel, m_maxRefLevelS;
  int m_RegridSteps, m_plot_interval, m_numGhost, m_nCompState, m_gridRedefineSteps, m_checkpointInterval;
  Real m_refineThreshold, m_tagBuffer, m_lambdafac, m_SrefineThresh;
  bool m_inviscid, m_isFunctionEBVTO, m_isVTSdefine, m_viscousOutput;

  //BL:: define memory container
  BoundaryLayer m_BL;

  //SFD:: define variables
  bool m_SFD;
  Real m_chi, m_Dl, m_h, m_Cchi, m_chi_exp;
  //data held in memory  
  Vector<LevelData<EBCellFAB>* > m_qbar;
  Vector<LevelData<EBCellFAB>* > m_source;

  //grid operators
  Vector<RefCountedPtr<EBQuadCFInterp> >    m_quadCFI;
  Vector<RefCountedPtr<EBCoarseAverage> >   m_aveOper;
  Vector<RefCountedPtr<EBCoarseAverage> >   m_aveSpac;
  Vector<RefCountedPtr<EBCoarseAverage> >   m_aveStat;
  Vector<DisjointBoxLayout> m_grids;
  Vector<EBISLayout>        m_ebisl;
  Vector<EBLevelGrid>       m_eblg;
  Vector<EBAMRNavier*>     m_NavierLevels;
  Vector<EBAMRSpecies*>    m_SpeciesLevels;
  Vector<LevelData<EBCellFAB>* > m_cellNData;

  // diffusion coefficients for NavierStokes object
  Vector<RefCountedPtr<LevelData<EBCellFAB> > >           m_aco;
  Vector<RefCountedPtr<LevelData<EBFluxFAB> > >           m_eta;
  Vector<RefCountedPtr<LevelData<EBFluxFAB> > >           m_lambda;
  Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >    m_etaIrreg;
  Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >    m_lambdaIrreg;
  Vector<RefCountedPtr<LevelData<EBCellFAB> > >           m_rcv;

  //wall data
  Vector< RefCountedPtr<LevelData<BaseIVFAB<Real> > > > m_rhosRHS;

  //output streams
  ofstream     m_out ;

  //additional memory for species calculations
  //@@Spec@@//
  int m_nSpec;  // number of species to be read from input
  int m_finestLevel, m_finestLevelNS;

  //Potential solution
  ElectricPotential m_ElectricPotential;


  Vector<LevelData<EBCellFAB>* > m_Velo;
  Vector<LevelData<EBCellFAB>* > m_Temp;
  Vector<LevelData<EBCellFAB>* > m_phi;
  Vector<LevelData<EBFluxFAB>* > m_gradPhi;
  Vector<LevelData<EBCellFAB>* > m_cellGradPhi;
  Vector< LevelData<BaseIVFAB<Real>>* > m_IVGradPhi;
  Vector<LevelData<EBCellFAB>* > m_FSrc;
  Vector<LevelData<BaseIVFAB<Real> >* > m_rhos;

  Vector<EBPWLFineInterp* > m_ebFineInterps;

  //conversion factors
  Real m_Tinf,m_Vel0,m_Src0,m_SrcV;
  
  /***------------------Useless Stuff Below-------------------------***/
  
  //we disallow copy construction and assignment of objects
  //which contain pointered data
  EBAMRPlasma(const EBAMRPlasma& a_input)
  {
    MayDay::Error("invalid constructor");
  }

  //we disallow copy construction and assignment of objects
  //which contain pointered data
  void operator=(const EBAMRPlasma& a_input)
  {
    MayDay::Error("invalid constructor");
  }

};

#endif
