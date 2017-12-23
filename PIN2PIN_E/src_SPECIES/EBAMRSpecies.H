#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBAMRSPECIES_H_
#define _EBAMRSPECIES_H_


#include "AMRLevel.H"
#include "EBCellFAB.H"
#include "BaseIVFAB.H"
#include "LevelData.H"
#include "EBLevelRedist.H"
#include "EBCoarToCoarRedist.H"
#include "EBCoarToFineRedist.H"
#include "EBFineToCoarRedist.H"
#include "EBCoarseAverage.H"
#include "EBPWLFineInterp.H"
#include "EBFluxRegister.H"
#include "EBLevelSpecies.H"
#include "Box.H"
#include "IntVectSet.H"
#include "Vector.H"
#include "DisjointBoxLayout.H"
#include "EBSpeciesFluxOpFactory.H"
#include "EBSpeciesFluxOp.H"
#include "BiCGStabSolver.H"
#include "EBLevelTGA.H"
#include "EBPatchSpecies.H"

#include "NamespaceHeader.H"

///interface for load balancing
typedef int(*LoadBalanceFunc)(Vector<int>&         a_procs,
                               const Vector<Box>&   a_boxes,
                               const ProblemDomain& a_domain,
                               bool a_verbose);

/// AMRLevel for advection-diffusion-reaction species integration
/**
 */
class EBAMRSpecies : public AMRLevel
{
public:

  /// Default constructor
  EBAMRSpecies();

/*********/
  ~EBAMRSpecies();
  void setPhysics(PlasmaPhysics*  a_PlasmaPhysics);
  void SetVelTempPhi(const Vector< LevelData<EBCellFAB>* >& a_Velo,
				   const Vector< LevelData<EBCellFAB>* >& a_Temp,
				   const Vector< LevelData<EBCellFAB>* >& a_phi);
  void SetVelTempPhi(const Vector< LevelData<EBCellFAB>* >& a_Velo,
		     const Vector< LevelData<EBCellFAB>* >& a_Temp,
		     const Vector< LevelData<EBCellFAB>* >& a_phi,
		     const Vector< LevelData<EBFluxFAB>* >& a_gradPhi);

  void SetVelTempPhi(const Vector< LevelData<EBCellFAB>* >& a_Velo,
		     const Vector< LevelData<EBCellFAB>* >& a_Temp,
		     const Vector< LevelData<EBCellFAB>* >& a_phi,
		     const Vector< LevelData<EBFluxFAB>* >& a_gradPhi,
		     const Vector< LevelData<EBCellFAB>* >& a_cellGradPhi);
  
  void SetVelTempPhi(const Vector< LevelData<EBCellFAB>* >& a_Velo,
		     const Vector< LevelData<EBCellFAB>* >& a_Temp,
		     const Vector< LevelData<EBCellFAB>* >& a_phi,
		     const Vector< LevelData<EBFluxFAB>* >& a_gradPhi,
		     const Vector< LevelData<EBCellFAB>* >& a_cellGradPhi,
		     const Vector< LevelData<BaseIVFAB<Real>>* >& a_IVGradPhi);

  
  void ZeroVelSetTempPhi(const Vector< LevelData<EBCellFAB>* >& a_Velo,
			 const Vector< LevelData<EBCellFAB>* >& a_Temp,
			 const Vector< LevelData<EBCellFAB>* >& a_phi,
			 const Vector< LevelData<EBFluxFAB>* >& a_gradPhi,
			 const Vector< LevelData<EBCellFAB>* >& a_cellGradPhi,
			 const Vector< LevelData<BaseIVFAB<Real>>* >& a_IVGradPhi);
/*******/
  void define(const Real&                 a_cfl,
	      const RealVect&             a_domainLength,
	      const Real&                 a_refineThresh,
	      const int&                  a_tagBufferSize,
	      const Real&                 a_initialDtMultiplier,
	      const bool&                 a_useLimiting);
/*******/
  void define(AMRLevel*            a_coarserLevelPtr,
	      const ProblemDomain& a_problemDomain,
	      int                  a_level,
	      int                  a_refRatio);
  void defineRefinementBoxes();
  void getIVBoxLimits(IntVect& a_loInt, IntVect& a_hiInt, 
		      RealVect a_loReal,RealVect a_hiReal, 
		      const Real& a_loAdd, const Real& a_hiAdd);
  
/*******/
  void defineSolvers(const bool& flag=true);
/*******/
  void getEBSFOFactory(RefCountedPtr<EBSpeciesFluxOpFactory>&                    a_factory,
		       const Vector<DisjointBoxLayout>&                           a_grids,
		       const Vector<EBISLayout>&                                  a_ebisl,
		       const Vector<EBLevelGrid>&                                 a_eblg,
		       Vector<int>&                                               a_refRat,
		       Real&                                                      a_lev0Dx,
		       Vector<ProblemDomain>&                                     a_domains,
		       const int&                                                 a_iSpec);
/*******/
  void getSpeciesBCFactories(RefCountedPtr<BaseDomainBCFactory>&         a_domBC,
			     RefCountedPtr<BaseEBBCFactory>&             a_ebBC,
			     const Vector<DisjointBoxLayout>&            a_grids,
			     const Vector<EBISLayout>&                   a_ebisl,
			     const int&                                  a_iSpec);
/*******/
  Real advance();
/*********/
  Real diffusiveAdvance(LevelData<EBCellFAB>& a_source);
/*********/
  void makeDiffusiveSource(LevelData<EBCellFAB>& a_diffusiveSrc);
  void makeFluidSource(LevelData<EBCellFAB>& a_fluidSrc);
/*********/
  void  implicitDiffusion(EBFluxRegister*& finerFRPtr, EBFluxRegister*& coarserFRPtr, 
			  LevelData<EBCellFAB>& oldUCoarse, LevelData<EBCellFAB>& stateUCoarse,
			  Real& tCoarserOld, Real& tCoarserNew, LevelData<EBCellFAB>& a_src);
/*********/
  void makeAddReactiveSource(LevelData<EBCellFAB>& a_diffusiveSrc, const Real& a_dt=0);
/*********/
  void makefluxGphi(LevelData<EBCellFAB>& a_fluxGphi);
/*********/
  void getHierarchyAndGrids(Vector<EBAMRSpecies*>&        a_hierarchy,
			    Vector<DisjointBoxLayout>&   a_grids,
			    Vector<EBISLayout>&          a_ebisl,
			    Vector<EBLevelGrid>&         a_eblg,
			    Vector<int>&                 a_refRat,
			    ProblemDomain&               a_lev0Dom,
			    Real&                        a_lev0Dx,
			    Vector<ProblemDomain>&       a_domains,
			    Vector<Real>&                a_dxs);
/*********/
  void getCoarseFineDataPointers(LevelData<EBCellFAB>** a_coarserDataOldPtr, LevelData<EBCellFAB>** a_coarserDataNewPtr,
				 LevelData<EBCellFAB>** a_finerDataOldPtr, LevelData<EBCellFAB>** a_finerDataNewPtr,
				 Real& a_tCoarserOld,  Real& a_tCoarserNew, Real& a_tFinerOld,  Real& a_tFinerNew);
/*********/
  void getCoarseFineDataRegisters(LevelData<EBCellFAB>** a_coarserDataOldPtr, LevelData<EBCellFAB>** a_coarserDataNewPtr,
				  EBFluxRegister** a_coarFR, EBFluxRegister** a_fineFR,
				  Real& a_tCoarserOld,  Real& a_tCoarserNew, Real& a_tFinerOld,  Real& a_tFinerNew);
/*********/
  // gogetters functions
  LevelData<EBCellFAB>* getStatePtr();
  LevelData<EBCellFAB>* getBarPtr();
  Vector< RefCountedPtr<LevelData<BaseIVFAB<Real> > > > getRhosRHS();
/*********/
  void interpolateInTime(LevelData<EBCellFAB>&       a_U,
			 const LevelData<EBCellFAB>& a_UOld,
			 const LevelData<EBCellFAB>& a_UNew,
			 const DisjointBoxLayout&    a_grids,
			 const Real&                 a_time,
			 const Real&                 a_told,
			 const Real&                 a_tnew);
/*********/
  void tagCells(IntVectSet& a_tags);
/*********/
  void regridWithDBL(const DisjointBoxLayout& a_DBL);
/*********/
  void postTimeStep();
/*********/
  void extrapolateZeroArea();
  /*********/
  void doImplicitReflux();
/*******/
  void setSolverCoef(Real a_alpha, Real a_beta);
/*******/
  void initialGrid(const Vector<Box>& a_new_grids);
  void initialGrid(const DisjointBoxLayout& a_DBL);
/*******/
  void levelSetup();
/*******/
  void initialData();
/*******/
void postInitialize();
/*******/
Real computeDt();
/*******/
Real computeInitialDt();
/*******/
EBAMRSpecies* getCoarserLevel() const;
/*******/
EBAMRSpecies* getFinerLevel() const;
/*******/
void syncWithFineLevel();
/*******/
void fillAdvectionVelocity();
/*******/
void getElectronEnergy(LevelData<EBCellFAB>&  a_temp, int a_tComp);
/*******/
void setDiffCoeff(Vector<EBAMRSpecies*>&      a_hierarchy,
		  Vector<DisjointBoxLayout>&   a_grids,
		  Vector<EBISLayout>&          a_ebisl,
		  Vector<EBLevelGrid>&         a_eblg,
		  Vector<int>&                 a_refRat,
		  ProblemDomain&               a_lev0Dom,
		  Vector<Real>&                a_dxs,
		  Vector<ProblemDomain>&       a_domains);
  void junkDiffusion();
/*******/
  LevelData<EBCellFAB>& getStateNew();
/*******/
  void doSmushing(bool a_doSmushing);
/*******/
  void doRZCoords(bool a_doRZCoords);
/*******/
  void hasSourceTerm(bool a_hasSourceTerm);
/*******/
  void useMassRedistribution(bool a_useMassRedist);
/*******/
  void patchSpecies(const EBPatchSpeciesFactory* const a_ebPatchSpeciesFactory);
/*******/
  void refinementThreshold(Real a_refineThresh);
/*******/
  void tagBufferSize(int a_tagBufferSize);
/*******/
  void domainLength(RealVect a_domainLength);
/*******/
  void redistRadius(int a_redistRad);

  void EBcopyLevelData(LevelData<EBCellFAB>&a_copiedTo, const LevelData<EBCellFAB>&a_copyFrom, 
		       const Interval& a_fromInterv,const EBISLayout& a_ebisl);
  void EBcopyLevelData(LevelData<EBCellFAB>&a_copiedTo, const Interval& a_toInterv, const LevelData<EBCellFAB>&a_copyFrom, const Interval& a_fromInterv, const bool doClone=true);
  void EBcopyLevelData(LevelData<EBFluxFAB>&a_copiedTo, const Interval& a_toInterv, const LevelData<EBFluxFAB>&a_copyFrom, const Interval& a_fromInterv, const bool doClone=true);
  void EBcopyLevelData(LevelData<BaseIVFAB<Real> >&a_copiedTo, const Interval& a_toInterv, const LevelData<BaseIVFAB<Real> >&a_copyFrom, const Interval& a_fromInterv, const bool doClone=true);

  void setDt(const Real& a_dt);
  bool PlotReactiveSource(const int& a_nplot=1,const Real& a_Dt=1.0, const bool&  a_onlyfluid = false);

  void multiplyCoeffsByRadius(const Vector<Real>&           a_dxs,
			      const Vector<ProblemDomain>&  a_domains);

  //Unused stuff
  void tagCellsInit(IntVectSet& a_tags);
  void regrid(const Vector<Box>& a_new_grids);
  void writeCheckpointHeader(HDF5Handle& a_handle) const;
  void readCheckpointHeader(HDF5Handle& a_handle);
  void writeCheckpointLevel(HDF5Handle& a_handle) const;
  void readCheckpointLevel(HDF5Handle& a_handle);
  void readCheckpointData(HDF5Handle& a_handle);
  void writePlotHeader(HDF5Handle& a_handle) const;
  void writePlotHeaderOld(HDF5Handle& a_handle) const;
  void writePlotLevel(HDF5Handle& a_handle) const;
  void writePlotLevelOld(HDF5Handle& a_handle) const;
  void outputHDF(const Vector<LevelData<EBCellFAB>* >& a_Ptr, const string& a_name);
  void outputIVFAB(const Vector<LevelData<BaseIVFAB<Real> >* >& a_Ptr, const string& a_name, const int& a_base=0);
  static void setLoadBalance(LoadBalanceFunc a_func)
  {
    s_loadBalance      = a_func;
    s_isLoadBalanceSet = true;
  }
  void cellGrad(LevelData<EBCellFAB>&        a_gradPhi,
		const LevelData<EBCellFAB>&  a_phi,
		const int&                   a_comp);
  void cellGradTens(LevelData<EBCellFAB>&        a_gradPhi,
		    const LevelData<EBCellFAB>&  a_phi,
		    const int&                   a_nComp);
  void SharfetterGummel(LevelData<EBFluxFAB>&        a_fluxMobData,
			LevelData<EBFluxFAB>&        a_fluxDiffData,
			const Real&		       a_dx);

  
/*****************************/
  Real setImplicitTimeStep() const
  {
    return m_PlasmaPhysics->setImplicitTimeStep();
    
  }
  //*******************************************************************************************************//
  static LoadBalanceFunc   s_loadBalance;
  static bool              s_isLoadBalanceSet;
  static bool              s_isOnlyPlasma;
  static int               s_NewPlotFile;

  // Conserved state, U, at old and new time
  LevelData<EBCellFAB> m_stateOld,m_stateNew;
  //Cell Velocity, Temperature  and edge velocity
  Vector<LevelData<EBCellFAB>* > m_Velo;
  Vector<LevelData<EBCellFAB>* > m_Temp;
  Vector<LevelData<EBCellFAB>* > m_phi;
  Vector<LevelData<EBFluxFAB>* > m_gradPhi;
  Vector<LevelData<EBCellFAB>* > m_cellGradPhi;
  Vector<LevelData<BaseIVFAB<Real>>* > m_IVGradPhi;
  LevelData<EBFluxFAB>           m_AdvectiveVelocity;
  LevelData<EBFluxFAB> m_fluxMobData;

  Vector<string> m_stateNames;
  bool m_isDefined;
  Real                 m_cfl;
  Real                 m_refineThresh;
  Real                 m_initialDtMultiplier;
  bool                 m_useLimiting;
  RealVect             m_domainLength;
  RealVect             m_dx;
  bool                 m_useInject;

  // averaging from coarse to fine level
  EBCoarseAverage m_ebCoarseAverage;

  // interpolation from fine to coarse level
  EBPWLFineInterp m_ebFineInterp;

  Real m_dtNew, m_dtExtern;
  int m_nGhost, m_nlevels;
  EBLevelSpecies m_levelSpecies;
  EBFluxRegister m_ebFluxRegister;

  //redistribution mass
  LevelData<BaseIVFAB<Real> > m_massDiff;
  LayoutData<IntVectSet> m_sets;

  EBCoarToFineRedist m_ebCoarToFineRedist;
  EBCoarToCoarRedist m_ebCoarToCoarRedist;
  EBFineToCoarRedist m_ebFineToCoarRedist;

  // tag buffer size
  int m_tagBufferSize;
  int m_redistRad;
  int m_tagAllIrreg;

  bool m_hasCoarser;
  bool m_hasFiner;
  bool m_doImplicitReflux;
  bool m_hasDiffusion;

  // EB stuff added by LM 11/28/13
  bool m_useMassRedist,m_doSmushing,m_doRZCoords,m_doRZCoordsCANC, m_hasSourceTerm, m_isPrintResiduals;
  DisjointBoxLayout m_grids;
  EBISLayout m_ebisl;
  Box m_domainBox;

  //Patch and species pointers
  const EBPatchSpeciesFactory* m_ebPatchSpeciesFactory;
  EBPatchSpecies* m_ebPatchSpecies;
  PlasmaPhysics* m_PlasmaPhysics;
  

  

  // Species Flux Added Luca Massa 11/23/12
  int m_nComp, m_nSpecies, m_plot, m_ghost_Dk, m_steps;

  // for plotting
  RealVect m_origin,m_aspect;
  


  /*------------------------------Class Static Variables     ------------------------*/
  // diffusion coefficients
  static Vector<Vector<RefCountedPtr<LevelData<EBCellFAB> > >  >       s_aco;
  static Vector<Vector<RefCountedPtr<LevelData<EBFluxFAB> > >  >       s_Dk;
  static Vector<Vector<RefCountedPtr<LevelData<EBFluxFAB> > >  >       s_vecW;
  static Vector<Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > > > s_DkIrreg;
  static Vector<Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > > > s_vecWIrreg;
  static Real s_alpha,s_beta;

  // solver to do level  diffusion and flux register interaction
  static Vector<RefCountedPtr<EBLevelCrankNicolson> >                 s_diffuseLevTGA;
  static Vector<RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB> > > > s_diffuseAMRMG;
  static Vector<RefCountedPtr<EBSpeciesFluxOpFactory > >              s_diffuseOpFact;
  static BiCGStabSolver<LevelData<EBCellFAB> >                        s_botSolver;

  // operators
  static Vector<Vector<EBSpeciesFluxOp *>  >                          s_EBAMROps;

  // databased flux boundary conditions
  static Vector< RefCountedPtr<LevelData<BaseIVFAB<Real> > > > s_data;
  // data for the rho RHS
  static Vector< RefCountedPtr<LevelData<BaseIVFAB<Real> > > > s_rhosRHS;
  /*------------------------------End of Class Static Variables     ------------------------*/


  //sfd stuff
  bool m_SFD;
  Real m_Dl,m_chi,m_chiImp,m_chiDt;
  // bar vector at new time
  LevelData<EBCellFAB> m_qbar;

  // for confining the chemistry
  Box m_chemBox;
  

private:

  // Disallowed for all the usual reasons
  void operator=(const EBAMRSpecies&);
  EBAMRSpecies(const EBAMRSpecies&);
};

#include "NamespaceFooter.H"

#endif
