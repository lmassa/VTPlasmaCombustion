#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBAMRNAVIER_H_
#define _EBAMRNAVIER_H_

#include "AMRLevel.H"
#include "EBAMRGodunov.H"
#include "PlasmaPhysics.H"
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
#include "EBLevelGodunov.H"
#include "Box.H"
#include "IntVectSet.H"
#include "Vector.H"
#include "DisjointBoxLayout.H"
#include "EBViscousTensorOpFactory.H"
#include "EBConductivityOpFactory.H"
#include "BiCGStabSolver.H"
#include "EBLevelTGA.H"
#include "EBPatchPolytropic.H"
#include "NamespaceHeader.H"

///interface for load balancing
typedef int(*LoadBalanceFunc)(Vector<int>&         a_procs,
                               const Vector<Box>&   a_boxes,
                               const ProblemDomain& a_domain,
                               bool a_verbose);


class ViscousParameters
{
public:
  int           ebBcType;
  int           domBcType;
  int           ebBcTypeT;
  Real          ebBCValueT;
  int           orderEB;
};

///
/**
 */
class EBAMRNavier : public AMRLevel
{
public:

  EBAMRNavier();

  virtual ~EBAMRNavier();

  void sumConserved(Real& a_sumcons,
                    const int& a_ivar) const;

  virtual void define(AMRLevel*  a_coarser_level_ptr,
              const Box& a_problem_domain,
              int        a_level,
              int        a_ref_ratio);

  virtual void define(AMRLevel*  a_coarser_level_ptr,
              const ProblemDomain& a_problem_domain,
              int        a_level,
              int        a_ref_ratio);

  /// advance by one timestep
  virtual Real advance();
  void residual(LevelData<EBCellFAB>&     a_updated);

  /// things to do after a timestep
  virtual void postTimeStep();

  /// setup fine->coarse stuff
  /**
     stuff that needs to be setup from the finer
     level.  A bunch of objects depend upon the layouts
     from both levels and the finer level changes more
     often from regrid so this needs to be called from the finer
     level
  */
  virtual void syncWithFineLevel();

  void dumpDebug();
  void dumpDebug(const string& a_debstring);

  /// create tags
  virtual void tagCells(IntVectSet& a_tags) ;

  /// create tags at initialization
  virtual void tagCellsInit(IntVectSet& a_tags) ;

  /// regrid
  virtual void regrid(const Vector<Box>& a_new_grids);

  /// initialize grids
  virtual void initialGrid(const Vector<Box>& a_new_grids);

  /// initialize data
  virtual void initialData();

  /// things to do after initialization
  virtual void postInitialize();

  /// object for physics integration
  virtual void patchGodunov(const EBPatchGodunovFactory* const a_ebpatchGodunovFactory);

  //Get injection object
  void setSymmFlag(const int & a_symmetric);
  void setViscPntrs(Vector<RefCountedPtr<LevelData<EBCellFAB> > >*        a_aco,
		    Vector<RefCountedPtr<LevelData<EBFluxFAB> > >*        a_eta,
		    Vector<RefCountedPtr<LevelData<EBFluxFAB> > >*        a_lambda,
		    Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >* a_etaIrreg,
		    Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >* a_lambdaIrreg,
		    Vector<RefCountedPtr<LevelData<EBCellFAB> > >*        a_rcv);
  void defineRefinementBoxes();

#ifdef CH_USE_HDF5
  /// write checkpoint header
  virtual void writeCheckpointHeader(HDF5Handle& a_handle) const;

  /// write checkpoint data for this level
  virtual void writeCheckpointLevel(HDF5Handle& a_handle) const;

  /// read checkpoint header
  virtual void readCheckpointHeader(HDF5Handle& a_handle);

  /// read checkpoint data for this level
  virtual void readCheckpointLevel(HDF5Handle& a_handle);

  /// write plotfile header
  virtual void writePlotHeader(HDF5Handle& a_handle) const;

  void writePlotHeaderOld(HDF5Handle& a_handle) const;

  /// write plot file data for this level
  virtual void writePlotLevel(HDF5Handle& a_handle) const;

  void writePlotLevelOld(HDF5Handle& a_handle) const;  
#endif
  void checkdata(const string& a_msg = string(""));
  ///
  virtual Real getDt() const;

  /// compute dt
  virtual Real computeDt();

  /// compute dt with initial data
  virtual Real computeInitialDt();

  /// CFL
  virtual void CFL(Real a_cfl);

  ///
  virtual void doSmushing(bool a_doSmushing);

  ///
  virtual   void doRZCoords(bool a_doRZCoords);

  ///
  virtual   void hasSourceTerm(bool a_hasSourceTerm);

  ///redistribution radius
  virtual   void redistRadius(int a_redistRad);
  virtual   int getRedistRadius();

  /// longest length of the domain in real space
  virtual   void domainLength(RealVect a_domainLength);

  /// refinement threshold
  virtual   void refinementThreshold(Real a_refine_thresh);

  /// tag buffer size
  virtual   void tagBufferSize(int a_tag_buffer_size);

  virtual   void useMassRedistribution(bool a_useMassRedist);

  virtual void levelSetup();

  const EBISLayout& getEBISLayout() const;
  void getVelTemp(LevelData<EBCellFAB>&  a_vel,LevelData<EBCellFAB>&  a_temp, const bool& a_removeSelfSimilar =false);
  void getVelTemp(Vector<LevelData<EBCellFAB>* >  a_vel, Vector<LevelData<EBCellFAB>* >   a_temp, Vector<DisjointBoxLayout>&   a_grids, Vector<EBPWLFineInterp* > a_ebFineInterps);
  virtual void fillConsAndPrim(LevelData<EBCellFAB>& a_data) const;

  void tagAll(bool a_tagAll);

  // state vector at old time
  LevelData<EBCellFAB>& getStateOld();

  // state vector at new time
  LevelData<EBCellFAB>& getStateNew();

  void setQbar(LevelData<EBCellFAB>*  a_qbar,Real a_chi);
  void setSource(LevelData<EBCellFAB>*  a_src);
  void setSource(Vector<LevelData<EBCellFAB>* >  a_src);

  static int s_NewPlotFile;

  static void setLoadBalance(LoadBalanceFunc a_func)
  {
    s_loadBalance      = a_func;
    s_isLoadBalanceSet = true;
  }

  //lucanote: I've unprotected these two members
  EBAMRNavier* getCoarserLevel() const;
  EBAMRNavier* getFinerLevel() const;

  //vioscus routines
  void defineSolvers();
  void getHierarchyAndGrids(Vector<EBAMRNavier*>&        a_hierarchy,
			    Vector<DisjointBoxLayout>&   a_grids,
			    Vector<EBISLayout>&          a_ebisl,
			    Vector<EBLevelGrid>&         a_eblg,
			    Vector<int>&                 a_refRat,
			    ProblemDomain&               a_lev0Dom,
			    Real&                        a_lev0Dx,
			    Vector<ProblemDomain>&       a_domains,
			    Vector<Real>&                a_dxs);
  
  void  getEBVTOFactory(RefCountedPtr<EBViscousTensorOpFactory>&                    a_factory,
			const Vector<DisjointBoxLayout>&                            a_grids,
			const Vector<EBISLayout>&                                   a_ebisl,
			const Vector<EBLevelGrid>&                                  a_eblg,
			Vector<int>&                                                a_refRat,
			Real&                                                       a_lev0Dx);
  
  void  getViscousBCFactories(RefCountedPtr<BaseDomainBCFactory>&                  a_domBC,
			      RefCountedPtr<BaseEBBCFactory>&                      a_ebbc,
			      const Vector<DisjointBoxLayout>&                     a_grids,
			      const Vector<EBISLayout>&                            a_ebisl);
  void getConductivityFactory(RefCountedPtr<EBConductivityOpFactory>&                    a_factory,
			      const Vector<DisjointBoxLayout>&                            a_grids,
			      const Vector<EBISLayout>&                                   a_ebisl,
			      const Vector<EBLevelGrid>&                                  a_eblg,
			      Vector<int>&                                                a_refRat,
			      Real&                                                       a_lev0Dx,
			      Vector<ProblemDomain>&                                      a_domains);
  void getConductivityBCFactories(RefCountedPtr<BaseDomainBCFactory>&         a_domBC,
				  RefCountedPtr<BaseEBBCFactory>&             a_ebbc,
				  const Vector<DisjointBoxLayout>&            a_grids,
				  const Vector<EBISLayout>&                   a_ebisl);
  
  void makeDiffusiveSource(LevelData<EBCellFAB>& a_diffusiveSrc, const bool& incFlag=false);
  void implicitDiffusion(EBFluxRegister*& finerFRPtr, EBFluxRegister*& coarserFRPtr, 
			 const LevelData<EBCellFAB>& oldUCoarse, const LevelData<EBCellFAB>& stateUCoarse,
			 Real& tCoarserOld, Real& tCoarserNew);
  void getCoarseFineDataPointers(LevelData<EBCellFAB>** a_coarserDataOldPtr, LevelData<EBCellFAB>** a_coarserDataNewPtr,
				 LevelData<EBCellFAB>** a_finerDataOldPtr, LevelData<EBCellFAB>** a_finerDataNewPtr,
				 Real& a_tCoarserOld,  Real& a_tCoarserNew, Real& a_tFinerOld,  Real& a_tFinerNew);

  void consToVelTemp(LevelData<EBCellFAB>& a_V, LevelData<EBCellFAB>& a_T, const LevelData<EBCellFAB>& a_consState);
  void consToVelTemp(EBCellFAB& a_V, EBCellFAB& a_T, const EBCellFAB& a_consState, const Box& a_box);
  void updateConserved(LevelData<EBCellFAB>& a_consState, const LevelData<EBCellFAB>& a_V, const LevelData<EBCellFAB>& a_T);
  void consToDiffusiveProps(LevelData<EBCellFAB>& a_R, LevelData<EBCellFAB>& a_E, LevelData<EBCellFAB>& a_RCV, const LevelData<EBCellFAB>& a_consState, const ProblemDomain& a_domain);
  void interpolateInTime(LevelData<EBCellFAB>& a_U,const LevelData<EBCellFAB>& a_UOld,const LevelData<EBCellFAB>& a_UNew,
			 const DisjointBoxLayout& a_grids, const Real& a_time, const Real& a_told, const Real& a_tnew);
  void setViscCoeff(Vector<RefCountedPtr<LevelData<EBCellFAB> > >&        a_aco,
		    Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&        a_eta,
		    Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&        a_lambda,
		    Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >& a_etaIrreg,
		    Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >& a_lambdaIrreg,
		    Vector<RefCountedPtr<LevelData<EBCellFAB> > >&        a_rcv,
		    Vector<EBAMRNavier*>&        a_hierarchy,
		    Vector<DisjointBoxLayout>&   a_grids,
		    Vector<EBISLayout>&          a_ebisl,
		    Vector<EBLevelGrid>&         a_eblg,
		    Vector<int>&                 a_refRat,
		    ProblemDomain&               a_lev0Dom,
		    Vector<Real>&                a_dxs,
		    Vector<ProblemDomain>&       a_domains);
  void addDKEDt(LevelData<EBCellFAB>& a_srcT,const LevelData<EBCellFAB>& a_updateV, const LevelData<EBCellFAB>& a_initialV, const LevelData<EBCellFAB>& a_consState);
  void computeVorticity(LevelData<EBCellFAB>& a_vort, const LevelData<EBCellFAB>& a_velo);
  void getVorticity(LevelData<EBCellFAB>& a_vort);
  void PlotVorticity(const int& a_nplot);
  void PlotVelocity(const int& a_nplot, const Real& a_Vel0 = 1, const bool& a_removeSelfSimilar=false);
  void copyLevel(LevelData<EBCellFAB>& a_dest, const LevelData<EBCellFAB>& a_src, const Interval& destComps, const Interval& srcComps);

  
  // for confining the geom
  void getIVBoxLimits(IntVect& a_loInt, IntVect& a_hiInt, RealVect a_loReal,RealVect a_hiReal, const Real& a_loAdd =0, const Real& a_hiAdd=0);
  Box m_vortBox;
protected:

  static LoadBalanceFunc   s_loadBalance;
  static bool              s_isLoadBalanceSet;
  bool m_tagAll;
  bool m_useMassRedist;
  Box m_domainBox;
  Real calculateMass() const;


  DisjointBoxLayout m_grids;
  EBISLayout m_ebisl;

  // state vector at old time
  LevelData<EBCellFAB> m_stateOld;

  // state vector at new time
  LevelData<EBCellFAB> m_stateNew;


  // number of components of m_state
  int m_nComp;

  int m_nGhost ;
  // names of components
  Vector<string> m_stateNames;
  Vector<string> m_primNames;

  // grid spacing
  RealVect m_origin;
  RealVect m_dx;
  RealVect m_aspect;

  // domain length
  RealVect m_domainLength;

  // cfl number
  Real m_cfl;

  // new time step
  Real m_dtNew;

  // averaging from coarse to fine level
  EBCoarseAverage m_ebCoarseAverage;

  // interpolation from fine to coarse level
  EBPWLFineInterp m_ebFineInterp;

  // level solver
  EBLevelGodunov m_ebLevelGodunov;

  // flux register
  EBFluxRegister m_ebFluxRegister;

  // patch integrator
  EBPatchGodunov* m_ebPatchGodunov;
  EBPatchPolytropic* m_ebPolytropic;

  // patch integrator
  const EBPatchGodunovFactory* m_ebPatchGodunovFactory;

  // refinement threshold for gradient
  Real m_refineThresh;

  //redistribution mass
  LevelData<BaseIVFAB<Real> > m_massDiff;
  LayoutData<IntVectSet> m_sets;

  RedistStencil m_redStencil;
  EBCoarToFineRedist m_ebCoarToFineRedist;
  EBCoarToCoarRedist m_ebCoarToCoarRedist;
  EBFineToCoarRedist m_ebFineToCoarRedist;

  // tag buffer size
  int m_tagBufferSize;
  int m_redistRad;

  bool m_doRZCoords;
  bool m_hasSourceTerm;
  bool m_doSmushing;
  bool m_hasCoarser;
  bool m_hasFiner, m_hasFinerLevel;
  bool m_isDefined;


  //output streams
  std::ofstream     m_out;
  Real m_NAN;
  int m_dumps;
  
  //symmetric flag
  int m_symmetric;

  //Viscous instance variables
  Real m_alpha, m_gamma, m_cv, m_lambdafac,m_Tinf;
  bool m_hasDiffusion, m_explicitDiffusion, m_useInject, m_isFunctionEBVTO;  
  ViscousParameters m_params;
  Vector<RefCountedPtr<LevelData<EBCellFAB> > >*         m_aco;
  Vector<RefCountedPtr<LevelData<EBFluxFAB> > >*         m_eta;
  Vector<RefCountedPtr<LevelData<EBFluxFAB> > >*         m_lambda;
  Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >*  m_etaIrreg;
  Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >*  m_lambdaIrreg;
  Vector<RefCountedPtr<LevelData<EBCellFAB> > >*         m_rcv;

  //Viscous class variables
  static RefCountedPtr<EBLevelBackwardEuler>                 s_viscLevTGA;
  static RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB> > > s_viscAMRMG;
  static RefCountedPtr<EBViscousTensorOpFactory >            s_viscOpFact;

  static RefCountedPtr<EBLevelBackwardEuler>                 s_condLevTGA;
  static RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB> > > s_condAMRMG;
  static RefCountedPtr<EBConductivityOpFactory  >            s_condOpFact;

  static BiCGStabSolver<LevelData< EBCellFAB> >              s_botSolver;
  static Vector<EBViscousTensorOp *> s_EBAMROps;
  static Vector<EBConductivityOp  *> s_EBAMROpsT;
  static Real s_beta, s_betaCond;

  LevelData<EBCellFAB> m_initialV;
  LevelData<EBCellFAB> m_initialT;
  LevelData<EBCellFAB> m_srcV;
  LevelData<EBCellFAB> m_srcT;
  LevelData<EBCellFAB> m_srcF;
  LevelData<EBCellFAB>* m_srcE;
  LevelData<EBCellFAB> m_updateT;
  LevelData<EBCellFAB> m_updateV;
  LevelData<EBCellFAB> m_newVCoarse;
  LevelData<EBCellFAB> m_newTCoarse;
  LevelData<EBCellFAB> m_oldVCoarse;
  LevelData<EBCellFAB> m_oldTCoarse;



  //sfd stuff
  bool m_SFD;
  Real m_Dl,m_chi;
  // bar vector at new time
  LevelData<EBCellFAB> m_qbar;
  
  int m_plot;
  void PlotReactiveSource();
  void outputHDF(const Vector<LevelData<EBCellFAB>* >& a_Ptr, const string& a_name);

private:
  //disallowed for all the usual reasons
  void operator=(const EBAMRNavier& a_input)
  {
    MayDay::Error("invalid operator");
  }
  EBAMRNavier(const EBAMRNavier& a_input)
  {
    MayDay::Error("invalid operator");
  }

};

#include "NamespaceFooter.H"
#endif
