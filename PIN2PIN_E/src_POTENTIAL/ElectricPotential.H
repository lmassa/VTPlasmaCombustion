#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _ELECTRICPOTENTIAL_H_
#define _ELECTRICPOTENTIAL_H_


#include "ParmParse.H"
#include "PolyGeom.H"
#include "GeometryShop.H"
#include "PlaneIF.H"
#include "BRMeshRefine.H"
#include "PlasmaPhysics.H"
#include "AMRMultiGrid.H"
#include "EBAMRPoissonOpRZFactory.H"
#include "BiCGStabSolver.H"
#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "EBCellFactory.H"
#include "EBLevelDataOps.H"
#include "PoissonDomMixBcFunc.H"
#include "PoissonDomMixBcFunc.H"
#include "PoissonEBMixBcFunc.H"
#include "ElectrodesEBDirBcFunc.H"
#include "IonizationWaveEBDirBcFunc.H"

#include "MixedPoissonDomainBC.H"
#include "MixedPoissonEBBC.H"
#include "BaseMixBCValue.H"
#include "CH_HDF5.H"
#include "MultilevelLinearOp.H"
#include "NamespaceHeader.H"

/**
 */
//LM if you have already define Poisson Parameters in PlasmaPhysics, you should not define it here


class PoissonParameters
{
public:
  IntVect       nCells;
  int           maxGridSize;
  int           whichGeom;
  int           blockFactor;
  int           bufferSize;
  Real          fillRatio;
  int           ebBcType;
  int           domBcType;
  int           maxLevel;
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

  Vector<Vector<int> >  domMixBc;
  Real					domBcValue;
};


class ElectricPotential
{
public:
  /// Constructor
  /**
   */
  ElectricPotential();
  /// Define
  void defineGeometry();
  void setupForGivenDBLRun(PlasmaPhysics* a_PlasmaPhysics, const Vector<DisjointBoxLayout>& a_DBL);
  void setupForRestart(PlasmaPhysics* a_PlasmaPhysics, const Vector<DisjointBoxLayout>& a_DBL, HDF5Handle& a_handle);
  void  defineConductivityCoef(Vector<RefCountedPtr<LevelData<EBCellFAB> > >&         a_aco,
			       Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&         a_bco,
			       Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >&  a_bcoIrreg,
			       const Vector<DisjointBoxLayout>&                       a_grids,
			       const Vector<EBISLayout>&                              a_ebisl,
			       const PoissonParameters&                               a_params);
  
  void setConductivityCoefs(LevelData<EBCellFAB >&             a_aco,
			    LevelData<EBFluxFAB>&              a_bco,
			    LevelData<BaseIVFAB<Real> >   &    a_bcoIrreg,
			    const  Real&                       a_dx,
			    const  PoissonParameters&          a_params,
			    const int                     &    a_ilev=0);

  /// Destructor
  /**
   */
  ~ElectricPotential();


  Real solve(const Vector<LevelData<EBCellFAB>* > a_cellNData, const bool& a_evaluate_residual=false);
  void updateRhos(Vector<LevelData<BaseIVFAB<Real> >* > a_rhos);
  void zeroRhos();

  Vector<LevelData<EBCellFAB>* > getPhi();
  Vector<LevelData<EBFluxFAB>* > getGradPhi();
  Vector<LevelData<EBCellFAB>* > getCellGradPhi();
  Vector<LevelData<BaseIVFAB<Real>>* > getIVGradPhi();
  void getEBAMRPFactory(RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > >&   a_levelOpFactory,
			const Vector<DisjointBoxLayout>&        a_grids,
			const Vector< EBISLayout >&             a_ebisl,
			const PoissonParameters&                a_params,
			const int&                              a_numPreCondIters,
			const int&                              a_relaxType,
			const Real&                             a_time,
			const Real&                             a_alpha,
			const Real&                             a_beta);
  void getBCFactories(RefCountedPtr<BaseDomainBCFactory>& a_baseDomainBCFactory,
		      RefCountedPtr<BaseEBBCFactory>&     a_baseEBBCFactory,
		      const PoissonParameters&            a_params);
  
  void getConductivityBCFactories(RefCountedPtr<BaseDomainBCFactory>&                     a_domBC,
				  RefCountedPtr<BaseEBBCFactory>&                         a_ebBC,
				  const Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&    a_bcoef,
				  const Vector<DisjointBoxLayout>&                        a_grids,
				  const Vector<EBISLayout>&                               a_ebisl,
				  const PoissonParameters&                                a_params);
  void getPoissonParameters(PoissonParameters&  a_params, bool a_forceSingleLevel=false);


  void definePoissonGeometry(const PoissonParameters&  a_params);

  void defineSolver();
  void unDefineSolver(){m_isSolverDefined = false;};

  void tagCells(IntVectSet& a_tags, int a_ilev);

  void regridWithDBL(const Vector<DisjointBoxLayout>& a_DBL);

  void setTime(Real a_time);

  void getGradients();

  void cellGrad(LevelData<EBCellFAB>&         a_gradPhi,
		const LevelData<EBCellFAB>&   a_phi,
		const int&                    a_comp,
		const DisjointBoxLayout &     a_grids,
		const EBISLayout &            a_ebisl,
		const ProblemDomain &         a_domain,
		const RealVect &              a_dx);
  

  void outputHDF(const Vector<LevelData<EBCellFAB>* >& a_Ptr, const string& a_name);
  
  void writeCheckpointFile(HDF5Handle& a_handle) const;

  Real getElectrodePotential() const;

  Real getYEB() const;
  bool doTagging() const;
  void zeroField(Vector<LevelData<EBCellFAB>* > a_cellData);
  void zeroSurf(Vector<LevelData<BaseIVFAB<Real> >* > a_cellData);
  void multiplyRHSByRadius();
  void testRHS();
  void extrpolateSolToCentroids();
  void evaluateIVGrad(LevelData<EBCellFAB> &        a_cellData,
		      LevelData<EBCellFAB> &        a_gradData,
		      LevelData<BaseIVFAB<Real>> &  a_IVData,
		      const DisjointBoxLayout &     a_grids,
		      const EBISLayout &            a_ebisl,
		      const ProblemDomain &         a_domain,
		      const RealVect &              a_dx,
		      const int &                   a_ilev);
  void testSOL(Vector< LevelData<EBCellFAB>* > a_resid) const;
  
protected:
  
  int m_nGhost;

  PoissonParameters m_params;
  PlasmaPhysics* m_PlasmaPhysics;
  int m_nComp;


  Vector< LevelData<EBCellFAB>* > m_phi;
  Vector< LevelData<EBFluxFAB>* > m_gradPhi;
  Vector< LevelData<EBCellFAB>* > m_rhs;
  Vector< LevelData<EBCellFAB>* > m_cellGradPhi;
  Vector< LevelData<BaseIVFAB<Real>>* > m_IVGradPhi;

  Real m_alpha, m_beta;
  int m_numPreCondIters,m_relaxtype;
  BiCGStabSolver<LevelData<EBCellFAB> > m_bottomSolver;
  AMRMultiGrid<LevelData<EBCellFAB> > m_solver;
  
  int m_numSmooths;
  Real m_hang, m_eps, m_epsg, m_epsd;
  int m_numMG, m_iterMax;

  
  Vector<DisjointBoxLayout> m_grids;
  Vector<EBISLayout>        m_ebisl;
  Vector<EBLevelGrid>       m_eblg;
  Vector<ProblemDomain>     m_domain;
  Vector<Real>              m_dx;
  Vector<RealVect>          m_dxv;

  //regrid stuff
  int m_tagBufferSize;
  Real m_refineThresh;
  int m_verbosity;

  Real m_time;
  int m_call;

  //luca:: solver stuff
  bool m_isSolverDefined;
  bool m_isDebug, m_isTestResid, m_doRZCoords;
  bool m_isKappaWeight;
  // multi level operator
  MultilevelLinearOp<EBCellFAB >* m_moppa;

  //to get the max volt on the electrode
  Real m_omegaEB, m_voltageEB;
  
  //geometry to be passed outside
  Real m_yEB;
  bool m_isGeometryDefined;
  
};

#include "NamespaceFooter.H"
#endif
