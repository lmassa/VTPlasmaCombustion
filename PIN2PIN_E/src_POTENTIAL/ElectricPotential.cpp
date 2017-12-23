#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <algorithm>    // std::min
#include "ElectricPotential.H"
#include "ElectricPotentialF_F.H"
#include "EBAMRSpeciesF_F.H"
#include "SlabService.H"
#include "PlaneService.H"
#include "PlaneIF.H"
#include "IntersectionIF.H"
#include "UnionIF.H"
#include "MultiSphereIF.H"
#include "GeometryShop.H"
#include "NeumannPoissonDomainBC.H"
#include "DirichletPoissonDomainBC.H"
#include "MixedPoissonDomainBC.H"
#include "MixedConductivityDomainBC.H"
#include "BaseEBBC.H"
#include "DirichletPoissonEBBC.H"
#include "NeumannPoissonEBBC.H"
#include "EBPWLFillPatch.H"
#include "EBPWLFineInterp.H"
#include "LayoutData.H"
#include "BoxLayoutData.H"
#include "EBCellFAB.H"
#include "EBAMRDataOps.H"
#include "EBLevelCCProjector.H" // for averaging
#include "DirichletConductivityDomainBC.H"
#include   "NeumannConductivityDomainBC.H"
#include   "MixedConductivityDomainBC.H"
#include "DirichletConductivityEBBC.H"
#include   "NeumannConductivityEBBC.H"
//debug//#include "MultilevelLinearOp.H"
#include "NamespaceHeader.H"

/***/
//constructor
ElectricPotential::ElectricPotential()
{

  m_verbosity = 1;
  m_nComp = 1;
  m_time = 0;
  m_call = 0; 
  m_isSolverDefined = false;
  m_isGeometryDefined = false;
  m_isDebug = false;
  m_isTestResid = false;
  // LM I think we shall kappaweight it here.
  m_isKappaWeight = true;
  m_doRZCoords = true;

  //input deck
  ParmParse pp;
  m_relaxtype=1;
  pp.get("mg_relax_type", m_relaxtype);
  /*
  if (m_relaxtype !=0)
    {
      bool doLazyRelax;
      pp.get("do_lazy_relax", doLazyRelax);
      EBAMRPoissonOpRZ::doLazyRelax(doLazyRelax);
    }
  Real slowFactor=-1;
  pp.query("slow_factor", slowFactor);
  if(slowFactor > 0)EBAMRPoissonOpRZ::setSlowSafetyFactor(slowFactor);
  */

  m_numPreCondIters = 40;
  pp.get("num_pre_cond_iters",m_numPreCondIters);

  pp.get("eps_gas", m_epsg);
  pp.get("eps_diel", m_epsd);

  pp.get("mg_num_smooths", m_numSmooths);
  pp.get("mg_hang", m_hang);
  pp.get("mg_eps",  m_eps);
  pp.get("mg_num_cycles", m_numMG);
  pp.get("mg_iter_max",   m_iterMax);

  pp.get("alpha", m_alpha);
  pp.get("beta",  m_beta);

  //refinement
  pp.get("electric_refine_thresh", m_refineThresh);
  //query three times to get one of the three, preferably the last
  m_tagBufferSize=1;
  pp.query("vorticity_tag_buffer_size", m_tagBufferSize);
  pp.query("species_tag_buffer_size", m_tagBufferSize);
  pp.query("electric_tag_buffer_size", m_tagBufferSize);

  //added to print the elec potential value
  Real frequency=0;
  m_omegaEB = 0;m_voltageEB=0;
  pp.query("eb_dir_bc_value", m_voltageEB);
  pp.query("frequency", frequency);
  m_omegaEB = 2*M_PI*frequency;

  //verbosity
  pp.query("verbosity",m_verbosity);
  
  if(m_verbosity >=5) m_moppa=NULL;
  
}
void ElectricPotential::defineGeometry()
{
  getPoissonParameters(m_params);

  definePoissonGeometry(m_params);

  //allocate data holders
  m_phi.resize(m_params.numLevels);
  m_rhs.resize(m_params.numLevels);
  m_gradPhi.resize(m_params.numLevels);
  m_cellGradPhi.resize(m_params.numLevels);
  m_IVGradPhi.resize(m_params.numLevels,NULL);
  for (int ilev = 0; ilev <= m_params.maxLevel; ilev++)
    {
      m_phi[ilev]   = new LevelData<EBCellFAB>();
      m_rhs[ilev]   = new LevelData<EBCellFAB>();
      m_gradPhi[ilev]   = new LevelData<EBFluxFAB>();
      m_cellGradPhi[ilev]   = new LevelData<EBCellFAB>();
      m_IVGradPhi[ilev]   = new LevelData<BaseIVFAB<Real> >();
    }

  m_domain.resize(m_params.numLevels);
  m_dx.resize(m_params.numLevels);
  m_dxv.resize(m_params.numLevels);
  m_domain[0] = m_params.coarsestDomain;
  m_dx[0] = m_params.coarsestDx[0];
  for (int ilev = 1; ilev < m_params.numLevels; ++ilev)
    {
      m_domain[ilev] = refine(m_domain[ilev-1],m_params.refRatio[ilev-1]);
      m_dx[ilev] = m_dx[ilev-1]/m_params.refRatio[ilev-1];
    }
  for (int ilev = 0; ilev < m_params.numLevels; ++ilev)m_dxv[ilev]= m_dx[ilev]* RealVect::Unit;

  m_isGeometryDefined=true;
}

/********************************************************/
void ElectricPotential::defineSolver()
{


  ParmParse pp;

  const EBIndexSpace* const ebisPtr = Electric_EBIS::instance();

  // this part was taken from PoissonUtilities::defineSolver
  RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > > operatorFactory;
  int time = 0e0;
  getEBAMRPFactory(operatorFactory, m_grids,  m_ebisl, m_params, m_numPreCondIters, m_relaxtype,
                   time, m_alpha,  m_beta);
  ProblemDomain coarsestDomain(m_params.coarsestDomain);


  m_bottomSolver.m_verbosity = 0;
  m_bottomSolver.m_numRestarts=20;
  m_bottomSolver.m_imax=100;
  m_bottomSolver.m_recounts=100;
  m_bottomSolver.m_eps=1e-8; //this is redefined in the mutligrid solver

  int numLevels=m_grids.size();
  int maxDepth =1;
  pp.query("max_depth", maxDepth);
  m_solver.setMaxDepth(maxDepth);
  m_solver.define(coarsestDomain, *operatorFactory, &m_bottomSolver, numLevels);
  Real normThresh = 1.0e-30;
  m_solver.setSolverParameters(m_numSmooths, m_numSmooths, m_numSmooths,
                               m_numMG, m_iterMax, m_eps, m_hang, normThresh);
  m_solver.m_verbosity = m_verbosity;
  int imin =  max(floor(m_iterMax/2.0),5.0);
  pp.query("imin", imin);
  m_solver.m_imin =imin;
  //for non-perfect convergence set the bottom solver cushion
  if(m_eps > 1e-8) m_solver.setBottomSolverEpsCushion(m_eps);

  
  if(m_isTestResid){
    Vector<ProblemDomain> domains(numLevels);
    for(int k=0;k < numLevels;k++)domains[k]=m_domain[k];
    if(m_moppa != NULL)delete m_moppa;
    m_moppa = new MultilevelLinearOp< EBCellFAB >();
    m_moppa->define(m_grids, m_params.refRatio, domains, m_dxv, operatorFactory, 0);
  }
  m_isSolverDefined = true;
}
/******************************************/
void ElectricPotential::setupForGivenDBLRun(PlasmaPhysics* a_PlasmaPhysics, const Vector<DisjointBoxLayout>& a_DBL)
{
  m_PlasmaPhysics = a_PlasmaPhysics;

  const EBIndexSpace* const ebisPtr = Electric_EBIS::instance();
  m_grids.resize(m_params.numLevels);
  m_ebisl.resize(m_params.numLevels);
  m_eblg.resize(m_params.numLevels);
    
  // fill the grid and Layout
  m_nGhost = 4;
  int maxloop = std::min(m_params.numLevels,(int)a_DBL.size());
  pout () << "ElectricPotential::setupForGivenDBLRun- Max number of Levels are " << maxloop << endl;
  for (int ilev = 0; ilev < maxloop; ++ilev)
    {
      m_grids[ilev] = a_DBL[ilev]; 
      ebisPtr->fillEBISLayout(m_ebisl[ilev], m_grids[ilev], m_domain[ilev], m_nGhost);
      m_eblg[ilev] = EBLevelGrid(m_grids[ilev], m_domain[ilev], m_nGhost, ebisPtr);
      EBCellFactory factory(m_ebisl[ilev]);
      m_phi[ilev]->define(m_grids[ilev], m_nComp, m_params.ghostPhi, factory);
      m_rhs[ilev]->define(m_grids[ilev], m_nComp, m_params.ghostRHS, factory);
      EBLevelDataOps::setToZero(*m_phi[ilev]);
      EBLevelDataOps::setToZero(*m_rhs[ilev]);
    }
}
void ElectricPotential::setupForRestart(PlasmaPhysics* a_PlasmaPhysics, const Vector<DisjointBoxLayout>& a_DBL, HDF5Handle& a_handle)
{

  m_PlasmaPhysics = a_PlasmaPhysics;
  int numLevels =  a_DBL.size();

  const EBIndexSpace* const ebisPtr = Electric_EBIS::instance();
  m_grids.resize(numLevels);
  m_ebisl.resize(numLevels);
  m_eblg.resize(numLevels);
  m_phi.resize(numLevels);
  m_rhs.resize(numLevels);
  m_gradPhi.resize(numLevels);
  m_cellGradPhi.resize(numLevels);
  m_IVGradPhi.resize(numLevels);
    
  // fill the grid and Layout
  m_nGhost = 4;
  for (int ilev = 0; ilev < numLevels; ++ilev)
    {
      m_grids[ilev] = a_DBL[ilev]; 
      ebisPtr->fillEBISLayout(m_ebisl[ilev], m_grids[ilev], m_domain[ilev], m_nGhost);
      m_eblg[ilev] = EBLevelGrid(m_grids[ilev], m_domain[ilev], m_nGhost, ebisPtr);
      //m_eblg[ilev] = EBLevelGrid(m_grids[ilev], m_ebisl[ilev], m_domain[ilev]);
      EBCellFactory factory(m_ebisl[ilev]);
      m_phi[ilev]->define(m_grids[ilev], m_nComp, m_params.ghostPhi, factory);
      m_rhs[ilev]->define(m_grids[ilev], m_nComp, m_params.ghostRHS, factory);
      EBLevelDataOps::setToZero(*m_phi[ilev]);
      EBLevelDataOps::setToZero(*m_rhs[ilev]);


      // Setup the level string
      char levelStr[20];
      sprintf(levelStr,"%d",ilev);
      const std::string label = std::string("SPlevel_") + levelStr;
      a_handle.setGroup(label);
      int dataStatusNew = read<EBCellFAB>(a_handle,*m_phi[ilev], "phi", m_grids[ilev], Interval(), false); 
      if ((dataStatusNew != 0) )
	{
	  MayDay::Error("file does not contain state data");
	}
      m_phi[ilev]->exchange();
    }
  
  defineSolver();

  getGradients();
   
}
/***/
ElectricPotential::
~ElectricPotential()
{
}
/********************************************************/
Real ElectricPotential::solve(const Vector<LevelData<EBCellFAB>* > a_cellNData, const bool& a_evaluate_residual)
{

  if(!m_isSolverDefined) defineSolver();

  //setup the RHS
  int numLevels =  a_cellNData.size();
  EBAMRDataOps::setToZero(m_rhs);
  if(!m_isDebug){
    for (int ilev = 0; ilev < numLevels; ++ilev)
      {   
	if(m_verbosity > 3)pout() << ilev << " Evaluating RHS in ElectricPotential " << endl;   
	m_PlasmaPhysics->ElectricPotentialRHS(*m_rhs[ilev],*a_cellNData[ilev], m_dx[ilev]);
	if(m_isKappaWeight) EBLevelDataOps::kappaWeight(*m_rhs[ilev]);
	if(m_verbosity > 3)pout() << ilev << " Evaluated RHS in ElectricPotential " << endl; 
      }
  }
  else
    {
      testRHS();
    }
  
  if(m_doRZCoords) multiplyRHSByRadius();
  
  m_solver.solve(m_phi, m_rhs, m_params.maxLevel, 0, false);
  Real retVal = m_solver.m_residNorm;
  //LM Check:: this is necessary for verification
  // Should I do it also when doing standard computations????
  if (m_isDebug) extrpolateSolToCentroids();
  std::ofstream ofs;
  if (m_isDebug)
    {
      ofs.open ("g.dat", std::ofstream::out | std::ofstream::trunc);
    }
  
  if(m_verbosity >=10)
    {
      pout () << " Residuals " << retVal << endl;
      EBAMRDataOps::setCoveredVal(m_phi,10);
      outputHDF(m_phi,string("phi"));
      Vector< LevelData<EBCellFAB>* > resid(numLevels,NULL);
      for (int ilev = 0; ilev <= m_params.maxLevel; ilev++)
	{
	  EBCellFactory factory(m_ebisl[ilev]);
	  resid[ilev] = new LevelData<EBCellFAB>(m_grids[ilev], m_nComp, m_params.ghostRHS, factory);
	}
      if(m_isDebug)testSOL(resid);
      else if (m_isTestResid)
	m_moppa->residual(resid,m_phi,m_rhs,false);// homog last
      outputHDF(resid,string("residC"));
      for (int ilev = 0; ilev <= m_params.maxLevel; ilev++)delete resid[ilev];
    }
  
  

  getGradients();
  
  if (m_isDebug)
    {
      ofs.close();
    }

  return retVal;
}

/********************************************************/
void ElectricPotential::updateRhos(Vector<LevelData<BaseIVFAB<Real> >* > a_rhos)
{ 
  zeroSurf(a_rhos);
  // here I tested it makes no difference, the MG operators are likely to be homog
  Vector< AMRLevelOp<LevelData<EBCellFAB> >* > ops = m_solver.getAMROperators();
  //Vector< MGLevelOp<LevelData<EBCellFAB> >* > ops = m_solver.getAllOperators();
  for (int iop = 0; iop < ops.size(); iop++)  
    {
      AMRLevelOp<LevelData<EBCellFAB> >* opA = &(*ops[iop]);
      //MGLevelOp<LevelData<EBCellFAB> >* opA = &(*ops[iop]);
      EBAMRPoissonOpRZ* op = dynamic_cast<EBAMRPoissonOpRZ* >(opA);
      //EBAMRPoissonOpRZ& op = *(static_cast <EBAMRPoissonOpRZ> (opA));
      op->setEBBCSource(a_rhos, m_domain);
      if(m_verbosity > 3)pout() << iop << " Updated rhos in ElectricPotential " << endl;
    }

  if(m_isTestResid){
    Vector<RefCountedPtr<AMRLevelOp<LevelData<EBCellFAB> > > > opx = m_moppa->m_vectOperators;
    for (int iop = 0; iop < opx.size(); iop++)  
      {
	AMRLevelOp<LevelData<EBCellFAB> >* opA = &(*opx[iop]);
	EBAMRPoissonOpRZ* op = dynamic_cast<EBAMRPoissonOpRZ* >(opA);
	//EBAMRPoissonOpRZ& op = *(static_cast <EBAMRPoissonOpRZ> (opA));
	op->setEBBCSource(a_rhos, m_domain);
	if(m_verbosity > 3)pout() << iop << " Updated rhos in ElectricPotential " << endl;
      }
  }

}
/********************************************************/
void ElectricPotential::zeroRhos()
{ 
  // here I tested it makes no difference, the MG operators are likely to be homog
  Vector< AMRLevelOp<LevelData<EBCellFAB> >* > ops = m_solver.getAMROperators();
  //Vector< MGLevelOp<LevelData<EBCellFAB> >* > ops = m_solver.getAllOperators();
  for (int iop = 0; iop < ops.size(); iop++)  
    {
      AMRLevelOp<LevelData<EBCellFAB> >* opA = &(*ops[iop]);
      //MGLevelOp<LevelData<EBCellFAB> >* opA = &(*ops[iop]);
      EBAMRPoissonOpRZ* op = dynamic_cast<EBAMRPoissonOpRZ* >(opA);
      //EBAMRPoissonOpRZ& op = *(static_cast <EBAMRPoissonOpRZ> (opA));
      op->zeroEBBCSource();
    }


}

Vector<LevelData<EBCellFAB>* >
ElectricPotential::getPhi()
{
  return m_phi;
}

Vector<LevelData<EBFluxFAB>* >
ElectricPotential::getGradPhi()
{
  return m_gradPhi;
}
Vector<LevelData<EBCellFAB>* >
ElectricPotential::getCellGradPhi()
{
  return m_cellGradPhi;
}
Vector<LevelData<BaseIVFAB<Real>>* >
ElectricPotential::getIVGradPhi()
{
  return m_IVGradPhi;
}

/****************/
//LM untouched
void ElectricPotential:: 
getEBAMRPFactory(RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > >&   a_levelOpFactory,
                 const Vector<DisjointBoxLayout>&        a_grids,
                 const Vector< EBISLayout >&             a_ebisl,
                 const PoissonParameters&                a_params,
                 const int&                              a_numPreCondIters,
                 const int&                              a_relaxType,
                 const Real&                             a_time,
                 const Real&                             a_alpha,
                 const Real&                             a_beta)
{
  CH_TIME("PoissonUtilities::getEBAMRPFactory");  
  const EBIndexSpace* const ebisPtr = Electric_EBIS::instance();
  RefCountedPtr<BaseDomainBCFactory> baseDomainBCFactory;
  RefCountedPtr<BaseEBBCFactory>         baseEBBCFactory;


  RealVect vectDx = RealVect::Unit;
  vectDx *= a_params.coarsestDx;
  int numLevels = -1;
  Vector<RefCountedPtr<EBQuadCFInterp> > quadCFI(a_grids.size());
  ProblemDomain levDom =  a_params.coarsestDomain;
  ProblemDomain coarDom;
  for (int ilev = 0; ilev < a_grids.size(); ilev++)
    {
      if (ilev > 0)
        {
          int nref = a_params.refRatio[ilev-1];
	  pout () << "Doing Level ....... " << ilev << endl;
          quadCFI[ilev] = RefCountedPtr<EBQuadCFInterp>(new EBQuadCFInterp(a_grids[ilev],
                                                                           a_grids[ilev-1],
                                                                           a_ebisl[ilev],
                                                                           a_ebisl[ilev-1],
                                                                           coarDom, nref, m_nComp,
                                                                           *m_eblg[ilev].getCFIVS(),
									   ebisPtr));
          coarDom.refine(a_params.refRatio[ilev]);
        }
      coarDom = levDom;
      levDom.refine(a_params.refRatio[ilev]);
    }
  /*a_levelOpFactory = RefCountedPtr<EBAMRPoissonOpRZFactory>(new EBAMRPoissonOpRZFactory(m_eblg, a_params.refRatio, quadCFI, //
                                                                                    vectDx, RealVect::Zero,
                                                                                    a_numPreCondIters,a_relaxType,
                                                                                    baseDomainBCFactory, baseEBBCFactory,
                                                                                    a_alpha, a_beta, a_time,a_params.ghostPhi, a_params.ghostRHS,
                                                                                    numLevels));*/
  Vector<RefCountedPtr<LevelData<EBCellFAB> > >           aco;
  Vector<RefCountedPtr<LevelData<EBFluxFAB> > >           bco;
  Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >    bcoIrreg;
  defineConductivityCoef(aco, bco, bcoIrreg, a_grids, a_ebisl, a_params);
  
  getConductivityBCFactories(baseDomainBCFactory, baseEBBCFactory, bco, a_grids, a_ebisl, a_params);
  a_levelOpFactory = RefCountedPtr<EBAMRPoissonOpRZFactory>(new EBAMRPoissonOpRZFactory(m_eblg, quadCFI, a_alpha, a_beta,
											aco, bco, bcoIrreg,//
											a_params.coarsestDx[0], a_params.refRatio,
											baseDomainBCFactory, baseEBBCFactory,
											a_params.ghostPhi, a_params.ghostRHS, a_relaxType));
}

void ElectricPotential:: 
defineConductivityCoef(Vector<RefCountedPtr<LevelData<EBCellFAB> > >&           a_aco,
                       Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&           a_bco,
                       Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >&    a_bcoIrreg,
                       const Vector<DisjointBoxLayout>&                         a_grids,
                       const Vector<EBISLayout>&                                a_ebisl,
                       const PoissonParameters&                                 a_params)
{
  CH_TIME("PoissonUtilities::defineConductivityCoef");
  a_aco.resize(        a_params.numLevels);
  a_bco.resize(        a_params.numLevels);
  a_bcoIrreg.resize(   a_params.numLevels);
  Real dxLev = a_params.coarsestDx[0];
  ProblemDomain domLev = a_params.coarsestDomain;
  for (int ilev = 0; ilev < a_params.numLevels; ilev++)
    {
      LayoutData<IntVectSet> irregSets(a_grids[ilev]);
      int nghost = 1;
      for (DataIterator dit = a_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
          //has to correspond to number of ghost cells
          Box grownBox = grow(a_grids[ilev].get(dit()), nghost);
          grownBox &= domLev;
          irregSets[dit()] = a_ebisl[ilev][dit()].getIrregIVS(grownBox);
        }
      EBFluxFactory        ebfluxfact(a_ebisl[ilev]);
      EBCellFactory        ebcellfact(a_ebisl[ilev]);
      BaseIVFactory<Real>  baseivfact(a_ebisl[ilev], irregSets);

      a_aco[ilev]         = RefCountedPtr<LevelData<EBCellFAB       > >(new LevelData<EBCellFAB       >(a_grids[ilev], 1, nghost*IntVect::Unit, ebcellfact));
      a_bco[ilev]         = RefCountedPtr<LevelData<EBFluxFAB       > >(new LevelData<EBFluxFAB       >(a_grids[ilev], 1, nghost*IntVect::Unit, ebfluxfact));
      a_bcoIrreg[ilev]    = RefCountedPtr<LevelData<BaseIVFAB<Real> > >(new LevelData<BaseIVFAB<Real> >(a_grids[ilev], 1, nghost*IntVect::Unit, baseivfact));
      setConductivityCoefs(*a_aco[ilev], *a_bco[ilev], *a_bcoIrreg[ilev], dxLev, a_params, ilev);
      dxLev /=      a_params.refRatio[ilev];
      domLev.refine(a_params.refRatio[ilev]);
  }
}


void ElectricPotential::
setConductivityCoefs(LevelData<EBCellFAB>          &    a_aco,
		     LevelData<EBFluxFAB>          &    a_bco,
		     LevelData<BaseIVFAB<Real> >   &    a_bcoIrreg,
		     const Real                    &    a_dx,
		     const  PoissonParameters      &    a_params,
		     const int                     &    a_ilev)
{
  int icomp = 0; 
  Real acoVal=0, bcoVal=1;
  RealVect vectDx = a_dx*(RealVect::Unit);
  for (DataIterator dit = a_aco.dataIterator(); dit.ok(); ++dit)
    {
      a_aco     [dit()].setVal(acoVal);
      a_bco     [dit()].setVal(bcoVal);
      a_bcoIrreg[dit()].setVal(bcoVal);
    }
  if (true)
    {
      RealVect cylinderAxis = BASISREALV(0);

      for (DataIterator dit = a_bco.dataIterator(); dit.ok(); ++dit)
        {
	  for (int idir = 0; idir < SpaceDim; idir++)
	    {
	      FaceStop::WhichFaces stopCritGrid =  FaceStop::SurroundingWithBoundary;
	      EBFaceFAB& dataEBFAB = a_bco[dit()][idir];
	      const Box& region = dataEBFAB.getCellRegion();
	      Box boxS = a_bco.box(dit());
	      IntVectSet ivsBox(region);
	      const EBISBox& ebisBox = dataEBFAB.getEBISBox();
	      for (FaceIterator faceit(ivsBox, ebisBox.getEBGraph(), idir, stopCritGrid); faceit.ok(); ++faceit)
		{
		  const FaceIndex& face = faceit();
		  RealVect centroid = ebisBox.centroid(face);
		  centroid *= a_dx;
		  RealVect faceloc = EBArith::getFaceLocation(face, vectDx, centroid);
		  Real dist = getDistanceFromAxis(faceloc, cylinderAxis, RealVect::Zero);
		  dataEBFAB(face,icomp) = m_doRZCoords? dist : 1;
		}
	    }
	  
	  EBCellFAB& dataEBFAB = a_aco[dit()];
	  const EBISBox& ebisBox = dataEBFAB.getEBISBox();
	  const Box& region = dataEBFAB.getRegion();
	  Box grownBox = grow(region, 0);
	  grownBox &= m_domain[a_ilev];
	  IntVectSet ivsIrreg = ebisBox.getIrregIVS(grownBox);
	  for (VoFIterator vofit(ivsIrreg, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
	    {
	      RealVect centroid = ebisBox.bndryCentroid(vofit());
	      centroid *= a_dx;
	      RealVect vofloc = EBArith::getVofLocation(vofit(), vectDx, centroid);
	      Real dist = getDistanceFromAxis(vofloc, cylinderAxis, RealVect::Zero);
	      a_bcoIrreg[dit()](vofit(), icomp) = m_doRZCoords? dist : 1;
	    }

	}
    }
  
  a_aco.exchange(Interval(0,0));
  a_bco.exchange(Interval(0,0));
  a_bcoIrreg.exchange(Interval(0,0));

}

void ElectricPotential::getBCFactories(RefCountedPtr<BaseDomainBCFactory>& a_baseDomainBCFactory,
				       RefCountedPtr<BaseEBBCFactory>&     a_baseEBBCFactory,
				       const PoissonParameters&            a_params)
{
  CH_TIME("PoissonUtilities::getBCFactories");
  ParmParse pp;
  //[MD] for mixed Dirichlet/Neumann BCs
  if (a_params.domBcType == 100) 
    {
      pout() << "mixed Dirichlet/Neumann bcs on domain" << endl;
      Real domDirBcValue, domNeuBcValue, EbDirBcValue, omega, frequency;
      pp.get("dom_dir_bc_value", domDirBcValue);
      pp.get("dom_neu_bc_value", domNeuBcValue);
      pp.get("Potential", EbDirBcValue);
      pp.get("frequency", frequency);
      omega = 2*M_PI*frequency;

      int max_level = a_params.maxLevel;
      RealVect fineDx = a_params.coarsestDx;
      for (int ilev = 0; ilev < max_level; ilev++)
	{
	  fineDx /= a_params.refRatio[ilev];
	}	
      vector<int> plane_BC_hi(SpaceDim);
      vector<Real> domYInt(3,0);
      pp.getarr("plane_BC_hi",plane_BC_hi,0,SpaceDim);
      domYInt[0] = fineDx[1];
      domYInt[1] = plane_BC_hi[1]*fineDx[1];
      domYInt[2] = a_params.domainLength[1];
		
      int orderEB;
      pp.get("order_ebbc", orderEB);
      // pp.query("dom_order_eb",orderEB);
      Vector<Vector<int> >  a_domMixBc(2, vector<int>(SpaceDim));
      Vector<int> dom_mixbc_v(2*SpaceDim, 0); 
      pp.getarr("dom_mixbc",dom_mixbc_v,0,2*SpaceDim);
	  
      //this is defined with solver, check bool onlyHomog=false;
      for (int idir = 0; idir < SpaceDim; idir++)
	{
	  //1st index for face side (lo or hi), second for face direction (x or y)
	  a_domMixBc[0][idir] = dom_mixbc_v[2*idir];
	  a_domMixBc[1][idir] = dom_mixbc_v[2*idir+1];

	}
      
      vector<Real> domDxInt(2, 0); 
      vector<Real> EbDxInt(4, 0); 	
      PoissonDomMixBcFunc* DomDbdBcPtr = new PoissonDomMixBcFunc();
      DomDbdBcPtr->define(domDxInt, domYInt, domDirBcValue, domNeuBcValue, EbDxInt, EbDirBcValue, omega);
      RefCountedPtr<BaseMixBCValue> baseDomDbdBcPtr = RefCountedPtr<BaseMixBCValue>(DomDbdBcPtr);

      MixedPoissonDomainBCFactory* domainBCFactory = new MixedPoissonDomainBCFactory();     
      domainBCFactory->setEBOrder(orderEB);
      domainBCFactory->setArguments(a_domMixBc, domDirBcValue, domNeuBcValue, baseDomDbdBcPtr);
      domainBCFactory->setTime(&m_time);
      a_baseDomainBCFactory = RefCountedPtr<BaseDomainBCFactory>(domainBCFactory);
    }
  else if (a_params.domBcType == 0)
    {
      pout() << "constant Neumann bcs on domain" << endl;
      Real domBCValue;
      pp.get("dom_neu_bc_value", domBCValue);

      NeumannPoissonDomainBCFactory* domainBCFactory = new NeumannPoissonDomainBCFactory();
      domainBCFactory->setValue(domBCValue);

      a_baseDomainBCFactory = RefCountedPtr<BaseDomainBCFactory>(domainBCFactory);
    }
  else if (a_params.domBcType == 1)
    {
      pout() << "constant Dirichlet bcs on domain" << endl;
      DirichletPoissonDomainBCFactory* domainBCFactory = new DirichletPoissonDomainBCFactory();
      Real domBCValue;
      pp.get("domain_dir_bc_value", domBCValue);
      domainBCFactory->setValue(domBCValue);

      a_baseDomainBCFactory = RefCountedPtr<BaseDomainBCFactory>(domainBCFactory);
    }

  else
    {
      MayDay::Error("Unknown domain BC type");
    }

  //[MD] for mixed Dirichlet/Neumann EB BCs
  if (a_params.ebBcType == 100) 
    {
      int EbBcSel;
      Real EbDirBcValue=0, EbNeuBcValue=0, epsHi, epsLo, src, omega, frequency;
      pp.get("eb_dir_bc_value", EbDirBcValue);
      pp.get("frequency", frequency);
      omega = 2*M_PI*frequency;
      //MDCancel 06/13: only do this for DC case, for AC starting sine is used and below is not necessary
      //if(EbDirBcValue>0 && omega==0)
      //EbDirBcValue = Min(100 + (EbDirBcValue-100)/8*m_call,EbDirBcValue);
      //else if(EbDirBcValue<0 && omega==0)
      //EbDirBcValue = -Min(100 + (abs(EbDirBcValue)-100)/8*m_call,abs(EbDirBcValue));

	
      pout() << "mixed Dirichlet/Neumann bcs on EB,  EbDirBcValue =" << EbDirBcValue << endl;

      pp.get("eb_bc_sel", EbBcSel);
      //MDChange: now defined in constructor
      //pp.get("eps_gas", epsHi);
      //pp.get("eps_diel", epsLo);
      
      //pp.get("src", src);
      // defined by solver, redundancy, check 
      // bool onlyHomog=false;
      int orderEB;
      pp.get("order_ebbc", orderEB);
      vector<Real> EbDxInt(4, 0); 	 
      pp.getarr("eb_d_xint",EbDxInt,0,4);
 	  
      int max_level = a_params.maxLevel;
      RealVect fineDx = a_params.coarsestDx;
      for (int ilev = 0; ilev < max_level; ilev++)
	{
	  fineDx /= a_params.refRatio[ilev];
	}

      PoissonEBMixBcFunc* EbDbdBcPtr = new PoissonEBMixBcFunc();
      EbDbdBcPtr->define(EbDxInt, EbDirBcValue, EbNeuBcValue, omega);
      RefCountedPtr<BaseMixBCValue> baseEbDbdBcPtr = RefCountedPtr<BaseMixBCValue>(EbDbdBcPtr);

      MixedPoissonEBBCFactory* ebBC = new MixedPoissonEBBCFactory();     
      ebBC->setOrder(orderEB);
      ebBC->setTime(&m_time);
      //MD below m_epsd/g - also need ebBV->setArguments(rhos) with rhos at ilevPtr or something like that
      ebBC->setArguments(EbDirBcValue, EbNeuBcValue, baseEbDbdBcPtr, EbBcSel, fineDx, m_epsd, m_epsg);

      a_baseEBBCFactory = RefCountedPtr<BaseEBBCFactory>(ebBC);
    }
  else if (a_params.ebBcType == 0)
    {
      Real domBCValue;
      pp.get("eb_bc_value", domBCValue);
      NeumannPoissonEBBCFactory* ebBC = new NeumannPoissonEBBCFactory();
      ebBC->setValue(domBCValue);

      a_baseEBBCFactory =RefCountedPtr<BaseEBBCFactory>( ebBC);
    }
  else if (a_params.ebBcType == 1)
    {
      
      Real Potential=0, omega, frequency, sphereDistance, xloc;
      pp.get("Potential", Potential);
      pp.get("frequency", frequency);
      omega = 2*M_PI*frequency;
      pp.get("OSUsphere_distance",sphereDistance);
      xloc = sphereDistance/2.0;

      if(m_isDebug) Potential=0;
      pout() << " Dirichlet bcs on EB,  Potential =" << Potential << " xloc =" << xloc  << endl;
      
      int orderEB;
      pp.get("order_ebbc", orderEB);

      
 	  
      int max_level = a_params.maxLevel;
      RealVect fineDx = a_params.coarsestDx;
      for (int ilev = 0; ilev < max_level; ilev++)
	{
	  fineDx /= a_params.refRatio[ilev];
	}

      IonizationWaveEBDirBcFunc* EbDbdBcPtr = new IonizationWaveEBDirBcFunc();
      EbDbdBcPtr->define(Potential, omega, xloc, &m_time);
      RefCountedPtr<BaseBCValue> baseEbDbdBcPtr = RefCountedPtr<BaseBCValue>(EbDbdBcPtr);

      DirichletPoissonEBBCFactory* ebBC = new DirichletPoissonEBBCFactory();     
      ebBC->setOrder(orderEB);
      ebBC->setFunction(baseEbDbdBcPtr);
      a_baseEBBCFactory = RefCountedPtr<BaseEBBCFactory>(ebBC);
      
    }
  else
    {
      MayDay::Error("Unknown EB BC type");
    }
}

void ElectricPotential::getConductivityBCFactories(RefCountedPtr<BaseDomainBCFactory>&                     a_baseDomainBCFactory,
						   RefCountedPtr<BaseEBBCFactory>&                         a_baseEBBCFactory,
						   const Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&    a_bcoef,
						   const Vector<DisjointBoxLayout>&                        a_grids,
						   const Vector<EBISLayout>&                               a_ebisl,
						   const PoissonParameters&                                a_params)
{
  CH_TIME("ElectricPotential::getConductivityBCFactories");
  ParmParse pp;
  //[MD] for mixed Dirichlet/Neumann BCs
  if (a_params.domBcType == 100) 
    {
      
     pout() << "mixed Dirichlet/Neumann bcs on domain" << endl;
      Real domDirBcValue=0, domNeuBcValue=0, EbDirBcValue=0, omega, frequency;
      pp.get("dom_dir_bc_value", domDirBcValue);
      pp.get("dom_neu_bc_value", domNeuBcValue);
      pp.get("Potential", EbDirBcValue);
      pp.get("frequency", frequency);
      omega = 2*M_PI*frequency;

		
      int orderEB;
      pp.get("order_ebbc", orderEB);
      // pp.query("dom_order_eb",orderEB);
      Vector<Vector<int> >  a_domMixBc(2, vector<int>(SpaceDim));
      Vector<int> dom_mixbc_v(2*SpaceDim, 0); 
      pp.getarr("dom_mixbc",dom_mixbc_v,0,2*SpaceDim);
	  
      //this is defined with solver, check bool onlyHomog=false;
      for (int idir = 0; idir < SpaceDim; idir++)
	{
	  //1st index for face side (lo or hi), second for face direction (x or y)
	  a_domMixBc[0][idir] = dom_mixbc_v[2*idir];
	  a_domMixBc[1][idir] = dom_mixbc_v[2*idir+1];

	}
      
      int max_level = a_params.maxLevel;
      RealVect fineDx = a_params.coarsestDx;
      for (int ilev = 0; ilev < max_level; ilev++)
	{
	  fineDx /= a_params.refRatio[ilev];
	}
      vector<int> plane_BC_hi(SpaceDim);
      vector<Real> domYInt(3,0);
      vector<Real> domDxInt(2, 0); 
      vector<Real> EbDxInt(4, 0); 	
      PoissonDomMixBcFunc* DomDbdBcPtr = new PoissonDomMixBcFunc();
      DomDbdBcPtr->define(domDxInt, domYInt, domDirBcValue, domNeuBcValue, EbDxInt, EbDirBcValue, omega);
      RefCountedPtr<BaseMixBCValue> baseDomDbdBcPtr = RefCountedPtr<BaseMixBCValue>(DomDbdBcPtr);

      MixedConductivityDomainBCFactory* domainBCFactory = new MixedConductivityDomainBCFactory();     
      domainBCFactory->setEBOrder(orderEB);
      domainBCFactory->setArguments(a_domMixBc, domDirBcValue, domNeuBcValue, baseDomDbdBcPtr);
      domainBCFactory->setTime(&m_time);
      a_baseDomainBCFactory = RefCountedPtr<BaseDomainBCFactory>(domainBCFactory);
    }
  else if (a_params.domBcType == 0)
    {
      pout() << "constant Neumann bcs on domain" << endl;
      Real domBCValue;
      pp.get("dom_neu_bc_value", domBCValue);

      NeumannConductivityDomainBCFactory* domainBCFactory = new NeumannConductivityDomainBCFactory();
      domainBCFactory->setValue(domBCValue);

      a_baseDomainBCFactory = RefCountedPtr<BaseDomainBCFactory>(domainBCFactory);
    }
  else if (a_params.domBcType == 1)
    {
      pout() << "constant Dirichlet bcs on domain" << endl;
      DirichletConductivityDomainBCFactory* domainBCFactory = new DirichletConductivityDomainBCFactory();
      Real domBCValue;
      pp.get("dom_dir_bc_value", domBCValue);
      domainBCFactory->setValue(domBCValue);

      a_baseDomainBCFactory = RefCountedPtr<BaseDomainBCFactory>(domainBCFactory);
    }

  else
    {
      MayDay::Error("Unknown domain BC type");
    }

  //[MD] for mixed Dirichlet/Neumann EB BCs
  if (a_params.ebBcType == 100) 
    {
     
      MayDay::Error("Not Implemented yet for RZ problems");
    }
  else if (a_params.ebBcType == 0)
    {
      Real domBCValue;
      pp.get("eb_bc_value", domBCValue);
      NeumannConductivityEBBCFactory* ebBC = new NeumannConductivityEBBCFactory();
      ebBC->setValue(domBCValue);

      a_baseEBBCFactory =RefCountedPtr<BaseEBBCFactory>( ebBC);
    }
  else if (a_params.ebBcType == 1)
    {
      
      Real Potential=0, omega, frequency, sphereDistance, xloc;
      pp.get("Potential", Potential);
      pp.get("frequency", frequency);
      omega = 2*M_PI*frequency;
      pp.get("OSUsphere_distance",sphereDistance);
      xloc = sphereDistance/2.0;

	
      if(m_isDebug) {
	Potential=0;
	pout() << "zeroing out Potential for debug" << endl;
      }
      pout() << " Dirichlet bcs on EB,  Potential =" << Potential << " xloc =" << xloc  << endl;
      
      int orderEB=1;
      pp.get("order_ebbc", orderEB);

      
 	  
      int max_level = a_params.maxLevel;
      RealVect fineDx = a_params.coarsestDx;
      for (int ilev = 0; ilev < max_level; ilev++)
	{
	  fineDx /= a_params.refRatio[ilev];
	}

      IonizationWaveEBDirBcFunc* EbDbdBcPtr = new IonizationWaveEBDirBcFunc();
      EbDbdBcPtr->define(Potential, omega, xloc, &m_time);
      RefCountedPtr<BaseBCValue> baseEbDbdBcPtr = RefCountedPtr<BaseBCValue>(EbDbdBcPtr);

      DirichletConductivityEBBCFactory* ebBC = new DirichletConductivityEBBCFactory();     
      ebBC->setOrder(orderEB);
      ebBC->setFunction(baseEbDbdBcPtr);
      a_baseEBBCFactory = RefCountedPtr<BaseEBBCFactory>(ebBC);
      
    }
  else
    {
      MayDay::Error("Unknown EB BC type");
    }
}
/********/
void ElectricPotential::getPoissonParameters(PoissonParameters&  a_params, bool a_forceSingleLevel)
{
  //MD not much to clean up, more to come after finalizing input format
  CH_TIME("PoissonUtilities::getPoissonParameters");
  ParmParse pp;

  pp.get("eb_bc_type",a_params.ebBcType);
  pp.get("domain_bc_type",a_params.domBcType);

  std::vector<int> nCellsArray(SpaceDim);
  pp.getarr("n_cells",nCellsArray,0,SpaceDim);

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_params.nCells[idir] = nCellsArray[idir];
      a_params.ghostPhi[idir] = 4;
      a_params.ghostRHS[idir] = 0;
      if (a_params.ebBcType == 0 || a_params.ebBcType == 2)
        {
	  a_params.ghostPhi[idir] = 1;
        }
    }
  if (pp.contains("ghostPhi") && pp.contains("ghostRhs"))
    {
      vector<int> ghost_phi_v(SpaceDim, 4);
      vector<int> ghost_rhs_v(SpaceDim, 4);
      pp.queryarr("ghostPhi",ghost_phi_v,0,SpaceDim);
      pp.queryarr("ghostRhs",ghost_rhs_v,0,SpaceDim);

      for (int idir = 0; idir < SpaceDim; idir++)
        {
          a_params.ghostPhi[idir] = ghost_phi_v[idir];
          a_params.ghostRHS[idir] = ghost_rhs_v[idir];
        }
    }
  if (a_forceSingleLevel)
    {
      a_params.maxLevel  = 0;
      a_params.numLevels = 1;
      a_params.refRatio.resize(1,2);
      a_params.blockFactor = 8;
      a_params.fillRatio = 1;
      a_params.bufferSize = 123;
    }
  else
    {
      pp.get("max_level_species", a_params.maxLevel);
      a_params.numLevels = a_params.maxLevel + 1;
      pp.getarr("ref_ratio",a_params.refRatio,0,a_params.numLevels);
      pp.get("block_factor",a_params.blockFactor);
      pp.get("fill_ratio",a_params.fillRatio);
      pp.get("grid_buffer_size",a_params.bufferSize);
    }
  IntVect lo = IntVect::Zero;
  IntVect hi = a_params.nCells;
  hi -= IntVect::Unit;

  std::vector<int> is_periodica(SpaceDim,0);
  bool is_periodic[SpaceDim];
  pp.queryarr("is_periodic", is_periodica,0, SpaceDim);
  // convert periodic from int->bool
  for (int dim=0; dim<SpaceDim; dim++)
    {
      is_periodic[dim] = (is_periodica[dim] == 1);
      if (is_periodic[dim])
        {
          pout() << "Using Periodic BCs in direction: " << dim << endl;
        }
    }

  a_params.coarsestDomain = ProblemDomain(lo, hi,is_periodic);

  std::vector<Real> dLArray(SpaceDim);
  pp.getarr("domain_length",dLArray,0,SpaceDim);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_params.domainLength[idir] = dLArray[idir];
    }

  pp.get("which_geom",   a_params.whichGeom);
  pp.get("max_grid_size",a_params.maxGridSize);

  //derived stuff
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_params.coarsestDx[idir] = a_params.domainLength[idir]/a_params.nCells[idir];
    }
  a_params.probLo = RealVect::Zero;
  a_params.probHi = RealVect::Zero;
  a_params.probHi += a_params.domainLength;

  a_params.noRefCorners = false;
  pp.query("no_ref_corners", a_params.noRefCorners);
  
}

/********/

void ElectricPotential:: definePoissonGeometry(const PoissonParameters&  a_params)
{
  CH_TIME("PoissonUtilities::definePoissonGeometry");
  int max_level = a_params.maxLevel;

  ProblemDomain finestDomain = a_params.coarsestDomain;
  for (int ilev = 0; ilev <  max_level; ilev++)
    {
      finestDomain.refine(a_params.refRatio[ilev]);
    }

  int verbosity = m_verbosity;
  RealVect origin = RealVect::Zero;

  RealVect coarseDx = a_params.coarsestDx;
  int ebMaxCoarsen = -1;
  RealVect fineDx = coarseDx;
  for (int ilev = 0; ilev < max_level; ilev++)
    {
      fineDx /= a_params.refRatio[ilev];
    }
  int whichgeom = a_params.whichGeom;
  int ebMaxSize = a_params.maxGridSize;
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  EBIndexSpace* electricPtr = Electric_EBIS::instance();

  ParmParse pp;
//OSU geometry
  int numSpheres=2;
  Real sphereRadius=1e99;
  Real sphereDistance = 1e99;
  pp.get("OSUsphere_radius",sphereRadius);
  pp.get("OSUsphere_distance",sphereDistance);
  Vector<Real>     radius(numSpheres);
  Vector<RealVect> center(numSpheres);
  center[0][0] = 0;
  center[1][0] = sphereDistance;
  for (int isphere = 0; isphere < numSpheres; isphere++)
     {
       radius[isphere] = sphereRadius;
     }
  if(whichgeom == 300){ //Pin2Pin
    bool inside = false;
    MultiSphereIF spheres(radius, center, inside);
    GeometryShop workshop(spheres,verbosity,fineDx);
    ebisPtr->define(finestDomain, origin, fineDx[0], workshop, ebMaxSize, ebMaxCoarsen);
    electricPtr->define(finestDomain, origin, fineDx[0], workshop, ebMaxSize, ebMaxCoarsen);
  }
  else if (whichgeom == 400){  //PlateToPlate
    bool inside = true;
    RealVect normal1(D_DECL(1.0,0.0,0.0));
    PlaneIF plane1(normal1,center[0] + normal1*radius[0]*(0),inside);
    RealVect normal2(D_DECL(-1.0,0.0,0.0));
    PlaneIF plane2(normal2,center[1] + normal2*radius[1]*(0),inside);
    IntersectionIF spheres(plane1,plane2);
    GeometryShop workshop(spheres,verbosity,fineDx);
    ebisPtr->define(finestDomain, origin, fineDx[0], workshop, ebMaxSize, ebMaxCoarsen);
    electricPtr->define(finestDomain, origin, fineDx[0], workshop, ebMaxSize, ebMaxCoarsen);
  }
    
}  

/****************************/
void ElectricPotential::tagCells(IntVectSet& a_tags, int a_ilev)
{
  CH_TIME("EBAMRSpecies::tagCells");
  if (m_verbosity >= 3)
    {
      pout() << "EBAMRSpecies::tagCells for level " << a_ilev << "Threshold " << m_refineThresh  << endl;
    }
  // LM tagging the gradinet of phi
  a_tags.makeEmpty();

  // If there is a coarser level interpolate undefined ghost cells
  //only interpolate the Electric field
  EBCellFactory factory(m_ebisl[a_ilev]);
  int nCons = 1;
  int nghost = m_nGhost;
  int phiIndex = 0;
  Interval consInterv(0, nCons-1);
  Interval intervDensity = consInterv;

  /*
    LevelData<EBCellFAB> consTemp(m_grids[a_ilev], nCons, (IntVect::Unit)*nghost, factory);
    m_phi[a_ilev]->copyTo(consInterv, consTemp, consInterv);
    if (a_ilev>0)
    {
    int refRatCrse = m_params.refRatio[a_ilev-1];
    EBPWLFillPatch patcher(m_grids[a_ilev],
    m_grids[a_ilev-1],
    m_ebisl[a_ilev],
    m_ebisl[a_ilev-1],
    m_domain[a_ilev-1].domainBox(),
    refRatCrse, m_nComp, nghost);
      
    Real coarTimeOld = 0.0;
    Real coarTimeNew = 1.0;
    Real fineTime    = 1.0;
    patcher.interpolate(consTemp,
    *m_phi[a_ilev-1],
    *m_phi[a_ilev-1],
    coarTimeOld,
    coarTimeNew,
    fineTime,
    intervDensity);
    }
    consTemp.exchange(intervDensity);*/

  IntVectSet localTags;
  // Compute undivided gradient
  for (DataIterator dit = m_grids[a_ilev].dataIterator(); dit.ok(); ++dit)
    {
      const Box& b = m_grids[a_ilev].get(dit());
      const EBISBox& ebisBox = m_ebisl[a_ilev][dit()];
      const EBCellFAB& gradFAB = (*m_cellGradPhi[a_ilev])[dit()];
      

      //pointwise op so just have to iterate over multivalued cells
      

      // Tag where gradient exceeds threshold

      IntVectSet ivsTot(b);
      for (VoFIterator vofit(ivsTot, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
	{
	  const VolIndex& vof = vofit();
	  const IntVect& iv = vof.gridIndex();
	  Real gradmag = 0.0;
	  
	  int tagGradComps = SpaceDim;
	  for (int icomp = 0; icomp < tagGradComps; icomp++)
	    {
	      Real gradDirVal = gradFAB(vof, icomp);
	      gradmag += gradDirVal*gradDirVal;
	    }
	  gradmag = sqrt(gradmag);
	  if (gradmag >= m_refineThresh) localTags |= iv;
	  
	  //IntVectSet irregIVS = ebisBox.getIrregIVS(b);
	  //localTags |= irregIVS;
	}
    }
  
  localTags.grow(m_tagBufferSize);
  
  // Need to do this in two steps unless a IntVectSet::operator &=
  // (ProblemDomain) operator is defined
  Box localTagsBox = localTags.minBox();
  localTagsBox &= m_domain[a_ilev];
  localTags &= localTagsBox;
  a_tags = localTags;
  
}

/***************************/
void ElectricPotential::regridWithDBL(const Vector<DisjointBoxLayout>& a_DBL)
{
  CH_TIME("ElectricPotential::regrid");

  int finest_level_old = m_grids.size()-1;
  
  m_params.numLevels = a_DBL.size();
  m_params.maxLevel = m_params.numLevels -1;

  m_grids.resize(m_params.numLevels);
  m_ebisl.resize(m_params.numLevels);
  m_eblg.resize(m_params.numLevels);
  m_phi.resize(m_params.numLevels,NULL);
  m_rhs.resize(m_params.numLevels,NULL);
  m_gradPhi.resize(m_params.numLevels,NULL);
  m_cellGradPhi.resize(m_params.numLevels,NULL);
  m_IVGradPhi.resize(m_params.numLevels,NULL);

  const EBIndexSpace* const ebisPtr = Electric_EBIS::instance();
  // save data for later copy
  //not using m_ebisl because it gets wiped later
  //the ebisl has to know about the fact that we really have
  //four ghost cells.
  int nGhostEBISL = 6;
  Interval interv(0,m_nComp-1);
  IntVect ivGhost = m_nGhost*IntVect::Unit;

  for (int ilev = 0; ilev < a_DBL.size(); ilev++)
    {
      LevelData<EBCellFAB> stateSaved;

      if(ilev <= finest_level_old)
	{
	  EBISLayout ebislOld;
	  ebisPtr->fillEBISLayout(ebislOld, m_grids[ilev], m_domain[ilev].domainBox(), nGhostEBISL);
	  EBCellFactory factoryOld(ebislOld);
	  stateSaved.define(m_grids[ilev], m_nComp, ivGhost, factoryOld);
	  m_phi[ilev]->copyTo(interv, stateSaved, interv);
	}

  

      m_grids[ilev] = a_DBL[ilev];
      ebisPtr->fillEBISLayout(m_ebisl[ilev], m_grids[ilev], m_domain[ilev].domainBox(), nGhostEBISL);
      m_eblg[ilev] = EBLevelGrid(m_grids[ilev], m_domain[ilev], nGhostEBISL, ebisPtr);
      //m_eblg[ilev] = EBLevelGrid(m_grids[ilev], m_ebisl[ilev], m_domain[ilev]);

      EBCellFactory factoryNew(m_ebisl[ilev]);
      // reshape state with new grids
      if(m_phi[ilev] == NULL) m_phi[ilev]   = new LevelData<EBCellFAB>();
      if(m_rhs[ilev] == NULL) m_rhs[ilev]   = new LevelData<EBCellFAB>();
      if(m_gradPhi[ilev] == NULL) m_gradPhi[ilev]   = new LevelData<EBFluxFAB>();
      if(m_cellGradPhi[ilev] == NULL) m_cellGradPhi[ilev]   = new LevelData<EBCellFAB>();
      if(m_IVGradPhi[ilev] == NULL)   m_IVGradPhi[ilev]     = new LevelData<BaseIVFAB<Real> >();
      m_phi[ilev]->define(m_grids[ilev],m_nComp,ivGhost, factoryNew);
      m_rhs[ilev]->define(m_grids[ilev],m_nComp,m_params.ghostRHS, factoryNew);


      // interpolate to coarser level
      if (ilev > 0 )
	{
	  EBPWLFineInterp ebInterpVec(m_grids[ ilev  ],
				      m_grids[ ilev-1],
				      m_ebisl[ ilev  ],
				      m_ebisl[ ilev-1],
				      m_domain[ilev-1],
				      m_params.refRatio[ilev-1],
				      1,
				      ebisPtr);
	  ebInterpVec.interpolate(*m_phi[ilev], *m_phi[ilev-1], interv);
	}
      
      // copy from old state
      if(ilev <= finest_level_old)
	stateSaved.copyTo(interv,*m_phi[ilev], interv);
    }
  
  EBAMRDataOps::quadCFInterpAll(m_phi,m_eblg,m_params.refRatio,ebisPtr);
  //EBAMRDataOps::pwlFillPatchAll(m_phi,m_eblg,m_params.refRatio,ebisPtr); //cant patch askMD if she can patch her Leveldatas
  EBAMRDataOps::averageDown(m_phi,m_eblg,m_params.refRatio);
  
  defineSolver();
 
  getGradients();
}
void ElectricPotential::outputHDF(const Vector<LevelData<EBCellFAB>* >& a_Ptr, const string& a_name)
{


  string  filename = a_name;
  std::ostringstream oss;
  oss << m_call;
  filename += oss.str();
  filename.append(".hdf5");

  Real dumReal =  1.0;

  //EBAMRDataOps::checkNANBOX(a_Ptr);
  long long totalPoints = 0;
  long long totalBoxes  = 0;
  for (int ilev = 0; ilev < m_params.numLevels; ilev++)
    {
      //m_grids[ilev].print();// print all boxes to cout
      long long pointsThisLevel = 0;
      for (LayoutIterator lit = m_grids[ilev].layoutIterator(); lit.ok(); ++lit)
	{
	  pointsThisLevel += m_grids[ilev][lit()].numPts();
	}
	  
      totalPoints += pointsThisLevel;
      totalBoxes += m_grids[ilev].size();
      pout() << "ElectricPotential:level[" << ilev
	     << "], number of boxes = " << m_grids[ilev].size()
	     << ", number of points = " << pointsThisLevel << endl;
    }
  pout() << "ElectricPotential:"
	 <<  "   total boxes = " << totalBoxes
	 <<  ", total points = " << totalPoints <<  endl;
  pout() << "Coarsest Domain Box " << m_params.coarsestDomain.domainBox() << endl;
  pout() << "Coarsest Dx " << m_dx[0] << endl;




  // calculate the number of boxes
  long long localPoints = 0;
  for (int ilev = 0; ilev < m_params.numLevels; ilev++)
    {
      long long pointsThisLevel = 0;
      for (DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
	  Box thisB = m_grids[ilev][dit()]  ;// & m_jetBox[ilev];
	  const EBISBox& ebisBox = (*a_Ptr[ilev])[dit()].getEBISBox();
          pointsThisLevel += ebisBox.getEBGraph().getMultiCells(thisB).numPts();
        }
      localPoints += pointsThisLevel;
    }
  long long multipoints = EBLevelDataOps::parallelSum(static_cast<int>(localPoints));
  pout() << a_name << " multipoints " << multipoints << endl;

  int pvars = a_Ptr[0]->nComp();
  char charstrrhs[100];
  Vector<string> namesrhs(pvars, string("var"));
  for (int idir = 0; idir < pvars; idir++)
    {
      sprintf(charstrrhs, "var%d",idir+1);
      namesrhs[idir] = charstrrhs;
    }
  bool replaceCovered = false;
  Vector<Real> coveredValues(pvars, -10.0);
  writeEBHDF5(filename, m_grids, a_Ptr, namesrhs,
	      m_params.coarsestDomain.domainBox(), m_dx[0], dumReal, m_time,
	      m_params.refRatio, m_params.numLevels,
	      replaceCovered, coveredValues);
  pout() << " End Output HDF for " << a_name << endl;
}
//**************
void ElectricPotential::setTime(Real a_time)
{
  m_time=a_time;
  if(m_verbosity > 2) pout() << " ElectricPotential::setTime "<< a_time << endl;
}
//**************
void ElectricPotential::getGradients()
{
  // run Interpolation and Exchanges on phi (lucacancel these should have already have been performed by the elliptic solver-delete next 5 lines for speed-)
  {const EBIndexSpace* const ebisPtr = Electric_EBIS::instance();
  EBAMRDataOps::exchangeAll(m_phi);
  EBAMRDataOps::averageDown(m_phi,m_eblg,m_params.refRatio);
  EBAMRDataOps::quadCFInterpAll(m_phi,m_eblg,m_params.refRatio,ebisPtr);
  //EBAMRDataOps::pwlFillPatchAll(m_phi,m_eblg,m_params.refRatio,ebisPtr);
  }

  // get the gradients
  for (int ilev = 0; ilev < m_grids.size(); ++ilev)
    {
      if(m_verbosity > 3)pout() << ilev << " Evaluating EPE gradients " << endl;
      const ProblemDomain& domain = m_eblg[ilev].getDomain();
      RealVect dxv = m_dx[ilev]*RealVect::Unit;
      EBFluxFactory fact(m_eblg[ilev].getEBISL());
      m_gradPhi[ilev]->define(m_grids[ilev], 1, 2*(IntVect::Unit), fact);
      LevelData<EBFluxFAB>& gradPhi = *m_gradPhi[ilev];
      EBLevelDataOps::setToZero(gradPhi);
      for (DataIterator dit = m_eblg[ilev].getDBL().dataIterator(); dit.ok(); ++dit)
	{ 
	  Box dblBox = m_eblg[ilev].getDBL()[dit()];
	  const EBISBox&  ebisBox = m_eblg[ilev].getEBISL()[dit()];
	  macGradient(gradPhi[dit()], (*m_phi[ilev])[dit()],
		      ebisBox, dblBox, domain, dxv);
	}
      ccpExtrapolateToDomainBoundaries(gradPhi,m_grids[ilev],m_ebisl[ilev],domain, dxv);
    }
  EBAMRDataOps::exchangeAll(m_gradPhi);
  EBAMRDataOps::averageDown(m_gradPhi,m_eblg,m_params.refRatio);
  EBAMRDataOps::setToZero(m_cellGradPhi);
  for (int ilev = 0; ilev < m_grids.size(); ++ilev)
    {    
      const ProblemDomain& domain = m_eblg[ilev].getDomain();
      RealVect dxv = m_dx[ilev]*RealVect::Unit;  
      EBCellFactory cfac(m_eblg[ilev].getEBISL());
      m_cellGradPhi[ilev]->define(m_grids[ilev], SpaceDim, m_params.ghostPhi, cfac);
      LevelData<EBFluxFAB>& gradPhi = *m_gradPhi[ilev];
      ccpAverageFaceToCells(*m_cellGradPhi[ilev],gradPhi,m_grids[ilev],m_ebisl[ilev],domain,dxv);
      //cellGrad(*m_cellGradPhi[ilev],*m_phi[ilev],0,m_grids[ilev],m_ebisl[ilev],domain,dxv);
    }

  if(m_verbosity > 3)pout() << " EPE exchanges and CFBC " << endl;
  const EBIndexSpace* const ebisPtr = Electric_EBIS::instance();
  EBAMRDataOps::exchangeAll(m_cellGradPhi);
  EBAMRDataOps::averageDown(m_cellGradPhi,m_eblg,m_params.refRatio);
  EBAMRDataOps::pwlFillPatchAll(m_cellGradPhi,m_eblg,m_params.refRatio,ebisPtr);
  EBAMRDataOps::quadCFInterpAll(m_cellGradPhi,m_eblg,m_params.refRatio,ebisPtr);

  // Evaluate IV
  for (int ilev = 0; ilev < m_grids.size(); ++ilev)
    {    
      const ProblemDomain& domain = m_eblg[ilev].getDomain();
      RealVect dxv = m_dx[ilev]*RealVect::Unit;
      LayoutData<IntVectSet> irregSets(m_grids[ilev]);
      for (DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit) 
	{
	  Box grownBox = grow(m_grids[ilev].get(dit()), m_params.ghostPhi); 
	  grownBox &= m_domain[ilev];
	  irregSets[dit()] = m_ebisl[ilev][dit()].getIrregIVS(grownBox);
	}
      BaseIVFactory<Real>  baseivfact(m_ebisl[ilev], irregSets);
      m_IVGradPhi[ilev]->define(m_grids[ilev], 1    , m_params.ghostPhi, baseivfact);

      LevelData<BaseIVFAB<Real>>& IVPhi = *m_IVGradPhi[ilev];
      evaluateIVGrad(*m_phi[ilev],*m_cellGradPhi[ilev],IVPhi,m_grids[ilev],m_ebisl[ilev],domain,dxv,ilev);
      //cellGrad(*m_cellGradPhi[ilev],*m_phi[ilev],0,m_grids[ilev],m_ebisl[ilev],domain,dxv);
      m_IVGradPhi[ilev]->exchange(Interval(0,0));
    }



}

void ElectricPotential::
cellGrad(LevelData<EBCellFAB>&         a_gradPhi,
         const LevelData<EBCellFAB>&   a_phi,
         const int&                    a_comp,
	 const DisjointBoxLayout &     a_grids,
	 const EBISLayout &            a_ebisl,
	 const ProblemDomain &         a_domain,
	 const RealVect &              a_dx)
{
  int ibox = 0;
  EBLevelDataOps::setToZero(a_gradPhi); 
  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit, ibox++)
    {
      EBCellFAB& gradPhi = a_gradPhi[dit()];
      const EBCellFAB& phi = a_phi[dit()];

      const Box& cellBox = a_grids.get(dit());
      const EBISBox& ebisBox = a_ebisl[dit()];
      IntVectSet ivsIrreg = ebisBox.getIrregIVS(cellBox);
      
      for (int derivDir = 0; derivDir < SpaceDim; derivDir++)
	{
	  Box loBox, hiBox, centerBox;
	  int hasLo, hasHi;
	  EBArith::loHiCenter(loBox, hasLo, hiBox, hasHi, centerBox, a_domain , cellBox, derivDir);

	  
          //
          BaseFab<Real>&       regGrad = gradPhi.getSingleValuedFAB();
          const BaseFab<Real>& regPhi  =     phi.getSingleValuedFAB();
          FORT_CELLGRADEBSPEC(CHF_FRA1(regGrad, derivDir),
                              CHF_CONST_FRA1(regPhi, a_comp),
                              CHF_CONST_REAL(a_dx[derivDir]),
                              CHF_BOX(loBox),
                              CHF_BOX(hiBox),
                              CHF_BOX(centerBox),
                              CHF_CONST_INT(hasLo),
                              CHF_CONST_INT(hasHi),
                              CHF_CONST_INT(derivDir));
	  
	  for (VoFIterator vofit(ivsIrreg, ebisBox.getEBGraph());   vofit.ok(); ++vofit)
	    {
	      const VolIndex& vof =  vofit();
	      Vector<FaceIndex> loFaces = ebisBox.getFaces(vof, derivDir, Side::Lo);
	      Vector<FaceIndex> hiFaces = ebisBox.getFaces(vof, derivDir, Side::Hi);
	      bool hasLoPt = ((loFaces.size() == 1)  && (!loFaces[0].isBoundary()));
	      bool hasHiPt = ((hiFaces.size() == 1)  && (!hiFaces[0].isBoundary()));
	      Real valLo=0, valHi=0;
	      if (hasLoPt) valLo = phi(loFaces[0].getVoF(Side::Lo), a_comp);
	      if (hasHiPt) valHi = phi(hiFaces[0].getVoF(Side::Hi), a_comp);
	      
	      Real valCe = phi(vof, a_comp);
	      
	      if (hasLoPt && hasHiPt)
		{
		  gradPhi(vof, derivDir) = (valHi - valLo)/(2.*a_dx[derivDir]);
		}
	      else if (hasHiPt)
		{
		  gradPhi(vof, derivDir) = (valHi - valCe)/(a_dx[derivDir]);
		}
	      else if (hasLoPt)
		{
		  gradPhi(vof, derivDir) = (valCe - valLo)/(a_dx[derivDir]);
		}
	      else
		{
		  gradPhi(vof, derivDir) = 0.0;
		}
	    }
	}
    }
  a_gradPhi.exchange();
}

/********************/
void ElectricPotential::writeCheckpointFile(HDF5Handle& a_handle) const
{
  
  for (int ilev = 0; ilev <= m_params.maxLevel; ilev++)
    { 
      // Setup the level string
      char levelStr[20];
      sprintf(levelStr,"%d",ilev);
      const std::string label = std::string("SPlevel_") + levelStr;

      a_handle.setGroup(label);
      write(a_handle,m_grids[ilev]);
      write(a_handle,*m_phi[ilev],"phi");
    }
}

/********************************************************/
Real ElectricPotential::getElectrodePotential() const
{ 
  if(m_omegaEB <=0)
    return m_voltageEB*cos(m_omegaEB*m_time);
  else
    return m_voltageEB*sin(m_omegaEB*m_time);

}
Real ElectricPotential::getYEB() const
{ 

  if(!m_isGeometryDefined)
    MayDay::Error("file does not contain state data");
  return m_yEB;
}
bool ElectricPotential::doTagging() const
{ 

  return m_refineThresh >0;
}
/**/ 
void ElectricPotential::zeroField(Vector<LevelData<EBCellFAB>* > a_cellData)
{ 

  //lucacancelthis block
  Real maxPot = getElectrodePotential()*1.2;
  Real maxVal=maxPot;Real minVal=-maxPot;int comp=0;
  Real maxphi, minphi; EBAMRDataOps::getMaxMin(maxphi, minphi, m_phi,0);//lucacancel
  maxphi = max(abs(maxphi),abs(minphi));

  if(abs(maxPot) > 1500 && maxphi > 1.2*abs(maxPot))
    {
      pout() << "Warning:: limiting potential" << endl;
      int numLevels = m_phi.size();
      for (int ilev = 0; ilev < numLevels; ilev++)
	{
	  RealVect dxv = m_dx[ilev]*RealVect::Unit;
	  for (DataIterator dit=m_phi[ilev]->dataIterator();dit.ok();++dit)
	    {
	      EBCellFAB& dataEBFAB = (*m_phi[ilev])[dit()];
	      const Box& region = (*a_cellData[ilev])[dit()].getRegion();
	      const IntVectSet ivsBox(region);
	      const EBISBox& ebisBox = m_ebisl[ilev][dit()];
	      //const EBISBox& ebisBox = dataEBFAB.getEBISBox();
	      for (VoFIterator vofit(ivsBox, ebisBox.getEBGraph());vofit.ok(); ++vofit)
		{
		  const VolIndex& vof = vofit();
		  RealVect xloc =  EBArith::getVofLocation(vof,dxv,RealVect::Zero) ;
		  Real val = dataEBFAB(vof,comp);
		  if(xloc[0] > 0.04 || val > maxVal || val < minVal)
		    {
		      for (int k = 0; k < a_cellData[ilev]->nComp(); ++k)
			{
			  (*a_cellData[ilev])[dit()](vof,k) = 0;

			}
		      
		    }
		}
	    }
	}
    }
}

/**/ 
void ElectricPotential::zeroSurf(Vector<LevelData<BaseIVFAB<Real> >* > a_cellData)
{ 

  int numLevels = a_cellData.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      RealVect dxv = m_dx[ilev]*RealVect::Unit;
      for (DataIterator dit=a_cellData[ilev]->dataIterator();dit.ok();++dit)
	{
	  const EBGraph&  ebgraph = m_ebisl[ilev][dit()].getEBGraph();
	  for (VoFIterator vofit((*a_cellData[ilev])[dit()].getIVS(), ebgraph); vofit.ok(); ++vofit)
	    {
	      VolIndex vof = vofit();
	      RealVect xloc =  EBArith::getVofLocation(vof,dxv,RealVect::Zero) ;
	      if(xloc[0] > m_params.domainLength[0]*0.9 || xloc[0] < 0.1* m_params.domainLength[0])
		for (int k = 0; k < a_cellData[ilev]->nComp(); ++k)
		  (*a_cellData[ilev])[dit()](vof,k) = 0;
	    }
	}
    }
      //}
}


//**************
void ElectricPotential::multiplyRHSByRadius()
{

  RealVect cylinderAxis = BASISREALV(0);
  int icomp = 0;
  for (int ilev = 0; ilev < m_rhs.size(); ilev++)
    {
      Real dx = m_dx[ilev];
      RealVect vectDx = dx*(RealVect::Unit);
      for (DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
	  
          EBCellFAB& rhsFAB = (*m_rhs[ilev])[dit()];
          const Box& region = rhsFAB.getRegion();
          IntVectSet ivsBox(region);
          const EBISBox& ebisBox = rhsFAB.getEBISBox();
          for (VoFIterator vofit(ivsBox, ebisBox.getEBGraph());vofit.ok(); ++vofit)
            {
              const VolIndex& vof = vofit();
	      RealVect centroid = ebisBox.centroid(vof);
	      centroid *= dx;
              RealVect vofloc = EBArith::getVofLocation(vofit(), vectDx, centroid);
              Real dist = getDistanceFromAxis(vofloc, cylinderAxis, RealVect::Zero);
	      rhsFAB(vof,icomp) *= dist;
            }
        }
      
    }

}

//**************
void ElectricPotential::testRHS()
{

  RealVect cylinderAxis = BASISREALV(0);
  int icomp = 0;
  for (int ilev = 0; ilev < m_rhs.size(); ilev++)
    {
      Real dx = m_dx[ilev];
      RealVect vectDx = dx*(RealVect::Unit);
      for (DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
	  
          EBCellFAB& rhsFAB = (*m_rhs[ilev])[dit()];
          const Box& region = rhsFAB.getRegion();
          IntVectSet ivsBox(region);
          const EBISBox& ebisBox = rhsFAB.getEBISBox();
          for (VoFIterator vofit(ivsBox, ebisBox.getEBGraph());vofit.ok(); ++vofit)
            {
              const VolIndex& vof = vofit();
	      RealVect centroid = ebisBox.centroid(vof);
	      centroid *= dx;
              RealVect vofloc = EBArith::getVofLocation(vofit(), vectDx, centroid);
	      Real x = vofloc[0];
              Real r = getDistanceFromAxis(vofloc, cylinderAxis, RealVect::Zero);
	      Real lap= -(20*M_PI*(-24 + 50*x + 25*(225*pow(r,4) + pow(x,2)*(23 + 25*(-2 + x)*x) + 5*pow(r,2)*(23 + 50*(-1 + x)*x)))*sin((5*M_PI*r)/4.)*sin(M_PI*x) +  r*cos((5*M_PI*r)/4.)*(-1600*M_PI*(-1 + 2*x)*(-1 + 25*pow(r,2) + 25*(-1 + x)*x)*cos(M_PI*x) + (-800*(69 + 250*pow(r,2) + 250*(-1 + x)*x) +  41*pow(M_PI,2)*(-1 + 25*pow(r,2) + 25*pow(x,2))* (24 + 25*pow(r,2) - 50*x + 25*pow(x,2)))*sin(M_PI*x)))/(10000.*r);
	      rhsFAB(vof,icomp) = lap;
            }
        }
      if (m_isKappaWeight) EBLevelDataOps::kappaWeight(*m_rhs[ilev]);
    }

}

//**************
void ElectricPotential::testSOL(Vector< LevelData<EBCellFAB>* > a_resid) const
{

  

  RealVect cylinderAxis = BASISREALV(0);
  int icomp = 0;
  for (int ilev = 0; ilev < a_resid.size(); ilev++)
    {
      Real dx = m_dx[ilev];
      RealVect vectDx = dx*(RealVect::Unit);
      for (DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
	  
          EBCellFAB& rhsFAB = (*a_resid[ilev])[dit()];
          EBCellFAB& phiFAB = (*m_phi[ilev])[dit()];
          const Box& region = rhsFAB.getRegion();
          IntVectSet ivsBox(region);
          const EBISBox& ebisBox = rhsFAB.getEBISBox();
          for (VoFIterator vofit(ivsBox, ebisBox.getEBGraph());vofit.ok(); ++vofit)
            {
              const VolIndex& vof = vofit();
	      RealVect centroid = ebisBox.centroid(vof);
	      centroid *= dx;
              RealVect vofloc = EBArith::getVofLocation(vofit(), vectDx, centroid);
	      Real x = vofloc[0];
              Real r = getDistanceFromAxis(vofloc, cylinderAxis, RealVect::Zero);
	      Real sol= (0.96 + pow(r,2) + (-2 + x)*x)*(-0.04 + pow(r,2) + pow(x,2))*cos((5*M_PI*r)/4.)*sin(M_PI*x);
	      rhsFAB(vof,icomp) = sol-phiFAB(vof,icomp);
            }
        }
    }

}


//**************
void ElectricPotential::extrpolateSolToCentroids()
{

  RealVect cylinderAxis = BASISREALV(0);
  int icomp = 0;
  for (int ilev = 0; ilev < m_phi.size(); ilev++)
    {
      Real dx = m_dx[ilev];
      RealVect vectDx = dx*(RealVect::Unit);
      for (DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
	  
          EBCellFAB& phiFAB = (*m_phi[ilev])[dit()];
          const Box& region = grow(m_grids[ilev].get(dit()),1) & m_domain[ilev];
          const EBISBox& ebisBox = phiFAB.getEBISBox();
          //IntVectSet ivsBox(region);
	  IntVectSet ivsBox = ebisBox.getIrregIVS(region);
          for (VoFIterator vofit(ivsBox, ebisBox.getEBGraph());vofit.ok(); ++vofit)
            {
              const VolIndex& vof = vofit();
	      //RealVect centroid = ebisBox.bndryCentroid(vof);
	      RealVect centroid = ebisBox.centroid(vof);
	      centroid *= dx;
	      VoFStencil extrapSten;
	      int order = EBArith::getExtrapolationStencil(extrapSten, centroid, vectDx, vofit(), ebisBox);
	      Real extrapval = applyVoFStencil(extrapSten, phiFAB, 0);
	      //pout() << extrapval-phiFAB(vof,icomp) << centroid << endl;
	      phiFAB(vof,icomp) = extrapval;
            }
        }
    }

}

/*****/
void ElectricPotential::
evaluateIVGrad(LevelData<EBCellFAB> &        a_cellData,
	       LevelData<EBCellFAB> &        a_gradData,
	       LevelData<BaseIVFAB<Real>> &  a_IVData,
	       const DisjointBoxLayout &     a_grids,
	       const EBISLayout &            a_ebisl,
	       const ProblemDomain &         a_domain,
	       const RealVect &              a_dxv,
	       const int &                   a_ilev)
{
  std::ofstream ofs;
  if (m_isDebug)
    {
      ofs.open ("g.dat", std::ofstream::out | std::ofstream::app);
    }
  
  int comp = 0;
  Vector< AMRLevelOp<LevelData<EBCellFAB> >* > ops = m_solver.getAMROperators();
  AMRLevelOp<LevelData<EBCellFAB> >* opA = &(*ops[a_ilev]);
  EBAMRPoissonOpRZ* op = dynamic_cast<EBAMRPoissonOpRZ* >(opA);
  LayoutData<BaseIVFAB<VoFStencil> >& poissSten = *(op->getFluxStencil(-1));
  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      const EBCellFAB& curPhi = a_cellData[dit()];
      const EBISBox& curEBISBox = a_ebisl[dit()];
      const EBGraph& ebgraph = curEBISBox.getEBGraph();
      const EBISBox& ebisBox = a_cellData[dit()].getEBISBox();  
      const Box grid = a_grids.get(dit());
      IntVectSet notRegular;
      notRegular |= curEBISBox.getIrregIVS  (grid);
      notRegular |= curEBISBox.getMultiCells(grid);
      
      BaseIVFAB<Real>& IVData = a_IVData[dit()];
      const BaseIVFAB<VoFStencil>& STData = poissSten[dit()];
      const BaseIVFAB<Real>& poissWeight = (op->getFluxWeight())[dit()];

       const IntVectSet& ivsIrreg = IVData.getIVS();
       const IntVectSet& ivs = STData.getIVS();
       //pout() << ivsIrreg << endl;
       //pout() << ivs << endl;
       //pout() << notRegular << endl;
      for (VoFIterator vofit(notRegular, ebgraph); vofit.ok(); ++vofit)
	{
	  VolIndex vof = vofit();
	  RealVect normal = ebisBox.normal(vof);
	  const IntVect& iv = vof.gridIndex();
	  
       
	  Real bValue =op->getEBValue(vof, dit(), m_params.probLo, curEBISBox, m_time);
          Real areaFrac = curEBISBox.bndryArea(vof);
	  if(areaFrac > 1e-16) areaFrac = a_dxv[0]/areaFrac;
	  VoFStencil curStencil = STData(vof, 0);
          curStencil *= areaFrac;
	  Real curWeight = poissWeight(vof, 0);
	  Real flux = curWeight * bValue;
	  for (int i = 0; i < curStencil.size(); i++)
	    {
	      const VolIndex& curVoF = curStencil.vof(i);
	      Real weight = curStencil.weight(i);
	      Real phiVal = curPhi(curVoF,comp);
	      flux += weight * phiVal;
	    }
	  
	  IVData(vof,0) = flux;
	  if(m_isDebug)
	    {
	      RealVect cylinderAxis = BASISREALV(0);
	      RealVect bndryCentroid = curEBISBox.bndryCentroid(vof);
	      bndryCentroid *= a_dxv[0];
	      RealVect vofloc = EBArith::getVofLocation(vof, a_dxv, bndryCentroid);
	      Real x = vofloc[0];
              Real r = getDistanceFromAxis(vofloc, cylinderAxis, RealVect::Zero);
	      Real Grad =0;
	      if(x< m_params.domainLength[0]/2)
		Grad =(-5*M_PI*r*(-1 + 25*pow(r,2) + 25*pow(x,2))*(25*pow(r,2) + (-6 + 5*x)*(-4 + 5*x))*sin((5*M_PI*r)/4.)*sin(M_PI*x) + 4*cos((5*M_PI*r)/4.)*(M_PI*x*(-1 + 25*pow(r,2) + 25*pow(x,2))*(24 + 25*pow(r,2) - 50*x + 25*pow(x,2))*cos(M_PI*x) + 50*(x + (pow(r,2) + pow(x,2))*(23 + 50*pow(r,2) + 25*x*(-3 + 2*x)))*sin(M_PI*x)))/(2500.*sqrt(pow(r,2) + pow(x,2)));
	      else
		Grad = (-5*M_PI*r*(-1 + 25*pow(r,2) + 25*pow(x,2))*(25*pow(r,2) + (-6 + 5*x)*(-4 + 5*x))*sin((5*M_PI*r)/4.)*sin(M_PI*x) + 4*cos((5*M_PI*r)/4.)*(M_PI*(-1 + x)*(-1 + 25*pow(r,2) + 25*pow(x,2))*(25*pow(r,2) + (-6 + 5*x)*(-4 + 5*x))*cos(M_PI*x) + 50*(50*pow(r,4) + (-1 + x)*(-1 + 2*x)*(-1 + 25*(-1 + x)*x) + pow(r,2)*(48 + 25*x*(-5 + 4*x)))*sin(M_PI*x)))/(2500.*sqrt(pow(r,2) + pow(-1 + x,2)));
	      //pout() << vofloc << " value " << IVData(vof,0) << " Grad " << Grad << endl;
	      ofs << a_ilev+ vofloc[0] << "  " << IVData(vof,0) << "  " << Grad  << "  " << bValue << endl;
	    }
  
	}
      
    }
}

#include "NamespaceFooter.H"
