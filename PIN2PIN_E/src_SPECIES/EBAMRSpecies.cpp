#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iomanip>

#include "parstream.H"

#include "LayoutIterator.H"
#include "CH_HDF5.H"
#include "SPMD.H"
#include "LoadBalance.H"
#include "BoxIterator.H"
#include "AMRIO.H"
#include "AMR.H"
#include "computeSum.H"
#include "computeNorm.H"

#include "EBAMRSpecies.H"
#include "PlasmaPhysics.H"
#include "ParmParse.H"
#include "EBAMRDataOps.H"
#include "EBLevelCCProjector.H" // for averaging
#include "EBAMRSpeciesF_F.H"
#include "EBPatchPolytropic.H"

#include "NamespaceHeader.H"
//EBLevelCrankNicolson is in lib/src/EBAMRElliptic/EBLevelTGA.H
Vector<RefCountedPtr<EBLevelCrankNicolson> >                  EBAMRSpecies::s_diffuseLevTGA;
Vector<RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB> > > >  EBAMRSpecies::s_diffuseAMRMG;
Vector<RefCountedPtr<EBSpeciesFluxOpFactory >  >              EBAMRSpecies::s_diffuseOpFact;


Vector<Vector<RefCountedPtr<LevelData<EBCellFAB> > >  >       EBAMRSpecies::s_aco;
Vector<Vector<RefCountedPtr<LevelData<EBFluxFAB> > >  >       EBAMRSpecies::s_Dk;
Vector<Vector<RefCountedPtr<LevelData<EBFluxFAB> > >  >       EBAMRSpecies::s_vecW;
Vector<Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > > > EBAMRSpecies::s_DkIrreg;
Vector<Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > > > EBAMRSpecies::s_vecWIrreg;

// the data for the flux based boundary conditions of the diffusive solvers
Vector< RefCountedPtr<LevelData<BaseIVFAB<Real> > > > EBAMRSpecies:: s_data;
Vector< RefCountedPtr<LevelData<BaseIVFAB<Real> > > > EBAMRSpecies:: s_rhosRHS;


BiCGStabSolver<LevelData< EBCellFAB> >                EBAMRSpecies::s_botSolver;
Vector<Vector<EBSpeciesFluxOp *>  >                   EBAMRSpecies::s_EBAMROps;

bool EBAMRSpecies::s_isLoadBalanceSet = false;
bool EBAMRSpecies::s_isOnlyPlasma = false;
LoadBalanceFunc EBAMRSpecies::s_loadBalance  = NULL;
int  EBAMRSpecies::s_NewPlotFile = 0;

// we should use dimensional transport coeffs, thus setting beta to -1 for LHS
Real EBAMRSpecies::s_beta = -1;
Real EBAMRSpecies::s_alpha = 0;
/***************************/
// Only and default constructor
EBAMRSpecies::EBAMRSpecies()
{
  m_cfl = 0.8;
  m_useMassRedist = true;
  m_doRZCoords = false;
  m_doRZCoordsCANC = true;
  m_hasSourceTerm = false;
  m_doSmushing = true;
  m_dx = RealVect::Unit;
  m_domainLength = RealVect::Unit;
  m_refineThresh = 5e-11;
  m_initial_dt_multiplier = 0.1;
  m_redistRad = 1;
  m_isDefined = false;
  m_useInject=false;
  m_doImplicitReflux = true;
  m_origin = RealVect::Zero;
  m_aspect= RealVect::Unit;
  m_plot = 0;
  m_ghost_Dk = 2;
  m_steps = 0;
  m_SFD = false;
  m_isPrintResiduals = false;
}
/*******/
void
EBAMRSpecies::
define(const Real&                 a_cfl,
       const RealVect&             a_domainLength,
       const Real&                 a_refineThresh,
       const int&                  a_tagBufferSize,
       const Real&                 a_initialDtMultiplier,
       const bool&                 a_useLimiting)
{
  m_isDefined = true;
  m_cfl = a_cfl;
  m_domainLength = a_domainLength;
  m_refineThresh = a_refineThresh;
  m_tagBufferSize = a_tagBufferSize;
  m_initialDtMultiplier = a_initialDtMultiplier;
  m_useLimiting = a_useLimiting;
  m_hasDiffusion     = true;
  
  // set EB flags hardwired for now
  m_useMassRedist=false;
  m_doSmushing=true;
  m_hasSourceTerm = true;

  /*
  ParmParse ppPlasma;
  m_SFD = ppPlasma.contains("sfd_chi") && ppPlasma.contains("sfd_Dl");
  //m_SFD = false; // cancel
  if(m_SFD)
    {
      ppPlasma.get("sfd_chi_exp", m_chi);
      ppPlasma.get("sfd_Dl", m_Dl);
    }
  */
  m_SFD = false;
  m_chi = 1e7;
  m_chiImp = 0e0;
  m_Dl = 1e-8;
  ParmParse ppPlasma;
  m_tagAllIrreg = -1;//negative number
  ppPlasma.query("tag_all_irreg", m_tagAllIrreg);


  ParmParse ppViscous("ns");
  m_doImplicitReflux = true;
  ppViscous.query("do_implicit_reflux", m_doImplicitReflux);
}
/********/
void EBAMRSpecies::define(AMRLevel*            a_coarserLevelPtr,
			  const ProblemDomain& a_problem_domain,
			  int                  a_level,
			  int                  a_refRatio)
{
  AMRLevel::define(a_coarserLevelPtr, a_problem_domain, a_level, a_refRatio);

  m_problem_domain = a_problem_domain;
  m_domainBox = m_problem_domain.domainBox();

  if (a_coarserLevelPtr != NULL)
    {
      EBAMRSpecies* amrg_ptr = dynamic_cast<EBAMRSpecies*>(a_coarserLevelPtr);

      if (amrg_ptr == NULL)
        {
          pout() << "EBAMRSpecies::define:cannot cast  to derived class"
                 << endl;
          MayDay::Error();
        }

      define(amrg_ptr->m_cfl,
	     amrg_ptr->m_domainLength,
	     amrg_ptr->m_refineThresh,
	     amrg_ptr->m_tagBufferSize,
	     amrg_ptr->m_initialDtMultiplier,
	     amrg_ptr->m_useLimiting);
    }

  for (int idir = 0; idir < SpaceDim; ++idir)
    {
      m_dx[idir] = m_domainLength[idir]/m_domainBox.size(idir);
    }
  m_nGhost = 4;

  //if (m_ebPatchSpecies != NULL)  delete m_ebPatchSpecies;

  m_ebPatchSpecies = m_ebPatchSpeciesFactory->create();
  m_ebPatchSpecies->define(m_problem_domain, m_dx);
  


  m_nComp = m_PlasmaPhysics->nComponents();
  m_stateNames  = Vector<string>(m_nComp, string("scalar"));


  //define chemistry box
  defineRefinementBoxes();
}

/*****************************/
void EBAMRSpecies::defineRefinementBoxes()
{
  // void
}
/********************************/
void EBAMRSpecies::getIVBoxLimits(IntVect& a_loInt, IntVect& a_hiInt, 
				 RealVect a_loReal,RealVect a_hiReal, 
				 const Real& a_loAdd, const Real& a_hiAdd)
{
  a_loReal += a_loAdd;
  a_hiReal += a_hiAdd;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_loInt[idir] = floor(a_loReal[idir]/m_dx[idir] - 0.5);
      a_hiInt[idir] =  ceil(a_hiReal[idir]/m_dx[idir] - 0.5); 
    }
}
/***************************/
void EBAMRSpecies::
setPhysics(PlasmaPhysics*  a_PlasmaPhysics)
{
  m_PlasmaPhysics = a_PlasmaPhysics;
}

/***************************/
void EBAMRSpecies::ZeroVelSetTempPhi(const Vector< LevelData<EBCellFAB>* >& a_Velo,
				     const Vector< LevelData<EBCellFAB>* >& a_Temp,
				     const Vector< LevelData<EBCellFAB>* >& a_phi,
				     const Vector< LevelData<EBFluxFAB>* >& a_gradPhi,
				     const Vector< LevelData<EBCellFAB>* >& a_cellGradPhi,
				     const Vector< LevelData<BaseIVFAB<Real>>* >& a_IVGradPhi)
{
  m_Velo = a_Velo;
  EBLevelDataOps::setToZero(*(m_Velo[m_level]));
  m_Temp = a_Temp;
  m_PlasmaPhysics->setConstantTemperature(*(m_Temp[m_level]), m_stateNew);
  
  m_phi = a_phi;
  m_gradPhi = a_gradPhi;
  m_cellGradPhi = a_cellGradPhi;
  m_IVGradPhi = a_IVGradPhi;
}
/***************************/
void EBAMRSpecies::SetVelTempPhi(const Vector< LevelData<EBCellFAB>* >& a_Velo,
				 const Vector< LevelData<EBCellFAB>* >& a_Temp,
				 const Vector< LevelData<EBCellFAB>* >& a_phi,
				 const Vector< LevelData<EBFluxFAB>* >& a_gradPhi,
				 const Vector< LevelData<EBCellFAB>* >& a_cellGradPhi,
				 const Vector< LevelData<BaseIVFAB<Real>>* >& a_IVGradPhi)
{  
  m_Velo = a_Velo;
  m_Temp = a_Temp;
  m_phi = a_phi;
  m_gradPhi = a_gradPhi;
  m_cellGradPhi = a_cellGradPhi;
  m_IVGradPhi = a_IVGradPhi;
}
/***************************/
void EBAMRSpecies::SetVelTempPhi(const Vector< LevelData<EBCellFAB>* >& a_Velo,
				 const Vector< LevelData<EBCellFAB>* >& a_Temp,
				 const Vector< LevelData<EBCellFAB>* >& a_phi,
				 const Vector< LevelData<EBFluxFAB>* >& a_gradPhi,
				 const Vector< LevelData<EBCellFAB>* >& a_cellGradPhi)
{  
  m_Velo = a_Velo;
  m_Temp = a_Temp;
  m_phi = a_phi;
  m_gradPhi = a_gradPhi;
  m_cellGradPhi = a_cellGradPhi;
}
/***************************/
void EBAMRSpecies::SetVelTempPhi(const Vector< LevelData<EBCellFAB>* >& a_Velo,
				 const Vector< LevelData<EBCellFAB>* >& a_Temp,
				 const Vector< LevelData<EBCellFAB>* >& a_phi,
				 const Vector< LevelData<EBFluxFAB>* >& a_gradPhi)
{  
  m_Velo = a_Velo;
  m_Temp = a_Temp;
  m_phi = a_phi;
  m_gradPhi = a_gradPhi;
}
/***************************/
void EBAMRSpecies::SetVelTempPhi(const Vector< LevelData<EBCellFAB>* >& a_Velo,
				 const Vector< LevelData<EBCellFAB>* >& a_Temp,
				 const Vector< LevelData<EBCellFAB>* >& a_phi)
{  
  m_Velo = a_Velo;
  m_Temp = a_Temp;
  m_phi = a_phi;
}
/*******/
void
EBAMRSpecies::
defineSolvers(const bool& flag)//flag is defaulted to true
{
  
  int numSmooth, numMG, maxIter, verby;
  Real eps, hang; 
  ParmParse ppViscous("ns");
  ppViscous.get("num_smooth", numSmooth);
  ppViscous.get("num_mg",     numMG);
  ppViscous.get("max_iterations", maxIter);
  ppViscous.get("tolerance", eps);
  ppViscous.get("verbosity", verby);
#ifdef CH_USE_FLOAT
  eps = sqrt(eps);
#endif
  ppViscous.get("hang",      hang);
  Real normThresh = 1.0e-15;


  s_botSolver.m_verbosity = 0;

  
  Vector<EBAMRSpecies*>        hierarchy;
  Vector<DisjointBoxLayout>    grids;
  Vector<EBISLayout>           ebisl;
  Vector<EBLevelGrid>          eblg;
  Vector<int>                  refRat;
  ProblemDomain                lev0Dom;
  Real                         lev0Dx;
  Vector<ProblemDomain>        domains;
  Vector<Real>                 dxs;
  getHierarchyAndGrids(hierarchy, grids, ebisl, eblg, refRat, lev0Dom, lev0Dx, domains, dxs);
  
  if(m_level >= m_phi.size()) return;
  //set on all levels
  if(flag) setDiffCoeff(hierarchy, grids, ebisl, eblg, refRat, lev0Dom, dxs, domains);

  //loop over the species (m_nComp)
  s_diffuseLevTGA.resize(m_nComp);
  s_diffuseAMRMG.resize(m_nComp);
  s_diffuseOpFact.resize(m_nComp);
  s_EBAMROps.resize(m_nComp);
  for(int iSpec =0;iSpec<m_nComp; iSpec++)
    {      
      getEBSFOFactory(s_diffuseOpFact[iSpec], grids, ebisl, eblg, refRat, lev0Dx, domains,iSpec);
      { // define the Species Flux Operator
	RefCountedPtr< AMRLevelOpFactory< LevelData< EBCellFAB > > > OpFact = s_diffuseOpFact[iSpec];
	
	s_diffuseAMRMG[iSpec] = RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB> > > (new AMRMultiGrid<LevelData<EBCellFAB> >());
	
	s_diffuseAMRMG[iSpec]->define(lev0Dom, *s_diffuseOpFact[iSpec], &s_botSolver, hierarchy.size());
	s_diffuseAMRMG[iSpec]->setSolverParameters(numSmooth, numSmooth, numSmooth,
						   numMG, maxIter, eps, hang, normThresh);
	s_diffuseAMRMG[iSpec]->m_verbosity = verby-1;
     
	
	s_diffuseLevTGA[iSpec] = RefCountedPtr<EBLevelCrankNicolson>
	  (new EBLevelCrankNicolson(grids, refRat, lev0Dom, OpFact, s_diffuseAMRMG[iSpec]));
	s_diffuseLevTGA[iSpec]->setEBLG(eblg);
	s_diffuseLevTGA[iSpec]->setBeta(s_beta);

	
	Vector<AMRLevelOp <LevelData<EBCellFAB> > *> AMROps = s_diffuseAMRMG[iSpec]->getAMROperators(); 
	s_EBAMROps[iSpec].resize(AMROps.size());
	for (int ilev = 0; ilev < AMROps.size(); ilev++)
	  {
	    s_EBAMROps[iSpec][ilev] =  (static_cast <EBSpeciesFluxOp*> (AMROps[ilev]));
	    s_EBAMROps[iSpec][ilev]->setAlphaAndBeta(s_alpha, s_beta);
	    s_EBAMROps[iSpec][ilev]->setLevel(ilev);
	  }

      }//end ebsfo AMRMG setup
    }
}

/**/
void EBAMRSpecies::
getEBSFOFactory(RefCountedPtr<EBSpeciesFluxOpFactory>&                      a_factory,
		const Vector<DisjointBoxLayout>&                            a_grids,
		const Vector<EBISLayout>&                                   a_ebisl,
		const Vector<EBLevelGrid>&                                  a_eblg,
		Vector<int>&                                                a_refRat,
		Real&                                                       a_lev0Dx,
		Vector<ProblemDomain>&                                      a_domains,
		const int&                                                  a_ispec)
{ 
  RefCountedPtr<BaseDomainBCFactory> domBC; // either zero Neumann or value Dirichlet
  RefCountedPtr<BaseEBBCFactory>     ebBC; //  zero Dirichlet
  getSpeciesBCFactories(domBC, ebBC, a_grids, a_ebisl, a_ispec);

  Vector<RefCountedPtr<EBQuadCFInterp> > quadCFI(a_grids.size());
  int nvar = 1;
  for (int ilev = 0; ilev < a_grids.size(); ilev++)
    {
      if (ilev > 0)
        {
          int nref = a_refRat[ilev-1];
          quadCFI[ilev] = RefCountedPtr<EBQuadCFInterp>(new EBQuadCFInterp(a_grids[ilev],
                                                                           a_grids[ilev-1],
                                                                           a_ebisl[ilev],
                                                                           a_ebisl[ilev-1],
                                                                           a_domains[ilev-1],
                                                                           nref, nvar,
                                                                           *a_eblg[ilev].getCFIVS()));
        }
    }

  ParmParse ppViscous("ns");
  int relaxType = 0; //jacobi
  ppViscous.get("relax_type", relaxType);

  // note the last passed value is the number of component of s_data
  a_factory = RefCountedPtr<EBSpeciesFluxOpFactory>
    (new EBSpeciesFluxOpFactory(a_eblg, quadCFI, s_alpha, s_beta, s_aco[a_ispec], s_Dk[a_ispec], s_vecW[a_ispec], s_DkIrreg[a_ispec], s_vecWIrreg[a_ispec], s_data,
				a_lev0Dx,  a_refRat, domBC, ebBC, m_nGhost*IntVect::Unit, m_nGhost*IntVect::Unit, relaxType, m_nComp+2));

  Vector<int> wallBCInfo = m_PlasmaPhysics->getWallBCinfo();
  //Note:: setData is in EBSpeciesFluxOpFactory.H
  a_factory->setData(s_data,a_ispec, wallBCInfo);
}
/********************/
void EBAMRSpecies::
getSpeciesBCFactories(RefCountedPtr<BaseDomainBCFactory>&         a_domBC,
		      RefCountedPtr<BaseEBBCFactory>&             a_ebBC,
		      const Vector<DisjointBoxLayout>&            a_grids,
		      const Vector<EBISLayout>&                   a_ebisl,
		      const int&                                  a_ispec)
{

  MixedSpeciesFluxDomainBCFactory* domainBCFactory = new MixedSpeciesFluxDomainBCFactory();
  domainBCFactory->setPhysics(m_PlasmaPhysics, a_ispec);
  a_domBC = RefCountedPtr<BaseDomainBCFactory>(domainBCFactory);


  Vector<int> wallBCInfo = m_PlasmaPhysics->getWallBCinfo();
  //NeumannSpeciesFluxEBBCFactory* ebBC = new NeumannSpeciesFluxEBBCFactory();
  MixedSpeciesFluxEBBCFactory* ebBC = new MixedSpeciesFluxEBBCFactory();
  a_ebBC = RefCountedPtr<BaseEBBCFactory>(ebBC);
}

/*******/
Real
EBAMRSpecies::
advance()
{
  if (s_verbosity >= 2) pout() << "EBAMRSpecies::advance " << m_level << endl;
  //These are debugging Hooks to print the AMR solution within substeps
  EBSpeciesFluxOp::s_step = AMR::s_step;
  for(int iSpec=0;iSpec<m_nComp;iSpec++) s_EBAMROps[iSpec][m_level]->increaseStep();

  // Copy the new to the old
  m_stateNew.copyTo(m_stateNew.interval(), m_stateOld, m_stateOld.interval());

  //SFD CFI exchange
  if (m_SFD) {
    Interval consInterv(0, m_nComp-1);
    if (m_hasCoarser)
      {
	const EBAMRSpecies* amrGodCoarserPtr = getCoarserLevel();
	int refRatCrse = amrGodCoarserPtr->refRatio();
	EBPWLFillPatch patcher(m_grids,
			       amrGodCoarserPtr->m_grids,
			       m_ebisl,
			       amrGodCoarserPtr->m_ebisl,
			       amrGodCoarserPtr->m_domainBox,
			       refRatCrse, m_nComp, m_nGhost);
	
	  Real coarTimeOld = 0.0;
	  Real coarTimeNew = 1.0;
	  Real fineTime    = 1.0;
	  patcher.interpolate(m_qbar,
			      amrGodCoarserPtr->m_qbar,
			      amrGodCoarserPtr->m_qbar,
			      coarTimeOld,
			      coarTimeNew,
			      fineTime,
			      consInterv);
      }
    m_qbar.exchange(consInterv);
  }
  
  //set pointers
  //set the advective velocity
  fillAdvectionVelocity();
  m_levelSpecies.setVelTempPtrs(&m_AdvectiveVelocity, m_Velo[m_level], m_Temp[m_level]);

  EBCellFactory factory(m_ebisl);
  LevelData<EBCellFAB> diffusiveSrc(m_grids, m_nComp, m_stateNew.ghostVect(),factory);  
  EBLevelDataOps::setToZero(diffusiveSrc); 
  makeAddReactiveSource(diffusiveSrc, m_dt/2);
  
  
  //add sfd source
  if (m_SFD)
    {
      int ibox = 0;
      Real coec = -m_chi; Real coeb = m_chi;
      for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit, ibox++)
	{
	  const Box& cellBox = m_grids.get(dit());
	  const EBISBox& ebisBox = m_ebisl[dit()];
	  if (!ebisBox.isAllCovered())
	    {
	      
	      EBCellFAB& consState = m_stateNew[dit()];
	      const Box& bigBox = consState.box();
	      
	      EBCellFAB& source = diffusiveSrc[dit()];
	      EBCellFAB& barState = m_qbar[dit()];
	      source.plus(consState,coec);
	      source.plus(barState,coeb);
            }
	}
    }

  //forced (mobility) and molecular diffusion
  Real new_dt = diffusiveAdvance(diffusiveSrc);

  makeAddReactiveSource(diffusiveSrc, m_dt/2);


  Interval interv(0, m_nComp-1);
  // increment redistribution register between this level
  // and next coarser level  
  if (m_hasCoarser)
    {
      for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
	{
	  m_ebFineToCoarRedist.increment(m_massDiff[dit()], dit(), interv);
	}
    }

  //initialize redistirbution register between this level
  // and next finer level
  //this includes re-redistribution registers
  if (m_hasFiner)
    {
      m_ebCoarToFineRedist.setToZero();
      m_ebCoarToCoarRedist.setToZero();
      for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
	{
	  BaseIVFAB<Real>& massDiffFAB = m_massDiff[dit()];
	  m_ebCoarToFineRedist.increment(massDiffFAB, dit(), interv);
	  m_ebCoarToCoarRedist.increment(massDiffFAB, dit(), interv);
	}
    }
  
  //SFD Update
  if (m_SFD)
    { 
      Real acoef = m_Dl/(m_dt+m_Dl);Real bcoef = m_dt/(m_dt+m_Dl); 
      EBLevelDataOps::axby(m_qbar,m_qbar, m_stateNew,  acoef, bcoef);
    }

   //lucacancel temporary fix (same as in EBAMRNavier)
  LevelData<EBCellFAB> consTemp(m_grids, m_nComp, m_nGhost*IntVect::Unit, factory);
  m_ebPatchSpecies->getEBPhysIBC()->initialize(consTemp, m_ebisl);
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      Box gridBox = m_stateNew[dit()].getRegion() & consTemp[dit()].getRegion();
      const EBGraph& ebgraph = m_ebisl[dit()].getEBGraph();
      IntVectSet ivsTot(gridBox);

      for (VoFIterator vofit(ivsTot, ebgraph); vofit.ok(); ++vofit)
	{
	  const VolIndex& vof = vofit();
	  const IntVect& iv = vof.gridIndex();
	  if (iv[0] <=0)
	    {
	      for (int icomp = 0; icomp < m_nComp; icomp++)
		m_stateNew[dit()](vof, icomp) = consTemp[dit()](vof, icomp);
	    }
	}
    }
  extrapolateZeroArea();

  // Update the time and store the new timestep
  m_time += m_dt;
  Real return_dt = m_cfl * new_dt;

  //save stable timestep to save computational effort
  m_dtNew = return_dt;

  return return_dt;
}
/*********/
Real
EBAMRSpecies::
diffusiveAdvance(LevelData<EBCellFAB>& a_source)
{ 

  if (s_verbosity >= 2) pout() << "EBAMRSpecies::diffusiveAdvance " << m_level << endl;
  //set up arguments to step
  //undefined lfr in case we need it
  EBFluxRegister lfr;
  //undefined leveldata in case we need it
  LevelData<EBCellFAB> ld;

  Real tCoarserOld, tCoarserNew, tFinerOld, tFinerNew;

  EBFluxRegister* coarFR = &lfr;
  EBFluxRegister* fineFR = &lfr;
  LevelData<EBCellFAB>* coarDataOld = &ld;
  LevelData<EBCellFAB>* coarDataNew = &ld;

  getCoarseFineDataRegisters(&coarDataOld, &coarDataNew, &coarFR,  &fineFR,
			    tCoarserOld, tCoarserNew, tFinerOld, tFinerNew);


  // Real new_dt = m_dt;
  // Advance the solve one timestep (note m_dt is not changedin step, thus we can use it in implicitDiffusion)
  Real new_dt = 1.0;
  if(!s_isOnlyPlasma) {
    new_dt = m_levelSpecies.step(m_stateNew,
				 m_massDiff,
				 a_source,
				 *fineFR,
				 *coarFR,
				 *coarDataOld,
				 *coarDataNew,
				 m_time,
				 tCoarserOld,
				 tCoarserNew,
				 m_dt);
  }
  

  if (m_hasDiffusion)
    {
      implicitDiffusion(fineFR,coarFR,*coarDataOld,*coarDataNew,tCoarserOld,tCoarserNew,a_source);
    }


  return new_dt;
}

/*********/
void
EBAMRSpecies::
makeDiffusiveSource(LevelData<EBCellFAB>& a_diffusiveSrc)
{
  if (s_verbosity >= 2) pout() << "EBAMRSpecies::makeDiffusiveSource" << m_level << endl;
  //assert the opertors are on the LHS
  CH_assert(s_beta <= 0);
    
  if (m_hasDiffusion)
    {

      //undefined leveldata in case we need it
      LevelData<EBCellFAB> ld;
      LevelData<EBCellFAB>* coarserDataOldPtr = &ld;
      LevelData<EBCellFAB>* coarserDataNewPtr = &ld;
      LevelData<EBCellFAB>* finerDataOldPtr = &ld;
      LevelData<EBCellFAB>* finerDataNewPtr = &ld;
      Real tCoarserOld, tCoarserNew, tFinerOld, tFinerNew;
      getCoarseFineDataPointers(&coarserDataOldPtr,&coarserDataNewPtr,
                                &finerDataOldPtr,&finerDataNewPtr,
                                tCoarserOld, tCoarserNew, tFinerOld, tFinerNew);
      LevelData<EBCellFAB>* stateUCoarse = new LevelData<EBCellFAB>();
      LevelData<EBCellFAB>* stateUFine = new LevelData<EBCellFAB>();

      EBCellFactory factory(m_ebisl);
      LevelData<EBCellFAB>& stateU  = m_stateOld;
      LevelData<EBCellFAB> lhs(m_grids, 1, m_stateOld.ghostVect(), factory);


      Real interpTime = m_time + 0.5*m_dt;
      if (m_hasCoarser)
	{
	  EBAMRSpecies* coarserPtr = getCoarserLevel();
	  EBCellFactory Cfactory(coarserPtr->m_ebisl);
	  stateUCoarse->define(coarserDataOldPtr->disjointBoxLayout(),coarserDataOldPtr->nComp(), coarserDataOldPtr->ghostVect(), Cfactory);
	  interpolateInTime(*stateUCoarse, *coarserDataOldPtr, *coarserDataNewPtr,coarserPtr->m_grids,interpTime, tCoarserOld, tCoarserNew);    
	  EBLevelDataOps::setCoveredVal(*stateUCoarse,0e0);  
	}
      if (m_hasFiner)
	{
	  EBAMRSpecies* finerPtr = getFinerLevel();
	  EBCellFactory Ffactory(finerPtr->m_ebisl);
	  stateUFine->define(finerDataOldPtr->disjointBoxLayout(),finerDataOldPtr->nComp(), finerDataOldPtr->ghostVect(), Ffactory);
	  interpolateInTime(*stateUFine, *finerDataOldPtr, *finerDataNewPtr,finerPtr->m_grids,interpTime, tFinerOld, tFinerNew);
	  EBLevelDataOps::setCoveredVal(*stateUFine,0e0);
	  
	}


      for(int iSpec=0;iSpec<m_nComp;iSpec++)
	{
	  Interval allInterv(iSpec,iSpec);
	  Interval thisInterv(0, 0);
	  // uisng the letter V in analogy with velocity solver
	  LevelData<EBCellFAB> stateV;
	  LevelData<EBCellFAB> stateVCoarse; 
	  LevelData<EBCellFAB> stateVFine;
	  //aliasing to solve species by species
	  /*
	  aliasLevelData<EBCellFAB>(stateV, &stateU, allInterv);
	  aliasLevelData<EBCellFAB>(stateVCoarse, stateUCoarse, allInterv);
	  aliasLevelData<EBCellFAB>(stateVFine,   stateUFine,   allInterv);
	  EBLevelDataOps::setCoveredVal(stateV,0e0);*/
	  EBcopyLevelData(stateV, stateU, allInterv, m_ebisl);
	  if (m_hasCoarser)
	    {
	      EBAMRSpecies* coarserPtr = getCoarserLevel();
	      EBcopyLevelData(stateVCoarse, *stateUCoarse, allInterv, coarserPtr->m_ebisl);
	    }
	  if (m_hasFiner)
	    {
	      EBAMRSpecies* finerPtr = getFinerLevel();
	      EBcopyLevelData(stateVFine,   *stateUFine,   allInterv, finerPtr->m_ebisl);
	    }
      
	  if (!m_hasCoarser && !m_hasFiner)
	    {
	      s_EBAMROps[iSpec][m_level]->applyOp (lhs, stateV, false);
	    }
	  else if (!m_hasCoarser)
	    {
	      s_EBAMROps[iSpec][m_level]->AMROperatorNC (lhs, stateVFine, stateV, false, s_EBAMROps[iSpec][m_level+1]);
	    }
	  else if(!m_hasFiner)
	    {
	      s_EBAMROps[iSpec][m_level]->AMROperatorNF (lhs, stateV, stateVCoarse, false);
	    }
	  else
	    {
	      s_EBAMROps[iSpec][m_level]->AMROperator (lhs, stateVFine, stateV, stateVCoarse, false, s_EBAMROps[iSpec][m_level+1]);
	    }

	  EBLevelDataOps::assign(a_diffusiveSrc,lhs, allInterv,thisInterv);
	}

      //divide by the volume fraction
      EBLevelDataOps::kappaDivide(a_diffusiveSrc);
      // bring on the RHS
      EBLevelDataOps::scale(a_diffusiveSrc,-1.0);

      a_diffusiveSrc.exchange();

      delete stateUCoarse; 
      delete stateUFine;

      
      /*pout() << "NAN Check diffusiveSrc " << m_level << endl;
      bool isNAN = EBLevelDataOps::checkNANINF(a_diffusiveSrc);
      if(isNAN) MayDay::Error("Checking diffusiveSrc.");
      pout() << "ByeBye";exit(0);*/
      
    }

}
/***********/
void 
EBAMRSpecies::implicitDiffusion(EBFluxRegister*& finerFRPtr, EBFluxRegister*& coarserFRPtr, 
			       LevelData<EBCellFAB>& oldUCoarse, LevelData<EBCellFAB>& stateUCoarse,
			       Real& tCoarserOld, Real& tCoarserNew, LevelData<EBCellFAB>& a_src)
{
  if (s_verbosity >= 2) pout() << "EBAMRSpecies::implicitDiffusion" << m_level << endl;
  // Advance Both Natural and Forced Diffusion Terms
  m_steps++;
  EBCellFactory factory(m_ebisl);
  LevelData<EBCellFAB> zeroSrc(m_grids, 1, m_stateNew.ghostVect(),factory); 
  EBLevelDataOps::setToZero(zeroSrc);
  EBFluxRegister* NfinerFRPtr=NULL; EBFluxRegister* NcoarserFRPtr=NULL;

  bool constTe = m_PlasmaPhysics->isConstantTe();
  bool constBg = m_PlasmaPhysics->isConstantBg();
  int maxSpec = constTe ? (m_nComp-1) : m_nComp;
  int minSpec = constBg ? (1) : 0;

  int k=m_PlasmaPhysics->findSpec(string("e-"));
  Real maxElec, minElec; EBLevelDataOps::getMaxMin(maxElec, minElec, m_stateNew,k);
  Real DTsign = 1e0;

  for(int iSpec=minSpec; iSpec<maxSpec; iSpec++)
    {
      //pout() << "Species " << iSpec << endl;
      Interval allInterv(iSpec,iSpec);
      Interval thisInterv(0, 0);
      //aliasing to solve species by species
      LevelData<EBCellFAB> state;
      LevelData<EBCellFAB> initial;
      LevelData<EBCellFAB> stateOld;
      LevelData<EBCellFAB> coarseState;
      LevelData<EBCellFAB> coarseStateOld;
      //LevelData<EBCellFAB> src;
      /*
      aliasLevelData<EBCellFAB>(state, &m_stateNew, allInterv);
      aliasLevelData<EBCellFAB>(stateOld, &m_stateOld, allInterv);
      aliasLevelData<EBCellFAB>(coarseState, &stateUCoarse, allInterv);
      aliasLevelData<EBCellFAB>(coarseStateOld, &oldUCoarse, allInterv); */
      EBcopyLevelData(state,    m_stateNew, allInterv, m_ebisl);
      EBcopyLevelData(stateOld, m_stateOld, allInterv, m_ebisl);
      EBcopyLevelData(initial,  m_stateNew, allInterv, m_ebisl);
      //EBcopyLevelData(src,  a_src, allInterv, m_ebisl);
      if (m_hasCoarser)
	{
	  EBAMRSpecies* coarserPtr = getCoarserLevel();
	  EBcopyLevelData(coarseState, stateUCoarse, allInterv, coarserPtr->m_ebisl);
	  EBcopyLevelData(coarseStateOld, oldUCoarse, allInterv, coarserPtr->m_ebisl);
	}
      //updateSoln is in ~/ChomboVT/lib/src/AMRElliptic/BaseLevelCrankNicolson.H
      s_diffuseLevTGA[iSpec]->updateSoln(state, initial, zeroSrc,
					 finerFRPtr, coarserFRPtr,
					 &coarseStateOld, &coarseState,
					 m_time, tCoarserOld, tCoarserNew,
					 DTsign*m_dt, m_level, false, false, iSpec);// DTsign >0 CN, DTsign <0 BWDEul
      EBcopyLevelData(m_stateNew, allInterv,  state,  thisInterv);
      int exitStatus = s_diffuseLevTGA[iSpec]->getExitStatus();
      if(exitStatus == 4 || exitStatus == 6)
	{
	  Real residNorm = s_diffuseLevTGA[iSpec]->getResidNorm();
	  if(residNorm > 1e-10)
	    {
	      pout() <<  m_steps << ", " << m_level <<   ", " << iSpec << " Failed Integration " << exitStatus << " residNorm " << residNorm << endl;
	      PlotReactiveSource(10000 + m_steps);
	      exit(0);
	    }
	}
      if(m_isPrintResiduals)
	{
	  LevelData<EBCellFAB> resid;
	  EBcopyLevelData(resid,    m_stateNew, allInterv, m_ebisl);
	  s_diffuseLevTGA[iSpec]->getResidual(resid, state, initial, zeroSrc,
					      &coarseStateOld, &coarseState,
					      m_time, tCoarserOld, tCoarserNew,
					      m_dt, m_level, false);
	  Real maxElec, minElec; EBLevelDataOps::getMaxMin(maxElec, minElec, resid,0);
	  pout() << "residual" << maxElec << ", " <<  minElec << endl;
	}
      // Uncomment the line below to plot the residual using PlotReactiveSource
      /* EBcopyLevelData(m_stateNew, allInterv,  resid,  thisInterv);*/
      
    }
  //PlotReactiveSource(1000+m_doRZCoordsCANC);
  //exit(0);
  //pout() << "NAN Check diffusiveSrc " << m_level << endl;
  /*bool isNAN = EBLevelDataOps::checkNANINF(m_stateNew);
  int tmp = isNAN;int allNAN;
  MPI_Allreduce(&tmp, &allNAN, 1, MPI_INT, MPI_SUM, Chombo_MPI::comm);
  if(allNAN>0) { pout() << "NAN in level " << m_level << endl;PlotReactiveSource(1000); MPI_Barrier(Chombo_MPI::comm);MayDay::Error("Checking m_stateNew after implicitDiffusion.");}
  //pout() << "ByeBye";exit(0);*/
}
/*********/
void
EBAMRSpecies::
makeAddReactiveSource(LevelData<EBCellFAB>& a_diffusiveSrc, const Real& a_dt)
{

  if (s_verbosity >= 2) pout() << "EBAMRSpecies::makeAddReactiveSource" << m_level << endl;
  if(a_dt <=0) return;

  // a_diffusiveSrc is zero coming in
  EBCellFactory factory(m_ebisl);
  //evaluate the cell grad ne
  LevelData<EBCellFAB> gradNe(m_grids, SpaceDim, m_stateNew.ghostVect(),factory);
  cellGrad(gradNe, m_stateNew, m_nComp-2);

  //this call zeros out the source terms where Tele > threshold
  LevelData<EBCellFAB> reactiveSrc(m_grids, m_nComp, m_stateNew.ghostVect(),factory);
  LevelData<EBCellFAB> eJHeat(m_grids, m_nComp, m_stateNew.ghostVect(),factory);
  EBLevelDataOps::setToZero(eJHeat);
  m_PlasmaPhysics->addJouleHeat(eJHeat, m_stateNew, *m_Temp[m_level], *m_cellGradPhi[m_level],*m_Velo[m_level], gradNe, m_grids, m_ebisl);
  LevelData<EBCellFAB>& mask=reactiveSrc;
  EBLevelDataOps::setToZero(mask);
  m_PlasmaPhysics->setMask(mask,a_diffusiveSrc,m_stateNew,m_grids,m_ebisl);
  m_PlasmaPhysics->plasmaSourcesODE(m_stateNew, *m_Temp[m_level], eJHeat, *m_cellGradPhi[m_level], m_ebisl,a_dt, m_level);
  

}
/*********/
void
EBAMRSpecies::
makeFluidSource(LevelData<EBCellFAB>& a_fluidSrc)
{
  if (s_verbosity >= 2) pout() << "EBAMRSpecies::makeFluidSource" << m_level << endl;

  EBLevelDataOps::setToZero(a_fluidSrc);
  m_PlasmaPhysics->fluidSources(a_fluidSrc,m_stateNew,*m_Temp[m_level],*m_cellGradPhi[m_level],*m_Velo[m_level], m_ebisl);
  
    //add the joule heating contribution
  EBCellFactory factory(m_ebisl);
  //evaluate the cell grad tensor
  int gComp = m_nComp-2;
  LevelData<EBCellFAB> gradN(m_grids, SpaceDim*gComp, m_stateNew.ghostVect(),factory);
  cellGradTens(gradN, m_stateNew, gComp);
  m_PlasmaPhysics->AddJouleHeatF(a_fluidSrc, m_stateNew, *m_Temp[m_level], *m_cellGradPhi[m_level], gradN, m_grids, m_ebisl);
  
}

/*********/
bool
EBAMRSpecies::
PlotReactiveSource(const int& a_nplot, const Real&  a_dt, const bool&  a_onlyfluid) 
{

  if (s_verbosity >= 2) pout() << "EBAMRSpecies::PlotReactiveSource" << m_level << endl;
      
  Vector<EBAMRSpecies*>        hierarchy;
  Vector<DisjointBoxLayout>    grids;
  Vector<EBISLayout>           ebisl;
  Vector<EBLevelGrid>          eblg;
  Vector<int>                  refRat;
  ProblemDomain                lev0Dom;
  Real                         lev0Dx;
  Vector<ProblemDomain>        domains;
  Vector<Real>                 dxs;
  getHierarchyAndGrids(hierarchy, grids, ebisl, eblg, refRat, lev0Dom, lev0Dx, domains, dxs);
  int numLevels = hierarchy.size();


  EBPatchPolytropic Ppatch;
  Interval momIntrv = Ppatch.momentumInterval();  
  int numC = Ppatch.numConserved();
  ParmParse pp;
  Real SrcV=1.0;
  // pp.get("pinf",SrcV);//uncomment to scale the source terms
    
  m_plot = a_nplot;
    
  int nreacr = m_PlasmaPhysics->nReactions();
  int nPlasmaCoeffs = m_PlasmaPhysics->nComponents(); 
  IntVect ivGhost = m_nGhost*IntVect::Unit;
  Vector<LevelData<EBCellFAB>* > srcdata(numLevels,NULL);
  Vector<LevelData<EBCellFAB>* > rhsdata(numLevels,NULL);
  Vector<LevelData<EBCellFAB>* > dtdata(numLevels,NULL);
  Vector<LevelData<EBCellFAB>* > fsrcdata(numLevels,NULL);
  Vector<LevelData<EBCellFAB>* > reacdata(numLevels,NULL);
  Vector<LevelData<EBCellFAB>* > state(numLevels);
  Vector<LevelData<EBCellFAB>* > stateOld(numLevels);
  Vector<LevelData<EBCellFAB>* > diffdata(numLevels,NULL);
  Vector<LevelData<EBCellFAB>* > mobdata(numLevels,NULL);
  Vector<LevelData<EBCellFAB>* > tdata(numLevels,NULL);
  Vector<LevelData<BaseIVFAB<Real> > * > diffIV(numLevels,NULL);
  Vector<LevelData<BaseIVFAB<Real> > * > mobIV(numLevels,NULL);
  Vector<LevelData<BaseIVFAB<Real> > * > wallIV(numLevels,NULL);
  for (int ilev = 0; ilev < numLevels; ++ilev)
    {
      EBCellFactory factory(ebisl[ilev]);
      srcdata[ilev] = new LevelData<EBCellFAB>();
      srcdata[ilev]->define(grids[ilev], m_nComp+3, ivGhost, factory); 
      rhsdata[ilev] = new LevelData<EBCellFAB>(grids[ilev], 1, ivGhost, factory);
      dtdata[ilev] = new LevelData<EBCellFAB>(grids[ilev], 1, ivGhost, factory);
      fsrcdata[ilev] = new LevelData<EBCellFAB>();
      fsrcdata[ilev]->define(grids[ilev], numC, ivGhost, factory);
      reacdata[ilev] = new LevelData<EBCellFAB>();
      reacdata[ilev]->define(grids[ilev], nreacr, ivGhost, factory);
      LevelData<EBCellFAB>& stateU = hierarchy[ilev]->m_stateNew;
      state[ilev] = &stateU;
      stateOld[ilev] = &hierarchy[ilev]->m_stateOld;
      LevelData<EBCellFAB> eJHeat(grids[ilev], m_nComp, ivGhost,factory);
      EBLevelDataOps::setToZero(eJHeat);

      //evaluate the cell grad ne
      int ie = m_PlasmaPhysics->findSpec(string("e-"));
      int iee = m_PlasmaPhysics->findSpec(string("ee"));
      LevelData<EBCellFAB> gradNe(grids[ilev], SpaceDim, stateU.ghostVect(),factory);
      hierarchy[ilev]->cellGrad(gradNe, stateU, ie);
      m_PlasmaPhysics->addJouleHeat(eJHeat, stateU, *m_Temp[ilev], *m_cellGradPhi[ilev], gradNe, *m_Velo[ilev], grids[ilev], ebisl[ilev]);
      // product of Joule/ne * ne:
      EBLevelDataOps::product(eJHeat,eJHeat,stateU,0,0,ie);
      EBLevelDataOps::setToZero(*srcdata[ilev]);
      m_PlasmaPhysics->plasmaSources(*srcdata[ilev], stateU, *m_Temp[ilev], *m_cellGradPhi[ilev], eJHeat, ebisl[ilev]);
      EBLevelDataOps::setToZero(*rhsdata[ilev]);
      m_PlasmaPhysics->ElectricPotentialRHS(*rhsdata[ilev], hierarchy[ilev]->m_stateNew, dxs[ilev]);
      EBLevelDataOps::setToZero(*dtdata[ilev]);
      m_PlasmaPhysics->maxDT(*dtdata[ilev], hierarchy[ilev]->m_stateNew, *m_Temp[ilev], a_dt);
      EBLevelDataOps::setToZero(*fsrcdata[ilev]);
      hierarchy[ilev]->makeFluidSource(*fsrcdata[ilev]);
      //hierarchy[ilev]->makefluxGphi(*srcdata[ilev]);
      m_PlasmaPhysics->reactionRates(*reacdata[ilev], stateU, *m_Temp[ilev],*m_cellGradPhi[ilev]);
      diffdata[ilev] = new LevelData<EBCellFAB>();
      diffdata[ilev]->define(grids[ilev], nPlasmaCoeffs, ivGhost, factory);
      // this version prints only the electron mobility
      mobdata[ilev] = new LevelData<EBCellFAB>();
      mobdata[ilev]->define(grids[ilev], SpaceDim, ivGhost, factory);
      LevelData<EBCellFAB> tbuff(grids[ilev], nPlasmaCoeffs, ivGhost, factory);
      // fill diffusion pointers
      RealVect dxv=dxs[ilev]*RealVect::Unit;
      for(int k=0; k<nPlasmaCoeffs; k++)
	{
	  Interval MDinterv(k,k);
	  Interval LMinterv(0,0);
	  ccpAverageFaceToCells(tbuff,*s_Dk [k][ilev],grids[ilev],ebisl[ilev],domains[ilev],dxv);
	  EBcopyLevelData(*diffdata[ilev],MDinterv,tbuff,LMinterv);
	}
      //int k=nPlasmaCoeffs-2;// for electrons
      int k=m_PlasmaPhysics->findSpec(string("e-"));
      ccpAverageFaceToCells(*mobdata[ilev],*s_vecW [k][ilev],grids[ilev],ebisl[ilev],domains[ilev],dxv);
      // EBLevelDataOps::incr(*mobdata[ilev],*m_Velo[ilev],-1e0);


      tdata[ilev] = new LevelData<EBCellFAB>();
      tdata[ilev]->define(grids[ilev], 1, ivGhost, factory); 
      // MD Note that only inside cells calculated below
      m_PlasmaPhysics->electronTemperature(*tdata[ilev], stateU, grids[ilev], ebisl[ilev]);
      
      LayoutData<IntVectSet> irregSets(grids[ilev]);
      
      for (DataIterator dit = grids[ilev].dataIterator(); dit.ok(); ++dit) 
	{
	  Box grownBox = grow(grids[ilev].get(dit()), m_ghost_Dk);
	  grownBox &= domains[ilev];
	  irregSets[dit()] = ebisl[ilev][dit()].getIrregIVS(grownBox);
	}
      BaseIVFactory<Real>  baseivfact(ebisl[ilev], irregSets);
      diffIV[ilev]    = new LevelData< BaseIVFAB<Real> >(grids[ilev], nPlasmaCoeffs, m_ghost_Dk*IntVect::Unit, baseivfact);
      mobIV[ilev]    = new LevelData< BaseIVFAB<Real> >(grids[ilev], nPlasmaCoeffs, m_ghost_Dk*IntVect::Unit, baseivfact);
      for(int k=0; k<nPlasmaCoeffs; k++)
	{
	  Interval MDinterv(k,k);
	  Interval LMinterv(0,0);
	  s_DkIrreg[k][ilev]->copyTo(LMinterv, *diffIV[ilev]  , MDinterv);
	  s_vecWIrreg[k][ilev]->copyTo(LMinterv, *mobIV[ilev] , MDinterv);
	}
      wallIV[ilev]    = &*s_data[ilev];	
      for (int k = momIntrv.begin(); k <= momIntrv.end(); k++)
	EBLevelDataOps::scale(*fsrcdata[ilev],1.0/SrcV,k);
    }
  //EBAMRDataOps::checkNANINF(srcdata);
 
  if(!a_onlyfluid)
    {
      outputHDF(m_cellGradPhi,string("cellGradPhi"));
      outputHDF(state,string("State"));
    }
  //
  //    outputHDF(reacdata,string("Rterms"));
  //    outputHDF(srcdata,string("Sterms"));
  if(a_nplot >= 100000)
    {
      outputHDF(fsrcdata,string("FSRC"));
      outputHDF(m_phi,string("phi"));
      outputHDF(m_Velo,  string("Vel"));
      outputHDF(diffdata,   string("Diff"));
      outputHDF(mobdata,   string("Mob"));
      outputHDF(stateOld,string("StateOld"));
      outputHDF(rhsdata,string("RHS"));
      outputHDF(dtdata,string("DT"));
      outputHDF(tdata,   string("Tele"));
      outputHDF(m_Temp,  string("Temp"));
    }
//    outputIVFAB(diffIV,   string("Diff"));
//    outputIVFAB(mobIV,   string("Mob"));
  // outputIVFAB(wallIV,   string("wall"));
  for (int ilev = 0; ilev < numLevels; ++ilev)
      {
	delete srcdata[ilev];
	delete rhsdata[ilev];
	delete dtdata[ilev];
	delete fsrcdata[ilev];
	delete reacdata[ilev];
	delete diffdata[ilev];
	delete mobdata[ilev];
	delete tdata[ilev];
	delete diffIV[ilev];
	delete mobIV[ilev];
      }
	
    return (true);
}

/*********/
void
EBAMRSpecies::
junkDiffusion()
{

  Vector<EBAMRSpecies*>        hierarchy;
  Vector<DisjointBoxLayout>    grids;
  Vector<EBISLayout>           ebisl;
  Vector<EBLevelGrid>          eblg;
  Vector<int>                  refRat;
  ProblemDomain                lev0Dom;
  Real                         lev0Dx;
  Vector<ProblemDomain>        domains;
  Vector<Real>                 dxs;
  getHierarchyAndGrids(hierarchy, grids, ebisl, eblg, refRat, lev0Dom, lev0Dx, domains, dxs);
  int numLevels = hierarchy.size();

  int nreacr = m_PlasmaPhysics->nReactions();
  int nPlasmaCoeffs = m_PlasmaPhysics->nComponents();
  IntVect ivGhost = m_nGhost*IntVect::Unit;
  Vector<LevelData<EBCellFAB>* > srcdata(numLevels,NULL);
  Vector<LevelData<EBCellFAB>* > reacdata(numLevels,NULL);
  Vector<LevelData<EBCellFAB>* > state(numLevels);
  Vector<LevelData<EBCellFAB>* > diffdata(numLevels,NULL);
  Vector<LevelData<EBCellFAB>* > mobdata(numLevels,NULL);
  Vector<LevelData<BaseIVFAB<Real> > * > diffIV(numLevels,NULL);
  Vector<LevelData<BaseIVFAB<Real> > * > mobIV(numLevels,NULL);
  Vector<LevelData<BaseIVFAB<Real> > * > wallIV(numLevels,NULL);
  for (int ilev = 0; ilev < numLevels; ++ilev)
    {
      RealVect dxv=dxs[ilev]*RealVect::Unit;
      EBCellFactory factory(ebisl[ilev]);
      LevelData<EBCellFAB> tbuff(grids[ilev], SpaceDim, ivGhost, factory);
      EBLevelDataOps::setToZero(tbuff);
      EBLevelDataOps::setVal(tbuff,-1e2,1);
      for(int k=0; k<nPlasmaCoeffs; k++)
	{
	  EBLevelDataOps::averageCellToFacesMAC(*s_vecW [k][ilev], tbuff, grids[ilev], ebisl[ilev], domains[ilev]);
	  ccpExtrapolateToDomainBoundaries(*s_vecW [k][ilev],grids[ilev], ebisl[ilev], domains[ilev], dxv);
	  EBLevelDataOps::setVal(*s_Dk [k][ilev],0e0);
	  for (DataIterator dit = grids[ilev].dataIterator(); dit.ok(); ++dit)
	    {
	      (  *s_DkIrreg[k][ilev])[dit()].setVal(0e0);
	      (*s_vecWIrreg[k][ilev])[dit()].setVal(0e2);
	      (        *s_data[ilev])[dit()].setVal(0e2);
	      (        *s_data[ilev])[dit()].setVal(nPlasmaCoeffs+1,0e2);
	    }
	}
    }
	
	
}
/*********/
void
EBAMRSpecies::
makefluxGphi(LevelData<EBCellFAB>& a_fluxGphi)
{
  //assert the opertors are on the LHS
  CH_assert(s_beta <= 0);
    

  //undefined leveldata in case we need it
  LevelData<EBCellFAB> ld;
  LevelData<EBCellFAB>* coarserDataOldPtr = &ld;
  LevelData<EBCellFAB>* coarserDataNewPtr = &ld;
  LevelData<EBCellFAB>* finerDataOldPtr = &ld;
  LevelData<EBCellFAB>* finerDataNewPtr = &ld;
  Real tCoarserOld, tCoarserNew, tFinerOld, tFinerNew;
  getCoarseFineDataPointers(&coarserDataOldPtr,&coarserDataNewPtr,
			    &finerDataOldPtr,&finerDataNewPtr,
			    tCoarserOld, tCoarserNew, tFinerOld, tFinerNew);

  LevelData<EBCellFAB>* stateUCoarse = new LevelData<EBCellFAB>();
  
  EBCellFactory factory(m_ebisl);
  LevelData<EBCellFAB>& stateU  = m_stateNew;


  Real interpTime = m_time + 0.5*m_dt;
  if (m_hasCoarser)
    {
      EBAMRSpecies* coarserPtr = getCoarserLevel();
      EBCellFactory Cfactory(coarserPtr->m_ebisl);
      stateUCoarse->define(coarserDataOldPtr->disjointBoxLayout(),coarserDataOldPtr->nComp(), coarserDataOldPtr->ghostVect(), Cfactory);
      interpolateInTime(*stateUCoarse, *coarserDataOldPtr, *coarserDataNewPtr,coarserPtr->m_grids,interpTime, tCoarserOld, tCoarserNew);
      EBLevelDataOps::setCoveredVal(*stateUCoarse,0e0);
      int ncomp = m_nComp;
      LayoutData<IntVectSet> cfivs;
      EBArith::defineCFIVS(cfivs, m_grids, m_problem_domain);
      EBQuadCFInterp quadCFIWithFine(m_grids,
				     coarserDataOldPtr->disjointBoxLayout(),
				     m_ebisl,
				     coarserPtr->m_ebisl,
				     coarserPtr->m_problem_domain,
				     m_ref_ratio,
				     ncomp, cfivs);

      Interval  intval(0,ncomp-1);
      quadCFIWithFine.interpolate(stateU, *stateUCoarse, intval);
      stateU.exchange(intval);
    }
  
  
  for(int iSpec=0;iSpec<m_nComp;iSpec++)
    {
      s_EBAMROps[iSpec][m_level]->getFluxGradPhi(a_fluxGphi,stateU,*m_gradPhi[m_level],m_AdvectiveVelocity,iSpec);      
    }

  delete stateUCoarse;

}
void
EBAMRSpecies::
getHierarchyAndGrids(Vector<EBAMRSpecies*>&        a_hierarchy,
                     Vector<DisjointBoxLayout>&   a_grids,
                     Vector<EBISLayout>&          a_ebisl,
		     Vector<EBLevelGrid>&         a_eblg,
                     Vector<int>&                 a_refRat,
                     ProblemDomain&               a_lev0Dom,
                     Real&                        a_lev0Dx,
		     Vector<ProblemDomain>&       a_domains,
                     Vector<Real>&                a_dxs)
{
  Vector<AMRLevel*> hierarchy = AMRLevel::getAMRLevelHierarchy();
  int nlevels = hierarchy.size();

  a_hierarchy.resize(nlevels);
  a_refRat.resize(   nlevels);
  a_grids.resize(    nlevels);
  a_ebisl.resize(    nlevels);
  a_eblg.resize(     nlevels);
  a_domains.resize(  nlevels);
  a_dxs.resize    (  nlevels);

  EBAMRSpecies* coarsestLevel = (EBAMRSpecies*)(hierarchy[0]);
  a_lev0Dx       = coarsestLevel->m_dx[0];
  a_lev0Dom      = coarsestLevel->m_problem_domain;


  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();

  for (int ilev = 0; ilev < nlevels; ilev++)
    {
      EBAMRSpecies* adLevel = (EBAMRSpecies*)(hierarchy[ilev]);

      a_hierarchy[ilev] = adLevel;
      a_grids [ilev] = adLevel->m_grids;
      a_refRat[ilev] = adLevel->m_ref_ratio;
      const ProblemDomain& domainLev = adLevel->problemDomain();
      ebisPtr->fillEBISLayout(a_ebisl[ilev], a_grids[ilev], domainLev, m_nGhost);
      a_eblg[ilev] = EBLevelGrid(a_grids[ilev], a_ebisl[ilev], domainLev);
      a_domains[ilev] = domainLev;
      a_dxs[ilev] = adLevel->m_dx[0];
    }
}
/********/
void
EBAMRSpecies::
getCoarseFineDataPointers(LevelData<EBCellFAB>** a_coarserDataOldPtr, LevelData<EBCellFAB>** a_coarserDataNewPtr,
			  LevelData<EBCellFAB>** a_finerDataOldPtr, LevelData<EBCellFAB>** a_finerDataNewPtr,
			  Real& a_tCoarserOld,  Real& a_tCoarserNew, Real& a_tFinerOld,  Real& a_tFinerNew)
{
  *a_coarserDataOldPtr = NULL;
  *a_coarserDataNewPtr = NULL;  
  *a_finerDataOldPtr = NULL;
  *a_finerDataNewPtr = NULL;

  a_tCoarserOld = 0.0;
  a_tCoarserNew = 0.0;
  a_tFinerOld = 0.0;
  a_tFinerNew = 0.0;

  // A coarser level exists
  if (m_hasCoarser)
    {
      EBAMRSpecies* coarserPtr = getCoarserLevel();

      *a_coarserDataOldPtr = &coarserPtr->m_stateOld;
      *a_coarserDataNewPtr = &coarserPtr->m_stateNew;

      a_tCoarserNew = Max(coarserPtr->m_time,0e0);
      a_tCoarserOld = Max(a_tCoarserNew - coarserPtr->m_dt,0e0);
    }

  // A finer level exists
  if (m_hasFiner)
    {
      EBAMRSpecies* finerPtr = getFinerLevel();

      *a_finerDataOldPtr = &finerPtr->m_stateOld;
      *a_finerDataNewPtr = &finerPtr->m_stateNew;

      a_tFinerNew = Max(finerPtr->m_time,0e0);
      a_tFinerOld = Max(a_tFinerNew - finerPtr->m_dt,0e0);
    }
}
/********/
void
EBAMRSpecies::
getCoarseFineDataRegisters(LevelData<EBCellFAB>** a_coarserDataOldPtr, LevelData<EBCellFAB>** a_coarserDataNewPtr,
			  EBFluxRegister** a_coarFR, EBFluxRegister** a_fineFR,
			  Real& a_tCoarserOld,  Real& a_tCoarserNew, Real& a_tFinerOld,  Real& a_tFinerNew)
{

  a_tCoarserOld = 0.0;
  a_tCoarserNew = 0.0;
  a_tFinerOld = 0.0;
  a_tFinerNew = 0.0;

  // A coarser level exists
  if (m_hasCoarser)
    {
      EBAMRSpecies* coarserPtr = getCoarserLevel();

      *a_coarserDataOldPtr = &coarserPtr->m_stateOld;
      *a_coarserDataNewPtr = &coarserPtr->m_stateNew;

      *a_coarFR = &coarserPtr->m_ebFluxRegister;

      a_tCoarserNew = Max(coarserPtr->m_time,0e0);
      a_tCoarserOld = Max(a_tCoarserNew - coarserPtr->m_dt,0e0);
    }

  // A finer level exists
  if (m_hasFiner)
    {      
      *a_fineFR = &m_ebFluxRegister;
    }
}
/********/
LevelData<EBCellFAB>*
EBAMRSpecies::getStatePtr()
{
  return (&m_stateNew);
}
Vector< RefCountedPtr<LevelData<BaseIVFAB<Real> > > >
EBAMRSpecies::getRhosRHS()
{
  return (s_rhosRHS);
}
/********/
LevelData<EBCellFAB>*
EBAMRSpecies::getBarPtr()
{
  if(m_SFD)
    return (&m_qbar);
  else
    return (&m_stateNew);
}
/********************************/
void
EBAMRSpecies::
interpolateInTime(LevelData<EBCellFAB>&       a_U,
		  const LevelData<EBCellFAB>& a_UOld,
		  const LevelData<EBCellFAB>& a_UNew,
		  const DisjointBoxLayout&    a_grids,
		  const Real&                 a_time,
		  const Real&                 a_told,
		  const Real&                 a_tnew)
{
  CH_assert(a_tnew >= a_told);
  Real newFac;
  if(a_tnew-a_told > 1.0e-12)
    newFac = (a_time-a_told)/(a_tnew-a_told);
  else
    newFac = 1.0;
  Real oldFac = 1.0-newFac;
  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      a_U[dit()].axby(a_UOld[dit()], a_UNew[dit()], oldFac, newFac);
    }
}
/****************************/
void EBAMRSpecies::tagCells(IntVectSet& a_tags)
{
  CH_TIME("EBAMRSpecies::tagCells");
 
  // This appends the Charge threshol in regridThresh[m_nComp-1]
  Vector<Real> SpecThres(m_nComp,-1);
  for(int iSpec=0; iSpec<m_nComp; iSpec++) SpecThres[iSpec] = m_PlasmaPhysics->regridThresh(iSpec);
  ///for(int iSpec=0; iSpec<m_nComp; iSpec++) pout() << "Threshold  " << iSpec << " " << SpecThres[iSpec] << " "; pout() << endl;

  a_tags.makeEmpty();
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      const EBCellFAB& specFAB = m_stateNew[dit()];
      Box gridBox = m_grids.get(dit());
      const EBGraph& ebgraph = m_ebisl[dit()].getEBGraph();
      IntVectSet ivsTot(gridBox);
      Box shrunkDom = m_problem_domain.domainBox(); //noshrinking here(refer to the EBAMViscous implementation for shrinking)
      ivsTot &= shrunkDom;

      if (m_refineThresh > 0)
	{
	  for (VoFIterator vofit(ivsTot, ebgraph); vofit.ok(); ++vofit)
	    {
	      const VolIndex& vof = vofit();
	      const IntVect& iv = vof.gridIndex();
	      for(int iSpec=0; iSpec<m_nComp-1; iSpec++)
		{
		  Real thresh = SpecThres[iSpec];
		  if (thresh > 0)
		    {
		      Real specmag = specFAB(vof, iSpec);
		      if(specmag >= thresh) a_tags |= iv;
		    }
		}
	      Real charge = m_PlasmaPhysics->gasCharge(specFAB,vof);
	      if(abs(charge) >= SpecThres[m_nComp-1]) a_tags |= iv;
	    } //end loop over vofs
	}

      //refine all irregular cells in the chemBox
      if(m_tagAllIrreg >= 0)
	{
	  if(m_level>= m_tagAllIrreg && !m_chemBox.isEmpty()) gridBox &=  m_chemBox;
	  IntVectSet irregIVS = ebgraph.getIrregCells(gridBox);
	  a_tags |= irregIVS;
	}
    }

  a_tags.grow(m_tagBufferSize);


  // Need to do this in two steps unless a IntVectSet::operator &=
  // (ProblemDomain) operator is defined
  Box localTagsBox = a_tags.minBox();
  localTagsBox &= m_problem_domain;
  a_tags &= localTagsBox;
}

/***************************/
void EBAMRSpecies::regridWithDBL(const DisjointBoxLayout& a_DBL)
{
  CH_TIME("EBAMRSpecies::regrid");
  if (s_verbosity >= 3)
    {
      pout() << " in EBAMRSpecies regrid for level " << m_level << endl;
    }

  const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
  // save data for later copy
  //not using m_ebisl because it gets wiped later
  //the ebisl has to know about the fact that we really have
  //four ghost cells.
  int nGhostEBISL = 6;
  Interval interv(0,m_nComp-1);
  LevelData<EBCellFAB> stateSaved;
  IntVect ivGhost = m_nGhost*IntVect::Unit;
  {
    CH_TIME("defines and copies and ebisls");
    EBISLayout ebislOld;
    ebisPtr->fillEBISLayout(ebislOld, m_grids, m_domainBox, nGhostEBISL);
    EBCellFactory factoryOld(ebislOld);
    stateSaved.define(m_grids, m_nComp, ivGhost, factoryOld);
    m_stateNew.copyTo(interv, stateSaved, interv);
  }
  

  m_grids = a_DBL;
  m_level_grids = a_DBL.boxArray();

  ebisPtr->fillEBISLayout(m_ebisl, m_grids, m_domainBox, nGhostEBISL);

  EBCellFactory factoryNew(m_ebisl);
  // reshape state with new grids
  m_stateNew.define(m_grids,m_nComp,ivGhost, factoryNew);
  m_stateOld.define(m_grids,m_nComp,ivGhost, factoryNew);
  if(m_SFD) m_qbar.define(m_grids,m_nComp, ivGhost, factoryNew);

  // set up data structures
  levelSetup();

  // interpolate from coarser level (fills this level (fine))
  if (m_hasCoarser)
    {
      EBAMRSpecies* coarPtr = getCoarserLevel();
      m_ebFineInterp.interpolate(m_stateNew, coarPtr->m_stateNew, interv);
    }

  // copy from old state
  stateSaved.copyTo(interv,m_stateNew, interv);
}

/*******/
void
EBAMRSpecies::
extrapolateZeroArea()
{
  if (s_verbosity >= 2)
    {
      pout() << "EBAMRSpecies::extrapolateZeroArea " << m_level << endl;
    }

  
  //extrapolate @ zero area cells
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      Box gridBox = m_stateNew[dit()].getRegion();
      const EBISBox& ebisBox = m_ebisl[dit()];
      const EBGraph& ebgraph = m_ebisl[dit()].getEBGraph();
      EBCellFAB& fab = m_stateNew[dit()];
      EBCellFAB& fabOld = m_stateOld[dit()];
      
      IntVectSet notRegular;
      notRegular |= ebisBox.getIrregIVS  (gridBox);
      notRegular |= ebisBox.getMultiCells(gridBox);

      for (VoFIterator vofit(notRegular,ebgraph); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
	  if( ebisBox.bndryArea(vof) < 1e-6)
	    {
	      VoFStencil extrapSten;
	      RealVect bndryCentroid = m_origin; //ebisBox.bndryCentroid(vof);
	      int order = EBArith::getExtrapolationStencil(extrapSten, bndryCentroid, m_dx, vof, ebisBox);
	      for (int icomp = 0; icomp < m_nComp; icomp++)
		{
		  //Real retval=0;
		  Real retval=1e99;
		  for (int isten = 0; isten < extrapSten.size(); isten++)
		    {
                      retval = min(retval, fab(extrapSten.vof(isten),icomp));
		      //retval += (extrapSten.weight(isten))*(fab(extrapSten.vof(isten), icomp));
		    }
		  fabOld(vof,icomp) = retval;
		}
	    }
	}
    }
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      Box gridBox = m_stateNew[dit()].getRegion();
      const EBISBox& ebisBox = m_ebisl[dit()];
      const EBGraph& ebgraph = m_ebisl[dit()].getEBGraph();
      EBCellFAB& fab = m_stateNew[dit()];
      EBCellFAB& fabOld = m_stateOld[dit()];
      IntVectSet notRegular;
      notRegular |= ebisBox.getIrregIVS  (gridBox);
      notRegular |= ebisBox.getMultiCells(gridBox);
      

      for (VoFIterator vofit(notRegular,ebgraph); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
	  if( ebisBox.bndryArea(vof) < 1e-6)
	    {
              //if( vof == VolIndex(IntVect(D_DECL(11,15,0)),0)) { pout() << " Level " << m_level << " Zero Area Point " << vof;for (int icomp = 0; icomp < m_nComp; icomp++)pout() << " " <<  fab(vof,icomp); pout () << endl;}
	      for (int icomp = 0; icomp < m_nComp; icomp++)  fab(vof,icomp) = fabOld(vof,icomp);
	    }
	}
    }
}
/*******/
void
EBAMRSpecies::
postTimeStep()
{
  if (s_verbosity >= 2)
    {
      pout() << "EBAMRSpecies::postTimeStep " << m_level << endl;
    }

 CH_TIME("EBAMRSpecies::postTimeStep");


  Interval interv(0, m_nComp-1);
  if (m_hasCoarser)
    {
      if (m_doSmushing)
        {
          //redistibute to coarser level
          EBAMRSpecies* coarPtr = getCoarserLevel();

          m_ebFineToCoarRedist.redistribute(coarPtr->m_stateNew, interv);
        }

      m_ebFineToCoarRedist.setToZero();
    }
  if (m_hasFiner)
    {
      EBAMRSpecies* finePtr = getFinerLevel();
      // reflux from finer level solution
      CH_assert(Abs(m_dx[0]-m_dx[1])<1.e-9);
      Real scale = -1.0/m_dx[0];

      if (m_doImplicitReflux)
        {
          doImplicitReflux();
        }
      else
        {
	  m_ebFluxRegister.reflux(m_stateNew, interv, scale);
	}

      finePtr->m_ebCoarseAverage.average(m_stateNew, finePtr->m_stateNew, interv);


      m_ebFluxRegister.incrementRedistRegister(m_ebCoarToFineRedist, interv, scale);

      m_ebFluxRegister.incrementRedistRegister(m_ebCoarToCoarRedist, interv, scale);

      if (m_doSmushing)
        {
          //redistibute to finer level
          m_ebCoarToFineRedist.redistribute(finePtr->m_stateNew, interv);

          //do the re redistirubtion
          m_ebCoarToCoarRedist.redistribute(m_stateNew, interv);
        }//m_doSmushing

      // average from finer level data
      finePtr->m_ebCoarseAverage.average(m_stateNew, finePtr->m_stateNew, interv);
      if(m_SFD) finePtr->m_ebCoarseAverage.average(m_qbar, finePtr->m_qbar, interv);

      m_ebCoarToFineRedist.setToZero();
      m_ebCoarToCoarRedist.setToZero();
    }


}

/*******/
void
EBAMRSpecies::
doImplicitReflux()
{
  // first find out what time coarser level is (if it exists)
  Real time_eps = 1.0e-10;
  // do multilevel operations if this is the coarsest level or if
  // coarser level is not at same time as this level. otherwise,
  // defer this until we get down to the coarsest level at this time.
  if (m_level == 0 || (abs(m_coarser_level_ptr->time() - m_time) > time_eps))
    {

      Vector<EBAMRSpecies*>        hierarchy;
      Vector<DisjointBoxLayout>    grids;
      Vector<EBISLayout>           ebisl;
      Vector<EBLevelGrid>          eblg;
      Vector<int>                  refRat;
      ProblemDomain                lev0Dom;
      Real                         lev0Dx;
      Vector<ProblemDomain>        domains;
      Vector<Real>                 dxs;
      getHierarchyAndGrids(hierarchy, grids, ebisl, eblg, refRat, lev0Dom, lev0Dx, domains, dxs);

      int finest_level = hierarchy.size()-1;

      // now do implicit refluxing
      // Vector of pointers to LevelData of FABS
      Vector<LevelData<EBCellFAB>* > refluxCor(finest_level+1, NULL);
      Vector<LevelData<EBCellFAB>* > refluxRHS(finest_level+1, NULL);
      // collectUN: AMR vector containing soln at new time
      Vector<LevelData<EBCellFAB>* > collectUN(finest_level+1, NULL);

      // loop over levels, allocate storage, set up for AMRSolve
      // if coarser level exists, define it as well for BCs.
      int startLev = Max(m_level-1, 0);
      int lbase = m_level;
      int lmax  = finest_level;
      Real alpha = s_alpha;
      Real beta  = -m_dt;
      // this resets the coeffients including eta, alpha, beta (for all species)
      setSolverCoef(alpha, beta);
      for (int lev = startLev; lev<= finest_level; lev++)
	{
	  // rhs has no ghost cells, correction does
	  EBCellFactory factory( ebisl[lev]);
	  // the number of ghost cells here has to matcht that in EBSpeciesFluxOpFactory
	  refluxRHS[lev]  = new LevelData<EBCellFAB>(grids[lev], 1, m_nGhost*IntVect::Unit,factory);
	  refluxCor[lev]  = new LevelData<EBCellFAB>(grids[lev], 1, m_nGhost*IntVect::Unit,factory);
	  collectUN[lev]  = &(hierarchy[lev]->m_stateNew);
	} // end loop over levels for setup.
      
      for(int iSpec =0;iSpec<m_nComp; iSpec++)
	{ 
	  Interval interv(iSpec,iSpec);
	  for (int lev = startLev; lev<= finest_level; lev++)
	    for (DataIterator dit = grids[lev].dataIterator(); dit.ok(); ++dit)
	      {
		(*refluxRHS[lev])[dit()].setVal(0.0);
		(*refluxCor[lev])[dit()].setVal(0.0);
	      }

	  // now loop over levels and set up RHS
	  // note that we don't need to look at finest level here,
	  // since it will have no finer level to get a reflux correction
	  // from.   Also this starts at m_level instead of startLev since,
	  // in the case m_level > 0, m_level-1 is only used for boundary conditions
	  for (int lev= m_level; lev < finest_level; lev++)
	    {
	      Real dxLev = dxs[lev];
	      Real refluxScale = 1.0/dxLev;
	      hierarchy[lev]->m_ebFluxRegister.reflux(*refluxRHS[lev], Interval(0,0), interv, refluxScale);
	    }


	  s_diffuseAMRMG[iSpec]->solve(refluxCor, refluxRHS, lmax, lbase);

	  for (int lev= m_level; lev <= finest_level; lev++)
	    {
	      EBLevelDataOps::incr( (*hierarchy[lev]).m_stateNew, *refluxCor[lev], 0, iSpec,1);
	    }
	}

      //remember that startLev can be different from m_level
      for (int lev = startLev; lev<= finest_level; lev++)
        {
          delete refluxRHS[lev];
          delete refluxCor[lev];
          refluxRHS[lev] = NULL;
          refluxCor[lev] = NULL;
        }
    } //end if times are in alignment
}
/*******/
void
EBAMRSpecies::
setSolverCoef(Real a_alpha, Real a_beta)
{
  // now set alpha and beta of all the operators in the solver.
  // this includes resetting the relaxation coefficient
  
  for(int iSpec =0;iSpec<m_nComp; iSpec++)
    {      
      Vector<MGLevelOp<LevelData<EBCellFAB> >* > ops = s_diffuseAMRMG[iSpec]->getAllOperators();
      for (int iop = 0; iop < ops.size(); iop++)
	{
	  EBSpeciesFluxOp* helmop = (EBSpeciesFluxOp*) ops[iop];
	  helmop->setAlphaAndBeta(a_alpha, a_beta);
	}
    }
}
/***************************/
void EBAMRSpecies::initialGrid(const Vector<Box>& a_new_grids)
{
  CH_TIME("EBAMRG::initialGrid");

  if (s_verbosity >= 3)
    {
      pout() << " in EBAMRSpecies initialGrid for level " << m_level << endl;
    }
  Vector<Box>& newGrids =   (Vector<Box>&) a_new_grids;
  mortonOrdering(newGrids);
  const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
  m_level_grids = a_new_grids;

  // load balance and create boxlayout
  Vector<int> proc_map;
  if (s_isLoadBalanceSet)
    {
      s_loadBalance(proc_map,a_new_grids, m_domainBox, false);
    }
  else
    {
      LoadBalance(proc_map,a_new_grids);
    }

  m_grids = DisjointBoxLayout(a_new_grids,proc_map);

  if (s_verbosity >= 3)
    {
      pout() << "EBAMRSpecies::initialgrid grids " << endl;
      dumpDBL(&m_grids);
    }

  //the ebisl has to know about the fact that we really have
  //four ghost cells.
  int nGhostEBISL = 6;
  if (s_verbosity >= 3)
    {
      pout() << " EBAMRSpecies::about to fill ebislayout  in EBAMRSpecies initialGrid for level " << m_level << endl;
    }

  ebisPtr->fillEBISLayout(m_ebisl, m_grids,  m_domainBox, nGhostEBISL);

  if (s_verbosity >= 3)
    {
      pout() << " done with filling ebislayout  in EBAMRSpecies initialGrid for level " << m_level << endl;
    }

  EBCellFactory factoryNew(m_ebisl);
  IntVect ivGhost = m_nGhost*IntVect::Unit;
  m_stateNew.define(m_grids,m_nComp, ivGhost, factoryNew);
  m_stateOld.define(m_grids,m_nComp, ivGhost, factoryNew);
  if(m_SFD) m_qbar.define(m_grids,m_nComp, ivGhost, factoryNew);

  // set up data structures
  levelSetup();
}
/*****/
void EBAMRSpecies::initialGrid(const DisjointBoxLayout& a_DBL)
{
  CH_TIME("EBAMRSpecies::initialGrid(DBL)");

  

  m_grids = a_DBL;
  m_level_grids = a_DBL.boxArray();


  //the ebisl has to know about the fact that we really have
  //four ghost cells.
  int nGhostEBISL = 6;
  if (s_verbosity >= 3)
    {
      pout() << " about to fill ebislayout  in EBAMRLinear initialGrid for level " << m_level << endl;
    }

  const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
  ebisPtr->fillEBISLayout(m_ebisl, m_grids,  m_domainBox, nGhostEBISL);

  if (s_verbosity >= 3)
    {
      pout() << " done with filling ebislayout  in EBAMRLinear initialGrid for level " << m_level << endl;
    }

  EBCellFactory factoryNew(m_ebisl);
  IntVect ivGhost = m_nGhost*IntVect::Unit;
  m_stateNew.define(m_grids,m_nComp, ivGhost, factoryNew);
  m_stateOld.define(m_grids,m_nComp, ivGhost, factoryNew);
  if(m_SFD) m_qbar.define(m_grids,m_nComp, ivGhost, factoryNew);

  // set up data structures
  levelSetup();
}
/*******/
void EBAMRSpecies::levelSetup()
{
  if (s_verbosity >= 3)
    {
      pout() << " in EBAMRSpecies levelSetup for level " << m_level << endl;
    }

  EBAMRSpecies* coarPtr = getCoarserLevel();
  EBAMRSpecies* finePtr = getFinerLevel();

  m_hasCoarser = (coarPtr != NULL);
  m_hasFiner   = (finePtr != NULL);

  int nRefCrse = -1;
  DisjointBoxLayout* crseGridsPtr = NULL;

  if (m_hasCoarser)
    {
      nRefCrse = m_coarser_level_ptr->refRatio();

      const DisjointBoxLayout& coarser_m_grids = coarPtr->m_grids;
      const EBISLayout& coarser_m_ebisl = coarPtr->m_ebisl;
      const Box& domainCoar = coarPtr->m_domainBox;

      m_ebCoarseAverage.define(m_grids,
			       coarser_m_grids,
			       m_ebisl,
			       coarser_m_ebisl,
			       domainCoar,
			       nRefCrse,
			       m_nComp, Chombo_EBIS::instance());
      m_ebFineInterp.define(m_grids,
			    coarser_m_grids,
			    m_ebisl,
			    coarser_m_ebisl,
			    domainCoar,
			    nRefCrse,
			    m_nComp);
      
      // Maintain levelSpecies
      m_levelSpecies.define(m_grids,
			    coarser_m_grids,
			    m_ebisl,
			    coarser_m_ebisl,
			    m_problem_domain,
			    nRefCrse,
			    m_dx,
			    m_useMassRedist,
			    m_doSmushing,
			    m_doRZCoords,
			    m_hasSourceTerm,
			    m_ebPatchSpeciesFactory,
			    m_hasCoarser,
			    m_hasFiner);
      
      //define fine to coarse redistribution object
      //for now set to volume weighting
      m_ebFineToCoarRedist.define(m_grids, coarser_m_grids,
				  m_ebisl, coarser_m_ebisl,
				  domainCoar, nRefCrse, m_nComp);
      m_ebFineToCoarRedist.setToZero();

      coarPtr->syncWithFineLevel();
    }
  else
    {
      m_levelSpecies.define(m_grids,
			    DisjointBoxLayout(),
			    m_ebisl,
			    EBISLayout(),
			    m_problem_domain,
			    m_ref_ratio,
			    m_dx,
			    m_useMassRedist,
			    m_doSmushing,
			    m_doRZCoords,
			    m_hasSourceTerm,
			    m_ebPatchSpeciesFactory,
			    m_hasCoarser,
			    m_hasFiner);
    }

  //set up mass redistribution array
  m_sets.define(m_grids);
  for (DataIterator dit = m_grids.dataIterator();
      dit.ok(); ++dit)
    {
      Box thisBox = m_grids.get(dit());
      m_sets[dit()] = m_ebisl[dit()].getIrregIVS(thisBox);
    }
  BaseIVFactory<Real> IVfactory(m_ebisl, m_sets);
  //the ghost cells on the mass redistribution array
  //are tied to the redistribution radius
  int nghost = 2*m_redistRad;
  IntVect ivGhost = nghost*IntVect::Unit;
  m_massDiff.define(m_grids, m_nComp, ivGhost, IVfactory);
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      m_massDiff[dit()].setVal(0.0);
    }

  //advection velocity  (This is a Mac advection velocity, thus only one component)
  EBFluxFactory        ebfluxfact(m_ebisl);
  m_AdvectiveVelocity.define(m_grids, 1, m_ghost_Dk*IntVect::Unit, ebfluxfact);

}

/*******/
void
EBAMRSpecies::
initialData()
{

  if (s_verbosity >= 3)
    {
      pout() << " in EBAMRSpecies initialData for level " << m_level << endl;
    }
  
  const EBPhysIBC* const ebphysIBCPtr = m_ebPatchSpecies->getEBPhysIBC();

  //initialize both new and old states to
  //be the same thing
  ebphysIBCPtr->initialize(m_stateNew, m_ebisl);
  ebphysIBCPtr->initialize(m_stateOld, m_ebisl);
  if(m_SFD) ebphysIBCPtr->initialize(m_qbar, m_ebisl);
}
/*******/
void
EBAMRSpecies::
postInitialize()
{
  if (s_verbosity >= 3)
    {
      pout() << "EBAMRSpecies::postInitialize " << m_level << endl;
    }

  // solve for post-regrid smoothing may eventually go here.
}
/*******/

Real
EBAMRSpecies::
computeDt()
{
  Real newDt;
  newDt = m_dtNew;

  return newDt;
}

/*******/
Real
EBAMRSpecies::
computeInitialDt()
{
  //quick return if the velocity is not initialized
  if(m_level >= m_Velo.size()) return 1e0;

  Real newDT = m_initial_dt_multiplier * m_dx[0]
    / m_levelSpecies.getMaxWaveSpeed(m_stateNew, m_Velo[m_level]);

  return newDT;
}

/*******/
EBAMRSpecies*
EBAMRSpecies::
getCoarserLevel() const
{
  EBAMRSpecies* amrADCoarserPtr = NULL;

  if (m_coarser_level_ptr != NULL)
    {
      amrADCoarserPtr = dynamic_cast<EBAMRSpecies*>(m_coarser_level_ptr);

      if (amrADCoarserPtr == NULL)
        {
          MayDay::Error("EBAMRSpecies::getCoarserLevel: dynamic cast failed");
        }
    }

  return amrADCoarserPtr;
}

/*******/
EBAMRSpecies*
EBAMRSpecies::
getFinerLevel() const
{
  EBAMRSpecies* amrADFinerPtr = NULL;

  if (m_finer_level_ptr != NULL)
    {
      amrADFinerPtr = dynamic_cast<EBAMRSpecies*>(m_finer_level_ptr);

      if (amrADFinerPtr == NULL)
        {
          MayDay::Error("EBAMRSpecies::getFinerLevel: dynamic cast failed");
        }
    }

  return amrADFinerPtr;
}
/***********/
void EBAMRSpecies::syncWithFineLevel()
{
  CH_TIME("EBAMRG::syncWithFineLevel");
  if (s_verbosity >= 3)
    {
      pout() << " in EBAMRSpecies syncWithFineLevel for level " << m_level << endl;
    }
  //stuff that needs to be setup from the finer
  //level.  A bunch of objects depend upon the layouts
  //from both levels and the finer level changes more
  //often from regrid so this needs to be called from the finer
  //level
  CH_assert(m_hasFiner);
  if (m_hasFiner)
    {
      EBAMRSpecies* finePtr = getFinerLevel();
      int nRefFine = refRatio();
      const DisjointBoxLayout& finer_m_grids = finePtr->m_grids;
      const EBISLayout& finer_m_ebisl = finePtr->m_ebisl;
      //define fine to coarse redistribution object
      //for now set to volume weighting
      m_ebCoarToFineRedist.define(finer_m_grids, m_grids,
                                  finer_m_ebisl,  m_ebisl,
                                  m_domainBox, nRefFine , m_nComp, 1, Chombo_EBIS::instance());
      //define coarse to coarse redistribution object
      m_ebCoarToCoarRedist.define(finer_m_grids, m_grids,
                                  finer_m_ebisl,  m_ebisl,
                                  m_domainBox, nRefFine , m_nComp);
      // maintain flux registers
      m_ebFluxRegister.define(finer_m_grids,
                              m_grids,
                              finer_m_ebisl,
                              m_ebisl,
                              m_domainBox,
                              nRefFine,
                              m_nComp, Chombo_EBIS::instance());

      //set all the registers to zero
      m_ebCoarToFineRedist.setToZero();
      m_ebCoarToCoarRedist.setToZero();
      m_ebFluxRegister.setToZero();
    }

}
/*******/
void
EBAMRSpecies::
fillAdvectionVelocity()
{
  
  EBLevelDataOps::averageCellToFacesMAC(m_AdvectiveVelocity, *m_Velo[m_level], m_grids, m_ebisl,
					m_problem_domain);
  ccpExtrapolateToDomainBoundaries(m_AdvectiveVelocity,m_grids, m_ebisl, m_problem_domain, m_dx);
  
}
/********/
void
EBAMRSpecies::getElectronEnergy(LevelData<EBCellFAB>&  a_temp, int a_tComp)
{
  int enComp = m_PlasmaPhysics->getElectronEnergyComponent();
  m_stateNew.copyTo(Interval(enComp,enComp),a_temp,Interval(a_tComp,a_tComp));
}
/********/
void EBAMRSpecies::setDiffCoeff(Vector<EBAMRSpecies*>&      a_hierarchy,
			       Vector<DisjointBoxLayout>&   a_grids,
			       Vector<EBISLayout>&          a_ebisl,
			       Vector<EBLevelGrid>&         a_eblg,
			       Vector<int>&                 a_refRat,
			       ProblemDomain&               a_lev0Dom,
			       Vector<Real>&                a_dxs,
			       Vector<ProblemDomain>&       a_domains)
{ 

  m_nlevels = a_hierarchy.size();
  s_vecW.resize(m_nComp);
  s_Dk.resize(m_nComp);
  s_aco.resize(m_nComp);
  s_DkIrreg.resize(m_nComp);
  s_vecWIrreg.resize(m_nComp);

  // legend::
  // s_data is the flux to the wall; i.e., if s_data = -100, there is a flux of 100 mole/sec/m^2 into the fluid
  // vecW is the +mu*grad phi; i.e., the charged particles move in the direction opposite to vecW
  // D_k is the diffusion coeffs, i.e. dn_k/dt = + D_k nabla n_k + vecW n_k

  s_data.resize(m_nlevels);
  s_rhosRHS.resize(m_nlevels);
  for(int k=0; k<m_nComp; k++)
    { 
      s_vecW[k].resize(m_nlevels);
      s_Dk[k].resize(m_nlevels);
      s_aco[k].resize(m_nlevels);
      s_DkIrreg[k].resize(m_nlevels);
      s_vecWIrreg[k].resize(m_nlevels);
    }

  // askMD:: I need here the number of neutral+ions+electron+electronEnergy  
  int nPlasmaCoeffs = m_PlasmaPhysics->nComponents();


  int nghost = m_nGhost; //for the conserved variables
  int nCons = m_nComp;   //note should be the same as nPlasmaCoeffs:: ie, I find the coefficients for all species
  Interval consInterv(0, nCons-1);
  for (int ilev = 0; ilev < m_nlevels; ilev++)
    {  
      Real dx = a_dxs[ilev];
      RealVect dxv = dx*RealVect::Unit;

      EBCellFactory factory(a_ebisl[ilev]);
      // askMD:: can this pass the whole conserved state or should I alias
      LevelData<EBCellFAB>& stateU  = a_hierarchy[ilev]->getStateNew();
      LevelData<EBCellFAB> cellNData(a_grids[ilev], nCons, (IntVect::Unit)*nghost, factory);
      LevelData<EBCellFAB>& cellphiData  = *m_phi[ilev];
      LevelData<EBCellFAB>& cellGphiData  = *m_cellGradPhi[ilev];
      LevelData<BaseIVFAB<Real>>& IVGphiData  = *m_IVGradPhi[ilev];
      LevelData<EBCellFAB>& cellTData  = *m_Temp[ilev];

      
      EBcopyLevelData(cellNData, consInterv, stateU, consInterv, false);
      //stateU.copyTo(consInterv, cellNData, consInterv);
      if (ilev >0)
	{
	  EBPWLFillPatch patcher(a_grids[ilev],
				 a_grids[ilev-1],
				 a_ebisl[ilev],
				 a_ebisl[ilev-1],
				 a_domains[ilev-1].domainBox(),
				 a_refRat[ilev-1], m_nComp, nghost);
      
	  Real coarTimeOld = 0.0;
	  Real coarTimeNew = 1.0;
	  Real fineTime    = 1.0;
	  patcher.interpolate(cellNData,
			      a_hierarchy[ilev-1]->m_stateOld,
			      a_hierarchy[ilev-1]->m_stateNew,
			      coarTimeOld,
			      coarTimeNew,
			      fineTime,
			      consInterv);
	}
      cellNData.exchange(consInterv);
      m_PlasmaPhysics->floor(cellNData);
      //EBLevelDataOps::setCoveredVal (cellNData, 0e0);
      EBLevelDataOps::setCoveredVal (cellTData, 0e0);
      EBLevelDataOps::setCoveredVal (cellphiData, 0e0);
      EBLevelDataOps::setCoveredVal (*m_gradPhi[ilev], 0e0);

      LayoutData<IntVectSet> irregSets(a_grids[ilev]);

      for (DataIterator dit = cellNData.dataIterator(); dit.ok(); ++dit) 
	{
	  Box grownBox = grow(a_grids[ilev].get(dit()), m_ghost_Dk);
	  grownBox &= a_domains[ilev];
	  irregSets[dit()] = a_ebisl[ilev][dit()].getIrregIVS(grownBox);
	}
      EBCellFactory        ebcellfact(a_ebisl[ilev]);
      EBFluxFactory        ebfluxfact(a_ebisl[ilev]);
      BaseIVFactory<Real>  baseivfact(a_ebisl[ilev], irregSets);

      LevelData<EBFluxFAB> fluxMobData(a_grids[ilev], nPlasmaCoeffs, m_ghost_Dk*IntVect::Unit, ebfluxfact);
      LevelData<EBFluxFAB> fluxDiffData(a_grids[ilev], nPlasmaCoeffs, m_ghost_Dk*IntVect::Unit, ebfluxfact);
      LevelData< BaseIVFAB<Real> > fluxMobIV(a_grids[ilev], nPlasmaCoeffs, m_ghost_Dk*IntVect::Unit, baseivfact);
      LevelData< BaseIVFAB<Real> > fluxDiffIV(a_grids[ilev], nPlasmaCoeffs, m_ghost_Dk*IntVect::Unit, baseivfact);

      s_data[ilev]    = RefCountedPtr< LevelData< BaseIVFAB<Real> > >(new LevelData< BaseIVFAB<Real> >(a_grids[ilev], nPlasmaCoeffs+2, m_ghost_Dk*IntVect::Unit, baseivfact));
      LevelData< BaseIVFAB<Real> >& wallData = *s_data[ilev];
      for (DataIterator dit = wallData.dataIterator(); dit.ok(); ++dit)
	{
	  wallData[dit()].setVal(0e0);
	}
      s_rhosRHS[ilev]    = RefCountedPtr< LevelData< BaseIVFAB<Real> > >(new LevelData< BaseIVFAB<Real> >(a_grids[ilev], 1, m_ghost_Dk*IntVect::Unit, baseivfact));
      LevelData< BaseIVFAB<Real> >& rhosRHS = *s_rhosRHS[ilev];
      for (DataIterator dit = rhosRHS.dataIterator(); dit.ok(); ++dit) rhosRHS[dit()].setVal(0e0);
      
      
      //gets the boxed FABS
      m_PlasmaPhysics->transportCoefficients(fluxMobData, fluxDiffData, cellNData, cellTData, cellphiData, *m_gradPhi[ilev],  a_grids[ilev], a_ebisl[ilev], a_domains[ilev],a_dxs[ilev]);
      //gets the IV FABS
      m_PlasmaPhysics->transportCoefficients(fluxMobIV, fluxDiffIV, wallData, rhosRHS, cellNData, cellTData, IVGphiData, cellGphiData, a_grids[ilev], a_ebisl[ilev], a_domains[ilev], a_dxs[ilev]);
      
      //change the sign of the mobility coeffs (kind of irrational, it would make more sense to calculate on the RHS)
      //EBLevelDataOps::scale(fluxMobData,-1.0); // PlasmaPhysics has them as -Gradphi*mu
      EBLevelDataOps::ABS(fluxDiffData);
      for (DataIterator dit = fluxMobIV.dataIterator(); dit.ok(); ++dit)// MD has them as -Gradphi*mu
	{
	  const EBGraph&  ebgraph = a_ebisl[ilev][dit()].getEBGraph();
	  for (VoFIterator vofit(irregSets[dit()], ebgraph); vofit.ok(); ++vofit)
	    {
	      VolIndex vof = vofit();
	      for(int k=0; k<nPlasmaCoeffs; k++) fluxMobIV[dit()](vof,k) *= -1.0;
	      for(int k=0; k<nPlasmaCoeffs; k++) fluxDiffIV[dit()](vof,k) = abs(fluxDiffIV[dit()](vof,k));
	  
	    }
	}


      EBLevelDataOps::setCoveredVal (fluxMobData, 0e0);
      EBLevelDataOps::setCoveredVal (fluxDiffData, 0e0);
      wallData.exchange();
      rhosRHS.exchange();


      // fill diffusion pointers
      for(int k=0; k<nPlasmaCoeffs; k++)
	{	  
	  Interval MDinterv(k,k);
	  Interval LMinterv(0,0);
      
	  s_vecW[k][ilev]    = RefCountedPtr< LevelData<EBFluxFAB> >(new LevelData<EBFluxFAB>(a_grids[ilev], 1, m_ghost_Dk*IntVect::Unit, ebfluxfact));
	  s_Dk[k][ilev]      = RefCountedPtr< LevelData<EBFluxFAB> >(new LevelData<EBFluxFAB>(a_grids[ilev], 1, m_ghost_Dk*IntVect::Unit, ebfluxfact));
	  s_aco[k][ilev]     = RefCountedPtr< LevelData<EBCellFAB> >(new LevelData<EBCellFAB>(a_grids[ilev], 1, m_ghost_Dk*IntVect::Unit, ebcellfact));
	  s_DkIrreg[k][ilev] = RefCountedPtr< LevelData< BaseIVFAB<Real> > >(new LevelData< BaseIVFAB<Real> >(a_grids[ilev], 1, m_ghost_Dk*IntVect::Unit, baseivfact));
	  s_vecWIrreg[k][ilev]=RefCountedPtr< LevelData< BaseIVFAB<Real> > >(new LevelData< BaseIVFAB<Real> >(a_grids[ilev], 1, m_ghost_Dk*IntVect::Unit, baseivfact));

	  LevelData<EBCellFAB>&                aco = *s_aco[k][ilev];
	  LevelData<EBFluxFAB>&                Dk  = *s_Dk [k][ilev];
	  LevelData<EBFluxFAB>&               vecW = *s_vecW[k][ilev];
	  LevelData< BaseIVFAB<Real> >&  DkIrreg   = *s_DkIrreg[k][ilev];
	  LevelData< BaseIVFAB<Real> >&  vecWIrreg = *s_vecWIrreg[k][ilev];
      

	  EBLevelDataOps::setVal(Dk,0e0);  
	  EBLevelDataOps::setVal(vecW,0e0);
	  EBLevelDataOps::setVal(aco,1e0); 
	  for (DataIterator dit = Dk.dataIterator(); dit.ok(); ++dit)
	    {
	      (DkIrreg)  [dit()].setVal(0e0);
	      (vecWIrreg)[dit()].setVal(0e0);
	    }

	  //fluxMobData. copyTo(MDinterv, vecW, LMinterv);
	  //fluxDiffData.copyTo(MDinterv, Dk  , LMinterv);
	  //fluxMobIV. copyTo(MDinterv, vecWIrreg, LMinterv);
	  //fluxDiffIV.copyTo(MDinterv, DkIrreg  , LMinterv);
          if(k >=0){
	    EBcopyLevelData(vecW,LMinterv,fluxMobData, MDinterv);
	    EBcopyLevelData(Dk,  LMinterv,fluxDiffData,MDinterv);
	    EBcopyLevelData(vecWIrreg,LMinterv,fluxMobIV, MDinterv);
	    EBcopyLevelData(DkIrreg,  LMinterv,fluxDiffIV,MDinterv);
	    ccpExtrapolateToDomainBoundaries(vecW,a_grids[ilev], a_ebisl[ilev], a_domains[ilev], a_dxs[ilev]*RealVect::Unit);
	    ccpExtrapolateToDomainBoundaries(Dk,a_grids[ilev], a_ebisl[ilev], a_domains[ilev], a_dxs[ilev]*RealVect::Unit);
	  }
	  
	  // exchanges
	  //aco    .exchange(Interval(0,0));
	  //Dk     .exchange(Interval(0,0));
	  //vecW   .exchange(Interval(0,0));
	  //DkIrreg   .exchange(Interval(0,0));
	  //vecWIrreg .exchange(Interval(0,0));
	}
    }

  // Do this only if working w/radial symmatry
  if(m_doRZCoordsCANC) multiplyCoeffsByRadius( a_dxs, a_domains);
  
  Vector<RefCountedPtr<LevelData<EBFluxFAB> > > mobdata(m_nlevels); 
  Vector<RefCountedPtr<LevelData<EBFluxFAB> > > diffdata(m_nlevels); 
  Vector<RefCountedPtr<LevelData<EBCellFAB> > > acofdata(m_nlevels); 
  int k=m_PlasmaPhysics->findSpec(string("e-"));
  for (int ilev = 0; ilev < m_nlevels; ilev++)
    {  
      mobdata[ilev] = s_vecW[k][ilev];
      diffdata[ilev] = s_Dk[k][ilev];
      acofdata[ilev] = s_aco[k][ilev];
    }
  EBAMRDataOps::averageDown(mobdata,a_eblg,a_refRat);
  EBAMRDataOps::averageDown(diffdata,a_eblg,a_refRat);
  EBAMRDataOps::averageDown(acofdata,a_eblg,a_refRat);
  
  
}

//**************
//Necessary for RZ problems
void EBAMRSpecies::multiplyCoeffsByRadius(const Vector<Real>&           a_dxs,
					  const Vector<ProblemDomain>&  a_domains)
{

  //LM note: the distance function is in gasModel.H
  Real minRadius = m_dx[1]/1000;//to avoid division by zero
  int nPlasmaCoeffs = m_PlasmaPhysics->nComponents();
  RealVect cylinderAxis = BASISREALV(0);
  int icomp = 0;
  for (int ilev = 0; ilev < m_nlevels; ilev++)
    {
      for(int k=0; k<nPlasmaCoeffs; k++)
	{
      
	LevelData<EBCellFAB>&                aco = *s_aco[k][ilev];
	LevelData<EBFluxFAB>&                Dk  = *s_Dk [k][ilev];
	LevelData<EBFluxFAB>&               vecW = *s_vecW[k][ilev];
	LevelData< BaseIVFAB<Real> >&  DkIrreg   = *s_DkIrreg[k][ilev];
	LevelData< BaseIVFAB<Real> >&  vecWIrreg = *s_vecWIrreg[k][ilev];
	LevelData< BaseIVFAB<Real> >&  wallData  = *s_data[ilev];
	Real dx = a_dxs[ilev];
	RealVect vectDx = dx*(RealVect::Unit);
	
	
	for (DataIterator dit = Dk.dataIterator(); dit.ok(); ++dit)
	  {
	    //F-part
	    for (int idir = 0; idir < SpaceDim; idir++)
	      {
		FaceStop::WhichFaces stopCritGrid =  FaceStop::SurroundingWithBoundary;
		EBFaceFAB& dataEBFAB = Dk[dit()][idir];
		EBFaceFAB& dataEBFAB2 = vecW[dit()][idir];
		const Box& region = dataEBFAB.getCellRegion();
		IntVectSet ivsBox(region);
		const EBISBox& ebisBox = dataEBFAB.getEBISBox();
		for (FaceIterator faceit(ivsBox, ebisBox.getEBGraph(), idir, stopCritGrid); faceit.ok(); ++faceit)
		  {
		    const FaceIndex& face = faceit();
		    RealVect centroid = ebisBox.centroid(face);
		    centroid *= dx;
		    RealVect faceloc = EBArith::getFaceLocation(face, vectDx, centroid);
		    Real dist = Max(getDistanceFromAxis(faceloc, cylinderAxis, RealVect::Zero),minRadius);
		    dataEBFAB(face,icomp) *= dist;
		    dataEBFAB2(face,icomp) *= dist;
		  }
	      }// end flux part
	    EBCellFAB& dataEBFAB = aco[dit()];
	    const EBISBox& ebisBox = dataEBFAB.getEBISBox();
	    const Box& region = dataEBFAB.getRegion();
	    Box grownBox = grow(region, 0);
	    grownBox &= a_domains[ilev];
	    //I-part
	    IntVectSet ivsIrreg = ebisBox.getIrregIVS(grownBox);
	    for (VoFIterator vofit(ivsIrreg, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
	      {
		RealVect centroid = ebisBox.bndryCentroid(vofit());
		centroid *= dx;
		RealVect vofloc = EBArith::getVofLocation(vofit(), vectDx, centroid);
		Real dist = Max(getDistanceFromAxis(vofloc, cylinderAxis, RealVect::Zero), minRadius);
		DkIrreg[dit()](vofit(), icomp) *= dist;
		vecWIrreg[dit()](vofit(), icomp) *= dist;
		if(k==0) for (int jcomp = 0; jcomp < wallData.nComp(); jcomp++)wallData[dit()](vofit(), jcomp) *= dist;
	      }
	    //C-part
	    IntVectSet ivsBox(region);
	    for (VoFIterator vofit(ivsBox, ebisBox.getEBGraph());vofit.ok(); ++vofit)
            {
              const VolIndex& vof = vofit();
	      RealVect centroid = ebisBox.centroid(vof);
	      centroid *= dx;
              RealVect vofloc = EBArith::getVofLocation(vofit(), vectDx, centroid);
              Real dist = Max( getDistanceFromAxis(vofloc, cylinderAxis, RealVect::Zero), minRadius);
	      dataEBFAB(vof,icomp) *= dist;
            }

	  }

	aco.exchange(Interval(0,0));
	wallData.exchange();
      }
    }
  
}

void EBAMRSpecies::EBcopyLevelData(LevelData<EBCellFAB>&a_copiedTo, const LevelData<EBCellFAB>&a_copyFrom, 
				    const Interval& a_fromInterv,const EBISLayout& a_ebisl)
{
  EBCellFactory factory(a_ebisl);
  //a_copiedTo.define(*a_copyFrom,a_fromInterv,factory);
  Interval toInterv(0,a_fromInterv.size()-1);
  a_copiedTo.define(a_copyFrom.disjointBoxLayout(), a_fromInterv.size() , a_copyFrom.ghostVect(), factory);
  EBcopyLevelData(a_copiedTo, toInterv, a_copyFrom, a_fromInterv,true);
}
/****/
void EBAMRSpecies::EBcopyLevelData(LevelData<EBCellFAB>&a_copiedTo, const Interval& a_toInterv, const LevelData<EBCellFAB>&a_copyFrom, const Interval& a_fromInterv, const bool doClone)
{
  for (DataIterator dit = a_copyFrom.dataIterator(); dit.ok(); ++dit)
    {
      //if(doClone) a_copiedTo[dit()].define(a_copyFrom[dit()].getEBISBox(), a_copyFrom[dit()].getRegion(), a_toInterv.size());
      Box Bregion = a_copyFrom[dit()].getRegion(); //includes ghost cells
      a_copiedTo[dit()].copy(Bregion, a_toInterv, Bregion, a_copyFrom[dit()], a_fromInterv);
    }
  EBLevelDataOps::setCoveredVal(a_copiedTo,0e0);
}
void EBAMRSpecies::EBcopyLevelData(LevelData<EBFluxFAB>&a_copiedTo, const Interval& a_toInterv, const LevelData<EBFluxFAB>&a_copyFrom, const Interval& a_fromInterv, const bool doClone)
{
  for (DataIterator dit = a_copyFrom.dataIterator(); dit.ok(); ++dit)
    {
      //if(doClone) a_copiedTo[dit()].define(a_copyFrom[dit()].getEBISBox(), a_copyFrom[dit()].getRegion(), a_toInterv.size());
      Box Bregion = a_copyFrom[dit()].getRegion(); //includes ghost cells
      a_copiedTo[dit()].copy(Bregion, a_toInterv, Bregion, a_copyFrom[dit()], a_fromInterv);
    }
  EBLevelDataOps::setCoveredVal(a_copiedTo,0e0);
}
void EBAMRSpecies::EBcopyLevelData(LevelData<BaseIVFAB<Real> >&a_copiedTo, const Interval& a_toInterv, const LevelData<BaseIVFAB<Real> >&a_copyFrom, const Interval& a_fromInterv, const bool doClone)
{
  int compSize = a_fromInterv.size();
  for (DataIterator dit = a_copyFrom.dataIterator(); dit.ok(); ++dit)
    {
      const EBGraph&  ebgraph = a_copyFrom[dit()].getEBGraph();
      for (VoFIterator vofit(a_copyFrom[dit()].getIVS(), ebgraph); vofit.ok(); ++vofit)
	{
	  VolIndex vof = vofit();
	  for (int icomp = 0; icomp < compSize; icomp++)
            {
              int isrccomp = a_fromInterv.begin() + icomp;
              int idstcomp = a_toInterv.begin() + icomp;
              a_copiedTo[dit()](vof, idstcomp) = a_copyFrom[dit()](vof, isrccomp);
            }
	}
    }
}
/***************************/
LevelData<EBCellFAB>&
EBAMRSpecies::getStateNew()
{
  return m_stateNew;
}
/***************************/
void
EBAMRSpecies::
doSmushing(bool a_doSmushing)
{
  m_doSmushing = a_doSmushing;
}
/***************************/
void
EBAMRSpecies::
doRZCoords(bool a_doRZCoords)
{
  m_doRZCoords = a_doRZCoords;
}
/***************************/
void
EBAMRSpecies::
hasSourceTerm(bool a_hasSourceTerm)
{
  m_hasSourceTerm = a_hasSourceTerm;
}
/***************************/
void
EBAMRSpecies::
useMassRedistribution(bool a_useMassRedist)
{
  m_useMassRedist = a_useMassRedist;
}
/***************************/
void EBAMRSpecies::patchSpecies(const EBPatchSpeciesFactory* const a_ebPatchSpeciesFactory)
{
  m_ebPatchSpeciesFactory = a_ebPatchSpeciesFactory;
}
/*******/

/***************************/
void EBAMRSpecies::refinementThreshold(Real a_refineThresh)
{
  m_refineThresh = a_refineThresh;
}
/***************************/
void EBAMRSpecies::tagBufferSize(int a_tagBufferSize)
{
  m_tagBufferSize = a_tagBufferSize;
}
/***************************/
void EBAMRSpecies::domainLength(RealVect a_domainLength)
{
  m_domainLength = a_domainLength;
}
/***************************/
void EBAMRSpecies::redistRadius(int a_redistRad)
{
  m_redistRad = a_redistRad;
}
/***************************/
// Unused routines below
// need to make the class concrete, ie( Not abstract)
void EBAMRSpecies::tagCellsInit(IntVectSet& a_tags)
{
  if (s_verbosity >= 3)
    {
      pout() << "EBAMRSpecies::tagCellsInit for level " << m_level << endl;
    }

  // Since tags are calculated using only current time step data, use
  // the same tagging function for initialization and for regridding.
  tagCells(a_tags);
}
/***************************/
void EBAMRSpecies::regrid(const Vector<Box>& a_new_grids)
{
  CH_TIME("EBAMRSpecies::regrid");
  if (s_verbosity >= 3)
    {
      pout() << " in EBAMRSpecies regrid for level " << m_level << endl;
    }
  Vector<Box>& newGrids =   (Vector<Box>&) a_new_grids;
  {
    CH_TIME("mortonordering");
    mortonOrdering(newGrids);
  }

  const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
  // save data for later copy
  //not using m_ebisl because it gets wiped later
  //the ebisl has to know about the fact that we really have
  //four ghost cells.
  int nGhostEBISL = 6;
  Interval interv(0,m_nComp-1);
  LevelData<EBCellFAB> stateSaved;
  LevelData<EBCellFAB> stateSavedBar;
  IntVect ivGhost = m_nGhost*IntVect::Unit;
  {
    CH_TIME("defines and copies and ebisls");
    EBISLayout ebislOld;
    ebisPtr->fillEBISLayout(ebislOld, m_grids, m_domainBox, nGhostEBISL);
    EBCellFactory factoryOld(ebislOld);
    stateSaved.define(m_grids, m_nComp, ivGhost, factoryOld);
    m_stateNew.copyTo(interv, stateSaved, interv);
    if(m_SFD) 
      {
	stateSavedBar.define(m_grids, m_nComp, ivGhost, factoryOld);
	m_qbar.copyTo(interv, stateSavedBar, interv);
      }
  }
  //create grids and ebis layouts
  //the ebisl has to know about the fact that we really have
  //four ghost cells.
  m_level_grids = a_new_grids;
  Vector<int> proc_map;

  if (s_isLoadBalanceSet)
    {
      s_loadBalance(proc_map,a_new_grids, m_domainBox, false);
    }
  else
    {
      LoadBalance(proc_map,a_new_grids);
    }

  m_grids= DisjointBoxLayout(a_new_grids,proc_map);

  ebisPtr->fillEBISLayout(m_ebisl, m_grids, m_domainBox, nGhostEBISL);

  EBCellFactory factoryNew(m_ebisl);
  // reshape state with new grids
  m_stateNew.define(m_grids,m_nComp,ivGhost, factoryNew);
  m_stateOld.define(m_grids,m_nComp,ivGhost, factoryNew);
  if(m_SFD) m_qbar.define(m_grids,m_nComp,ivGhost, factoryNew);

  // set up data structures
  levelSetup(); //@regrid

  // interpolate to coarser level
  if (m_hasCoarser)
    {
      EBAMRSpecies* coarPtr = getCoarserLevel();
      m_ebFineInterp.interpolate(m_stateNew, coarPtr->m_stateNew, interv);
      if(m_SFD) m_ebFineInterp.interpolate(m_qbar, coarPtr->m_qbar, interv);
    }

  // copy from old state
  stateSaved.copyTo(interv,m_stateNew, interv);
  if(m_SFD) stateSavedBar.copyTo(interv,m_qbar, interv);
  // defineSolvers(); //defineSolvers @regrid (cant do this without regridding phi, should add phi to EBAMRSPECIES)
}/***************************/

/***************************/
void
EBAMRSpecies::setDt(const Real& a_dt)
{
  // for now this variable is not used
   m_dtExtern = a_dt;
}

void
EBAMRSpecies::
cellGrad(LevelData<EBCellFAB>&        a_gradPhi,
         const LevelData<EBCellFAB>&  a_phi,
         const int&                   a_comp)
{
  int ibox = 0;
  EBLevelDataOps::setToZero(a_gradPhi); 
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit, ibox++)
    {
      EBCellFAB& gradPhi = a_gradPhi[dit()];
      const EBCellFAB& phi = a_phi[dit()];

      const Box& cellBox = m_grids.get(dit());
      const EBISBox& ebisBox = m_ebisl[dit()];
      IntVectSet ivsIrreg = ebisBox.getIrregIVS(cellBox);
      
      for (int derivDir = 0; derivDir < SpaceDim; derivDir++)
	{
	  Box loBox, hiBox, centerBox;
	  int hasLo, hasHi;
	  EBArith::loHiCenter(loBox, hasLo, hiBox, hasHi, centerBox, m_problem_domain , cellBox, derivDir);

	  
          //
          BaseFab<Real>&       regGrad = gradPhi.getSingleValuedFAB();
          const BaseFab<Real>& regPhi  =     phi.getSingleValuedFAB();
          FORT_CELLGRADEBSPEC(CHF_FRA1(regGrad, derivDir),
                              CHF_CONST_FRA1(regPhi, a_comp),
                              CHF_CONST_REAL(m_dx[derivDir]),
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
		  gradPhi(vof, derivDir) = (valHi - valLo)/(2.*m_dx[derivDir]);
		}
	      else if (hasHiPt)
		{
		  gradPhi(vof, derivDir) = (valHi - valCe)/(m_dx[derivDir]);
		}
	      else if (hasLoPt)
		{
		  gradPhi(vof, derivDir) = (valCe - valLo)/(m_dx[derivDir]);
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
/************************/
void
EBAMRSpecies::
cellGradTens(LevelData<EBCellFAB>&        a_gradPhi,
	     const LevelData<EBCellFAB>&  a_phi,
	     const int&                   a_nComp)
{
  int ibox = 0;
  EBLevelDataOps::setToZero(a_gradPhi); 
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit, ibox++)
    {
      EBCellFAB& gradPhi = a_gradPhi[dit()];
      const EBCellFAB& phi = a_phi[dit()];

      const Box& cellBox = m_grids.get(dit());
      const EBISBox& ebisBox = m_ebisl[dit()];
      IntVectSet ivsIrreg = ebisBox.getIrregIVS(cellBox);
      
      for (int derivDir = 0; derivDir < SpaceDim; derivDir++)
	{
	  Box loBox, hiBox, centerBox;
	  int hasLo, hasHi;
	  EBArith::loHiCenter(loBox, hasLo, hiBox, hasHi, centerBox, m_problem_domain , cellBox, derivDir);

	  
          //
          BaseFab<Real>&       regGrad = gradPhi.getSingleValuedFAB();
          const BaseFab<Real>& regPhi  =     phi.getSingleValuedFAB();
	  
	  for (int icomp = 0; icomp < a_nComp; icomp++)
	    {
	      int gradComp = icomp*SpaceDim+derivDir;
	      FORT_CELLGRADEBSPEC(CHF_FRA1(regGrad, gradComp),
				  CHF_CONST_FRA1(regPhi, icomp),
				  CHF_CONST_REAL(m_dx[derivDir]),
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
		  if (hasLoPt) valLo = phi(loFaces[0].getVoF(Side::Lo), icomp);
		  if (hasHiPt) valHi = phi(hiFaces[0].getVoF(Side::Hi), icomp);
		  
		  Real valCe = phi(vof, icomp);
		  
		  if (hasLoPt && hasHiPt)
		    {
		      gradPhi(vof, gradComp) = (valHi - valLo)/(2.*m_dx[derivDir]);
		    }
		  else if (hasHiPt)
		    {
		      gradPhi(vof, gradComp) = (valHi - valCe)/(m_dx[derivDir]);
		    }
		  else if (hasLoPt)
		    {
		      gradPhi(vof, gradComp) = (valCe - valLo)/(m_dx[derivDir]);
		    }
		  else
		    {
		      gradPhi(vof, gradComp) = 0.0;
		    }
		}
	    }
	}
    }
  a_gradPhi.exchange();
}

void EBAMRSpecies::SharfetterGummel(LevelData<EBFluxFAB>&      a_fluxMobData,
				    LevelData<EBFluxFAB>&      a_fluxDiffData,
				    const Real&		       a_dx)
{ 
  int nspec = m_PlasmaPhysics->nComponents();
  for (DataIterator dit = a_fluxMobData.dataIterator(); dit.ok(); ++dit)
    {
      for (int idir = 0;idir<SpaceDim;idir++)
        {
          EBFaceFAB&       faceMobData = a_fluxMobData[dit()][idir];
          EBFaceFAB&       faceDiffData = a_fluxDiffData[dit()][idir];
	  BaseFab<Real>&   regFaceMobData =  faceMobData.getSingleValuedFAB();
	  BaseFab<Real>&   regFaceDiffData = faceDiffData.getSingleValuedFAB();
	  Box locBox = faceMobData.getRegion();
	  FORT_SHARFETTERGUMMEL_FD(CHF_FRA(regFaceMobData),
				CHF_FRA(regFaceDiffData),
				CHF_BOX(locBox),
				CHF_CONST_REAL(a_dx),
				CHF_CONST_INT(nspec));
        }
    }
  
}
/************************/
#ifdef CH_USE_HDF5
/***************************/
void EBAMRSpecies::writeCheckpointHeader(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
    {
      pout() << "EBAMRSpecies::writeCheckpointHeader" << endl;
    }

  //stuff in non-eb checkpoint header is already
  //set in the define function, such as
  // Setup the number of components
  // Setup the component names
  //so i will eliminate the middleman
}
/***************************/
void EBAMRSpecies::readCheckpointHeader(HDF5Handle& a_handle)
{
  if (s_verbosity >= 3)
    {
      pout() << "EBAMRSpecies::readCheckpointHeader" << endl;
    }

  //stuff in non-eb checkpoint header is already
  //set in the define function, such as
  // Setup the number of components
  // Setup the component names
  // So i will eliminate the middleman.
}
/***************************/
void EBAMRSpecies::writeCheckpointLevel(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
    {
      pout() << "EBAMRSpecies::writeCheckpointLevel" << endl;
    }

  //reasons for deviations from non-eb stuff
  // the tag buffer size is set by the factory.
  // dx is set in the define function
  // the problem domain is set in the define function

  // Setup the level string
  char levelStr[20];
  sprintf(levelStr,"%d",m_level);
  const std::string label = std::string("SPlevel_") + levelStr;

  a_handle.setGroup(label);

  // Setup the level header information
  HDF5HeaderData header;

  header.m_int ["ref_ratioS"]       = m_ref_ratio;
  header.m_real["dtS"]              = m_dt;

  if (s_verbosity >= 3)
    {
      pout() << "EBAMRSpecies::writeCheckpointLevel " << label << endl;
    }

  // Write the header for this level
  header.writeToFile(a_handle);

  write(a_handle,m_grids);
  write(a_handle,m_stateOld,"dataOldSpec");
  write(a_handle,m_stateNew,"dataNewSpec");
  if(m_SFD) 
    {
      write(a_handle,m_qbar,"dataBarSpec");
    }
}
/***************************/
/***************************/
void EBAMRSpecies::readCheckpointLevel(HDF5Handle& a_handle)
{
  CH_assert(m_isDefined);
  if (s_verbosity >= 3)
    {
      pout() << "EBAMRSpecies::readCheckpointLevel" << endl;
    }

  // Setup the level string
  char levelStr[20];
  sprintf(levelStr,"%d",m_level);
  const std::string label = std::string("SPlevel_") + levelStr;

  // Read the header for this level
  a_handle.setGroup(label);

  HDF5HeaderData header;
  header.readFromFile(a_handle);

  // Get the refinement ratio
  if (header.m_int.find("ref_ratioS") == header.m_int.end())
    {
      MayDay::Error("::readCheckpointLevel: file does not contain ref_ratio");
    }
  m_ref_ratio = header.m_int["ref_ratioS"];

  //reasons for deviations from non-eb stuff
  // the tag buffer size is set by the factory.
  // dx is set in the define function
  // the problem domain is set in the define function

  // Get dt
  if (header.m_real.find("dtS") == header.m_real.end())
    {
      MayDay::Error("readCheckpointLevel: file does not contain dt");
    }
  m_dt = header.m_real["dtS"];

  // Get the grids
  Vector<Box> vboxGrids;
  const int gridStatus = read(a_handle, vboxGrids);
  if (gridStatus != 0)
    {
      MayDay::Error("readCheckpointLevel: file has no grids");
    }

  Vector<int> proc_map;
  if (s_isLoadBalanceSet)
    {
      s_loadBalance(proc_map,vboxGrids, m_domainBox, false);
    }
  else
    {
      LoadBalance(proc_map,vboxGrids);
    }
  broadcast(proc_map, uniqueProc(SerialTask::compute));

  m_grids= DisjointBoxLayout(vboxGrids,proc_map);

  //this keeps the order of the AMRLevel m_level_grids
  //consistent with m_grids
  LayoutIterator lit = m_grids.layoutIterator();
  for (lit.begin(); lit.ok(); ++lit)
    {
      Box b = m_grids.get(lit());
      m_level_grids.push_back(b);
    }

  //the ebisl has to know about the fact that we really have
  //four ghost cells.
  int nGhostEBISL = 6;
  const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
  ebisPtr->fillEBISLayout(m_ebisl, m_grids,
                          m_domainBox, nGhostEBISL);
  EBCellFactory factoryNew(m_ebisl);
  //m_nghost is set in define function
  IntVect ivGhost = m_nGhost*IntVect::Unit;
  m_stateNew.define(m_grids,m_nComp, ivGhost, factoryNew);
  m_stateOld.define(m_grids,m_nComp, ivGhost, factoryNew);
  if(m_SFD) m_qbar.define(m_grids,m_nComp, ivGhost, factoryNew);

  readCheckpointData(a_handle);

}
void EBAMRSpecies::readCheckpointData(HDF5Handle& a_handle)
{
  // Setup the level string
  char levelStr[20];
  sprintf(levelStr,"%d",m_level);
  const std::string label = std::string("SPlevel_") + levelStr;

  // Read the header for this level
  a_handle.setGroup(label);

  //  Interval vars(0, m_nComp-1);
  //the false says to not redefine the data
  int dataStatusNew = read<EBCellFAB>(a_handle,
                                      m_stateNew,
                                      "dataNewSpec",
                                      m_grids,
                                      Interval(),
                                      false);

  int dataStatusOld = read<EBCellFAB>(a_handle,
                                      m_stateOld,
                                      "dataOldSpec",
                                      m_grids,
                                      Interval(),
                                      false);

  if ((dataStatusNew != 0) || (dataStatusOld != 0))
    {
      MayDay::Error("file does not contain state data");
    }
  // Set up data structures
  levelSetup(); //@readCheckpointLevel

  if(m_SFD) //read restart data
    {    
      int dataBar = 1;
      dataBar = read<EBCellFAB>(a_handle, m_qbar, "dataBarSpec", m_grids, Interval(), false);
      //if no SFD was present (not sure this works, it might issue a stop from read)
      if(dataBar!= 0) EBcopyLevelData(m_qbar, m_stateNew, m_qbar.interval(), m_ebisl);
      
      Interval consInterv(0, m_nComp-1);
      if (m_hasCoarser)
	{
	  const EBAMRSpecies* amrGodCoarserPtr = getCoarserLevel();
	  int refRatCrse = amrGodCoarserPtr->refRatio();
	  EBPWLFillPatch patcher(m_grids,
				 amrGodCoarserPtr->m_grids,
				 m_ebisl,
				 amrGodCoarserPtr->m_ebisl,
				 amrGodCoarserPtr->m_domainBox,
				 refRatCrse, m_nComp, m_nGhost);
      
	  Real coarTimeOld = 0.0;
	  Real coarTimeNew = 1.0;
	  Real fineTime    = 1.0;
	  patcher.interpolate(m_qbar,
			      amrGodCoarserPtr->m_qbar,
			      amrGodCoarserPtr->m_qbar,
			      coarTimeOld,
			      coarTimeNew,
			      fineTime,
			      consInterv);
	}
      m_qbar.exchange(consInterv);
    }
}
/***************************/
void EBAMRSpecies::writePlotHeader(HDF5Handle& a_handle) const
{
  if (s_NewPlotFile == 0)
    {
      writePlotHeaderOld(a_handle);
      return;
    }
  Vector<int> refRatios;
  const EBAMRSpecies* current = this;
  int nlevs = 0;
  while (current != NULL)
  {
    refRatios.push_back(current->refRatio());
    nlevs++;
    current = (const EBAMRSpecies*)(current-> m_finer_level_ptr);
  }

  headerEBFile(a_handle, nlevs, refRatios,
               m_problem_domain.domainBox(),
               m_origin, m_dx, m_aspect,
               m_stateNew.ghostVect());

  writeCellCenteredNames(a_handle, m_stateNames);
  a_handle.setGroup("/Expressions");
  HDF5HeaderData expressions;
  m_ebPatchSpecies->expressions(expressions);
  expressions.writeToFile(a_handle);

}
/***************************/
void EBAMRSpecies::writePlotHeaderOld(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
    {
      pout() << "EBAMRSpecies::writePlotHeader" << endl;
    }

  HDF5HeaderData header;
  int nCons = m_nComp;
  int consAndPrim = nCons ;
  int nCompTotal = nCons;

  Vector<string> names(nCompTotal);

  for (int i = 0; i < nCons; i++)
    {
      names[i] = m_stateNames[i];
    }

  //now output this into the hdf5 handle
  header.m_int["num_components"] = nCompTotal;
  // Setup the component names
  char compStr[30];
  for (int comp = 0; comp < nCompTotal; ++comp)
    {
      sprintf(compStr,"component_%d",comp);
      header.m_string[compStr] = names[comp];
    }

  // Write the header to the file
  header.writeToFile(a_handle);

  if (s_verbosity >= 4)
    {
      pout() << header << endl;
    }
}

void EBAMRSpecies::writePlotLevel(HDF5Handle& a_handle) const
{
  if (s_NewPlotFile == 0)
    {
      writePlotLevelOld(a_handle);
      return;
    }
  if (s_verbosity >= 3)
    {
      pout() << "EBAMRSpecies::writePlotLevel " << m_level<< endl;
    }
  writeCellCentered(a_handle, m_level, &m_stateNew);
}

/***************************/
void EBAMRSpecies::writePlotLevelOld(HDF5Handle& a_handle) const
{

  if (s_verbosity >= 3)
    {
      pout() << "EBAMRSpecies::writePlotLevel" << endl;
    }
  //fill output data including a bunch of geometric crud
  int nCons = m_nComp;
  int consAndPrim = nCons;
  int nCompTotal = nCons;

  Vector<Real> coveredValuesCons(nCons, -1.0);

  Vector<Real> coveredValues;
  coveredValues.append(coveredValuesCons);

  LevelData<FArrayBox> fabData(m_grids, nCompTotal, IntVect::Zero);

#ifdef CH_MPI
    MPI_Barrier(Chombo_MPI::comm);
#endif
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      const EBISBox& ebisbox = m_ebisl[dit()];
      const Box& grid = m_grids.get(dit());
      const EBCellFAB& consfab = m_stateNew[dit()];
      //cfivs, time and timestep fake and not used here
      const Real faket = 1.0;
      const IntVectSet emptyivs;


      FArrayBox& currentFab = fabData[dit()];

      // copy regular data
      currentFab.copy(consfab.getSingleValuedFAB(),0,0,nCons);


      // set special values
      // iterate through the current grid
      // NOTE:  this is probably an inefficient way to do this
      for (BoxIterator bit(grid); bit.ok(); ++bit)
        {
          const IntVect& iv = bit();
          // set special values for covered cells
          if (ebisbox.isCovered(iv))
            {
              for (int icomp = 0; icomp < consAndPrim; icomp++)
                {
                  Real cval = coveredValues[icomp];

                  currentFab(iv,icomp) = cval;
                }
            }
        }//end loop over cells
    }//end loop over grids

#ifdef CH_MPI
    MPI_Barrier(Chombo_MPI::comm);
#endif
  // Setup the level string
  char levelStr[20];
  sprintf(levelStr,"%d",m_level);
  const std::string label = std::string("level_") + levelStr;

  a_handle.setGroup(label);

  // Setup the level header information
  HDF5HeaderData header;

  header.m_int ["ref_ratio"]   = m_ref_ratio;
  header.m_real["dx"]          = m_dx[0];
  header.m_real["dt"]          = m_dt;
  header.m_real["time"]        = m_time;
  header.m_box ["prob_domain"] = m_problem_domain.domainBox();

  // Write the header for this level
  header.writeToFile(a_handle);

  if (s_verbosity >= 4)
    {
      pout() << header << endl;
    }

  // Write the data for this level
  write(a_handle,fabData.boxLayout());
  write(a_handle,fabData,"data");
}
void EBAMRSpecies::outputHDF(const Vector<LevelData<EBCellFAB>* >& a_Ptr, const string& a_name)
{
  string  filename = a_name;
  std::ostringstream oss;
  oss << m_plot + m_level*100;
  filename += oss.str();
  filename.append(".hdf5");
  Real dumReal =  1.0;

  Vector<EBAMRSpecies*>        hierarchy;
  Vector<DisjointBoxLayout>    grids;
  Vector<EBISLayout>           ebisl;
  Vector<EBLevelGrid>          eblg;
  Vector<int>                  refRat;
  ProblemDomain                lev0Dom;
  Real                         lev0Dx;
  Vector<ProblemDomain>        domains;
  Vector<Real>                 dxs;
  getHierarchyAndGrids(hierarchy, grids, ebisl, eblg, refRat, lev0Dom, lev0Dx, domains, dxs);
  int numLevels = hierarchy.size();
  
  int pvars = a_Ptr[0]->nComp();
  char charstrrhs[100];
  Vector<string> namesrhs(pvars, string("var"));
  if(pvars == SpaceDim) // vector in real space
    {
      namesrhs[0] = "xvar";
      namesrhs[1] = "yvar";
      if(SpaceDim == 3) namesrhs[2] = "zvar";
    }
  else
    {
      for (int idir = 0; idir < pvars; idir++)
	{
	  sprintf(charstrrhs, "var%d",idir+1);
	  namesrhs[idir] = charstrrhs;
	}
    }
  bool replaceCovered = true;
  Vector<Real> coveredValues(pvars, 0.0);
  writeEBHDF5noJunk(filename, grids, a_Ptr, namesrhs,
		    lev0Dom.domainBox(), lev0Dx, dumReal, m_time,
		    refRat, numLevels, replaceCovered, coveredValues);
}
#endif

void EBAMRSpecies::outputIVFAB(const Vector<LevelData<BaseIVFAB<Real> >* >& a_Ptr, const string& a_name, const int& a_base)
{
  string  filename = a_name;
  std::ostringstream oss;
  oss << m_plot;
  filename += oss.str();
  filename.append(".ivfab");
  ofstream  out;


	
  Vector<EBAMRSpecies*>        hierarchy;
  Vector<DisjointBoxLayout>    grids;
  Vector<EBISLayout>           ebisl;
  Vector<EBLevelGrid>          eblg;
  Vector<int>                  refRat;
  ProblemDomain                lev0Dom;
  Real                         lev0Dx;
  Vector<ProblemDomain>        domains;
  Vector<Real>                 dxs;
  getHierarchyAndGrids(hierarchy, grids, ebisl, eblg, refRat, lev0Dom, lev0Dx, domains, dxs);
  int numLevels = hierarchy.size();
	
  for (int ip=0; ip<numProc(); ++ip)
    {
      if(procID() == ip)
	{
	  if(ip==0)
	    out.open(filename.c_str(),ofstream::trunc);
	  else
	    out.open(filename.c_str(),ofstream::app);
	  for (int ilev = Min(a_base,numLevels-1); ilev < numLevels; ilev++)
	    {  
	      Real dx = dxs[ilev];
	      RealVect dxv = dx*RealVect::Unit;
	      int nComp = a_Ptr[ilev]->nComp();
	      for (DataIterator dit = a_Ptr[ilev]->dataIterator(); dit.ok(); ++dit)
		{
		  const EBGraph&  ebgraph = ebisl[ilev][dit()].getEBGraph();
		  for (VoFIterator vofit((*a_Ptr[ilev])[dit()].getIVS(), ebgraph); vofit.ok(); ++vofit)
		    {
		      VolIndex vof = vofit();
		      RealVect xloc =  EBArith::getVofLocation(vof,dxv,RealVect::Zero) ;
		      IntVect  iloc = vof.gridIndex();
		      out <<  iloc[0] ;
		      for(int k=1; k<SpaceDim;k++)out << ", " << iloc[k] ;
		      for(int k=0; k<SpaceDim;k++)out << ", " << xloc[k] ;
		      for(int k=0; k<nComp; k++)  out << ", " << (*a_Ptr[ilev])[dit()](vof,k);
		      out << endl;
			
		    }
		}
	    }
	  out.close();
	}
#ifdef CH_MPI
      MPI_Barrier(Chombo_MPI::comm);
#endif
    }

}
/****************/
EBAMRSpecies::
~EBAMRSpecies()
{
  if (m_ebPatchSpecies != NULL)  delete m_ebPatchSpecies;
}
#include "NamespaceFooter.H"
