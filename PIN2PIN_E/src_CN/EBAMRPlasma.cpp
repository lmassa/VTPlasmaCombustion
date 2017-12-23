#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <math.h>

#include "LevelData.H"
#include "EBCellFAB.H"
#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "EBCellFactory.H"
#include "EBViscousTensorOp.H"
#include "EBConductivityOp.H"
#include "EBLevelDataOps.H"
#include "BiCGStabSolver.H"
#include "EBLevelGrid.H"
#include "EBEllipticLoadBalance.H"


#include "ParmParse.H"
#include "EBAMRPlasma.H"
#include "EBPatchPolytropicF_F.H"
#include "EBArith.H"
#include "REAL.H"
#include "EBSpeciesIBCFactory.H"
#include "EBAMRSpecies.H"
#include "EBAMRSpeciesFactory.H"
#include "EBPatchPolytropic.H"
#include "EBAMRIO.H"
#include "EBArith.H"
#include "EBAMRDataOps.H"
#include "EBViscousTensorOp.H"
#include "EBViscousTensorOpFactory.H"
#include "DirichletPoissonEBBC.H"
#include   "NeumannPoissonEBBC.H"
#include "DirichletPoissonDomainBC.H"
#include   "NeumannPoissonDomainBC.H"
#include "DirichletViscousTensorEBBC.H"
#include   "NeumannViscousTensorEBBC.H"
#include "DirichletViscousTensorDomainBC.H"
#include   "NeumannViscousTensorDomainBC.H"
#include     "MixedViscousTensorDomainBC.H"
#include "DirichletConductivityDomainBC.H"
#include   "NeumannConductivityDomainBC.H"
#include "DirichletConductivityEBBC.H"
#include   "NeumannConductivityEBBC.H"
#include <iomanip>
#include <cmath>
#include <cstdio>
#include "UsingNamespace.H"
#include "memusage.H"
#include "memtrack.H"

/**********************/
EBAMRPlasma::EBAMRPlasma()
{

  CH_TIME("EBAMRPlasma::EBAMRPlasma");

  // Read parameters from file (some are in getPlasmaParameters)
  ParmParse ppgodunov;
  m_numGhost = 4;
  m_PlasmaPhysics.define();
  m_nSpec = m_PlasmaPhysics.nComponents();
  getPlasmaParameters(m_params);

  
  
  // Computational Domain
  ProblemDomain coarsestDomain = m_params.coarsestDomain;
  
  // Physical Domain
  vector<Real> length(SpaceDim);
  ppgodunov.getarr("domain_length",length,0,SpaceDim);
  RealVect domainLength, coarsestDx;
  for (int idir=0;idir<SpaceDim;idir++) domainLength[idir] = length[idir];
  std::vector<int> nCellsArray(SpaceDim);
  ppgodunov.getarr("n_cells",nCellsArray,0,SpaceDim);
  for (int idir = 0; idir < SpaceDim; idir++)
    coarsestDx[idir] = domainLength[idir]/nCellsArray[idir];
  

  // Geometry
  // this should call the EBIS, that's why it is at the beginning
  m_ElectricPotential.defineGeometry(); 



  // get the maxLevel first
  int maxLevel = 0;
  ppgodunov.get("max_level",m_maxRefLevel);
  m_maxRefLevelS = m_maxRefLevel;
  ppgodunov.query("max_level_species",m_maxRefLevelS);
  maxLevel = Max(m_maxRefLevel,m_maxRefLevelS);


  pout() << "Check...3";
  Real yEB = m_ElectricPotential.getYEB();
  m_PlasmaPhysics.setEBy(yEB);
  Real epsg, epsd;
  ppgodunov.get("eps_gas", epsg);
  ppgodunov.get("eps_diel", epsd);
  m_PlasmaPhysics.setEps(epsg, epsd);
//
 

  // 4 steps 
  // (1) define initial and boundary conditions
  // (2) define gamma-gas patch
  // (3) define inviscid solver
  // (4) define viscous solver

  //(1) define IBC
  m_gamma = 1.4;
  ppgodunov.get("gamma",m_gamma);
  Real MachInflow;
  ppgodunov.get("Mach_inflow",MachInflow);
  Real shockCenter;
  ppgodunov.get("shock_center",shockCenter);
  int inormal;
  ppgodunov.get("shock_normal",inormal);
  int ishockback;
  ppgodunov.get("shock_backward",ishockback);
  bool shockbackward = (ishockback == 1);
  //force these for now.
  bool doSmushing = true;
  bool doRZCoords = false;
  bool hasSourceTerm = false;
  bool useLimiting = true;

  pout() << "Check...4";
  m_Tinf = m_PlasmaPhysics.m_Tg; // TG hardwired cold gas Temp in PlasmaPhysics
  Real Rex,Prandtl,ReL,ReCH,pinf;
  ppgodunov.get("Rex",Rex);
  ppgodunov.get("Prandtl",   Prandtl);
  ppgodunov.get("pinf",   pinf);
  Real Runiv = 8.314e3; //J/kmole/K
  Real Rgas = Runiv/m_PlasmaPhysics.m_molWeight; // Joule/kg/K
  // ReL is Re for Lchar = 1m
  ReL = pinf*MachInflow*sqrt(m_gamma/(Rgas*m_Tinf))*(m_Tinf+79.4)/(1.4844e-6*pow(m_Tinf,1.5));
  ReCH = pinf*sqrt(m_gamma/(Rgas*m_Tinf))*(m_Tinf+79.4)/(1.4844e-6*pow(m_Tinf,1.5))/m_gamma;
  //the dimensional velo is the speed of sound at the inflow
  m_Vel0 = sqrt(m_gamma*Rgas*m_Tinf);
  // the dimensional source term is pinf*Vref
  m_Src0 = m_Vel0*pinf;
  m_SrcV = pinf;
  
    
  // Define the Physical IBC Factory;
  m_BL.define(MachInflow, m_gamma, m_Tinf, Prandtl, yEB, Rex, ReCH);
  m_bcGammaFactory.define(m_gamma, MachInflow, shockCenter, inormal, shockbackward, doRZCoords,&m_BL);

  pout() << " ReL = " << ReL << ", ReCH = " << ReCH << ", Rex= " <<  Rex << ", thetaMom= " << m_BL.m_thetaMom << ", X0= " << m_BL.getX0()  << endl;
  pout() << "non dimensional source " << m_Src0 << endl;

  //(2) define gamma-gas patch
  int ifourth, iflatten, iartvisc;
  ppgodunov.get("use_fourth_order_slopes", ifourth);
  ppgodunov.get("use_flattening"         , iflatten);
  ppgodunov.get("use_art_visc"           , iartvisc);
  bool useFourthOrderSlopes = (ifourth  ==1);
  bool useFlattening        = (iflatten ==1);
  bool useArtificialVisc    = (iartvisc ==1);
  

  //create patch integrator
  m_patchGammaFactory.define(&m_bcGammaFactory, m_gamma,
		      useFourthOrderSlopes, useFlattening, useArtificialVisc,
		      useLimiting, doRZCoords);

  //(3) define AMR inviscid object
  //get input line parameters for AMR factory
  // periodic BCs
  bool is_periodic[SpaceDim];
  for (int i = 0; i < SpaceDim; ++i)is_periodic[i] = false;
  int isperiodicz;
  ppgodunov.get("isperiodicz",isperiodicz);
  if(SpaceDim >=3 ) is_periodic[SpaceDim-1] = (isperiodicz == 1);
  Real initialCFL = 0.1;ppgodunov.get("initial_cfl",initialCFL);
  Real cfl = 0.8;ppgodunov.get("cfl",cfl);
  int redistRad = 0;ppgodunov.get("redist_radius",redistRad);
  Real refineThresh = 0.3;ppgodunov.get ("vorticity_refine_thresh",refineThresh);
  int tagBufferSize = 3;ppgodunov.get("vorticity_tag_buffer_size",tagBufferSize);
  int verbosity;ppgodunov.get("verbosity",verbosity);m_verbosity = verbosity;
  int iusemassredist;ppgodunov.get("use_mass_redist", iusemassredist);bool useMassRedist = (iusemassredist ==1);

  // setup problem domain
  bool symmetric = true;
  m_aco.resize(        maxLevel+1);
  m_eta.resize(        maxLevel+1);
  m_lambda.resize(     maxLevel+1);
  m_etaIrreg.resize(   maxLevel+1);
  m_lambdaIrreg.resize(maxLevel+1);
  m_rcv.resize(        maxLevel+1);
  EBAMRNavierFactory amrg_fact(initialCFL, cfl, redistRad,
			       domainLength, refineThresh, tagBufferSize,
			       verbosity, useMassRedist, doSmushing,
			       doRZCoords, hasSourceTerm, &m_patchGammaFactory, symmetric,
			       &m_aco, &m_eta, &m_lambda, &m_etaIrreg, &m_lambdaIrreg,
			       &m_rcv, false);

  //create AMR
  int num_read_levels = Max(maxLevel,1);
  std::vector<int> ref_ratios; // (num_read_levels,1);
  // note this requires a ref_ratio to be defined for the
  // finest level (even though it will never be used)
  ppgodunov.getarr("ref_ratio",ref_ratios,0,num_read_levels+1);
  
  m_amr.define(m_maxRefLevel,ref_ratios,coarsestDomain,&amrg_fact);

  // rather lengthy procedure to set-up the computation
  Real fixedDt = -1;
  ppgodunov.get("fixed_dt",fixedDt);
  bool isFixedDt = fixedDt>0;
  m_inviscid = false;
  ppgodunov.query("inviscid",m_inviscid);
  m_tagShrinkDomain =0;
  ppgodunov.query("tag_shrink_domain",m_tagShrinkDomain);
  m_tagVortComps = 1;  // tag only streamwise vorticity
  ppgodunov.query("tag_vorticity_comp",m_tagVortComps);
  m_tagVortComps = Min(m_tagVortComps,SpaceDim);
  m_RegridSteps = -1;
  ppgodunov.query("plasma_regrid_steps",m_RegridSteps);
  m_gridRedefineSteps = 1;
  ppgodunov.query("grid_redefine_steps",m_gridRedefineSteps);
  m_viscousOutput = false;
  ppgodunov.query("viscous_output",m_viscousOutput);

  int maxGridSize = 32;
  ppgodunov.get("max_grid_size",maxGridSize);
  int blockFactor = 1;
  ppgodunov.get("block_factor",blockFactor);
  Real fillRatio = 0.75;
  ppgodunov.get("fill_ratio",fillRatio);
  int gridBufferSize;
  ppgodunov.get("grid_buffer_size",gridBufferSize);


  if(isFixedDt) m_amr.fixedDt(fixedDt);
// set grid generation parameters
  m_amr.maxGridSize(maxGridSize);
  m_amr.blockFactor(blockFactor);
  m_amr.fillRatio(fillRatio);
  m_amr.gridBufferSize(gridBufferSize);

  m_checkpointInterval = 0;ppgodunov.get("checkpoint_interval",m_checkpointInterval);
  m_plot_interval = 0;
  ppgodunov.get("plot_interval",m_plot_interval);

  Real maxDtGrowth = 1.1;ppgodunov.get("max_dt_growth",maxDtGrowth);
  Real dtToleranceFactor = 1.1;ppgodunov.get("dt_tolerance_factor",dtToleranceFactor);


  //regrid Intervals (num_read_levels defined above)
  std::vector<int> regridIntervals; // (num_read_levels,1);
  ppgodunov.getarr("regrid_interval",regridIntervals,0,num_read_levels);

  // set output parameters in AMR
  m_amr.checkpointInterval(-1); //disable
  m_amr.plotInterval(m_plot_interval);
  m_amr.regridIntervals(regridIntervals);
  m_amr.maxDtGrow(maxDtGrowth);
  m_amr.dtToleranceFactor(dtToleranceFactor);
  m_amr.verbosity(m_verbosity);

  bool useSubcycling=true;
  if (ppgodunov.contains("use_subcycling"))
    {
      ppgodunov.get("use_subcycling", useSubcycling);
      if (!useSubcycling)
        {
          pout() << "SUBCYCLING IN TIME TURNED OFF!!!"  << endl;
        }
      m_amr.useSubcyclingInTime(useSubcycling);
    }
  std::string prefix;
  if (ppgodunov.contains("plot_prefix"))
    {
      ppgodunov.get("plot_prefix",prefix);
      m_amr.plotPrefix(prefix);
    }

  std::string sprefix;
  if (ppgodunov.contains("chk_prefix"))
    {
      ppgodunov.get("chk_prefix",sprefix);
      m_amr.checkpointPrefix(sprefix);
    }

  std::string restart_file;
  bool isRestart = false;
  if (!ppgodunov.contains("restart_file"))
    {
      if (!ppgodunov.contains("fixed_hierarchy"))
        {
          // initialize from scratch for AMR run
          // initialize hierarchy of levels
	  pout() << ">>>>>>>>>>>>>>>>>>>>>>" << endl;
          m_amr.setupForNewAMRRun();
	  pout() << "<<<<<<<<<<<<<<<<<<<" << endl;
        }
      else
        {
          MayDay::Error("fixed grid not implemented in EBAMRPlasma");
        }
    }
  else
    {
      ppgodunov.get("restart_file",restart_file);
      pout() << " restarting from file " << restart_file << endl;

#ifdef CH_USE_HDF5
      HDF5Handle handle(restart_file,HDF5Handle::OPEN_RDONLY);
      // read from checkpoint file
      m_amr.setupForRestart(handle);
      handle.close();
      isRestart = true;
#else
      MayDay::Error("amrGodunov restart only defined with hdf5");
#endif
    }


  //first create the IBC object
  m_bcSpeciesFactory.define(m_nSpec, &m_PlasmaPhysics);//calls create  from EBPatchGodunov::setEBPhysIBC
  //create patch integrator
  bool spec_use_limiting = true;
  bool spec_setMaxMin = true;
  Real spec_max_val = 1e10;
  Real spec_min_val = 0;   
  refineThresh = 5e-11;ppgodunov.query ("electron_refine_thresh",refineThresh);
  ppgodunov.query("species_tag_buffer_size",tagBufferSize);
  m_patchSpeciesFactory.define(&m_bcSpeciesFactory, &m_PlasmaPhysics, spec_use_limiting,
			       m_nSpec,spec_max_val,spec_min_val,spec_setMaxMin);
  
  bool spec_hasSourceTerm = true;
  EBAMRSpeciesFactory amrg_spec_fact(initialCFL, cfl, redistRad,
				     domainLength, refineThresh, tagBufferSize,
				     verbosity, useLimiting, doSmushing,
				     doRZCoords, spec_hasSourceTerm, 
				     &m_patchSpeciesFactory);
  m_amrS.define(m_maxRefLevelS,ref_ratios,coarsestDomain,&amrg_spec_fact);

  //turn-off internal regridding; this should be fine for few levels runs
  std::vector<int> regridIntervalsSpec(num_read_levels,-1);

  if(isFixedDt) m_amrS.fixedDt(fixedDt);
// set grid generation parameters
  m_amrS.maxGridSize(maxGridSize);
  m_amrS.blockFactor(blockFactor);
  m_amrS.fillRatio(fillRatio);
  m_amrS.gridBufferSize(gridBufferSize);
  m_amrS.checkpointInterval(-1); // disable
  m_amrS.plotInterval(m_plot_interval);
  m_amrS.regridIntervals(regridIntervalsSpec);
  m_amrS.maxDtGrow(maxDtGrowth);
  m_amrS.dtToleranceFactor(dtToleranceFactor);
  m_amrS.verbosity(m_verbosity);

  // do the restart on top of the intialization not to mess with the rigid AMR structure
  //move this block up here. Stopped working on the rhos, It should be dimensionaed before being restarted
  if(!isRestart)
    {
      pout() << "&&&&&&&&&&&&&&" << endl;
      m_amrS.setupForNewAMRRun();
      pout() << "####################" << endl;
    }
  else
    {
      HDF5Handle handle(restart_file,HDF5Handle::OPEN_RDONLY); 
      m_amrS.setupForRestartS(handle);
      handle.close(); 
      m_amrS.initialTime(currentTime());
      m_amrS.setStep(m_amr.getStep());
      m_amrS.setNewDt(m_amr.getOldDt()/m_Vel0);
    }
  m_amrS.useSubcyclingInTime(useSubcycling);
  m_amrS.plotPrefix(prefix);
  m_amrS.checkpointPrefix(sprefix);



  // (4) Define  Grids and level pointers
  definegrids();

  
  EBPatchPolytropic Ppatch;
  m_nCompState = Ppatch.numConserved();
  //define the vel and temp pointers
  m_Velo.resize(m_params.numLevels,NULL);
  m_Temp.resize(m_params.numLevels,NULL);
  m_FSrc.resize(m_params.numLevels,NULL);
  m_cellNData.resize(m_params.numLevels);
  m_rhos.resize(m_params.numLevels,NULL);
  for (int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
      EBCellFactory factory(m_ebisl[ilev]);
      m_Velo[ilev]   = new LevelData<EBCellFAB>();
      m_Temp[ilev]   = new LevelData<EBCellFAB>();
      m_FSrc[ilev]   = new LevelData<EBCellFAB>();
      m_rhos[ilev]   = new LevelData<BaseIVFAB<Real> >();
      m_Velo[ilev]->define(m_grids[ilev], SpaceDim, m_params.ghostPhi, factory);
      m_Temp[ilev]->define(m_grids[ilev], 2       , m_params.ghostPhi, factory);
      m_FSrc[ilev]->define(m_grids[ilev], m_nCompState    , m_params.ghostPhi, factory);
      LayoutData<IntVectSet> irregSets(m_grids[ilev]);
      m_ghost_Dk = 2;// m_ghost_Dk = 2;
      for(int k=0;k<SpaceDim;k++)m_ghost_Dk=Max(m_ghost_Dk,m_params.ghostPhi[k]);
      for (DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit) 
	{
	  Box grownBox = grow(m_grids[ilev].get(dit()), m_ghost_Dk); 
	  grownBox &= m_domain[ilev];
	  irregSets[dit()] = m_ebisl[ilev][dit()].getIrregIVS(grownBox);
	}
      BaseIVFactory<Real>  baseivfact(m_ebisl[ilev], irregSets);
      m_rhos[ilev]->define(m_grids[ilev], 1    , m_ghost_Dk*IntVect::Unit, baseivfact);
      for (DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit) 
	(*m_rhos[ilev])[dit()].setVal(0e0);
    }


  //get the velocity only becasue the electron temeprature data have not yet been defined
  EBAMRDataOps::setVal(m_Velo,1.0);
  EBAMRDataOps::setVal(m_Temp,1.0);
  EBAMRDataOps::setVal(m_FSrc,0.0);
  EBAMRDataOps::setCoveredVal(m_FSrc,0.0);
  m_NavierLevels[0]->getVelTemp(m_Velo,m_Temp,m_grids,m_ebFineInterps);
  EBAMRDataOps::scale(m_Velo,m_Vel0);
  EBAMRDataOps::scale(m_Temp,m_Tinf,0);// scale only gas temperature
  updateGhost(m_Velo,10);
  updateGhost(m_Temp,10);


  //define the potential
  if(!isRestart)
    m_ElectricPotential.setupForGivenDBLRun(&m_PlasmaPhysics, m_grids);
  else
    {
      HDF5Handle handle(restart_file,HDF5Handle::OPEN_RDONLY);
      // read from checkpoint file for phi
      m_ElectricPotential.setupForRestart(&m_PlasmaPhysics, m_grids, handle); 
      //restart rhos
      for (int level = 0; level <= m_finestLevel; ++level)
	{
	  // Setup the level string
	  char levelStr[20];
	  sprintf(levelStr,"%d",level);
	  const std::string label = std::string("SPlevel_") + levelStr;
	  handle.setGroup(label);
	  int dataRhos = read<BaseIVFAB<Real> >(handle,*m_rhos[level],"rhos",m_grids[level],Interval(),false);
	  if (dataRhos != 0)  MayDay::Error("rhosfile does not contain state data");
	}
      handle.close();
    }    
  for (int ilev = 0; ilev <= m_finestLevel; ++ilev)m_rhos[ilev]->exchange(Interval(0,0));
  m_phi = m_ElectricPotential.getPhi();
  m_gradPhi = m_ElectricPotential.getGradPhi();
  m_cellGradPhi = m_ElectricPotential.getCellGradPhi();
  m_IVGradPhi = m_ElectricPotential.getIVGradPhi();

  // register the pointers with plasma
  for (int ilev = 0; ilev <= m_finestLevel; ilev++) 
    m_SpeciesLevels[ilev]->SetVelTempPhi(m_Velo,m_Temp,m_phi,m_gradPhi,m_cellGradPhi,m_IVGradPhi);

    

  // get the electron energy and species
  for (int ilev = 0; ilev <= m_amrS.getFinestLevel(); ilev++)
      m_SpeciesLevels[ilev]->getElectronEnergy(*m_Temp[ilev],1);
  updateGhost(m_Temp,10);
  getSpecies();
  

  m_NavierLevels[0]->PlotVelocity(1, 1, false);
  
  // intialize the electrostatic potential to the initial species
  //Solve Electric Poisson Equation
  m_ElectricPotential.setTime(currentTime());
  if(! isRestart)
    {
      m_ElectricPotential.solve(m_cellNData);
    }
  m_phi = m_ElectricPotential.getPhi();
  m_gradPhi = m_ElectricPotential.getGradPhi();
  m_cellGradPhi = m_ElectricPotential.getCellGradPhi();
  getVelTemp();
  for (int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
      m_SpeciesLevels[ilev]->SetVelTempPhi(m_Velo,m_Temp,m_phi);
    }
  m_NavierLevels[m_finestLevelNS]-> defineSolvers();
  m_SpeciesLevels[m_finestLevel]->defineSolvers();
   
  bool isNAN2 = EBAMRDataOps::checkNANINF(m_cellGradPhi);
  if(isNAN2)MayDay::Error("EBAMRPlasma::(Constructor) checkNANINF m_cellGradPhi."); 

  //Regreid the initial configuration
  bool doInitialRegrid=false;ppgodunov.query("do_initial_regrid",doInitialRegrid);
  if(isRestart && (m_amrS.getFinestLevel() < maxLevel || doInitialRegrid))
    {
      m_amr.regrid(0);
      regrid();
    }
  else
    printGridSize();

  //refinethresh for Source
  m_SrefineThresh = 2.5e2;
  ppgodunov.query ("source_refine_thresh",m_SrefineThresh);

  stdOut(0);
  
  m_SpeciesLevels[0]->PlotReactiveSource(0);//lucacancel
  //m_SpeciesLevels[0]->m_plot = m_amr.getFinestLevel();
  //m_SpeciesLevels[0]->outputIVFAB(m_rhos, string("rhos"), 10000);

  pout () << "Ending Plasma Constructor" << endl;
}
//**************
void EBAMRPlasma::run(Real a_max_time, int a_max_step, const char* a_inFile)
{

  EBPatchPolytropic Ppatch;
  Interval momIntrv = Ppatch.momentumInterval();
  int enrgIndx = Ppatch.energyIndexC();

  //setup the plot frequency
  Real frequency, omega, nplot_per_period, dtp, ndtp, dtstep;
  ParmParse pp;	 
  pp.get("frequency", frequency);
  pp.get("nplot_per_period", nplot_per_period);
  omega = 2*M_PI*frequency;
  if(omega > 0)
    dtp = 8.0*atan(1.0)/(omega*nplot_per_period);
  else
    dtp = 8.0*atan(1.0)/(5e4*nplot_per_period);

  // to perform only fluid simulations
  bool onlyfluid = false;
  pp.query("only_fluid", onlyfluid);
  bool plotVorticity = false;
  pp.query("plot_vorticity", plotVorticity);  
  int nblownup=0;
  Real EPresid=0; //the residual of the Electric potentail equation

  int nSpecSteps=1;
  pp.query("species_to_fluid_steps", nSpecSteps);
  if(nSpecSteps > 1) pout () << "Using " << nSpecSteps << " per each fluid step." << endl;
  
  
  for (int steps = 0 ; m_amr.ok(a_max_time,a_max_step); m_amr.add(),m_amrS.add(),++steps)
    {          
      
      //CheckPoint
      if (steps > 0 && m_checkpointInterval > 0 && (m_amr.getStep() % m_checkpointInterval) == 0) writeCheckpointFile(a_inFile);

      // preps
      m_phi = m_ElectricPotential.getPhi();
      m_gradPhi = m_ElectricPotential.getGradPhi();
      m_cellGradPhi = m_ElectricPotential.getCellGradPhi();
      m_IVGradPhi = m_ElectricPotential.getIVGradPhi();
      // usng the updated temperature (mu) redefine the NS solvers (not strictly necessary, the temperature changes slowly and this is called by regrid)
      m_NavierLevels[m_finestLevelNS]->defineSolvers();

      // fluid source
      for (int ilev = 0; ilev <= m_finestLevel; ilev++)
	{
	  if(!onlyfluid)
	    {
	      m_SpeciesLevels[ilev]->makeFluidSource(*m_FSrc[ilev]);
	      //scale by ,1/srcV 1/src0 to get non-dimensional value
	      EBLevelDataOps::scale(*m_FSrc[ilev],1.0/m_Src0,enrgIndx);
	      for (int k = momIntrv.begin(); k <= momIntrv.end(); k++)
		EBLevelDataOps::scale(*m_FSrc[ilev],1.0/m_SrcV,k);
	    }
	}      
      if(!onlyfluid) EBAMRDataOps::averageDown(m_FSrc,m_eblg,m_params.refRatio);
      m_NavierLevels[0]->setSource(m_FSrc);

      
      //NS solution
      m_amr.oneStep();
      definegrids();
      getVelTemp();

      //Synchronize
      Real dtSpecies = m_amr.getOldDt()/m_Vel0;
      Real tSpecies = m_amr.getCurrentTime()/m_Vel0;
      m_amrS.setNewDt(dtSpecies/nSpecSteps);
      m_ElectricPotential.setTime(tSpecies);

      if(!onlyfluid)
	{
	  //Plasma Solution
	  for (int ilev = 0; ilev <= m_finestLevel; ilev++) 
	    m_SpeciesLevels[ilev]->SetVelTempPhi(m_Velo,m_Temp,m_phi,m_gradPhi,m_cellGradPhi,m_IVGradPhi); 
	  m_SpeciesLevels[m_finestLevel]->defineSolvers();
	  for (int istep = 1; istep <= nSpecSteps; istep++) m_amrS.oneStep();
	  getSpecies();
	  
	  //Charge density
	  for (int ilev = 0; ilev <= m_finestLevel; ilev++) 
	    {
	      if(m_verbosity > 3)pout() << ilev << " Updating rhos via Plasma Physics " << endl;
	      m_PlasmaPhysics.rhos(*m_rhos[ilev], *m_rhosRHS[ilev], m_grids[ilev], m_ebisl[ilev], m_domain[ilev], 2, dtSpecies);
	    }
	  
	  //Solve Electric Poisson Equation
	  if ( (m_amr.getStep() % 1) == 0)
	    {
	      m_ElectricPotential.updateRhos(m_rhos);
	      EPresid=m_ElectricPotential.solve(m_cellNData,true);
	    }
	}

      stdOut(EPresid);
     
      // plot the species only solution space
      ndtp = tSpecies/dtp;
      dtstep = (ndtp - floor(ndtp))*dtp;
      if (Min(tSpecies, dtstep) <= dtSpecies)
	{
	  int nplot = floor(ndtp + 0.5);
	  m_SpeciesLevels[0]->PlotReactiveSource(nplot, dtSpecies,onlyfluid);
	  if(plotVorticity) m_NavierLevels[0]->PlotVorticity(nplot);
	  bool removeSelfSimilar = true;
	  m_NavierLevels[0]->PlotVelocity(nplot, m_Vel0, removeSelfSimilar);//m_NavierLevels[0]->PlotVelocity(nplot);//
	  //if (!onlyfluid) m_SpeciesLevels[0]->outputIVFAB(m_rhos, string("rhos"), 10000);
	}

      //Regrid
      if (!onlyfluid && m_amr.getStep() > 0 && m_RegridSteps >0 && (m_amr.getStep() % m_RegridSteps) == 0) regrid();
      
    }


}

//**************
void EBAMRPlasma::runPlasmaOnly(Real a_max_time, int a_max_step, const char* a_inFile)
{


  EBAMRSpecies::s_isOnlyPlasma = true;
  Real EPresid=0;
  Real dtSpecies = m_SpeciesLevels[m_finestLevel]->setImplicitTimeStep();//CheckA
  for (int steps = 0 ; m_amrS.ok(a_max_time,a_max_step); m_amrS.add(),++steps)
    {          
      
      // preps
      grabPotentials();
      definegrids();
      getVelTemp(false);


      //Plasma Solution
      m_amrS.setNewDt(dtSpecies);
      for (int ilev = 0; ilev <= m_finestLevel; ilev++) m_SpeciesLevels[ilev]->ZeroVelSetTempPhi(m_Velo,m_Temp, m_phi,m_gradPhi,m_cellGradPhi,m_IVGradPhi); 
      m_SpeciesLevels[m_finestLevel]->defineSolvers();
      m_amrS.oneStep();

      
      Real tSpecies = m_amrS.getUpdatedTime();
      m_ElectricPotential.setTime(tSpecies);
      getSpecies();  // grabs  m_cellNData
      EPresid=m_ElectricPotential.solve(m_cellNData,true);

	  
      stdOut(EPresid);
      
      if((steps % 25) == 0) m_SpeciesLevels[0]->PlotReactiveSource(steps, dtSpecies,false);

      if ( m_RegridSteps >0 && (m_amrS.getStep() % m_RegridSteps) == 0) regrid(true);

      
    }

  pout() << "Done in EBAMRPlasma::runPlasmaOnly" << endl;


}
/*******************************/
void EBAMRPlasma::definegrids()
{
  CH_TIME("EBAMRPlasma::definegrids");


  //the assumptions in this "Grids" versio of the code is that the max number of levels in the fluids and species grids are the same
  m_finestLevel = m_amrS.getFinestLevel();
  m_finestLevelNS = m_amr.getFinestLevel();
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  Vector<AMRLevel*> AMRLevels =   m_amr.getAMRLevels();
  Vector<AMRLevel*> AMRLevelsS = m_amrS.getAMRLevels();
  int AMRsize = AMRLevelsS.size();
  m_params.numLevels = Min(m_finestLevel+1,AMRsize);

  int numLevels = m_params.numLevels;
  m_NavierLevels.resize(numLevels);
  m_SpeciesLevels.resize(numLevels);
  m_grids.resize(numLevels);
  m_ebisl.resize(numLevels);
  m_eblg.resize(numLevels);
  m_domain.resize(m_params.maxLevel+1); //slightly larger
  m_dx.resize(m_params.maxLevel+1);
  m_ebFineInterps.resize(numLevels);


  m_domain[0] = m_params.coarsestDomain;
  m_dx[0]     = m_params.coarsestDx[0];
  for (int ilev = 1; ilev <= m_params.maxLevel; ++ilev)
    {
      m_domain[ilev] = refine(m_domain[ilev-1],m_params.refRatio[ilev-1]);
      m_dx[ilev] = m_dx[ilev-1]/Real(m_params.refRatio[ilev-1]);
    }

  // get the DBLS and EBIS
  for (int ilev = 0; ilev < numLevels; ++ilev)
    {
      if(ilev <= m_finestLevelNS) m_NavierLevels[ilev] =  (static_cast <EBAMRNavier*> (AMRLevels[ilev]));
      m_SpeciesLevels[ilev] = (static_cast <EBAMRSpecies*> (AMRLevelsS[ilev]));
      m_grids[ilev] = m_SpeciesLevels[ilev]->getStateNew().disjointBoxLayout();
      m_ebFineInterps[ilev] = &(m_SpeciesLevels[ilev]->m_ebFineInterp);
      ebisPtr->fillEBISLayout(m_ebisl[ilev], m_grids[ilev], m_domain[ilev], m_numGhost);
      m_eblg[ilev] = EBLevelGrid(m_grids[ilev], m_ebisl[ilev], m_domain[ilev]);
    }
}
//*******************************************
void EBAMRPlasma::regrid(const bool a_isPlasmaOnly)
{
  int numLevels = m_params.numLevels;

  int top_level = Min(m_amrS.getFinestLevel(), m_maxRefLevelS-1);
  if (top_level>=0)
    {
      int finest_level_old = m_amrS.getFinestLevel();
      Vector<LevelData<EBCellFAB>* > rhosData;
      Vector<DisjointBoxLayout> oldDBL;
      if(!a_isPlasmaOnly) regridRhos(rhosData,oldDBL,0); //for now plasma-only without surface accumulation
      Vector<IntVectSet> tagsVect(top_level+1);
      // here we have three modules together NS, Species, Electric
      // tag each individually and take the union of the tags
      for (int ilev = 0; ilev <= top_level; ilev++)
	{
	  IntVectSet SPtags, EFtags, SRCtags;
	  if(a_isPlasmaOnly)
	    {
	      m_SpeciesLevels[ilev]->tagCells(tagsVect[ilev]);
	    }
	  else
	    {
	      tagCells(SRCtags,ilev);
	      m_SpeciesLevels[ilev]->tagCells(SPtags);
	      if(!m_ElectricPotential.doTagging())
		tagsVect[ilev] = SRCtags | SPtags;
	      else
		{
		  m_ElectricPotential.tagCells(EFtags,ilev);
		  tagsVect[ilev] = SRCtags | SPtags | EFtags;
		}
	    }
	  //if(ilev == top_level) pout() << "EBAMRPlasma::regrid: tags[" << ilev << "]: " << SRCtags << endl;
	}
      // after the section below the 3 modules should have the same grid
      m_amrS.regrid(tagsVect, 0);
      definegrids(); //updats m_params.numLevels
      //CH_assert(m_amrS.getFinestLevel() == m_amr.getFinestLevel()); //not necessary, just to start
      m_Velo.resize(m_params.numLevels,NULL);
      m_Temp.resize(m_params.numLevels,NULL);
      m_FSrc.resize(m_params.numLevels,NULL);
      m_cellNData.resize(m_params.numLevels);
      m_rhos.resize(m_params.numLevels,NULL);
      for (int ilev = 0; ilev <= m_finestLevel; ilev++)
	{
	  EBCellFactory factory(m_ebisl[ilev]);	  
	  if(m_Velo[ilev]==NULL)m_Velo[ilev]   = new LevelData<EBCellFAB>();
	  if(m_Temp[ilev]==NULL)m_Temp[ilev]   = new LevelData<EBCellFAB>();
	  if(m_FSrc[ilev]==NULL)m_FSrc[ilev]   = new LevelData<EBCellFAB>();
	  m_Velo[ilev]->define(m_grids[ilev], SpaceDim, m_params.ghostPhi, factory);
	  m_Temp[ilev]->define(m_grids[ilev], 2       , m_params.ghostPhi, factory);
	  m_FSrc[ilev]->define(m_grids[ilev],m_nCompState, m_params.ghostPhi, factory);
	}
      EBAMRDataOps::setVal(m_FSrc,0.0);
      EBAMRDataOps::setCoveredVal(m_FSrc,0.0);

      m_ElectricPotential.regridWithDBL(m_grids);
      if(!a_isPlasmaOnly) regridRhos(rhosData,oldDBL,1);
      if(!a_isPlasmaOnly) regridRhos(rhosData,oldDBL,2);


      getVelTemp();  //@regrid
      m_phi = m_ElectricPotential.getPhi();
      m_gradPhi = m_ElectricPotential.getGradPhi();
      m_cellGradPhi = m_ElectricPotential.getCellGradPhi();
      m_IVGradPhi = m_ElectricPotential.getIVGradPhi();
      if(a_isPlasmaOnly)
	for (int ilev = 0; ilev <= m_finestLevel; ilev++)
	  m_SpeciesLevels[ilev]->ZeroVelSetTempPhi(m_Velo,m_Temp, m_phi,m_gradPhi,m_cellGradPhi,m_IVGradPhi);
      else
	for (int ilev = 0; ilev <= m_finestLevel; ilev++)
	  m_SpeciesLevels[ilev]->SetVelTempPhi    (m_Velo,m_Temp,m_phi,m_gradPhi,m_cellGradPhi,m_IVGradPhi); 
      m_SpeciesLevels[m_finestLevel]->defineSolvers();      
      getSpecies();  //@EBAMRPlasma::regrid

      if(!a_isPlasmaOnly) m_ElectricPotential.updateRhos(m_rhos);
      Real EPresid=m_ElectricPotential.solve(m_cellNData,true);
      pout() << " EPresid after regridding " << EPresid << endl;
    }
  printGridSize();

}
//*******************************************
void EBAMRPlasma::printGridSize()
{
  // ouput the (new) number of boxes
  long long totalPoints = 0;
  for (int ilev = 0; ilev < m_params.numLevels; ilev++)
    {
      long long pointsThisLevel = 0;
      for (LayoutIterator lit = m_grids[ilev].layoutIterator(); lit.ok(); ++lit)
        {
          pointsThisLevel += m_grids[ilev][lit()].numPts();
        }
      totalPoints += pointsThisLevel;
      pout() << "regridStep level# " << ilev
	     <<  " points on this level = " << pointsThisLevel <<  endl;
    }
  pout() << "regridStep:"
         <<  " total points = " << totalPoints <<  endl;
}
/*********/
/***************************/
void EBAMRPlasma::regridRhos(Vector<LevelData<EBCellFAB>* >& a_rhosData, Vector<DisjointBoxLayout>& a_oldDBL, const int& flag)
{
  CH_TIME("EBAMRPlasma::regrid");

  // this is a very inefficient way to the rhos regridding: Because I didnt want to
  // rewrite an interpolation algorithm I copy the IVFAB onto a EBCellFAb and use that class' patch interpolator

  const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
  int numComp=1;
  if(flag == 0)
    {
      int numLevels = m_grids.size();
      a_oldDBL.resize(numLevels);
      a_rhosData.resize(numLevels,NULL);
      IntVect ivGhost = m_ghost_Dk*IntVect::Unit;
      for (int ilev = 0; ilev < m_grids.size(); ilev++)
	{
	  EBCellFactory factory(m_ebisl[ilev]);
	  a_rhosData[ilev] = new LevelData<EBCellFAB>(m_grids[ilev], 1, ivGhost, factory);
	  EBLevelDataOps::setToZero(*a_rhosData[ilev]);
	  for (DataIterator dit = m_rhos[ilev]->dataIterator(); dit.ok(); ++dit)
	    {
	      BaseIVFAB<Real>& data = (*m_rhos[ilev])[dit()];
	      const EBGraph& ebgraph = data.getEBGraph();
	      const IntVectSet& ivsIrreg = data.getIVS(); 
	      for (VoFIterator vofit(ivsIrreg, ebgraph); vofit.ok(); ++vofit)
		{
		  VolIndex vof = vofit();
		  for (int i=0; i<numComp; i++)
		    (*a_rhosData[ilev])[dit()](vof, i) = data(vof, i);
		  IntVect gi = vof.gridIndex();
		  gi.shift(1,1);
		  VolIndex vofUP(gi,0);
		  for (int i=0; i<numComp; i++)
		    (*a_rhosData[ilev])[dit()](vofUP, i) = data(vof, i);
		}
	    }
	  a_oldDBL[ilev] = m_grids[ilev];
	}
      updateGhost(a_rhosData,100);
    }

  if(flag == 1)
    {
      
      a_rhosData.resize(m_params.numLevels,NULL);
  
      // save data for later copy
      //not using m_ebisl because it gets wiped later
      //the ebisl has to know about the fact that we really have
      //four ghost cells.
      int nGhostEBISL = 6;
      Interval interv(0,1-1);
      IntVect ivGhost = m_ghost_Dk*IntVect::Unit;
      
      for (int ilev = 0; ilev < m_grids.size(); ilev++)
	{
	  LevelData<EBCellFAB> stateSaved;
	  if(ilev < a_oldDBL.size())
	    {
	      EBISLayout ebislOld;
	      ebisPtr->fillEBISLayout(ebislOld, a_oldDBL[ilev], m_domain[ilev].domainBox(), nGhostEBISL);
	      EBCellFactory factoryOld(ebislOld);
	      stateSaved.define(a_oldDBL[ilev], 1, ivGhost, factoryOld);
	      a_rhosData[ilev]->copyTo(interv, stateSaved, interv);
	    }

	  EBCellFactory factoryNew(m_ebisl[ilev]);
	  if(a_rhosData[ilev] == NULL) a_rhosData[ilev]   = new LevelData<EBCellFAB>();
	  a_rhosData[ilev]->define(m_grids[ilev],1,ivGhost, factoryNew);

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
	      ebInterpVec.interpolate(*a_rhosData[ilev], *a_rhosData[ilev-1], interv);
	    }
      
	  // copy from old state
	  if(ilev < a_oldDBL.size())
	    stateSaved.copyTo(interv,*a_rhosData[ilev], interv);
	}
      
      updateGhost(a_rhosData,100);
    }

  if(flag == 2)
    {
      for (int ilev = 0; ilev <= m_finestLevel; ilev++)
	{
	  LayoutData<IntVectSet> irregSets(m_grids[ilev]);
	  for (DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit) 
	    {
	      Box grownBox = grow(m_grids[ilev].get(dit()), m_ghost_Dk); 
	      grownBox &= m_domain[ilev];
	      irregSets[dit()] = m_ebisl[ilev][dit()].getIrregIVS(grownBox);
	    }
	  if(m_rhos[ilev]==NULL)m_rhos[ilev]   = new LevelData<BaseIVFAB<Real> >();
	  BaseIVFactory<Real>  baseivfact(m_ebisl[ilev], irregSets);
	  m_rhos[ilev]->define(m_grids[ilev], 1    , m_ghost_Dk*IntVect::Unit, baseivfact);
	  for (DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit) 
	    {
	      BaseIVFAB<Real>& data = (*m_rhos[ilev])[dit()];
	      data.setVal(0e0);
	      const EBGraph& ebgraph = data.getEBGraph();
	      const IntVectSet& ivsIrreg = data.getIVS(); 
	      for (VoFIterator vofit(ivsIrreg, ebgraph); vofit.ok(); ++vofit)
		{
		  VolIndex vof = vofit();
		  for (int i=0; i<numComp; i++)
		    data(vof, i)=(*a_rhosData[ilev])[dit()](vof, i);
		}
	    }
	  delete a_rhosData[ilev];
	  m_rhos[ilev]->exchange(Interval(0,0));
	}
    }


}
/********************************/
void EBAMRPlasma::
getVelTemp(const bool& a_dimensionalScaling)
{
  //Multilevel call to NavierLevels in order to interpolate
  m_NavierLevels[0]->getVelTemp(m_Velo,m_Temp,m_grids,m_ebFineInterps);
  for (int ilev = 0; ilev <= m_amrS.getFinestLevel(); ilev++)
    {
      m_SpeciesLevels[ilev]->getElectronEnergy(*m_Temp[ilev],1);
    }
  // dimensional scaling
  if( a_dimensionalScaling )
    {
      EBAMRDataOps::scale(m_Velo,m_Vel0);
      EBAMRDataOps::scale(m_Temp,m_Tinf,0);// scale only gas temperature
    }
  updateGhost(m_Velo,10);
  updateGhost(m_Temp,10);
  /*
  Real maxVelo, minVelo; EBAMRDataOps::getMaxMin(maxVelo, minVelo, m_Velo,1);
  pout() << " Max Min Vert velocity " << maxVelo << ", " << minVelo << endl;
  */
}
/*************************/
void EBAMRPlasma::
getSpecies()
{

  // askMD:: the species are the state vector of m_amrS

  for (int ilev = 0; ilev <= m_amrS.getFinestLevel(); ilev++)
    {
      m_cellNData[ilev] = m_SpeciesLevels[ilev]->getStatePtr();
    }  
  updateGhost(m_cellNData,1e5);
  
  int ie = m_PlasmaPhysics.findSpec(string("e-"));
  Real maxElec, minElec; EBAMRDataOps::getMaxMin(maxElec, minElec, m_cellNData,ie);
  if(maxElec > m_PlasmaPhysics.regridThresh(ie)*1000)
    {
      m_SpeciesLevels[0]->PlotReactiveSource(2000);
      m_SpeciesLevels[0]->outputIVFAB(m_rhos, string("rhos"), 2000);
      pout() << "Stopping with max min elec concentration " << maxElec << ", " << minElec << endl;
      exit(0);
    }
  int ib = 0;
  Real maxBack, minBack; EBAMRDataOps::getMaxMin(maxBack, minBack, m_cellNData,ib);
  if(maxBack >20)
    {
      m_SpeciesLevels[0]->PlotReactiveSource(2000);
      m_SpeciesLevels[0]->outputIVFAB(m_rhos, string("rhos"), 2000);
      pout() << "Stopping with max min back concentration " << maxBack << ", " << minBack << endl;
      exit(0);
    }

  m_PlasmaPhysics.floor(m_cellNData);

  //hook to get the wall fluxes
  m_rhosRHS=m_SpeciesLevels[0]->getRhosRHS();
}
//*******************************************
void EBAMRPlasma::
getPlasmaParameters(PlasmaParameters&  a_params, bool a_forceSingleLevel)
{
  CH_TIME("EBAMRPlasma::getPlasmaParameters");


  // NS input only for LM code
  ParmParse ppViscous("ns");;
  // hardwire ebBcType
  a_params.ebBcType = 1;  // value is read thorugh parmparse
  ppViscous.get("domain_bc_type",a_params.domBcType);
  //temperature Boundary conditions
  ppViscous.get("eb_bc_typeT",a_params.ebBcTypeT);
  ppViscous.get("eb_bc_valueT",a_params.ebBCValueT);
  ppViscous.get("order_ebbc", a_params.orderEB);

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_params.ghostPhi[idir] = 4;
      a_params.ghostRHS[idir] = 0;
      if (a_params.ebBcType == 0 || a_params.ebBcType == 2)
        {
          a_params.ghostPhi[idir] = 1;
        }
    }
  if (ppViscous.contains("ghostPhi") && ppViscous.contains("ghostRhs"))
    {
      vector<int> ghost_phi_v(SpaceDim, 4);
      vector<int> ghost_rhs_v(SpaceDim, 4);
      ppViscous.queryarr("ghostPhi",ghost_phi_v,0,SpaceDim);
      ppViscous.queryarr("ghostRhs",ghost_rhs_v,0,SpaceDim);

      for (int idir = 0; idir < SpaceDim; idir++)
        {
          a_params.ghostPhi[idir] = ghost_phi_v[idir];
          a_params.ghostRHS[idir] = ghost_rhs_v[idir];
        }
    }

  // general input same in MD and LM codes
  ParmParse ppPlasma;


  std::vector<int> nCellsArray(SpaceDim);
  ppPlasma.getarr("n_cells",nCellsArray,0,SpaceDim);
  for (int idir = 0; idir < SpaceDim; idir++)
    a_params.nCells[idir] = nCellsArray[idir];

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
      
      int maxLevel = 0;
      ppPlasma.get("max_level",maxLevel);
      if(ppPlasma.contains("max_level_species"))
	{
	  int maxLevelSpecies = 0;
	  ppPlasma.get("max_level_species",maxLevelSpecies);
	  maxLevel = max(maxLevel,maxLevelSpecies);
	}
      a_params.maxLevel = maxLevel;
      a_params.numLevels = a_params.maxLevel + 1;
      ppPlasma.getarr("ref_ratio",a_params.refRatio,0,a_params.numLevels);
      ppPlasma.get("block_factor",a_params.blockFactor);
      ppPlasma.get("fill_ratio",a_params.fillRatio);
      ppPlasma.get("grid_buffer_size",a_params.bufferSize);
    }
  a_params.finestLevel=a_params.maxLevel;
  IntVect lo = IntVect::Zero;
  IntVect hi = a_params.nCells;
  hi -= IntVect::Unit; //cells index start from 0

  std::vector<int> is_periodica(SpaceDim,0);
  bool is_periodic[SpaceDim];
  ppPlasma.queryarr("is_periodic", is_periodica,0, SpaceDim);
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
  ppPlasma.getarr("domain_length",dLArray,0,SpaceDim);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_params.domainLength[idir] = dLArray[idir];
    }
  ppPlasma.get("max_grid_size",a_params.maxGridSize);

  //derived stuff
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_params.coarsestDx[idir] = a_params.domainLength[idir]/a_params.nCells[idir];
    }
  a_params.probLo = RealVect::Zero;
  a_params.probHi = RealVect::Zero;
  a_params.probHi += a_params.domainLength;

  a_params.noRefCorners = false;
  ppPlasma.query("no_ref_corners", a_params.noRefCorners);
  
}
void EBAMRPlasma::tagCells(IntVectSet& a_tags, const int& a_lev)
{
  CH_TIME("EBAMRPlasma::tagCells");
 
  Real SrefineThresh = 2.5e2;
  int StagBufferSize = 2;
  EBPatchPolytropic Ppatch;
  Interval momIntrv = Ppatch.momentumInterval();
  a_tags.makeEmpty();

  m_SpeciesLevels[a_lev]->makeFluidSource(*m_FSrc[a_lev]);
  for (DataIterator dit = m_grids[a_lev].dataIterator(); dit.ok(); ++dit)
    {
      const EBCellFAB& srcFAB = (*m_FSrc[a_lev])[dit()];
      Box gridBox = m_grids[a_lev].get(dit());
      const EBGraph& ebgraph = m_ebisl[a_lev][dit()].getEBGraph();
      IntVectSet ivsTot(gridBox);

      for (VoFIterator vofit(ivsTot, ebgraph); vofit.ok(); ++vofit)
	{
	  const VolIndex& vof = vofit();
	  const IntVect& iv = vof.gridIndex();
	  Real srcmag = 0.0;
	  
	  int tagSrcComps = SpaceDim;
	  for (int icomp = 0; icomp < tagSrcComps; icomp++)
	    {
	      int jcomp = momIntrv.begin() + icomp;
	      Real srcDirVal = srcFAB(vof, jcomp);
	      srcmag += srcDirVal*srcDirVal;
	    }
	  srcmag = sqrt(srcmag);
	  if (srcmag >= SrefineThresh)
	    {
	      a_tags |= iv;
	    }
	} //end loop over vofs
    }

  a_tags.grow(StagBufferSize);


  // Need to do this in two steps unless a IntVectSet::operator &=
  // (ProblemDomain) operator is defined
  Box localTagsBox = a_tags.minBox();
  localTagsBox &= m_domain[a_lev];
  a_tags &= localTagsBox;
}

/***********************/
void
EBAMRPlasma::
updateGhost(Vector<LevelData<EBCellFAB>* >&  a_phi, const Real& a_limit)
{
  CH_TIME("EBAMRPlasma::updateGhost");

  for (int ilev = 0; ilev < m_params.numLevels; ilev++)
    {
      LevelData<EBCellFAB>& ldata = *a_phi[ilev];
      int pvars = ldata.nComp();
      for (DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
	{
	  for (int icomp = 0; icomp < pvars; icomp++)
	    ldata[dit()].setCoveredCellVal(0.0,icomp);
          Box cbox = ldata[dit()].getRegion();
	  Box bbox = m_grids[ilev].get(dit());
	  IntVectSet ivs(cbox);
	  ivs -= bbox; // only ghost
	  for (VoFIterator vofit(ivs, m_ebisl[ilev][dit()].getEBGraph()); vofit.ok(); ++vofit)
	    {
	      const VolIndex& vof = vofit();
	      for (int icomp = 0; icomp < pvars; icomp++)
		{
		  if ( ! (abs(ldata[dit()](vof,icomp)) < a_limit))
		    {
		      ldata[dit()](vof,icomp) = 1.0;		      
		    }
		}
	      
	    }
	}
    }
  EBAMRDataOps::quadCFInterpAll(a_phi,m_eblg,m_params.refRatio);
  EBAMRDataOps::pwlFillPatchAll(a_phi,m_eblg,m_params.refRatio);
  EBAMRDataOps::averageDown(a_phi,m_eblg,m_params.refRatio);
  //averageDown(a_phi);
}

//-----------------------------------------------------------------------
void EBAMRPlasma::writeCheckpointFile(const char* a_inFile) const
{

  // copying from AMR::run, this is done before check pointing
  for (int level = 0; level < m_NavierLevels.size(); ++level)
    {
      m_NavierLevels[level]->time(m_amr.getCurrentTime());
    } 

  m_amr.writeCheckpointFile();

  string iter_str = m_amr.checkPointName();
  HDF5Handle handle(iter_str.c_str(), HDF5Handle::OPEN_RDWR);

  m_amrS.writeCheckpointData(handle);
  m_ElectricPotential.writeCheckpointFile(handle);

  
  for (int level = 0; level <= m_finestLevel; ++level)
    {
      // Setup the level string
      char levelStr[20];
      sprintf(levelStr,"%d",level);
      const std::string label = std::string("SPlevel_") + levelStr;
      handle.setGroup(label);
      write(handle,*m_rhos[level],"rhos");
    }

  handle.close();

  if(procID() == 0)
    {
      ofstream inputdeck;
      inputdeck.open(a_inFile,ios::app);
      inputdeck << "restart_file = " << m_amr.checkPointName() << endl;
      inputdeck.close();
    }

}
//-----------------------------------------------------------------------
Real EBAMRPlasma::currentTime() const
{
  return (m_amr.getCurrentTime()/m_Vel0);
}
//-----------------------------------------------------------------------
void EBAMRPlasma::stdOut(const Real& a_EPresid) const
{
      Real maxphi, minphi; EBAMRDataOps::getMaxMin(maxphi, minphi, m_phi,0);//lucacancel
      maxphi = max(abs(maxphi),abs(minphi));
      /*Real maxV = abs(m_ElectricPotential.getElectrodePotential());
      if(maxphi > 20*maxV)
	{
	  m_SpeciesLevels[0]->PlotReactiveSource(2000);
	  pout() << "Stopping with max min phi value " << maxphi << ", " << minphi <<   "; EPr= " << a_EPresid << endl;
	  exit(0);
	  }*/
      pout () << m_amr.getStep() << "; V_eb= " << m_ElectricPotential.getElectrodePotential() << "; EPr= " << a_EPresid
	      <<  "; phi_m= " << maxphi << "; ";
      EBAMRDataOps::checkThisData(m_cellNData, string(" m_cellNData "));
}

void EBAMRPlasma::conclude() const
{
  m_amr.conclude();
}
/****************/

