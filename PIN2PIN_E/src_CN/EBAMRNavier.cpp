#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cmath>

#include "parstream.H"
#include "ParmParse.H"
#include "PolyGeom.H"

#include "DebugOut.H"
#include "Box.H"
#include "Vector.H"
#include "IntVectSet.H"
#include "EBCellFAB.H"
#include "DisjointBoxLayout.H"
#include "LayoutIterator.H"
#include "CH_HDF5.H"
#include "SPMD.H"
#include "LoadBalance.H"
#include "ProblemDomain.H"
#include "BoxIterator.H"
#include "EBAMRIO.H"
#include "EBAMRDataOps.H"

#include "AMRLevel.H"
#include "EBAMRNavier.H"
#include "EBCellFactory.H"
#include "BaseIVFactory.H"
#include "VoFIterator.H"
#include "EBPatchGodunovF_F.H"
#include "EBPatchPolytropicF_F.H" 
#include "EBAMRNavierF_F.H" 
#include "EBIndexSpace.H"
#include "EBArith.H"
#include "EBFastFR.H"
#include "EBLoHiCenter.H"
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
#include "NamespaceHeader.H"

int  EBAMRNavier::s_NewPlotFile = 0;
bool EBAMRNavier::s_isLoadBalanceSet = false;
LoadBalanceFunc EBAMRNavier::s_loadBalance  = NULL;
IntVect ivdebamrN(D_DECL(16, 11, 9));

//define a ste of class variables
RefCountedPtr<EBLevelBackwardEuler>                   EBAMRNavier::s_viscLevTGA = RefCountedPtr<EBLevelBackwardEuler>();
RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB> > >   EBAMRNavier::s_viscAMRMG  = RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB> > >();
RefCountedPtr<EBViscousTensorOpFactory >              EBAMRNavier::s_viscOpFact = RefCountedPtr<EBViscousTensorOpFactory >();

RefCountedPtr<EBLevelBackwardEuler>                   EBAMRNavier::s_condLevTGA = RefCountedPtr<EBLevelBackwardEuler>();
RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB> > >   EBAMRNavier::s_condAMRMG  = RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB> > >();
RefCountedPtr<EBConductivityOpFactory >               EBAMRNavier::s_condOpFact = RefCountedPtr<EBConductivityOpFactory >();

BiCGStabSolver<LevelData< EBCellFAB> >                EBAMRNavier::s_botSolver;
Vector<EBViscousTensorOp *>                           EBAMRNavier::s_EBAMROps;
Vector<EBConductivityOp  *>                           EBAMRNavier::s_EBAMROpsT;
Real  EBAMRNavier::s_beta = 1e0;
Real  EBAMRNavier::s_betaCond = 0e0;
/***************************/
EBAMRNavier::EBAMRNavier()
{
  m_cfl = 0.8;
  m_tagAll = false;
  m_useMassRedist = true;
  m_doRZCoords = false;
  m_hasSourceTerm = false;
  m_doSmushing = true;
  m_origin = RealVect::Zero;
  m_dx = RealVect::Unit;
  m_aspect = RealVect::Unit;
  m_domainLength = RealVect::Unit;
  m_refineThresh = 0.2;
  m_initial_dt_multiplier = 0.1;
  m_ebPatchGodunov = NULL;
  m_srcE = NULL;
  m_redistRad = 1;
  m_isDefined = false;
  m_useInject=false;
  m_SFD = false;
  m_plot=0;
}
/***************************/
EBAMRNavier::~EBAMRNavier()
{
  if (m_ebPatchGodunov != NULL)
    delete m_ebPatchGodunov;
}
/***************************/
void
EBAMRNavier::
tagAll(bool a_tagAll)
{
  m_tagAll = a_tagAll;
}
/***************************/
void
EBAMRNavier::
doSmushing(bool a_doSmushing)
{
  m_doSmushing = a_doSmushing;
}
/***************************/
void
EBAMRNavier::
doRZCoords(bool a_doRZCoords)
{
  m_doRZCoords = a_doRZCoords;
}
/***************************/
void
EBAMRNavier::
hasSourceTerm(bool a_hasSourceTerm)
{
  m_hasSourceTerm = a_hasSourceTerm;
}
/***************************/
void
EBAMRNavier::
useMassRedistribution(bool a_useMassRedist)
{
  m_useMassRedist = a_useMassRedist;
}
/***************************/
void EBAMRNavier::redistRadius(int a_redistRad)
{
  m_redistRad = a_redistRad;
}
/***************************/
int EBAMRNavier::getRedistRadius()
{
  return(m_redistRad);
}
/***************************/
void EBAMRNavier::define(AMRLevel*  a_coarser_level_ptr,
                          const Box& a_problem_domain,
                          int        a_level,
                          int        a_ref_ratio)
{
  MayDay::Error("EBAMRNavier::define -\n\tShould never be called with a Box for a problem domain");
}
/***************************/
void EBAMRNavier::define(AMRLevel*  a_coarser_level_ptr,
                          const ProblemDomain& a_problem_domain,
                          int        a_level,
                          int        a_ref_ratio)
{
  CH_TIME("EBAMRNavier::define");
  if (s_verbosity >= 3)
    {
      pout() << "EBAMRNavier::define, level=" << a_level << endl;
    }

  m_isDefined = true;
  AMRLevel::define(a_coarser_level_ptr, a_problem_domain, a_level, a_ref_ratio);
  m_domainBox = m_problem_domain.domainBox();

  if (a_coarser_level_ptr != NULL)
    {
      EBAMRNavier* amrg_ptr = dynamic_cast<EBAMRNavier*>(a_coarser_level_ptr);
      if (amrg_ptr == NULL)
        {
          pout() << "EBAMRG::define:cannot cast  to derived class"  << endl;
          MayDay::Error();
        }

      m_cfl = amrg_ptr->m_cfl;
      m_domainLength = amrg_ptr->m_domainLength;
      m_refineThresh = amrg_ptr->m_refineThresh;
      m_tagBufferSize = amrg_ptr->m_tagBufferSize;
    }
  for (int idir = 0; idir < SpaceDim; ++idir)
    {
      m_dx[idir] = m_domainLength[idir]/m_domainBox.size(idir);
    }
  m_nGhost = 4;

  if (m_ebPatchGodunov != NULL)
    delete m_ebPatchGodunov;

  m_ebPatchGodunov = m_ebPatchGodunovFactory->create();
  m_ebPatchGodunov->define(m_problem_domain, m_dx);
  m_ebPolytropic = (static_cast <EBPatchPolytropic*> (m_ebPatchGodunov));

  m_nComp = m_ebPatchGodunov->numConserved();
  m_stateNames = m_ebPatchGodunov->stateNames();
  m_primNames = m_ebPatchGodunov->primNames();

  //definition of the viscous part
  ParmParse ppPlasma;
  m_nGhost = 4;
  m_hasDiffusion=true;      ppPlasma.query("has_diffusion",   m_hasDiffusion);
  m_explicitDiffusion=false;ppPlasma.query("explicit_diffusion",m_explicitDiffusion);
  
  //get the Reynolds and Prandtl numbers
  Real ReCH, ReL, MachInflow, Prandtl, Tinf, pinf;
  ppPlasma.get("Prandtl",    Prandtl);
  ppPlasma.get("gamma",      m_gamma);
  ppPlasma.get("Mach_inflow",MachInflow);
  ppPlasma.get("pinf",    pinf);

  
  PlasmaPhysics PlasmaPhysicsT;
  PlasmaPhysicsT.define();
  Real Runiv = 8.314e3; //J/kmole/K
  Real Rgas = Runiv/PlasmaPhysicsT.m_molWeight; // Joule/kg/K
  Tinf = PlasmaPhysicsT.m_Tg;
  // ReL is Re for Lchar = 1m, mu_inf for Helium 
  ReL = pinf*MachInflow*sqrt(m_gamma/(Rgas*Tinf))*(Tinf+79.4)/(1.4844e-6*pow(Tinf,1.5));
  ReCH = ReL/m_gamma/MachInflow;
  s_beta = -1.0/ReCH;
  m_alpha = 1e0;
  Real coe1 = 1.0/(Prandtl*(m_gamma-1.0));
  m_cv = 1.0/m_gamma/(m_gamma-1.0);
  s_betaCond = s_beta*coe1;
  m_lambdafac = -2./3.;
  //MDCheck if ok to use in file 
  m_Tinf = Tinf;

  pout() << m_level << "navier level; ReL " << ReL << ", Tinf " << m_Tinf  << endl;

  m_SFD = ppPlasma.contains("sfd_chi") && ppPlasma.contains("sfd_Dl");
  //m_SFD = false; // cancel
  if(m_SFD)
    {
      ppPlasma.get("sfd_chi", m_chi);
      ppPlasma.get("sfd_Dl", m_Dl);
    }

  //definition of refinement boxes
  defineRefinementBoxes();
}

/***************************/
Real EBAMRNavier::advance()
{
  CH_TIME("EBAMRNavier::advance");
  EBPatchGodunov::s_whichLev = m_level;
  dumpDebug(string("going into advance"));
  if (s_verbosity >= 3)
    {
      pout() << " in EBAMRNavier advance for level =" << m_level << ", with dt = " << m_dt << endl;
    }
  {
    CH_TIME("copy new to old");
    copyLevel(m_stateOld,m_stateNew,m_stateOld.interval(), m_stateNew.interval());
    //m_stateNew.copyTo(m_stateNew.interval(), m_stateOld, m_stateOld.interval());
  }

  Real new_dt = 0.0;
  //set up arguments to step
  //undefined lfr in case we need it
  EBFluxRegister lfr;
  //undefined leveldata in case we need it
  const LevelData<EBCellFAB> ld;

  //set arguments to dummy arguments and
  //then fix if they are available
  EBFluxRegister* coarFR = &lfr;
  EBFluxRegister* fineFR = &lfr;
  const LevelData<EBCellFAB>* coarDataOld = &ld;
  const LevelData<EBCellFAB>* coarDataNew = &ld;

  Real tCoarserOld = 0.0;
  Real tCoarserNew = 0.0;
  if (m_hasCoarser)
    {
      EBAMRNavier* coarPtr = getCoarserLevel();
      //recall that my flux register goes between my
      //level and the next finer level
      coarFR = &coarPtr->m_ebFluxRegister;
      coarDataOld = &coarPtr->m_stateOld;
      coarDataNew = &coarPtr->m_stateNew;
      tCoarserNew = coarPtr->m_time;
      tCoarserOld = tCoarserNew - coarPtr->m_dt;
      //time should never be greater than the newest coarse
      //time.  time might be very slightly smaller than
      //tCoarserOld because of the above subtraction.
      Real eps = 1.0e-10;
      if ((m_time > tCoarserNew) || (m_time < (tCoarserOld - eps)))
        {
          MayDay::Error("out of bounds time input to AMRGodunov");
        }
      //correct for said floating-point nastiness
      m_time = Min(Max(m_time, tCoarserOld), tCoarserNew) ;
    }
  if (m_hasFiner)
    {
      //recall that my flux register goes between my
      //level and the next finer level
      fineFR = &m_ebFluxRegister;
    }

#ifndef NDEBUG
  if (!m_hasCoarser && (s_verbosity > 1))
    {
      Real summass;
      int densityIndex = m_ebPatchGodunov->densityIndex();
      sumConserved(summass,  densityIndex);

      pout() << "sum mass = " << summass << endl;
    }
#endif

  {

    if (m_SFD) 
      {
	Interval consInterv(0, m_nComp-1);
	if (m_hasCoarser)
	  {
	    const EBAMRNavier* amrGodCoarserPtr = getCoarserLevel();
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
	m_ebLevelGodunov.getQbar(&m_qbar, m_chi);
    }
    //source term
    LevelData<EBCellFAB> diffusiveSrc;
    if(m_hasDiffusion && m_explicitDiffusion)
      {
	//PlotReactiveSource();//lucacancel
	EBCellFactory factory(m_ebisl);
	diffusiveSrc.define(m_grids, m_nComp, m_stateNew.ghostVect(),factory);
	makeDiffusiveSource(diffusiveSrc);
	//add plasma source
	EBLevelDataOps::incr(diffusiveSrc,m_srcF,1.0);
	m_ebLevelGodunov.setESource(&diffusiveSrc);
      }
    CH_TIME("levelgodunov step");
    EBPatchGodunov::setCurLevel(m_level);
    new_dt = m_ebLevelGodunov.step(m_stateNew,
                                   m_massDiff,
                                   *fineFR,
                                   *coarFR,
                                   *coarDataOld,
                                   *coarDataNew,
                                   m_time,
                                   tCoarserOld,
                                   tCoarserNew,
                                   m_dt);
  }
  
  //SFD Update
  if (m_SFD)
    { 
      Real acoef = m_Dl/(m_dt+m_Dl);Real bcoef = m_dt/(m_dt+m_Dl); 
      EBLevelDataOps::axby(m_qbar,m_qbar, m_stateNew,  acoef, bcoef);
    }

  if (m_hasDiffusion && !m_explicitDiffusion)
    {
      implicitDiffusion(fineFR,coarFR,*coarDataOld,*coarDataNew,tCoarserOld,tCoarserNew);      
    }
  dumpDebug(string("after levelgodunov"));

  if ((s_verbosity > 2))
    {
      pout() << "for this proc max  wave speed = " << EBPatchGodunov::getMaxWaveSpeed()
             << " at cell =  " << EBPatchGodunov::getMaxWaveSpeedIV() 
	     << "new_dt is " << new_dt << endl;
    }
  // flux register manipulation happens in levelgodunov
  // level redistribution happens in levelgodunov

  Interval interv(0, m_nComp-1);
  // increment redistribution register between this level
  // and next coarser level
  {
    CH_TIME("coarse-fine guano");
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
  }


  //lucacancel temporary
  EBCellFactory factory(m_ebisl);
  LevelData<EBCellFAB> consTemp(m_grids, m_nComp, m_nGhost*IntVect::Unit, factory);
  m_ebPatchGodunov->getEBPhysIBC()->initialize(consTemp, m_ebisl);
  
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

  //deal with time and time step
  m_time += m_dt;
  Real return_dt = m_cfl * new_dt;

  //save stable timestep to save computational effort
  m_dtNew = return_dt;

  dumpDebug(string("leaving advance   "));
  return return_dt;
}
/***************************/
void EBAMRNavier::postTimeStep()
{
  CH_TIME("EBAMRNavier::postTimeStep");
  if (s_verbosity >= 3)
    {
      pout() << " in EBAMRNavier postTimeStep for level " << m_level << endl;
    }
  Interval interv(0, m_nComp-1);
  if (m_hasCoarser)
    {
      if (m_doSmushing)
        {
          //redistibute to coarser level
          EBAMRNavier* coarPtr = getCoarserLevel();
          //if use mass weighting, need to
          //fix weights of redistribution object
          //EBFineToCoarRedist::s_verbose = false;
          //if (!m_hasFiner && m_hasCoarser)
          //  {
          //    EBFineToCoarRedist::s_verbose = true;
          //    EBFineToCoarRedist::s_ivdebug = ivdebamrN;
          //  }

          if (m_useMassRedist)
            {
              int densevar = m_ebPatchGodunov->densityIndex();
              coarPtr->m_stateNew.exchange(interv);

              m_ebFineToCoarRedist.resetWeights(coarPtr->m_stateNew, densevar);
            }

          m_ebFineToCoarRedist.redistribute(coarPtr->m_stateNew, interv);

          coarPtr->dumpDebug(string("after finetocoar  "));
        }// m_doSmushing

      m_ebFineToCoarRedist.setToZero();
    }
  if (m_hasFiner)
    {
      EBAMRNavier* finePtr = getFinerLevel();
      // reflux from finer level solution
      CH_assert(Abs(m_dx[0]-m_dx[1])<1.e-9);
      Real scale = -1.0/m_dx[0];

      m_ebFluxRegister.reflux(m_stateNew, interv, scale);

      dumpDebug(string("after reflux      "));

      //do averaging from finer level
      finePtr->m_ebCoarseAverage.average(m_stateNew, finePtr->m_stateNew, interv);


      //the flux register must modify the redistribution registers


      m_ebFluxRegister.incrementRedistRegister(m_ebCoarToFineRedist,
                                               interv, scale);

      m_ebFluxRegister.incrementRedistRegister(m_ebCoarToCoarRedist,
                                               interv, scale);

      if (m_doSmushing)
        {
          //if use mass weighting, need to
          //fix weights of redistribution object
          if (m_useMassRedist)
            {
              int densevar = m_ebPatchGodunov->densityIndex();
              m_stateNew.exchange(interv);
              dumpDebug(string("before ctofresetw "));
              m_ebCoarToFineRedist.resetWeights(m_stateNew, densevar);
              dumpDebug(string("before ctocresetw "));
              m_ebCoarToCoarRedist.resetWeights(m_stateNew, densevar);
              dumpDebug(string("after  ctocresetw "));
            }
          //redistibute to finer level
          m_ebCoarToFineRedist.redistribute(finePtr->m_stateNew, interv);

          finePtr->dumpDebug(string("after coartofine  "));

          //do the re redistirubtion
          m_ebCoarToCoarRedist.redistribute(m_stateNew, interv);

          dumpDebug(string("after coartocoar  "));
        }//m_doSmushing

      // average from finer level data
      finePtr->m_ebCoarseAverage.average(m_stateNew, finePtr->m_stateNew, interv);
      if(m_SFD) finePtr->m_ebCoarseAverage.average(m_qbar, finePtr->m_qbar, interv);


      dumpDebug(string("after average     "));
      m_ebCoarToFineRedist.setToZero();
      m_ebCoarToCoarRedist.setToZero();
    } //m_hasfiner
}
/****************************/
void EBAMRNavier::tagCells(IntVectSet& a_tags)
{
  CH_TIME("EBAMRNavier::tagCells");
  if (s_verbosity >= 3)
  {
  pout() << "EBAMRNavier::tagCells for level " << m_level << "Threshold " << m_refineThresh  << endl;
      }

          // If there is a coarser level interpolate undefined ghost cells
  int densityIndex = m_ebPatchGodunov->densityIndex();
  Interval momIntrv = m_ebPolytropic->momentumInterval();
  Interval intervDensMom(Min(densityIndex,momIntrv.begin()), Max(densityIndex,momIntrv.end()));
  EBCellFactory factory(m_ebisl);
  int nghost = m_nGhost;
  int nCons = m_ebPatchGodunov->numConserved();
  LevelData<EBCellFAB> consTemp(m_grids, nCons, (IntVect::Unit)*nghost, factory);
  Interval consInterv(0, nCons-1);
  m_stateNew.copyTo(consInterv, consTemp, consInterv);
  if (m_hasCoarser)
    {
      const EBAMRNavier* amrGodCoarserPtr = getCoarserLevel();
      int refRatCrse = amrGodCoarserPtr->refRatio();
      EBPWLFillPatch patcher(m_grids,
			     amrGodCoarserPtr->m_grids,
			     m_ebisl,
			     amrGodCoarserPtr->m_ebisl,
			     amrGodCoarserPtr->m_domainBox,
			     refRatCrse, m_nComp, nghost);
      
      Real coarTimeOld = 0.0;
      Real coarTimeNew = 1.0;
      Real fineTime    = 0.0;
      patcher.interpolate(consTemp,
			  amrGodCoarserPtr->m_stateOld,
			  amrGodCoarserPtr->m_stateNew,
			  coarTimeOld,
			  coarTimeNew,
			  fineTime,
			  intervDensMom);
    }
  consTemp.exchange(intervDensMom);
  consToVelTemp(m_initialV,m_initialT,consTemp);
  computeVorticity(m_updateV, m_initialV);
  a_tags.makeEmpty();
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      const EBCellFAB& vortFAB = m_updateV[dit()];
      Box gridBox = m_grids.get(dit());
      const EBGraph& ebgraph = m_ebisl[dit()].getEBGraph();
      IntVectSet ivsTot(gridBox);
      Box shrunkDom = m_problem_domain.domainBox(); //noshrinking here(refer to the EBAMViscous implementation for shrinking)
      ivsTot &= shrunkDom;
      if(!m_vortBox.isEmpty()) ivsTot &=  m_vortBox;

      if (m_refineThresh > 0)
	{
	  for (VoFIterator vofit(ivsTot, ebgraph); vofit.ok(); ++vofit)
	    {
	      const VolIndex& vof = vofit();
	      const IntVect& iv = vof.gridIndex();
	      Real vortmag = 0.0;

	      int tagVortComps = SpaceDim;
	      for (int icomp = 0; icomp < tagVortComps; icomp++)
		{
		  Real vortDirVal = vortFAB(vof, icomp);
		  vortmag += vortDirVal*vortDirVal;
		}
	      vortmag = sqrt(vortmag);
	      if (vortmag >= m_refineThresh)
		{
		  a_tags |= iv;
		}
	    } //end loop over vofs
	}

        IntVectSet irregIVS = ebgraph.getIrregCells(gridBox);
        a_tags |= irregIVS;
    }

  a_tags.grow(m_tagBufferSize);

  // Need to do this in two steps unless a IntVectSet::operator &=
  // (ProblemDomain) operator is defined
  Box localTagsBox = a_tags.minBox();
  localTagsBox &= m_problem_domain;
  a_tags &= localTagsBox;
}
/***************************/
void EBAMRNavier::tagCellsInit(IntVectSet& a_tags)
{
  if (s_verbosity >= 3)
    {
      pout() << "EBAMRNavier::tagCellsInit for level " << m_level << endl;
    }

  // Since tags are calculated using only current time step data, use
  // the same tagging function for initialization and for regridding.
  tagCells(a_tags);
}
/***************************/
void EBAMRNavier::regrid(const Vector<Box>& a_new_grids)
{
  CH_TIME("EBAMRNavier::regrid");
  if (s_verbosity >= 3)
    {
      pout() << " in EBAMRNavier regrid for level " << m_level << endl;
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

  ebisPtr->fillEBISLayout(m_ebisl, m_grids, m_domainBox, nGhostEBISL); //@regrid

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
      EBAMRNavier* coarPtr = getCoarserLevel();
      m_ebFineInterp.interpolate(m_stateNew, coarPtr->m_stateNew, interv);
      if(m_SFD) m_ebFineInterp.interpolate(m_qbar, coarPtr->m_qbar, interv);
    }

  // copy from old state
  stateSaved.copyTo(interv,m_stateNew, interv);
  if(m_SFD) stateSavedBar.copyTo(interv,m_qbar, interv);
  defineSolvers();
}
/***************************/
void EBAMRNavier::initialGrid(const Vector<Box>& a_new_grids)
{
  CH_TIME("EBAMRG::initialGrid");

  if (s_verbosity >= 3)
    {
      pout() << " in EBAMRNavier initialGrid for level " << m_level << endl;
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
  if (s_verbosity >= 3)
    {
      pout() << " just loadbalanced " << m_level << endl;
    }

  m_grids = DisjointBoxLayout(a_new_grids,proc_map);

  if (s_verbosity >= 3)
    {
      pout() << "EBAMRNavier::initialgrid grids " << endl;
      dumpDBL(&m_grids);
    }

  //the ebisl has to know about the fact that we really have
  //four ghost cells.
  int nGhostEBISL = 6;
  if (s_verbosity >= 3)
    {
      pout() << " EBAMRNavier::about to fill ebislayout  in EBAMRSpecies initialGrid for level " << m_level << endl;
    }

  ebisPtr->fillEBISLayout(m_ebisl, m_grids,  m_domainBox, nGhostEBISL); //@initialGrid


  EBCellFactory factoryNew(m_ebisl);
  IntVect ivGhost = m_nGhost*IntVect::Unit;
  m_stateNew.define(m_grids,m_nComp, ivGhost, factoryNew);
  m_stateOld.define(m_grids,m_nComp, ivGhost, factoryNew);
  if(m_SFD) m_qbar.define(m_grids,m_nComp, ivGhost, factoryNew);

  // set up data structures
  levelSetup();  //@initialgrid
}
/***************************/
void EBAMRNavier::initialData()
{
  CH_TIME("EBAMRNavier::initialData");

  if (s_verbosity >= 3)
    {
      pout() << " in EBAMRNavier initialData for level " << m_level << endl;
    }

  const EBPhysIBC* const ebphysIBCPtr = m_ebPatchGodunov->getEBPhysIBC();

  //initialize both new and old states to
  //be the same thing
  ebphysIBCPtr->initialize(m_stateNew, m_ebisl);
  ebphysIBCPtr->initialize(m_stateOld, m_ebisl);  
  if(m_SFD) ebphysIBCPtr->initialize(m_qbar, m_ebisl);
}
/***************************/
void EBAMRNavier::postInitialize()
{
  if (s_verbosity >= 3)
    {
      pout() << "EBAMRSpecies::postInitialize " << m_level << endl;
    }
}
/***************************/
void EBAMRNavier::syncWithFineLevel()
{
  CH_TIME("EBAMRG::syncWithFineLevel");
  if (s_verbosity >= 3)
    {
      pout() << " in EBAMRNavier syncWithFineLevel for level " << m_level << endl;
    }
  //stuff that needs to be setup from the finer
  //level.  A bunch of objects depend upon the layouts
  //from both levels and the finer level changes more
  //often from regrid so this needs to be called from the finer
  //level
  CH_assert(m_hasFiner);
  if (m_hasFiner)
    {
      EBAMRNavier* finePtr = getFinerLevel();
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
/***************************/
void EBAMRNavier::patchGodunov(const EBPatchGodunovFactory* const a_ebPatchGodunovFactory)
{
  m_ebPatchGodunovFactory = a_ebPatchGodunovFactory;
}
/***************************/
Real EBAMRNavier::computeDt()
{
  Real newDt;
  newDt = m_dtNew;

  return newDt;
}
/***************************/
Real EBAMRNavier::computeInitialDt()
{
  Real maxwavespeed =  m_ebLevelGodunov.getMaxWaveSpeed(m_stateNew);
  CH_assert(Abs(m_dx[0]-m_dx[1])<1.e-9);
  Real newDT = m_initial_dt_multiplier * m_dx[0] /maxwavespeed;
  m_dt = newDT;

  return newDT;
}
/***************************/
void EBAMRNavier::CFL(Real a_cfl)
{
  m_cfl = a_cfl;
}
/***************************/
void EBAMRNavier::domainLength(RealVect a_domainLength)
{
  m_domainLength = a_domainLength;
}
/***************************/
void EBAMRNavier::refinementThreshold(Real a_refineThresh)
{
  m_refineThresh = a_refineThresh;
}
/***************************/
void EBAMRNavier::tagBufferSize(int a_tagBufferSize)
{
  m_tagBufferSize = a_tagBufferSize;
}

/***************************/
void
EBAMRNavier::sumConserved(Real& a_sumcons,
                           const int& a_ivar) const
{
  CH_TIME("EBAMRG::sumConserved");
  Real sumgrid = 0;
  for (DataIterator dit= m_grids.dataIterator(); dit.ok(); ++dit)
    {
      const EBCellFAB& state = m_stateNew[dit()];
      Box thisBox = m_grids.get(dit());
      IntVectSet uberIVS(thisBox);
      const EBISBox& ebisBox = m_ebisl[dit()];
      for (VoFIterator vofit(uberIVS, ebisBox.getEBGraph());
          vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          Real consVal = state(vof, a_ivar);
          Real volFrac = ebisBox.volFrac(vof);
          Real volume = volFrac;
          if (m_doRZCoords)
            {
              Real cellVol, kvol;
              EBArith::getKVolRZ(kvol, cellVol, ebisBox, m_dx[0], vof);
              volume = cellVol*kvol;
            }
          sumgrid += consVal*volume;
        }
    }

  Vector<Real> all_sum;
  gather(all_sum,sumgrid,uniqueProc(SerialTask::compute));
  Real sumallgrid = 0.;
  if (procID() == uniqueProc(SerialTask::compute))
    {
      for (int i = 0; i < all_sum.size(); ++i)
        {
          sumallgrid += all_sum[i];
        }
    }
  broadcast(sumallgrid,uniqueProc(SerialTask::compute));
  a_sumcons = sumallgrid;
}
/***************************/
EBAMRNavier*
EBAMRNavier::getCoarserLevel() const
{
  EBAMRNavier* retval = NULL;
  if (m_coarser_level_ptr != NULL)
    {
      retval = dynamic_cast <EBAMRNavier*> (m_coarser_level_ptr);

      if (retval == NULL)
        {
          pout() << "EBAMRG::getCoarserLevel: dynamic cast failed"
                 << endl;
          MayDay::Error();
        }
    }
  return retval;
}

/***************************/
EBAMRNavier*
EBAMRNavier::getFinerLevel() const
{
  EBAMRNavier* retval = NULL;
  if (m_finer_level_ptr != NULL)
    {
      retval = dynamic_cast <EBAMRNavier*> (m_finer_level_ptr);

      if (retval == NULL)
        {
          pout() << "EBAMRG::getFinerLevel: dynamic cast failed"
                 << endl;
          MayDay::Error();
        }
    }
  return retval;
}
/***************************/
void
EBAMRNavier::levelSetup()
{
  CH_TIME("EBAMRG::levelSetup");
  if (s_verbosity >= 3)
    {
      pout() << " in EBAMRNavier levelSetup for level " << m_level << endl;
    }

  EBAMRNavier* coarPtr = getCoarserLevel();
  EBAMRNavier* finePtr = getFinerLevel();

  m_hasCoarser = (coarPtr != NULL);
  m_hasFiner   = (finePtr != NULL);
  //same as above but check that the finer grid has points
  m_hasFinerLevel   = hasFinerLevel();

  if (m_hasCoarser)
    {
      int nRefCrse = m_coarser_level_ptr->refRatio();
      const DisjointBoxLayout& coarser_m_grids = coarPtr->m_grids;
      const EBISLayout& coarser_m_ebisl = coarPtr->m_ebisl;
      const Box& domainCoar = coarPtr->m_domainBox;

      {
        CH_TIME("ave_interp_defs");
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
      }

      // maintain levelgodunov
      m_ebLevelGodunov.define(m_grids,
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
                              m_ebPatchGodunovFactory,
                              m_hasCoarser,
                              m_hasFiner);

      //define fine to coarse redistribution object
      //for now set to volume weighting
      {
        CH_TIME("fineToCoar_defs");
        m_ebFineToCoarRedist.define(m_grids, coarser_m_grids, m_ebisl, coarser_m_ebisl,
	   domainCoar, nRefCrse, m_nComp); // (full) definition
	//m_ebFineToCoarRedist.define(m_eblg, coarser_m_eblg, nRefCrse, m_nComp, 1);	
        m_ebFineToCoarRedist.setToZero();
      }

      coarPtr->syncWithFineLevel();
    }
  else
    {
      m_ebLevelGodunov.define(m_grids,
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
                              m_ebPatchGodunovFactory,
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


  //define diffusive data-containers
  EBCellFactory factory(m_ebisl);
  m_initialV.define(m_grids, SpaceDim, m_stateNew.ghostVect(), factory);
  m_initialT.define(m_grids, 1       , m_stateNew.ghostVect(), factory);
  m_srcV.define(m_grids, SpaceDim, m_stateNew.ghostVect(), factory);
  m_srcT.define(m_grids, 1       , m_stateNew.ghostVect(), factory);
  m_updateV.define(m_grids, SpaceDim, m_stateNew.ghostVect(), factory);
  m_updateT.define(m_grids, 1       , m_stateNew.ghostVect(), factory);
  EBLevelDataOps::setCoveredVal(m_initialV,0e0);
  EBLevelDataOps::setCoveredVal(m_initialT,0e0);
  if (m_hasCoarser)
    {
      EBAMRNavier* coarserPtr = getCoarserLevel();
      EBCellFactory Cfactory(coarserPtr->m_ebisl);
      m_newVCoarse.define(coarserPtr->m_grids, SpaceDim, m_stateNew.ghostVect(), Cfactory);
      m_newTCoarse.define(coarserPtr->m_grids, 1       , m_stateNew.ghostVect(), Cfactory);
      m_oldVCoarse.define(coarserPtr->m_grids, SpaceDim, m_stateNew.ghostVect(), Cfactory);
      m_oldTCoarse.define(coarserPtr->m_grids, 1       , m_stateNew.ghostVect(), Cfactory);
      EBLevelDataOps::setCoveredVal(m_newVCoarse,0e0);
      EBLevelDataOps::setCoveredVal(m_oldVCoarse,0e0);
      EBLevelDataOps::setCoveredVal(m_newTCoarse,0e0);
      EBLevelDataOps::setCoveredVal(m_oldTCoarse,0e0);
    }

  //define Fluid Source container
  m_srcF.define(m_grids, m_nComp, m_stateNew.ghostVect(), factory);
  EBLevelDataOps::setToZero(m_srcF);

}
/***************************/
LevelData<EBCellFAB>&
EBAMRNavier::getStateNew()
{
  return m_stateNew;
}
/***************************/
LevelData<EBCellFAB>&
EBAMRNavier::getStateOld()
{
  return m_stateOld;
}
/***************************/
Real
EBAMRNavier::getDt() const
{
  return m_dt;
}
/***************************/
const EBISLayout&
EBAMRNavier::getEBISLayout() const
{
  return m_ebisl;
}

/***************************/
void EBAMRNavier::getVelTemp(LevelData<EBCellFAB>&  a_vel,LevelData<EBCellFAB>&  a_temp, const bool& a_removeSelfSimilar)
{
  
  consToVelTemp(a_vel,a_temp,m_stateNew);
  if(a_removeSelfSimilar)
    {
      EBCellFactory factory(m_ebisl);
      LevelData<EBCellFAB> consTemp(m_grids, m_nComp, m_nGhost*IntVect::Unit, factory);
      m_ebPatchGodunov->getEBPhysIBC()->initialize(consTemp, m_ebisl);
      consToVelTemp(m_updateV,m_updateT,consTemp);
      EBLevelDataOps::incr(a_vel,m_updateV,-1.0);
      EBLevelDataOps::incr(a_temp,m_updateT,0,0,1,-1.0);
    }

}

/***************************/
void EBAMRNavier::getVelTemp(Vector<LevelData<EBCellFAB>* >  a_vel, Vector<LevelData<EBCellFAB>* >   a_temp, Vector<DisjointBoxLayout>&   a_grids, Vector<EBPWLFineInterp* > a_ebFineInterps)
{

  Vector<EBAMRNavier*>        hierarchy;
  Vector<DisjointBoxLayout>    grids;
  Vector<EBISLayout>           ebisl;
  Vector<EBLevelGrid>          eblg;
  Vector<int>                  refRat;
  ProblemDomain                lev0Dom;
  Real                         lev0Dx;
  Vector<ProblemDomain>        domains;
  Vector<Real>                 dxs;
  getHierarchyAndGrids(hierarchy, grids, ebisl, eblg, refRat, lev0Dom, lev0Dx, domains, dxs);
  Interval vInterv(0,SpaceDim-1);
  Interval tInterv(0,0);

  //
  int numLevels = hierarchy.size();
  int numLevelsSpec = a_grids.size();
  //CH_assert(numLevelsSpec >= numLevels);
  for (int ilev = 0; ilev < numLevelsSpec; ++ilev)
    {
  
      if(ilev < numLevels) consToVelTemp(hierarchy[ilev]->m_initialV,hierarchy[ilev]->m_initialT,hierarchy[ilev]->m_stateNew);
      if (ilev > 0 ) 
	{
	  a_ebFineInterps[ilev]->interpolate(*a_vel[ilev],  *a_vel[ilev-1],  vInterv);
	  a_ebFineInterps[ilev]->interpolate(*a_temp[ilev], *a_temp[ilev-1], tInterv);
	}
      if(ilev < numLevels) 
	{
	  hierarchy[ilev]->m_initialV.copyTo(vInterv,*a_vel[ilev], vInterv);
	  hierarchy[ilev]->m_initialT.copyTo(tInterv,*a_temp[ilev], tInterv);
	}
    }

}
/***************************/
void EBAMRNavier::setQbar(LevelData<EBCellFAB>*  a_qbar,Real a_chi)
{
  m_ebLevelGodunov.getQbar(a_qbar, a_chi);
}
void EBAMRNavier::setSource(LevelData<EBCellFAB>*  a_src)
{
  if (s_verbosity >= 3)
    {
      pout() << "EBAMRNavier::setSource1, level=" << m_level << endl;
    }

  if(m_hasDiffusion && m_explicitDiffusion)
    m_srcE = a_src;
  else
    m_ebLevelGodunov.setESource(a_src);
}

void EBAMRNavier::setSource(Vector<LevelData<EBCellFAB>* >  a_src)
{
  if (s_verbosity >= 3)
    {
      pout() << "EBAMRNavier::setSource2 level=" << m_level << endl;
    }
  Vector<EBAMRNavier*>        hierarchy;
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
  int numLevelsSpec = a_src.size();
  CH_assert(numLevelsSpec >= numLevels);
  Interval interv(0,m_nComp-1);
  for (int ilev = 0; ilev < numLevels; ++ilev)
    {
      if(ilev>0) hierarchy[ilev]->m_ebFineInterp.interpolate(hierarchy[ilev]->m_srcF, hierarchy[ilev-1]->m_srcF, interv);
      a_src[ilev]->copyTo(interv,hierarchy[ilev]->m_srcF, interv);  
    }


  if(! m_hasDiffusion || !m_explicitDiffusion)
    {
      for (int ilev = 0; ilev < numLevels; ++ilev)
	hierarchy[ilev]->m_ebLevelGodunov.setESource(&(hierarchy[ilev]->m_srcF));
    }


  if (s_verbosity >= 3)
    {
      pout() << "EBAMRNavier::setSource2 (end) level=" << m_level << endl;
    }

}

/***************************/
void EBAMRNavier::fillConsAndPrim(LevelData<EBCellFAB>& a_data) const
{
  CH_TIME("EBAMRG::fillConsAndPrim");
  if (s_verbosity >= 3)
    {
      pout() << "EBAMRNavier::fillConsAndPrim" << endl;
    }
  //fill output data including a bunch of geometric crud
  int nCons = m_ebPatchGodunov->numConserved();
  int nPrim = m_ebPatchGodunov->numPrimitives() ;
  int consAndPrim = nCons + nPrim;

  EBCellFactory ebcellfact(m_ebisl);
  a_data.define(m_grids, consAndPrim, IntVect::Zero, ebcellfact);

  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      const EBISBox& ebisbox = m_ebisl[dit()];
      const Box& grid = m_grids.get(dit());
      const EBCellFAB& consfab = m_stateNew[dit()];
      EBCellFAB primfab(ebisbox, grid, nPrim);
      //cfivs, time and timestep fake and not used here
      Real faket = 1.0;
      int logflag = 0;
      IntVectSet emptyivs;
      m_ebPatchGodunov->setValidBox(grid, ebisbox, emptyivs, faket, faket);
      m_ebPatchGodunov->consToPrim(primfab, consfab, grid,logflag);

      EBCellFAB& outputfab = a_data[dit()];

      Interval consIntervSrc(0, nCons-1);
      Interval consIntervDst(0, nCons-1);
      Interval primIntervSrc(0, nPrim-1);
      Interval primIntervDst(nCons, consAndPrim-1);

      // copy regular data
      outputfab.copy(grid, consIntervDst,  grid, consfab, consIntervSrc);
      outputfab.copy(grid, primIntervDst,  grid, primfab, primIntervSrc);
      Real coveredVal = -10.0;
      for (int ivar = 0; ivar < consAndPrim; ivar++)
        {
          outputfab.setInvalidData(coveredVal, ivar);
        }

    }//end loop over grids
}
/***************************/
#ifdef CH_USE_HDF5
/***************************/
void EBAMRNavier::writeCheckpointHeader(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
    {
      pout() << "EBAMRNavier::writeCheckpointHeader" << endl;
    }

  //stuff in non-eb checkpoint header is already
  //set in the define function, such as
  // Setup the number of components
  // Setup the component names
  //so i will eliminate the middleman
}
/***************************/
void EBAMRNavier::readCheckpointHeader(HDF5Handle& a_handle)
{
  if (s_verbosity >= 3)
    {
      pout() << "EBAMRNavier::readCheckpointHeader" << endl;
    }

  //stuff in non-eb checkpoint header is already
  //set in the define function, such as
  // Setup the number of components
  // Setup the component names
  // So i will eliminate the middleman.
}
/***************************/
void EBAMRNavier::writeCheckpointLevel(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
    {
      pout() << "EBAMRNavier::writeCheckpointLevel" << endl;
    }

  //reasons for deviations from non-eb stuff
  // the tag buffer size is set by the factory.
  // dx is set in the define function
  // the problem domain is set in the define function

  // Setup the level string
  char levelStr[20];
  sprintf(levelStr,"%d",m_level);
  const std::string label = std::string("level_") + levelStr;

  a_handle.setGroup(label);

  // Setup the level header information
  HDF5HeaderData header;

  header.m_int ["ref_ratio"]       = m_ref_ratio;
  header.m_real["dt"]              = m_dt;

  if (s_verbosity >= 3)
    {
      pout() << "EBAMRNavier::writeCheckpointLevel " << label << endl;
    }

  // Write the header for this level
  header.writeToFile(a_handle);

  write(a_handle,m_grids,"NSgrid");
  write(a_handle,m_stateOld,"dataOld");
  write(a_handle,m_stateNew,"dataNew");
}
/***************************/
/***************************/
void EBAMRNavier::readCheckpointLevel(HDF5Handle& a_handle)
{
  CH_assert(m_isDefined);
  if (s_verbosity >= 3)
    {
      pout() << "EBAMRNavier::readCheckpointLevel" << endl;
    }

  // Setup the level string
  char levelStr[20];
  sprintf(levelStr,"%d",m_level);
  const std::string label = std::string("level_") + levelStr;

  // Read the header for this level
  a_handle.setGroup(label);

  HDF5HeaderData header;
  header.readFromFile(a_handle);

  // Get the refinement ratio
  if (header.m_int.find("ref_ratio") == header.m_int.end())
    {
      MayDay::Error("::readCheckpointLevel: file does not contain ref_ratio");
    }
  m_ref_ratio = header.m_int["ref_ratio"];

  //reasons for deviations from non-eb stuff
  // the tag buffer size is set by the factory.
  // dx is set in the define function
  // the problem domain is set in the define function

  // Get dt
  if (header.m_real.find("dt") == header.m_real.end())
    {
      MayDay::Error("readCheckpointLevel: file does not contain dt");
    }
  m_dt = header.m_real["dt"];

  // Get the grids
  Vector<Box> vboxGrids;
  const int gridStatus = read(a_handle, vboxGrids,"NSgrid");
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
  ebisPtr->fillEBISLayout(m_ebisl, m_grids, m_domainBox, nGhostEBISL); //@readCheckpointLevel
  EBCellFactory factoryNew(m_ebisl);
  //m_nghost is set in define function
  IntVect ivGhost = m_nGhost*IntVect::Unit;
  m_stateNew.define(m_grids,m_nComp, ivGhost, factoryNew);
  m_stateOld.define(m_grids,m_nComp, ivGhost, factoryNew);
  EBLevelDataOps::setToZero(m_stateNew);
  EBLevelDataOps::setToZero(m_stateOld);

  //  Interval vars(0, m_nComp-1);
  //the false says to not redefine the data
  int dataStatusNew = read<EBCellFAB>(a_handle,
                                      m_stateNew,
                                      "dataNew",
                                      m_grids,
                                      Interval(),
                                      false);

  int dataStatusOld = read<EBCellFAB>(a_handle,
                                      m_stateOld,
                                      "dataOld",
                                      m_grids,
                                      Interval(),
                                      false);

  if ((dataStatusNew != 0) || (dataStatusOld != 0))
    {
      MayDay::Error("file does not contain state data");
    }
  // Set up data structures
  levelSetup(); //@readCheckpointLevel


  if(m_SFD)
    {
      m_qbar.define(m_grids,m_nComp, ivGhost, factoryNew);
      m_stateNew.copyTo(m_stateNew.interval(), m_qbar, m_qbar.interval());
      //for (DataIterator dit = m_stateNew.dataIterator(); dit.ok(); ++dit)
      //{
      //	  m_qbar[dit()].clone(m_stateNew[dit()]);
      //	  Box region = m_stateNew[dit()].getRegion(); //includes ghost cells
      //	  m_qbar[dit()].copy(region, m_qbar.interval(), region, m_stateNew[dit()], m_stateNew.interval());
      //	}
      
      Interval consInterv(0, m_nComp-1);
      if (m_hasCoarser)
	{
	  const EBAMRNavier* amrGodCoarserPtr = getCoarserLevel();
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
			      amrGodCoarserPtr->m_stateOld,
			      amrGodCoarserPtr->m_stateNew,
			      coarTimeOld,
			      coarTimeNew,
			      fineTime,
			      consInterv);
	}
      m_qbar.exchange(consInterv);
    }
}
/***************************/
void EBAMRNavier::writePlotHeader(HDF5Handle& a_handle) const
{
  if (s_NewPlotFile == 0)
    {
      writePlotHeaderOld(a_handle);
      return;
    }
  Vector<int> refRatios;
  const EBAMRNavier* current = this;
  int nlevs = 0;
  while (current != NULL)
  {
    refRatios.push_back(current->refRatio());
    nlevs++;
    current = (const EBAMRNavier*)(current-> m_finer_level_ptr);
  }

  headerEBFile(a_handle, nlevs, refRatios,
               m_problem_domain.domainBox(),
               m_origin, m_dx, m_aspect,
               m_stateNew.ghostVect());

  writeCellCenteredNames(a_handle, m_stateNames);
  a_handle.setGroup("/Expressions");
  HDF5HeaderData expressions;
  m_ebPatchGodunov->expressions(expressions);
  expressions.writeToFile(a_handle);

}
/***************************/
void EBAMRNavier::writePlotHeaderOld(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
    {
      pout() << "EBAMRNavier::writePlotHeader" << endl;
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

void EBAMRNavier::writePlotLevel(HDF5Handle& a_handle) const
{
  if (s_NewPlotFile == 0)
    {
      writePlotLevelOld(a_handle);
      return;
    }
  if (s_verbosity >= 3)
    {
      pout() << "EBAMRNavier::writePlotLevel " << m_level<< endl;
    }
  writeCellCentered(a_handle, m_level, &m_stateNew);
}

/***************************/
void EBAMRNavier::writePlotLevelOld(HDF5Handle& a_handle) const
{

  if (s_verbosity >= 3)
    {
      pout() << "EBAMRNavier::writePlotLevel" << endl;
    }
  //fill output data including a bunch of geometric crud
  int nCons = m_ebPatchGodunov->numConserved();
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

      // m_ebPatchGodunov->setValidBox(grid, ebisbox, emptyivs, faket, faket);

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
#endif
/*******/
/*                             Diffusive Part               */
/*******/
/*****************************/
void EBAMRNavier:: setSymmFlag(const int & a_symmetric)
{
  m_symmetric = a_symmetric;
}
/*****************************/
void EBAMRNavier:: setViscPntrs(Vector<RefCountedPtr<LevelData<EBCellFAB> > >*        a_aco,
				Vector<RefCountedPtr<LevelData<EBFluxFAB> > >*        a_eta,
				Vector<RefCountedPtr<LevelData<EBFluxFAB> > >*        a_lambda,
				Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >* a_etaIrreg,
				Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >* a_lambdaIrreg,
				Vector<RefCountedPtr<LevelData<EBCellFAB> > >*        a_rcv)
{
  m_aco                 = a_aco;
  m_eta                 = a_eta;
  m_lambda              = a_lambda;
  m_etaIrreg            = a_etaIrreg;
  m_lambdaIrreg         = a_lambdaIrreg;
  m_rcv                 = a_rcv;
}
/*****************************/
void EBAMRNavier::defineRefinementBoxes()
{
  
}
void EBAMRNavier::defineSolvers()
{
  if (m_hasDiffusion)
    {
      CH_TIME("EBAMRNavier::defineSolvers");


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
      Real normThresh = 1.0e-30;

      // hardwire ebBcType
      m_params.ebBcType = 1;  // value is read thorugh parmparse
      ppViscous.get("domain_bc_type",m_params.domBcType);
      //temperature Boundary conditions
      ppViscous.get("eb_bc_typeT",m_params.ebBcTypeT);
      ppViscous.get("eb_bc_valueT",m_params.ebBCValueT);
      ppViscous.get("order_ebbc", m_params.orderEB);
      

      Vector<EBAMRNavier*>            hierarchy;
      Vector<DisjointBoxLayout>       grids;
      Vector<EBISLayout>              ebisl;
      Vector<EBLevelGrid>             eblg;
      Vector<int>                     refRat;
      ProblemDomain                   lev0Dom;
      Real                            lev0Dx;
      Vector<ProblemDomain>           domains;
      Vector<Real>                    dxs;
      getHierarchyAndGrids(hierarchy, grids, ebisl, eblg, refRat, lev0Dom, lev0Dx, domains, dxs);

      //set on all levels:: try if it works when defining it on only 1 level
      setViscCoeff(*m_aco, *m_eta, *m_lambda, *m_etaIrreg, *m_lambdaIrreg, *m_rcv, hierarchy, grids, ebisl, eblg, refRat, lev0Dom, dxs, domains);

      getEBVTOFactory(s_viscOpFact, grids, ebisl, eblg, refRat, lev0Dx);
      getConductivityFactory(s_condOpFact, grids, ebisl, eblg, refRat, lev0Dx, domains);


      s_botSolver.m_verbosity = 0;

      { // the EBVTO operator
	RefCountedPtr< AMRLevelOpFactory< LevelData< EBCellFAB > > > OpFact = s_viscOpFact;
	
	s_viscAMRMG = RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB> > > (new AMRMultiGrid<LevelData<EBCellFAB> >());
	
	s_viscAMRMG->define(lev0Dom, *s_viscOpFact, &s_botSolver, hierarchy.size());
	s_viscAMRMG->setSolverParameters(numSmooth, numSmooth, numSmooth,
					 numMG, maxIter, eps, hang, normThresh);
	s_viscAMRMG->m_verbosity = verby-1;
	
	s_viscLevTGA = RefCountedPtr<EBLevelBackwardEuler>
	  (new EBLevelBackwardEuler(grids, refRat, lev0Dom, OpFact, s_viscAMRMG));
	s_viscLevTGA->setEBLG(eblg);
	s_viscLevTGA->setBeta(s_beta);
	
	Vector<AMRLevelOp <LevelData<EBCellFAB> > *> AMROps = s_viscAMRMG->getAMROperators(); 
	s_EBAMROps.resize(AMROps.size());
	for (int ilev = 0; ilev < AMROps.size(); ilev++)
	  {
	    s_EBAMROps[ilev] =  (static_cast <EBViscousTensorOp*> (AMROps[ilev]));
	    s_EBAMROps[ilev]->setAlphaAndBeta(m_alpha, s_beta);
	  }
      }//end ebvto AMRMG setup

      {//same operations for the conductivity operator
	RefCountedPtr< AMRLevelOpFactory< LevelData< EBCellFAB > > > OpFact = s_condOpFact;
	
	s_condAMRMG = RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB> > > (new AMRMultiGrid<LevelData<EBCellFAB> >());
	
	s_condAMRMG->define(lev0Dom, *s_condOpFact, &s_botSolver, hierarchy.size());
	s_condAMRMG->setSolverParameters(numSmooth, numSmooth, numSmooth,
					 numMG, maxIter, eps, hang, normThresh);
	s_condAMRMG->m_verbosity = verby-1;
	
	s_condLevTGA = RefCountedPtr<EBLevelBackwardEuler>
	  (new EBLevelBackwardEuler(grids, refRat, lev0Dom, OpFact, s_condAMRMG));
	s_condLevTGA->setEBLG(eblg);
	s_condLevTGA->setBeta(s_betaCond);
	
	Vector<AMRLevelOp <LevelData<EBCellFAB> > *> AMROps = s_condAMRMG->getAMROperators(); 
	s_EBAMROpsT.resize(AMROps.size());
	for (int ilev = 0; ilev < AMROps.size(); ilev++)
	  {
	    s_EBAMROpsT[ilev] =  (static_cast <EBConductivityOp*> (AMROps[ilev]));
	    s_EBAMROpsT[ilev]->setAlphaAndBeta(m_alpha, s_betaCond);
	  }
      }//end conductivity AMRMG setup
    }

}

/*******/
void
EBAMRNavier::
getHierarchyAndGrids(Vector<EBAMRNavier*>&        a_hierarchy,
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

  EBAMRNavier* coarsestLevel = (EBAMRNavier*)(hierarchy[0]);
  a_lev0Dx       = coarsestLevel->m_dx[0];
  a_lev0Dom      = coarsestLevel->m_problem_domain;


  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();

  for (int ilev = 0; ilev < nlevels; ilev++)
    {
      EBAMRNavier* adLevel = (EBAMRNavier*)(hierarchy[ilev]);

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

/**/
void EBAMRNavier::
getEBVTOFactory(RefCountedPtr<EBViscousTensorOpFactory>&                    a_factory,
                const Vector<DisjointBoxLayout>&                            a_grids,
                const Vector<EBISLayout>&                                   a_ebisl,
                const Vector<EBLevelGrid>&                                  a_eblg,
		Vector<int>&                                                a_refRat,
		Real&                                                       a_lev0Dx)
{ 
  RefCountedPtr<BaseDomainBCFactory> domBC; // either zero Neumann or value Dirichlet
  RefCountedPtr<BaseEBBCFactory>     ebBC; //  zero Dirichlet
  getViscousBCFactories(domBC, ebBC, a_grids, a_ebisl);


  a_factory = RefCountedPtr<EBViscousTensorOpFactory>
    (new EBViscousTensorOpFactory(a_eblg, m_alpha, s_beta, *m_aco, *m_eta, *m_lambda, *m_etaIrreg, *m_lambdaIrreg,
                                  a_lev0Dx, a_refRat, domBC, ebBC,
                                  m_nGhost*IntVect::Unit, m_nGhost*IntVect::Unit, -1, false));
}
/**/
void EBAMRNavier::
getViscousBCFactories(RefCountedPtr<BaseDomainBCFactory>&         a_domBC,
		      RefCountedPtr<BaseEBBCFactory>&             a_ebBC,
		      const Vector<DisjointBoxLayout>&            a_grids,
		      const Vector<EBISLayout>&                   a_ebisl)
{

  ParmParse ppViscousBC("ns");
  if (m_params.domBcType == 0)
    {
      Real domBCValue =0;  

      NeumannViscousTensorDomainBCFactory* domainBCFactory = new NeumannViscousTensorDomainBCFactory();
      domainBCFactory->setValue(domBCValue);

      a_domBC = RefCountedPtr<BaseDomainBCFactory>(domainBCFactory);
    }
  else if (m_params.domBcType == 1)
    {
      DirichletViscousTensorDomainBCFactory* domainBCFactory = new DirichletViscousTensorDomainBCFactory();
      Real domBCValue;
      ppViscousBC.get("domain_bc_value", domBCValue);
      domainBCFactory->setValue(domBCValue);

      a_domBC = RefCountedPtr<BaseDomainBCFactory>(domainBCFactory);
    }
  else if (m_params.domBcType == 2)
    {
      Real domBCValue =0;  // hardwire these two
      Vector<Vector<Vector<int> > > isDirichlet;
      isDirichlet.resize(2);
      for (int j = 0; j < 2; j++)
	{
	  isDirichlet[j].resize(SpaceDim);
	  for (int i = 0; i < SpaceDim; i++) {isDirichlet[j][i].resize(SpaceDim,0);}
	}
      if (SpaceDim >= 3) isDirichlet[1][2][2] = 1;
      
      MixedViscousTensorDomainBCFactory* domainBCFactory = new MixedViscousTensorDomainBCFactory();
      domainBCFactory->setValue(domBCValue);
      domainBCFactory->setType(isDirichlet);

      a_domBC = RefCountedPtr<BaseDomainBCFactory>(domainBCFactory);
    }
  else
    {
      MayDay::Error("Unknown domain BC type");
    }


  Real ebBCValue = 0;
  DirichletViscousTensorEBBCFactory* ebBC = new DirichletViscousTensorEBBCFactory();

  ebBC->setValue(ebBCValue);
  m_isFunctionEBVTO = false;

  a_ebBC = RefCountedPtr<BaseEBBCFactory>(ebBC);
}
/**/
void EBAMRNavier::
getConductivityFactory(RefCountedPtr<EBConductivityOpFactory>&                    a_factory,
		       const Vector<DisjointBoxLayout>&                            a_grids,
		       const Vector<EBISLayout>&                                   a_ebisl,
		       const Vector<EBLevelGrid>&                                  a_eblg,
		       Vector<int>&                                                a_refRat,
		       Real&                                                       a_lev0Dx,
		       Vector<ProblemDomain>&                                      a_domains)
{ 
  RefCountedPtr<BaseDomainBCFactory> domBC; // either zero Neumann or value Dirichlet
  RefCountedPtr<BaseEBBCFactory>     ebBC; //  zero Dirichlet
  getConductivityBCFactories(domBC, ebBC, a_grids, a_ebisl);
  
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
  int relaxType = 1;
  ppViscous.get("relax_type", relaxType);
  a_factory = RefCountedPtr<EBConductivityOpFactory>
    (new EBConductivityOpFactory(a_eblg, quadCFI, m_alpha, s_betaCond, *m_rcv, *m_eta, *m_etaIrreg,
                                 a_lev0Dx,  a_refRat, domBC, ebBC,
                                 m_nGhost*IntVect::Unit, m_nGhost*IntVect::Unit, relaxType));
}

void EBAMRNavier::
getConductivityBCFactories(RefCountedPtr<BaseDomainBCFactory>&         a_domBC,
			   RefCountedPtr<BaseEBBCFactory>&             a_ebBC,
			   const Vector<DisjointBoxLayout>&            a_grids,
			   const Vector<EBISLayout>&                   a_ebisl)
{
  CH_TIME("EBAMRViscous::getConductivityBCFactories");


  Real domBCValueT = 0;
  NeumannConductivityDomainBCFactory* domainBCFactory = new NeumannConductivityDomainBCFactory();
  domainBCFactory->setValue(domBCValueT);
  a_domBC = RefCountedPtr<BaseDomainBCFactory>(domainBCFactory);

  // Note:: I did not implement mixed boundary conditions so if using injection force to use dirichlet
  //if (m_isFunctionEBVTO && !m_params.ebBcTypeT == 1)
  //  MayDay::Error("did not implement mixed boundary conditionsso use Dirichlet T with injection");

  if (m_params.ebBcTypeT == 0)
    {
      NeumannConductivityEBBCFactory* ebBC = new NeumannConductivityEBBCFactory();
      ebBC->setValue(m_params.ebBCValueT);

      a_ebBC =RefCountedPtr<BaseEBBCFactory>( ebBC);
    }
  else if (m_params.ebBcTypeT == 1)
    {
      MayDay::Error("Dirichlet temperature not implemented");
    }
  else
    {
      MayDay::Error("Unknown EB BC type Temperature");
    }
}
void
EBAMRNavier::
makeDiffusiveSource(LevelData<EBCellFAB>& a_diffusiveSrc, const bool& incFlag)
//incflag defaults to false
{
  CH_TIME("EBAMRNavier::makeDiffusiveSource");
  
  if(!incFlag)
    for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
      a_diffusiveSrc[dit()].setVal(0);

  if (m_hasDiffusion)
    {
      Interval momIntrv = m_ebPolytropic->momentumInterval();
      Interval velIntrv = Interval(0, SpaceDim-1);
      Interval nrgIntrv = Interval(m_ebPolytropic->energyIndexC(), m_ebPolytropic->energyIndexC());

      //note beta is -1/Re
      for (int ilev = 0; ilev < s_EBAMROps.size(); ilev++) 
	{
	  s_EBAMROps[ilev]->setAlphaAndBeta(0., s_beta);
	  s_EBAMROpsT[ilev]->setAlphaAndBeta(0., s_betaCond);
	}


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
      EBCellFactory factory(m_ebisl);
      LevelData<EBCellFAB>* stateVCoarse;
      LevelData<EBCellFAB>* stateTCoarse;
      LevelData<EBCellFAB>* stateVFine;
      LevelData<EBCellFAB>* stateTFine;

      //note using m_stateNew rather than Old
      LevelData<EBCellFAB>& stateU  = m_stateNew;
      LevelData<EBCellFAB> stateV(m_grids, SpaceDim, m_stateOld.ghostVect(), factory);
      LevelData<EBCellFAB> stateT(m_grids, 1       , m_stateOld.ghostVect(), factory);
      LevelData<EBCellFAB> lhs(m_grids, SpaceDim, m_stateOld.ghostVect(), factory);
      LevelData<EBCellFAB> lhsT(m_grids, 1      , m_stateOld.ghostVect(), factory);
      

      Real interpTime = m_time + 0.5*m_dt;
      if (m_hasCoarser)
	{
	  Real CinterpTime = Min(interpTime,tCoarserNew);
	  EBAMRNavier* coarserPtr = getCoarserLevel();
	  EBCellFactory Cfactory(coarserPtr->m_ebisl);
	  LevelData<EBCellFAB>* stateUCoarse = new LevelData<EBCellFAB>(coarserDataOldPtr->disjointBoxLayout(),coarserDataOldPtr->nComp(), coarserDataOldPtr->ghostVect(), Cfactory);
	  interpolateInTime(*stateUCoarse, *coarserDataOldPtr, *coarserDataNewPtr,coarserPtr->m_grids,CinterpTime, tCoarserOld, tCoarserNew);
	  //new stuff:: patcher
	  int refRatCrse = coarserPtr->refRatio();
	  EBPWLFillPatch patcher(m_grids,
				 coarserPtr->m_grids,
				 m_ebisl,
				 coarserPtr->m_ebisl,
				 coarserPtr->m_domainBox,
				 refRatCrse, m_nComp, m_nGhost);
	  patcher.interpolate(stateU,
			      *coarserDataOldPtr, 
			      *coarserDataNewPtr,
			      tCoarserOld,
			      tCoarserNew,
			      CinterpTime,
			      Interval(0, m_nComp-1));
      //get vel and T
	  stateVCoarse = new LevelData<EBCellFAB>(stateUCoarse->disjointBoxLayout(),SpaceDim, stateUCoarse->ghostVect(), Cfactory);
	  stateTCoarse = new LevelData<EBCellFAB>(stateUCoarse->disjointBoxLayout(),1       , stateUCoarse->ghostVect(), Cfactory);
	  consToVelTemp(*stateVCoarse,*stateTCoarse,*stateUCoarse);
	  delete stateUCoarse;
	  EBLevelDataOps::setCoveredVal(*stateVCoarse,0e0);
	  EBLevelDataOps::setCoveredVal(*stateTCoarse,0e0);
	}
      if (m_hasFinerLevel)
	{
	  Real FinterpTime = Min(interpTime,tFinerNew);
	  EBAMRNavier* finerPtr = getFinerLevel();
	  EBCellFactory Ffactory(finerPtr->m_ebisl);
	  LevelData<EBCellFAB>* stateUFine = new LevelData<EBCellFAB>(finerDataOldPtr->disjointBoxLayout(),finerDataOldPtr->nComp(), finerDataOldPtr->ghostVect(), Ffactory);
	  interpolateInTime(*stateUFine, *finerDataOldPtr, *finerDataNewPtr,finerPtr->m_grids,FinterpTime, tFinerOld, tFinerNew);
	  stateVFine = new LevelData<EBCellFAB>(stateUFine->disjointBoxLayout(),SpaceDim, stateUFine->ghostVect(), Ffactory);
	  stateTFine = new LevelData<EBCellFAB>(stateUFine->disjointBoxLayout(),1       , stateUFine->ghostVect(), Ffactory);
	  consToVelTemp(*stateVFine,*stateTFine,*stateUFine);
	  delete stateUFine;
	  EBLevelDataOps::setCoveredVal(*stateVFine,0e0);
	  EBLevelDataOps::setCoveredVal(*stateTFine,0e0);
	}

      consToVelTemp(stateV,stateT,stateU);
      EBLevelDataOps::setCoveredVal(stateV,0e0);
      EBLevelDataOps::setCoveredVal(stateT,0e0);
      EBLevelDataOps::setVal(lhs,0e0);
      EBLevelDataOps::setVal(lhsT,0e0);


      if (!m_hasCoarser && !m_hasFinerLevel)
        {
          s_EBAMROps[m_level]->applyOp (lhs, stateV, false);
	  s_EBAMROpsT[m_level]->applyOp(lhsT,stateT, false);
        }
      else if (!m_hasCoarser)
        {
          s_EBAMROps[m_level]->AMROperatorNC (lhs, *stateVFine, stateV, false, s_EBAMROps[m_level+1]);
          s_EBAMROpsT[m_level]->AMROperatorNC(lhsT,*stateTFine, stateT, false, s_EBAMROpsT[m_level+1]);
        }
      else if(!m_hasFinerLevel) //note the difference btwn m_hasFiner<< ", " << m_hasFinerLevel (checks grid>0);
        {
          s_EBAMROps[m_level]->AMROperatorNF (lhs, stateV, *stateVCoarse, false);
          s_EBAMROpsT[m_level]->AMROperatorNF(lhsT,stateT, *stateTCoarse, false);
        }
      else
	{
	  s_EBAMROps[m_level]->AMROperator (lhs, *stateVFine, stateV, *stateVCoarse, false, s_EBAMROps[m_level+1]);
	  s_EBAMROpsT[m_level]->AMROperator(lhsT,*stateTFine, stateT, *stateTCoarse, false, s_EBAMROpsT[m_level+1]);
	}


      if(1){
	LevelData<EBCellFAB> dsu(m_grids, 1, m_stateOld.ghostVect(), factory);
	s_EBAMROps[m_level]->getKappaDivSigmaU(dsu,stateV,stateVCoarse,m_level);
	EBLevelDataOps::incr(lhsT,dsu,s_beta); // scale by -1/Re and add to the LHS
      }


      if(!incFlag)
	{
	  EBLevelDataOps::assign(a_diffusiveSrc,lhs, momIntrv,velIntrv);
	  EBLevelDataOps::assign(a_diffusiveSrc,lhsT,nrgIntrv,lhsT.interval());
	  //divide by the volume fraction
	  EBLevelDataOps::kappaDivide(a_diffusiveSrc);
	  // bring on the RHS
	  EBLevelDataOps::scale(a_diffusiveSrc,-1.0);
	}
      else
	{
	  EBLevelDataOps::kappaDivide(lhs);
	  EBLevelDataOps::kappaDivide(lhsT);
	  EBLevelDataOps::incr(a_diffusiveSrc,lhs, momIntrv.begin(),0,SpaceDim,-1.0);
	  EBLevelDataOps::incr(a_diffusiveSrc,lhsT,nrgIntrv.begin(),0,1,-1.0);
	}
	  

      /// Over the coarse-fine interface, the diffusive source is set
      /// to zero. At fine-fine interfaces, it is filled in by
      /// neighboring boxes.
      a_diffusiveSrc.exchange();

      if(m_hasFinerLevel) {delete stateVFine;delete stateTFine;}
      if(m_hasCoarser) {delete stateVCoarse;delete stateTCoarse;}

      
      /*pout() << "NAN Check diffusiveSrc " << m_level << endl;
      bool isNAN = EBLevelDataOps::checkNANINF(a_diffusiveSrc);
      if(isNAN) MayDay::Error("Checking diffusiveSrc.");
      pout() << "ByeBye";exit(0);*/
      
    }
}
/******************************************/
void 
EBAMRNavier::implicitDiffusion(EBFluxRegister*& finerFRPtr, EBFluxRegister*& coarserFRPtr, 
			       const LevelData<EBCellFAB>& oldUCoarse, const LevelData<EBCellFAB>& stateUCoarse,
			       Real& tCoarserOld, Real& tCoarserNew)
{

  Interval momIntrv = m_ebPolytropic->momentumInterval();
  Interval nrgIntrv = Interval(m_ebPolytropic->energyIndexC(), m_ebPolytropic->energyIndexC());
 
  Interval intV(0, SpaceDim-1);
  Interval intT(0, 0);
  EBFluxRegister* NfinerFRPtr=NULL; EBFluxRegister* NcoarserFRPtr=NULL;

  LevelData<EBCellFAB>& stateU  = m_stateNew;
  consToVelTemp(m_initialV,m_initialT,stateU);
  EBLevelDataOps::setCoveredVal(m_initialV,0e0);
  EBLevelDataOps::setCoveredVal(m_initialT,0e0);

  /*
  copyLevel(m_srcV,m_stateNew,intV,momIntrv);//srcV and srcT are on the RHS
  copyLevel(m_srcT,m_stateNew,intT,nrgIntrv);
  for (DataIterator dit=m_srcV.dataIterator(); dit.ok(); ++dit)
    {
      m_srcV[dit()].minus(stateU[dit()],momIntrv.begin(),intV.begin(),intV.size());
      m_srcV[dit()] /= m_dt;
      m_srcT[dit()].minus(stateU[dit()],nrgIntrv.begin(),intT.begin(),intT.size());
      m_srcT[dit()] /= m_dt;
    }
  */

  EBLevelDataOps::setToZero(m_srcV);
  EBLevelDataOps::setToZero(m_srcT);    

  // stateV = Velocity(m_stateOld); updateV  = Velocity(m_stateOld); dV = (Momentum(m_stateNew)-Momentum(m_stateOld))/dt
  // same for stateT with energy and temperature

  if (m_hasCoarser)
    {
      consToVelTemp(m_newVCoarse,m_newTCoarse,stateUCoarse);
      consToVelTemp(m_oldVCoarse,m_oldTCoarse,oldUCoarse);
      EBLevelDataOps::setCoveredVal(m_newVCoarse,0e0);
      EBLevelDataOps::setCoveredVal(m_newTCoarse,0e0);
      EBLevelDataOps::setCoveredVal(m_oldVCoarse,0e0);
      EBLevelDataOps::setCoveredVal(m_oldTCoarse,0e0);
    }
  
  //! Integrates the helmholtz equation represented by this object, placing
  //! the new solution in \a a_phiNew.
  //! \param a_phiNew The new solution (the value of phi at time n + 1) will
  //!                 be stored here.
  //! \param a_phiOld The old solution (the value of phi at time n).
  //! \param a_src The source term on the right hand side of the Helmholtz
  //!              equation.
  //! \param a_flux This will store the flux computed at the current grid
  //!               level during the solution of the Helmholtz equation.
  //! \param a_fineFluxRegPtr A pointer to the flux register representing the
  //!                         finer grid level adjacent to this one, or NULL
  //!                         if there is no finer grid level.
  //! \param a_crseFluxRegPtr A pointer to the flux register representing the
  //!                         coarser grid level adjacent to this one, or NULL
  //!                         if there is no coarser grid level.
  //! \param a_oldTime The time at the beginning of the integration step at
  //!                  the current grid level.
  //! \param a_crseOldTime The time at the beginning of the integration step
  //!                      at the coarser adjacent grid level. This parameter
  //!                      is ignored if there is no coarser grid level.
  //! \param a_crseNewTime The time at the end of the integration step
  //!                      at the coarser adjacent grid level. This parameter
  //!                      is ignored if there is no coarser grid level.
  //! \param a_dt The size of the integration step at the current grid level.
  //! \param a_level The current grid level.
  //! \param a_zeroPhi If set to true, \a a_phiNew will be set to zero before
  //!                  the integration takes place. Otherwise, a_phiNew is
  //!                  assumed to be an initial estimate for the solution in
  //!                  the iterative linear solve.
  //! \param a_fluxStartComponent An index identifying the component at which
  //!                             flux data begins within \a a_fineFluxRegPtr
  //!                             and \a a_crseFluxRegPtr.

  //! updateSoln (LevelData< EBCellFAB > &a_phiNew, LevelData< EBCellFAB > &a_phiOld, LevelData< EBCellFAB > &a_src, EBFluxRegister *a_fineFluxRegPtr, EBFluxRegister *a_crseFluxRegPtr, const LevelData< EBCellFAB > *a_crsePhiOldPtr, const LevelData< EBCellFAB > *a_crsePhiNewPtr, Real a_oldTime, Real a_crseOldTime, Real a_crseNewTime, Real a_dt, int a_level, bool a_zeroPhi=true, bool a_rhsAlreadyKappaWeighted=false, int a_fluxStartComponent=0)
  copyLevel(m_updateV, m_initialV, intV, intV);
  s_viscLevTGA->updateSoln(m_updateV, m_initialV, m_srcV,
			   NfinerFRPtr, NcoarserFRPtr,  //lucacancel delete N
			   &m_oldVCoarse, &m_newVCoarse,
			   m_time, tCoarserOld, tCoarserNew,
			   m_dt, m_level, false, false, momIntrv.begin());

  LevelData<EBCellFAB>& dsu = m_updateT;
  EBLevelDataOps::setToZero(dsu);
  s_EBAMROps[m_level]->getKappaDivSigmaU(dsu,m_updateV,&m_newVCoarse,m_level);
  EBLevelDataOps::incr(m_srcT,dsu,-s_beta); // scale by 1/Re and add to srcT (RHS)

  addDKEDt(m_srcT,m_updateV, m_initialV,m_stateNew);
  copyLevel(m_updateT, m_initialT, intT, intT);
  s_condLevTGA->updateSoln(m_updateT, m_initialT, m_srcT,
			   NfinerFRPtr, NcoarserFRPtr, //lucacancel delete N
			   &m_oldTCoarse, &m_newTCoarse,
			   m_time, tCoarserOld, tCoarserNew,
			   m_dt, m_level, false, false, nrgIntrv.begin());
  updateConserved(m_stateNew,m_updateV,m_updateT);
}
/********/
void
EBAMRNavier::
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
      EBAMRNavier* coarserPtr = getCoarserLevel();

      *a_coarserDataOldPtr = &coarserPtr->m_stateOld;
      *a_coarserDataNewPtr = &coarserPtr->m_stateNew;

      a_tCoarserNew = Max(coarserPtr->m_time,0e0);
      a_tCoarserOld = Max(a_tCoarserNew - coarserPtr->m_dt,0e0);
    }

  // A finer level exists
  if (m_hasFiner)
    {
      EBAMRNavier* finerPtr = getFinerLevel();

      *a_finerDataOldPtr = &finerPtr->m_stateOld;
      *a_finerDataNewPtr = &finerPtr->m_stateNew;

      a_tFinerNew = Max(finerPtr->m_time,0e0);
      a_tFinerOld = Max(a_tFinerNew - finerPtr->m_dt,0e0);
    }
}
/*
void 
EBAMRNavier::consToVel(LevelData<EBCellFAB>& pa_V, const LevelData<EBCellFAB>& pa_consState, const LevelData<EBCellFAB>& a_consState)
{
  for (DataIterator dit = a_consState.dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB& ebV1     = pa_V[dit()];
      const EBCellFAB& ebConst1 = pa_consState[dit()];
      const EBCellFAB& ebConst = a_consState[dit()];

      Box intBox = ebV1.getRegion() & ebConst1.getRegion();
      consToVel(ebV1, ebConst1, ebConst, intBox);
    }
}
void 
EBAMRNavier::consToVel(LevelData<EBCellFAB>& pa_V, LevelData<EBCellFAB>& a_V, const LevelData<EBCellFAB>& pa_consState, const LevelData<EBCellFAB>& a_consState)
{
  for (DataIterator dit = a_consState.dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB& ebV1     = pa_V[dit()];
      EBCellFAB& ebV      =  a_V[dit()];
      const EBCellFAB& ebConst1 = pa_consState[dit()];
      const EBCellFAB& ebConst  =  a_consState[dit()];

      Box intBox = ebV1.getRegion() & ebConst1.getRegion();
      consToVel(ebV1, ebV, ebConst1, ebConst, intBox);
    }
}
void 
EBAMRNavier::consToVel(EBCellFAB& pa_V, const EBCellFAB& pa_consState,  const EBCellFAB& a_consState, const Box& a_box)
{
  BaseFab<Real>&         regV1 =    pa_V.getSingleValuedFAB();
  const BaseFab<Real>&   regCons1 = pa_consState.getSingleValuedFAB();
  const BaseFab<Real>&   regCons =  a_consState.getSingleValuedFAB();


  FORT_PGETVELOCITY(CHF_FRA(regV1),CHF_CONST_FRA(regCons1),CHF_CONST_FRA(regCons),CHF_BOX(a_box)); 

  const EBISBox& ebisBox = pa_consState.getEBISBox();
  IntVectSet ivsMulti = ebisBox.getMultiCells(a_box);
  if(!ivsMulti.isEmpty())
    {
      Interval momIntrv = m_ebPolytropic->momentumInterval();
      Interval velIntrv = Interval(0, SpaceDim-1);
      int densIndx = m_ebPolytropic->densityIndexC();
    
      for (VoFIterator vofit(ivsMulti, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
	{
	  const VolIndex& vof = vofit();

	  for (int ivar = 0; ivar <= velIntrv.end(); ivar++)
	    {
	      Real dens = a_consState(vof,densIndx);
	      Real dens2 = dens*dens;
	      int  mvar = momIntrv.begin()+ivar;
	      pa_V(vof,ivar) = pa_consState(vof,mvar)/dens - 
		a_consState(vof,mvar)/dens2 * pa_consState(vof,densIndx) ;	
	    }  
	}
    }
}
void 
EBAMRNavier::consToVel(EBCellFAB& pa_V, EBCellFAB& a_V, const EBCellFAB& pa_consState,  const EBCellFAB& a_consState, const Box& a_box)
{
  BaseFab<Real>&         regV1 =    pa_V.getSingleValuedFAB();
  BaseFab<Real>&         regV  =     a_V.getSingleValuedFAB();
  const BaseFab<Real>&   regCons1 = pa_consState.getSingleValuedFAB();
  const BaseFab<Real>&   regCons =  a_consState.getSingleValuedFAB();


  FORT_PGETVELOCITYBOTH(CHF_FRA(regV1),CHF_FRA(regV),CHF_CONST_FRA(regCons1),CHF_CONST_FRA(regCons),CHF_BOX(a_box)); 

  const EBISBox& ebisBox = pa_consState.getEBISBox();
  IntVectSet ivsMulti = ebisBox.getMultiCells(a_box);
  if(!ivsMulti.isEmpty())
    {
      Interval momIntrv = m_ebPolytropic->momentumInterval();
      Interval velIntrv = Interval(0, SpaceDim-1);
      int densIndx = m_ebPolytropic->densityIndexC();
    
      for (VoFIterator vofit(ivsMulti, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
	{
	  const VolIndex& vof = vofit();

	  for (int ivar = 0; ivar <= velIntrv.end(); ivar++)
	    {
	      Real dens = a_consState(vof,densIndx);
	      Real dens2 = dens*dens;
	      int  mvar = momIntrv.begin()+ivar;
	      pa_V(vof,ivar) = pa_consState(vof,mvar)/dens - 
		a_consState(vof,mvar)/dens2 * pa_consState(vof,densIndx) ;
	      a_V(vof,ivar) = a_consState(vof,mvar)/dens;	
	    }  
	}
    }
}
*/
void 
EBAMRNavier::consToVelTemp(LevelData<EBCellFAB>& a_V, LevelData<EBCellFAB>& a_T, const LevelData<EBCellFAB>& a_consState)
{
  
  EBLevelDataOps::setToZero(a_V);
  EBLevelDataOps::setToZero(a_T);
  for (DataIterator dit = a_consState.dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB& ebV     = a_V[dit()];
      EBCellFAB& ebT     = a_T[dit()];
      const EBCellFAB& ebConst = a_consState[dit()];

      Box intBox = ebV.getRegion() & ebConst.getRegion();
      consToVelTemp(ebV, ebT, ebConst, intBox);
    }
}
void 
EBAMRNavier::consToVelTemp(EBCellFAB& a_V, EBCellFAB& a_T, const EBCellFAB& a_consState, const Box& a_box)
{
  BaseFab<Real>&         regV =    a_V.getSingleValuedFAB();
  BaseFab<Real>&         regT =    a_T.getSingleValuedFAB();
  const BaseFab<Real>&   regCons =  a_consState.getSingleValuedFAB();


  FORT_GETVELNTEMP(CHF_FRA(regV),CHF_FRA1(regT,0),CHF_CONST_FRA(regCons),CHF_BOX(a_box),CHF_CONST_REAL(m_cv)); 

  const EBISBox& ebisBox = a_consState.getEBISBox();
  IntVectSet ivsMulti = ebisBox.getMultiCells(a_box);
  if(!ivsMulti.isEmpty())
    {
      Interval momIntrv = m_ebPolytropic->momentumInterval();
      Interval velIntrv = Interval(0, SpaceDim-1);
      int densIndx = m_ebPolytropic->densityIndexC();
      int nrgIndx  = m_ebPolytropic->energyIndexC();
    
      for (VoFIterator vofit(ivsMulti, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
	{
	  const VolIndex& vof = vofit();

	  Real dens = a_consState(vof,densIndx);
	  Real kinetic = 0e0;
	  for (int ivar = 0; ivar <= velIntrv.end(); ivar++)
	    {
	      int mvar = momIntrv.begin()+ivar;
	      Real vel = a_consState(vof,mvar)/dens;
	      a_V(vof,ivar) = vel;
	      kinetic += vel*vel;
	    }
	  kinetic /= 2e0;
	  a_T(vof,0) = (a_consState(vof,nrgIndx)/dens  - kinetic)/m_cv;
	}
    }
}
//
void 
EBAMRNavier::consToDiffusiveProps(LevelData<EBCellFAB>& a_R, LevelData<EBCellFAB>& a_E, LevelData<EBCellFAB>& a_RCV, const LevelData<EBCellFAB>& a_consState, const ProblemDomain& a_domain)
{
  for (DataIterator dit = a_consState.dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB& ebR     = a_R[dit()];
      EBCellFAB& ebE     = a_E[dit()];
      EBCellFAB& ebRCV   = a_RCV[dit()];
      const EBCellFAB& ebConst = a_consState[dit()];

      Box intBox = ebR.getRegion() & ebConst.getRegion();
      // Box intBox = grow(a_consState.box(dit()),1); intBox &= a_domain;
      
      
      BaseFab<Real>&       regR =    ebR.getSingleValuedFAB();
      BaseFab<Real>&       regRCV =  ebRCV.getSingleValuedFAB();
      BaseFab<Real>&       regE =    ebE.getSingleValuedFAB();
      const BaseFab<Real>& regCons = ebConst.getSingleValuedFAB();


      FORT_GETDENSETARCV(CHF_FRA1(regR,0),CHF_FRA1(regE,0),CHF_FRA1(regRCV,0),CHF_CONST_FRA(regCons),
			  CHF_BOX(intBox),CHF_CONST_REAL(m_cv),CHF_CONST_REAL(m_Tinf));

      const EBISBox& ebisBox = ebConst.getEBISBox();
      IntVectSet ivsMulti = ebisBox.getMultiCells(intBox);
      if(!ivsMulti.isEmpty())
	{
	  Interval momIntrv = m_ebPolytropic->momentumInterval();
	  Interval velIntrv = Interval(0, SpaceDim-1);
	  int densIndx = m_ebPolytropic->densityIndexC();
	  int nrgIndx  = m_ebPolytropic->energyIndexC();
	  Real rho = 1; // not used
	  Real etainf;
	  FORT_POINTAIRVISCOSITY(CHF_REAL(etainf),CHF_CONST_REAL(rho),CHF_CONST_REAL(m_Tinf));
    
	  for (VoFIterator vofit(ivsMulti, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
	    {
	      const VolIndex& vof = vofit();
	      
	      Real dens = ebConst(vof,densIndx);
	      Real kinetic = 0e0;
	      for (int ivar = 0; ivar <= velIntrv.end(); ivar++)
		{
		  int mvar = momIntrv.begin()+ivar;
		  Real vel = ebConst(vof,mvar)/dens;
		  kinetic += vel*vel;
		}
	      kinetic /= 2e0;
	      Real Temp = (ebConst(vof,nrgIndx)/dens  - kinetic)/m_cv;
	      Real Tdim = Temp*m_Tinf;
	      Real etaDim;
	      FORT_POINTAIRVISCOSITY(CHF_REAL(etaDim),CHF_CONST_REAL(rho),CHF_CONST_REAL(Tdim));
	      ebR(vof,0) = dens;
	      ebE(vof,0) = Min(Max(etaDim/etainf,0.33),3.0);
	      ebRCV(vof,0) = dens*m_cv;
	    }
	}
    }
}
void 
EBAMRNavier::updateConserved(LevelData<EBCellFAB>& a_consState, const LevelData<EBCellFAB>& a_V, const LevelData<EBCellFAB>& a_T)
{
  for (DataIterator dit = a_consState.dataIterator(); dit.ok(); ++dit)
    {
      const EBCellFAB& ebV     = a_V[dit()];
      const EBCellFAB& ebT     = a_T[dit()];
      EBCellFAB& ebConst = a_consState[dit()];

      Box intBox = ebV.getRegion() & ebConst.getRegion();


      const BaseFab<Real>& regV =   ebV.getSingleValuedFAB();
      const BaseFab<Real>& regT =   ebT.getSingleValuedFAB();
      BaseFab<Real>&       regCons= ebConst.getSingleValuedFAB();


      FORT_UPDATECONS(CHF_FRA(regCons),CHF_CONST_FRA(regV),CHF_CONST_FRA1(regT,0),CHF_BOX(intBox),CHF_CONST_REAL(m_cv)); 

      const EBISBox& ebisBox = ebConst.getEBISBox();
      IntVectSet ivsMulti = ebisBox.getMultiCells(intBox);
      if(!ivsMulti.isEmpty())
	{
	  Interval momIntrv = m_ebPolytropic->momentumInterval();
	  Interval velIntrv = Interval(0, SpaceDim-1);
	  int densIndx = m_ebPolytropic->densityIndexC();
	  int nrgIndx  = m_ebPolytropic->energyIndexC();
    
	  for (VoFIterator vofit(ivsMulti, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
	    {
	      const VolIndex& vof = vofit();

	      Real dens = ebConst(vof,densIndx);
	      Real kinetic = 0e0;
	      for (int ivar = 0; ivar <= velIntrv.end(); ivar++)
		{
		  int mvar = momIntrv.begin()+ivar;
		  Real vel = ebV(vof,ivar);
		  ebConst(vof,ivar) = dens*vel;
		  kinetic += vel*vel;
		}
	      kinetic /= 2e0;
	      ebConst(vof,nrgIndx) = dens*(m_cv*ebT(vof,0)+kinetic);
	    }
	}
      
    }
}

void EBAMRNavier::setViscCoeff(Vector<RefCountedPtr<LevelData<EBCellFAB> > >&        a_aco,
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
			       Vector<ProblemDomain>&       a_domains)
{ 
  int nlevels = a_hierarchy.size();
  int ghost_eta = 4;
  for (int ilev = 0; ilev < nlevels; ilev++)
    {
      LevelData<EBCellFAB>& stateU  = a_hierarchy[ilev]->getStateNew();
      
      EBCellFactory factory(a_ebisl[ilev]);
      int nghost = m_nGhost;
      int nCons = m_ebPatchGodunov->numConserved();
      Interval consInterv(0, nCons-1);
      LevelData<EBCellFAB> consTemp(a_grids[ilev], nCons, (IntVect::Unit)*nghost, factory);
      stateU.copyTo(consInterv, consTemp, consInterv);
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
	  patcher.interpolate(consTemp,
			      a_hierarchy[ilev-1]->m_stateOld,
			      a_hierarchy[ilev-1]->m_stateNew,
			      coarTimeOld,
			      coarTimeNew,
			      fineTime,
			      consInterv);
	}
      consTemp.exchange(consInterv);

      LayoutData<IntVectSet> irregSets(a_grids[ilev]);

      for (DataIterator dit = consTemp.dataIterator(); dit.ok(); ++dit) 
	{
	  Box grownBox = grow(a_grids[ilev].get(dit()), ghost_eta);
	  grownBox &= a_domains[ilev];
	  irregSets[dit()] = a_ebisl[ilev][dit()].getIrregIVS(grownBox);
	}
      EBCellFactory        ebcellfact(a_ebisl[ilev]);
      EBFluxFactory        ebfluxfact(a_ebisl[ilev]);
      BaseIVFactory<Real>  baseivfact(a_ebisl[ilev], irregSets);
      
      a_lambda[ilev]= RefCountedPtr< LevelData<EBFluxFAB> >(new LevelData<EBFluxFAB>(a_grids[ilev], 1, ghost_eta*IntVect::Unit, ebfluxfact));
      a_eta[ilev]   = RefCountedPtr< LevelData<EBFluxFAB> >(new LevelData<EBFluxFAB>(a_grids[ilev], 1, ghost_eta*IntVect::Unit, ebfluxfact));
      a_aco[ilev]   = RefCountedPtr< LevelData<EBCellFAB> >(new LevelData<EBCellFAB>(a_grids[ilev], 1, ghost_eta*IntVect::Unit, ebcellfact));
      a_rcv[ilev]   = RefCountedPtr< LevelData<EBCellFAB> >(new LevelData<EBCellFAB>(a_grids[ilev], 1, ghost_eta*IntVect::Unit, ebcellfact));
      a_etaIrreg[ilev]    = RefCountedPtr< LevelData< BaseIVFAB<Real> > >(new LevelData< BaseIVFAB<Real> >(a_grids[ilev], 1, ghost_eta*IntVect::Unit, baseivfact));
      a_lambdaIrreg[ilev] = RefCountedPtr< LevelData< BaseIVFAB<Real> > >(new LevelData< BaseIVFAB<Real> >(a_grids[ilev], 1, ghost_eta*IntVect::Unit, baseivfact));
      
      LevelData<EBCellFAB>&    aco =  *a_aco[ilev];
      LevelData<EBFluxFAB>&    eta =  *a_eta[ilev];
      LevelData<EBFluxFAB>&    lambda=*a_lambda[ilev];
      LevelData<EBCellFAB>&    rcv =  *a_rcv[ilev];
      LevelData< BaseIVFAB<Real> >&  etaIrreg    = *a_etaIrreg[ilev];
      LevelData< BaseIVFAB<Real> >&  lambdaIrreg = *a_lambdaIrreg[ilev];
      LevelData<EBCellFAB> etacell(a_grids[ilev], 1, ghost_eta*IntVect::Unit, ebcellfact);

      EBLevelDataOps::setVal(etacell,1e0);
      EBLevelDataOps::setVal(eta,1e0);
      EBLevelDataOps::setVal(lambda,0e0);
      EBLevelDataOps::setVal(aco,m_gamma);
      EBLevelDataOps::setVal(rcv,m_gamma*m_cv);
      for (DataIterator dit = eta.dataIterator(); dit.ok(); ++dit)
	{
	  etaIrreg[dit()].setVal(1e0);
	  lambdaIrreg[dit()].setVal(0e0);
	}
      
      consToDiffusiveProps(aco,etacell,rcv,consTemp,a_domains[ilev]);
      //average eta at the cell interface
      EBLevelDataOps::averageCellToFace(eta, etacell,a_grids[ilev],a_ebisl[ilev],a_domains[ilev], 0, 0, 1, true);
      EBLevelDataOps::assign(lambda,eta);
      EBLevelDataOps::scale(lambda,m_lambdafac);

      Real dx = a_dxs[ilev];
      RealVect dxv = dx*RealVect::Unit;
      for (DataIterator dit = eta.dataIterator(); dit.ok(); ++dit)
	{
	  Box grownBox = grow(a_grids[ilev].get(dit()), 1);
	  grownBox &= a_domains[ilev];
	  const EBISBox& ebisBox = eta[dit()].getEBISBox();
	  IntVectSet ivsIrreg = ebisBox.getIrregIVS(grownBox);
	  for (VoFIterator vofit(ivsIrreg, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
	    {
	      RealVect bndryCentroid = ebisBox.bndryCentroid(vofit());
	      bndryCentroid *= dx;
	      VoFStencil extrapSten;
	      int order = EBArith::getExtrapolationStencil(extrapSten, bndryCentroid, dxv, vofit(), ebisBox);
	  
	      Real extrapval = applyVoFStencil(extrapSten, etacell[dit()], 0);
	      etaIrreg[dit()](vofit(), 0) = extrapval;
	      lambdaIrreg[dit()](vofit(), 0) = m_lambdafac*extrapval;
	    }
	}
      aco        .exchange(Interval(0,0));
      eta        .exchange(Interval(0,0));
      lambda     .exchange(Interval(0,0));
      etaIrreg   .exchange(Interval(0,0));
      lambdaIrreg.exchange(Interval(0,0));
      rcv        .exchange(Interval(0,0));

    }
}
void
EBAMRNavier::
addDKEDt(LevelData<EBCellFAB>& a_srcT,const LevelData<EBCellFAB>& a_updateV, const LevelData<EBCellFAB>& a_initialV, const LevelData<EBCellFAB>& a_consState)
{
  for (DataIterator dit = a_srcT.dataIterator(); dit.ok(); ++dit)
    {
      const EBCellFAB& ebV1 = a_updateV[dit()];
      const EBCellFAB& ebV0 = a_initialV[dit()];
      const EBCellFAB& ebC  = a_consState[dit()];
      EBCellFAB& ebS = a_srcT[dit()];

      Box intBox = ebS.getRegion() & ebV0.getRegion();


      const BaseFab<Real>& regV1 = ebV1.getSingleValuedFAB();
      const BaseFab<Real>& regV0 = ebV0.getSingleValuedFAB();
      BaseFab<Real>&       regS  = ebS.getSingleValuedFAB();
      const BaseFab<Real>& regC  = ebC.getSingleValuedFAB();


      FORT_INCRDEKDT(CHF_FRA1(regS,0),CHF_CONST_FRA(regV1),CHF_CONST_FRA(regV0),CHF_CONST_FRA(regC),CHF_CONST_REAL(m_dt),CHF_BOX(intBox)); 

      const EBISBox& ebisBox = ebC.getEBISBox();
      IntVectSet ivsMulti = ebisBox.getMultiCells(intBox);
      if(!ivsMulti.isEmpty())
	{
	  int densIndx = m_ebPolytropic->densityIndexC();
	  Interval velIntrv = Interval(0, SpaceDim-1);
    
	  for (VoFIterator vofit(ivsMulti, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
	    {
	      const VolIndex& vof = vofit();

	      Real dens = ebC(vof,densIndx);
	      Real kinetic1 = 0e0;Real kinetic0 = 0e0;
	      for (int ivar = 0; ivar <= velIntrv.end(); ivar++)
		{
		  Real vel = ebV0(vof,ivar);
		  kinetic0 += vel*vel;
		  vel = ebV1(vof,ivar);
		  kinetic1 += vel*vel;
		}
	      ebS(vof,0) = ebS(vof,0) -dens*(kinetic1-kinetic0)/2e0/m_dt;
	    }
	}
      
    }
}
/*****************/
/****************************/
void EBAMRNavier::getVorticity(LevelData<EBCellFAB>& a_vort)
{

          // If there is a coarser level interpolate undefined ghost cells
  int densityIndex = m_ebPatchGodunov->densityIndex();
  Interval momIntrv = m_ebPolytropic->momentumInterval();
  Interval intervDensMom(Min(densityIndex,momIntrv.begin()), Max(densityIndex,momIntrv.end()));
  EBCellFactory factory(m_ebisl);
  int nghost = m_nGhost;
  int nCons = m_ebPatchGodunov->numConserved();
  LevelData<EBCellFAB> consTemp(m_grids, nCons, (IntVect::Unit)*nghost, factory);
  Interval consInterv(0, nCons-1);
  m_stateNew.copyTo(consInterv, consTemp, consInterv);
  if (m_hasCoarser)
    {
      const EBAMRNavier* amrGodCoarserPtr = getCoarserLevel();
      int refRatCrse = amrGodCoarserPtr->refRatio();
      EBPWLFillPatch patcher(m_grids,
			     amrGodCoarserPtr->m_grids,
			     m_ebisl,
			     amrGodCoarserPtr->m_ebisl,
			     amrGodCoarserPtr->m_domainBox,
			     refRatCrse, m_nComp, nghost);
      
      Real coarTimeOld = 0.0;
      Real coarTimeNew = 1.0;
      Real fineTime    = 0.0;
      patcher.interpolate(consTemp,
			  amrGodCoarserPtr->m_stateOld,
			  amrGodCoarserPtr->m_stateNew,
			  coarTimeOld,
			  coarTimeNew,
			  fineTime,
			  intervDensMom);
    }
  consTemp.exchange(intervDensMom);
  consToVelTemp(m_initialV,m_initialT,consTemp);
  computeVorticity(a_vort, m_initialV);
}
void
EBAMRNavier::
computeVorticity(LevelData<EBCellFAB>& a_vort, const LevelData<EBCellFAB>& a_velo)

{
  CH_TIME("EBAMRViscous::computeVorticity");

  const DisjointBoxLayout& grid = m_grids;
  const EBISLayout&  ebis = m_ebisl;
  const ProblemDomain& domain = m_problem_domain;
  const Real& dxLev = m_dx[0];
  EBCellFactory ebcellfact(ebis);
  //define the data holder
  //vorticity is a vector in 3d, scalar in 2d
  int ncomp = 1;
  if (SpaceDim==3)
    {
      ncomp = SpaceDim;
    }


  //do actual computation
  for (DataIterator dit = grid.dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB&  vortFAB = a_vort[dit()];
      const EBCellFAB&  veloFAB = a_velo[dit()];
      const EBGraph& ebgraph =   ebis[dit()].getEBGraph();
      vortFAB.setVal(0.);

      BaseFab<Real>&  regVort = vortFAB.getSingleValuedFAB();
      const BaseFab<Real>&  regVelo = veloFAB.getSingleValuedFAB();
      Box interiorBox = grid.get(dit());
      interiorBox.grow(1);
      interiorBox &= domain;
      interiorBox.grow(-1);
      int vortIndex;
#if CH_SPACEDIM==3
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          vortIndex = idir;
#else
          vortIndex = 0;
          int idir = 3;
#endif

          //compute on all cells as regular
          FORT_COMPUTEVORT(CHF_FRA1(regVort, vortIndex),
                           CHF_CONST_FRA(regVelo),
                           CHF_BOX(interiorBox),
                           CHF_CONST_REAL(dxLev),
                           CHF_CONST_INT(idir));;

          //do cells on domain boundary as if they were irregular
          IntVectSet ivsIrreg(grid.get(dit()));
          ivsIrreg -= interiorBox;
          ivsIrreg |= ebgraph.getIrregCells(interiorBox);
          int diffDirVec[2];
          if (SpaceDim==2 || idir ==2)
            {
              diffDirVec[0] = 0;
              diffDirVec[1] = 1;
            }
          else if (idir == 0)
            {
              diffDirVec[0] = 1;
              diffDirVec[1] = 2;
            }
          else if (idir==1)
            {
              diffDirVec[0] = 2;
              diffDirVec[1] = 0;
            }
          else
            {
              MayDay::Error("missed a case");
            }

          for (VoFIterator vofit(ivsIrreg, ebgraph); vofit.ok(); ++vofit)
            {
              const VolIndex vof = vofit();
              Real vortValue = 0.0;

              //taking the derivative d(udvdir)/d(xdiffdir)
              for (int idiff = 0; idiff < 2; idiff++)
                {

                  Real signDiff = 0;
                  int diffDir= -1;
                  int velComp = -1;
                  if (idiff == 0)
                    {
                      signDiff = 1.0;
                      diffDir  = diffDirVec[0];
                      velComp  = diffDirVec[1];
                    }
                  else if (idiff == 1)
                    {
                      signDiff = -1.0;
                      diffDir  = diffDirVec[1];
                      velComp  = diffDirVec[0];
                    }
                  else
                    {
                      MayDay::Error("missed a case");
                    }

                  Vector<FaceIndex> hiFaces =  ebgraph.getFaces(vof, diffDir, Side::Hi);
                  Vector<FaceIndex> loFaces =  ebgraph.getFaces(vof, diffDir, Side::Lo);

                  bool hasHi = (hiFaces.size() == 1) && (!hiFaces[0].isBoundary());
                  bool hasLo = (loFaces.size() == 1) && (!loFaces[0].isBoundary());
                  Real diffValue = 0.0;
                  if (hasHi && hasLo)
                    {
                      const VolIndex& vofHi = hiFaces[0].getVoF(Side::Hi);
                      const VolIndex& vofLo = loFaces[0].getVoF(Side::Lo);
                      Real hiValue = veloFAB(vofHi, velComp);
                      Real loValue = veloFAB(vofLo, velComp);
                      diffValue =0.5*(hiValue - loValue)/dxLev;
                    }
                  else if (hasHi)
                    {
                      const VolIndex& vofHi = hiFaces[0].getVoF(Side::Hi);
                      Real hiValue = veloFAB(vofHi, velComp);
                      Real loValue = veloFAB(vof,   velComp);
                      diffValue =(hiValue - loValue)/dxLev;
                    }
                  else if (hasLo)
                    {
                      const VolIndex& vofLo = loFaces[0].getVoF(Side::Lo);
                      Real hiValue = veloFAB(vof  , velComp);
                      Real loValue = veloFAB(vofLo, velComp);
                      diffValue =  (hiValue - loValue)/dxLev;
                    }
                  else
                    {
                      diffValue = 0.0;
                    }
                  vortValue += signDiff*diffValue;
                } //end loop over idiff

              vortFAB(vof, vortIndex) = vortValue;

            } //end loop over irregular vofs
#if CH_SPACEDIM==3
        } //end loop over vort components in 3d
#endif

    }
}
void
EBAMRNavier::
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
/********************************/
void 
EBAMRNavier::copyLevel(LevelData<EBCellFAB>& a_dest, const LevelData<EBCellFAB>& a_src, const Interval& destComps, const Interval& srcComps)
{
  for (DataIterator dit = a_dest.dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB& dest     = a_dest[dit()];
      const EBCellFAB& src      = a_src[dit()];

      Box intBox = dest.getRegion() & src.getRegion();
      dest.copy(intBox, destComps,   intBox,  src,   srcComps);
    }
}      
/*********/
void EBAMRNavier::residual(LevelData<EBCellFAB>&     a_updated)
{
  EBPatchGodunov::s_whichLev = m_level;
  
  EBLevelDataOps::setToZero(a_updated);
  copyLevel(a_updated,m_stateNew,a_updated.interval(),m_stateNew.interval());


  //set up arguments to step
  //undefined lfr in case we need it
  EBFluxRegister lfr;
  //undefined leveldata in case we need it
  const LevelData<EBCellFAB> ld;

  //set arguments to dummy arguments and
  //then fix if they are available
  EBFluxRegister* coarFR = &lfr;
  EBFluxRegister* fineFR = &lfr;
  const LevelData<EBCellFAB>* coarDataOld = &ld;
  const LevelData<EBCellFAB>* coarDataNew = &ld;

  Real tCoarserOld = 0.0;
  Real tCoarserNew = 0.0;
  if (m_hasCoarser)
    {
      EBAMRNavier* coarPtr = getCoarserLevel();
      //recall that my flux register goes between my
      //level and the next finer level
      coarFR = &coarPtr->m_ebFluxRegister;
      coarDataOld = &coarPtr->m_stateOld;
      coarDataNew = &coarPtr->m_stateNew;
      tCoarserNew = coarPtr->m_time;
      tCoarserOld = tCoarserNew - coarPtr->m_dt;
      //time should never be greater than the newest coarse
      //time.  time might be very slightly smaller than
      //tCoarserOld because of the above subtraction.
      Real eps = 1.0e-10;
      //correct for said floating-point nastiness
      m_time = Min(Max(m_time, tCoarserOld),tCoarserNew-eps);
    }
  if (m_hasFiner)
    {
      //recall that my flux register goes between my
      //level and the next finer level
      fineFR = &m_ebFluxRegister;
    }

  //source term
  LevelData<EBCellFAB> diffusiveSrc;
  if(m_hasDiffusion && m_explicitDiffusion) ////lucacancel testresidual
    {
      EBCellFactory factory(m_ebisl);
      diffusiveSrc.define(m_grids, m_nComp, m_stateNew.ghostVect(),factory);
      makeDiffusiveSource(diffusiveSrc);
      m_ebLevelGodunov.setESource(&diffusiveSrc);
      ///copyLevel(a_updated,diffusiveSrc,a_updated.interval(), diffusiveSrc.interval()); //lucacancel testresidual
      //return; //lucacancel testresidual
    }
    
    EBPatchGodunov::setCurLevel(m_level);
    m_ebLevelGodunov.residual(a_updated,
			      m_stateNew,
			      m_massDiff,
			      *fineFR,
			      *coarFR,
			      *coarDataOld,
			      *coarDataNew,
			      m_time,
			      tCoarserOld,
			      tCoarserNew,
			      m_dt);
    
  
}
/*********/
void
EBAMRNavier::
PlotReactiveSource()
{


  Vector<EBAMRNavier*>        hierarchy;
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

  //if (numLevels == 1 || (m_level == numLevels-1 && (abs(m_coarser_level_ptr->time() - m_time) <= m_dt)))
  {
    m_plot ++;
        
    int numC = m_nComp;
    IntVect ivGhost = m_nGhost*IntVect::Unit;
    Vector<LevelData<EBCellFAB>* > srcdata(numLevels,NULL);
    Vector<LevelData<EBCellFAB>* > resdata(numLevels,NULL);
    Vector<LevelData<EBCellFAB>* > zesdata(numLevels,NULL);
    for (int ilev = 0; ilev < numLevels; ++ilev)
      {
	EBCellFactory factory(ebisl[ilev]);
	srcdata[ilev] = new LevelData<EBCellFAB>();
	srcdata[ilev]->define(grids[ilev], m_nComp, ivGhost, factory);
	EBLevelDataOps::setToZero(*srcdata[ilev]);
	hierarchy[ilev]->makeDiffusiveSource(*srcdata[ilev]);
	resdata[ilev] = new LevelData<EBCellFAB>();
	resdata[ilev]->define(grids[ilev], m_nComp, ivGhost, factory);
	EBLevelDataOps::setToZero(*resdata[ilev]);
	hierarchy[ilev]->residual(*resdata[ilev]);


	/*
	for (DataIterator dit = grids[ilev].dataIterator(); dit.ok(); ++dit)
	{
	  EBCellFAB& Ndata = (*srcdata[ilev])[dit()];
	  Box box = Ndata.getRegion(); 
	  const EBISBox& ebisBox = ebisl[ilev][dit()];
	  IntVectSet ivs(box);
	  for (VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
	    {
	      const VolIndex& vof = vofit();
	      if(Ndata(vof,2) > 0.1) pout() << vof << ", " << ilev << ", TAGG, " << Ndata(vof,2) << endl;
	    }
	
	}
	*/

		
      }
    //the zes stuff should be cmntd-out and used only to chekc
    // the continuity residual by comparing Resid to the case with
    // one component zeroed out
    int zComp = 2;
    for (int ilev = 0; ilev < numLevels; ++ilev)
      {
	EBLevelDataOps::setVal(hierarchy[ilev]->m_stateNew,0e0,zComp);
	EBLevelDataOps::setVal(hierarchy[ilev]->m_stateOld,0e0,zComp);
      }
    for (int ilev = 0; ilev < numLevels; ++ilev)
      {
	EBCellFactory factory(ebisl[ilev]);
	zesdata[ilev] = new LevelData<EBCellFAB>();
	zesdata[ilev]->define(grids[ilev], m_nComp, ivGhost, factory);
	EBLevelDataOps::setToZero(*zesdata[ilev]);
	hierarchy[ilev]->residual(*zesdata[ilev]);
	EBAMRDataOps::checkNANINF(srcdata);	
      }

 
    outputHDF(srcdata,string("Sterms"));
    outputHDF(resdata,string("Resid"));
    outputHDF(zesdata,string("Zesid"));
    for (int ilev = 0; ilev < numLevels; ++ilev)
      {
	delete srcdata[ilev];
	delete resdata[ilev];
	delete zesdata[ilev];
      }
	
    pout() << "Exiting Program from PlotReactiveSource" << endl;
    exit(0);
  }
}

void
EBAMRNavier::
PlotVorticity(const int& a_nplot)
{


  Vector<EBAMRNavier*>        hierarchy;
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

  //if (numLevels == 1 || (m_level == numLevels-1 && (abs(m_coarser_level_ptr->time() - m_time) <= m_dt)))
  {
    m_plot = a_nplot;
        
    int numC = m_nComp;
    IntVect ivGhost = m_nGhost*IntVect::Unit;
    Vector<LevelData<EBCellFAB>* > vortdata(numLevels,NULL);
    for (int ilev = 0; ilev < numLevels; ++ilev)
      {
	EBCellFactory factory(ebisl[ilev]);
	vortdata[ilev] = new LevelData<EBCellFAB>(grids[ilev], SpaceDim, ivGhost, factory);
	EBLevelDataOps::setToZero(*vortdata[ilev]);
	hierarchy[ilev]->getVorticity(*vortdata[ilev]);
		
      }

 
    outputHDF(vortdata,string("Vort"));
    for (int ilev = 0; ilev < numLevels; ++ilev)
      {
	delete vortdata[ilev];
      }
  }
}
void
EBAMRNavier::
PlotVelocity(const int& a_nplot, const Real& a_Vel0, const bool& a_removeSelfSimilar)
{


  Vector<EBAMRNavier*>        hierarchy;
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

  //if (numLevels == 1 || (m_level == numLevels-1 && (abs(m_coarser_level_ptr->time() - m_time) <= m_dt)))
  {
    m_plot = a_nplot;
        
    int numC = m_nComp;
    IntVect ivGhost = m_nGhost*IntVect::Unit;
    Vector<LevelData<EBCellFAB>* > velodata(numLevels,NULL), tempdata(numLevels,NULL);
    for (int ilev = 0; ilev < numLevels; ++ilev)
      {
	EBCellFactory factory(ebisl[ilev]);
	velodata[ilev] = new LevelData<EBCellFAB>(grids[ilev], SpaceDim, ivGhost, factory);
	tempdata[ilev] = new LevelData<EBCellFAB>(grids[ilev], 1, ivGhost, factory);
	EBLevelDataOps::setToZero(*velodata[ilev]);
	EBLevelDataOps::setToZero(*tempdata[ilev]);
	hierarchy[ilev]->getVelTemp(*velodata[ilev],*tempdata[ilev],a_removeSelfSimilar);		
      }

 
    EBAMRDataOps::scale(velodata,a_Vel0);
    outputHDF(velodata,string("Vel"));
    for (int ilev = 0; ilev < numLevels; ++ilev)
      {
	delete velodata[ilev];
	delete tempdata[ilev];
      }
  }
}

void EBAMRNavier::outputHDF(const Vector<LevelData<EBCellFAB>* >& a_Ptr, const string& a_name)
{
  string  filename = a_name;
  std::ostringstream oss;
  oss << m_plot;
  filename += oss.str();
  filename.append(".hdf5");
  Real dumReal =  1.0;

  Vector<EBAMRNavier*>        hierarchy;
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
  writeEBHDF5(filename, grids, a_Ptr, namesrhs,
	      lev0Dom.domainBox(), lev0Dx, dumReal, m_time,
	      refRat, numLevels,
	      replaceCovered, coveredValues);
}
/********************************/
void EBAMRNavier::getIVBoxLimits(IntVect& a_loInt, IntVect& a_hiInt, 
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
void
EBAMRNavier::dumpDebug()
{
  dumpDebug("arg");
}

void
EBAMRNavier::dumpDebug(const string& a_debstring)
{
  int ilev = m_level;
}
#include "NamespaceFooter.H"
