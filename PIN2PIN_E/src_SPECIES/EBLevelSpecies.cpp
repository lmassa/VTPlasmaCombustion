#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "DebugOut.H"
#include "EBLevelSpecies.H"
#include "BaseIVFactory.H"
#include "BaseIFFactory.H"
#include "BaseIFFAB.H"
#include "EBFluxFAB.H"
#include "FaceIterator.H"
#include "REAL.H"
#include "EBCellFactory.H"
#include "FaceIterator.H"
#include <cstdio>
#include "EBAMRIO.H"
#include "EBArith.H"
#include "NamespaceHeader.H"
int EBLevelSpecies::s_timestep = 0;
/*****************************/
/*****************************/
EBLevelSpecies::EBLevelSpecies()
{
  m_isDefined = false;
  m_ebPatchSpecies = NULL;
  m_SFD = false;
}
/*****************************/
/*****************************/
void
EBLevelSpecies::define(const DisjointBoxLayout&      a_thisDBL,
                       const DisjointBoxLayout&      a_coarDBL,
                       const EBISLayout&             a_thisEBISL,
                       const EBISLayout&             a_coarEBISL,
                       const ProblemDomain&          a_domain,
                       const int&                    a_nRefine,
                       const RealVect&               a_dx,
                       const bool&                   a_useMassRedist,
                       const bool&                   a_doSmushing,
                       const bool&                   a_doRZCoords,
                       const bool&                   a_hasSourceTerm,
                       const EBPatchSpeciesFactory*  const a_patchSpecies,
                       const bool&                   a_hasCoarser,
                       const bool&                   a_hasFiner)
{
  CH_TIME("EBLevelSpecies::define");
  CH_assert(a_dx[0] > 0.0);
  CH_assert(a_nRefine > 0);

  m_isDefined = true;
  m_thisGrids = a_thisDBL;
  m_thisEBISL = a_thisEBISL;
  m_refRatCrse= a_nRefine;
  m_dx = a_dx;
  m_domain = a_domain;
  m_hasCoarser= a_hasCoarser;
  m_hasFiner= a_hasFiner;
  m_doSmushing = a_doSmushing;
  m_doRZCoords = a_doRZCoords;
  m_hasSourceTerm = a_hasSourceTerm;
  m_useMassRedist = false;//lucanote hardwired
  m_useFlattening = false; //lucanote hardwired
  if (m_doRZCoords && !m_hasSourceTerm)
    {
      MayDay::Error("LG: RZ implies need of a source term--inconsistent inputs");
    }
  if (m_hasCoarser)
    {
      m_coarGrids = a_coarDBL;
      m_coarEBISL = a_coarEBISL;
    }

  //lucacancel // delete seem to be able to delete this pointer becasue of the Plasma Physics stuff //not threadsafe
  //if (m_ebPatchSpecies != NULL) delete m_ebPatchSpecies;

  m_ebPatchSpecies = a_patchSpecies->create();
  m_ebPatchSpecies->define(m_domain, m_dx);

  // Determing the number of ghost cells necessary here
  m_nGhost = 4;
  m_nCons = m_ebPatchSpecies->numConserved();
  m_nFlux = m_ebPatchSpecies->numFluxes();

  //define redistribution object for this level
  //for now set to volume weighting
  m_ebLevelRedist.define(m_thisGrids, m_thisEBISL, m_domain, m_nCons);

  if (m_hasCoarser)
    {
      CH_TIME("coarse_stuff_defs");
      ProblemDomain domainCrse = coarsen(m_domain, m_refRatCrse);

      //patcher is defined with the number of conserved vars.
      m_patcher.define(m_thisGrids, m_coarGrids,
                       m_thisEBISL, m_coarEBISL,
                       domainCrse, m_refRatCrse, m_nCons,
                       m_nGhost);
    }

  m_irregSetsSmall.define(m_thisGrids);
  for (DataIterator dit = m_thisGrids.dataIterator();
      dit.ok(); ++dit)
    {
      const EBISBox& ebisBox = m_thisEBISL[dit()];
      if (!ebisBox.isAllCovered())
        {
          const Box&     thisBox = m_thisGrids.get(dit());
          m_irregSetsSmall[dit()] = ebisBox.getIrregIVS(thisBox);

        }
    }
  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      CH_TIME("flux_interpolant_defs");
      EBArith::defineFluxInterpolant(m_fluxInterpolants[faceDir],
                                     m_irregSetsGrown  [faceDir],
                                     m_thisGrids, m_thisEBISL,
                                     m_domain, m_nFlux, faceDir);
    }
  {
    CH_TIME("EBIrregFlux_defs");
    BaseIVFactory<Real> cellFactorySmall(m_thisEBISL, m_irregSetsSmall);
    m_nonConsDivergence.define(m_thisGrids, m_nCons,
                               IntVect::Zero, cellFactorySmall);
    m_ebIrregFaceFlux.define(m_thisGrids, m_nCons,
                             IntVect::Zero, cellFactorySmall);
  }

  {
    CH_TIME("coarse-fine_ivs_defs");
    //ebpatch needs coarse-fine IVS to know where to drop order for interpolation
    m_cfIVS.define(m_thisGrids);
    for (DataIterator ditOuter = m_thisGrids.dataIterator();
        ditOuter.ok(); ++ditOuter)
      {
        const EBISBox& ebisBox = m_thisEBISL[ditOuter()];
        if (!ebisBox.isAllCovered())
          {
            // make a grown box of the grid and then subtract off each grid in the
            //domain to make the coarse-fine IVS
            // CFIVS = grow(b, 1) - sum(bi)
            Box grownBox = grow(m_thisGrids.get(ditOuter()), 1);
            grownBox &= m_domain;
            IntVectSet complementIVS(grownBox);

            for (LayoutIterator litInner = m_thisGrids.layoutIterator();
                litInner.ok(); ++litInner)
              {
                Box innerBox = m_thisGrids.get(litInner());
                complementIVS -= innerBox;
              }

            m_cfIVS[ditOuter()] = complementIVS;
          }
      }
  }

  {
    CH_TIME("flattening_defs");
    //create temp data with the correct number of ghost cells
    IntVect ivGhost = m_nGhost*IntVect::Unit;
    EBCellFactory factory(m_thisEBISL);
    IntVect flatGhostIV = 3*IntVect::Unit;
    m_flattening.define(m_thisGrids, 1, flatGhostIV, factory);
  }
}

/*****************************/
void EBLevelSpecies::setVelTempPtrs(LevelData<EBFluxFAB>* a_faceVelocity,
				   LevelData<EBCellFAB>* a_cellVelocity,
				   LevelData<EBCellFAB>* a_cellTemperature)
{
  m_faceVelocity=a_faceVelocity;
  m_cellVelocity=a_cellVelocity;
  m_cellTemperature=a_cellTemperature;
}
/*****************************/
void
EBLevelSpecies::
floorConserved(LevelData<EBCellFAB>&         a_consState,
               Real a_time, Real a_dt)
{
  for (DataIterator dit = m_thisGrids.dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB& consState = a_consState[dit()];
      const IntVectSet& cfivs = m_cfIVS[dit()];
      const EBISBox& ebisBox = m_thisEBISL[dit()];
      const Box& cellBox = m_thisGrids.get(dit());
      m_ebPatchSpecies->setValidBox(cellBox, ebisBox, cfivs, a_time, a_dt);
      m_ebPatchSpecies->floorConserved(consState, cellBox);
    }
}
/*****************************/
void
EBLevelSpecies::
fillConsState(LevelData<EBCellFAB>&         a_consState,
              const LevelData<EBCellFAB>&   a_consStateCoarseOld,
              const LevelData<EBCellFAB>&   a_consStateCoarseNew,
              const Real&                   a_time,
              const Real&                   a_coarTimeOld,
              const Real&                   a_coarTimeNew)
{
  Interval consInterv(0, m_nCons-1);
  Interval fluxInterv(0, m_nFlux-1);

  if (m_hasCoarser)
    {
      m_patcher.interpolate(a_consState,
                            a_consStateCoarseOld,
                            a_consStateCoarseNew,
                            a_coarTimeOld,
                            a_coarTimeNew,
                            Max(a_time,a_coarTimeOld),
                            consInterv);
    }

  // Exchange all the data between grids
  a_consState.exchange(consInterv);
}
/***************/
void
EBLevelSpecies::
computeFlattening(Real a_time, Real a_dt,
                  LevelData<EBCellFAB>&         a_consState)
{
  CH_TIME("eblevelgodunov::compute_flattening");
  //compute flattening coefficients. this saves a ghost cell. woo hoo.
  int ibox = 0;
  bool verbose = false;
  for (DataIterator dit = m_thisGrids.dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB& consState = a_consState[dit()];
      const EBISBox& ebisBox = m_thisEBISL[dit()];
      if (!ebisBox.isAllCovered())
        {
          const Box& cellBox = m_thisGrids.get(dit());
          const IntVectSet& cfivs = m_cfIVS[dit()];
          m_ebPatchSpecies->setValidBox(cellBox, ebisBox, cfivs, a_time, a_dt);
          m_ebPatchSpecies->setCoveredConsVals(consState);
          int nPrim  = m_ebPatchSpecies->numPrimitives();
          EBCellFAB primState(ebisBox, consState.getRegion(), nPrim);
          int logflag = 0;
          //debug
          //verbose = true;
          //end debug
          m_ebPatchSpecies->consToPrim(primState, consState, consState.getRegion(), logflag, verbose);
          EBCellFAB& flatteningFAB = m_flattening[dit()];
          //this will set the stuff over the coarse-fine interface
          flatteningFAB.setVal(1.);
          if (m_ebPatchSpecies->usesFlattening())
            {
              m_ebPatchSpecies->computeFlattening(flatteningFAB,
                                                  primState,
                                                  cellBox);
            }
          ibox++;
        }
    }
  Interval zerointerv(0,0);
  m_flattening.exchange(zerointerv);
}
/***************/
Real
EBLevelSpecies::
doIrregularUpdate(EBFluxRegister&               a_fineFluxRegister,
                  EBFluxRegister&               a_coarFluxRegister,
                  LevelData<BaseIVFAB<Real> >&  a_massDiff,
                  Real a_time, Real a_dt,
                  LevelData<EBCellFAB>&         a_consState)
{
  //now do the irregular update and the max wave speed
  Real maxWaveSpeed = 1.0e-12;
  int ibox = 0;
  Interval consInterv(0, m_nCons-1);
  Interval fluxInterv(0, m_nFlux-1);
  for (DataIterator dit = m_thisGrids.dataIterator(); dit.ok(); ++dit, ibox++)
    {
      const Box& cellBox = m_thisGrids.get(dit());
      const EBISBox& ebisBox = m_thisEBISL[dit()];
      if (!ebisBox.isAllCovered())
        {
          const IntVectSet& cfivs = m_cfIVS[dit()];

          EBCellFAB& consState = a_consState[dit()];
          BaseIVFAB<Real>& redMass = a_massDiff[dit()];
	  
	  EBFluxFAB& advVel = (*m_faceVelocity)[dit];
	  EBCellFAB& cellVel= (*m_cellVelocity)[dit];
	  EBCellFAB& cellTemp=(*m_cellTemperature)[dit];
	  m_ebPatchSpecies->setVelocities(cellVel, advVel);
          m_ebPatchSpecies->setValidBox(cellBox, ebisBox, cfivs, a_time, a_dt);

          BaseIFFAB<Real> centroidFlux[SpaceDim];
          const BaseIFFAB<Real>* interpolantGrid[SpaceDim];
          const IntVectSet& ivsIrregSmall = m_irregSetsSmall[dit()];
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              const BaseIFFAB<Real>& interpol = m_fluxInterpolants[idir][dit()];
              interpolantGrid[idir] = &interpol;
              BaseIFFAB<Real>& fluxDir= centroidFlux[idir];
              fluxDir.define(ivsIrregSmall, ebisBox.getEBGraph(), idir, m_nFlux);
            }
          m_ebPatchSpecies->interpolateFluxToCentroids(centroidFlux, interpolantGrid, ivsIrregSmall);

          Real maxWaveSpeedGrid = 0.0;
          //update the state and interpolate the flux
          const BaseIVFAB<Real>& nonConsDiv = m_nonConsDivergence[dit()];
          const BaseIVFAB<Real>& ebIrregFlux = m_ebIrregFaceFlux[dit()];
          m_ebPatchSpecies->irregularUpdate(consState,
                                            maxWaveSpeedGrid, redMass,
                                            centroidFlux, ebIrregFlux, nonConsDiv,
                                            cellBox, ivsIrregSmall);

          m_ebLevelRedist.increment(redMass, dit(), consInterv);
          maxWaveSpeed = Max(maxWaveSpeed, maxWaveSpeedGrid);

          //do fluxregister mambo
          /*
            Coarse flux register is flux register with the coarse level.
            Fine flux register is the flux register with the fine level.
            To the finer level FR, this level is the coarse level.
            To the coarser level FR, this level is the fine level.
          */
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              Real scale = a_dt;

              BaseIFFAB<Real> fluxRegFlux;
              if (m_doRZCoords)
                {
                  //gather fluxes into flux-register compatible form
                  fluxRegFlux.define(ivsIrregSmall, ebisBox.getEBGraph(), idir, m_nCons);
                  m_ebPatchSpecies->assembleFluxIrr(fluxRegFlux, centroidFlux[idir],
                                                    idir,  cellBox, ivsIrregSmall);
                }

              if (m_hasFiner)
                {
                  a_fineFluxRegister.incrementCoarseIrregular(centroidFlux[idir],
                                                              scale,dit(),
                                                              consInterv, idir);
                }

              if (m_hasCoarser)
                {
                  for (SideIterator sit; sit.ok(); ++sit)
                    {
                      a_coarFluxRegister.incrementFineIrregular(centroidFlux[idir],
                                                                scale, dit(),
                                                                consInterv, idir,sit());
                    }
                }
            }
        } //!ebisBox.isAllCovered()
    }// end of loop over grids.
  a_consState.exchange(consInterv);
  return maxWaveSpeed;
}
/*****************************/
Real
EBLevelSpecies::
step(LevelData<EBCellFAB>&         a_consState,
     LevelData<BaseIVFAB<Real> >&  a_massDiff,
     LevelData<EBCellFAB>&         a_source,
     EBFluxRegister&               a_fineFluxRegister,
     EBFluxRegister&               a_coarFluxRegister,
     const LevelData<EBCellFAB>&   a_consStateCoarseOld,
     const LevelData<EBCellFAB>&   a_consStateCoarseNew,
     const Real&                   a_time,
     const Real&                   a_coarTimeOld,
     const Real&                   a_coarTimeNew,
     const Real&                   a_dt)
{
  CH_TIME("levelgodunov_step");
  Interval consInterv(0, m_nCons-1);
  Interval fluxInterv(0, m_nFlux-1);
  CH_assert(isDefined());
  CH_assert(a_consState.disjointBoxLayout() == m_thisGrids);
  m_hasSourceTerm = a_source.isDefined();
  EBPatchGodunov::s_conservativeSource  = true;
  EBPatchGodunov::s_curComp  = 1;
  

  {
    CH_TIME("fillConsState");
    //cf interpolation and exchange
    fillConsState(a_consState, a_consStateCoarseOld, a_consStateCoarseNew, a_time, a_coarTimeOld, a_coarTimeNew);
  }

  {
    CH_TIME("early_coarse_fine");
    m_ebLevelRedist.setToZero();
    // clear flux registers with fine level
    //remember this level is the coarse level to the fine FR
    if (m_hasFiner)
      {
        a_fineFluxRegister.setToZero();
      }
  }

  if(m_useFlattening) //should be always false
  {
    CH_TIME("computeFlattening");
    //compute flattening coefficients. this saves a ghost cell. woo hoo.
    computeFlattening(a_time, a_dt, a_consState);
  }

  {
    CH_TIME("doRegularUpdate");
    //this includes copying flux into flux interpolant and updating
    //regular grids and incrementing flux registers.
    doRegularUpdate(a_fineFluxRegister, a_coarFluxRegister, a_time, a_dt, a_consState, a_source);
  }

  Real maxWaveSpeed;
  {
    CH_TIME("doIrregularUpdate");
    //this does irregular update and deals with flux registers.
    //also computes the mass increment
    maxWaveSpeed = doIrregularUpdate(a_fineFluxRegister, a_coarFluxRegister, a_massDiff, a_time, a_dt, a_consState);
  }

  if (m_doSmushing)
    {
      CH_TIME("smushing");
      if (m_useMassRedist)
        {
          //if use mass weighting, need to
          //fix weights of redistribution object
          int densityIndex = m_ebPatchSpecies->densityIndex();
          a_consState.exchange(consInterv);
          m_ebLevelRedist.resetWeights(a_consState, densityIndex);
        }

      //do fluxregister samba
      //redistribute mass at this level.  Redistribution to other
      //levels (including reredist) is handled in EBAMRSpecies::postTimeStep
      m_ebLevelRedist.redistribute(a_consState, consInterv);
    }

  {
    CH_TIME("floors");
    floorConserved(a_consState, a_time, a_dt);
  }

  //max wave speed already gets broadcast
  maxWaveSpeed = getMaxWaveSpeed(a_consState);

  // Find the minimum of dt's over this level
  Real dtNew = m_dx[0] / maxWaveSpeed;

  //return the maximum stable time step
  return dtNew;
}
/***************/
void
EBLevelSpecies::
doRegularUpdate(EBFluxRegister&               a_fineFluxRegister,
                EBFluxRegister&               a_coarFluxRegister,
                Real                          a_time, 
		Real                          a_dt,
                LevelData<EBCellFAB>&         a_consState, 
		const LevelData<EBCellFAB>&   a_source)
{
  bool verbose = false;
  int ibox = 0;
  Interval consInterv(0, m_nCons-1);
  Interval fluxInterv(0, m_nFlux-1);
  if (m_SFD)
    {
      floorConserved(*m_qbar, a_time, a_dt);
    }
  for (DataIterator dit = m_thisGrids.dataIterator(); dit.ok(); ++dit, ibox++)
    {
      const Box& cellBox = m_thisGrids.get(dit());
      const EBISBox& ebisBox = m_thisEBISL[dit()];
      if (!ebisBox.isAllCovered())
        {
          const IntVectSet& cfivs = m_cfIVS[dit()];

          EBCellFAB& consState = a_consState[dit()];

	  EBFluxFAB& advVel = (*m_faceVelocity)[dit];
	  EBCellFAB& cellVel= (*m_cellVelocity)[dit];
	  EBCellFAB& cellTemp=(*m_cellTemperature)[dit];
	  m_ebPatchSpecies->setVelocities(cellVel, advVel);
          m_ebPatchSpecies->setValidBox(cellBox, ebisBox, cfivs, a_time, a_dt);


          if (m_SFD)
            {
	      MayDay::Error("Luca: Not dealing with SFD yet");
            }

          EBFluxFAB flux(ebisBox, cellBox, m_nFlux);
          BaseIVFAB<Real>& nonConsDiv    = m_nonConsDivergence[dit()];
          BaseIVFAB<Real>& ebIrregFlux   = m_ebIrregFaceFlux[dit()];
          flux.setVal(7.89);
          ebIrregFlux.setVal(7.89);
          const IntVectSet& ivsIrreg     = m_irregSetsSmall[dit()];
          const EBCellFAB& flatteningFAB = m_flattening[dit()];

          m_ebPatchSpecies->regularUpdate(consState, flux, ebIrregFlux,
                                          nonConsDiv,flatteningFAB,
                                          a_source[dit()], cellBox, ivsIrreg,
                                          dit(),verbose);

          //do fluxregister cha-cha
          /*
            Coarse flux register is flux register with the coarse level.
            Fine flux register is the flux register with the fine level.
            To the finer level FR, this level is the coarse level.
            To the coarser level FR, this level is the fine level.
          */
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              Real scale = a_dt;

              EBFaceFAB fluxRegFlux;
              if (m_doRZCoords)
                {
                  //gather fluxes into flux-register compatible form
                  fluxRegFlux.define(ebisBox, cellBox, idir, m_nCons);
                  m_ebPatchSpecies->assembleFluxReg(fluxRegFlux, flux[idir],
                                                    idir, cellBox);
                }
              if (m_hasFiner)
                {
                  a_fineFluxRegister.incrementCoarseRegular(flux[idir], scale,dit(),
                                                            consInterv, idir);
                }

              if (m_hasCoarser)
                {
                  for (SideIterator sit; sit.ok(); ++sit)
                    {
                      a_coarFluxRegister.incrementFineRegular(flux[idir],scale, dit(),
                                                              consInterv, idir,sit());
                    }
                }
            }

          //copy fluxes into sparse interpolant
          for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
            {
              IntVectSet ivsIrregGrown = m_irregSetsGrown[faceDir][dit()];
              ivsIrregGrown &= cellBox;
              FaceStop::WhichFaces stopCrit = FaceStop::SurroundingWithBoundary;

              BaseIFFAB<Real>& interpol = m_fluxInterpolants[faceDir][dit()];
              interpol.setVal(7.7777e7);
              EBFaceFAB& fluxDir = flux[faceDir];
              for (FaceIterator faceit(ivsIrregGrown, ebisBox.getEBGraph(),
                                      faceDir, stopCrit);
                  faceit.ok(); ++faceit)
                {
                  for (int ivar = 0; ivar < m_nFlux; ivar++)
                    {
                      interpol(faceit(), ivar) = fluxDir(faceit(), ivar);
                    }
                }

            }
        }
    }

  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      m_fluxInterpolants[faceDir].exchange(fluxInterv);
    }
}
/*****************************/
Real EBLevelSpecies::getMaxWaveSpeed(const LevelData<EBCellFAB>& a_state,
				     LevelData<EBCellFAB>* a_cellVelocity)
{
  m_cellVelocity=a_cellVelocity;
  return getMaxWaveSpeed(a_state);
}
/*****************************/
Real EBLevelSpecies::getMaxWaveSpeed(const LevelData<EBCellFAB>& a_state)
{
  CH_assert(a_state.disjointBoxLayout() == m_thisGrids);
  Real speed = 0.0;
  for (DataIterator dit = m_thisGrids.dataIterator(); dit.ok(); ++dit)
    {
      const Box& validBox   = m_thisGrids.get(dit());
      const EBISBox& ebisBox= m_thisEBISL[dit()];
      if (!ebisBox.isAllCovered())
        {
          //place holder not used maxWaveSpeed calc
          Real time = 0.0;
          const IntVectSet& cfivs = m_cfIVS[dit()];
	  
	  EBCellFAB& cellVel= (*m_cellVelocity)[dit];
	  m_ebPatchSpecies->setVelocities(cellVel);
          m_ebPatchSpecies->setValidBox(validBox, ebisBox, cfivs, time, time);
          Real speedOverBox = m_ebPatchSpecies->getMaxWaveSpeed(a_state[dit()],
                                                                validBox);
          speed = Max(speed,speedOverBox);
        }
    }

  // gather speed
  Vector<Real> all_speeds;
  gather(all_speeds,speed,uniqueProc(SerialTask::compute));
  if (procID() == uniqueProc(SerialTask::compute))
    {
      speed = all_speeds[0];
      for (int i = 1; i < all_speeds.size (); ++i)
        {
          speed = Max(speed,all_speeds[i]);
        }
    }
  broadcast(speed,uniqueProc(SerialTask::compute));

  return speed;
}
/*****************************/
bool
EBLevelSpecies::isDefined() const
{
  return m_isDefined;
}
/*****************************/
/*****************************/
//lucaadd
void EBLevelSpecies::getQbar(LevelData<EBCellFAB>*  a_qbar,Real a_chi)
{
  m_SFD = true;
  m_chi = a_chi;
  m_qbar = a_qbar;
}
/*****************************/
/*****************************/
EBLevelSpecies::~EBLevelSpecies()
{
  if (m_ebPatchSpecies != NULL)
    delete m_ebPatchSpecies;
}
/*****************************/
/*****************************/
#include "NamespaceFooter.H"
