#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "LoadBalance.H"
#include "EBArith.H"

#include "EBSpeciesFluxOp.H"
#include "EBArith.H"

#include "CH_Timer.H"
#include "EBSpeciesFluxOpFactory.H"
#include "EBCoarseAverage.H"
#include "NamespaceHeader.H"


//-----------------------------------------------------------------------
void
coarsenSFO(LevelData<EBCellFAB>                 & a_acoefCoar,
             LevelData<EBFluxFAB>               & a_DkCoar,
             LevelData<EBFluxFAB>               & a_vecWCoar,
             LevelData<BaseIVFAB<Real> >        & a_DkCoarIrreg,
             LevelData<BaseIVFAB<Real> >        & a_vecWCoarIrreg,
             LevelData<BaseIVFAB<Real> >        & a_dataCoar,
             const EBLevelGrid                  & a_eblgFine,
             const EBLevelGrid                  & a_eblgCoar,
             const LevelData<EBCellFAB>         & a_acoefFine,
             const LevelData<EBFluxFAB>         & a_DkFine,
             const LevelData<EBFluxFAB>         & a_vecWFine,
             const LevelData<BaseIVFAB<Real> >  & a_DkFineIrreg,
             const LevelData<BaseIVFAB<Real> >  & a_vecWFineIrreg,
             const LevelData<BaseIVFAB<Real> >  & a_dataFine,
             const int                          & a_refToDepth,
             const int                          & a_nSpecies)
{
  CH_assert(a_refToDepth > 0);

  Interval interv(0, 0);
  Interval intSpec(0, a_nSpecies-1);
  if (a_refToDepth == 1)
    {
      a_acoefFine. copyTo(interv,  a_acoefCoar, interv);
      a_DkFine.   copyTo(interv,    a_DkCoar, interv);
      a_vecWFine.copyTo(interv, a_vecWCoar, interv);
      a_DkFineIrreg.  copyTo(interv, a_DkCoarIrreg,   interv);
      a_vecWFineIrreg.copyTo(interv, a_vecWCoarIrreg, interv);
      a_dataFine.copyTo(intSpec, a_dataCoar, intSpec);
    }
  else
    {
      EBCoarseAverage averageOp(a_eblgFine.getDBL(),    a_eblgCoar.getDBL(),
                                a_eblgFine.getEBISL(),  a_eblgCoar.getEBISL(),
                                a_eblgCoar.getDomain(), a_refToDepth, 1,
                                a_eblgCoar.getEBIS());
      EBCoarseAverage averageSpec(a_eblgFine.getDBL(),    a_eblgCoar.getDBL(),
				  a_eblgFine.getEBISL(),  a_eblgCoar.getEBISL(),
				  a_eblgCoar.getDomain(), a_refToDepth, a_nSpecies,
				  a_eblgCoar.getEBIS());

      //MayDay::Warning("might want to figure out what harmonic averaging is in this context");
      averageOp.average( a_acoefCoar     ,  a_acoefFine     , interv);
      averageOp.average(   a_DkCoar     ,    a_DkFine     , interv);
      averageOp.average(   a_DkCoarIrreg,    a_DkFineIrreg, interv);
      averageOp.average(a_vecWCoar     , a_vecWFine     , interv);
      averageOp.average(a_vecWCoarIrreg, a_vecWFineIrreg, interv);
      averageSpec.average(a_dataCoar, a_dataFine, intSpec);

    }
  a_acoefCoar.exchange(interv);
  a_DkCoar.exchange(interv);
  a_vecWCoar.exchange(interv);
  a_DkCoarIrreg.exchange(interv);
  a_vecWCoarIrreg.exchange(interv);
  a_dataCoar.exchange(intSpec);
}
//-----------------------------------------------------------------------
EBSpeciesFluxOpFactory::
EBSpeciesFluxOpFactory(const Vector<EBLevelGrid>&                                  a_eblgs,
                        const Vector<RefCountedPtr<EBQuadCFInterp> >&               a_quadCFI,
                        const Real&                                                 a_alpha,
                        const Real&                                                 a_beta,
                        const Vector<RefCountedPtr<LevelData<EBCellFAB> > >&        a_acoef,
                         const Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&        a_Dk,
                         const Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&        a_vecW,
                         const Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >& a_DkIrreg,
                         const Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >& a_vecWIrreg,
                         const Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >& a_data,
                        const Real&                                                 a_dxCoarse,
                        const Vector<int>&                                          a_refRatio,
                        const RefCountedPtr<BaseDomainBCFactory>&                   a_domainBCFactory,
                        const RefCountedPtr<BaseEBBCFactory>    &                   a_ebBCFactory,
                        const IntVect&                                              a_ghostCellsPhi,
                        const IntVect&                                              a_ghostCellsRhs,
                        const int &                                                 a_relaxType,
                        const int &                                                 a_nSpecies,
                        int a_numLevels)
{
  CH_assert(a_eblgs.size() <= a_refRatio.size());
  m_dataBased = false;
  m_relaxType = a_relaxType;
  m_quadCFI = a_quadCFI;
  m_ghostCellsRhs = a_ghostCellsRhs;
  m_ghostCellsPhi = a_ghostCellsPhi;
  m_acoef         = a_acoef;
  m_acoef0.resize(a_acoef.size());
  m_acoef1.resize(a_acoef.size());
  m_Dk           = a_Dk;
  m_vecW        = a_vecW;
  m_vecWIrreg   = a_vecWIrreg;
  m_DkIrreg      = a_DkIrreg;
  m_data   = a_data;
  m_alpha         = a_alpha;
  m_beta          = a_beta;
  m_nSpecies = a_nSpecies;
  if (a_numLevels > 0)
  {
    m_numLevels = a_numLevels;
  }
  else
  {
    m_numLevels = a_eblgs.size();
  }

  m_domainBCFactory = a_domainBCFactory;
  m_ebBCFactory     = a_ebBCFactory;

  m_eblgs           = a_eblgs;
  m_refRatio        = a_refRatio;
  m_dx.resize(m_numLevels);

  m_dx[0] = a_dxCoarse;
  for (int ilev = 1; ilev < m_numLevels; ilev++)
  {
    m_dx[ilev] = m_dx[ilev-1];
    m_dx[ilev] /= m_refRatio[ilev-1];
  }
  m_alpha = a_alpha;
  m_beta = a_beta;

  m_eblgsMG.resize(m_numLevels);
  m_acoefMG.resize(m_numLevels);
  m_DkMG.resize(m_numLevels);
  m_vecWMG.resize(m_numLevels);
  m_DkIrregMG.resize(m_numLevels);
  m_vecWIrregMG.resize(m_numLevels);
  m_dataMG.resize(m_numLevels);
  m_hasMGObjects.resize(m_numLevels);
  for (int ilev = 0; ilev < m_numLevels; ilev++)
  {
    if ((ilev==0) || (m_refRatio[ilev] > 2))
    {
      m_hasMGObjects[ilev] = true;

      int mgRef = 2;
      m_eblgsMG      [ilev].resize(0);
      m_acoefMG      [ilev].resize(0);
      m_DkMG        [ilev].resize(0);
      m_vecWMG     [ilev].resize(0);
      m_DkIrregMG   [ilev].resize(0);
      m_vecWIrregMG[ilev].resize(0);
      m_dataMG[ilev].resize(0);

      m_eblgsMG      [ilev].push_back(m_eblgs      [ilev]);
      m_acoefMG      [ilev].push_back(m_acoef      [ilev]);
      m_DkMG        [ilev].push_back(m_Dk        [ilev]);
      m_vecWMG     [ilev].push_back(m_vecW     [ilev]);
      m_DkIrregMG   [ilev].push_back(m_DkIrreg   [ilev]);
      m_vecWIrregMG[ilev].push_back(m_vecWIrreg[ilev]);
      m_dataMG[ilev].push_back(m_data[ilev]);

      bool hasCoarser = true;
      while (hasCoarser)
      {
        int imgsize = m_eblgsMG[ilev].size();
        const EBLevelGrid& eblgFine=  m_eblgsMG[ilev][imgsize-1];
        DisjointBoxLayout dblCoarMG;
        ProblemDomain  domainCoarMG;
        int maxBoxSize = 32;
        bool dumbool;
        hasCoarser = EBAMRPoissonOp::getCoarserLayouts(dblCoarMG,
                                                       domainCoarMG,
                                                       eblgFine.getDBL(),
                                                       eblgFine.getEBISL(),
                                                       eblgFine.getDomain(),
                                                       mgRef,
                                                       eblgFine.getEBIS(),
                                                       maxBoxSize, dumbool);

        if (hasCoarser)
        {
          m_eblgsMG[ilev].push_back(EBLevelGrid(dblCoarMG, domainCoarMG, 4, eblgFine.getEBIS()));
          int img = m_eblgsMG[ilev].size() - 1;
          const EBLevelGrid& eblgCoar = m_eblgsMG[ilev][img  ];
          const EBLevelGrid& eblgFine = m_eblgsMG[ilev][img-1];

          int nghost = 1;
          LayoutData<IntVectSet> irregSets(eblgCoar.getDBL());
          for (DataIterator dit = eblgCoar.getDBL().dataIterator(); dit.ok(); ++dit)
          {
            Box grownBox = grow(eblgCoar.getDBL().get(dit()), nghost);
            grownBox &= domainCoarMG;
            irregSets[dit()] = eblgCoar.getEBISL()[dit()].getIrregIVS(grownBox);
          }
          EBFluxFactory       ebfluxfact(eblgCoar.getEBISL());
          EBCellFactory       ebcellfact(eblgCoar.getEBISL());
          BaseIVFactory<Real> baseivfact(eblgCoar.getEBISL(), irregSets);

          RefCountedPtr<LevelData<BaseIVFAB<Real> > >    DkIrregCoar( new LevelData<BaseIVFAB<Real> >(eblgCoar.getDBL(), 1, nghost*IntVect::Unit, baseivfact) );
          RefCountedPtr<LevelData<BaseIVFAB<Real> > > vecWIrregCoar( new LevelData<BaseIVFAB<Real> >(eblgCoar.getDBL(), 1, nghost*IntVect::Unit, baseivfact) );
          RefCountedPtr<LevelData<EBCellFAB> >              acoefCoar( new LevelData<EBCellFAB>       (eblgCoar.getDBL(), 1, nghost*IntVect::Unit, ebcellfact) );
          RefCountedPtr<LevelData<EBFluxFAB> >                DkCoar( new LevelData<EBFluxFAB>       (eblgCoar.getDBL(), 1, nghost*IntVect::Unit, ebfluxfact) );
          RefCountedPtr<LevelData<EBFluxFAB> >             vecWCoar( new LevelData<EBFluxFAB>       (eblgCoar.getDBL(), 1, nghost*IntVect::Unit, ebfluxfact) );
          RefCountedPtr<LevelData<BaseIVFAB<Real> > > dataCoar( new LevelData<BaseIVFAB<Real> >(eblgCoar.getDBL(), m_nSpecies, nghost*IntVect::Unit, baseivfact) );

          coarsenSFO(*acoefCoar, *DkCoar, *vecWCoar, *DkIrregCoar, *vecWIrregCoar, *dataCoar, eblgFine, eblgCoar,
		     *m_acoefMG[ilev][img-1], *m_DkMG[ilev][img-1], *m_vecWMG[ilev][img-1],
		     *m_DkIrregMG[ilev][img-1], *m_vecWIrregMG[ilev][img-1], *m_dataMG[ilev][img-1], mgRef,m_nSpecies);

          m_acoefMG      [ilev].push_back(acoefCoar);
          m_DkMG        [ilev].push_back( DkCoar);
          m_vecWMG     [ilev].push_back(  vecWCoar);
          m_DkIrregMG   [ilev].push_back( DkIrregCoar);
          m_vecWIrregMG[ilev].push_back(  vecWIrregCoar);
          m_dataMG[ilev].push_back(       dataCoar);
        }
      }

    }
    else
    {
      m_hasMGObjects[ilev] = false;
    }
  }
}

/****/
EBSpeciesFluxOpFactory::~EBSpeciesFluxOpFactory()
{
}
/****/
int
EBSpeciesFluxOpFactory::
refToFiner(const ProblemDomain& a_domain) const
{
  int retval = -1;
  bool found = false;
  for (int ilev = 0; ilev < m_eblgs.size(); ilev++)
    {
      if (m_eblgs[ilev].getDomain() == a_domain)
        {
          retval = m_refRatio[ilev];
          found = true;
        }
    }
  if (!found)
    {
      MayDay::Error("Domain not found in AMR hierarchy");
    }
  return retval;
}


//-----------------------------------------------------------------------
EBSpeciesFluxOp*
EBSpeciesFluxOpFactory::
MGnewOp(const ProblemDomain& a_domainFine,
        int                  a_depth,
        bool                 a_homoOnly)
{
  //find out if there is a real starting point here.
  int ref=-1;
  bool found = false;

  RefCountedPtr<LevelData<EBCellFAB> >               acoef, acoef0, acoef1;
  RefCountedPtr<LevelData<EBFluxFAB> >                 Dk;
  RefCountedPtr<LevelData<EBFluxFAB> >              vecW;
  RefCountedPtr<LevelData<BaseIVFAB<Real> > >     DkIrreg;
  RefCountedPtr<LevelData<BaseIVFAB<Real> > >  vecWIrreg;
  RefCountedPtr<LevelData<BaseIVFAB<Real> > >  data;

  for (int ilev = 0; ilev < m_numLevels && !found; ilev++)
  {
    if (a_domainFine == m_eblgs[ilev].getDomain())
    {
      found = true;
      ref = ilev ;
    }
  }
  if (!found)
  {
    MayDay::Error("No corresponding AMRLevel to starting point of MGnewOp");
  }

  //multigrid operator.  coarn by two from depth  no
  //coarse or finer to worry about.
  Real          dxMGLevel;
  EBLevelGrid eblgMGLevel;
  EBLevelGrid eblgCoarMG;
  RefCountedPtr<EBQuadCFInterp> quadCFI; //only defined if on an amr level
  Real dxCoar = 1.0;
  dxCoar *= -1.0;
  int refToDepth = 1;
  if (ref > 0)
  {
    dxCoar = m_dx[ref-1];
  }
  bool hasCoarMGObjects = false;
  if (a_depth == 0)
  {
    eblgMGLevel    = m_eblgs[ref];

    acoef = m_acoef[ref];
    if (m_acoef0[ref] != NULL)
    {
      CH_assert(m_acoef1[ref] != NULL);
      acoef0 = m_acoef0[ref];
      acoef1 = m_acoef1[ref];
    }

    Dk        =  m_Dk[ref];
    vecW      =  m_vecW[ref];
    DkIrreg    = m_DkIrreg[ref];
    vecWIrreg =  m_vecWIrreg[ref];
    data      =  m_data[ref];
    dxMGLevel      = m_dx[ref];
    quadCFI        = m_quadCFI[ref];

    hasCoarMGObjects = m_hasMGObjects[ref];
    if (hasCoarMGObjects)
    {
      eblgCoarMG = m_eblgsMG[ref][1];
    }
  }
  else
  {
    int icoar = 1;
    for (int idep = 0; idep < a_depth; idep++)
    {
      icoar *= 2;
    }
    refToDepth = icoar;
    const ProblemDomain domainFine = m_eblgs[ref].getDomain();
    ProblemDomain domainBoxMGLevel = coarsen(domainFine, icoar);
    bool foundMGLevel = false;
    int numMGLevels = m_eblgsMG[ref].size();
    for (int img = 0; img < numMGLevels; img++)
    {
      if (m_eblgsMG[ref][img].getDomain() == domainBoxMGLevel)
      {
        eblgMGLevel = m_eblgsMG[ref][img];
        acoef = m_acoefMG[ref][img];
        Dk       =  m_DkMG[ref][img];
        vecW    =   m_vecWMG[ref][img];
        DkIrreg =   m_DkIrregMG[ref][img];
        vecWIrreg = m_vecWIrregMG[ref][img];
        data = m_dataMG[ref][img];
        foundMGLevel = true;

        hasCoarMGObjects = ((img+1) < (numMGLevels));
        if (hasCoarMGObjects)
        {
          eblgCoarMG = m_eblgsMG[ref][img+1];
        }
        break;
      }
    }
    bool coarsenable = foundMGLevel;

    dxMGLevel = m_dx[ref];
    dxMGLevel *= Real(icoar);

    if (!coarsenable)
    {
      //not coarsenable.
      //return null
      return NULL;
    }
  }

  
  m_ebBCFactory->setData(data, m_ispec, m_wallBCInfo);
  SpeciesFluxBaseEBBC*     speciesEBBC  = (SpeciesFluxBaseEBBC*)     m_ebBCFactory->create(eblgMGLevel.getDomain(), eblgMGLevel.getEBISL(),
      dxMGLevel*RealVect::Unit, &m_ghostCellsPhi, &m_ghostCellsRhs);
  if (m_dataBased)
    {
      speciesEBBC->setData(data, m_ispec, m_wallBCInfo);
    }
  SpeciesFluxBaseDomainBC* speciesDomBC = (SpeciesFluxBaseDomainBC*) m_domainBCFactory->create(eblgMGLevel.getDomain(), eblgMGLevel.getEBISL(),
      dxMGLevel*RealVect::Unit);

  RefCountedPtr<SpeciesFluxBaseEBBC>      ebbc(speciesEBBC);
  RefCountedPtr<SpeciesFluxBaseDomainBC> dombc(speciesDomBC);
  //no need for finer or coarser amr levels here
  bool hasFine   = false;
  bool hasCoar   = false;
  int bogRef = 2;
  //hook for optimization.
  //need to put this in  ala EBAMRPoissonOp to get it.
  bool layoutChanged = true;
  //
  EBSpeciesFluxOp* newOp = NULL;
  if (acoef0 == NULL)
  {
    // Time-independent a coefficient.
    newOp = new EBSpeciesFluxOp(EBLevelGrid(), eblgMGLevel, EBLevelGrid(), eblgCoarMG, quadCFI,
                                 dombc, ebbc, dxMGLevel, dxCoar, bogRef, bogRef, hasFine, hasCoar,
                                 hasCoarMGObjects, layoutChanged, m_alpha, m_beta,
                                 acoef, Dk, vecW, DkIrreg, vecWIrreg, m_ghostCellsPhi, m_ghostCellsRhs, m_relaxType, m_ispec);
  }
  else
  {
    // Time-dependent a coefficient.
    MayDay::Error("Time dependent Species Flux Operator Not Implemented");
  }

  return newOp;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
EBSpeciesFluxOp*
EBSpeciesFluxOpFactory::
AMRnewOp(const ProblemDomain& a_domainFine)
{
  //figure out which level we are at.
  int ref=-1;
  bool found = false;
  EBLevelGrid eblgFine, eblgCoar;
  for (int ilev = 0; ilev < m_numLevels && !found; ilev++)
    {
      if (a_domainFine == m_eblgs[ilev].getDomain())
        {
          found = true;
          ref = ilev ;
        }
    }
  if (!found)
    {
      MayDay::Error("No corresponding AMRLevel to starting point of AMRnewOp");
    }

  int refToFiner   = 2;
  int refToCoarser = 2;
  Real dxCoar = -1.0;
  if (ref > 0)
    {
      eblgCoar = m_eblgs[ref-1];
      dxCoar = m_dx[ref-1];
      refToCoarser= m_refRatio[ref-1];
    }
  if (ref < m_numLevels-1)
    {
      eblgFine = m_eblgs[ref+1];
      refToFiner = m_refRatio[ref];
    }
  //creates coarse and finer info and bcs and all that
  EBLevelGrid      eblgMGLevel = m_eblgs[ref];
  Real               dxMGLevel =    m_dx[ref];

  bool hasCoarMGObjects = m_hasMGObjects[ref];
  EBLevelGrid eblgCoarMG;
  if (hasCoarMGObjects)
    {
      eblgCoarMG = m_eblgsMG[ref][1];
    }
  m_ebBCFactory->setData(m_data[ref], m_ispec, m_wallBCInfo);
  SpeciesFluxBaseEBBC*   speciesEBBC  = (SpeciesFluxBaseEBBC*)     m_ebBCFactory->create(m_eblgs[ref].getDomain(),
											 m_eblgs[ref].getEBISL(), dxMGLevel*RealVect::Unit,
											 &m_ghostCellsPhi, &m_ghostCellsRhs);
  if (m_dataBased)
    {
      speciesEBBC->setData(m_data[ref], m_ispec, m_wallBCInfo);
    }

  SpeciesFluxBaseDomainBC* speciesDomBC = (SpeciesFluxBaseDomainBC*) m_domainBCFactory->create(m_eblgs[ref].getDomain(),
											    m_eblgs[ref].getEBISL(), dxMGLevel*RealVect::Unit);
  RefCountedPtr<SpeciesFluxBaseEBBC>      ebbc(speciesEBBC);
  RefCountedPtr<SpeciesFluxBaseDomainBC> dombc(speciesDomBC);

  bool hasFine = (ref < (m_eblgs.size()-1));
  bool hasCoar = (ref > 0);
  //optimization hook.  need to store the result out of EBArith::getCoarserLayoutss
  bool layoutChanged = true;

  EBSpeciesFluxOp* newOp = NULL;
  if (m_acoef0[ref] == NULL)
  {
    // Time-independent a coefficient.
    newOp = new EBSpeciesFluxOp(eblgFine, eblgMGLevel, eblgCoar, eblgCoarMG, m_quadCFI[ref],
                                 dombc, ebbc,  dxMGLevel,dxCoar, refToFiner, refToCoarser,
                                 hasFine, hasCoar, hasCoarMGObjects,  layoutChanged,
                                 m_alpha, m_beta, m_acoef[ref], m_Dk[ref], m_vecW[ref], m_DkIrreg[ref], m_vecWIrreg[ref],
                                 m_ghostCellsPhi, m_ghostCellsRhs, m_relaxType, m_ispec);
  }
  else
  {
    // Time-dependent a coefficient.
    MayDay::Error("Time dependent Species Flux Operator Not Implemented");
  }

  return newOp;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
EBSpeciesFluxOpFactory::
reclaim(MGLevelOp<LevelData<EBCellFAB> >* a_reclaim)
{
  delete a_reclaim;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
EBSpeciesFluxOpFactory::
AMRreclaim(EBSpeciesFluxOp* a_reclaim)
{
  delete a_reclaim;
}
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"
