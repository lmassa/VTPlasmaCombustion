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
#include "EBAMRPoissonOp.H"

#include "EBSpeciesFluxOp.H"
#include "EBQuadCFInterp.H"
#include "EBLevelCCProjector.H" // for averaging

#include "EBConductivityOpF_F.H"
#include "EBSpeciesFluxOpF_F.H"
#include "InterpF_F.H"
#include "EBAMRPoissonOpF_F.H"
#include "BCFunc.H"
#include "CH_Timer.H"
#include "BCFunc.H"
#include "EBLevelGrid.H"
#include "EBAlias.H"
#include "EBCoarseAverage.H"
#include "ParmParse.H"
#include "NamespaceHeader.H"
//IntVect EBSpeciesFluxOp::s_ivDebug = IntVect(D_DECL(111, 124, 3));
bool EBSpeciesFluxOp::s_turnOffBCs = false; //REALLY needs to default to false
int EBSpeciesFluxOp::s_step = -1;

//-----------------------------------------------------------------------
EBSpeciesFluxOp::
EBSpeciesFluxOp(const EBLevelGrid &                                  a_eblgFine,
		const EBLevelGrid &                                  a_eblg,
		const EBLevelGrid &                                  a_eblgCoar,
		const EBLevelGrid &                                  a_eblgCoarMG,
		const RefCountedPtr<EBQuadCFInterp>&                 a_quadCFI,
		const RefCountedPtr<SpeciesFluxBaseDomainBC>&        a_domainBC,
		const RefCountedPtr<SpeciesFluxBaseEBBC>&            a_ebBC,
		const Real    &                                      a_dx,
		const Real    &                                      a_dxCoar,
		const int&                                           a_refToFine,
		const int&                                           a_refToCoar,
		const bool&                                          a_hasFine,
		const bool&                                          a_hasCoar,
		const bool&                                          a_hasMGObjects,
		const bool&                                          a_layoutChanged,
		const Real&                                          a_alpha,
		const Real&                                          a_beta,
		const RefCountedPtr<LevelData<EBCellFAB> >&          a_acoef,
		const RefCountedPtr<LevelData<EBFluxFAB> >&          a_Dk,
		const RefCountedPtr<LevelData<EBFluxFAB> >&          a_vecW,
		const RefCountedPtr<LevelData<BaseIVFAB<Real> > >&   a_DkIrreg,
		const RefCountedPtr<LevelData<BaseIVFAB<Real> > >&   a_vecWIrreg,
		const IntVect&                                       a_ghostCellsPhi,
		const IntVect&                                       a_ghostCellsRHS,
		const int&                                           a_relaxType,
		const int&                                           a_ispec):
  LevelTGAHelmOp<LevelData<EBCellFAB>, EBFluxFAB>(false), // is time-independent
  m_ispec(a_ispec),
  m_relaxType(a_relaxType),
  m_ghostCellsPhi(a_ghostCellsPhi),
  m_ghostCellsRHS(a_ghostCellsRHS),
  m_quadCFIWithCoar(a_quadCFI),
  m_eblg(a_eblg),
  m_eblgFine(),
  m_eblgCoar(),
  m_eblgCoarMG(),
  m_eblgCoarsenedFine(),
  m_domainBC(a_domainBC),
  m_ebBC(a_ebBC),
  m_dxFine(),
  m_dx(a_dx),
  m_dxCoar(a_dxCoar),
  m_acoef(a_acoef),
  m_acoef0(), // Not used since acoef is time-independent
  m_acoef1(), // Not used since acoef is time-independent
  //m_Dk(a_Dk),
  //m_vecW(a_vecW),
  m_DkIrreg(a_DkIrreg),
  m_vecWIrreg(a_vecWIrreg),
  m_alpha(a_alpha),
  m_beta(a_beta),
  m_alphaDiagWeight(),
  m_betaDiagWeight(),
  m_refToFine(a_hasFine ? a_refToFine : 1),
  m_refToCoar(a_hasCoar ? a_refToCoar : 1),
  m_hasFine(a_hasFine),
  m_hasInterpAve(false),
  m_hasCoar(a_hasCoar),
  m_ebAverage(),
  m_ebInterp(),
  m_opEBStencil(),
  m_relCoef(),
  m_vofIterIrreg(),
  m_vofIterMulti(),
  m_vofIterDomLo(),
  m_vofIterDomHi(),
  m_loCFIVS(),
  m_hiCFIVS(),
  m_fastFR(),
  m_hasMGObjects(a_hasMGObjects),
  m_layoutChanged(a_layoutChanged),
  m_ebAverageMG(),
  m_ebInterpMG(),
  m_dblCoarMG(),
  m_ebislCoarMG(),
  m_domainCoarMG(),
  m_colors()
{
  CH_TIME("EBSpeciesFluxOp::SpeciesFluxOp");
  int ncomp = 1;

  //transform the Mobility and diffusion FluxFABS in SharfetterGummel fluxes
  EBFluxFactory ebfluxfact(m_eblg.getEBISL());
  m_vecW = RefCountedPtr< LevelData<EBFluxFAB> >(new LevelData<EBFluxFAB>(m_eblg.getDBL(), 1, a_vecW->ghostVect(), ebfluxfact));
  m_Dk = RefCountedPtr< LevelData<EBFluxFAB> >(new LevelData<EBFluxFAB>(m_eblg.getDBL(), 1, a_Dk->ghostVect(), ebfluxfact));
  LevelData<EBFluxFAB>& mDk = *m_Dk;
  LevelData<EBFluxFAB>& mVw = *m_vecW;
  const LevelData<EBFluxFAB>& aDk = *a_Dk;
  const LevelData<EBFluxFAB>& aVw = *a_vecW;
  for (DataIterator dit = aVw.dataIterator(); dit.ok(); ++dit)
    {
      for (int idir = 0;idir<SpaceDim;idir++)
        {
          EBFaceFAB&       faceMobData = mVw[dit()][idir];
          EBFaceFAB&       faceDiffData = mDk[dit()][idir];
	  BaseFab<Real>&   regFaceMobData =  faceMobData.getSingleValuedFAB();
	  BaseFab<Real>&   regFaceDiffData = faceDiffData.getSingleValuedFAB();
	  const  BaseFab<Real>&   aregFaceMobData =  aVw[dit()][idir].getSingleValuedFAB();
	  const  BaseFab<Real>&   aregFaceDiffData =  aDk[dit()][idir].getSingleValuedFAB();
	  Box locBox = faceMobData.getRegion();
	  FORT_SG1(CHF_FRA1(regFaceMobData,0),
		     CHF_FRA1(regFaceDiffData,0),
		     CHF_CONST_FRA1(aregFaceMobData,0),
		     CHF_CONST_FRA1(aregFaceDiffData,0),
		     CHF_BOX(locBox),
		     CHF_CONST_REAL(a_dx));
        }
    }

  if (m_hasFine)
    {
      m_eblgFine       = a_eblgFine;
      m_dxFine         = m_dx/a_refToFine;
    }

  EBCellFactory fact(m_eblg.getEBISL());
  m_relCoef.define(m_eblg.getDBL(), 1, IntVect::Zero, fact);
  if (m_hasCoar)
    {
      m_eblgCoar       = a_eblgCoar;

      for (int idir = 0; idir < SpaceDim; idir++)
        {
          m_loCFIVS[idir].define(m_eblg.getDBL());
          m_hiCFIVS[idir].define(m_eblg.getDBL());

          for (DataIterator dit = m_eblg.getDBL().dataIterator();  dit.ok(); ++dit)
            {
              m_loCFIVS[idir][dit()].define(m_eblg.getDomain(), m_eblg.getDBL().get(dit()),
                                            m_eblg.getDBL(), idir,Side::Lo);
              m_hiCFIVS[idir][dit()].define(m_eblg.getDomain(), m_eblg.getDBL().get(dit()),
                                            m_eblg.getDBL(), idir,Side::Hi);
            }
        }

      //if this fails, then the AMR grids violate proper nesting.
      ProblemDomain domainCoarsenedFine;
      DisjointBoxLayout dblCoarsenedFine;

      int maxBoxSize = 32;
      bool dumbool;
      bool hasCoarser = EBAMRPoissonOp::getCoarserLayouts(dblCoarsenedFine,
                                                          domainCoarsenedFine,
                                                          m_eblg.getDBL(),
                                                          m_eblg.getEBISL(),
                                                          m_eblg.getDomain(),
                                                          m_refToCoar,
                                                          m_eblg.getEBIS(),
                                                          maxBoxSize, dumbool);

      //should follow from coarsenable
      if (hasCoarser)
        {
          m_eblgCoarsenedFine = EBLevelGrid(dblCoarsenedFine, domainCoarsenedFine, 4, m_eblg.getEBIS());
          m_hasInterpAve = true;
          m_ebInterp.define( m_eblg.getDBL(),     m_eblgCoar.getDBL(),
                             m_eblg.getEBISL(), m_eblgCoar.getEBISL(),
                             domainCoarsenedFine, m_refToCoar, ncomp, m_eblg.getEBIS(),
                             a_ghostCellsPhi);
          m_ebAverage.define(m_eblg.getDBL(),     m_eblgCoarsenedFine.getDBL(),
                             m_eblg.getEBISL(), m_eblgCoarsenedFine.getEBISL(),
                             domainCoarsenedFine, m_refToCoar, ncomp, m_eblg.getEBIS(),
                             a_ghostCellsRHS);
        }
    }

  //special mg objects for when we do not have
  //a coarser level or when the refinement to coarse
  //is greater than two
  //flag for when we need special MG objects
  if (m_hasMGObjects)
    {
      int mgRef = 2;
      m_eblgCoarMG = a_eblgCoarMG;

      m_ebInterpMG.define( m_eblg.getDBL(),     m_eblgCoarMG.getDBL(),
                           m_eblg.getEBISL(), m_eblgCoarMG.getEBISL(),
                           m_eblgCoarMG.getDomain(), mgRef, ncomp, m_eblg.getEBIS(),
                           a_ghostCellsPhi);
      m_ebAverageMG.define(m_eblg.getDBL(),     m_eblgCoarMG.getDBL(),
                           m_eblg.getEBISL(), m_eblgCoarMG.getEBISL(),
                           m_eblgCoarMG.getDomain() , mgRef, ncomp, m_eblg.getEBIS(),
                           a_ghostCellsRHS);

    }

  //define stencils for the operator
  defineStencils();

}

//-----------------------------------------------------------------------
EBSpeciesFluxOp::
~EBSpeciesFluxOp()
{
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
fillGrad(const LevelData<EBCellFAB>& a_phi)
{
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
finerOperatorChanged(const MGLevelOp<LevelData<EBCellFAB> >& a_operator,
                     int a_coarseningFactor)
{
  const EBSpeciesFluxOp& op =
    dynamic_cast<const EBSpeciesFluxOp&>(a_operator);

  // Perform multigrid coarsening on the operator data.
  Interval interv(0, 0); // All data is scalar.
  EBLevelGrid eblgCoar = m_eblg;
  EBLevelGrid eblgFine = op.m_eblg;
  LevelData<EBCellFAB>& acoefCoar = *m_acoef;
  const LevelData<EBCellFAB>& acoefFine = *(op.m_acoef);
  LevelData<EBFluxFAB>& DkCoar = *m_Dk;
  const LevelData<EBFluxFAB>& DkFine = *(op.m_Dk);
  LevelData<EBFluxFAB>& vecWCoar = *m_vecW;
  const LevelData<EBFluxFAB>& vecWFine = *(op.m_vecW);
  LevelData<BaseIVFAB<Real> >& DkCoarIrreg = *m_DkIrreg;
  const LevelData<BaseIVFAB<Real> >& DkFineIrreg = *(op.m_DkIrreg);
  LevelData<BaseIVFAB<Real> >& vecWCoarIrreg = *m_vecWIrreg;
  const LevelData<BaseIVFAB<Real> >& vecWFineIrreg = *(op.m_vecWIrreg);
  if (a_coarseningFactor != 1)
    {
      EBCoarseAverage averageOp(eblgFine.getDBL(), eblgCoar.getDBL(),
                                eblgFine.getEBISL(), eblgCoar.getEBISL(),
                                eblgCoar.getDomain(), a_coarseningFactor, 1,
                                eblgCoar.getEBIS());

      //MayDay::Warning("might want to figure out what harmonic averaging is in this context");
      averageOp.average(acoefCoar, acoefFine, interv);
      averageOp.average(DkCoar, DkFine, interv);
      averageOp.average(DkCoarIrreg, DkFineIrreg, interv);
      averageOp.average(vecWCoar, vecWFine, interv);
      averageOp.average(vecWCoarIrreg, vecWFineIrreg, interv);
    }

  // Handle parallel domain ghost elements.
  acoefCoar.exchange(interv);
  DkCoar.exchange(interv);
  DkCoarIrreg.exchange(interv);
  vecWCoar.exchange(interv);
  vecWCoarIrreg.exchange(interv);

  // Recompute the relaxation coefficient for the operator.
  calculateAlphaWeight();
  calculateRelaxationCoefficient();

  // Notify any observers of this change.
  notifyObserversOfChange();
}
//-----------------------------------------------------------------------
Real
EBSpeciesFluxOp::
getSafety()
{
  Real safety = 1.0;
  return safety;
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
resetACoefficient(RefCountedPtr<LevelData<EBCellFAB> >& a_acoef)
{
  // This cannot be used if the a coefficient is time-dependent!
  CH_assert(m_acoef0 == NULL);
  CH_assert(m_acoef1 == NULL);

  m_acoef = a_acoef;
  calculateAlphaWeight();
  calculateRelaxationCoefficient();
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
calculateAlphaWeight()
{
  for (DataIterator dit = m_eblg.getDBL().dataIterator();  dit.ok(); ++dit)
    {
      VoFIterator& vofit = m_vofIterIrreg[dit()];
      for (vofit.reset(); vofit.ok(); ++vofit)
        {
          const VolIndex& VoF = vofit();
          Real volFrac = m_eblg.getEBISL()[dit()].volFrac(VoF);
          Real alphaWeight = (*m_acoef)[dit()](VoF, 0);
          alphaWeight *= volFrac;

          m_alphaDiagWeight[dit()](VoF, 0) = alphaWeight;
        }
    }
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
getDivFStencil(VoFStencil&      a_vofStencil,
               const VolIndex&  a_vof,
               const DataIndex& a_dit)
{
  CH_TIME("EBSpeciesFluxOp::getDivFStencil");
  const EBISBox& ebisBox = m_eblg.getEBISL()[a_dit];
  a_vofStencil.clear();
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      for (SideIterator sit; sit.ok(); ++sit)
        {
          int isign = sign(sit());
          Vector<FaceIndex> faces = ebisBox.getFaces(a_vof, idir, sit());
          for (int iface = 0; iface < faces.size(); iface++)
            {
              VoFStencil fluxStencil;
              getFluxStencil(fluxStencil, faces[iface], a_dit);
              Real areaFrac = ebisBox.areaFrac(faces[iface]);
              fluxStencil *= Real(isign)*areaFrac/m_dx;
              a_vofStencil += fluxStencil;
            }
        }
    }
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
getFluxStencil(VoFStencil&      a_fluxStencil,
               const FaceIndex& a_face,
               const DataIndex& a_dit)
{
  /// stencil for flux computation.   the truly ugly part of this computation
  /// beta and eta are multiplied in here

  CH_TIME("EBSpeciesFluxOp::getFluxStencil");
  //need to do this by interpolating to centroids
  //so get the stencil at each face center and add with
  //interpolation weights
  FaceStencil interpSten = EBArith::getInterpStencil(a_face,
                                                     (*m_eblg.getCFIVS())[a_dit],
                                                     m_eblg.getEBISL()[a_dit],
                                                     m_eblg.getDomain());

  a_fluxStencil.clear();
  for (int isten = 0; isten < interpSten.size(); isten++)
    {
      const FaceIndex& face = interpSten.face(isten);
      const Real&    weight = interpSten.weight(isten);
      VoFStencil faceCentSten;
      getFaceCenteredFluxStencil(faceCentSten, face, a_dit);
      faceCentSten *= weight;
      a_fluxStencil += faceCentSten;
    }
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
getFaceCenteredFluxStencil(VoFStencil&      a_fluxStencil,
                           const FaceIndex& a_face,
                           const DataIndex& a_dit)
{
  CH_TIME("EBSpeciesFluxOp::getFaceCenteredFluxStencil");
  //face centered gradient is just a centered diff
  int faceDir= a_face.direction();
  a_fluxStencil.clear();
  //LM no beta here?
  if (!a_face.isBoundary())
    {
      Real DkFace = (*m_Dk)[a_dit][faceDir](a_face, 0);
      Real vecWFace = (*m_vecW)[a_dit][faceDir](a_face, 0);
      a_fluxStencil.add(a_face.getVoF(Side::Hi),  DkFace/m_dx, 0);
      a_fluxStencil.add(a_face.getVoF(Side::Lo), -vecWFace/m_dx, 0);
    }
  else
    {
      //the boundary condition handles this one.
    }
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
setAlphaAndBeta(const Real& a_alpha,
                const Real& a_beta)
{
  CH_TIME("EBSpeciesFluxOp::setAlphaAndBeta");
  m_alpha = a_alpha;
  m_beta  = a_beta;
  calculateAlphaWeight(); //need to do this because the a coef has probably been changed under us
  calculateRelaxationCoefficient();
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
setTime(Real a_oldTime, Real a_mu, Real a_dt)
{
  // This only affects a coefficients that are time-dependent.
  if (m_acoef0 != NULL)
    {
      CH_assert(m_acoef1 != NULL); // All or nothing!

      // The a coefficient is linearly interpolated as
      // acoef = acoef0 + mu * (acoef1 - acoef0)
      //       = (1 - mu) * acoef0 + mu * acoef1.
      for (DataIterator dit = m_eblg.getDBL().dataIterator();  dit.ok(); ++dit)
        {
          (*m_acoef)[dit()].axby((*m_acoef0)[dit()], (*m_acoef1)[dit()],
                                 1.0 - a_mu, a_mu);
        }

      // Notify any observers that the operator's state has changed.
      notifyObserversOfChange();

      // Redefine the stencils.
      defineStencils();
    }
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
kappaScale(LevelData<EBCellFAB> & a_rhs)
{
  CH_TIME("EBSpeciesFluxOp::kappaScale");
  EBLevelDataOps::kappaWeight(a_rhs);
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
diagonalScale(LevelData<EBCellFAB> & a_rhs,
              bool a_kappaWeighted)
{

  CH_TIME("EBSpeciesFluxOp::diagonalScale");
  //  dumpLevelPoint(a_rhs, string("EBSpeciesFluxOp: diagonalScale: phi coming in = "));
  if (a_kappaWeighted)
    EBLevelDataOps::kappaWeight(a_rhs);

  //  dumpLevelPoint(a_rhs, string("EBSpeciesFluxOp: diagonalScale: kappa*phi = "));

  //also have to weight by the coefficient
  for (DataIterator dit = m_eblg.getDBL().dataIterator();  dit.ok(); ++dit)
    {
      a_rhs[dit()] *= (*m_acoef)[dit()];
    }
  //  dumpLevelPoint(a_rhs, string("EBSpeciesFluxOp: diagonalScale: acoef*kappa*phi = "));

}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
divideByIdentityCoef(LevelData<EBCellFAB> & a_rhs)
{

  CH_TIME("EBSpeciesFluxOp::divideByIdentityCoef");

  for (DataIterator dit = m_eblg.getDBL().dataIterator();  dit.ok(); ++dit)
    {
      a_rhs[dit()] /= (*m_acoef)[dit()];
    }
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
calculateRelaxationCoefficient()
{
  CH_TIME("ebsfo::calculateRelCoef");
  // define regular relaxation coefficent
  Real safety = getSafety();
  int ncomp = 1;
  for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      const Box& grid = m_eblg.getDBL().get(dit());

      //lucaNote assume time-independent operators
      // For time-independent acoef, initialize lambda = alpha * acoef.
      const EBCellFAB& acofab = (*m_acoef)[dit()];
      m_relCoef[dit()].setVal(0.);
      m_relCoef[dit()].plus(acofab, 0, 0, 1);
      m_relCoef[dit()]*= m_alpha;

      // Compute the relaxation coefficient in regular cells.
      BaseFab<Real>& regRel = m_relCoef[dit()].getSingleValuedFAB();
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          BaseFab<Real>& regDk = (*m_Dk)[dit()][idir].getSingleValuedFAB();
          BaseFab<Real>& regVW = (*m_vecW)[dit()][idir].getSingleValuedFAB();
          FORT_DECRINVRELCOEFEBSFO(CHF_FRA1(regRel,0),
                                  CHF_FRA1(regDk,0),
                                  CHF_FRA1(regVW,0),
                                  CHF_CONST_REAL(m_beta),
                                  CHF_BOX(grid),
                                  CHF_REAL(m_dx),
                                  CHF_INT(idir));
        }

      //now invert so lambda = stable lambda for variable coef lapl
      //(according to phil, this is the correct lambda)
      FORT_INVERTLAMBDAEBCO(CHF_FRA1(regRel,0), CHF_REAL(safety), CHF_BOX(grid));

      // Now go over the irregular cells.
      VoFIterator& vofit = m_vofIterIrreg[dit()];
      for (vofit.reset(); vofit.ok(); ++vofit)
        {
          const VolIndex& VoF = vofit();

          Real alphaWeight = m_alphaDiagWeight[dit()](VoF, 0);
          Real  betaWeight =  m_betaDiagWeight[dit()](VoF, 0);
          alphaWeight *= m_alpha;
          betaWeight  *= m_beta;

          Real diagWeight  = alphaWeight + betaWeight;

          // Set the irregular relaxation coefficients.
          if (Abs(diagWeight) > 1.0e-9)
            {
              m_relCoef[dit()](VoF, 0) = safety/diagWeight;
            }
          else
            {
              m_relCoef[dit()](VoF, 0) = 0.;
            }
        }
    }
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
defineStencils()
{
  CH_TIME("EBSpeciesFluxOp::defineStencils");
  // create ebstencil for irregular applyOp
  m_opEBStencil.define(m_eblg.getDBL());
  // create vofstencils for applyOp and

  Real fakeBeta = 1;
  m_domainBC->setCoef(m_eblg,   fakeBeta,      m_Dk,      m_vecW);
  m_ebBC->setCoef(    m_eblg,   fakeBeta,      m_DkIrreg, m_vecWIrreg);

  Real dxScale = 1.0/m_dx;
  m_ebBC->define((*m_eblg.getCFIVS()), dxScale); //has to happen AFTER coefs are set
  LayoutData<BaseIVFAB<VoFStencil> >* fluxStencil = m_ebBC->getFluxStencil(0);

  m_vofIterIrreg.define(     m_eblg.getDBL()); // vofiterator cache
  m_vofIterMulti.define(     m_eblg.getDBL()); // vofiterator cache
  m_alphaDiagWeight.define(  m_eblg.getDBL());
  m_betaDiagWeight.define(   m_eblg.getDBL());
  Box sideBoxLo[SpaceDim];
  Box sideBoxHi[SpaceDim];
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      Box domainBox = m_eblg.getDomain().domainBox();
      sideBoxLo[idir] = adjCellLo(domainBox, idir, 1);
      sideBoxLo[idir].shift(idir,  1);
      sideBoxHi[idir] = adjCellHi(domainBox, idir, 1);
      sideBoxHi[idir].shift(idir, -1);
      m_vofIterDomLo[idir].define( m_eblg.getDBL()); // vofiterator cache for domain lo
      m_vofIterDomHi[idir].define( m_eblg.getDBL()); // vofiterator cache for domain hi
    }
  EBArith::getMultiColors(m_colors);

  for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      const Box& curBox = m_eblg.getDBL().get(dit());
      const EBISBox& ebisBox = m_eblg.getEBISL()[dit()];
      const EBGraph& ebgraph = ebisBox.getEBGraph();

      IntVectSet irregIVS = ebisBox.getIrregIVS(curBox);
      IntVectSet multiIVS = ebisBox.getMultiCells(curBox);

      BaseIVFAB<VoFStencil> opStencil(irregIVS,ebgraph, 1);

      //cache the vofIterators
      m_alphaDiagWeight[dit()].define(irregIVS,ebisBox.getEBGraph(), 1);
      m_betaDiagWeight [dit()].define(irregIVS,ebisBox.getEBGraph(), 1);
      m_vofIterIrreg   [dit()].define(irregIVS,ebisBox.getEBGraph());
      m_vofIterMulti   [dit()].define(multiIVS,ebisBox.getEBGraph());

      for (int idir = 0; idir < SpaceDim; idir++)
        {
          IntVectSet loIrreg = irregIVS;
          IntVectSet hiIrreg = irregIVS;
          loIrreg &= sideBoxLo[idir];
          hiIrreg &= sideBoxHi[idir];
          m_vofIterDomLo[idir][dit()].define(loIrreg,ebisBox.getEBGraph());
          m_vofIterDomHi[idir][dit()].define(hiIrreg,ebisBox.getEBGraph());
        }

      VoFIterator& vofit = m_vofIterIrreg[dit()];
      for (vofit.reset(); vofit.ok(); ++vofit)
        {
          const VolIndex& VoF = vofit();

          VoFStencil& curStencil = opStencil(VoF,0);

          //Dk is included here in the flux consistent
          //with the regular
	  //LMnote:beta is not in there (altough It claims otherwise), so sten is on the RHS
	  // also it does all the directions
          getDivFStencil(curStencil,VoF, dit());
          if (fluxStencil != NULL)
            {
              BaseIVFAB<VoFStencil>& fluxStencilBaseIVFAB = (*fluxStencil)[dit()];
              //this fills the stencil with the gradient*beta*Dk
              VoFStencil  fluxStencilPt = fluxStencilBaseIVFAB(VoF,0);
              curStencil += fluxStencilPt;
            }
          Real betaWeight = EBArith::getDiagWeight(curStencil, VoF);
	  // LM: The betaWeight block is horrible, it shouldbe handled by domBC
	  // here it adds the contribution to n_i from faces on the dom boundary 
	  // which are not accounted for by getDivFStencil
          const IntVect& iv = VoF.gridIndex();
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              Box loSide = bdryLo(m_eblg.getDomain(),idir);
              loSide.shiftHalf(idir,1);
              Real adjust = 0;
              if (loSide.contains(iv))
                {
                  Real faceAreaFrac = 0.0;
                  Real weightedAreaFrac = 0;
                  Vector<FaceIndex> faces = ebisBox.getFaces(VoF,idir,Side::Lo);
                  for (int i = 0; i < faces.size(); i++)
                    {
                      weightedAreaFrac = ebisBox.areaFrac(faces[i]);
                      weightedAreaFrac *= (*m_Dk)[dit()][idir](faces[i],0);
                      faceAreaFrac +=  weightedAreaFrac;
                    }
		  if(faces.size() > 0)
		    adjust += -weightedAreaFrac /(m_dx*m_dx);
                }
              Box hiSide = bdryHi(m_eblg.getDomain(),idir);
              hiSide.shiftHalf(idir,-1);
              if (hiSide.contains(iv))
                {
                  Real faceAreaFrac = 0.0;
                  Real weightedAreaFrac = 0;
                  Vector<FaceIndex> faces = ebisBox.getFaces(VoF,idir,Side::Hi);
                  for (int i = 0; i < faces.size(); i++)
                    {
                      weightedAreaFrac = ebisBox.areaFrac(faces[i]);
                      weightedAreaFrac *= (*m_vecW)[dit()][idir](faces[i],0);
                      faceAreaFrac +=  weightedAreaFrac;
                    }
		  if(faces.size() > 0)
		    adjust += -weightedAreaFrac /(m_dx*m_dx) ;
                }
              betaWeight += adjust;
            }

          //add in identity term
          Real volFrac = ebisBox.volFrac(VoF);
          Real alphaWeight = (*m_acoef)[dit()](VoF, 0);
          alphaWeight *= volFrac;

          m_alphaDiagWeight[dit()](VoF, 0) = alphaWeight;
          m_betaDiagWeight[dit()](VoF, 0)  = betaWeight;
        }

      //Operator ebstencil
      m_opEBStencil[dit()] = RefCountedPtr<EBStencil>
        (new EBStencil(m_vofIterIrreg[dit()].getVector(), opStencil, m_eblg.getDBL().get(dit()),
                       m_eblg.getEBISL()[dit()], m_ghostCellsPhi, m_ghostCellsRHS, 0, true));
    }//dit
  calculateAlphaWeight();
  calculateRelaxationCoefficient();

  if (m_hasFine)
    {

      int ncomp = 1;
      m_fastFR.define(m_eblgFine, m_eblg, m_refToFine, ncomp);
      m_hasEBCF = m_fastFR.hasEBCF();
    }
  defineEBCFStencils();
  defineColorStencils(sideBoxLo, sideBoxHi);
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
defineColorStencils(Box a_sideBoxLo[SpaceDim],
                    Box a_sideBoxHi[SpaceDim])
{
  LayoutData<BaseIVFAB<VoFStencil> >* fluxStencil = m_ebBC->getFluxStencil(0);
  //define the stencils and iterators specific to gsrb.
  for (int icolor=0; icolor < m_colors.size(); ++icolor)
    {
      m_colorEBStencil[icolor].define(m_eblg.getDBL());
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          m_vofItIrregColorDomLo[icolor][idir].define( m_eblg.getDBL());
          m_vofItIrregColorDomHi[icolor][idir].define( m_eblg.getDBL());
          m_cacheEBxDomainFluxLo[icolor][idir].define( m_eblg.getDBL());
          m_cacheEBxDomainFluxHi[icolor][idir].define( m_eblg.getDBL());
        }
      for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
        {

          const EBISBox& curEBISBox = m_eblg.getEBISL()[dit()];
          const EBGraph& curEBGraph = curEBISBox.getEBGraph();
          Box dblBox( m_eblg.getDBL().get(dit()) );

          IntVectSet ivsColor(DenseIntVectSet(dblBox, false));

          VoFIterator& vofit = m_vofIterIrreg[dit()];
          for (vofit.reset(); vofit.ok(); ++vofit)
            {
              const VolIndex& vof = vofit();
              const IntVect& iv = vof.gridIndex();

              bool doThisVoF = true;
              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  if (iv[idir] % 2 != m_colors[icolor][idir])
                    {
                      doThisVoF = false;
                      break;
                    }
                }

              if (doThisVoF)
                {
                  ivsColor |= iv;
                }
            }

          BaseIVFAB<VoFStencil> colorStencilBaseIVFAB(ivsColor, curEBGraph, 1);

          for (int idir = 0; idir < SpaceDim; idir++)
            {
              IntVectSet loIrregColor = ivsColor;
              IntVectSet hiIrregColor = ivsColor;
              loIrregColor &= a_sideBoxLo[idir];
              hiIrregColor &= a_sideBoxHi[idir];
              m_vofItIrregColorDomLo[icolor][idir][dit()].define(loIrregColor,curEBGraph);
              m_vofItIrregColorDomHi[icolor][idir][dit()].define(hiIrregColor,curEBGraph);
              m_cacheEBxDomainFluxLo[icolor][idir][dit()].define(loIrregColor,curEBGraph,1);
              m_cacheEBxDomainFluxHi[icolor][idir][dit()].define(hiIrregColor,curEBGraph,1);
            }

          VoFIterator vofitcolor(ivsColor, curEBGraph);
          int ivof = 0;
          for (vofitcolor.reset(); vofitcolor.ok(); ++vofitcolor)
            {
              const VolIndex& vof = vofitcolor();

              VoFStencil& colorStencil =   colorStencilBaseIVFAB(vof,0);
              getDivFStencil(colorStencil,vof, dit());

              if (fluxStencil != NULL)
                {
                  BaseIVFAB<VoFStencil>& ebFluxStencilBaseIVFAB = (*fluxStencil)[dit()];
                  //this fills the stencil with the gradient
                  const VoFStencil&  ebFluxStencilPt = ebFluxStencilBaseIVFAB(vof,0);
                  //if the stencil returns empty, this means that our
                  //geometry is underresolved and we just set the
                  //stencil to zero.   This might not work.
                  if (ebFluxStencilPt.size() != 0)
                    {
                      colorStencil += ebFluxStencilPt;
                    }
                }
              ivof++;
            }

          const Vector<VolIndex>& srcVofs = vofitcolor.getVector();

          m_colorEBStencil[icolor][dit()]  = RefCountedPtr<EBStencil>
            (new EBStencil(srcVofs, colorStencilBaseIVFAB,  m_eblg.getDBL().get(dit()),
                           m_eblg.getEBISL()[dit()], m_ghostCellsPhi, m_ghostCellsRHS, 0, true));
        }//dit
    }//color
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
defineEBCFStencils()
{
  ///this routine is ugly and complicated.
  //I will attempt to comment it but I fear it is a lost cause
  //because the algorithm is so arcane.
  //We are attempting to only do stuff at the very specific
  //points where there is an embedded boundary crossing a coarse-fine
  //interface.   We happen to know that EBFastFR already has done this
  //choreography and we want to leverage it.

  //EBFastFR has data structures in it that serve as buffers and so on
  //that we will (thankfully) be able to leave alone.   We are only
  //going to access the data structures wherein it identifies which
  //coarse cells are associated with the coarse-fine interface
  // where the EB crosses and use this list to build up  which faces
  // need to be cal

  //important factoid: beta gets multiplied in at the last minute
  //(on evaluation) because it can change as diffusion solvers progress.
  if (m_hasFine && m_hasEBCF)
    {
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          for (SideIterator sit; sit.ok(); ++sit)
            {
              //coarse fine stuff is between me and next finer level
              //fine stuff lives over m_eblgfine
              //coar stuff lives over m_eblg
              int index = m_fastFR.index(idir, sit());
              m_stencilCoar[index].define(m_eblg.getDBL());
              m_faceitCoar [index].define(m_eblg.getDBL());

              for (DataIterator dit =      m_eblg.getDBL().dataIterator();  dit.ok(); ++dit)
                {
                  Vector<FaceIndex>& facesEBCFCoar =  m_faceitCoar[index][dit()];
                  Vector<VoFStencil>& stencEBCFCoar= m_stencilCoar[index][dit()];
                  Vector<VoFIterator>& vofitlist = m_fastFR.getVoFItCoar(dit(), idir, sit());
                  //first build up the list of the faces
                  for (int ivofit = 0; ivofit < vofitlist.size(); ivofit++)
                    {
                      VoFIterator& vofit = vofitlist[ivofit];
                      for (vofit.reset(); vofit.ok(); ++vofit)
                        {
                          //on the coarse side of the CF interface we are
                          //looking in the flip direction because we look
                          //back at the interface
                          Vector<FaceIndex> facespt = m_eblg.getEBISL()[dit()].getFaces(vofit(), idir, flip(sit()));
                          facesEBCFCoar.append(facespt);
                        }
                    }

                  stencEBCFCoar.resize(facesEBCFCoar.size());
                  for (int iface = 0; iface < stencEBCFCoar.size(); iface++)
                    {
                      IntVectSet cfivs; //does not apply here
                      getFluxStencil(stencEBCFCoar[iface], facesEBCFCoar[iface], dit());
                    }
                }
            }
        }
    }
}

//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
residual(LevelData<EBCellFAB>&       a_residual,
         const LevelData<EBCellFAB>& a_phi,
         const LevelData<EBCellFAB>& a_rhs,
         bool                        a_homogeneousPhysBC)
{
  CH_TIME("EBSpeciesFluxOp::residual");
  //this is a multigrid operator so only homogeneous CF BC
  //and null coar level
  CH_assert(a_residual.ghostVect() == m_ghostCellsRHS);
  CH_assert(a_phi.ghostVect() == m_ghostCellsPhi);
  applyOp(a_residual,a_phi,NULL, a_homogeneousPhysBC, true);
  axby(a_residual,a_residual,a_rhs,-1.0, 1.0);
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
applyOpNoBoundary(LevelData<EBCellFAB>&        a_opPhi,
                  const LevelData<EBCellFAB>&  a_phi)
{
  s_turnOffBCs = true;
  applyOp(a_opPhi, a_phi, true);
  s_turnOffBCs = false;
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
applyOp(LevelData<EBCellFAB>&             a_opPhi,
        const LevelData<EBCellFAB>&       a_phi,
        bool                              a_homogeneousPhysBC)
{
  //homogeneous CFBCs because that is all we can do.
  applyOp(a_opPhi, a_phi, NULL, a_homogeneousPhysBC, true);
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
applyOp(LevelData<EBCellFAB>&                    a_lhs,
        const LevelData<EBCellFAB>&              a_phi,
        const LevelData<EBCellFAB>* const        a_phiCoar,
        const bool&                              a_homogeneousPhysBC,
        const bool&                              a_homogeneousCFBC)
{
  CH_TIME("ebsfo::applyOp");
  LevelData<EBCellFAB>& phi = const_cast<LevelData<EBCellFAB>&>(a_phi);
  if (m_hasCoar && (!s_turnOffBCs))
    {
      applyCFBCs(phi, a_phiCoar, a_homogeneousCFBC);
    }
  phi.exchange(phi.interval());

  EBLevelDataOps::setToZero(a_lhs);
  incr( a_lhs, a_phi, m_alpha); //this multiplies by alpha
  for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      a_lhs[dit()].mult((*m_acoef)[dit()], 0, 0, 1);

      Box loBox[SpaceDim],hiBox[SpaceDim];
      int hasLo[SpaceDim],hasHi[SpaceDim];
      EBCellFAB      & phi = (EBCellFAB&)(a_phi[dit()]);
      EBCellFAB      & lph = a_lhs[dit()];
      //phi.setCoveredCellVal(0.0, 0);

      const BaseFab<Real>  & phiFAB = phi.getSingleValuedFAB();
      BaseFab<Real>        & lphFAB = lph.getSingleValuedFAB();
      Box dblBox = m_eblg.getDBL()[dit()];
      int nComps = 1;
      Box curPhiBox = phiFAB.box();

      if (!s_turnOffBCs)
        {
          incrOpRegularAllDirs(loBox, hiBox, hasLo, hasHi,
			       dblBox, curPhiBox, nComps,
			       lphFAB,
			       phiFAB,
			       a_homogeneousPhysBC,
			       dit()); 
        }
      else
        {
          //the all dirs code is wrong for no bcs = true
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              incrOpRegularDir(a_lhs[dit()], a_phi[dit()], a_homogeneousPhysBC, idir, dit());
            }
        }

      applyOpIrregular(a_lhs[dit()], a_phi[dit()], a_homogeneousPhysBC, dit());
    }
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
applyCFBCs(LevelData<EBCellFAB>&             a_phi,
           const LevelData<EBCellFAB>* const a_phiCoar,
           bool a_homogeneousCFBC)
{
  CH_TIMERS("EBSpeciesFluxOp::applyCFBCs");
  CH_TIMER("inhomogeneous_cfbcs_define",t1);
  CH_TIMER("inhomogeneous_cfbcs_execute",t3);
  CH_TIMER("homogeneous_cfbs",t2);
  CH_assert(a_phi.nComp() == 1);

  if (m_hasCoar)
    {
      if (!a_homogeneousCFBC)
        {
          CH_START(t1);
          if (a_phiCoar==NULL)
            {
              MayDay::Error("cannot enforce inhomogeneous CFBCs with NULL coar");
            }
          //define coarse fine interpolation object on the fly
          //because most operators do not need it
          CH_assert(a_phiCoar->nComp() == 1);
          CH_STOP(t1);

          CH_START(t3);
          Interval interv(0,0);
          m_quadCFIWithCoar->interpolate(a_phi, *a_phiCoar, interv);

          CH_STOP(t3);

        }
      else
        {
          CH_START(t2);
          applyHomogeneousCFBCs(a_phi);
          CH_STOP(t2);
        }
    }
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
applyHomogeneousCFBCs(LevelData<EBCellFAB>&   a_phi)
{
  CH_TIME("EBSpeciesFluxOp::applyHomogeneousCFBCs");
  CH_assert(a_phi.nComp() == 1);
  CH_assert( a_phi.ghostVect() >= IntVect::Unit);
  for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          for (SideIterator sit; sit.ok(); sit.next())
            {
              applyHomogeneousCFBCs(a_phi[dit()],dit(),idir,sit());
            }
        }
    }
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
applyHomogeneousCFBCs(EBCellFAB&            a_phi,
                      const DataIndex&      a_datInd,
                      int                   a_idir,
                      Side::LoHiSide        a_hiorlo)
{
  if (m_hasCoar)
    {
      CH_TIMERS("EBSpeciesFluxOp::applyHomogeneousCFBCs2");
      CH_TIMER("packed_applyHomogeneousCFBCs",t1);
      CH_TIMER("unpacked_applyHomogeneousCFBCs",t2);
      CH_assert((a_idir >= 0) && (a_idir  < SpaceDim));
      CH_assert((a_hiorlo == Side::Lo )||(a_hiorlo == Side::Hi ));
      CH_assert(a_phi.nComp() == 1);
      int ivar = 0;

      const CFIVS* cfivsPtr = NULL;

      if (a_hiorlo == Side::Lo)
        {
          cfivsPtr = &m_loCFIVS[a_idir][a_datInd];
        }
      else
        {
          cfivsPtr = &m_hiCFIVS[a_idir][a_datInd];
        }

      const IntVectSet& interpIVS = cfivsPtr->getFineIVS();
      if (cfivsPtr->isPacked() )
        {
          CH_START(t1);
          const int ihiorlo = sign(a_hiorlo);
          FORT_INTERPHOMO(CHF_FRA(a_phi.getSingleValuedFAB()),
                          CHF_BOX(cfivsPtr->packedBox()),
                          CHF_CONST_REAL(m_dx),
                          CHF_CONST_REAL(m_dxCoar),
                          CHF_CONST_INT(a_idir),
                          CHF_CONST_INT(ihiorlo));

          CH_STOP(t1);
        }
      else
        {
          if (!interpIVS.isEmpty())
            {
              CH_START(t2);
              Real halfdxcoar = m_dxCoar/2.0;
              Real halfdxfine = m_dx/2.0;
              Real xg = halfdxcoar -   halfdxfine;
              Real xc = halfdxcoar +   halfdxfine;
              Real xf = halfdxcoar + 3*halfdxfine;
              Real hf = m_dx;
              Real denom = xf*xc*hf;

              const EBISBox&  ebisBox = m_eblg.getEBISL()[a_datInd];
              const EBGraph&  ebgraph = m_eblg.getEBISL()[a_datInd].getEBGraph();
              for (VoFIterator vofit(interpIVS, ebgraph); vofit.ok(); ++vofit)
                {
                  const VolIndex& VoFGhost = vofit();

                  IntVect ivGhost  = VoFGhost.gridIndex();
                  IntVect ivClose =  ivGhost;
                  IntVect ivFar   =  ivGhost;

                  Vector<VolIndex> farVoFs;
                  Vector<VolIndex> closeVoFs = ebisBox.getVoFs(VoFGhost,
                                                               a_idir,
                                                               flip(a_hiorlo),
                                                               1);
                  bool hasClose = (closeVoFs.size() > 0);
                  bool hasFar = false;
                  Real phic = 0.0;
                  Real phif = 0.0;
                  if (hasClose)
                    {
                      const int& numClose = closeVoFs.size();
                      for (int iVof=0;iVof<numClose;iVof++)
                        {
                          const VolIndex& vofClose = closeVoFs[iVof];
                          phic += a_phi(vofClose,0);
                        }
                      phic /= Real(numClose);

                      farVoFs = ebisBox.getVoFs(VoFGhost,
                                                a_idir,
                                                flip(a_hiorlo),
                                                2);
                      hasFar   = (farVoFs.size()   > 0);
                      if (hasFar)
                        {
                          const int& numFar = farVoFs.size();
                          for (int iVof=0;iVof<numFar;iVof++)
                            {
                              const VolIndex& vofFar = farVoFs[iVof];
                              phif += a_phi(vofFar,0);
                            }
                          phif /= Real(numFar);
                        }
                    }

                  Real phiGhost;
                  if (hasClose && hasFar)
                    {
                      // quadratic interpolation  phi = ax^2 + bx + c
                      Real A = (phif*xc - phic*xf)/denom;
                      Real B = (phic*hf*xf - phif*xc*xc + phic*xf*xc)/denom;

                      phiGhost = A*xg*xg + B*xg;
                    }
                  else if (hasClose)
                    {
                      //linear interpolation
                      Real slope =  phic/xc;
                      phiGhost   =  slope*xg;
                    }
                  else
                    {
                      phiGhost = 0.0; //nothing to interpolate from
                    }
                  a_phi(VoFGhost, ivar) = phiGhost;
                }
              CH_STOP(t2);
            }
        }
    }
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
incrOpRegularDir(EBCellFAB&             a_lhs,
                 const EBCellFAB&       a_phi,
                 const bool&            a_homogeneous,
                 const int&             a_dir,
                 const DataIndex&       a_datInd)
{
  CH_TIME("ebsfo::incrOpReg");
  const Box& grid = m_eblg.getDBL()[a_datInd];
  Box domainFaces = m_eblg.getDomain().domainBox();
  domainFaces.surroundingNodes(a_dir);
  Box interiorFaces = grid;
  interiorFaces.surroundingNodes(a_dir);
  interiorFaces.grow(a_dir, 1);
  interiorFaces &=  domainFaces;
  interiorFaces.grow( a_dir, -1);

  //do flux difference for interior points
  FArrayBox interiorFlux(interiorFaces, 1);
  const FArrayBox& phi  = (FArrayBox&)(a_phi.getSingleValuedFAB());
  getFlux(interiorFlux, phi,  interiorFaces, a_dir, m_dx, a_datInd);

  Box loBox, hiBox, centerBox;
  int hasLo, hasHi;
  EBArith::loHiCenter(loBox, hasLo, hiBox, hasHi, centerBox, m_eblg.getDomain(),grid, a_dir);

  //do the low high center thing
  BaseFab<Real>& reglhs        = a_lhs.getSingleValuedFAB();
  Box dummyBox(IntVect::Zero, IntVect::Unit);
  FArrayBox domainFluxLo(dummyBox, 1);
  FArrayBox domainFluxHi(dummyBox, 1);

  RealVect origin = RealVect::Zero;
  Real time = 0.0;
  RealVect dxVect = m_dx*RealVect::Unit;
  if (hasLo==1)
    {
      Box loBoxFace = loBox;
      loBoxFace.shiftHalf(a_dir, -1);
      domainFluxLo.resize(loBoxFace, 1);
      if (!s_turnOffBCs)
        {
          //using EBAMRPoissonOp BCs which require a cell centered box.
          domainFluxLo.shiftHalf(a_dir, 1);
          m_domainBC->getFaceFlux(domainFluxLo,phi,origin,dxVect,a_dir,Side::Lo,a_datInd,time,a_homogeneous);
          domainFluxLo *= m_beta;
          domainFluxLo.shiftHalf(a_dir,-1);
        }
      else
        {
          //extrapolate to domain flux if there is no bc
          for (BoxIterator boxit(loBoxFace); boxit.ok(); ++boxit)
            {
              const IntVect& iv = boxit();
              IntVect ivn= iv;
              ivn[a_dir]++;
              domainFluxLo(iv, 0) = interiorFlux(ivn, 0);
            }
        }
    }
  if (hasHi==1)
    {
      Box hiBoxFace = hiBox;
      hiBoxFace.shiftHalf(a_dir, 1);
      domainFluxHi.resize(hiBoxFace, 1);
      if (!s_turnOffBCs)
        {
          //using EBAMRPoissonOp BCs which require a cell centered box.
          domainFluxHi.shiftHalf(a_dir, -1);
          m_domainBC->getFaceFlux(domainFluxHi,phi,origin,dxVect,a_dir,Side::Hi,a_datInd,time,a_homogeneous);
          domainFluxHi *= m_beta;
          domainFluxHi.shiftHalf(a_dir,  1);
        }
      else
        {
          //extrapolate to domain flux if there is no bc
          for (BoxIterator boxit(hiBoxFace); boxit.ok(); ++boxit)
            {
              const IntVect& iv = boxit();
              IntVect ivn= iv;
              ivn[a_dir]--;
              domainFluxHi(iv, 0) = interiorFlux(ivn, 0);
            }
        }
    }
  Real unity = 1.0;  //beta already in the flux
  FORT_INCRAPPLYEBSFO(CHF_FRA1(reglhs,0),
                     CHF_CONST_FRA1(interiorFlux, 0),
                     CHF_CONST_FRA1(domainFluxLo, 0),
                     CHF_CONST_FRA1(domainFluxHi, 0),
                     CHF_CONST_REAL(unity),
                     CHF_CONST_REAL(m_dx),
                     CHF_BOX(loBox),
                     CHF_BOX(hiBox),
                     CHF_BOX(centerBox),
                     CHF_CONST_INT(hasLo),
                     CHF_CONST_INT(hasHi),
                     CHF_CONST_INT(a_dir));
}
//-----------------------------------------------------------------------
void // Called by applyOp
EBSpeciesFluxOp::
incrOpRegularAllDirs(Box * a_loBox,
                     Box * a_hiBox,
                     int * a_hasLo,
                     int * a_hasHi,
                     Box & a_curDblBox,
                     Box & a_curPhiBox,
                     int a_nComps,
                     BaseFab<Real> & a_curOpPhiFAB,
                     const BaseFab<Real> & a_curPhiFAB,
                     bool a_homogeneousPhysBC,
                     const DataIndex& a_dit)
{
  CH_TIME("EBSpeciesFluxOp::incrOpRegularAllDirs");
  CH_assert(m_domainBC != NULL);

  //need to monkey with the ghost cells to account for boundary conditions
  if (!s_turnOffBCs)
    {
      BaseFab<Real>& phiFAB = (BaseFab<Real>&) a_curPhiFAB;
      applyDomainFlux(a_loBox, a_hiBox, a_hasLo, a_hasHi,
                      a_curDblBox, a_nComps, phiFAB,
                      a_homogeneousPhysBC, a_dit);
    }

  for (int comp = 0; comp<a_nComps; comp++)
    {
      //data ptr fusses if it is truly zero size
      BaseFab<Real> dummy(Box(IntVect::Zero, IntVect::Zero), 1);

      BaseFab<Real>* bc[3];  BaseFab<Real>* Wc[3];
      //need three coeffs because this has to work in 3d
      //this is my klunky way to make the call dimension-independent
      for (int iloc = 0; iloc < 3; iloc++)
        {
          if (iloc >= SpaceDim)
            {
              bc[iloc]= &dummy;
              Wc[iloc]= &dummy;
            }
          else
            {
              bc[iloc] = &((*m_Dk  )[a_dit][iloc].getSingleValuedFAB());
              Wc[iloc] = &((*m_vecW)[a_dit][iloc].getSingleValuedFAB());
            }
        }
      FORT_SPECIESFLUXINPLACE(CHF_FRA1(a_curOpPhiFAB,comp),
                               CHF_CONST_FRA1(a_curPhiFAB,comp),
                               CHF_CONST_FRA1((*bc[0]),comp),
                               CHF_CONST_FRA1((*bc[1]),comp),
                               CHF_CONST_FRA1((*bc[2]),comp),
                               CHF_CONST_FRA1((*Wc[0]),comp),
                               CHF_CONST_FRA1((*Wc[1]),comp),
                               CHF_CONST_FRA1((*Wc[2]),comp),
                               CHF_CONST_REAL(m_beta),
                               CHF_CONST_REAL(m_dx),
                               CHF_BOX(a_curDblBox));
    }
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
applyOpIrregular(EBCellFAB&             a_lhs,
                 const EBCellFAB&       a_phi,
                 const bool&            a_homogeneous,
                 const DataIndex&       a_datInd)
{
  CH_TIME("ebsfo::applyOpIrr");

//  bool stopHere = false;
//  if (a_lhs.box().contains(IntVect(26, 1)))
//    {
//      stopHere = true;
//    }

  //  dumpFABPoint(a_lhs, a_datInd, string("EBSpeciesFluxOp::applyopirr before apply lhs="));
  RealVect vectDx = m_dx*RealVect::Unit;
  //  m_opEBStencil[a_datInd]->apply(a_lhs, a_phi,
  //  m_alphaDiagWeight[a_datInd], m_alpha, m_beta, false, s_ivDebug,
  //  EBCellFAB::s_verbose);
  m_opEBStencil[a_datInd]->apply(a_lhs, a_phi, m_alphaDiagWeight[a_datInd], m_alpha, m_beta, false);

  
  Real levelFlag = m_level; //using the time variable to pass the level for debug

  if (!a_homogeneous)
    {
      const Real factor = m_beta/m_dx; //Dk handled within applyEBFlux
      m_ebBC->applyEBFlux(a_lhs, a_phi, m_vofIterIrreg[a_datInd], (*m_eblg.getCFIVS()),
                          a_datInd, RealVect::Zero, vectDx, factor,
                          a_homogeneous, levelFlag);
    }

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      int comp = 0;
      for (m_vofIterDomLo[idir][a_datInd].reset(); m_vofIterDomLo[idir][a_datInd].ok();  ++m_vofIterDomLo[idir][a_datInd])
        {
          Real flux;
          const VolIndex& vof = m_vofIterDomLo[idir][a_datInd]();
          m_domainBC->getFaceFlux(flux,vof,comp,a_phi,
                                  RealVect::Zero,vectDx,idir,Side::Lo, a_datInd, 0.0,
                                  a_homogeneous); //@applyOpIrregular
          //area gets multiplied in by bc operator
          a_lhs(vof,comp) -= flux*m_beta/m_dx;
        }
      for (m_vofIterDomHi[idir][a_datInd].reset(); m_vofIterDomHi[idir][a_datInd].ok();  ++m_vofIterDomHi[idir][a_datInd])
        {
          Real flux;
          const VolIndex& vof = m_vofIterDomHi[idir][a_datInd]();
          m_domainBC->getFaceFlux(flux,vof,comp,a_phi,
                                  RealVect::Zero,vectDx,idir,Side::Hi,a_datInd,0.0,
                                  a_homogeneous); //@applyOpIrregular
          //area gets multiplied in by bc operator
          a_lhs(vof,comp) += flux*m_beta/m_dx;
        }
    }
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
applyDomainFlux(Box * a_loBox,
                Box * a_hiBox,
                int * a_hasLo,
                int * a_hasHi,
                Box & a_dblBox,
                int a_nComps,
                BaseFab<Real> & a_phiFAB,
                bool a_homogeneousPhysBC,
                const DataIndex& a_dit)
{
  CH_TIME("EBSpeciesFluxOp::applyDomainFlux");
  CH_assert(m_domainBC != NULL);

  for (int idir=0; idir<SpaceDim; idir++)
    {

      EBArith::loHi(a_loBox[idir], a_hasLo[idir],
                    a_hiBox[idir], a_hasHi[idir],
                    m_eblg.getDomain(),a_dblBox, idir);

      for (int comp = 0; comp<a_nComps; comp++)
        {

          if (a_hasLo[idir] == 1 )
            {
              Box lbox=a_loBox[idir];
              lbox.shift(idir,-1);
              FArrayBox loFaceFlux(a_loBox[idir],a_nComps);
              int side = -1;
              Real time = 0;
              m_domainBC->getFaceFlux(loFaceFlux,a_phiFAB,RealVect::Zero,m_dx*RealVect::Unit,idir,Side::Lo,a_dit,time,a_homogeneousPhysBC);//@applyDomainFlux

              BaseFab<Real>& bc = ((*m_Dk)[a_dit][idir].getSingleValuedFAB());
              BaseFab<Real>& Wc = ((*m_vecW)[a_dit][idir].getSingleValuedFAB());
              //again, following the odd convention of EBAMRPoissonOp
              //(because I am reusing its BC classes),
              //the input flux here is CELL centered and the input box
              //is the box adjacent to the domain boundary on the valid side.
              //because I am not insane (yet) I will just shift the flux's box
              //over and multiply by the appropriate coefficient
              bc.shiftHalf(idir, 1);
              Wc.shiftHalf(idir, 1);
              FORT_EBSFOREGAPPLYDOMAINFLUX(CHF_FRA1(a_phiFAB,comp),
                                          CHF_CONST_FRA1(loFaceFlux,comp),
                                          CHF_CONST_FRA1(bc,comp),
                                          CHF_CONST_FRA1(Wc,comp),
                                          CHF_CONST_REAL(m_dx),
                                          CHF_CONST_INT(side),
                                          CHF_CONST_INT(idir),
                                          CHF_BOX(lbox));
              bc.shiftHalf(idir, -1);
              Wc.shiftHalf(idir, -1);
            }

          if (a_hasHi[idir] == 1)
            {
              Box hbox=a_hiBox[idir];
              hbox.shift(idir,1);
              FArrayBox hiFaceFlux(a_hiBox[idir],a_nComps);
              int side = 1;
              Real time = 0;
              m_domainBC->getFaceFlux(hiFaceFlux,a_phiFAB,RealVect::Zero,m_dx*RealVect::Unit,idir,Side::Hi,a_dit,time,a_homogeneousPhysBC);

              BaseFab<Real>& bc = ((*m_Dk)  [a_dit][idir].getSingleValuedFAB());
              BaseFab<Real>& Wc = ((*m_vecW)[a_dit][idir].getSingleValuedFAB());
              //again, following the odd convention of EBAMRPoissonOp
              //(because I am reusing its BC classes),
              //the input flux here is CELL centered and the input box
              //is the box adjacent to the domain boundary on the valid side.
              //because I am not insane (yet) I will just shift the flux's box
              //over and multiply by the appropriate coefficient
              bc.shiftHalf(idir, -1);
              Wc.shiftHalf(idir, -1);
              FORT_EBSFOREGAPPLYDOMAINFLUX(CHF_FRA1(a_phiFAB,comp),
                                          CHF_CONST_FRA1(hiFaceFlux,comp),
                                          CHF_CONST_FRA1(bc,comp),
                                          CHF_CONST_FRA1(Wc,comp),
                                          CHF_CONST_REAL(m_dx),
                                          CHF_CONST_INT(side),
                                          CHF_CONST_INT(idir),
                                          CHF_BOX(hbox));
              bc.shiftHalf(idir, 1);
              Wc.shiftHalf(idir, 1);

            }

        }
    }
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
preCond(LevelData<EBCellFAB>&       a_lhs,
        const LevelData<EBCellFAB>& a_rhs)
{
  CH_TIME("EBSpeciesFluxOp::preCond");
  EBLevelDataOps::assign(a_lhs, a_rhs);
  EBLevelDataOps::scale(a_lhs, m_relCoef);

  relax(a_lhs, a_rhs, 40);
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
relax(LevelData<EBCellFAB>&       a_phi,
      const LevelData<EBCellFAB>& a_rhs,
      int                         a_iterations)
{
  CH_TIME("ebsfo::relax");
  if (m_relaxType == 0)
    {
      relaxPoiJac(a_phi, a_rhs, a_iterations);
    }
  else if (m_relaxType == 1)
    {
      relaxGauSai(a_phi, a_rhs, a_iterations);
    }
  else if (m_relaxType == 2)
    {
      relaxGSRBFast(a_phi, a_rhs, a_iterations);
    }
  else
    {
      MayDay::Error("EBSpeciesFluxOp::bogus relaxtype");
    }
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
relaxPoiJac(LevelData<EBCellFAB>&       a_phi,
            const LevelData<EBCellFAB>& a_rhs,
            int                         a_iterations)
{
  CH_TIME("EBSpeciesFluxOp::relaxPoiJac");

  CH_assert(a_phi.ghostVect() == m_ghostCellsPhi);
  CH_assert(a_rhs.ghostVect() == m_ghostCellsRHS);

  CH_assert(a_phi.nComp() == 1);
  CH_assert(a_rhs.nComp() == 1);

  LevelData<EBCellFAB> lphi;
  create(lphi, a_rhs);
  for (int whichIter =0; whichIter < a_iterations; whichIter++)
    {
      applyHomogeneousCFBCs(a_phi);

      //after this lphi = L(phi)
      //this call contains bcs and exchange
      applyOp(  lphi,  a_phi, true);

      for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
        {
          lphi[dit()] -=     a_rhs[dit()];
          lphi[dit()] *= m_relCoef[dit()];
          //this is a safety factor because pt jacobi needs a smaller relaxation param
          lphi[dit()] *= -0.5;
          a_phi[dit()] += lphi[dit()];
        }
    }
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
relaxGSRBFast(LevelData<EBCellFAB>&       a_phi,
              const LevelData<EBCellFAB>& a_rhs,
              int                         a_iterations)
{
  CH_TIME("EBSpeciesFluxOp::relaxGSRBFast");

  CH_assert(a_phi.ghostVect() == m_ghostCellsPhi);
  CH_assert(a_rhs.ghostVect() == m_ghostCellsRHS);
  CH_assert(a_phi.nComp() == 1);
  CH_assert(a_rhs.nComp() == 1);


  for (int whichIter =0; whichIter < a_iterations; whichIter++)
    {


      //this is a multigrid operator so only homogeneous CF BC and null coar level
      CH_assert(a_rhs.ghostVect()    == m_ghostCellsRHS);
      CH_assert(a_phi.ghostVect()    == m_ghostCellsPhi);

      const DisjointBoxLayout& dbl = a_phi.disjointBoxLayout();

      int nComps = a_phi.nComp();
      int ibox = 0;
      for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit)
        {
          Box dblBox(m_eblg.getDBL().get(dit()));
          EBCellFAB& phi = a_phi[dit()];
          BaseFab<Real>& phiFAB       = phi.getSingleValuedFAB();

          Box loBox[SpaceDim],hiBox[SpaceDim];
          int hasLo[SpaceDim],hasHi[SpaceDim];

          {
            CH_TIME("EBSpeciesFluxOp::levelGSRB::applyDomainFlux");
            applyDomainFlux(loBox, hiBox, hasLo, hasHi,
                            dblBox, nComps, phiFAB,
                            true, dit());
          }
          ibox++;
        }

      // do first red, then black passes
      for (int redBlack =0; redBlack <= 1; redBlack++)
        {
          CH_TIME("EBSpeciesFluxOp::levelGSRB::Compute");

          a_phi.exchange();

          if (m_hasCoar)
            {
              CH_TIME("EBSpeciesFluxOp::levelGSRB::homogeneousCFInterp");
              applyCFBCs(a_phi, NULL, true);
            }
          ibox = 0;
          for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit)
            {
              EBCellFAB& phifab = a_phi[dit()];
              const EBCellFAB& rhsfab = a_rhs[dit()];

              //cache phi
              for (int c = 0; c < m_colors.size()/2; ++c)
                {
                  m_colorEBStencil[m_colors.size()/2*redBlack+c][dit()]->cachePhi(phifab);
                }

              //reg cells
              const Box& region = dbl.get(dit());
              Box dblBox(m_eblg.getDBL().get(dit()));
              //dummy has to be real because basefab::dataPtr is retarded
              BaseFab<Real> dummy(Box(IntVect::Zero, IntVect::Zero), 1);

              BaseFab<Real>      & reguPhi =      (a_phi[dit()]).getSingleValuedFAB();
              const BaseFab<Real>& reguRHS =     (a_rhs[dit()] ).getSingleValuedFAB();
              const BaseFab<Real>& relCoef = (m_relCoef[dit()] ).getSingleValuedFAB();
              const BaseFab<Real>& regACoe =((*m_acoef)[dit()] ).getSingleValuedFAB();
              const BaseFab<Real>* regBCoe[3];const BaseFab<Real>*  regWCoe[3];
              //need three coeffs because this has to work in 3d
              //this is my klunky way to make the call dimension-independent
              for (int iloc = 0; iloc < 3; iloc++)
                {
                  if (iloc >= SpaceDim)
                    {
                      regBCoe[iloc]= &dummy;
                      regWCoe[iloc]= &dummy;
                    }
                  else
                    {
                      regBCoe[iloc] = &((*m_Dk  )[dit()][iloc].getSingleValuedFAB());
                      regWCoe[iloc] = &((*m_vecW)[dit()][iloc].getSingleValuedFAB());
                    }
                }


              for (int comp = 0; comp < a_phi.nComp(); comp++)
                {
                  FORT_SPECIESFLUXGSRB(CHF_FRA1(        reguPhi,    comp),
                                        CHF_CONST_FRA1(  reguRHS,    comp),
                                        CHF_CONST_FRA1(  relCoef,    comp),
                                        CHF_CONST_FRA1(  regACoe,    comp),
                                        CHF_CONST_FRA1((*regBCoe[0]),comp),
                                        CHF_CONST_FRA1((*regBCoe[1]),comp),
                                        CHF_CONST_FRA1((*regBCoe[2]),comp),
                                        CHF_CONST_FRA1((*regWCoe[0]),comp),
                                        CHF_CONST_FRA1((*regWCoe[1]),comp),
                                        CHF_CONST_FRA1((*regWCoe[2]),comp),
                                        CHF_CONST_REAL(m_alpha),
                                        CHF_CONST_REAL(m_beta),
                                        CHF_CONST_REAL(m_dx),
                                        CHF_BOX(region),
                                        CHF_CONST_INT(redBlack));
                }

              //uncache phi
              for (int c = 0; c < m_colors.size()/2; ++c)
                {
                  m_colorEBStencil[m_colors.size()/2*redBlack+c][dit()]->uncachePhi(phifab);
                }

              for (int c = 0; c < m_colors.size()/2; ++c)
                {
                  GSColorAllIrregular(phifab, rhsfab, m_colors.size()/2*redBlack+c, dit());
                }
              ibox++;
            }
        }
    }
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
GSColorAllIrregular(EBCellFAB&                   a_phi,
                    const EBCellFAB&             a_rhs,
                    const int&                   a_icolor,
                    const DataIndex&             a_dit)
{
  CH_TIME("EBConductivyOp::GSColorAllIrregular");

  int comp = 0;
  Real time = 0;
  RealVect vdx = m_dx*RealVect::Unit;
  const BaseIVFAB<Real>& curAlphaWeight = m_alphaDiagWeight[a_dit];
  const BaseIVFAB<Real>& curBetaWeight  = m_betaDiagWeight[a_dit];
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      CH_TIME("domain cross eb stuff cache");
      for (m_vofItIrregColorDomLo[a_icolor][idir][a_dit].reset(); m_vofItIrregColorDomLo[a_icolor][idir][a_dit].ok();  ++m_vofItIrregColorDomLo[a_icolor][idir][a_dit])
        {
          const VolIndex& vof = m_vofItIrregColorDomLo[a_icolor][idir][a_dit]();
          Real flux;
          m_domainBC->getFaceFlux(flux,vof,comp,a_phi,
                                  RealVect::Zero,vdx,idir,Side::Lo, a_dit, time, true);

          m_cacheEBxDomainFluxLo[a_icolor][idir][a_dit](vof, comp) = flux;
        }

      for (m_vofItIrregColorDomHi[a_icolor][idir][a_dit].reset(); m_vofItIrregColorDomHi[a_icolor][idir][a_dit].ok();  ++m_vofItIrregColorDomHi[a_icolor][idir][a_dit])
        {
          const VolIndex& vof = m_vofItIrregColorDomHi[a_icolor][idir][a_dit]();
          Real flux;
          m_domainBC->getFaceFlux(flux,vof,comp,a_phi,
                                  RealVect::Zero,vdx,idir,Side::Hi,a_dit,time, true);

          m_cacheEBxDomainFluxHi[a_icolor][idir][a_dit](vof, comp) = flux;
        }
    }
  {
    CH_TIME("color ebstencil bit");
    //phi = (I-lambda*L)phiOld
    Real safety = getSafety();
    m_colorEBStencil[a_icolor][a_dit]->relax(a_phi, a_rhs, curAlphaWeight, curBetaWeight, m_alpha, m_beta, safety);
  }


  //apply domain bcs to (I-lambda*L)phi (already done in colorStencil += fluxStencil, and hom only here))
  //apply domain bcs to (I-lambda*L)phi
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      CH_TIME("domain cross eb stuff uncache");
      for (m_vofItIrregColorDomLo[a_icolor][idir][a_dit].reset(); m_vofItIrregColorDomLo[a_icolor][idir][a_dit].ok();  ++m_vofItIrregColorDomLo[a_icolor][idir][a_dit])
        {
          const VolIndex& vof = m_vofItIrregColorDomLo[a_icolor][idir][a_dit]();
          Real weightIrreg = m_alpha*curAlphaWeight(vof,0) + m_beta*curBetaWeight(vof,0);
          Real lambda = 0.0;
          if (Abs(weightIrreg) > 1.e-12)
            {
              lambda = 1./weightIrreg;
            }
          a_phi(vof,comp) += lambda * m_cacheEBxDomainFluxLo[a_icolor][idir][a_dit](vof, comp) * m_beta/m_dx;
        }

      for (m_vofItIrregColorDomHi[a_icolor][idir][a_dit].reset(); m_vofItIrregColorDomHi[a_icolor][idir][a_dit].ok();  ++m_vofItIrregColorDomHi[a_icolor][idir][a_dit])
        {
          const VolIndex& vof = m_vofItIrregColorDomHi[a_icolor][idir][a_dit]();
          Real weightIrreg = m_alpha*curAlphaWeight(vof,0) + m_beta*curBetaWeight(vof,0);
          Real lambda = 0.0;
          if (Abs(weightIrreg) > 1.e-12)
            {
              lambda = 1./weightIrreg;
            }
          a_phi(vof,comp) -= lambda * m_cacheEBxDomainFluxHi[a_icolor][idir][a_dit](vof, comp) * m_beta/m_dx;
        }
    }
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
relaxGauSai(LevelData<EBCellFAB>&       a_phi,
            const LevelData<EBCellFAB>& a_rhs,
            int                         a_iterations)
{
  CH_TIME("EBSpeciesFluxOp::relaxGauSai");

  CH_assert(a_phi.ghostVect() == m_ghostCellsPhi);
  CH_assert(a_rhs.ghostVect() == m_ghostCellsRHS);

  CH_assert(a_phi.nComp() == 1);
  CH_assert(a_rhs.nComp() == 1);

  LevelData<EBCellFAB> lphi;
  create(lphi, a_rhs);
  for (int whichIter =0; whichIter < a_iterations; whichIter++)
    {
      for (int icolor = 0; icolor < m_colors.size(); icolor++)
        {
          applyHomogeneousCFBCs(a_phi);

          //after this lphi = L(phi)
          //this call contains bcs and exchange
          applyOp(  lphi,  a_phi, true);
          gsrbColor(a_phi, lphi, a_rhs, m_colors[icolor]);
        }
    }
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
gsrbColor(LevelData<EBCellFAB>&       a_phi,
          const LevelData<EBCellFAB>& a_lph,
          const LevelData<EBCellFAB>& a_rhs,
          const IntVect&              a_color)
{
  CH_TIME("EBSpeciesFluxOp::gsrbColor");

  const DisjointBoxLayout& dbl = a_phi.disjointBoxLayout();
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
    {
      Box dblBox  = dbl.get(dit());
      BaseFab<Real>&       regPhi =     a_phi[dit()].getSingleValuedFAB();
      const BaseFab<Real>& regLph =     a_lph[dit()].getSingleValuedFAB();
      const BaseFab<Real>& regRhs =     a_rhs[dit()].getSingleValuedFAB();
      IntVect loIV = dblBox.smallEnd();
      IntVect hiIV = dblBox.bigEnd();

      for (int idir = 0; idir < SpaceDim; idir++)
        {
          if (loIV[idir] % 2 != a_color[idir])
            {
              loIV[idir]++;
            }
        }

      const BaseFab<Real>& regRel = m_relCoef[dit()].getSingleValuedFAB();
      if (loIV <= hiIV)
        {
          Box coloredBox(loIV, hiIV);
          FORT_GSRBEBSFO(CHF_FRA1(regPhi,0),
                        CHF_CONST_FRA1(regLph,0),
                        CHF_CONST_FRA1(regRhs,0),
                        CHF_CONST_FRA1(regRel,0),
                        CHF_BOX(coloredBox));
        }

      for (m_vofIterMulti[dit()].reset(); m_vofIterMulti[dit()].ok(); ++m_vofIterMulti[dit()])
        {
          const VolIndex& vof = m_vofIterMulti[dit()]();
          const IntVect& iv = vof.gridIndex();

          bool doThisVoF = true;
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              if (iv[idir] % 2 != a_color[idir])
                {
                  doThisVoF = false;
                  break;
                }
            }

          if (doThisVoF)
            {
              Real lph    = a_lph[dit()](vof, 0);
              Real rhs    = a_rhs[dit()](vof, 0);
              Real resid  = rhs - lph;
              Real lambda = m_relCoef[dit()](vof, 0);
              a_phi[dit()](vof, 0) += lambda*resid;
            }
        }
    }
}
//-----------------------------------------------------------------------
void EBSpeciesFluxOp::
AMRResidual(LevelData<EBCellFAB>&       a_residual,
            const LevelData<EBCellFAB>& a_phiFine,
            const LevelData<EBCellFAB>& a_phi,
            const LevelData<EBCellFAB>& a_phiCoar,
            const LevelData<EBCellFAB>& a_rhs,
            bool a_homogeneousPhysBC,
            AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp)
{
  CH_TIMERS("EBSpeciesFluxOp::AMRResidual");
  CH_TIMER("AMROperator", t1);
  CH_TIMER("axby", t2);
  CH_assert(a_residual.ghostVect() == m_ghostCellsRHS);
  CH_assert(a_rhs.ghostVect() == m_ghostCellsRHS);
  CH_assert(a_residual.nComp() == 1);
  CH_assert(a_phi.nComp() == 1);
  CH_assert(a_rhs.nComp() == 1);

  CH_START(t1);
  AMROperator(a_residual, a_phiFine, a_phi, a_phiCoar,
              a_homogeneousPhysBC, a_finerOp);
  CH_STOP(t1);

  //  dumpLevelPoint(a_residual, string("EBSpeciesFluxOp: AMRResidual: lphi = "));
  //  dumpLevelPoint(a_rhs,      string("EBSpeciesFluxOp: AMRResidual: rhs = "));
  //multiply by -1 so a_residual now holds -L(phi)
  //add in rhs so a_residual = rhs - L(phi)
  CH_START(t2);
  axby(a_residual,a_residual,a_rhs,-1.0, 1.0);
  CH_STOP(t2);
  //  dumpLevelPoint(a_residual, string("EBSpeciesFluxOp: AMRResidual: resid = "));
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
AMRUpdateResidual(LevelData<EBCellFAB>&       a_residual,
                  const LevelData<EBCellFAB>& a_correction,
                  const LevelData<EBCellFAB>& a_coarCorrection)
{
  CH_TIME("EBSpeciesFluxOp::AMRUpdateResidual");
  CH_assert(a_residual.ghostVect() == m_ghostCellsRHS);
  CH_assert(a_correction.ghostVect() == m_ghostCellsPhi);
  CH_assert(a_coarCorrection.ghostVect() == m_ghostCellsPhi);

  LevelData<EBCellFAB> lcorr;
  bool homogeneousPhys = true;
  bool homogeneousCF   = false;

  EBCellFactory ebcellfactTL(m_eblg.getEBISL());
  IntVect ghostVec = a_residual.ghostVect();

  lcorr.define(m_eblg.getDBL(), 1, ghostVec, ebcellfactTL);

  applyOp(lcorr, a_correction, &a_coarCorrection, homogeneousPhys, homogeneousCF);

  incr(a_residual, lcorr, -1);
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
AMRResidualNF(LevelData<EBCellFAB>&       a_residual,
              const LevelData<EBCellFAB>& a_phi,
              const LevelData<EBCellFAB>& a_phiCoar,
              const LevelData<EBCellFAB>& a_rhs,
              bool a_homogeneousPhysBC)
{
  CH_TIME("ebsfo::amrresNF");
  CH_assert(a_residual.ghostVect() == m_ghostCellsRHS);
  CH_assert(a_rhs.ghostVect() == m_ghostCellsRHS);

  AMROperatorNF(a_residual, a_phi, a_phiCoar, a_homogeneousPhysBC);
  axby(a_residual,a_residual,a_rhs,-1.0, 1.0);
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
AMRResidualNC(LevelData<EBCellFAB>&       a_residual,
              const LevelData<EBCellFAB>& a_phiFine,
              const LevelData<EBCellFAB>& a_phi,
              const LevelData<EBCellFAB>& a_rhs,
              bool a_homogeneousPhysBC,
              AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp)
{
  //dummy. there is no coarse when this is called
  CH_assert(a_residual.ghostVect() == m_ghostCellsRHS);
  CH_assert(a_rhs.ghostVect() == m_ghostCellsRHS);
  LevelData<EBCellFAB> phiC;
  AMRResidual(a_residual, a_phiFine, a_phi, phiC, a_rhs, a_homogeneousPhysBC, a_finerOp);
}
//-----------------------------------------------------------------------
void EBSpeciesFluxOp::
AMROperator(LevelData<EBCellFAB>&       a_LofPhi,
            const LevelData<EBCellFAB>& a_phiFine,
            const LevelData<EBCellFAB>& a_phi,
            const LevelData<EBCellFAB>& a_phiCoar,
            bool a_homogeneousPhysBC,
            AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp)
{
  CH_TIMERS("EBSpeciesFluxOp::AMROperator");
  CH_TIMER("applyOp", t1);
  CH_TIMER("reflux", t2);
  CH_assert(a_LofPhi.ghostVect() == m_ghostCellsRHS);
  CH_assert(a_LofPhi.nComp() == 1);
  CH_assert(a_phi.nComp() == 1);

  //apply the operator between this and the next coarser level.
  CH_START(t1);
  applyOp(a_LofPhi, a_phi, &a_phiCoar,  a_homogeneousPhysBC, false);
  CH_STOP(t1);

  //  dumpLevelPoint(a_LofPhi, string("EBSpeciesFluxOp: AMROperator: pre-reflux lphi = "));
  //now reflux to enforce flux-matching from finer levels
  if (m_hasFine)
    {
      CH_assert(a_finerOp != NULL);
      CH_START(t2);

      reflux(a_LofPhi, a_phiFine, a_phi, a_finerOp);

      CH_STOP(t2);
    }
  //  dumpLevelPoint(a_LofPhi, string("EBSpeciesFluxOp: AMROperator: post-reflux lphi = "));
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
AMROperatorNF(LevelData<EBCellFAB>&       a_LofPhi,
              const LevelData<EBCellFAB>& a_phi,
              const LevelData<EBCellFAB>& a_phiCoar,
              bool a_homogeneousPhysBC)
{
  CH_assert(a_LofPhi.ghostVect() == m_ghostCellsRHS);

  applyOp(a_LofPhi,a_phi, &a_phiCoar,  a_homogeneousPhysBC, false);

}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
AMROperatorNC(LevelData<EBCellFAB>&       a_LofPhi,
              const LevelData<EBCellFAB>& a_phiFine,
              const LevelData<EBCellFAB>& a_phi,
              bool a_homogeneousPhysBC,
              AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp)
{
  CH_TIME("ebsfo::amrOpNC");
  //dummy. there is no coarse when this is called
  CH_assert(a_LofPhi.ghostVect() == m_ghostCellsRHS);
  LevelData<EBCellFAB> phiC;
  AMROperator(a_LofPhi, a_phiFine, a_phi, phiC,
              a_homogeneousPhysBC, a_finerOp);
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
reflux(LevelData<EBCellFAB>& a_residual,
       const LevelData<EBCellFAB>& a_phiFine,
       const LevelData<EBCellFAB>& a_phi,
       AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp)
{
  CH_TIMERS("EBSpeciesFluxOp::fastReflux");
  CH_TIMER("setToZero",t2);
  CH_TIMER("incrementCoar",t3);
  CH_TIMER("incrementFine",t4);
  CH_TIMER("reflux_from_reg",t5);
  Interval interv(0,0);

  CH_START(t2);
  m_fastFR.setToZero();
  CH_STOP(t2);
  CH_START(t3);
  incrementFRCoar(m_fastFR, a_phiFine, a_phi);
  CH_STOP(t3);

  CH_START(t4);
  incrementFRFine(m_fastFR, a_phiFine, a_phi, a_finerOp);
  CH_STOP(t4);
  CH_START(t5);

  Real scale = 1.0/m_dx;
  m_fastFR.reflux(a_residual, interv, scale);

  CH_STOP(t5);
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
incrementFRCoar(EBFastFR&             a_fluxReg,
                const LevelData<EBCellFAB>& a_phiFine,
                const LevelData<EBCellFAB>& a_phi)
{
  CH_TIME("EBSpeciesFluxOp::incrementFRCoar");
  CH_assert(a_phiFine.nComp() == 1);
  CH_assert(a_phi.nComp() == 1);

  int ncomp = 1;
  Interval interv(0,0);

  for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      const EBCellFAB& coarfab = a_phi[dit()];
      const EBISBox& ebisBox = m_eblg.getEBISL()[dit()];
      const Box&  box = m_eblg.getDBL().get(dit());
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          //no boundary faces here.

          Box ghostedBox = box;
          ghostedBox.grow(1);
          ghostedBox.grow(idir,-1);
          ghostedBox &= m_eblg.getDomain();

          EBFaceFAB coarflux(ebisBox, ghostedBox, idir, ncomp);

          //old way
          //getFlux(coarflux, coarfab, ghostedBox, box, m_eblg.getDomain(), ebisBox, m_dx, dit(), idir);

          // new way
          getFluxRegOnly(coarflux, coarfab, ghostedBox, m_dx, dit(), idir);
          for (SideIterator sit; sit.ok(); ++sit)
            {
              Vector<FaceIndex>*  faceit=NULL;
              Vector<VoFStencil>* stencil=NULL;
              int index = EBFastFR::index(idir, sit());
              if (m_hasEBCF)
                {
                  faceit  = &( m_faceitCoar[index][dit()]);
                  stencil = &(m_stencilCoar[index][dit()]);
                }
              getFluxEBCF(coarflux, coarfab, ghostedBox, *faceit, *stencil);
            }

          //          dumpFlux(coarflux, idir,  string("incrementFRCoar: flux = "));
          Real scale = 1.0; //beta and Dk already in flux
          for (SideIterator sit; sit.ok(); ++sit)
            {
              a_fluxReg.incrementCoarseBoth(coarflux, scale, dit(), interv, idir, sit());
            }
        }
    }
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
incrementFRFine(EBFastFR&             a_fluxReg,
                const LevelData<EBCellFAB>& a_phiFine,
                const LevelData<EBCellFAB>& a_phi,
                AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp)
{
  CH_TIME("EBSpeciesFluxOp::incrementFRFine");
  CH_assert(a_phiFine.nComp() == 1);
  CH_assert(a_phi.nComp() == 1);
  CH_assert(m_hasFine);
  int ncomp = 1;
  Interval interv(0,0);
  EBSpeciesFluxOp& finerEBAMROp = (EBSpeciesFluxOp& )(*a_finerOp);

  //ghost cells of phiFine need to be filled
  LevelData<EBCellFAB>& phiFine = (LevelData<EBCellFAB>&) a_phiFine;
  finerEBAMROp.m_quadCFIWithCoar->interpolate(phiFine, a_phi, interv);
  phiFine.exchange(interv);

  DataIterator ditf = a_phiFine.dataIterator();
  for (ditf.reset(); ditf.ok(); ++ditf)
    {
      const Box&     boxFine = m_eblgFine.getDBL().get(ditf());
      const EBISBox& ebisBoxFine = m_eblgFine.getEBISL()[ditf()];
      const EBCellFAB& phiFine = a_phiFine[ditf()];

      for (int idir = 0; idir < SpaceDim; idir++)
        {
          for (SideIterator sit; sit.ok(); sit.next())
            {
              Box fabBox = adjCellBox(boxFine, idir, sit(), 1);
              fabBox.shift(idir, -sign(sit()));

              Box ghostedBox = fabBox;
              ghostedBox.grow(1);
              ghostedBox.grow(idir,-1);
              ghostedBox &= m_eblgFine.getDomain();

              EBFaceFAB fluxFine(ebisBoxFine, ghostedBox, idir, ncomp);
              finerEBAMROp.getFlux(fluxFine, phiFine, ghostedBox, fabBox, m_eblgFine.getDomain(),
                                   ebisBoxFine, m_dxFine, ditf(), idir);

              Real scale = 1.0; //beta and Dk already in flux

              a_fluxReg.incrementFineBoth(fluxFine, scale, ditf(), interv, idir, sit());
            }
        }
    }
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
getFlux(FArrayBox&                    a_flux,
        const FArrayBox&              a_phi,
        const Box&                    a_faceBox,
        const int&                    a_idir,
        const Real&                   a_dx,
        const DataIndex&              a_datInd)
{
  CH_TIME("ebsfo::getflux3");
  const EBFaceFAB& bcoebff = (*m_Dk)[a_datInd][a_idir];
  const EBFaceFAB& Wcoebff = (*m_vecW)[a_datInd][a_idir];
  const FArrayBox& regBCo = (const FArrayBox&)(bcoebff.getSingleValuedFAB());
  const FArrayBox& regWCo = (const FArrayBox&)(Wcoebff.getSingleValuedFAB());
  FORT_GETFLUXEBSFO(CHF_FRA1(a_flux,0),
                   CHF_CONST_FRA1(regBCo,0),
                   CHF_CONST_FRA1(regWCo,0),
                   CHF_CONST_FRA1(a_phi, 0),
                   CHF_BOX(a_faceBox),
                   CHF_CONST_REAL(a_dx),
                   CHF_CONST_REAL(m_beta),
                   CHF_CONST_INT(a_idir));

  a_flux *= m_beta;
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
getFlux(EBFluxFAB&                    a_flux,
        const EBCellFAB&              a_data,
        const Box&                    a_grid,
        const DataIndex&              a_dit,
	const int                     a_icomp)
{
  CH_TIME("ebsfo::getflux1");
  a_flux.setVal(0.);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      Box ghostedBox = a_grid;
      ghostedBox.grow(1);
      ghostedBox.grow(idir,-1);
      ghostedBox &= m_eblg.getDomain();


      getFlux(a_flux[idir], a_data, ghostedBox, a_grid,
              m_eblg.getDomain(), m_eblg.getEBISL()[a_dit], m_dx, a_dit, idir, a_icomp);

    }
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
getFlux(EBFluxFAB&                    a_flux,
        const LevelData<EBCellFAB>&   a_data,
        const Box&                    a_grid,
        const DataIndex&              a_dit,
	const int                     a_icomp)
{
  CH_TIME("ebsfo::getflux1");
  a_flux.define(m_eblg.getEBISL()[a_dit], a_grid, 1);
  a_flux.setVal(0.);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      Box ghostedBox = a_grid;
      ghostedBox.grow(1);
      ghostedBox.grow(idir,-1);
      ghostedBox &= m_eblg.getDomain();


      getFlux(a_flux[idir], a_data[a_dit], ghostedBox, a_grid,
              m_eblg.getDomain(), m_eblg.getEBISL()[a_dit], m_dx, a_dit, idir, a_icomp);

    }
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
getFlux(EBFaceFAB&                    a_fluxCentroid,
        const EBCellFAB&              a_phi,
        const Box&                    a_ghostedBox,
        const Box&                    a_fabBox,
        const ProblemDomain&          a_domain,
        const EBISBox&                a_ebisBox,
        const Real&                   a_dx,
        const DataIndex&              a_datInd,
        const int&                    a_idir,
	const int                     a_icomp)
{
  CH_TIME("ebsfo::getFlux2");
  //has some extra cells so...
  a_fluxCentroid.setVal(0.);
  Box cellBox = a_ghostedBox;
  //want only interior faces
  cellBox.grow(a_idir, 1);
  cellBox &= a_domain;
  cellBox.grow(a_idir,-1);

  Box faceBox = surroundingNodes(cellBox, a_idir);
  EBFaceFAB fluxCenter(a_ebisBox, a_ghostedBox, a_idir,1);

  //make a EBFaceFAB (including ghost cells) that will hold centered gradients
  BaseFab<Real>& regFlux  = fluxCenter.getSingleValuedFAB();
  const BaseFab<Real>& regPhi = a_phi.getSingleValuedFAB();
  const EBFaceFAB& bcoebff = (*m_Dk)[a_datInd][a_idir];
  const FArrayBox& regBCo = (const FArrayBox&)(bcoebff.getSingleValuedFAB());
  const EBFaceFAB& Wcoebff = (*m_vecW)[a_datInd][a_idir];
  const FArrayBox& regWCo = (const FArrayBox&)(Wcoebff.getSingleValuedFAB());

  FORT_GETFLUXEBSFO(CHF_FRA1(regFlux,0),
                   CHF_CONST_FRA1(regBCo,0),
                   CHF_CONST_FRA1(regWCo,0),
                   CHF_CONST_FRA1(regPhi, a_icomp),
                   CHF_BOX(faceBox),
                   CHF_CONST_REAL(a_dx),
                   CHF_CONST_REAL(m_beta),
                   CHF_CONST_INT(a_idir));


  a_fluxCentroid.copy(fluxCenter);

  // VolIndex vofdeblo1(IntVect(D_DECL(223, 247,0)), 0);
  // VolIndex vofdebhi1(IntVect(D_DECL(223, 248,0)), 0);
  // VolIndex vofdeblo2(IntVect(D_DECL(224, 247,0)), 0);
  // VolIndex vofdebhi2(IntVect(D_DECL(224, 248,0)), 0);
  // FaceIndex facedeb1(vofdeblo1, vofdebhi1, 1);
  // FaceIndex facedeb2(vofdeblo2, vofdebhi2, 1);

  IntVectSet ivsCell = a_ebisBox.getIrregIVS(cellBox);
  if (!ivsCell.isEmpty())
    {
      FaceStop::WhichFaces stopCrit = FaceStop::SurroundingNoBoundary;

      for (FaceIterator faceit(ivsCell, a_ebisBox.getEBGraph(), a_idir,stopCrit);
          faceit.ok(); ++faceit)
        {
          const FaceIndex& face = faceit();
	  //lucanote:: assuming that hi is where the normal points to
          Real phiHi = a_phi(face.getVoF(Side::Hi), a_icomp);
          Real phiLo = a_phi(face.getVoF(Side::Lo), a_icomp);
          Real fluxFace = bcoebff(face, 0)*(phiHi - phiLo)/a_dx + 
	    Wcoebff(face, 0)*(phiHi + phiLo)/2e0;
          //          if (EBCellFAB::s_verbose && ((face==facedeb1) || (face==facedeb2)))
          //            {
          //              pout() << "EBSpeciesFluxOp::getFlux at "<< face ;
          //              pout() << ", phiHi, phiLo, flux = " << phiHi << ", " << phiLo << ", "<< fluxFace << endl;
          //            }
          fluxCenter(face, 0) = fluxFace;
        }
      //interpolate from face centers to face centroids
      Box cellBox = a_fluxCentroid.getCellRegion();
      EBArith::interpolateFluxToCentroids(a_fluxCentroid,
                                          fluxCenter,
                                          a_fabBox,
                                          a_ebisBox,
                                          a_domain,
                                          a_idir);
    }

  a_fluxCentroid *= m_beta;
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
getFluxEBCF(EBFaceFAB&                    a_flux,
            const EBCellFAB&              a_phi,
            const Box&                    a_ghostedBox,
            Vector<FaceIndex>&            a_faceitEBCF,
            Vector<VoFStencil>&           a_stenEBCF)
{
  CH_TIME("EBSpeciesFluxOp::getFluxEBCF");

  //only do the evil stuff if you have a coarse-fine /  EB crossing situation

  if (m_hasEBCF)
    {
      CH_TIME("EBCF stuff");
      for (int iface = 0; iface < a_faceitEBCF.size(); iface++)
        {
          const FaceIndex& face =     a_faceitEBCF[iface];
          const VoFStencil& stencil   = a_stenEBCF[iface];
          Real fluxval = 0;
          for (int isten = 0; isten < stencil.size(); isten++)
            {
              fluxval += stencil.weight(isten)*(a_phi(stencil.vof(isten), 0));
            }
          //note the last minute beta
          a_flux(face, 0) = m_beta*fluxval;
        }
    }
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
getFluxRegOnly(EBFaceFAB&                    a_fluxCentroid,
               const EBCellFAB&              a_phi,
               const Box&                    a_ghostedBox,
               const Real&                   a_dx,
               const DataIndex&              a_datInd,
               const int&                    a_idir)
{
  CH_TIME("ebsfo::getFluxRegOnly");
  const ProblemDomain& domain = m_eblg.getDomain();

  //has some extra cells so...
  a_fluxCentroid.setVal(0.);
  int ncomp = a_phi.nComp();
  CH_assert(ncomp == a_fluxCentroid.nComp());
  Box cellBox = a_ghostedBox;
  //want only interior faces
  cellBox.grow(a_idir, 1);
  cellBox &= domain;
  cellBox.grow(a_idir,-1);

  Box faceBox = surroundingNodes(cellBox, a_idir);

  //make a EBFaceFAB (including ghost cells) that will hold centered gradients
  BaseFab<Real>& regFlux = a_fluxCentroid.getSingleValuedFAB();
  const BaseFab<Real>& regPhi = a_phi.getSingleValuedFAB();
  const EBFaceFAB& bcoebff = (*m_Dk)[a_datInd][a_idir];
  const FArrayBox& regBCo = (const FArrayBox&)(bcoebff.getSingleValuedFAB());
  const EBFaceFAB& Wcoebff = (*m_vecW)[a_datInd][a_idir];
  const FArrayBox& regWCo = (const FArrayBox&)(Wcoebff.getSingleValuedFAB());

  FORT_GETFLUXEBSFO(CHF_FRA1(regFlux,0),
                   CHF_CONST_FRA1(regBCo,0),
                   CHF_CONST_FRA1(regWCo,0),
                   CHF_CONST_FRA1(regPhi, 0),
                   CHF_BOX(faceBox),
                   CHF_CONST_REAL(a_dx),
                   CHF_CONST_REAL(m_beta),
                   CHF_CONST_INT(a_idir));

  a_fluxCentroid *= m_beta;
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
AMRRestrict(LevelData<EBCellFAB>&       a_resCoar,
            const LevelData<EBCellFAB>& a_residual,
            const LevelData<EBCellFAB>& a_correction,
            const LevelData<EBCellFAB>& a_coarCorrection,
	    bool a_skip_res)
{
  CH_TIME("EBSpeciesFluxOp::AMRRestrict");
  CH_assert(a_residual.ghostVect() == m_ghostCellsRHS);
  CH_assert(a_correction.ghostVect() == m_ghostCellsPhi);
  CH_assert(a_coarCorrection.ghostVect() == m_ghostCellsPhi);

  CH_assert(a_residual.nComp() == 1);
  CH_assert(a_resCoar.nComp() == 1);
  CH_assert(a_correction.nComp() == 1);

  LevelData<EBCellFAB> resThisLevel;
  bool homogeneousPhys = true;
  bool homogeneousCF =   false;

  EBCellFactory ebcellfactTL(m_eblg.getEBISL());
  IntVect ghostVec = a_residual.ghostVect();

  resThisLevel.define(m_eblg.getDBL(), 1, ghostVec, ebcellfactTL);
  EBLevelDataOps::setVal(resThisLevel, 0.0);

  //API says that we must average(a_residual - L(correction, coarCorrection))
  applyOp(resThisLevel, a_correction, &a_coarCorrection, homogeneousPhys, homogeneousCF);
  incr(resThisLevel, a_residual, -1.0);
  scale(resThisLevel,-1.0);

  //use our nifty averaging operator
  Interval variables(0, 0);
  CH_assert(m_hasInterpAve);
  m_ebAverage.average(a_resCoar, resThisLevel, variables);

}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
AMRProlong(LevelData<EBCellFAB>&       a_correction,
           const LevelData<EBCellFAB>& a_coarCorrection)
{
  CH_TIME("EBSpeciesFluxOp::AMRProlong");
  //use cached interpolation object
  Interval variables(0, 0);
  CH_assert(m_hasInterpAve);
  m_ebInterp.pwcInterp(a_correction, a_coarCorrection, variables);
}

//-----------------------------------------------------------------------
void EBSpeciesFluxOp::
restrictResidual(LevelData<EBCellFAB>&       a_resCoar,
                 LevelData<EBCellFAB>&       a_phiThisLevel,
                 const LevelData<EBCellFAB>& a_rhsThisLevel)
{
  CH_TIME("EBSpeciesFluxOp::restrictResidual");

  CH_assert(a_resCoar.nComp() == 1);
  CH_assert(a_phiThisLevel.nComp() == 1);
  CH_assert(a_rhsThisLevel.nComp() == 1);

  LevelData<EBCellFAB> resThisLevel;
  bool homogeneous = true;

  EBCellFactory ebcellfactTL(m_eblg.getEBISL());
  IntVect ghostVec = a_rhsThisLevel.ghostVect();

  resThisLevel.define(m_eblg.getDBL(), 1, ghostVec, ebcellfactTL);

  // Get the residual on the fine grid
  residual(resThisLevel,a_phiThisLevel,a_rhsThisLevel,homogeneous);

  // now use our nifty averaging operator
  Interval variables(0, 0);
  if (m_layoutChanged)
    {
      m_ebAverageMG.average(a_resCoar, resThisLevel, variables);
    }
  else
    {
      m_ebAverageMG.averageMG(a_resCoar, resThisLevel, variables);
    }
}
//-----------------------------------------------------------------------
void EBSpeciesFluxOp::
prolongIncrement(LevelData<EBCellFAB>&       a_phiThisLevel,
                 const LevelData<EBCellFAB>& a_correctCoar)
{
  CH_TIME("EBSpeciesFluxOp::prolongIncrement");
  Interval vars(0, 0);
  if (m_layoutChanged)
    {
      m_ebInterpMG.pwcInterp(a_phiThisLevel, a_correctCoar, vars);
    }
  else
    {
      m_ebInterpMG.pwcInterpMG(a_phiThisLevel, a_correctCoar, vars);
    }
}
//-----------------------------------------------------------------------
Real
EBSpeciesFluxOp::
AMRNorm(const LevelData<EBCellFAB>& a_coarResid,
        const LevelData<EBCellFAB>& a_fineResid,
        const int& a_refRat,
        const int& a_ord)

{
  // compute norm over all cells on coarse not covered by finer
  CH_TIME("EBSpeciesFluxOp::AMRNorm");
  MayDay::Error("never called");
  //return norm of temp
  return norm(a_coarResid, a_ord);
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
create(LevelData<EBCellFAB>&       a_lhs,
       const LevelData<EBCellFAB>& a_rhs)
{
  CH_TIME("ebsfo::create");
  int ncomp = a_rhs.nComp();
  EBCellFactory ebcellfact(m_eblg.getEBISL());
  a_lhs.define(m_eblg.getDBL(), ncomp, a_rhs.ghostVect(), ebcellfact);
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
createCoarsened(LevelData<EBCellFAB>&       a_lhs,
                const LevelData<EBCellFAB>& a_rhs,
                const int &                 a_refRat)
{
  CH_TIME("ebsfo::createCoar");
  int ncomp = a_rhs.nComp();
  IntVect ghostVect = a_rhs.ghostVect();

  CH_assert(m_eblg.getDBL().coarsenable(a_refRat));

  //fill ebislayout
  DisjointBoxLayout dblCoarsenedFine;
  coarsen(dblCoarsenedFine, m_eblg.getDBL(), a_refRat);

  EBISLayout ebislCoarsenedFine;
  IntVect ghostVec = a_rhs.ghostVect();

  ProblemDomain coarDom = coarsen(m_eblg.getDomain(), a_refRat);
  m_eblg.getEBIS()->fillEBISLayout(ebislCoarsenedFine, dblCoarsenedFine, coarDom , ghostVec[0]);
  if (m_refToCoar > 1)
    {
      ebislCoarsenedFine.setMaxRefinementRatio(m_refToCoar, m_eblg.getEBIS());
    }

  //create coarsened data
  EBCellFactory ebcellfactCoarsenedFine(ebislCoarsenedFine);
  a_lhs.define(dblCoarsenedFine, ncomp,ghostVec, ebcellfactCoarsenedFine);
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
createCoarser(LevelData<EBCellFAB>&       a_coar,
              const LevelData<EBCellFAB>& a_fine,
              bool                        a_ghosted)
{
  CH_TIME("ebsfo::createCoarser");
  CH_assert(a_fine.nComp() == 1);
  const DisjointBoxLayout& dbl = m_eblgCoarMG.getDBL();
  ProblemDomain coarDom = coarsen(m_eblg.getDomain(), 2);

  int nghost = a_fine.ghostVect()[0];
  EBISLayout coarEBISL;

  const EBIndexSpace* const ebisPtr = m_eblg.getEBIS();
  ebisPtr->fillEBISLayout(coarEBISL,
                          dbl, coarDom, nghost);

  EBCellFactory ebcellfact(coarEBISL);
  a_coar.define(dbl, 1,a_fine.ghostVect(),ebcellfact);
}
//-----------------------------------------------------------------------
int EBSpeciesFluxOp::
refToCoarser()
{
  return m_refToCoar;
}
//-----------------------------------------------------------------------
int EBSpeciesFluxOp::
refToFiner()
{
  return m_refToFine;
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
assign(LevelData<EBCellFAB>&       a_lhs,
       const LevelData<EBCellFAB>& a_rhs)
{
  CH_TIME("ebsfo::assign");
  EBLevelDataOps::assign(a_lhs,a_rhs);
}
//-----------------------------------------------------------------------
Real
EBSpeciesFluxOp::
dotProduct(const LevelData<EBCellFAB>& a_1,
           const LevelData<EBCellFAB>& a_2)
{
  CH_TIME("ebsfo::dotProd");
  ProblemDomain domain;
  Real volume;

  return EBLevelDataOps::kappaDotProduct(volume,a_1,a_2,EBLEVELDATAOPS_ALLVOFS,domain);
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
incr(LevelData<EBCellFAB>&       a_lhs,
     const LevelData<EBCellFAB>& a_x,
     Real                        a_scale)
{
  CH_TIME("ebsfo::incr");
  EBLevelDataOps::incr(a_lhs,a_x,a_scale);
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
axby(LevelData<EBCellFAB>&       a_lhs,
     const LevelData<EBCellFAB>& a_x,
     const LevelData<EBCellFAB>& a_y,
     Real                        a_a,
     Real                        a_b)
{
  CH_TIME("ebsfo::axby");
  EBLevelDataOps::axby(a_lhs,a_x,a_y,a_a,a_b);
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
scale(LevelData<EBCellFAB>& a_lhs,
      const Real&           a_scale)
{
  CH_TIME("ebsfo::scale");
  EBLevelDataOps::scale(a_lhs,a_scale);
}
//-----------------------------------------------------------------------
Real
EBSpeciesFluxOp::
norm(const LevelData<EBCellFAB>& a_rhs,
     int                         a_ord)
{
  CH_TIME("ebsfo::norm");
  //  ParmParse pp;
  //  Real thresh;
  //  pp.get("dump_threshold", thresh);
  //  pout() << endl << "EBSpeciesFluxOp::norm: res =" << endl;
  //  dumpEBLevelThresh(&a_rhs, thresh);
  //  pout() << "after dump " << endl;

  Real maxNorm =  EBAMRPoissonOp::staticMaxNorm(a_rhs, m_eblg);
#ifdef CH_MPI
  Real tmp = 1.;
  int result = MPI_Allreduce(&maxNorm, &tmp, 1, MPI_CH_REAL,
                             MPI_MAX, Chombo_MPI::comm);
  if (result != MPI_SUCCESS)
    { //bark!!!
      MayDay::Error("sorry, but I had a communcation error on norm");
    }
  maxNorm = tmp;
#endif
  return maxNorm;
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
setToZero(LevelData<EBCellFAB>& a_lhs)
{
  CH_TIME("ebsfo::setToZero");
  EBLevelDataOps::setToZero(a_lhs);
}
//-----------------------------------------------------------------------
void
EBSpeciesFluxOp::
setVal(LevelData<EBCellFAB>& a_lhs, const Real& a_value)
{
  CH_TIME("ebsfo::setVal");
  EBLevelDataOps::setVal(a_lhs, a_value);
}
//-----------------------------------------------------------------------

void EBSpeciesFluxOp::
dumpLevelPoint(const LevelData<EBCellFAB>& a_res, const string& a_blab)
{
  for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      dumpFABPoint(a_res[dit()], dit(), a_blab);
    }
}
//-----------------------------------------------------------------------

void
EBSpeciesFluxOp::
dumpFABPoint(const EBCellFAB&       a_lhs,
             const DataIndex&       a_datInd,
             const string&          a_blab)
{
  //  if (EBCellFAB::s_verbose)
  //    {
  //      if (m_eblg.getDBL().get(a_datInd).contains(s_ivDebug))
  //        {
  //          pout() << a_blab << endl;
  //          Vector<VolIndex> vofs = m_eblg.getEBISL()[a_datInd].getVoFs(s_ivDebug);
  //          for (int ivof = 0; ivof < vofs.size(); ivof++)
  //            {
  //              pout() << "vof = " << vofs[ivof].gridIndex() << ", " << vofs[ivof].cellIndex() << "--  data =" ;
  //              for (int ivar = 0; ivar < a_lhs.nComp(); ivar++)
  //                {
  //                  pout() << a_lhs(vofs[ivof], ivar) << "  ";
  //                }
  //              pout() << endl;
  //            }
  //        }
  //    }
}
//-----------------------------------------------------------------------
void EBSpeciesFluxOp::getFluxGradPhi(LevelData<EBCellFAB>&       a_fluxGphi,
				     const LevelData<EBCellFAB>& a_cellN,
				     const LevelData<EBFluxFAB>& a_gradPhi,
				     const LevelData<EBFluxFAB>&  a_velo,
				     const int                   a_ispec)
{
  int ispec;
  if(a_ispec == -1) 
    ispec = m_ispec; //default input
  else
    ispec = a_ispec;
  CH_assert(ispec == m_ispec);

  //Note!; this version does not include te advective velocity!!
  // a_velo is not used
  pout() << "a_velo not included routine incomplete"<<endl;


  const ProblemDomain& domain = m_eblg.getDomain();
  EBFluxFactory fact(m_eblg.getEBISL());
  LevelData<EBFluxFAB>   flux  (m_eblg.getDBL(), 1, a_gradPhi.ghostVect(), fact);

  

  RealVect dxv = m_dx*RealVect::Unit;
  EBLevelDataOps::setVal(a_fluxGphi,0e0,m_ispec);
  for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    { 
      Box dblBox = m_eblg.getDBL()[dit()];
      getFlux(flux[dit()], a_cellN[dit()], dblBox,  dit(), ispec);
      flux[dit()] *= a_gradPhi[dit()];
    }

  ccpExtrapolateToDomainBoundaries(flux,m_eblg.getDBL(),m_eblg.getEBISL(),domain, dxv);
  ccpAverageFaceToCellsSum(a_fluxGphi,flux,m_eblg.getDBL(),m_eblg.getEBISL(),domain, dxv);

}


/*****/
void
EBSpeciesFluxOp::
ccpAverageFaceToCellsSum(LevelData<EBCellFAB> &        a_cellData,
			 const LevelData<EBFluxFAB> &  a_macData,
			 const DisjointBoxLayout &     a_grids,
			 const EBISLayout &            a_ebisl,
			 const ProblemDomain &         a_domain,
			 const RealVect &              a_dx)
{
  LevelData<EBFluxFAB>& nonConstFlux = (LevelData<EBFluxFAB>&) a_macData;
  nonConstFlux.exchange();
  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      const EBGraph & ebgraph = a_ebisl[dit()].getEBGraph();
      const Box grid = a_grids.get(dit());
      const EBFluxFAB &  macData =  a_macData[dit()];
      EBCellFAB &       cellData = a_cellData[dit()];

      ccpAverageFaceToCellsSum(cellData, macData, ebgraph, grid, a_domain, a_dx);
    }
}
/*****/
void
EBSpeciesFluxOp::
ccpAverageFaceToCellsSum(EBCellFAB &             a_cellData,
			 const EBFluxFAB &       a_fluxData,
			 const EBGraph &         a_ebGraph,
			 const Box &             a_grid,
			 const ProblemDomain &   a_domain,
			 const RealVect &        a_dx)
{
  int icell = m_ispec;
  int icomp = 0;
  a_cellData.setVal(m_ispec,0e0);
  IntVectSet ivsIrreg = a_ebGraph.getIrregCells(a_grid);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      const EBFaceFAB & faceData = a_fluxData[idir];
      const BaseFab<Real> & regFaceData =   faceData.getSingleValuedFAB();
      BaseFab<Real> &       regCellData = a_cellData.getSingleValuedFAB();
      FORT_INCRAVEFACETOCELL( CHF_FRA1(regCellData, icell),
                             CHF_CONST_FRA1(regFaceData, icomp),
                             CHF_CONST_INT(idir),
                             CHF_BOX(a_grid));

      for (VoFIterator vofit(ivsIrreg, a_ebGraph); vofit.ok(); ++vofit)
        {
          const VolIndex & vof = vofit();
          int numFaces = 0;
          Real cellVal = 0.;
          for (SideIterator sit; sit.ok(); ++sit)
            {
              Vector<FaceIndex> faces = a_ebGraph.getFaces(vof, idir, sit());
              //if we have faces, then use them.  otherwise, need to extrapolate to covered
              if (faces.size() > 0)
                {
                  for (int iface = 0; iface < faces.size(); iface++)
                    {
                      cellVal += faceData(faces[iface], icomp);
                      numFaces++;
                    }
                }
              else
                {
                  const EBISBox&  ebisBox = a_cellData.getEBISBox();
                  Real extrapValCov = ccpGetCoveredExtrapValue(vof, idir, sit(),
                                                               faceData,
                                                               ebisBox,
                                                               a_grid,
                                                               a_domain,
                                                               a_dx,
                                                               icomp);


                  //add in covered face value as if it were a normal value
                  cellVal += extrapValCov;
                  numFaces++;
                }
            }
          if (numFaces > 1)
            {
              cellVal /= Real(numFaces);
            }
          a_cellData(vof, icell) += cellVal;
        }
    } //end loop over directions
}

void
EBSpeciesFluxOp::
checkFAB(const EBCellFAB&  a_data)
{  
  int nCons = Min(a_data.nComp(),5);  // to avoid checking all the primitive comps  
  Box box = a_data.getRegion();
  const EBISBox& ebisBox = a_data.getEBISBox();
  const EBGraph& ebgraph = ebisBox.getEBGraph();
  
  IntVectSet ivs(box);
  for (VoFIterator vofit(ivs, ebgraph);vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      IntVect gi = vof.gridIndex();
      for (int k=0;k<nCons;k++) pout() <<  "(" << a_data(vof, k)<< ", " <<ebisBox.isIrregular(gi) << " ), ";
    }
}
void
EBSpeciesFluxOp::
checkLevel(const LevelData<EBCellFAB>&  a_data, const string& a_msg)
{  
  if(!m_out.is_open()){ 
    char charstr[100];
    sprintf(charstr, "poisson%d.dbg", procID());
    m_out.open(charstr);m_dumps=0;
    m_out.precision(9);
  }
  const LevelData<EBCellFAB>&  AD = a_data;
  int nCons = Min(AD.nComp(),7);  // to avoid checking all the primitive comps  
  for (DataIterator dit = AD.dataIterator(); dit.ok(); ++dit)
    {
      Box box = a_data.box(dit());
      const EBISBox& ebisBox = AD[dit()].getEBISBox();
      IntVectSet ivs(box);
      for (VoFIterator vofit(ivs, ebisBox.getEBGraph());vofit.ok(); ++vofit)
	{
	  const VolIndex& vof = vofit();
	  IntVect gi = vof.gridIndex();
	  
	  if(ebisBox.isIrregular(gi))
	    {
	      pout()  << a_msg << vof;
	      for (int k=0;k<nCons;k++) pout()  << AD[dit()](vof, k)<< ", ";	  
	      pout() << endl;
	      m_out  << a_msg  << ", "<< m_dumps << ", " << gi[0] << ", "<< gi[1] << ", ";
	      for (int k=0;k<nCons;k++) m_out << AD[dit()](vof, k)<< ", ";	  
	      m_out << endl;
	    }
	}
    }
  if(m_dumps > 700){
    m_out.close();
    exit(0);}
  else
    m_dumps ++;
}
void
EBSpeciesFluxOp::
checkLevel(const EBCellFAB&  a_data, const string& a_msg)
{  
  if(!m_out.is_open()){ 
    char charstr[100];
    sprintf(charstr, "speciesFlux%d.dbg", procID());
    m_out.open(charstr);m_dumps=0;
    m_out.precision(9);
  }
  int nCons = Min(a_data.nComp(),7);  // to avoid checking all the primitive comps  
  Box box = a_data.getRegion();
  const EBISBox& ebisBox = a_data.getEBISBox();
  const EBGraph& ebgraph = ebisBox.getEBGraph();
  if(ebgraph.hasIrregular()) 
    {
      m_out  << a_msg  << ", " << box << endl;
    }
  else
    return;
  IntVectSet ivs(box);
  for (VoFIterator vofit(ivs, ebgraph);vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      IntVect gi = vof.gridIndex();
      if(ebisBox.isIrregular(gi))
	{
	  pout()  << a_msg << vof;
	  for (int k=0;k<nCons;k++) pout()  << a_data(vof, k)<< ", ";	  
	  pout() << endl;
	  m_out  << a_msg  << ", "<< m_dumps << ", " << gi[0] << ", "<< gi[1] << ", ";
	  for (int k=0;k<nCons;k++) m_out << a_data(vof, k)<< ", ";	  
	  m_out << endl;
	}
    }
  if(m_dumps > 700){
    m_out.close();
    exit(0);}
  else
    m_dumps ++;
}


void
EBSpeciesFluxOp::
checkLevel(const BaseIVFAB<Real>&  a_data, const string& a_msg)
{  
  if(!m_out.is_open()){ 
    char charstr[100];
    sprintf(charstr, "speciesFlux%d.dbg", procID());
    m_out.open(charstr);m_dumps=0;
    m_out.precision(9);
  }
  int nCons = Min(a_data.nComp(),7);  // to avoid checking all the primitive comps 
  const EBGraph& ebgraph = a_data.getEBGraph();
  if(ebgraph.hasIrregular()) 
    {
      m_out  << a_msg  << ", " << endl;
    }
  else
    return;
  IntVectSet ivs = a_data.getIVS();
  for (VoFIterator vofit(ivs, ebgraph);vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      IntVect gi = vof.gridIndex();
      pout()  << a_msg << vof;
      for (int k=0;k<nCons;k++) pout()  << a_data(vof, k)<< ", ";	  
      pout() << endl;
      m_out  << a_msg  << ", "<< m_dumps << ", " << gi[0] << ", "<< gi[1] << ", ";
      for (int k=0;k<nCons;k++) m_out << a_data(vof, k)<< ", ";	  
      m_out << endl;
    }
  if(m_dumps > 700){
    m_out.close();
    exit(0);}
  else
    m_dumps ++;
}

bool
EBSpeciesFluxOp::
checkNANFAB(const EBCellFAB&  a_data, const string& a_msg)
{ 
  Real vNAN = 1e399;
  bool dataIsNANINF = false;

  int nCons = Min(a_data.nComp(),5);  // to avoid checking all the primitive comps  
  Box box = a_data.getRegion();
  const EBISBox& ebisBox = a_data.getEBISBox();
  const EBGraph& ebgraph = ebisBox.getEBGraph();

  //if(~ebgraph.hasIrregular())  return dataIsNANINF;

  IntVectSet ivs(box);
  for (VoFIterator vofit(ivs, ebgraph);vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      IntVect gi = vof.gridIndex();
      bool dout = false;
      for (int k=0;k<nCons;k++) dout = dout || std::isnan(a_data(vof, k));// ! (abs(a_data(vof, k))<vNAN);
      if (dout){
	pout()  << a_msg   << ", " << gi[0] << ", "<< gi[1] << ", ";
	for (int k=0;k<nCons;k++) pout() << a_data(vof, k)<< ", ";	  
	pout() << endl;
	dataIsNANINF = true;}
    }
  if(dataIsNANINF) MayDay::Error(a_msg.c_str());
  return dataIsNANINF;
}
#include "NamespaceFooter.H"