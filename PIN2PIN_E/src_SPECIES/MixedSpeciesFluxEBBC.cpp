#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "BoxIterator.H"

#include "MixedSpeciesFluxEBBC.H"
#include "EBStencil.H"
#include "IonizationWaveEBDirBcFunc.H"
#include "NamespaceHeader.H"
/*****************/

/*****************/
MixedSpeciesFluxEBBC::MixedSpeciesFluxEBBC(const ProblemDomain& a_domain,
					   const EBISLayout&    a_layout,
					   const RealVect&      a_dx,
					   const IntVect*       a_ghostCellsPhi /*=0*/,
					   const IntVect*       a_ghostCellsRhs /*=0*/)
  :m_bc(a_domain, a_layout, a_dx,a_ghostCellsPhi,a_ghostCellsRhs)
{
  m_domain = a_domain;
  m_layout = a_layout;
  m_dx = a_dx;
  m_isDefined = false;
  m_dataBased = false;
  m_order=2;
  m_time = 0;

  //create a new ionization wave function (this is a bit hardwired)
  ParmParse pp;
  Real sphereDistance;
  pp.get("OSUsphere_distance",sphereDistance);
  Real xloc = sphereDistance/2.0;
  IonizationWaveEBDirBcFunc* EBBCfun = new IonizationWaveEBDirBcFunc();
  EBBCfun->define(1.0, 0, xloc, &m_time);
  m_func = RefCountedPtr<BaseBCValue>(EBBCfun);
  m_bc.setFunction(m_func);
  m_bc.setValue(0); //only Homogeneous Dirichlet BC for now
  m_bc.setOrder(m_order); 

}

void MixedSpeciesFluxEBBC::setOrder(int a_order)
{
  CH_assert(a_order >= 1 && a_order <= 2);

  if (m_order != a_order)
    {
      m_isDefined = false;
    }
  m_order = a_order;
}

//Called by EBSpeciesFluxOp::defineStencils() with args (*m_eblg.getCFIVS()),  1.0/m_dx;
void MixedSpeciesFluxEBBC::define(const LayoutData<IntVectSet>& a_cfivs,
                                  const Real&                   a_factor)
{

// Note:: No beta in here
  m_bc.define(a_cfivs, a_factor);
  LayoutData<BaseIVFAB<VoFStencil> >& poissSten = *(m_bc.getFluxStencil(0));

  const DisjointBoxLayout& dbl = m_layout.getDisjointLayout();
  LayoutData<VoFIterator > vofItIrreg;
  vofItIrreg.define(dbl); // vofiterator cache

  m_fluxStencil.define(dbl);
  m_fluxWeight.define(dbl);

  //make the Dirichlet stencils
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
    {
      const Box& curBox = dbl[dit()];
      const EBISBox& curEBISBox = m_layout[dit()];
      const EBGraph& curEBGraph = curEBISBox.getEBGraph();
      const IntVectSet& cfivsThisBox = a_cfivs[dit()];

      IntVectSet notRegular;
      int nComps = 1;

      notRegular |= curEBISBox.getIrregIVS  (curBox);
      notRegular |= curEBISBox.getMultiCells(curBox);

      vofItIrreg[dit()].define(notRegular,curEBISBox.getEBGraph());

      BaseIVFAB<VoFStencil>& curStencilBaseIVFAB = m_fluxStencil[dit()];
      BaseIVFAB<Real>&       curWeightBaseIVFAB  = m_fluxWeight[dit()];

      curStencilBaseIVFAB.define(notRegular,curEBGraph,nComps);
      curWeightBaseIVFAB.define(notRegular,curEBGraph,nComps);

      for (VoFIterator vofit(notRegular,curEBGraph); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();

          VoFStencil& curStencil = curStencilBaseIVFAB(vof,0);
          Real areaFrac = curEBISBox.bndryArea(vof);
          Real&       curWeight  = curWeightBaseIVFAB(vof,0);

	  	  
          bool isDirBC = isDirichlet(vof, curEBISBox);
	  if ((*m_data)[dit].getIVS().contains(vof.gridIndex()) && (!isDirBC) )
	    {
	      RealVect bndryCentroid = curEBISBox.bndryCentroid(vof);
	      bndryCentroid *= m_dx;
	      //int order = EBArith::getExtrapolationStencil(curStencil, bndryCentroid, m_dx, vof, curEBISBox);
	      int order = EBArith::getFirstOrderExtrapolationStencil(curStencil, bndryCentroid, m_dx, vof, curEBISBox,-1,NULL,0);
	      curWeight = 0e0;
	      // this fakeBeta makes the flux into the fluid
	      Real fakeBeta = -1;
	      Real specFac = (*m_data)[dit](vof, m_ispec)*fakeBeta;
	      //Need to magnify weight in fluxStencil with areafrac*factor, factor = 1/dx;
	      //Pass the factor in because m_dx[0] in here may not be same as in EBAMRPoissonOp
	      curStencil *= areaFrac*a_factor*specFac;
	    }
	  if(isDirBC)
	    {
	      //here use as flux -V*ne+D*\pop{ne}{n}
	      Real fakeBeta = -1;
	      Real specFac = areaFrac*a_factor * (*m_data)[dit](vof, m_ispec)*fakeBeta;
              curStencil = poissSten[dit()](vofit(), 0);
              Real bcoeFactor = (*m_Dk)[dit()](vofit(), 0);
	      curStencil *= bcoeFactor;
	      curStencil.add(vof,specFac);
	    }
	}
    }

  m_isDefined = true;
}
/*****************/
MixedSpeciesFluxEBBC::~MixedSpeciesFluxEBBC()
{
}
void MixedSpeciesFluxEBBC::applyEBFlux(EBCellFAB&                     a_lphi,
                                          const EBCellFAB&              a_phi,
                                          VoFIterator&                  a_vofit,
                                          const LayoutData<IntVectSet>& a_cfivs,
                                          const DataIndex&              a_dit,
                                          const RealVect&               a_probLo,
                                          const RealVect&               a_dx,
                                          const Real&                   a_factor,
                                          const bool&                   a_useHomogeneous,
                                          const Real&                   a_time)
{
  CH_TIME("MixedSpeciesFluxEBBC::applyEBFlux");
  CH_assert(a_lphi.nComp() == 1 );
  CH_assert(a_phi.nComp()  == 1);

  //NOTE :: a_factor is beta/dx

  const EBISBox&   ebisBox = a_phi.getEBISBox();
  for (a_vofit.reset(); a_vofit.ok(); ++a_vofit)
    {
      Real flux = 0.0;
      const VolIndex& vof = a_vofit();
      

      
      bool isDirBC = isDirichlet(vof, ebisBox);

      if(isDirBC)
	{
	  const BaseIVFAB<Real>& poissWeight = (m_bc.getFluxWeight())[a_dit];
	  Real poissWeightPt = poissWeight(vof, 0);
	  Real value = 0; // only homog for now
	  const Real& bcoef    = (*m_Dk)[a_dit](vof,0);
	  Real flux = poissWeightPt*value*bcoef;
	}
      else
	{
	  if (m_dataBased && m_wallBCInfo.size() > m_ispec)
	    {
	      if ( m_wallBCInfo[m_ispec] >=0 && (*m_data)[a_dit].getIVS().contains(vof.gridIndex()) )
		{
		  int nspec = m_wallBCInfo[m_ispec];
		  if(nspec < (*m_data).nComp()) flux = (*m_data)[a_dit](vof, nspec);
		}
	      //               
	    }
	  else if (m_bc.m_isFunction)
	    {
	      MayDay::Error("MixedSpeciesFluxEBBC:: not implemented");
	    }
	  else
	    {
	      if (m_bc.m_onlyHomogeneous)
		{
		  MayDay::Error("MixedSpeciesFluxEBBC::getFaceFlux called with undefined inhomogeneous BC");
		}
	      
	      flux = m_bc.m_value;
	    }
	}

      const Real& areaFrac = ebisBox.bndryArea(vof);
      flux *= areaFrac;
      // the real beta is in a_factor; this *fake* beta makes the flux into the fluid
      Real fakeBeta = -1;
      flux *= fakeBeta;
      a_lphi(vof,0) += flux * a_factor;
    }
}

/*****************/
void MixedSpeciesFluxEBBC::getEBFlux(Real&                         a_flux,
                                        const VolIndex&               a_vof,
                                        const LevelData<EBCellFAB>&   a_phi,
                                        const LayoutData<IntVectSet>& a_cfivs,
                                        const DataIndex&              a_dit,
                                        const RealVect&               a_probLo,
                                        const RealVect&               a_dx,
                                        const bool&                   a_useHomogeneous,
                                        const Real&                   a_time,
                                        const pair<int,Real>*         a_cacheHint )
{
  m_bc.getEBFlux(a_flux,
                 a_vof,
                 a_phi,
                 a_cfivs,
                 a_dit,
                 a_probLo,
                 a_dx,
                 a_useHomogeneous,
                 a_time,
                 a_cacheHint );
  Real bcoef = (*m_Dk)[a_dit](a_vof,0);
  a_flux *= bcoef;

}
/*****************/
void MixedSpeciesFluxEBBC::setValue(Real a_value)
{
  m_bc.setValue(a_value);
}
/*****************/
bool MixedSpeciesFluxEBBC::isDirichlet(const VolIndex& a_vof, const EBISBox& a_EBISBox) const
{
  RealVect centroid = a_EBISBox.bndryCentroid(a_vof);
  centroid *= m_dx;
  RealVect vectDx = m_dx*RealVect::Unit;
  RealVect vofloc = EBArith::getVofLocation(a_vof, vectDx, centroid);
  //this should apply only to electrons
  if(m_wallBCInfo[m_ispec] >=0 && vofloc[1] > -m_dx[0]/2)
    {
      return m_func->value(vofloc, RealVect::Unit, m_time, m_ispec) > 0.5;
    }
  else return false;
}
/*****************/
void MixedSpeciesFluxEBBC::setFunction(RefCountedPtr<BaseBCValue> a_func)
{
  m_bc.setFunction(a_func);
}

/*****************/
MixedSpeciesFluxEBBCFactory::MixedSpeciesFluxEBBCFactory()
{
  m_value = 12345.6789;
  m_func = RefCountedPtr<BaseBCValue>();
  m_dataBased = false;
  m_onlyHomogeneous = true;
  m_isFunction = false;
}
/*****************/
MixedSpeciesFluxEBBCFactory::~MixedSpeciesFluxEBBCFactory()
{
}
/*****************/
void MixedSpeciesFluxEBBCFactory::setValue(Real a_value)
{
  m_value = a_value;
  m_func = RefCountedPtr<BaseBCValue>();

  m_onlyHomogeneous = false;
  m_isFunction = false;
}
/*****************/
void MixedSpeciesFluxEBBCFactory::setFunction(RefCountedPtr<BaseBCValue> a_func)
{
  m_value = 12345.6789;
  m_func = a_func;

  m_onlyHomogeneous = false;
  m_isFunction = true;
}
/*****************/
MixedSpeciesFluxEBBC* MixedSpeciesFluxEBBCFactory::create(const ProblemDomain& a_domain,
                                                                const EBISLayout&    a_layout,
                                                                const RealVect&      a_dx,
                                                                const IntVect*       a_ghostCellsPhi /*=0*/,
                                                                const IntVect*       a_ghostCellsRhs /*=0*/)
{
  MixedSpeciesFluxEBBC* fresh = new MixedSpeciesFluxEBBC(a_domain,a_layout,a_dx,a_ghostCellsPhi,a_ghostCellsRhs);


  if (!m_onlyHomogeneous)
    {
      if (m_dataBased)
        {
          fresh->setData(m_data, m_ispec, m_wallBCInfo);
        }
      else if (!m_isFunction)
        {
          fresh->setValue(m_value);
        }
      else
        {
          fresh->setFunction(m_func);
        }
    }

  return fresh;
}
#include "NamespaceFooter.H"
