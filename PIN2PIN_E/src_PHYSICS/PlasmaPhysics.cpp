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
#include "PlasmaPhysics.H"
#include "PlasmaPhysicsF_F.H"
#include "BaseEBCellFAB.H"
#include "EBCellFAB.H"
#include "FArrayBox.H"
#include "LoHiSide.H"
#include "FaceIterator.H"
#include "VoFIterator.H"
#include "BoxIterator.H"
#include "CH_Timer.H"
#include "EBArith.H"
#include "EBIndexSpace.H"
#include "EBLevelDataOps.H"
#include "CH_assert.H"
#include "Stencils.H"
#include "LevelData.H"
#include "EBISLayout.H"
#include "EBPatchPolytropic.H"
#include "ParmParse.H"
#include "Helium.H"
#include "ConstCHIMU.H"

#include "NamespaceHeader.H"

/***/

//Global variables
int PlasmaPhysics::cvode_NEQ;
plasmaParameters PlasmaPhysics::m_ppars;
Vector<Real> PlasmaPhysics::m_minData;
Vector<Real> PlasmaPhysics::m_maxData;
string PlasmaPhysics:: m_bgSpec;
Vector<Vector<Real> > PlasmaPhysics::m_rrates;
Vector<Vector<Real> > PlasmaPhysics::m_transp;
Vector<Vector<Real> > PlasmaPhysics::m_reactant;
Vector<Vector<Real> > PlasmaPhysics::m_netReac;
Vector<Real> PlasmaPhysics::m_actReac;
Vector<Real> PlasmaPhysics::m_esecemcoef;
Vector<Real> PlasmaPhysics::m_esecemen;
Vector<int> PlasmaPhysics::m_species;
Vector<int> PlasmaPhysics::m_reactions;
int PlasmaPhysics::m_tgComp;
int PlasmaPhysics::m_eeComp;
int PlasmaPhysics::m_bgindex;
int PlasmaPhysics::m_o2index;
int PlasmaPhysics::m_eindex;
Real PlasmaPhysics::m_GasTemperature;
Real PlasmaPhysics::m_GasPressure;
Real PlasmaPhysics::m_odeSource;
Real PlasmaPhysics::m_odeEfield;
Real PlasmaPhysics::m_nBg;
Real PlasmaPhysics::m_Te;
Real PlasmaPhysics::m_Temax;
Real PlasmaPhysics::m_Tg;
bool PlasmaPhysics::m_print = false;
bool PlasmaPhysics::m_DC = true;
int PlasmaPhysics::m_constTe;
int PlasmaPhysics::m_constBg;
VolIndex PlasmaPhysics::m_vof;
Real PlasmaPhysics::m_tfac;

/***/
// Contructor
PlasmaPhysics::PlasmaPhysics()
{
  m_isDefined = false;
  {
    //Helium* myGas = new Helium();
    ConstCHIMU* myGas = new ConstCHIMU();
    m_gasModel = RefCountedPtr<gasModel>(myGas);
  }
  m_gasModel->setPlasmaParameters(m_ppars, m_minData, m_maxData,  m_bgSpec, m_charSpec);
  
  m_constTe = m_ppars.isConstTe;
  m_constBg = m_ppars.isConstBg;
	
}
void
PlasmaPhysics::define()
{ 
  //getPlasmaParameters();
  //bg stands for background gas (e.g. N2 or He, here He)
  m_bgindex = findSpec(m_bgSpec);  
  m_o2index = -100;
  m_eindex = findSpec(string("e-"));
  m_tgComp = 0;
  m_eeComp = 1;

  
  m_species = m_ppars.species; //This excludes the EE
  m_reactions = m_ppars.reactions;
  m_rrates = m_ppars.reacRateCoef;
  m_reactant = m_ppars.reactant;
  m_transp = m_ppars.transpCoef;
  m_netReac = m_ppars.netReac;
  m_actReac = m_ppars.actReac;
  m_isDefined = true;
  Real Torr2Pa = 133.3224;
  /***********/
  // enforcing constant Te and/or background gas concentration
  m_GasPressure = Torr2Pa*300; //Roy POP06 300 torr
  ParmParse pp;
  pp.query("pinf", m_GasPressure);
  // enforcing constant Te and/or background gas concentration  
  pp.query("constTe", m_constTe);
  pp.query("constBg", m_constBg);
  m_DC = false;//for wall BC
  pp.query("DC", m_DC); //for wall BC
  // if frequency <= 0 set it to DC
  if(pp.contains("frequency"))
    {
      Real frequency;pp.get("frequency", frequency); 
      if(frequency<=0)m_DC = true;
    }
  Real eVT = 11604.505;//1 eV-temp (hardwired move it to the gas model)
  m_Te = eVT;
  m_Temax = 15.0* eVT; //not passed to Fortran, update temax there if changed here
  // for IC and BC and source term
  m_Tg = 300.0; //[K]
  pp.query("tinf", m_Tg);
  m_molWeight = 4.0; //kg/kMole of the background gas
  m_useThermalVelocity = (m_constTe==0);
  //convert EE -> Te
  m_tfac =  2.0/(3.0*m_ppars.kb*m_ppars.NA);
  /***********/

  Real bgfrac = 1.0; // Background gas
  m_nBg = m_GasPressure*bgfrac/(m_ppars.kb*m_Tg)/m_ppars.NA; //[mol.m^-3]
  m_nO2 = m_GasPressure*(1.0-bgfrac)/(m_ppars.kb*m_Tg)/m_ppars.NA;
  m_maxData[0] = 2*m_nBg;

  //Wall BC info
  setWallBCinfo();

  //cvode stuff
  m_odeflag = true;
  cvode_NEQ = nComponents();
  cvode_mem = NULL;
  cvode_y=cvode_abstol=NULL;
  cvode_A = NULL;
  cvode_LS = NULL;
  cvode_y = N_VNew_Serial(cvode_NEQ);
  cvode_abstol = N_VNew_Serial(cvode_NEQ);
  /* Create dense SUNMatrix for use in linear solves */
  cvode_A = SUNDenseMatrix(cvode_NEQ, cvode_NEQ);
  /* Create dense SUNLinearSolver object for use by CVode */
  pout() << "cvode_LS allocation ";
  cvode_LS = SUNDenseLinearSolver(cvode_y, cvode_A);
  

  Real ATOL=1e-13;
  if(m_ppars.isDimensionless) ATOL=1e-6;
  for(int k=1;k<=cvode_NEQ;k++)
    {
      Ith(cvode_abstol,k) = ATOL;
      Ith(cvode_y,k) = ATOL;
    }
  Ith(cvode_y,cvode_NEQ) = 1e-3;
  NV_Ith_S(cvode_abstol,m_bgindex) = 1e-3;
  //MDCancel for air: NV_Ith_S(cvode_abstol,m_o2index) = 1e-3;
 
  cvode_reltol = 1e-6;
  realtype T0 = 0e0;

  
  /* Call CVodeCreate to create the solver memory and specify the 
   * Backward Differentiation Formula and the use of a Newton iteration */
  if(m_ppars.isStiff)
    cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  else
    cvode_mem = CVodeCreate(CV_ADAMS, CV_FUNCTIONAL);
  
  
  int flag = CVodeInit(cvode_mem, cvode_f, T0, cvode_y);
  flag = CVodeSVtolerances(cvode_mem, cvode_reltol, cvode_abstol);

  if(m_ppars.isStiff){
    /* Call CVDlsSetLinearSolver to attach the matrix and linear solver to CVode */
    flag = CVDlsSetLinearSolver(cvode_mem, cvode_LS, cvode_A);
    flag = CVDlsSetJacFn(cvode_mem, NULL);
  }
  
  //flag = CVodeSetMaxOrd(cvode_mem, 2);
  flag = CVodeSetMaxNumSteps(cvode_mem, 10000);
  //LM locating the top electrode and its potential


  //LM Refinement Thresholds
  m_threshSet = false;
  
}
/***/
PlasmaPhysics::
~PlasmaPhysics()
{
}

void PlasmaPhysics::setEBy(const Real& a_yEB)
{
  m_yEB = a_yEB;
}

void PlasmaPhysics::setEps(	const Real& a_epsg,
				const Real& a_epsd)
{
  
  if(m_ppars.isDimensionless)
    {
      m_epsg = m_ppars.q*m_ppars.NA*abs(m_ppars.chargeNr[m_eindex]);
    }
  else
    {
      m_epsg = a_epsg;
      m_epsd = a_epsd;
    }
}

int PlasmaPhysics::nReactions()
{
  return m_reactions.size();
}

int PlasmaPhysics::nSpecies()
{
  return m_species.size();
}

plasmaParameters PlasmaPhysics::plasmaPars()
{
  return m_ppars;
}

/*******************
Initialize Flow
******/
void PlasmaPhysics::initialize(LevelData<EBCellFAB>& a_conState,
			       const EBISLayout&     a_ebisl,
			       const Real&           a_dx,
			       const RealVect&       a_domainLength)
{
  // concentration in [m^-3]
  Real ne = 1.0e4/m_ppars.NA;
  Real electronEnergy;
  if (m_constTe == 1)
    electronEnergy = 3.0/2.0*ne*m_ppars.NA*m_ppars.kb*m_Te;
  else
    electronEnergy = 3.0/2.0*ne*m_ppars.NA*m_ppars.kb*300.0;
  m_minData[2+1]= electronEnergy/1e5;



  EBLevelDataOps::setToZero(a_conState);
  EBLevelDataOps::setVal(a_conState, m_nBg, m_bgindex);
  //MDCancel for air: EBLevelDataOps::setVal(a_conState, m_nO2, m_o2index);
  EBLevelDataOps::setVal(a_conState, m_ppars.initNe, m_eindex);
  EBLevelDataOps::setVal(a_conState,electronEnergy,m_eindex+1);

  // run the gas specific intialization
  m_gasModel->initialize(a_conState, m_ppars, a_ebisl, a_dx, a_domainLength);
  
}
  

/*******************
Box BCs
******/
Real PlasmaPhysics::BCvalue(int&           a_ispec,
			    const int&     a_idir,
			    const Side::LoHiSide& a_side)							
{
  Real ne = 0;
  Real electronEnergy = 3.0/2.0*ne*m_ppars.kb*m_ppars.NA*300.0;

  if(m_ppars.isConstMobDiff)
    {
      if (m_ppars.isNeutral[a_ispec]) return m_nBg;
      if(a_ispec == m_eindex) return (a_side == Side::Lo)?  m_ppars.initNe : 0.0;
      if(a_ispec == m_ppars.iIndex) return (a_side == Side::Lo)?  m_ppars.initNe : 0.0;
      return 0;
    }
  else
    {
      if(a_ispec == m_bgindex)
	return m_nBg;
      else if (a_ispec == m_o2index)
	return m_nO2;
      else if (a_ispec == m_eindex)
	return ne;
      else if (a_ispec == m_eindex+1)
	return electronEnergy;
      else
	return 0;
    }
}

/*************************
 Plasma Source Terms
******/

void PlasmaPhysics::plasmaSources(LevelData<EBCellFAB>&   	a_src,
				  const LevelData<EBCellFAB>& 	a_NData,
				  const LevelData<EBCellFAB>& 	a_teData,
				  const LevelData<EBCellFAB>& 	a_GphiData,
				  const LevelData<EBCellFAB>& 	a_fGphiData,
				  const EBISLayout&	        a_ebisl)

{
  // a_src contains all species + last component for electron energy source

  CH_assert(a_src.nComp() >= nSpecies()+1+3);//last 3 comp. for separate Joule, Elastic and Inelastic terms (plotting)
  EBCellFactory factory(a_ebisl);
  LevelData<EBCellFAB> lsen(a_src.disjointBoxLayout(), 4 , a_src.ghostVect(), factory);
  LevelData<EBCellFAB> lrr(a_src.disjointBoxLayout(), nReactions() , a_src.ghostVect(), factory);

  for (DataIterator dit = a_src.dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB& src = a_src[dit()];
      EBCellFAB& sen = lsen[dit()];
      EBCellFAB& rr = lrr[dit()];

      const EBCellFAB& NData = a_NData[dit()];
      const EBCellFAB& teData = a_teData[dit()];
      const EBCellFAB& fGphiData = a_fGphiData[dit()];
      const EBCellFAB& GphiData = a_GphiData[dit()];

      reactionRates(rr, NData, teData, GphiData);
      plasmaSourceEnergy(sen, rr, NData, teData, fGphiData);
      plasmaSourceTerms(src, rr, sen);
    }
}
 
//-------------------
     
void PlasmaPhysics::plasmaSourceTerms(EBCellFAB& 	a_src,
				      const EBCellFAB&	a_rr,
				      const EBCellFAB&	a_sen)

{
  const EBISBox&   ebis   = a_src.getEBISBox();  
  const EBISBox&   rrEbis   = a_rr.getEBISBox(); 
  const EBISBox&   senEbis   = a_sen.getEBISBox(); 
  //MDProb CH_assert(ebis == rrEbis);
		
  Box locRegion = a_src.getRegion() & a_rr.getRegion() & a_sen.getRegion();

  if (!locRegion.isEmpty())
    {
      BaseFab<Real>& srcRegFAB = a_src.getSingleValuedFAB();
      const BaseFab<Real>& rrRegFAB = a_rr.getSingleValuedFAB();
      const BaseFab<Real>& senRegFAB = a_sen.getSingleValuedFAB();

      //bool sameRegBox = ((srcRegFAB.box() == rrRegFAB.box()) && (srcRegFAB.box() == senRegFAB.box()));
      //bool nvofTest   = false; no Vof test here
      //bool regionTest = false;
      //MDnote: add vof check and pointer calculation later
      BaseIVFAB<Real>& srcIrrBFAB = a_src.getMultiValuedFAB();
      const BaseIVFAB<Real>& rrIrrBFAB = a_rr.getMultiValuedFAB();
      const BaseIVFAB<Real>& senIrrBFAB = a_sen.getMultiValuedFAB();
      //if (locRegion == a_src.getRegion())
      //	regionTest = true;
      for (int s = 0; s < nSpecies(); s++)
	{
	  if ((s != m_bgindex) ||  ((s = m_bgindex) && (m_constBg == 0)))
	    {
	      Vector<Real> par = m_netReac[m_species[s]];
			   				
	      FORT_PSOURCEFAB(CHF_FRA(srcRegFAB),
			      CHF_CONST_FRA(rrRegFAB),
			      CHF_CONST_FRA(senRegFAB),
			      CHF_CONST_VR(par),
			      CHF_CONST_INT(s),
			      CHF_BOX(locRegion));
     	
	      // Multi-Cells	
	      IntVectSet ivsMulti = a_rr.getMultiCells();
	      ivsMulti &= a_sen.getMultiCells();
	      ivsMulti &= a_src.getMultiCells();  // opt might be to take this out
	      ivsMulti &= locRegion;
	      IVSIterator ivsit(ivsMulti);			
	      for (ivsit.reset(); ivsit.ok(); ++ivsit)
		{
		  const IntVect& iv = ivsit();
		  Vector<VolIndex> vofs = ebis.getVoFs(iv);
		  for (int ivof = 0; ivof < vofs.size(); ivof++)
		    {
		      const VolIndex& vof = vofs[ivof];
		      if (locRegion.contains(vof.gridIndex()))
			{   
			  srcIrrBFAB(vof, s) = 0;
			  for (int r = 0; r < nReactions(); r++)
			    {													
			      srcIrrBFAB(vof, s) += rrIrrBFAB(vof, r) * par[r];
			    }
			  if (s == m_eindex)
			    {
			      for (int i = 0; i < 4; i++)
				srcIrrBFAB(vof, nSpecies()+i) = senIrrBFAB(vof, i);
			    }
			}			 
		    }
		}
	    }
	}
    }
}

//-----------------
// Note: JData includes the electron Joule heating calculated in a prior call of AddJouleHeat 
// Only non-zero if m_constTe=0 
void PlasmaPhysics::plasmaSourceEnergy(EBCellFAB&			a_sen,
				       const EBCellFAB&	a_rr,
				       const EBCellFAB& 	a_NData,
				       const EBCellFAB&	a_teData,
				       const EBCellFAB& 	a_JData)
						
{
  const EBISBox&   ebis   = a_sen.getEBISBox();  
  const EBISBox&   rrDataEbis   = a_rr.getEBISBox(); 
  const EBISBox&   NDataEbis   = a_NData.getEBISBox(); 
  const EBISBox&   teDataEbis   = a_teData.getEBISBox(); 
  const EBISBox&   JDataEbis   = a_JData.getEBISBox(); 
   	
  Box locRegion = a_sen.getRegion() & a_rr.getRegion() & a_teData.getRegion() & a_NData.getRegion() & a_JData.getRegion();
  
  a_sen.setVal(0e0);
  if (!locRegion.isEmpty() && (m_constTe == 0))
    //if (!locRegion.isEmpty())
    // Variable Te, evaluate electron source term)
    {   
      BaseFab<Real>& senRegFAB = a_sen.getSingleValuedFAB();
      const BaseFab<Real>& rrRegFAB = a_rr.getSingleValuedFAB();
      const BaseFab<Real>& NRegFAB = a_NData.getSingleValuedFAB();
      const BaseFab<Real>& teRegFAB = a_teData.getSingleValuedFAB();
      const BaseFab<Real>& JRegFAB = a_JData.getSingleValuedFAB();
   
      BaseIVFAB<Real>& senIrrBFAB = a_sen.getMultiValuedFAB();
      const BaseIVFAB<Real>& rrIrrBFAB = a_rr.getMultiValuedFAB();
      const BaseIVFAB<Real>& NIrrBFAB = a_NData.getMultiValuedFAB();
      const BaseIVFAB<Real>& teIrrBFAB = a_teData.getMultiValuedFAB();
      const BaseIVFAB<Real>& JIrrBFAB = a_JData.getMultiValuedFAB();
	   
      //if (locRegion == a_sen.getRegion())
      //regionTest = true;
      Real buffel, buffin, tel, teg, exsbg, exso2, een, lte, ve, Ng; 
      Vector<Real> epar = m_ppars.eInelEn;
      Vector<Real> mpar = m_ppars.mass;
      Vector<Real> sigebg = m_ppars.sigebg;  
      Vector<Real> sigeo2 = m_ppars.sigeo2;   
      Real etaT = m_ppars.etaT;
 
      FORT_PESOURCEFAB(CHF_FRA(senRegFAB),
		       CHF_CONST_FRA(rrRegFAB),
		       CHF_CONST_FRA(NRegFAB),
		       CHF_CONST_FRA(teRegFAB),
		       CHF_CONST_FRA(JRegFAB),
		       CHF_CONST_VR(epar),
		       CHF_CONST_VR(mpar),
		       CHF_CONST_VR(sigebg),  
		       CHF_CONST_VR(sigeo2), 
		       CHF_CONST_VI(m_species),
		       CHF_CONST_REAL(m_nBg),
		       CHF_CONST_REAL(m_Te),
		       CHF_CONST_INT(m_tgComp),
		       CHF_CONST_INT(m_eindex),
		       CHF_CONST_INT(m_bgindex),
		       CHF_CONST_INT(m_o2index),  
		       CHF_CONST_INT(m_constTe),
		       CHF_CONST_INT(m_constBg),
		       CHF_BOX(locRegion));
       
      //Multi-cells
      IntVectSet ivsMulti = a_sen.getMultiCells();
      ivsMulti &= a_rr.getMultiCells();
      // opt might be to take this out
      ivsMulti &= a_NData.getMultiCells(); 
      ivsMulti &= a_teData.getMultiCells(); 
      ivsMulti &= a_JData.getMultiCells(); 
      ivsMulti &= locRegion;
	
      IVSIterator ivsit(ivsMulti);			
      for (ivsit.reset(); ivsit.ok(); ++ivsit)
	{
	  const IntVect& iv = ivsit();
	  Vector<VolIndex> vofs = ebis.getVoFs(iv);
	  for (int ivof = 0; ivof < vofs.size(); ivof++)
	    {
	      const VolIndex& vof = vofs[ivof];
	      if (locRegion.contains(vof.gridIndex()))
		{   
		  //if (m_constTe == 0)	
		  //	{
		  een = teIrrBFAB(vof, m_eeComp)/(NIrrBFAB(vof, m_eindex)*m_ppars.NA);
		  tel=convertEe2Te(teIrrBFAB(vof, m_eeComp), NIrrBFAB(vof, m_eindex));
		  /*
		    }
		    else
		    {
		    tel = m_Te;
		    een = convertTe2Een(tel);		  
		    }
		  */
		  teg = teIrrBFAB(vof, m_tgComp);  
		
		  if (m_constBg == 0)
		    Ng = NIrrBFAB(vof, m_bgindex);
		  else
		    Ng = m_nBg;
			 
		  //Inelastic collision losses		 
		  buffin = 0;		
		  for (int r = 0; r < nReactions(); r++)												
		    {
		      if (een/m_ppars.q > abs(epar[r]) )
			buffin -= rrIrrBFAB(vof, r) * epar[r];  
		    }
		  buffin *= m_ppars.q * m_ppars.NA;  
		  //Elastic collision losses
		  /*	
			exsbg = 1.0;  
			exso2 = 1.0;  
			for (int l=0; l<3; l++)
			{
			exsbg *= exp(sigebg[l*3]*exp(-pow((log(een)-sigebg[3*l+1])/sigebg[3*l+2],2)));  
			exso2 *= exp(sigeo2[l*3]*exp(-pow((log(een)-sigeo2[3*l+1])/sigeo2[3*l+2],2)));  
			}
		  */
		  // MDCancel Helium only
		  lte = log10(tel);
		  exsbg = sigebg[0];
		  for (int j = 1; j < 8; j++)	
		    exsbg = exsbg + sigebg[j] * pow(lte,j);		
		  exsbg = pow(10.0,exsbg);		
		  ve = sqrt(2.0 * m_ppars.kb * m_ppars.NA * abs(tel) / mpar[m_species[m_eindex]]);
		  buffel = 3.0 * m_ppars.kb * pow(m_ppars.NA, 2.0) * NIrrBFAB(vof, m_eindex) * mpar[m_species[m_eindex]] * (tel-teg) * ve 
		    * ( Ng * exsbg / mpar[m_species[m_bgindex]]);  
		  //MDCancel for air: * ( NIrrBFAB(vof, m_bgindex) * exsbg / mpar[m_species[m_bgindex]]) + NIrrBFAB(vof, m_o2index) * exso2 / mpar[m_species[m_o2index]] );  
		  if (m_constTe == 0)
		    //Joule is already included as is in the input array JIrrBFAB
		    senIrrBFAB(vof, 0) = JIrrBFAB(vof, m_eindex) - buffel - buffin ; //total e- energy source term
		  // For plotting:
		  senIrrBFAB(vof, 1) = JIrrBFAB(vof, m_eindex); //Joule energy gain
		  senIrrBFAB(vof, 2) = buffel; //Elastic loss	
		  senIrrBFAB(vof, 3) = buffin; //Inelastic loss 			    							
		}			 
	    }
	}	
    }
}

/****************************************
 Fluid Source Terms (non-incl. Joule)
*************************/

void PlasmaPhysics::fluidSources(LevelData<EBCellFAB>&   	a_src,
				 const LevelData<EBCellFAB>& 	a_NData,
				 const LevelData<EBCellFAB>& 	a_teData,
				 const LevelData<EBCellFAB>& 	a_GphiData,
				 const LevelData<EBCellFAB>& 	a_Velo,
				 const EBISLayout&		a_ebisl)

{
  // a_src contains fx and fy body forces as components [0->SpaceDim-1] and gas energy source term as last component [SpaceDim]
  EBCellFactory factory(a_ebisl);
  LevelData<EBCellFAB> lrr(a_src.disjointBoxLayout(), nReactions() , a_src.ghostVect(), factory);

  bool isNAN = EBLevelDataOps::checkNANINF(a_GphiData);
  if(isNAN) MayDay::Error("Checking GphiData.");
  //pout() << "ByeBye";exit(0);

  for (DataIterator dit = a_src.dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB& src = a_src[dit()];
      EBCellFAB& rr = lrr[dit()];

      const EBCellFAB& NData = a_NData[dit()];
      const EBCellFAB& teData = a_teData[dit()];
      const EBCellFAB& GphiData = a_GphiData[dit()];
      const EBCellFAB& Velo = a_Velo[dit()];

      reactionRates(rr, NData, teData, GphiData);
      fluidSourceTerms(src, rr, NData, teData, GphiData, Velo);
    }
}

//-----------------
// Note: fluidSourceTerms does not include Joule heating, which is added in a subsequent call of AddJouleHeatF
//
  
void PlasmaPhysics::fluidSourceTerms(EBCellFAB&		a_src,
				     const EBCellFAB&	a_rr,
				     const EBCellFAB& 	a_NData,
				     const EBCellFAB&	a_teData,
				     const EBCellFAB& 	a_GphiData,
				     const EBCellFAB& 	a_Velo)
						
{
  const EBISBox&   ebis   = a_src.getEBISBox();  
  const EBISBox&   rrDataEbis   = a_rr.getEBISBox(); 
  const EBISBox&   NDataEbis   = a_NData.getEBISBox(); 
  const EBISBox&   teDataEbis   = a_teData.getEBISBox(); 
  const EBISBox&   GphiDataEbis   = a_GphiData.getEBISBox(); 

  EBPatchPolytropic Ppatch;
  Interval momIntrv = Ppatch.momentumInterval();
  int enrgIndx = Ppatch.energyIndexC();
   	
  Box locRegion = a_src.getRegion() & a_rr.getRegion() & a_teData.getRegion() & a_NData.getRegion() ;
 
  if (!locRegion.isEmpty())
    {
      a_src.setVal(0e0);
	 
      BaseFab<Real>& srcRegFAB = a_src.getSingleValuedFAB();
      const BaseFab<Real>& rrRegFAB = a_rr.getSingleValuedFAB();
      const BaseFab<Real>& NRegFAB = a_NData.getSingleValuedFAB();
      const BaseFab<Real>& teRegFAB = a_teData.getSingleValuedFAB();
      const BaseFab<Real>& GphiRegFAB = a_GphiData.getSingleValuedFAB();
      const BaseFab<Real>& VeloRegFAB = a_Velo.getSingleValuedFAB();
   
      BaseIVFAB<Real>& srcIrrBFAB = a_src.getMultiValuedFAB();
      const BaseIVFAB<Real>& rrIrrBFAB = a_rr.getMultiValuedFAB();
      const BaseIVFAB<Real>& NIrrBFAB = a_NData.getMultiValuedFAB();
      const BaseIVFAB<Real>& teIrrBFAB = a_teData.getMultiValuedFAB();
      const BaseIVFAB<Real>& GphiIrrBFAB = a_GphiData.getMultiValuedFAB();
      const BaseIVFAB<Real>& VeloIrrBFAB = a_Velo.getMultiValuedFAB();
	   
      //if (locRegion == a_src.getRegion())
      //	regionTest = true;
      int jj;
      Real buff, tel, teg, exsbg, exso2, een, charge, lte, ve, Ng;  
      Vector<Real> gpar = m_ppars.gInelEn;
      Vector<Real> zpar = m_ppars.chargeNr;
      Vector<Real> mpar = m_ppars.mass;
      Vector<Real> sigebg = m_ppars.sigebg;  
      Vector<Real> sigeo2 = m_ppars.sigeo2;  
      Real etaT = m_ppars.etaT;
      int mombeg = momIntrv.begin();
 
      FORT_FSOURCEFAB(CHF_FRA(srcRegFAB),
		      CHF_CONST_FRA(rrRegFAB),
		      CHF_CONST_FRA(NRegFAB),
		      CHF_CONST_FRA(teRegFAB),
		      CHF_CONST_FRA(GphiRegFAB),
		      CHF_CONST_FRA(VeloRegFAB),
		      CHF_CONST_VR(gpar),
		      CHF_CONST_VR(zpar),
		      CHF_CONST_VR(mpar),
		      CHF_CONST_VR(sigebg),  
		      CHF_CONST_VR(sigeo2),  
		      CHF_CONST_VI(m_species),
		      CHF_CONST_REAL(m_nBg),
		      CHF_CONST_REAL(m_Te),
		      CHF_CONST_REAL(etaT),
		      CHF_CONST_INT(m_tgComp),
		      CHF_CONST_INT(m_eindex),
		      CHF_CONST_INT(m_bgindex),
		      CHF_CONST_INT(m_o2index),  
		      CHF_CONST_INT(m_constBg), 
		      CHF_CONST_INT(m_constTe), 
		      CHF_CONST_INT(mombeg),  
		      CHF_CONST_INT(enrgIndx),  
		      CHF_BOX(locRegion));
	
      //Multi-cells
      IntVectSet ivsMulti = a_src.getMultiCells();
      ivsMulti &= a_rr.getMultiCells();
      // opt might be to take this out
      ivsMulti &= a_NData.getMultiCells(); 
      ivsMulti &= a_teData.getMultiCells(); 
      ivsMulti &= locRegion;
	
      IVSIterator ivsit(ivsMulti);			
      for (ivsit.reset(); ivsit.ok(); ++ivsit)
	{
	  const IntVect& iv = ivsit();
	  Vector<VolIndex> vofs = ebis.getVoFs(iv);
	  for (int ivof = 0; ivof < vofs.size(); ivof++)
	    {
	      const VolIndex& vof = vofs[ivof];
	      if (locRegion.contains(vof.gridIndex()))
		{   			

		  //MDCancel 
		  if (m_constTe == 0)
		    {
		      for (int r = 0; r < nReactions(); r++)
			{			
			  //MD 01/13/13 minus sign (* negative gpar components)										
			  srcIrrBFAB(vof, enrgIndx) -= rrIrrBFAB(vof, r) * gpar[r];
			}
		      srcIrrBFAB(vof, enrgIndx) *= m_ppars.q * m_ppars.NA;  
		      if (m_constTe == 0)	
			{
			  een = teIrrBFAB(vof, m_eeComp)/(NIrrBFAB(vof, m_eindex)*m_ppars.NA);
			  tel = convertEe2Te(teIrrBFAB(vof, m_eeComp), NIrrBFAB(vof, m_eindex));
			}
		      else
			{
			  tel = m_Te;
			  een = convertTe2Een(tel);		  
			}
		      teg = teIrrBFAB(vof, m_tgComp);  
		
		      if (m_constBg == 0)
			Ng = NIrrBFAB(vof, m_bgindex);
		      else
			Ng = m_nBg;
		      /*
			exsbg = 1.0;  
			exso2 = 1.0;  
			for (int l=0; l<3; l++)
			{
			exsbg *= exp(sigebg[l*3]*exp(-pow((log(een)-sigebg[3*l+1])/sigebg[3*l+2],2)));  
			exso2 *= exp(sigeo2[l*3]*exp(-pow((log(een)-sigeo2[3*l+1])/sigeo2[3*l+2],2)));  
			}
		      */
		      lte = log10(tel);
		      exsbg = sigebg[0];
		      for (int j = 1; j < 8; j++)
			exsbg = exsbg + sigebg[j] * pow(lte,j);
		      exsbg = pow(10.0,exsbg);
		
		      ve = sqrt(2.0 * m_ppars.kb * m_ppars.NA * abs(tel) / mpar[m_species[m_eindex]]);		
		      buff = 3.0 * m_ppars.kb * pow(m_ppars.NA, 2.0) * NIrrBFAB(vof, m_eindex) * mpar[m_species[m_eindex]] * (tel-teg) * ve 
			* ( Ng * exsbg / mpar[m_species[m_bgindex]]);  
		      //MDCancel for air: * ( NIrrBFAB(vof, m_bgindex) * exsbg / mpar[m_species[m_bgindex]]) + NIrrBFAB(vof, m_o2index) * exso2 / mpar[m_species[m_o2index]] );  
		
		      srcIrrBFAB(vof, enrgIndx) += buff;	
		    }	//end if m_constTe
		 							
		  // Momentum source and Advective loss Q*E*V
		  charge = gasCharge(NIrrBFAB, vof);

		  for (int j = 0; j < SpaceDim; j++)
		    //MD note: mombeg in 2D currently=1 ==> jj=1,2 and enrgIndx=3 (in output array for plot, look for indexes 1->3 not 0->2)
		    {		    
		      jj = mombeg + j;
		      srcIrrBFAB(vof, jj) = -GphiIrrBFAB(vof, j) * charge;
		      if (m_constTe == 0)
			srcIrrBFAB(vof, enrgIndx) += srcIrrBFAB(vof, jj) * VeloIrrBFAB(vof, j);
		    }		   			 
		}
	    }	
	}
    }
}
/*******************
 Reaction Rates
******/

void PlasmaPhysics::reactionRates(LevelData<EBCellFAB>&   	a_rr,
				  const LevelData<EBCellFAB>& 	a_NData,
				  const LevelData<EBCellFAB>& 	a_teData,
				  const LevelData<EBCellFAB>& 	a_GphiData)

{
  CH_assert(a_rr.nComp() >= nReactions());
  CH_assert(a_NData.nComp() >= nSpecies()); 
  //CH_assert(a_fGphiData.nComp() >= nSpecies()); 
  for (DataIterator dit = a_rr.dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB& rr = a_rr[dit()];
      const EBCellFAB& teData = a_teData[dit()];
      const EBCellFAB& NData = a_NData[dit()];
      const EBCellFAB& GphiData =	a_GphiData[dit()];
      //const EBCellFAB& fGphiData =	a_fGphiData[dit()];
      //CH_assert(rr.isDefined() && teData.isDefined());
      reactionRates(rr, NData, teData, GphiData);
    }
}
      
void PlasmaPhysics::reactionRates(EBCellFAB&            a_rr, 
				  const EBCellFAB& 	a_NData, 
				  const EBCellFAB& 	a_teData,
				  const EBCellFAB& 	a_GphiData)

{
  const EBISBox&   ebis   = a_rr.getEBISBox(); 
  const EBISBox&   teDataEbis   = a_teData.getEBISBox();
  const EBISBox&   NDataEbis   = a_NData.getEBISBox();  
  const EBISBox&   GphiDataEbis   = a_GphiData.getEBISBox();
  //const EBISBox&   fGphiDataEbis   = a_fGphiData.getEBISBox();

  
  Real Torr2Pa = 133.3224;
  Real cA=4.4e2/Torr2Pa;
  Real cB=14.0*pow(1e2/Torr2Pa,0.4);
  int eIndex=m_ppars.eIndex;
  int iIndex=m_ppars.iIndex;
		
  Box locRegion = a_rr.getRegion() & a_teData.getRegion() & a_NData.getRegion();

  if (!locRegion.isEmpty())
    {
      int porder = 9;
      Real arg, tvar, nel, mue;

      a_rr.setVal(0e0);
      const Box& region = a_rr.getRegion();
      const IntVectSet ivsBox(region);

      
      for (VoFIterator vofit(ivsBox, ebis.getEBGraph());vofit.ok(); ++vofit)
	{
	  const VolIndex& vof = vofit();
	  
	  if (locRegion.contains(vof.gridIndex()))
	    {
	      nel = a_NData(vof, m_eindex);

	      
	      if(m_ppars.isConstMobDiff)
		{

		  Real Emag=0; for (int j = 0; j < SpaceDim; j++)  Emag += pow(a_GphiData(vof, j),2); Emag = sqrt(Emag);
		  Real mve=m_ppars.constMobEle;
		  Real alpha=m_ppars.constAlpha;
		  Real nu = alpha*mve*Emag;
		  a_rr(vof, 0) = nu;  //One Single Reaction
		}
	      else
		{
		  for (int r = 0; r < nReactions(); r++)
		    {
		      const Vector<Real>& coef = m_rrates[r];
		      const Vector<Real>& react = m_reactant[r];
		      const Vector<Real>& transpe = m_transp[m_eindex];
		      
		      int destComp = r;
		      int teDataComp = (coef[0] > 0.1)? m_tgComp : m_eeComp;
		      
		      tvar = a_teData(vof, teDataComp);
		      if (coef[0]==-1) 
			{
			  tvar = (m_constTe == 0)  ?   convertEe2Te(a_teData(vof,teDataComp),nel) :  m_Te;
			    
			  if ( tvar > coef[1]) 
			    {
			      arg = coef[2]; for (int k=1; k< porder+1; k++) arg = arg + coef[k+2]/pow(tvar,k);
			      a_rr(vof, destComp) = exp(arg); // otherwise zero per initialization
			    }	
			}
		      else if (coef[0]==0)
			{
			  tvar = convertEe2Te(a_teData(vof,teDataComp),nel);
			  a_rr(vof, destComp) = coef[1]*pow(tvar,coef[2]);
			}				    				
		      else 
			a_rr(vof, destComp) = coef[1]*pow(tvar,coef[2])*exp(coef[3]/tvar);
		  
		      Real Ng;
		      if (m_constBg == 0)
			Ng = a_NData(vof, m_bgindex);// MDCancel for air + a_NData(vof, m_o2index);	
		      else 
			Ng = m_nBg;
		  
		      // Ionization
		      if (r == 0)
			{
			  // local E-field amplitude
			  Real avE=0; for (int j = 0; j < SpaceDim; j++)  avE += pow(a_GphiData(vof, j),2); avE = sqrt(avE);

			  // electron mobility
			  Real lte = log(tvar);
			  Real arg = transpe[1]; for (int k=1; k< 8+1; k++) arg = arg + transpe[k+1]*pow(lte,k);
			  mue = exp(arg)/(Ng*m_ppars.NA);
				
			  // ionization rate custom correction factor for Helium Roy POP06 system
			  Real EiPcmt=avE/m_GasPressure; //[V/cm/Torr]
			  Real alp=cA*exp(-cB/pow(EiPcmt,0.4))*m_GasPressure;
				
			  a_rr(vof, destComp) = alp*mue*avE;
			}
		      //Recombination
		      else if (r == 1)
			{
			  a_rr(vof, destComp) = 1.12e-7+2.2e-27*1e-6*(Ng*m_ppars.NA);
			}
		      // end of Helium adder
		      for ( int s = 0; s < nSpecies(); s++) 				
			a_rr(vof, destComp)*= pow(m_ppars.NA * a_NData(vof, s),react[m_species[s]]);
		  
		      a_rr(vof, destComp) *= pow(Ng * m_ppars.NA,react[nSpecies()]) / m_ppars.NA;
		    }
		}
	    } //if (locRegion.contains(vof.gridIndex()))
	} //for (VoFIterator vofit(ivsBox, ebis.getEBGraph());vofit.ok(); ++vofit)

    }// if (!locRegion.isEmpty())
}

/*******************
 Poisson RHS
******/

void PlasmaPhysics::ElectricPotentialRHS(LevelData<EBCellFAB>&			a_rhs,
					 const LevelData<EBCellFAB>& 	a_NData,
					 const Real&						a_dx)

{
  CH_assert(a_NData.nComp() >= nSpecies()); 
  Real eps;
  int yEBindex = m_yEB/a_dx;
  //Initialize
  EBLevelDataOps::setVal(a_rhs, 0.0);
  Vector<Real> coef = m_ppars.chargeNr;

  for (DataIterator dit = a_rhs.dataIterator(); dit.ok(); ++dit)
    {
      //DataIndex d = dit();
      EBCellFAB& rhs = a_rhs[dit()];	
      const EBCellFAB& NData = a_NData[dit()];
      //CH_assert(rhs.isDefined() && NData.isDefined);
      const EBISBox&   ebis   = rhs.getEBISBox(); 
      const EBISBox&   NDataEbis   = NData.getEBISBox(); 
      //MDProb: CH_assert(ebis == NDataEbis);

      Box locRegion = NData.getRegion() & rhs.getRegion();
      if (!locRegion.isEmpty())
	{			
	  BaseFab<Real>& rhsRegFAB = rhs.getSingleValuedFAB();
	  const BaseFab<Real>& NDataRegFAB = NData.getSingleValuedFAB();
	  //bool sameRegBox = (rhsRegFAB.box() == NDataRegFAB.box());
	  // bool nvofTest   = false; no vof test done here, add later with pointer calc
	  bool regionTest = false;
	  BaseIVFAB<Real>& rhsIrrBFAB = rhs.getMultiValuedFAB();
	  const BaseIVFAB<Real>& NDataIrrBFAB = NData.getMultiValuedFAB();				
							
	  //if (locRegion == NData.getRegion() && locRegion == rhs.getRegion())
	  //  regionTest = true;
				
	  for (int s = 0; s < nSpecies(); s++)
	    {			
	      FORT_RHSFAB(CHF_FRA1(rhsRegFAB, 0),
			  CHF_CONST_FRA(NDataRegFAB),
			  CHF_CONST_VR(coef),
			  CHF_CONST_VR(m_species),
			  CHF_CONST_REAL(m_epsg),
			  CHF_CONST_REAL(m_epsd),
			  CHF_CONST_INT(yEBindex),
			  CHF_CONST_INT(s),
			  CHF_BOX(locRegion));
				
	      // Multi-Cell
	      IntVectSet ivsMulti = NData.getMultiCells();
	      ivsMulti &= rhs.getMultiCells();  // opt might be to take this out
	      ivsMulti &= locRegion;
	      IVSIterator ivsit(ivsMulti);		
	
	      for (ivsit.reset(); ivsit.ok(); ++ivsit)
		{
		  const IntVect& iv = ivsit();
		  Vector<VolIndex> vofs = ebis.getVoFs(iv);
		  for (int ivof = 0; ivof < vofs.size(); ivof++)
		    {
		      const VolIndex& vof = vofs[ivof];
		      if (locRegion.contains(vof.gridIndex()))
			{   
			  const IntVect& iv = vof.gridIndex();
			  if(iv[1] < yEBindex)
			    eps = m_epsd;
			  else
			    eps = m_epsg;	
			  rhsIrrBFAB(vof, 0) -= m_ppars.q*coef[m_species[s]]*NDataIrrBFAB(vof, s) * m_ppars.NA /eps;														
			}
		    }
		}		
	    }	
	}
    } 
}


void PlasmaPhysics::printCoefficients(LevelData<EBFluxFAB>&     a_fluxMobData,
				      const DisjointBoxLayout&  a_grids,
				      const int&                a_ispec,
				      const int&                a_flag)
{
  int idir = 1;
  int nspec = Min(a_ispec,a_fluxMobData.nComp()-1);
  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      EBFluxFAB&       fluxMobData = a_fluxMobData[dit()];
      EBFaceFAB&       faceMobData = fluxMobData[idir];
      BaseFab<Real>&       regFaceMobData = faceMobData.getSingleValuedFAB();
      const Box&       dblBox = a_grids.get(dit());
      Box cellBox = dblBox;
      Box faceBox = cellBox;
      faceBox.surroundingNodes(idir);
      FORT_PRINT(CHF_FRA1(regFaceMobData, nspec),
		 CHF_CONST_INT(idir),
		 CHF_BOX(faceBox),
		 CHF_CONST_INT(a_flag));
    }
}
//******************************************//
void PlasmaPhysics::maxDT(LevelData<EBCellFAB>&         a_cellSData,
			  const LevelData<EBCellFAB>&   a_cellNData,
			  const LevelData<EBCellFAB>&   a_cellTData,
			  const Real&                   a_dt)
{
  int isource = 0;
  int s = findSpec(string("e-"));
  int tcomp;
  if(m_transp[s][0] == 0)
    tcomp = m_tgComp;
  else
    tcomp = m_eeComp;

  //if evaluating DT_limit/DT pass eps/dt instead than eps
  Real epsEff = m_epsg/a_dt;
  
  for (DataIterator dit = a_cellSData.dataIterator(); dit.ok(); ++dit)
    {
      const EBCellFAB& cellNData = a_cellNData[dit()];
      const EBCellFAB& cellTData = a_cellTData[dit()];
      EBCellFAB& cellSData = a_cellSData[dit()];
      BaseFab<Real>&       regSData = cellSData.getSingleValuedFAB();
      const BaseFab<Real>& regNData = cellNData.getSingleValuedFAB();
      const BaseFab<Real>& regTData = cellTData.getSingleValuedFAB();
      Box cellBox = cellSData.getRegion(); //a_grids.get(dit());
      FORT_MAXDT(CHF_FRA1(regSData, 0),
		 CHF_CONST_FRA1(regTData,tcomp),
		 CHF_CONST_FRA(regNData),
		 CHF_CONST_REAL(epsEff),
		 CHF_BOX(cellBox),
		 CHF_CONST_VR(m_transp[s]),
		 CHF_CONST_REAL(m_nBg),
		 CHF_CONST_REAL(m_Te),
		 CHF_CONST_INT(m_eindex),
		 CHF_CONST_INT(m_bgindex),
		 CHF_CONST_INT(m_o2index),
		 CHF_CONST_INT(m_constBg),
		 CHF_CONST_INT(m_constTe));      
    }
}

/*****************************
 Transport Coefficients
******/
//new 08/15/13 added from LM
/*************************/      
void PlasmaPhysics::transportCoefficients(LevelData<EBFluxFAB>&        a_fluxMobData,
					  LevelData<EBFluxFAB>&        a_fluxDiffData,
					  const LevelData<EBCellFAB>&  a_cellNData,
					  const LevelData<EBCellFAB>&  a_cellTData,
					  const LevelData<EBCellFAB>&  a_cellphiData,
					  const LevelData<EBFluxFAB>&  a_GphiData,
					  const DisjointBoxLayout&     a_grids,
					  const EBISLayout&            a_ebisl,
					  const ProblemDomain&         a_domain,
					  const Real&		       a_dx)	
{
  // commented out for now for testing of selected species subsets 
  // CH_assert((a_fluxMobData.nComp() >= nSpecies()) && (a_fluxDiffData.nComp() >= nSpecies()) );

  IntVect ghostCellphiIn = a_cellphiData.ghostVect();
  IntVect ghostCellTIn = a_cellTData.ghostVect();
  IntVect ghostCellNIn = a_cellNData.ghostVect();
  IntVect ghostFluxMobIn = a_fluxMobData.ghostVect();
  IntVect ghostFluxDiffIn = a_fluxDiffData.ghostVect();

  
  int porder = m_ppars.porder;  
  int nspec = nSpecies();
  bool electron = false;
  

  CH_assert(ghostCellphiIn[0]==ghostCellphiIn[1]);
  CH_assert(ghostCellphiIn[0]==ghostCellphiIn[SpaceDim-1]);
  CH_assert(ghostCellTIn[0]==ghostCellTIn[1]);
  CH_assert(ghostCellTIn[0]==ghostCellTIn[SpaceDim-1]);
  CH_assert(ghostCellNIn[0]==ghostCellNIn[1]);
  CH_assert(ghostCellNIn[0]==ghostCellNIn[SpaceDim-1]);

  CH_assert((ghostCellphiIn[0]==ghostCellTIn[0]) && (ghostCellphiIn[0]==ghostCellNIn[0]));

  CH_assert(ghostFluxMobIn[0]==ghostFluxMobIn[1]);
  CH_assert(ghostFluxMobIn[0]==ghostFluxMobIn[SpaceDim-1]);
  CH_assert(ghostFluxDiffIn[0]==ghostFluxDiffIn[1]);
  CH_assert(ghostFluxDiffIn[0]==ghostFluxDiffIn[SpaceDim-1]);

  CH_assert(ghostFluxMobIn[0]==ghostFluxDiffIn[0]);

  int ghostFluxTan   = ghostCellphiIn[0];
  ghostFluxTan = std::min(ghostFluxMobIn[0],ghostFluxTan);
  CH_assert(ghostFluxMobIn[0] >= ghostFluxTan);

  //Extrapolate temeprature and number density to the faces
  EBFluxFactory fluxFact(a_ebisl);
  LevelData<EBFluxFAB>  tempFace(a_grids, a_cellTData.nComp(),   IntVect::Zero, fluxFact);
  LevelData<EBFluxFAB>  NFace(a_grids, a_cellNData.nComp(),   IntVect::Zero, fluxFact);
  EBLevelDataOps::averageCellToFacesAll(tempFace, a_cellTData, a_grids, a_ebisl, a_domain);
  EBLevelDataOps::averageCellToFacesAll(NFace, a_cellNData, a_grids, a_ebisl, a_domain);

  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      const EBCellFAB& cellphiData = a_cellphiData[dit()];
      const EBGraph&   ebgraph = a_ebisl[dit()].getEBGraph();
      const Box&       dblBox = a_grids.get(dit());
      EBFluxFAB&       fluxMobData = a_fluxMobData[dit()];
      EBFluxFAB&       fluxDiffData = a_fluxDiffData[dit()];

      for (int idir = 0; idir < SpaceDim; idir++)
        {
	  EBFaceFAB&       faceMobData = fluxMobData[idir];
	  EBFaceFAB&       faceDiffData = fluxDiffData[idir];
	  const EBFaceFAB& faceTempData = tempFace[dit()][idir];
	  const EBFaceFAB& faceNData =       NFace[dit()][idir];
	  const EBFaceFAB& faceGphiData = a_GphiData[dit()][idir];
	  FaceStop::WhichFaces stopCritGrid =  FaceStop::SurroundingWithBoundary;
	  const Box& region = faceTempData.getCellRegion();
	  IntVectSet ivsBox(region);
	  const EBISBox& ebisBox = faceMobData.getEBISBox();
	   for (FaceIterator faceit(ivsBox, ebisBox.getEBGraph(), idir, stopCritGrid); faceit.ok(); ++faceit)
	     {
	       int TvarComp;
	       Real Eface, Ngface, Tface, ENface, nsface, neface, mvar, mob, zsign;
	       Vector<Real > transp;
	       
	       const FaceIndex& face = faceit();
	       Eface = -faceGphiData(face,0);
	       if (m_constBg == 0)
		 Ngface = faceNData(face, m_bgindex);   //MDCancel for air: + regCellNData(ivCell, m_o2index);
	       else
		 Ngface = m_nBg;
	       neface = faceNData(face, m_eindex);  
	       ENface = abs(Eface) / (Ngface * m_ppars.NA); //reduced field
	       if (ENface < 1.0e-21) ENface = 1.0e-21;  //Mah ??
	       for (int s = 0; s < nspec; s++)
		 {
		   zsign = sgn(m_ppars.chargeNr[m_species[s]]);
		   if(m_ppars.isConstMobDiff)
		     {
		       electron = (s == m_eindex);
		       if(electron){
			 faceMobData(face, s) =  m_ppars.constMobEle* zsign*Eface;
			 faceDiffData(face, s) = m_ppars.constDiffEle;
		       }
		       else if (m_ppars.isNeutral[s]) {
			 faceMobData(face, s) =  0;
			 faceDiffData(face, s) = m_ppars.constDiffIon;
		       }
		       else
			 {
			 faceMobData(face, s) =  m_ppars.constMobIon* zsign*Eface;
			 faceDiffData(face, s) = m_ppars.constDiffIon;
			 }
		       
		     }
		   else
		     //Denison's Model
		     {
		       //set the temperature to be used in the rate terms
		       transp = m_transp[s];
		       if(transp[0] == 0)
			 {
			   electron = false;
			   Tface = faceTempData(face, m_tgComp);
			   mvar = log(ENface);
			 }
		       else
			 {
			   Tface =convertEe2Te(faceTempData(face, m_eeComp), neface,  m_constTe); 
			   
			   electron = true;
			   mvar = log(Tface);
			 }		
		       if (transp[14] == 0)
			 {
			   // Ions or Electrons, calculate both mobility and diffusion:				
			   Real mob = transp[1];			    
			   for (int k=1; k< porder+1; k++)mob = mob + transp[k+1]*pow(mvar,k);						
			   mob = exp(mob)/(Ngface*m_ppars.NA);
			   
			   faceMobData(face, s) = mob * Eface * zsign;
			   faceDiffData(face, s) = mob * Tface * transp[11];
			   if (electron)
			     {
			       faceMobData(face, nspec) = 5.0/3.0 * mob * Eface * zsign;
			       faceDiffData(face, nspec) = - 5.0/3.0 * mob * Tface * transp[11];
			     }
			 }
		       else			
			 {
			   // Neutral, no mobility 
			   faceDiffData(face, s) = pow(Tface,transp[12]) * transp[11] / (Ngface*m_ppars.NA);
			 }
		     }
		 }
	     }
        }
    }
  EBLevelDataOps::exchangeAll(a_fluxMobData);
  EBLevelDataOps::exchangeAll(a_fluxDiffData);
}


void PlasmaPhysics::transportCoefficients(LevelData< BaseIVFAB<Real> >&   	a_fluxMobData,
					  LevelData< BaseIVFAB<Real> >&   	a_fluxDiffData,
					  LevelData< BaseIVFAB<Real> >&   	a_fluxWallData,
					  LevelData< BaseIVFAB<Real> >&   	a_rhosRHS,
					  const LevelData<EBCellFAB>&           a_cellNData,
					  const LevelData<EBCellFAB>&           a_cellTData,
					  const LevelData<BaseIVFAB<Real>>&     a_IVGPhiData,
					  const LevelData<EBCellFAB>&           a_cellGPhiData,
					  const DisjointBoxLayout&              a_grids,
					  const EBISLayout&                     a_ebisl,
					  const ProblemDomain&                  a_domain,
					  const Real&                           a_dx)
											
{
  int nspec = nSpecies();
  CH_assert((a_fluxMobData.nComp() >= nspec+1) && (a_fluxDiffData.nComp() >= nspec+1) && (a_fluxWallData.nComp() >= nspec+1+2));
  Real  flux, buffe, buffbg, buffo2,  pi, Eface, Tface, Ngface, ENface, nsface, neface, nbg, no2, mvar, mob, zsign, fluxwall, thermVelFact, arg, alpe;
  pi = 4.0*atan(1.0);
  int sel, porder;
  porder = m_ppars.porder;  
  RealVect dxv = a_dx*RealVect::Unit;
  Vector<Real > transp;
  bool electron;
  bool debug = false;

  Vector<Real > rhosFact(nspec,0);
  for (int s = 0; s < nspec; s++) rhosFact[s] = m_ppars.chargeNr[m_species[s]]*m_ppars.q * m_ppars.NA;

  // all the derivatives are outgoing, the normal is pointing out of the domain
  // for consistency with current species solver, keeping  species flux/N = flux/N * (-ey) = - sign(z) * mu * E * ey , E=-grad  phi
  // so that species flux is the outgoing component for negative charges in a positive E-field and vice-versa

  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      const EBGraph&    ebgraph = a_ebisl[dit()].getEBGraph();
      const EBISBox& ebisBox = a_cellNData[dit()].getEBISBox();  

      const EBCellFAB& cellTData = a_cellTData[dit()];
      const EBCellFAB& cellNData = a_cellNData[dit()];
      const EBCellFAB& cellGPhiData = a_cellGPhiData[dit()];
      const BaseIVFAB<Real>& IVGPhiData = a_IVGPhiData[dit()];
      BaseIVFAB<Real>& fluxMobData = a_fluxMobData[dit()];
      BaseIVFAB<Real>& fluxDiffData = a_fluxDiffData[dit()];
      BaseIVFAB<Real>& fluxWallData = a_fluxWallData[dit()];
      BaseIVFAB<Real>& rhosRHS = a_rhosRHS[dit()]; //it was init to 0 in EBAMRSpecies::setDiffCoeff
	
      IntVectSet ivsIrreg = fluxMobData.getIVS();
      for (VoFIterator vofit(ivsIrreg, ebgraph); vofit.ok(); ++vofit)
	{
	  VolIndex vof = vofit();
	  RealVect bndryCentroid = ebisBox.bndryCentroid(vof);
	  RealVect normal = ebisBox.normal(vof);
	  const IntVect& iv = vof.gridIndex();   
	
	  bndryCentroid *= a_dx;
	  VoFStencil extrapSten, normDSten, zeroSten;
	  Real weight;
	  int order = EBArith::getExtrapolationStencil(extrapSten, bndryCentroid, dxv, vof, ebisBox);
	  //int order = EBArith::getFirstOrderExtrapolationStencil(extrapSten, bndryCentroid, dxv, vof, ebisBox);
	  EBArith::getLeastSquaresGradSten(normDSten, weight, vof, ebisBox, dxv, a_domain);
	  		
	  // For DC applyVoFStencil is to the outside of the domain into the wall ==> reversing sign to get E=-grad phi
	  if(isElectrode(vof,a_dx,bndryCentroid))
	    Eface =  -IVGPhiData(vof,0);  //Points in the solid
	  else
	    Eface = -normal.dotProduct(applyVoFStencilVect(extrapSten, cellGPhiData)); //Luca :: Check the direction of the normal
	  if (m_constBg == 0)
	    Ngface = applyVoFStencilComp(extrapSten, a_cellNData[dit()], m_bgindex); //MDCancel for air:  + applyVoFStencilComp(extrapSten, a_cellNData[dit()], m_o2index);
	  else
	    Ngface = m_nBg;
	  
	  ENface = abs(Eface) / (Ngface*m_ppars.NA);
	  neface = applyVoFStencilComp(extrapSten, a_cellNData[dit()], m_eindex); 
	  if(neface < m_minData[m_eindex]) neface =  Max(a_cellNData[dit()](vof,m_eindex),m_minData[m_eindex]);
	    
	  Tface=convertEe2Te(applyVoFStencilComp(extrapSten, a_cellTData[dit()], m_eeComp), neface, m_constTe);
	  
	  sel = -1;
	  buffe = 0;
	  buffbg = 0;
	  buffo2 = 0;

	  
	  for (int s = 0; s < nspec; s++) fluxWallData(vof, s) =  0e0;
	  nbg = 0;
	  for (int s = 0; s < nspec; s++)
	    {	
	      int S = m_species[s];
	      transp = m_transp[S]; 
	      thermVelFact = 0.25*sqrt(8.0*m_ppars.kb*m_ppars.NA*Tface/(pi*m_ppars.mass[S]));
	      zsign = sgn(m_ppars.chargeNr[S]);
	      
	      nsface = applyVoFStencilComp(extrapSten, a_cellNData[dit()], s);
	      if(nsface < m_minData[s]) nsface =  Max(a_cellNData[dit()](vof,s),m_minData[s]);	
	      
	      if(m_ppars.isConstMobDiff)
		{
		  electron = (s == m_eindex);
		  if(electron){
		    fluxMobData(vof, s) =  m_ppars.constMobEle* zsign*Eface;
		    fluxDiffData(vof, s) = m_ppars.constDiffEle;
		  }
		  else if (m_ppars.isNeutral[s]) {
		    fluxMobData(vof, s) =  0;
		    fluxDiffData(vof, s) = m_ppars.constDiffIon;
		  }
		  else
		    {
		      fluxMobData(vof, s) =  m_ppars.constMobIon* zsign*Eface;
		      fluxDiffData(vof, s) = m_ppars.constDiffIon;
		    }
		  // assume the flux to be the mobility times the concentration out of the fluid
		  if(fluxMobData(vof, s) > 0.0)
		    {
		      fluxWallData(vof, s) += fluxMobData(vof, s);
		      if (s != m_eindex)
			fluxWallData(vof, m_wallBCInfo[m_eindex] ) -= m_ppars.secEmElec[S]*fluxMobData(vof, s)*nsface;
		    }
		  else
		    fluxWallData(vof, s) =0;
	
		  
		}
	      else
		{
		  //Denison's Model
		  if(transp[0] == 0)
		    { 
		      electron = false;
		      Tface = Max(applyVoFStencilComp(extrapSten, a_cellTData[dit()], m_tgComp),3e2);
		      mvar = log(ENface);
		    }
		  else
		    {
		      
		      electron = true;
		      sel = s; 
		      mvar = log(Tface);
		    }	
		  
		  // excited neutrals and ions recombine if sticking = 1, or are reflected otherwise, here assume sticking of 1
		  // (one could implement intermediate states by using surface reaction rates)  
		  fluxWallData(vof, s) = 0;
		  Real WallFluxMob = 0;
		  // April 13 
		  if (s != m_bgindex && m_useThermalVelocity) 
		    fluxWallData(vof, s) += thermVelFact;
		  
		  if (transp[14] == 0)
		    {
		      // Ions or Electrons, calculate mobility and diffusion:			
		      arg = transp[1];
		      for (int k=1; k< porder+1; k++) arg += transp[k+1]*pow(mvar,k);
		      mob = exp(arg)/(Ngface*m_ppars.NA);
		      // flux into the wall = - sign(z) * mu * E * ey , E=-grad phi is physical one
		      fluxMobData(vof, s) = - Eface * zsign * mob;
		      fluxDiffData(vof, s) = mob * Tface * transp[11];
		      int wallnorm = -1; 
		      // wallnorm * Eface > 0 positive for ion+ to the wall, or < 0 for e- to the wall
		      WallFluxMob = Max(wallnorm * Eface * zsign * mob, 0e0);
		      fluxWallData(vof, s) += WallFluxMob;
		      
		      if (electron)
			// energy terms
			{
			  fluxMobData(vof, nspec) =  -5.0/3.0 * mob * Eface* zsign;
			  fluxDiffData(vof, nspec) = -5.0/3.0 * mob * Tface * transp[11]; 
			  fluxWallData(vof, nspec) =  4e0/3e0 * thermVelFact;
			}
		      else
			{
			  // surface recombination => neutral flux into domain, only reflect into the domain positive current to the wall
			  // MDCancel for air buffo2 -= m_ppars.stickFac2o2[s] * Min(fluxWallData(vof,s),0e0) * nsface;
			  buffbg -= m_ppars.stickFac2bg[s] * Max(fluxWallData(vof,s),0e0) * nsface;			   				
			}
		    }
		  else			
		    {
		      // Neutrals, diffusion only 
		      fluxDiffData(vof, s) = pow(Tface,transp[12]) * transp[11] / (Ngface*m_ppars.NA);
		      if (s != m_bgindex && s!= m_o2index)
			{
			  // only reflect into the domain positive current to the wall
			  // MDCancel for air buffo2 -= m_ppars.stickFac2o2[s] * Min(fluxWallData(vof,s),0e0) * nsface;
			  buffbg -= m_ppars.stickFac2bg[s] * Max(fluxWallData(vof,s),0e0) * nsface;		  		
			}
		      else 
			{
			  if (s == m_bgindex)
			    nbg = Ngface;
			  //MDCancel for air
			  //else if (s == m_o2index)
			  //	no2 = nsface;
			}
		    }
		  // Secondary emission from ions and excited neutrals
		  // only add if Eface=-gradPhi<0 so that Eface * wallnorm is > 0
		  if (s != m_eindex)
		    {
		      if (Eface < 0) alpe=1; else alpe=0;
		      if (isElectrode(vof,a_dx,bndryCentroid))
			{
			  // only reflect into the domain positive flux to the wall 
			  buffe -= m_ppars.secEmElec[S] * alpe * Max(fluxWallData(vof,s),0e0) * nsface;
			}
		      else
			buffe -= m_ppars.secEmDiel[S] * alpe * Max(fluxWallData(vof,s),0e0) * nsface;
		    }
		  if(std::isnan(buffe) || std::isinf(buffe))
		    {
		      pout() << "PlasmaPhysics::transportCoefficients (irreg): s= " << s << ", Eface= "<< Eface <<  ", mob= "<< mob <<  ", buffe= "<< buffe <<  ", fluxwalldata= "<<fluxWallData(vof,s) <<  ", nsface= "<<nsface << endl;
		      pout() << "PlasmaPhysics::transportCoefficients (irreg): S= " << S <<  ", thermVelFact= " << thermVelFact << ", m_ppars.mass[S]= " << m_ppars.mass[S] <<  ", Tface= " << Tface <<  ", pi= " << pi << ", m_ppars.kb= "<< m_ppars.kb <<  ", m_ppars.NA= "<< m_ppars.NA << endl;
		    }
		}//Denison Model
	      //LMadd: make the right hand side of the rhos equation
	      rhosRHS(vof,0) += fluxWallData(vof,s)*nsface*rhosFact[s];
	    }// s
	  
	  // Add the buffers from surface recombinations and secondary emission
	  if(sel >0) // denison's model
	    {
	      fluxWallData(vof, nspec+1) = buffe;
	      fluxWallData(vof, nspec+2) = buffe * m_ppars.NA * m_ppars.q * m_ppars.secEmEnergy[m_species[sel]] ;
	      if ((nbg > 0) && (m_constBg ==0)) fluxWallData(vof, m_bgindex) += buffbg / nbg;
	    }
	  
	  if(isElectrode(vof,a_dx,bndryCentroid)) rhosRHS(vof,0) =0;						
	}// vofiterator
    }// dataiterator
	  
  a_fluxMobData.exchange(Interval(0,a_fluxMobData.nComp()-1));
  a_fluxDiffData.exchange(Interval(0,a_fluxDiffData.nComp()-1));
  a_fluxWallData.exchange();
}

/*******************
 Rhos ODE RHS
******/
void PlasmaPhysics::rhos(LevelData< BaseIVFAB<Real> >&          a_rhos,
			 const LevelData< BaseIVFAB<Real> >&   	a_rhosRHS,
			 const DisjointBoxLayout&      		a_grids,
			 const EBISLayout&             		a_ebisl,
			 const ProblemDomain&          		a_domain,
			 const int&          			a_ghost,			
			 const Real&          			a_dt)
											
{
  Real rhs;
  CH_assert(a_rhosRHS.nComp() == 1);
  CH_assert(a_dt > 0);

  // a_rhosRHS is a buffer filled by this->transportCoefficients, please check in there for flux convention and units
  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      Box grownBox = grow(a_grids.get(dit()), a_ghost);
      grownBox &= a_domain;
      const EBGraph&    ebgraph = a_ebisl[dit()].getEBGraph(); 
      BaseIVFAB<Real>& rhos = a_rhos[dit()];
      const BaseIVFAB<Real>& rhosRHS = a_rhosRHS[dit()];
      IntVectSet ivsIrreg = ebgraph.getIrregCells(grownBox);

      for (VoFIterator vofit(ivsIrreg, ebgraph); vofit.ok(); ++vofit)
	{
	  VolIndex vof = vofit();
	  rhos(vof, 0) += a_dt*rhosRHS(vof, 0);	
	  
	  				
	}// vofiterator
    }// dataiterator

  a_rhos.exchange(Interval(0,0));
}


/****************/

Real PlasmaPhysics::applyVoFStencilComp(const VoFStencil& a_sten, const EBCellFAB& a_fab, const int& a_comp)
{
  Real retval = 0.;
  for (int isten = 0; isten < a_sten.size(); isten++)
    {
      retval += (a_sten.weight(isten))*(a_fab((a_sten.vof(isten)), a_comp));
    }
  return retval;
}
RealVect PlasmaPhysics::applyVoFStencilVect(const VoFStencil& a_sten, const EBCellFAB& a_fab)
{
  RealVect retval=RealVect::Zero;
  for (int icomp = 0; icomp < SpaceDim; icomp++)
    {
      for (int isten = 0; isten < a_sten.size(); isten++)
	{
	  retval[icomp] += (a_sten.weight(isten))*(a_fab((a_sten.vof(isten)), icomp));
	}
    }
  return retval;
}

/*********************************************************************
 Electron Energy source Term: Joule heating contribution (only)
 *****************************************/

// Routine adds only the leading part of the joule heating to the source term, called before adding other sources by call of plasmaSources

void PlasmaPhysics::addJouleHeat(LevelData<EBCellFAB>&         a_cellSData,
				 const LevelData<EBCellFAB>&   a_cellNData,
				 const LevelData<EBCellFAB>&   a_cellTData,
				 const LevelData<EBCellFAB>&   a_cellEData,
				 const LevelData<EBCellFAB>&   a_cellGNData,
				 const LevelData<EBCellFAB>&   a_cellVelData,
				 const DisjointBoxLayout&      	a_grids,
				 const EBISLayout&             	a_ebisl)
											
{
  // add it to the electron energy source term
  
  int isource = findSpec(string("ee"));
  int s = findSpec(string("e-"));
  int tcomp;
  if(m_transp[s][0] == 0)
    tcomp = m_tgComp;
  else
    tcomp = m_eeComp;
  
  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      //EBs
      const EBCellFAB& ebE = a_cellEData[dit()];
      const EBCellFAB& ebT = a_cellTData[dit()];
      const EBCellFAB& ebN = a_cellNData[dit()];
      const EBCellFAB& ebGN = a_cellGNData[dit()];
      const EBCellFAB& ebVel = a_cellVelData[dit()];
      EBCellFAB&       ebS = a_cellSData[dit()];

      //ebS.setVal(0e0);
      if (m_constTe == 0)
	{
	  //Basefabs
	  const BaseFab<Real>& regE = ebE.getSingleValuedFAB();
	  const BaseFab<Real>& regT = ebT.getSingleValuedFAB();
	  const BaseFab<Real>& regN = ebN.getSingleValuedFAB();
	  const BaseFab<Real>& regGN = ebGN.getSingleValuedFAB();
	  const BaseFab<Real>& regVel = ebVel.getSingleValuedFAB();
	  BaseFab<Real>&       regS = ebS.getSingleValuedFAB();

	  Box cellBox = ebS.getRegion() & ebT.getRegion() & ebE.getRegion(); //a_grids.get(dit());
        
	  FORT_JOULEHEAT(CHF_FRA1(regS, isource),
			 CHF_CONST_FRA(regE),
			 CHF_CONST_FRA1(regT,tcomp),
			 CHF_CONST_FRA(regN),
			 CHF_CONST_FRA(regGN),
			 CHF_CONST_FRA(regVel),
			 CHF_BOX(cellBox),
			 CHF_CONST_VR(m_transp[s]),
			 CHF_CONST_REAL(m_nBg),
			 CHF_CONST_REAL(m_Te),
			 CHF_CONST_INT(m_eindex),
			 CHF_CONST_INT(m_bgindex),
			 CHF_CONST_INT(m_o2index),
			 CHF_CONST_INT(m_constBg),
			 CHF_CONST_INT(m_constTe));
      
	  const EBISBox& ebisBox = ebS.getEBISBox();
	  //Cancel fix up irregular faces
	  const EBGraph&    ebgraph = a_ebisl[dit()].getEBGraph();
	  IntVectSet ivsIrreg = ebgraph.getIrregCells(cellBox);
	  for (VoFIterator vofit(ivsIrreg, ebgraph); vofit.ok(); ++vofit)
	    {
	      VolIndex vof = vofit();
	      ebS(vof,isource) = 0e0;
	    }
	}//MDCancel
    }
}
/************************/
//fluid Joule heating
/*****************/
// Called after fluidSources, do not zero the input ebcellfab

void PlasmaPhysics::AddJouleHeatF(LevelData<EBCellFAB>&         a_cellSData,
				  const LevelData<EBCellFAB>&   a_cellNData,
				  const LevelData<EBCellFAB>&   a_cellTData,
				  const LevelData<EBCellFAB>&   a_cellEData,
				  const LevelData<EBCellFAB>&   a_cellGNData,
				  const DisjointBoxLayout&      a_grids,
				  const EBISLayout&             a_ebisl)
											
{
  // add it to the electron energy source term
  EBPatchPolytropic Ppatch;
  int enrgIndx = Ppatch.energyIndexC();
  //the enthalpy equation 
  int isource = enrgIndx;
  Real etaT = m_ppars.etaT;
  
  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      //EBs
      const EBCellFAB& ebE = a_cellEData[dit()];
      const EBCellFAB& ebT = a_cellTData[dit()];
      const EBCellFAB& ebN = a_cellNData[dit()];
      const EBCellFAB& ebGN = a_cellGNData[dit()];
      EBCellFAB&       ebS = a_cellSData[dit()];

      //MDCancel test
      if (m_constTe == 0)
	{ 
	  //Basefabs
	  const BaseFab<Real>& regE = ebE.getSingleValuedFAB();
	  const BaseFab<Real>& regT = ebT.getSingleValuedFAB();
	  const BaseFab<Real>& regN = ebN.getSingleValuedFAB();
	  const BaseFab<Real>& regGN = ebGN.getSingleValuedFAB();
	  BaseFab<Real>&       regS = ebS.getSingleValuedFAB();

	  Box cellBox = ebS.getRegion() & ebT.getRegion() & ebE.getRegion(); //a_grids.get(dit());
      
	  for (int s = 0; s < nSpecies()-1; s++) // skip the electron
	    {
	      int tcomp;
	      if(m_transp[s][0] == 0)
		tcomp = m_tgComp;
	      else
		tcomp = m_eeComp;
	      Real zsign = sgn(m_ppars.chargeNr[m_species[s]]);
      
	      FORT_JOULEHEATFLUID(CHF_FRA1(regS, isource),
				  CHF_CONST_FRA(regE),
				  CHF_CONST_FRA1(regT,tcomp),
				  CHF_CONST_FRA(regN),
				  CHF_CONST_FRA(regGN),
				  CHF_BOX(cellBox),
				  CHF_CONST_VR(m_transp[s]),
				  CHF_CONST_INT(m_eindex),
				  CHF_CONST_INT(m_bgindex),
				  CHF_CONST_INT(m_o2index),
				  CHF_CONST_INT(m_constBg),
				  CHF_CONST_INT(m_constTe),
				  CHF_CONST_INT(s),
				  CHF_CONST_REAL(zsign),
				  CHF_CONST_REAL(m_nBg),
				  CHF_CONST_REAL(m_Te),
				  CHF_CONST_REAL(etaT));
	    }
     
	}//MDcancel
    }
}

/******************/
//now unused, in the future for multivals

void PlasmaPhysics::electronTemperature(LevelData<EBCellFAB>&         a_cellTData,
					const LevelData<EBCellFAB>&   a_cellNData,
					const DisjointBoxLayout&      	a_grids,
					const EBISLayout&             	a_ebisl) 

										
{
	 
  int snrg = findSpec(string("ee"));
  int sden = findSpec(string("e-"));
  
  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      //EBs
      const EBCellFAB& ebN = a_cellNData[dit()];
      EBCellFAB&       ebT = a_cellTData[dit()];
      ebT.setVal(300e0);
      //Basefabs
      const BaseFab<Real>& regN = ebN.getSingleValuedFAB();
      BaseFab<Real>& regT = ebT.getSingleValuedFAB();

      Box cellBox = ebN.getRegion() & ebT.getRegion(); //a_grids.get(dit());    
      
      FORT_ELECTEMP(CHF_FRA1(regT, 0),
		    CHF_CONST_FRA1(regN,sden),
		    CHF_CONST_FRA1(regN,snrg),
		    CHF_CONST_REAL(m_Te),
		    CHF_CONST_INT(m_constTe),
		    CHF_BOX(cellBox));
      
    }
}

/*****************/

int PlasmaPhysics::findSpec(const string& a_spec)
{
  int retval=-1;
  for(int k=0;k< m_charSpec.size(); k++)
    {
      if(a_spec.compare(m_charSpec[k]) == 0)
	{
	  retval = k;
	  break;
	}
    }
  return retval;
}

/*****************/

// calculates electron source terms 

void PlasmaPhysics::plasmaSourcesODE(LevelData<EBCellFAB>& 	  a_NData,
				     const LevelData<EBCellFAB>&  a_teData,
				     const LevelData<EBCellFAB>&  a_fGphiData,
				     const LevelData<EBCellFAB>&  a_GphiData,				
				     const EBISLayout&	          a_ebisl,
				     const Real&                  a_dt,
				     const int&                   a_level)
{
  int flag=0; // the CVODE flag
  for (DataIterator dit = a_NData.dataIterator(); dit.ok(); ++dit)
    {

      EBCellFAB& NData = a_NData[dit()];
      const EBCellFAB& teData = a_teData[dit()];
      const EBCellFAB& fGphiData = a_fGphiData[dit()]; //not used becasue the joule terms are evaluated separately
      const EBCellFAB& GphiData = a_GphiData[dit()];

      Box box = NData.getRegion(); // a_NData.box(dit()); //to exclude ghost cells
      const EBISBox& ebisBox = a_ebisl[dit()];
      IntVectSet ivs(box);
      for (VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
	{
	  VolIndex vof = vofit();
	  //setup variables
	  m_GasTemperature = teData(vof,0);
	  m_odeSource = a_fGphiData[dit()](vof,cvode_NEQ-1);
	  //MDCancel: temporarily for Helium 
	  Real avE=0;				
	  for (int j = 0; j < SpaceDim; j++)
	    {
	      avE += pow(GphiData(vof, j),2);
	    }
	  m_odeEfield = sqrt(avE);
	
	  //m_print = vof == VolIndex(IntVect(D_DECL(242,12,0)),0) && a_level == 1 && a_NData.box(dit()).contains(vof.gridIndex()); // given vof, level and not a ghost point
	  m_vof = vof;
	  realtype *local_y;
	  local_y = NV_DATA_S(cvode_y);
	  for(int k=0; k< cvode_NEQ; k++)  local_y[k] = NData(vof,k);
	  realtype T0 = 0e0;
	  realtype t=0;
	  /* cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
	     flag = CVodeInit(cvode_mem, cvode_f, T0, cvode_y);
	     flag = CVodeSVtolerances(cvode_mem, cvode_reltol, cvode_abstol);
	     flag = CVDense(cvode_mem, cvode_NEQ);
	     flag = CVDlsSetDenseJacFn(cvode_mem, NULL); */
	  bool redux = false;
	  if(m_odeflag) //ODE
	    {
	      flag = CVodeReInit(cvode_mem, T0, cvode_y);
	      flag = CVode(cvode_mem, a_dt, cvode_y, &t, CV_NORMAL);
	      /* CVodeFree(&cvode_mem); */
	      // failure check
	      if(flag < 0 || m_print) 
		{
		  pout() << "Level " << a_level << " VOF " << vof << ", in " << a_NData.box(dit()).contains(vof.gridIndex()) << ", E " << m_odeEfield << ", isIrregular? " << ebisBox.isIrregular(vof.gridIndex()) << " Boundary Area " << ebisBox.bndryArea(vof);
		  for(int k=0; k< cvode_NEQ; k++) pout() << ", N("<< k << ")=" << NData(vof,k); pout() << endl;
		  pout() <<flag << ", end, t= " << t << vof;for(int k=0; k< cvode_NEQ; k++) pout() << ", " << local_y[k]; pout() << endl;
		  for(int k=0; k< cvode_NEQ; k++)  local_y[k] = abs(NData(vof,k)/2.0);
		  redux = false;
		  //if(flag < 0) MayDay::Error("PlasmaPhysics::NAN in plasmaSourcesODE");
		}
	    }
	  if(!m_odeflag || redux)
	    {  // FWD EUler
	      flag=cvode_f(t, cvode_y, cvode_abstol);
	      realtype *local_ydot;
	      local_ydot = NV_DATA_S(cvode_abstol);
	      for(int k=0; k< cvode_NEQ; k++)  local_y[k] += local_ydot[k]*a_dt;
	    }
	  
	  for(int k=0; k< cvode_NEQ; k++)  NData(vof,k) = Min(Max(local_y[k],m_minData[k]),m_maxData[k]);
	  
	}
    }
  a_NData.exchange(Interval(0,a_NData.nComp()-1));
}

/**************************
 Source Terms for cvode
**************/

int PlasmaPhysics::cvode_f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  //extract cvode data
  int numY=nSpecies();
  realtype *local_y;
  local_y = NV_DATA_S(y);
  realtype *local_ydot;
  local_ydot = NV_DATA_S(ydot);
  Vector<Real> abs_y(numY,0); 
  for (int s = 0; s < numY; s++) abs_y[s]=abs(local_y[s]);
  
  int porder = 9;
  Real  arg=0e0, tvar, mue, mui, alp, Et, cA, cB, lte, Torr2Pa;
  Vector<Real> epar = m_ppars.eInelEn;
  Vector<Real> mpar = m_ppars.mass;
  Vector<Real> sigebg = m_ppars.sigebg; //changed wrt dec31
  Vector<Real> sigeo2 = m_ppars.sigeo2;
  Real etaT = m_ppars.etaT;

  //for electron energy source term
  int se = cvode_NEQ-1;
  CH_assert(se == numY);
  local_ydot[se] = 0; 


  Vector<Real> coef, react, transp;
  Vector <Real> rates(nReactions(),0e0);
  Real nel = local_y[m_eindex];
  //MDCancel un-comment for air : 
  //Real no2 = local_y[m_o2index];
  Real nbg = local_y[m_bgindex];
  Real elEnrgy = local_y[cvode_NEQ-1];
  Real tel, een, eV;
  if (m_constTe == 0)
    {
      tel = PlasmaPhysics::convertEe2Te(elEnrgy,nel); //K
      een = elEnrgy/(nel*m_ppars.NA); //J
    }
  else
    {
      tel = m_Te;
      een = PlasmaPhysics::convertTe2Een(tel);
    }
  eV = een/m_ppars.q;
  Real teg = m_GasTemperature;
  Real Ng;
  if (m_constBg == 0)
    Ng = local_y[m_bgindex];//MDCancel un-comment for air : + local_y[m_o2index];
  else
    Ng = m_nBg;
  Vector<Real> transpe = m_transp[m_eindex];
  Real rion, recb;

  /*-----------------
    Reaction Rates 
    ------*/

  if(m_ppars.isConstMobDiff)
    {
      Real mve=m_ppars.constMobEle;
      Real alpha=m_ppars.constAlpha;
      Real nu = alpha*mve*m_odeEfield;
      rates[0] = abs(nel*nu);
      goto speciesSourceTerms;
    }
  for (int r = 0; r < nReactions(); r++)
    {
      coef = m_rrates[r];
      react = m_reactant[r];
      int destComp = r;
      
      if (coef[0]==-1)
	{
	  tvar = abs(tel);
	  if ( tvar > coef[1]) 
	    {
	      arg = coef[2];						    	
	      for (int m=1; m <= porder; m++) arg = arg + coef[m+2]/pow(tvar,m);								
	      rates[r] = exp(arg);//note initialize to 0 so if false rates=0
	    }
	}
      else if (coef[0]==0)
	{
	  tvar = abs(tel);
	  if(nel>m_minData[m_eindex]) 
	    rates[r] = coef[1]*pow(tvar,coef[2]);
	  else
	    rates[r] = 0; //MDCancel
	}				    				
      else 
	{
	  tvar = teg;
	  rates[r] = coef[1]*pow(tvar,coef[2])*exp(coef[3]/tvar);
	}

      // MDCancel below per Roy POP06 for Helium 	   
      // Ionization
      if (r == 0)
	{				
	  // electron mobility
	  transp = m_transp[m_eindex];
	  lte = log(tvar);		 	
	  arg = transpe[1];			 	
	  for (int k=1; k< 8+1; k++)
	    arg = arg + transpe[k+1]*pow(lte,k);
	  mue = exp(arg)/(Ng*m_ppars.NA);

				
	  // ionization rate custom correction factor for Helium Roy POP06 system
	  Torr2Pa = 133.3224;
	  Et= m_odeEfield/m_GasPressure; //[V/cm/Torr]
	  cA=4.4e2/Torr2Pa;
	  cB=14.0*pow(1e2/Torr2Pa,0.4);
	  alp=cA*exp(-cB/pow(Et,0.4))*m_GasPressure;
				
	  rates[r] = alp*mue*m_odeEfield;	
	  if (std::isnan(rates[r])) 
	    // MDCancel short-cut per Roy POP 12
	    {
	      mue=m_ppars.q/m_ppars.mass[m_species[m_eindex]]/1.0e12;  
	      rates[r] = alp*mue*m_odeEfield;	
	    }						
	  rion = rates[r];		
	}
      //Recombination
      if (r == 1)
	{
	  rates[r] = 1.12e-7+2.2e-27*1e-6*(Ng*m_ppars.NA);
	  recb = rates[r];
	}
      // end of Helium specific ionization rate
      if(std::isnan(rates[r]) || std::isinf(rates[r]))
	pout() << "PlasmaPhysics::cvode_f [SI]: r= " << r <<  ",m_vof= "<< m_vof << ", rates[r]= "<< rates[r] << ", mue= "<< mue << ", alp= "<< alp 
	       << ", tel= "<< tel << ", tvar= "<< tvar << ", nel= "<< nel << ", m_odeEfield= "<< m_odeEfield << ", een= "<< een << endl;
		
      for ( int s = 0; s < numY; s++)
	{					
	  rates[r] *= pow(abs_y[s]*m_ppars.NA,react[m_species[s]]);	
	}
	
      rates[r] *= pow(Ng*m_ppars.NA,react[nSpecies()])/m_ppars.NA; 
    }//end r it

  /*-----------------------
    Species source terms
    ----------*/
 speciesSourceTerms:
  for (int s = 0; s < numY; s++)
    {
      Vector<Real> par = m_netReac[m_species[s]];
      
      local_ydot[s] = 0;
      for (int r = 0; r < nReactions(); r++)
	{
	  //if (eV > abs(epar[r])) 		
	  local_ydot[s] += rates[r] * par[r];
	}
    }
    
  /*
    if(m_print)
    pout() << "PlasmaPhysics::cvode_f: tel= "<< tel << ", Vof " << m_vof 
    << ", local_ydot[0]= "<< local_ydot[0] <<  ", local_y[0]= "<< local_y[0] 
    << ", local_ydot[1]= "<< local_ydot[1] <<  ", local_y[1]= "<< local_y[1] 
    << ", local_ydot[2]= "<< local_ydot[2] <<  ", local_y[2]= "<< local_y[2]
    << ", rates[0]= "<< rates[0] <<  ", rates[1]= "<< rates[1]
    << ", alp= "<< alp <<  ", m_odeEfield= "<< m_odeEfield 
    << ", rion= "<< rion <<  ", recb= "<< recb <<  ", Ng= "<< Ng 
    << ", nel= "<< nel << ", eV= "<< een/m_ppars.q << ", mue= "<< mue 
    << endl;
  */
	
  /*---------------------------------
    Electron energy source term 
    ------------------*/

  Real buff = 0; Real rateInel = 0; Real rateJ = 0;
  local_ydot[se] = 0;

  if(m_ppars.isConstMobDiff) goto enforcePositivity;
  
 energySourceTerm:
  if(m_constTe == 0)
    {
      // Evaluate Se, variable Te:


      if(nel >0 && elEnrgy > 0e0)
	// Inelastic collision energy rate
	{
	  Real ecxs = 1.0;
	  for (int r = 0; r < nReactions(); r++)  
	    {
	      if (eV > abs(epar[r])) 
		local_ydot[se] -= rates[r] * epar[r];
	    }
	  local_ydot[se] *= m_ppars.q * m_ppars.NA;
	  rateInel = local_ydot[se];

	  // Elastic collision energy rate
	  /*       
		   Real exsbg = 1.0;  
		   Real exso2 = 1.0;  
		   for (int l=0; l<3; l++)
		   {
		   exsbg *= exp(sigebg[l*3]*exp(-pow((log(een)-sigebg[3*l+1])/sigebg[3*l+2],2)));  
		   exso2 *= exp(sigeo2[l*3]*exp(-pow((log(een)-sigeo2[3*l+1])/sigeo2[3*l+2],2)));  
		   }
	  */
  	  // MDCancel Helium only
	  Real lte = log10(tel);
	  Real exsbg = sigebg[0];
	  for (int j = 1; j <8; j++)  
	    exsbg = exsbg + sigebg[j] * pow(lte,j);
	  exsbg = pow(10.0,exsbg);
	
	  Real ve = sqrt(2.0 * m_ppars.kb * m_ppars.NA * abs(tel) / mpar[m_species[m_eindex]]) ;
	  buff = 3.0 * m_ppars.kb * nel * pow(m_ppars.NA, 2.0) * mpar[m_species[m_eindex]] * (tel-teg) * ve
	    * (  Ng*exsbg / mpar[m_species[m_bgindex]] );
	  //MDCancel for air * (  nbg*exsbg / mpar[m_species[m_bgindex]] +  no2*exso2 / mpar[m_species[m_o2index]] );
	  local_ydot[se] -= buff;
	}

      // Joules for electrons
      local_ydot[se] += m_odeSource*nel;
  
      
    }// end if m_constTe = 0

 enforcePositivity:
  for (int s = 0; s < numY; s++)
    {
      if(local_y[s] < 0) local_ydot[s] = abs(local_ydot[s]);
    }

  //if(m_print) pout () << "PlasmaPhysics::cvode_f, losses: odesrc=" << m_odeSource << ",elastic/el " << -buff/nel  << ",inel/el=" << rateInel /nel << " ,tel=" << tel << ", nel=" << nel << endl;
  /*  commented out because the check is expensive
  bool stopit = false;
  for (int r = 0; r < nReactions(); r++)stopit = stopit || std::isnan(rates[r]);
  //for(int k=0; k< cvode_NEQ; k++) stopit = stopit || (abs(local_y[k]) > 1e2);
  if (stopit || m_print)
    {
      pout() << t<< ", E " << m_odeEfield<< ", Tel " << tel;for(int k=0; k< cvode_NEQ; k++) pout() << ",  " << local_y[k]; pout() << endl;
      //pout() << t;for(int k=0; k< cvode_NEQ; k++) pout() << "; " << local_ydot[k]; pout() << endl;
      //pout() << t;for(int k=0;k<nReactions();k++) pout() << ":  " << rates[k]; pout() << endl;
      //pout() << "PlasmaPhysics::cvode_f, odeEfield= " << m_odeEfield << ", Tel=  "<< tel<< endl;
      //if(stopit) pout() << "Stopping for Vof " << m_vof << endl;
      //if(stopit) MayDay::Error("PlasmaPhysics::NAN in rate");
    }
  */
   
  return(0);
}

/*****************/

void PlasmaPhysics::zeroFarfield(LevelData<EBCellFAB>& 	  a_NData,
				 const Box&                   a_chemBox) 

// unused for now, see EBAMRSpecies plot

{

  for (DataIterator dit = a_NData.dataIterator(); dit.ok(); ++dit)
    {

      EBCellFAB& NData = a_NData[dit()];

      Box box = a_NData.box(dit()); //a_NData.getRegion(); //to include ghost cells
      const EBISBox& ebisBox = NData.getEBISBox();
      IntVectSet ivs(box);
      for (VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
	{
	  VolIndex vof = vofit();
	  if(!a_chemBox.contains(vof.gridIndex()))
	    {
	      for(int k=0; k< cvode_NEQ; k++)
		{
		  if(k != m_bgindex && k != m_o2index) NData(vof,k) = m_minData[k];
		  if(k == m_bgindex) NData(vof,k) = m_nBg;
		  if(k == m_o2index) NData(vof,k) = m_nO2;
		}
	    }
	}
    }
}

/*****************/

void PlasmaPhysics::resetNeutrals(LevelData<EBCellFAB>& a_conState)
{
  EBLevelDataOps::setVal(a_conState, m_nBg, m_bgindex);
  //MDCancel un-comment for air
  //EBLevelDataOps::setVal(a_conState, m_nO2, m_o2index);
}

void PlasmaPhysics::setConstantTemperature(LevelData<EBCellFAB>& a_Temp, const LevelData<EBCellFAB>&   a_cellNData) const
{
  EBLevelDataOps::setVal(a_Temp, m_Tg, m_tgComp);
  Real factor = m_Te/m_tfac;
  EBLevelDataOps::setVal(a_Temp, factor, m_eeComp);
  EBLevelDataOps::scale (a_Temp, a_cellNData, m_eeComp, m_eindex,1);
}

void PlasmaPhysics::zeroDomainBoundaries(LevelData<EBCellFAB>& 	  a_NData,
					 const ProblemDomain &    a_domain)

{
  for (DataIterator dit = a_NData.dataIterator(); dit.ok(); ++dit)
    {

      EBCellFAB& NData = a_NData[dit()];

      Box agrid = a_NData.box(dit()); //a_NData.getRegion(); //to include ghost cells
      const EBISBox& ebisBox = NData.getEBISBox();
      
      for (int idir = 0; idir < SpaceDim; idir++)
        {
	  for (SideIterator sit; sit.ok(); ++sit)
	    {
	      IntVect ivSideGrid =   agrid.sideEnd(sit());
	      IntVect ivSideDom  = a_domain.domainBox().sideEnd(sit());
	      //check to see if we actually are at domain boundary
	      if (ivSideGrid[idir] == ivSideDom[idir])
		{
		  //create box and ivs adjacent to domain boundary
		  //and interior to domain
		  Box sideBox = adjCellBox(agrid, idir, sit(), 1);
		  int ishift = -sign(sit());
		  sideBox.shift(idir, ishift);
		  IntVectSet ivs(sideBox);
		  for (VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
		    {
		      VolIndex vof = vofit();
		      for(int k=0; k< cvode_NEQ; k++)
			{
			  if(k != m_bgindex && k != m_o2index) NData(vof,k) = m_minData[k];
			  if(k == m_bgindex) NData(vof,k) = m_nBg;
			  if(k == m_o2index) NData(vof,k) = m_nO2;
			}
		    }
		}
	    }
	}
    }
}


/*****************/

void PlasmaPhysics::setMask(LevelData<EBCellFAB>&         a_mask,
			    LevelData<EBCellFAB>&   	  a_src,
			    const LevelData<EBCellFAB>&   a_cellNData,
			    const DisjointBoxLayout&      a_grids,
			    const EBISLayout&             a_ebisl) //now unused, in the future for multivals
											
{
  
  int sden = findSpec(string("e-"));
  int snrg = findSpec(string("ee"));
  
  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {

      //EBs
      EBCellFAB&       ebM = a_mask[dit()];
      EBCellFAB&       ebS = a_src[dit()];
      const EBCellFAB& ebN = a_cellNData[dit()];

      //Basefabs
      BaseFab<Real>& regM = ebM.getSingleValuedFAB();
      BaseFab<Real>& regS = ebS.getSingleValuedFAB();
      const BaseFab<Real>& regN = ebN.getSingleValuedFAB();

      Box cellBox = ebN.getRegion() & ebM.getRegion(); //a_grids.get(dit());
      
      ebM.setVal(0e0);

      FORT_MASK(CHF_FRA1(regM, 0),
		CHF_FRA(regS),
		CHF_CONST_FRA1(regN,sden),
		CHF_CONST_FRA1(regN,snrg),
		CHF_BOX(cellBox));
      
    }

}

/*****************/

//MDCancel direct discharge 2 electrodes
bool PlasmaPhysics::isElectrode(const VolIndex& a_vof, const Real& a_dx, const RealVect& a_centroid)
{
  return true;
}

// LM a routine to add only the leading part of the joule heating to the source term
void PlasmaPhysics::floor(LevelData<EBCellFAB>&   a_data) 		
{
  EBLevelDataOps::setCoveredVal(a_data,m_bgindex,m_nBg);
  for (DataIterator dit=a_data.dataIterator();dit.ok();++dit)
    {
      EBCellFAB& dataEBFAB = a_data[dit()];
      const Box& region = dataEBFAB.getRegion();
      const IntVectSet ivsBox(region);
      const EBISBox& ebisBox = dataEBFAB.getEBISBox();
      for (VoFIterator vofit(ivsBox, ebisBox.getEBGraph());vofit.ok(); ++vofit)
	{
	  const VolIndex& vof = vofit();
	  
	  for(int s=0;s< nComponents();s++)
	    {
	      Real maxVal = m_maxData[s];
	      Real minVal = m_minData[s];
	      Real& val = dataEBFAB(vof,s);
	      val = Max(minVal,val);
	      val = Min(maxVal,val);
	    }
	}
    }
}

//New from LM (added 08/15/13)

void PlasmaPhysics::floor(Vector<LevelData<EBCellFAB>* >& a_data) 		
{
  int numLevels = a_data.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {     
      for(int s=0;s< nComponents();s++)
	{
	  // the BG cannot be zero othrws Fortran produces NANs
	  if(s == m_bgindex)
	    EBLevelDataOps::setCoveredVal(*a_data[ilev],s,m_nBg);
	  else
	    EBLevelDataOps::setCoveredVal(*a_data[ilev],s,m_minData[s]);
	}
      for (DataIterator dit=a_data[ilev]->dataIterator();dit.ok();++dit)
	{
	  EBCellFAB& dataEBFAB =  (*a_data[ilev])[dit()];
	  const Box& region = dataEBFAB.getRegion();
	  const IntVectSet ivsBox(region);
	  const EBISBox& ebisBox = dataEBFAB.getEBISBox();
	  for (VoFIterator vofit(ivsBox, ebisBox.getEBGraph());vofit.ok(); ++vofit)
	    {
	      const VolIndex& vof = vofit();
	      
	      for(int s=0;s< nComponents();s++)
		{
		  Real maxVal = m_maxData[s];
		  Real minVal = m_minData[s];
		  Real& val = dataEBFAB(vof,s);
		  val = Max(minVal,val);
		  val = Min(maxVal,val);
		}
	    }
	}
    }
}

Real  PlasmaPhysics::regridThresh(int a_spec)
  { 

    if(!m_threshSet)
      {
	ParmParse pp;
	m_cationsThresh = -1;m_electronsThresh=-1;m_anionsThresh=-1;
	pp.query ("electron_refine_thresh",m_electronsThresh);
	pp.query ("cations_refine_thresh",m_cationsThresh);
	pp.query ("anions_refine_thresh",m_anionsThresh);
	m_threshSet=true;
	pout() << "Setting RefineThresh to " << m_electronsThresh << ", " << m_cationsThresh<< ", " << m_anionsThresh << endl;
      }

    Real retval = -100;
    if(m_ppars.isDimensionless)
      {
	if(a_spec >= nSpecies())
	  return m_ppars.charToSpecThrFac*m_electronsThresh;  //Charge Threshold
	else if(m_ppars.isNeutral[a_spec])
	  return 1e99;
	else
	  return m_electronsThresh;
      }
    
    int ispec = min(a_spec,nSpecies()-1);
    unsigned found = m_charSpec[ispec].find_last_of("+-");
    if(found <= m_charSpec[ispec].length())
      {
	string lastString = m_charSpec[ispec].substr(found,1);
	if(lastString.compare("+") == 0)
	  retval=m_cationsThresh;
	else if (ispec == m_eindex)
	  retval=m_electronsThresh;
	else if (lastString.compare("-") == 0)
	  retval=m_anionsThresh;
      }
    if(a_spec >= nSpecies()) retval *= m_ppars.charToSpecThrFac;

    return retval;
  }
#include "NamespaceFooter.H"
