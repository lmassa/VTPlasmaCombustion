#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _PLASMAPHYSICS_H_
#define _PLASMAPHYSICS_H_

#include "IntVect.H"
#include "LevelData.H"
#include "RefCountedPtr.H"
#include "BaseIVFactory.H"
#include "EBCellFAB.H"
#include "EBCellFactory.H"
#include "EBFluxFAB.H"
#include "EBFluxFactory.H"
#include "EBISLayout.H"

#include "gasModel.H"


// LM Stiff ODE solver
#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <cvode/cvode_direct.h>        /* access to CVDls interface            */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */

#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ */


#include "NamespaceHeader.H"




/// A class containing the Plasma physics parameters and routines
/**
 */

//**
class PlasmaPhysics
{
public:
  /// Constr
  PlasmaPhysics();

  ~PlasmaPhysics();

  void define();


  void initialize(LevelData<EBCellFAB>& a_conState,
		  const EBISLayout&     a_ebisl,
		  const Real&           a_dx,
		  const RealVect&       a_domainLength);

  void setEBy(const Real& a_yEB);

  void setEps(	const Real& a_epsg,
				const Real& a_epsd);
						
  Real BCvalue(int&           a_ispec,
	       const int&     a_idir,
	       const Side::LoHiSide& a_side=Side::Lo);

  void plasmaSources(LevelData<EBCellFAB>&   		a_src,
                     const LevelData<EBCellFAB>& 	a_NData,
                     const LevelData<EBCellFAB>& 	a_teData,
                     const LevelData<EBCellFAB>& 	a_fGphiData,
                     const LevelData<EBCellFAB>& 	a_GphiData,
                     const EBISLayout&			    a_ebisl);
																
  void plasmaSourceTerms(EBCellFAB&   	a_src,
			 const EBCellFAB&			a_rr,
			 const EBCellFAB&			a_sen);
																

  void plasmaSourceEnergy(EBCellFAB&		a_sen,
                          const EBCellFAB&	a_rr,
                          const EBCellFAB&	a_NData,
			  			const EBCellFAB& 	a_teData,
                          const EBCellFAB& 	a_fGphiData);


  void fluidSources(LevelData<EBCellFAB>&   a_src,
		    const LevelData<EBCellFAB>& 	a_NData,
		    const LevelData<EBCellFAB>& 	a_teData,
		    const LevelData<EBCellFAB>& 	a_GphiData,
		    const LevelData<EBCellFAB>& 	a_Velo,
		    const EBISLayout&				a_ebisl);
  
															
  void fluidSourceTerms(EBCellFAB&	        a_src,
			const EBCellFAB&	a_rr,
			const EBCellFAB&	a_NData,
			const EBCellFAB& 	a_teData,
			const EBCellFAB& 	a_GphiData,
			const EBCellFAB& 	a_Velo);
  
  void reactionRates(LevelData<EBCellFAB>&  			a_rr,
		     const LevelData<EBCellFAB>& 		a_NData,
		     const LevelData<EBCellFAB>& 		a_teData,
			 const LevelData<EBCellFAB>& 	a_GphiData);
  
  void reactionRates(	EBCellFAB& 			a_rr, 
						const EBCellFAB&	a_NData,
						const EBCellFAB& 	a_teData,
						const EBCellFAB& 	a_GphiData);

  void ElectricPotentialRHS( 	LevelData<EBCellFAB>&   			a_rr,
						        const LevelData<EBCellFAB>& 		a_NData,
								const Real&							a_dx);

  void transportCoefficients(LevelData<EBFluxFAB>&              a_fluxMobData,
			     LevelData<EBFluxFAB>&      	a_fluxDiffData,
			     const LevelData<EBCellFAB>&   	a_cellNData,
			     const LevelData<EBCellFAB>&   	a_cellTData,
			     const LevelData<EBCellFAB>&        a_cellphiData,
			     const LevelData<EBFluxFAB>&        a_GphiData,
			     const DisjointBoxLayout&      	a_grids,
			     const EBISLayout&             	a_ebisl,
			     const ProblemDomain&          	a_domain,
			     const Real&			a_dx);

  void printCoefficients(LevelData<EBFluxFAB>&     a_fluxMobData,
			 const DisjointBoxLayout&  a_grids,
			 const int&                a_ispec,
			 const int&                a_flag);

  //IV version
  void transportCoefficients(	LevelData< BaseIVFAB<Real> >&     a_fluxMobData,
				LevelData< BaseIVFAB<Real> >&     a_fluxDiffData,
				LevelData< BaseIVFAB<Real> >&     a_fluxWallData,
				LevelData< BaseIVFAB<Real> >&     a_rhosRHS,
				const LevelData<EBCellFAB>&       a_cellNData,
				const LevelData<EBCellFAB>&       a_cellTData,
				const LevelData<BaseIVFAB<Real>>& a_IVGPhiData,
				const LevelData<EBCellFAB>&       a_cellGPhiData,
				const DisjointBoxLayout&          a_grids,
				const EBISLayout&                 a_ebisl,
				const ProblemDomain&              a_domain,
				const Real&                       a_dx);
  
  
  void rhos(LevelData< BaseIVFAB<Real> >&   		a_rhos,
	    const LevelData< BaseIVFAB<Real> >&   	a_fluxWallData,
	    const DisjointBoxLayout&   			a_grids,
	    const EBISLayout&          			a_ebisl,
	    const ProblemDomain&       			a_domain,
	    const int&          			a_ghost,
	    const Real&          			a_dt);
  
  static int nSpecies();
  static int nReactions();
  plasmaParameters plasmaPars();

  // LM121812 luca's routines
  int getElectronEnergyComponent() const
  {
    return (m_eindex + 1);
  }
  // LM121812 luca's routines
  int getElectronComponent() const
  {
    return (m_eindex);
  }

  Real setImplicitTimeStep() const
  {
    return m_gasModel->setImplicitTimeStep();
  }

  
  Real gasCharge(const EBCellFAB& a_specFAB,  const VolIndex& a_vof) const
  {
    		
    Real Q = 0;
    for (int s = 0; s <= m_eindex; s++) Q += m_ppars.chargeNr[m_species[s]] * a_specFAB(a_vof, s);
    if(!m_ppars.isDimensionless)Q*=(m_ppars.q * m_ppars.NA);
    return Q;
  }
  
  Real gasCharge(const BaseIVFAB<Real>& a_specFAB,  const VolIndex& a_vof) const
  {
    		
    Real Q = 0;
    for (int s = 0; s <= m_eindex; s++) Q += m_ppars.chargeNr[m_species[s]] * a_specFAB(a_vof, s);
    if(!m_ppars.isDimensionless)Q*=(m_ppars.q * m_ppars.NA);
    return Q;
  }

  
  Vector<int>& getWallBCinfo()
  {
    return m_wallBCInfo;
  }
  
  void setWallBCinfo()
  {
    m_wallBCInfo.resize(nComponents(),-1); //includes the EE
    m_gasModel->setWallBCinfo(m_wallBCInfo);
  }

  
  int nComponents() const
  {
    return (m_species.size()+1);
  }

  static Real convertEe2Te (const Real& ee, const Real& ne, const int flag = 0)
  {
    return (flag == 0)? Min(ee/ne*m_tfac,m_Temax) : m_Te;
  }
  
  static Real convertTe2Een (const Real& Te) //een is per unit nel
  {
    return Te / m_tfac / m_ppars.NA;
  }
  
  bool isDirichlet(int a_idir,const Side::LoHiSide& a_side)
  {
	return (a_idir == 1);
	// return (a_idir == 1 && a_side == Side::Lo);
	//return (true);
	//return (a_idir == 1 || sign(a_side)>0);
	//return (a_idir == 1 && sign(a_side)>0);
  }

  bool isConstantTe(){ return (m_constTe == 1);}
  bool isConstantBg(){ return (m_constBg == 1);}
  //  LM's routine (not used preparing for future)
  Real m_cationsThresh,m_electronsThresh,m_anionsThresh;
  bool m_threshSet;
  Real regridThresh(int a_spec);

  Real applyVoFStencilComp(const VoFStencil& a_sten, const EBCellFAB& a_fab, const int& a_comp);

  RealVect applyVoFStencilVect(const VoFStencil& a_sten, const EBCellFAB& a_fab);



  //LM routine for Joule heating
  void addJouleHeat(LevelData<EBCellFAB>&         a_cellSData,
		    const LevelData<EBCellFAB>&   a_cellNData,
		    const LevelData<EBCellFAB>&   a_cellTData,
		    const LevelData<EBCellFAB>&   a_cellEData,
		    const LevelData<EBCellFAB>&   a_cellGNata,
		    const LevelData<EBCellFAB>&   a_cellVelData,
		    const DisjointBoxLayout&      a_grids,
		    const EBISLayout&             a_ebisl);
  //LM routine for max DT
  void maxDT(LevelData<EBCellFAB>&         a_cellSData,
	     const LevelData<EBCellFAB>&   a_cellNData,
	     const LevelData<EBCellFAB>&   a_cellTData,
	     const Real&                   a_dt=1.0);

  void AddJouleHeatF(LevelData<EBCellFAB>&         a_cellSData,
			       const LevelData<EBCellFAB>&   a_cellNData,
			       const LevelData<EBCellFAB>&   a_cellTData,
			       const LevelData<EBCellFAB>&   a_cellEData,
			       const LevelData<EBCellFAB>&   a_cellGNData,
			       const DisjointBoxLayout&      a_grids,
				  const EBISLayout&             a_ebisl);

  void electronTemperature(LevelData<EBCellFAB>&         a_cellTData,
			   const LevelData<EBCellFAB>&   a_cellNData,
			   const DisjointBoxLayout&      	a_grids,
			   const EBISLayout&             	a_ebisl);

  void floor(LevelData<EBCellFAB>&   a_data);

//added from LM 08/15/13
  void floor(Vector<LevelData<EBCellFAB>* >& a_data);

  // LM templates
  template <typename T> static int sgn(T val) 
  {
    return (T(0) < val) - (val < T(0));
  }
  template<typename T, size_t N> T * end(T (&ra)[N]) 
  {
    return ra + N;
  }

  //LM species
  Vector<string> m_charSpec;
  int findSpec(const string& a_spec);
  
  //LM cvode stuff

  realtype cvode_reltol;
  N_Vector cvode_y, cvode_abstol;
  SUNMatrix cvode_A;
  SUNLinearSolver cvode_LS;

  void *cvode_mem;

  void plasmaSourcesODE(LevelData<EBCellFAB>& 	     a_NData,
			const LevelData<EBCellFAB>&  a_teData,
			const LevelData<EBCellFAB>&  a_fGphiData,
			const LevelData<EBCellFAB>&  a_GphiData,
			const EBISLayout&	     a_ebisl,
			const Real&                  a_dt,
			const int&                   a_level=0);
  void zeroFarfield(LevelData<EBCellFAB>& 	  a_NData,
		    const Box&                   a_chemBox);
  void zeroDomainBoundaries(LevelData<EBCellFAB>& 	  a_NData,
			    const ProblemDomain &    a_domain);
  void resetNeutrals(LevelData<EBCellFAB>& a_conState);

//update LM 08/15/13
  static int cvode_f(realtype t, N_Vector y, N_Vector ydot, void *user_data=NULL);


  void setMask(LevelData<EBCellFAB>&         a_mask,
	       LevelData<EBCellFAB>&   		a_src,
	       const LevelData<EBCellFAB>&   a_cellNData,
	       const DisjointBoxLayout&      	a_grids,
	       const EBISLayout&             	a_ebisl);

  void setConstantTemperature( LevelData<EBCellFAB>& a_Temp, const LevelData<EBCellFAB>&   a_cellNData) const;

  
  //end of cvode stuff

  //LM unprotected variables
  Real  m_molWeight;
  bool m_useThermalVelocity;
  //LM electrode location::
  Vector<Real> m_EbDxInt;

  RefCountedPtr<gasModel> m_gasModel;
  Vector<int> m_wallBCInfo;

  //static variables
  static int cvode_NEQ;
  static plasmaParameters m_ppars;
  static Vector<Vector<Real> > m_rrates;
  static Vector<Vector<Real> > m_transp;
  static Vector<Vector<Real> > m_reactant;
  static Vector<Vector<Real> > m_netReac;
  static Vector<Real> m_actReac;
  static Vector<Real> m_esecemcoef;
  static Vector<Real> m_esecemen;
  static Vector<int> m_species;
  static Vector<int> m_reactions;
  static int m_tgComp;
  static int m_eeComp;
  static int m_bgindex;
  static int m_o2index;
  static int m_eindex;
  static Real m_GasPressure;
  static Real m_GasTemperature;
  static Real m_odeSource;
  static Real m_odeEfield;
  static Real m_nBg;
  static Real m_Te;
  static Real m_Temax;
  static Real m_Tg;
  static Vector<Real> m_minData;
  static Vector<Real> m_maxData;
  static bool m_print;
  static bool m_DC;
  static int m_constTe;
  static int m_constBg;
  static string m_bgSpec;
  static VolIndex m_vof;
  static Real m_tfac;

protected:
   
  int m_yEB;
  Real m_epsg;
  Real m_epsd;
  bool m_isDefined;
  Real m_nO2;//m_nBg

  Real m_EbDirBcValue;
  bool isElectrode(const VolIndex& a_vof, const Real& a_dx, const RealVect& a_centroid);

//added from LM 08/15/13
  //cvode vars
  bool m_odeflag;
private:

  //disallowed 
  void operator=(const PlasmaPhysics& a_input)
  {
    MayDay::Error("invalid operator");
  }
  PlasmaPhysics(const PlasmaPhysics& a_input)
  {
    MayDay::Error("invalid operator");
  }
 
};
#include "NamespaceFooter.H"
#endif
