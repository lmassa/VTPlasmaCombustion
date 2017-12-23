#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBPATCHSPECIES_H_
#define _EBPATCHSPECIES_H_
#include "EBPatchGodunov.H"
#include "EBSpeciesMacros.H"
#include "EBSpeciesIBC.H"
#include "Stencils.H"
#include "CoveredFaceFAB.H"
#include "BaseIFFAB.H"

#define TOL 1.e-12
///
/**
   An EBPatchAdvect for species gas dynamics.
*/
class EBPatchSpecies : public EBPatchGodunov
{
public: 
  ///
  /**
   */
  EBPatchSpecies();

  ///
  /**
   */
  virtual ~EBPatchSpecies();
  /*****************************/
  void setVelocities(const EBCellFAB& a_normalVel,
		     const EBFluxFAB& a_advectionVel); 
  void setVelocities(const EBCellFAB& a_normalVel); 

  void setSpeciesNumber(int a_nSpec);
  int getSpeciesNumber();

/*****************************/
void setValidBox(const Box&        a_validBox,
            const EBISBox&    a_ebisBox,
            const IntVectSet& a_coarseFineIVS,
            const Real&       a_cumulativeTime,
		 const Real&       a_timeStep);
/******/
void slope(EBCellFAB&       a_slopePrim,
	   EBCellFAB&       a_slopeNLim,
	   const EBCellFAB& a_primState,
	   const EBCellFAB& a_flattening,
	   const int&       a_dir,
	   const Box&       a_box,
	   bool a_doAgg);
/******/
  void doSecondOrderSlopes(EBCellFAB&       a_delta2W,
			   EBCellFAB&       a_deltaWL,
			   EBCellFAB&       a_deltaWR,
			   EBCellFAB&       a_deltaWC,
			   const EBCellFAB& a_primState,
			   const int&       a_dir,
			   const Box&       a_box);

/******/
void extrapToCoveredFaces(BaseIVFAB<Real>&        a_extendedPrim,
			  const EBFaceFAB&        a_primFace,
			  const EBCellFAB&        a_primState,
			  const Vector<VolIndex>& a_coveredFaces,
			  const int&              a_faceDir,
			  const Side::LoHiSide&   a_sd,
			  const Box&              a_box);
/********/
void pointExtrapToCovered2D(Vector<Real>&           a_extrapVal,
			    const EBFaceFAB&        a_primFace,
			    const EBCellFAB&        a_primState,
			    const int&              a_faceDir,
			    const VolIndex&         a_vof,
			    const RealVect&         a_normal,
			    const Side::LoHiSide&   a_sd,
			    const int&              a_numPrim);
/********/
void pointExtrapToCovered3D(Vector<Real>&           a_extrapVal,
			    const EBFaceFAB&        a_primFace,
			    const EBCellFAB&        a_primState,
			    const int&              a_faceDir,
			    const VolIndex&         a_vof,
			    const RealVect&         a_normal,
			    const Side::LoHiSide&   a_sd,
			    const int&              a_numPrim);
/*****************************/
void normalPred(EBCellFAB&       a_rhoLo,
           EBCellFAB&       a_rhoHi,
           const EBCellFAB& a_rho,
           const EBCellFAB& a_dRho,
           const Real&      a_dtbydx,
           const int&       a_dir,
		const Box&       a_box);
/*******************/
void extrapolatePrim(EBFluxFAB&                       a_flux,
		     EBCellFAB&                       a_primState,
		     EBCellFAB                        a_slopePrim[SpaceDim],
		     EBCellFAB                        a_slopeNLim[SpaceDim],
		     Vector<BaseIVFAB<Real> * >&      a_coveredFluxMinu,
		     Vector<BaseIVFAB<Real> * >&      a_coveredFluxPlus,
		     const Vector<IntVectSet >&       a_coveredSetsMinu,
		     const Vector<IntVectSet >&       a_coveredSetsPlus,
		     const Vector<Vector<VolIndex> >& a_coveredFaceMinu,
		     const Vector<Vector<VolIndex> >& a_coveredFacePlus,
		     const EBCellFAB&                 a_flattening,
		     const EBCellFAB&                 a_consState,
		     const EBCellFAB&                 a_source,
		     const Box&                       a_box,
		     const DataIndex&                 a_dit,
		     bool                             a_verbose);
/***/
void riemann(EBFaceFAB&       a_flux,
	     const EBCellFAB& a_primLeft,
	     const EBCellFAB& a_primRigh,
	     const int&       a_faceDir,
	     const Box&       a_box);
/***/
///returns fluxes
/*****************************/
void riemann(BaseIVFAB<Real>&        a_coveredFlux,
	     const BaseIVFAB<Real>&  a_exteState,
	     const EBCellFAB&        a_primState,
	     const Vector<VolIndex>& a_vofset,
	     const int&              a_faceDir,
	     const Side::LoHiSide&   a_sd,
	     const Box&       a_box);
/******/
void getFlux(EBFaceFAB&       a_flux,
	     const EBFaceFAB& a_prim,
	     const int&       a_faceDir,
	     const Box&       a_box);

/*****************************/
void computeEBIrregFlux(BaseIVFAB<Real>&  a_ebIrregFlux,
			const EBCellFAB&  a_primState,
			const EBCellFAB   a_slopePrim[SpaceDim],
			const IntVectSet& a_irregIVS,
			const EBCellFAB&  a_source);
/*****************************/
void averageVelToCC(EBCellFAB&                        a_normalVel,
               const EBFluxFAB&                  a_advectionVel,
               const Vector<BaseIVFAB<Real> * >& a_coveredVeloLo,
               const Vector<BaseIVFAB<Real> * >& a_coveredVeloHi,
               const Vector<Vector<VolIndex> >&  a_coveredFaceLo,
               const Vector<Vector<VolIndex> >&  a_coveredFaceHi,
               const Box&                        a_box) const;


  Real getMaxWaveSpeed(const EBCellFAB& a_U,
		       const Box&       a_box);
  void floorConserved(EBCellFAB&  a_consState,  const Box&  a_box); 
  void floorConserved(BaseIVFAB<Real>& a_consState, const IntVectSet& a_ivsIrreg);
  void floorPrimitives(EBCellFAB&  a_primState, const Box&  a_box);  
  void floorPrimitives(BaseIVFAB<Real>& a_primState, const IntVectSet& a_ivsIrreg);
  

/******/
  virtual bool isDefined() const;
/******/
  bool usesFlattening() const;
/******/
  bool usesArtificialViscosity() const;
/******/
  bool usesFourthOrderSlopes() const;
/******/
void consToPrim(EBCellFAB&       a_primState,
           const EBCellFAB& a_consState,
           const Box&       a_box,
           int              a_logflag,
		bool             a_verbose);
/******/
void consToPrim(BaseIVFAB<Real>&  a_primState,
           const EBCellFAB&  a_consState,
		const IntVectSet& a_ivs);
/******/
  void consToPrim(BaseIVFAB<Real>&        a_primState,
		  const BaseIVFAB<Real>&  a_consState,
		  const IntVectSet&       a_ivs);
/******/
void primToCons(EBCellFAB&       a_consState,
		const EBCellFAB& a_primState,
		const Box&       a_box);
/******/
void primToCons(BaseIVFAB<Real>&       a_consState,
		const BaseIVFAB<Real>& a_primState,
		const IntVectSet&      a_ivs);
/*****************************/
  void useLimiting(bool a_limiting);
/*****************************/
/*****************************/
  Vector<string> stateNames();
/*****************************/
  Vector<string> primNames();
/*****************************/

  int numPrimitives() const;

  ///
  /**
     Returns number of components for flux variables.
  */
  int numFluxes() const;

///
/**
   Returns number of components for conserved variables.
*/
  int numConserved() const;

///
/**
 */
  int numSlopes() const;

  //inlined
  void setMaxMin(const Real& a_maxVal,
                 const Real& a_minVal)
  {
    m_isMaxMinSet = true;
    m_maxVal = a_maxVal;
    m_minVal = a_minVal;
  }
  EBCellFAB& getPrimState()
  {
    return m_primState;
  }

  //spurious stuff
 Interval velocityInterval() const
  {
    MayDay::Error("should not be called");
    Interval retval;
    return retval;
  }
  virtual int densityIndex() const
  {
    return 0;
  }
  int pressureIndex() const
  {
    return -1;
  }
  int bulkModulusIndex() const
  {
    MayDay::Error("should not be called");
    return -1;
  }
  Real artificialViscosityCoefficient() const
  {
    MayDay::Error("should not be called");
    return -1.0;
  }

  
  void setCoveredConsVals(EBCellFAB& a_consState);
  void getCoveredValuesCons(Vector<Real>& a_covValues);
  void getCoveredValuesPrim(Vector<Real>& a_covValues);

  /*virtual void regularUpdate(EBCellFAB&          a_consState,
		   EBFluxFAB&          a_flux,
		   BaseIVFAB<Real>&    a_ebIrregFlux,
		   BaseIVFAB<Real>&    a_nonConservativeDivergence,
		   const EBCellFAB&    a_flattening,
		   const EBCellFAB&    a_source,
		   const Box&          a_box,  
		   const IntVectSet&   a_ivsSmall, 
		   const DataIndex&    a_dit,
		   bool                a_verbose);*/
  
private:

  const EBCellFAB* m_normalVelPtr;
  const EBFluxFAB* m_advectionVelPtr;

  bool m_isVelSet;
  bool m_isMaxMinSet;
  Real m_maxVal;
  Real m_minVal;
  int m_nSpec;
  EBSpeciesIBC* m_specBc;
  //disallowed for all the usual reasons
  void operator=(const EBPatchSpecies& a_input)
  {
    MayDay::Error("invalid operator");
  }
  EBPatchSpecies(const EBPatchSpecies& a_input)
  {
    MayDay::Error("invalid operator");
  }

};
#endif
