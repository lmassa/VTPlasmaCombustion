#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "MixedSpeciesFluxDomainBC.H"
#include "MixedSpeciesFluxDomainBCF_F.H"
#include "EBArithF_F.H"
#include "NamespaceHeader.H"

/*************/
void 
SpeciesFluxBaseDomainBC::setPhysics(PlasmaPhysics*  a_PlasmaPhysics, const int&  a_ispec)
{
  m_PlasmaPhysics = a_PlasmaPhysics;
  m_ispec = a_ispec;  
  m_onlyHomogeneous = false;
}
/*****/
/*****/
MixedSpeciesFluxDomainBC::
MixedSpeciesFluxDomainBC()
{
    m_coefSet  = false;
}
/*****/
MixedSpeciesFluxDomainBC::
~MixedSpeciesFluxDomainBC()
{
}


void
MixedSpeciesFluxDomainBC::
getFaceFlux(BaseFab<Real>&        a_faceFlux,
            const BaseFab<Real>&  a_phi,
            const RealVect&       a_probLo,
            const RealVect&       a_dx,
            const int&            a_idir,
            const Side::LoHiSide& a_side,
            const DataIndex&      a_dit,
            const Real&           a_time,
            const bool&           a_useHomogeneous)
{ 

  if (m_onlyHomogeneous && !a_useHomogeneous)
    {
      MayDay::Error("MixedPoissonDomainBC::getFaceFlux called with undefined inhomogeneous BC");
    }


  CH_TIME("MixedPoissonDomainBC::getFaceFlux");

  //only zero flux if Neuman
  if(m_PlasmaPhysics->isDirichlet(a_idir,a_side)== false)
    {
      a_faceFlux.setVal(0e0);
      return;
    }

  //LM:: note the Neumann Boundary condition are in reality zero diffusion flux
  const Box& faceBox = a_faceFlux.box();

  int iside = -sign(a_side);  //LO:: +1; HI::-1

  const Box& box = a_faceFlux.box();
  BaseFab<Real>&   regD = (*m_Dk)  [a_dit][a_idir].getSingleValuedFAB();
  BaseFab<Real>&   regW = (*m_vecW)[a_dit][a_idir].getSingleValuedFAB();
  //shift these fabs to be aligned with flux
  // the flux is cell centered on the valid side so on the low side i in box is [0] both for the flux and phi
  // on te HI face i is [nx-1].
  //D,W are shifted next;For D and W on LO i -> face i-1/2, so to the left of the cell, on HI i-> i+1/2
  regD.shiftHalf(a_idir, iside); //if LO => iface between i and i-1
  regW.shiftHalf(a_idir, iside); // if HI=> iface between i+1 and i

  Real value;
  if (a_useHomogeneous)
    {
      value = 0.0;    
    }
  else
    {
      value = m_PlasmaPhysics->BCvalue(m_ispec,a_idir, a_side);
    }
  // so in the fortran: LO: i-> phi_0, Flux_left, D,W_left; HI: i->phi_(nx-1) Flux_right, D,W_right
  FORT_SETMIXEDFACEFLUX(CHF_FRA1(a_faceFlux,0),
			CHF_CONST_FRA1(a_phi,0),
			CHF_CONST_FRA1(regD,0),
			CHF_CONST_FRA1(regW,0),
			CHF_CONST_REAL(value),
			CHF_CONST_REALVECT(a_dx),
			CHF_CONST_INT(a_idir),
			CHF_CONST_INT(iside),
			CHF_BOX(box));  
  regD.shiftHalf(a_idir, -iside);
  regW.shiftHalf(a_idir, -iside);

}

/*****/
void
MixedSpeciesFluxDomainBC::
getFaceFlux(Real&                 a_faceFlux,
            const VolIndex&       a_vof,
            const int&            a_comp,
            const EBCellFAB&      a_phi,
            const RealVect&       a_probLo,
            const RealVect&       a_dx,
            const int&            a_idir,
            const Side::LoHiSide& a_side,
            const DataIndex&      a_dit,
            const Real&           a_time,
            const bool&           a_useHomogeneous)
{
  
  // simply use zeroflux BC
  //only zero flux if Neuman
  if(m_PlasmaPhysics->isDirichlet(a_idir,a_side)== false)
    {
      a_faceFlux=0e0;
      return;
    }

  Real imposedVal;
  if (a_useHomogeneous)
      imposedVal = 0.0;    
  else
      imposedVal = m_PlasmaPhysics->BCvalue(m_ispec,a_idir,a_side);

  const EBISBox& ebisBox = a_phi.getEBISBox();
  const ProblemDomain& domainBox = ebisBox.getDomain();
  Vector<FaceIndex> faces = ebisBox.getFaces(a_vof,a_idir,a_side);
  int iside = -sign(a_side);  //LO:: +1; HI::-1


  Real Dave = 0;
  Real Wave = 0;
  Real areaTot = 0;
  for (int iface = 0; iface < faces.size(); iface++)
    {
      Real areaFrac = ebisBox.areaFrac(faces[iface]);
      Real DFace  = (*m_Dk)[a_dit][a_idir](faces[iface], 0);
      Dave += areaFrac*DFace;
      Real WFace  = (*m_vecW)[a_dit][a_idir](faces[iface], 0);
      Wave += areaFrac*WFace;
      areaTot += areaFrac;
    }
  if (areaTot > 1.0e-8)
    {
      Dave /= areaTot;
      Wave /= areaTot;
    }

  Real faceVal = 0.0;
  if (faces.size() > 0)
    {
      if (faces.size()==1)
        {
          IntVectSet cfivs;
          FaceStencil faceSten = EBArith::getInterpStencil(faces[0], cfivs, ebisBox, domainBox);
          for (int isten=0; isten < faceSten.size(); isten++)
            {
              const Real& weight = faceSten.weight(isten);
              const FaceIndex& face = faceSten.face(isten);

              Real thisFaceVal;
	      getInsideVal(thisFaceVal,face,a_comp,a_phi,a_idir,
                             a_side,a_dit,a_time,false);
              faceVal += thisFaceVal*weight;
            }
          faceVal *= ebisBox.areaFrac(faces[0]);
        }
      else
        {
          MayDay::Error("MixedSpeciesFluxDomainBC::getHigherOrderFaceFlux has multi-valued faces (or could be 0 faces)");
        }
    }
  if(iside>0) a_faceFlux = (Dave*faceVal - Wave*imposedVal)/a_dx[a_idir];
  if(iside<0) a_faceFlux = (Dave*imposedVal - Wave*faceVal)/a_dx[a_idir];
  
}
void
MixedSpeciesFluxDomainBC::
getFaceGradPhi(Real&                 a_faceFlux,
               const FaceIndex&      a_face,
               const int&            a_comp,
               const EBCellFAB&      a_phi,
               const RealVect&       a_probLo,
               const RealVect&       a_dx,
               const int&            a_idir,
               const Side::LoHiSide& a_side,
               const DataIndex&      a_dit,
               const Real&           a_time,
               const bool&           a_useAreaFrac,
               const RealVect&       a_centroid,
               const bool&           a_useHomogeneous)
{
  int iside = -sign(a_side);
  const Real ihdx = 2.0 / a_dx[a_idir];
  const EBISBox& ebisBox = a_phi.getEBISBox();

  Real value = 0;
  if(!a_useHomogeneous) value = m_PlasmaPhysics->BCvalue(m_ispec,a_idir,a_side);
  
  const VolIndex& vof = a_face.getVoF(flip(a_side));

 

  a_faceFlux = iside * ihdx * (a_phi(vof,a_comp) - value);

  if (a_useAreaFrac)
    {
      MayDay::Error("MixedSpeciesFluxDomainBC::getFaceFlux -- useAreaFrac=TRUE");
      a_faceFlux *= ebisBox.areaFrac(a_face);
    }
}
/*****/
void
MixedSpeciesFluxDomainBC::
getFaceVel(Real&                 a_faceFlux,
           const FaceIndex&      a_face,
           const EBFluxFAB&      a_vel,
           const RealVect&       a_probLo,
           const RealVect&       a_dx,
           const int&            a_idir,
           const int&            a_icomp,
           const Real&           a_time,
           const Side::LoHiSide& a_side,
           const bool&           a_doDivFreeOutflow)
{
  a_faceFlux =  a_vel[a_idir](a_face, a_icomp);
}
void
MixedSpeciesFluxDomainBC::
getInsideVal(Real&                 a_faceVal,
	     const FaceIndex&      a_face,
	     const int&            a_comp,
	     const EBCellFAB&      a_phi,
	     const int&            a_idir,
	     const Side::LoHiSide& a_side,
	     const DataIndex&      a_dit,
	     const Real&           a_time,
	     const bool&           a_useAreaFrac)
{
  int iside = -sign(a_side);
  const EBISBox& ebisBox = a_phi.getEBISBox();
  
  const VolIndex& vof = a_face.getVoF(flip(a_side));
 

  a_faceVal = a_phi(vof,a_comp);

  if (a_useAreaFrac)
    {
      MayDay::Error("MixedSpeciesFluxDomainBC::getFaceFlux -- useAreaFrac=TRUE");
      a_faceVal *= ebisBox.areaFrac(a_face);
    }
}






/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/******/
MixedSpeciesFluxDomainBCFactory::
MixedSpeciesFluxDomainBCFactory()
{
  m_onlyHomogeneous = true;
}

/******/
MixedSpeciesFluxDomainBCFactory::
~MixedSpeciesFluxDomainBCFactory()
{
}
/*************/
void MixedSpeciesFluxDomainBCFactory::
setPhysics(PlasmaPhysics*  a_PlasmaPhysics, const int&  a_ispec)
{
  m_PlasmaPhysics = a_PlasmaPhysics;
  m_ispec = a_ispec;
  m_onlyHomogeneous = false;
}
/******/
MixedSpeciesFluxDomainBC*
MixedSpeciesFluxDomainBCFactory::
create(const ProblemDomain& a_domain,
       const EBISLayout&    a_layout,
       const RealVect&      a_dx)
{
  MixedSpeciesFluxDomainBC* newBC = new MixedSpeciesFluxDomainBC();
  if (!m_onlyHomogeneous)
    {
      newBC->setPhysics(m_PlasmaPhysics, m_ispec);
    }
  return newBC;
}
#include "NamespaceFooter.H"
