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

#include "EBSpeciesIBC.H"
#include "EBISLayout.H"
#include "EBLoHiCenter.H"
#include "ParmParse.H"
#include "VoFIterator.H"
#include "EBLGIntegrator.H"
#include "EBPatchSpeciesF_F.H"
#include "EBArith.H"
/****************************/
/****************************/
EBSpeciesIBC::
EBSpeciesIBC(const Real&     a_nSpec, PlasmaPhysics*  a_PlasmaPhysics)
  :EBPhysIBC()
{
  m_nSpec = a_nSpec;
  m_PlasmaPhysics=a_PlasmaPhysics;
  m_isDefined = false;
  m_isVelSet = false;
}
/****************************/
void
EBSpeciesIBC::define(const ProblemDomain&  a_domain,
		     const RealVect& a_dx)
{
  m_domain = a_domain;
  m_dx = a_dx[0];
  CH_assert(Abs(a_dx[0]-a_dx[1]) < 1.e-9);


  //need to know more stuff
  
  ParmParse pp;
  vector<Real> length(SpaceDim);
  pp.getarr("domain_length",length,0,SpaceDim);
  for (int idir=0;idir<SpaceDim;idir++) m_domainLength[idir] = length[idir];
  

  m_isDefined = true;
}
/****************************/
void
EBSpeciesIBC::
setVelocities(const EBCellFAB& a_normalVel,
              const EBFluxFAB& a_advectionVel)

{
  m_normalVelPtr    = &a_normalVel;
  m_advectionVelPtr = &a_advectionVel;
  m_isVelSet = true;
}
/****************************/
void
EBSpeciesIBC::
setBndrySlopes(EBCellFAB&       a_deltaPrim,
               const EBCellFAB& a_primState,
               const EBISBox&   a_ebisBox,
               const Box&       a_box,
               const int&       a_dir)
{
  // don't want to deal with the periodic case
  // CH_assert(!m_domain.isPeriodic(a_dir));  // luca's comment
  if (m_domain.isPeriodic(a_dir)) return;
 CH_assert(m_isDefined);

  Box loBox,hiBox,centerBox,domain;
  int hasLo,hasHi;
  Box slopeBox = a_deltaPrim.getRegion()&m_domain;
  int numSlop = a_deltaPrim.nComp();
  // Generate the domain boundary boxes, loBox and hiBox, if there are
  // domain boundarys there
  eblohicenter(loBox,hasLo,hiBox,hasHi,centerBox,domain,
             slopeBox,m_domain,a_dir);

  // Set the boundary slopes if necessary
  if ((hasLo != 0) || (hasHi != 0))
    {
      BaseFab<Real>& regDeltaPrim = a_deltaPrim.getSingleValuedFAB();
      const BaseFab<Real>& regPrimState = a_primState.getSingleValuedFAB();
      /**/
      FORT_SLOPEBCS(CHF_FRA(regDeltaPrim),
                    CHF_CONST_FRA(regPrimState),
                    CHF_CONST_REAL(m_dx),
                    CHF_CONST_INT(a_dir),
                    CHF_BOX(loBox),
                    CHF_CONST_INT(hasLo),
                    CHF_BOX(hiBox),
                    CHF_CONST_INT(hasHi));
      /**/
    }
  for (SideIterator sit; sit.ok(); ++sit)
    {
      bool doThisSide;
      Box thisBox;
      if (sit() == Side::Lo)
        {
          doThisSide = (hasLo != 0);
          thisBox = loBox;
        }
      else
        {
          doThisSide = (hasHi != 0);
          thisBox = hiBox;
        }
      if (doThisSide)
        {

          // the cells for which the regular calculation
          //are incorrect are the grown set of the multivalued
          //cells intersected with the appropriate bndry box
          //and added to the irregular cells on the boundary.
          Box boxGrown = thisBox;
          boxGrown.grow(a_dir, 1);
          boxGrown &= m_domain;
          IntVectSet ivs = a_ebisBox.getMultiCells(boxGrown);
          ivs &= loBox;
          IntVectSet ivsIrreg = a_ebisBox.getIrregIVS(thisBox);
          ivs |= ivsIrreg;
          for (VoFIterator vofit(ivs, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
            {
              const VolIndex& vof = vofit();
              //all slopes at boundary get set to zero
              //except normal velocity.  just following the
              //fortran here
              int inormVelVar = a_dir + QVELX;
              for (int ivar = 0; ivar < numSlop; ivar++)
                {
                  if (ivar != inormVelVar)
                    {
                      a_deltaPrim(vof, ivar) = 0.0;
                    }
                }

              //for normal velocity
              //do strangely limited slope
              //just lifted from the fortran.
              Side::LoHiSide otherSide = flip(sit());
              Vector<FaceIndex> faces =
                a_ebisBox.getFaces(vof, a_dir, otherSide);

              Real thisVel = a_primState(vof, inormVelVar);
              Real otherVel = 0.;
              for (int iface = 0; iface < faces.size(); iface++)
                {
                  VolIndex otherVoF = faces[iface].getVoF(otherSide);
                  otherVel += a_primState(otherVoF,inormVelVar);
                }
              if (faces.size() > 0)
                {
                  otherVel /= Real(faces.size());
                  Real slope;
                  if (sit() == Side::Lo)
                    {
                      slope = otherVel - thisVel;
                    }
                  else
                    {
                      slope = thisVel  - otherVel;
                    }
                  //trick to make sure state will not change sign
                  if (slope*thisVel < 0)
                    {
                      a_deltaPrim(vof, inormVelVar) = 0.0;
                    }
                  else
                    { //slope and vel same sign
                      Real rsign = 1.0;
                      if (thisVel < 0.0)
                        rsign = -1.0;
                      //told you this was odd.
                      Real dvmin = Min(Abs(slope), Abs((Real)2.*thisVel));
                      a_deltaPrim(vof, inormVelVar) = rsign*dvmin;
                    }
                }
              else //no faces on high side of low boundary vof
                {
                  a_deltaPrim(vof, inormVelVar) = 0.0;
                }
            } //end loop over irregular vofs at boundary
        }
    } //end loop over boundary sides
}
/****************************/
/****************************/
void
EBSpeciesIBC::
fluxBC(EBFluxFAB&            a_flux,
       const EBCellFAB&      a_primCenter,
       const EBCellFAB&      a_primExtrap,
       const Side::LoHiSide& a_side,
       const Real&           a_time,
       const EBISBox&        a_ebisBox,
       const DataIndex&      a_dit,
       const Box&            a_box,
       const Box&            a_faceBox,
       const int&            a_dir)
{
 CH_assert(m_isDefined);
 CH_assert(m_isVelSet);

  Box FBox = a_flux[a_dir].getSingleValuedFAB().box();
  Box cellBox = FBox;
  int numFlux = a_flux[a_dir].nComp();

  // Determine which side and thus shifting directions
  int isign = sign(a_side);
  cellBox.shiftHalf(a_dir,isign);

  // Is there a domain boundary next to this grid
  if (!m_domain.contains(cellBox))
    {
      cellBox &= m_domain;
      // Find the strip of cells next to the domain boundary
      Box boundaryBox = bdryBox(cellBox, a_dir, a_side, 1);

      // Shift things to all line up correctly
      boundaryBox.shiftHalf(a_dir,-isign);
      BaseFab<Real>& regFlux = a_flux[a_dir].getSingleValuedFAB();
      const BaseFab<Real>& regPrimExtrap = a_primExtrap.getSingleValuedFAB();
      BaseFab<Real>& regAdvect = const_cast<BaseFab<Real>&> ((*m_advectionVelPtr)[a_dir].getSingleValuedFAB());

      // Set the boundary fluxes
      bool inflow = (a_dir ==0 && a_side == Side::Lo) ;
      if (inflow)
	{
          regFlux.shiftHalf(a_dir,-isign);
          regAdvect.shiftHalf(a_dir,-isign);
	  for(int iSpec =0;iSpec<m_nSpec; iSpec++)
	    {
	      Real value = m_PlasmaPhysics->BCvalue(iSpec,a_dir);
	      FORT_SCALARFABMULT(CHF_FRA1(regFlux,iSpec),CHF_CONST_FRA1(regAdvect,0),
				 CHF_CONST_REAL(value),CHF_BOX(boundaryBox));
	    }
          regFlux.shiftHalf(a_dir,isign);
          regAdvect.shiftHalf(a_dir,isign);
          IntVectSet ivs = a_ebisBox.getMultiCells(boundaryBox);
          for (VoFIterator vofit(ivs, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
            {
              const VolIndex& vof = vofit();
              Vector<FaceIndex> bndryFaces = a_ebisBox.getFaces(vof, a_dir, a_side);
              for (int iface= 0; iface < bndryFaces.size(); iface++)
                {
                  const FaceIndex& face = bndryFaces[iface];
		  const Real velface  = (*m_advectionVelPtr)[a_dir](face,0);
		  for(int iSpec =0;iSpec<m_nSpec; iSpec++)
		    {
		      Real value = m_PlasmaPhysics->BCvalue(iSpec,a_dir);
                      a_flux[a_dir](face, iSpec) = velface*value;
		    }	      
		}
	    }	    
	  
	}
      else
	{
          regFlux.shiftHalf(a_dir,-isign);
          regAdvect.shiftHalf(a_dir,-isign);
	  
	  
	  FORT_SPECIESEXTRAPBC(CHF_FRA(regFlux),
			       CHF_CONST_FRA(regPrimExtrap),
			       CHF_CONST_FRA1(regAdvect,0),
			       CHF_CONST_INT(m_nSpec),
			       CHF_BOX(boundaryBox));

          // Shift returned fluxes to be face centered
          regFlux.shiftHalf(a_dir,isign);
          regAdvect.shiftHalf(a_dir,isign);

          IntVectSet ivs = a_ebisBox.getMultiCells(boundaryBox);
          for (VoFIterator vofit(ivs, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
            {
              const VolIndex& vof = vofit();
              Vector<FaceIndex> bndryFaces = a_ebisBox.getFaces(vof, a_dir, a_side);
              for (int iface= 0; iface < bndryFaces.size(); iface++)
                {
                  const FaceIndex& face = bndryFaces[iface];
		  const Real velface  = (*m_advectionVelPtr)[a_dir](face,0);
		  for(int iSpec =0;iSpec<m_nSpec; iSpec++)
		    {
		      Real value = a_primExtrap(vof, iSpec);
                      a_flux[a_dir](face, iSpec) = velface*value;
		    }
                }
            }
        } //end ifinflow
    } //!domain(contains(cellbox)
}
/****************************/
/****************************/
/****************************/
void
EBSpeciesIBC::
initialize(LevelData<EBCellFAB>& a_conState, const EBISLayout& a_ebisl) const
{
 CH_assert(m_isDefined);

 
 // do this in plasmaPhysics
 m_PlasmaPhysics->initialize(a_conState, a_ebisl, m_dx, m_domainLength);
}
/****************************/
/****************************/
EBSpeciesIBC::
~EBSpeciesIBC()
{
}
/****************************/
/****************************/
