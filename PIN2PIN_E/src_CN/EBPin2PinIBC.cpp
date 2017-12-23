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

#include "EBPin2PinIBC.H"
#include "EBISLayout.H"
#include "EBLoHiCenter.H"
#include "ParmParse.H"
#include "VoFIterator.H"
#include "EBLGIntegrator.H"
#include "EBPatchPolytropicF_F.H"
#include "EBArith.H"
/****************************/
/****************************/
EBPin2PinIBC::
EBPin2PinIBC(const Real&     a_gamma,
                 const Real&     a_ms,
                 const Real&     a_center,
                 const int&      a_shocknorm,
                 const bool&     a_shockbackward,
                 const bool&     a_doRZCoords,
		 const BoundaryLayer* const a_BLPtr)
  :EBPhysIBC()
{
  int ishockback = 0;
  if (a_shockbackward)
    {
      ishockback = 1;
    }
  /**/
  Real preshockpress = 1.0;
  ParmParse pp;
  if (pp.contains("preshockpress"))
    {
      pp.get("preshockpress", preshockpress);
    }
  Real preshockdense  = a_gamma;
  if (pp.contains("preshockdense"))
    {
      pp.get("preshockdense", preshockdense);
    }

  m_BLinit=true;pp.query("use_BL_initialization",   m_BLinit);

  FORT_SETPIN2PIN(CHF_CONST_REAL(a_gamma),
                      CHF_CONST_REAL(a_ms),
                      CHF_CONST_REAL(a_center),
                      CHF_CONST_REAL(preshockpress),
                      CHF_CONST_REAL(preshockdense),
                      CHF_CONST_INT(a_shocknorm),
                      CHF_CONST_INT(ishockback));
  /**/
  m_doRZCoords    = a_doRZCoords;
  m_gamma         = a_gamma;
  m_ms            = a_ms;
  m_center        = a_center;
  m_shocknorm     = a_shocknorm;
  m_shockbackward = a_shockbackward;
  m_BLPtr         = a_BLPtr;

  m_isFortranCommonSet = true;
  m_isDefined = false;
}
/****************************/
/****************************/
void
EBPin2PinIBC::
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
 CH_assert(m_isFortranCommonSet);

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
		      //dvmin = 0.0; // luca cancel 0.0
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
EBPin2PinIBC::
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
 CH_assert(m_isFortranCommonSet);
  // CH_assert(!m_domain.isPeriodic(a_dir));  // luca's comment
  if (m_domain.isPeriodic(a_dir)) return;

  Box FBox = a_flux[a_dir].getSingleValuedFAB().box();
  Box cellBox = FBox;
  int numFlux = a_flux[a_dir].nComp();

  EBCellFAB d_primExtrap(a_primExtrap.getEBISBox(),a_primExtrap.getRegion(),a_primExtrap.nComp());
  d_primExtrap.clone(a_primExtrap);
  d_primExtrap.copy(a_primExtrap.getRegion(),Interval(0,d_primExtrap.nComp()-1),d_primExtrap.getRegion(),a_primExtrap,Interval(0,a_primExtrap.nComp()-1));

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

      // Set the boundary fluxes
      bool inbackofshock =
        (!m_shockbackward && a_side == Side::Lo) ||
        ( m_shockbackward && a_side == Side::Hi);
      if (a_dir != m_shocknorm+1 || a_side == Side::Hi) // && inbackofshock // luca's change
        {
          regFlux.shiftHalf(a_dir,-isign);

          // Set the boundary fluxes
          /**/
	  if (a_dir == m_shocknorm && a_side == Side::Lo)
            {
	      const BaseFab<Real>& d_regPrimExtrap = d_primExtrap.getSingleValuedFAB();
	      const Vector<Vector<Vector<Real> > >& BLspline = m_BLPtr->exportCoeffs();
	      const Vector<Vector<Real> >& BLcons = m_BLPtr->exportCons();
	      const Vector<Real>& eta = m_BLPtr->exportEta();
	      Real x0 = m_BLPtr->getX0(); // note it assumes that theta is assigned at x=0
	      Real y0 = m_BLPtr->mp_boundLoc;
	      Real ReRoot2 = m_BLPtr->getReLRoot2();
	      IntVectSet ivs(boundaryBox);
	      for (VoFIterator vofit(ivs, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
		{
		  const VolIndex& vof = vofit();
		  Vector<FaceIndex> facesLo = a_ebisBox.getFaces(vof, a_dir, Side::Lo);
		  if(facesLo.size() <= 0) continue;
		  RealVect centroid = a_ebisBox.centroid(facesLo[0])*m_dx;
		  RealVect faceloc = EBArith::getFaceLocation(facesLo[0],m_dx*(RealVect::Unit),centroid);
		  Real xx = faceloc[0]- x0;
		  Real yy = faceloc[1]- y0;
		  Real ee = yy*ReRoot2/sqrt(faceloc[0] - x0);
		  // two dimensional assumption	
		  Vector<Real> cvars(CNUM,0.0);
		  Real Ke = 0;
		  Real Invrhoh = 1.0;
		  for (int i = 0; i < m_BLPtr->nvars(); i++)
		    {
		      int icomp =  m_BL2Patch[i];
		      Real& yeval = cvars[icomp];
		      FORT_SPLINE5(CHF_REAL(yeval),
				   CHF_CONST_REAL(ee),
				   CHF_CONST_VR(eta),
				   CHF_CONST_VR(BLcons[i]),
				   CHF_CONST_VR(BLspline[i][0]),
				   CHF_CONST_VR(BLspline[i][1]),
				   CHF_CONST_VR(BLspline[i][2]),
				   CHF_CONST_VR(BLspline[i][3]),
				   CHF_CONST_VR(BLspline[i][4]));
		      // messy because using conservative variables should use primitives
		      if (i == 0)
			Invrhoh = 0.5/yeval;
		      else if(i == 1)
			Ke = (yeval*yeval);
		      else if(i == 2)
			{
			  yeval = yeval/(ReRoot2*sqrt(xx));
			  Ke += (yeval*yeval);
			}
		      else if (i == 3)
			yeval += Ke *Invrhoh;
		      
		    }
		  int temp = 0;		  
		  Vector<Real> qgdnv(QNUM);
		  FORT_POINTCONS2PRM(CHF_VR(cvars),
				     CHF_VR(qgdnv),
				     CHF_CONST_INT(temp));
		  //int i = QVELY;
		  //if (m_dx < 2e-3 && yy < 1e-3) pout() << xx << ", " << yy << ", " <<d_primExtrap(vof, i) << ", " << qgdnv[i] << endl;
		  for (int i = 0; i < QNUM; i++)
		    {
		      if(i !=QPRES)  //|| yy > 5e-2
			d_primExtrap(vof, i) = qgdnv[i];
		    }		    
		  
		}
              FORT_DIRICHLETBC(CHF_FRA(regFlux),
                            CHF_CONST_FRA(d_regPrimExtrap),
                            CHF_CONST_REAL(m_dx),
                            CHF_CONST_INT(a_dir),
                            CHF_BOX(boundaryBox));
            }// a_dir == m_shocknorm && a_side == Side::Lo
	  else if (a_side == Side::Hi || a_dir==2) //a_side == Side::Hi || a_dir==2
	    FORT_OUTFLOW(CHF_FRA(regFlux), 
			 CHF_CONST_FRA(regPrimExtrap),
			 CHF_CONST_REAL(m_dx),
			 CHF_CONST_INT(a_dir),
			 CHF_BOX(boundaryBox));
          else if (m_doRZCoords)
            {
              FORT_EXTRAPBCRZ(CHF_FRA(regFlux),
                              CHF_CONST_FRA(regPrimExtrap),
                              CHF_CONST_REAL(m_dx),
                              CHF_CONST_INT(a_dir),
                              CHF_BOX(boundaryBox));
            }
          else
            {
              FORT_EXTRAPBC(CHF_FRA(regFlux),
                            CHF_CONST_FRA(regPrimExtrap),
                            CHF_CONST_REAL(m_dx),
                            CHF_CONST_INT(a_dir),
                            CHF_BOX(boundaryBox));
            }
          /**/

          // Shift returned fluxes to be face centered
          regFlux.shiftHalf(a_dir,isign);
          //now for the multivalued cells.  Since it is pointwise,
          //the regular calc is correct for all single-valued cells.
          IntVectSet ivs = a_ebisBox.getMultiCells(boundaryBox);
          for (VoFIterator vofit(ivs, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
            {
              const VolIndex& vof = vofit();
              Vector<Real> qgdnv(QNUM);
              Vector<Real> fluxv(FNUM);
              for (int ivar = 0; ivar < QNUM; ivar++)
                {
                  qgdnv[ivar] = d_primExtrap(vof, ivar);
                }
	      // luca cancel limit pressure
	      //qgdnv[QPRES] = Min(qgdnv[QPRES] , 1.50);
              Vector<FaceIndex> bndryFaces =
                a_ebisBox.getFaces(vof, a_dir, a_side);
              for (int iface= 0; iface < bndryFaces.size(); iface++)
                {
                  if (m_doRZCoords)
                    {
                      Real radius = Real(vof.gridIndex()[0]);
                      FORT_POINTGETFLUXRZ(CHF_VR(fluxv),
                                          CHF_VR(qgdnv),
                                          CHF_CONST_INT(a_dir),
                                          CHF_CONST_REAL(radius));
                    }
                  else
                    {
                      /**/
                      FORT_POINTGETFLUX(CHF_VR(fluxv),
                                        CHF_VR(qgdnv),
                                        CHF_CONST_INT(a_dir));
                    }
                  /**/
                  const FaceIndex& face = bndryFaces[iface];
                  for (int ivar = 0; ivar < FNUM; ivar++)
                    {
                      a_flux[a_dir](face, ivar) = fluxv[ivar];
                    }
                }
            }
        } // a_dir != m_shocknorm+1 || a_side == Side::Hi
      else
        {  // a_dir = 2 &&  a_side == Side::Lo
          regFlux.shiftHalf(a_dir,-isign);
          //solid wall bcs
          if (m_doRZCoords)
            {
              FORT_SOLIDBCRZ(CHF_FRA(regFlux),
                             CHF_CONST_FRA(regPrimExtrap),
                             CHF_CONST_INT(isign),
                             CHF_CONST_REAL(m_dx),
                             CHF_CONST_INT(a_dir),
                             CHF_BOX(boundaryBox));
            }
          else
            {
              FORT_SLIPWALLSOLIDBC(CHF_FRA(regFlux),
                                   CHF_CONST_FRA(regPrimExtrap),
                                   CHF_CONST_INT(isign),
                                   CHF_CONST_REAL(m_dx),
                                   CHF_CONST_INT(a_dir),
                                   CHF_BOX(boundaryBox));
            }
          // Shift returned fluxes to be face centered
          regFlux.shiftHalf(a_dir,isign);
          int inormMomVar = CMOMX + a_dir;
          int inormVelVar = QVELX + a_dir;
          //now for the multivalued cells.  Since it is pointwise,
          //the regular calc is correct for all single-valued cells.
          IntVectSet ivs = a_ebisBox.getMultiCells(boundaryBox);
          for (VoFIterator vofit(ivs, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
            {
              const VolIndex& vof = vofit();
              //set all fluxes except normal momentum  to zero.
              //set normal momemtum flux to sign(normal)*pressure
              Vector<FaceIndex> bndryFaces = a_ebisBox.getFaces(vof, a_dir, a_side);
              for (int iface= 0; iface < bndryFaces.size(); iface++)
                {
                  const FaceIndex& face = bndryFaces[iface];
                  //set all fluxes to zero then fix normal momentum
                  for (int ivar = 0; ivar < numFlux; ivar++)
                    {
                      a_flux[a_dir](face, ivar) = 0;
                    }
                  Real press = d_primExtrap(vof, QPRES);
                  Real dense = d_primExtrap(vof, QRHO);
                  Real unorm = d_primExtrap(vof, inormVelVar);
                  Real speed = sqrt(m_gamma*press/dense);

                  if (m_doRZCoords)
                    {
                      if (a_dir == 0)
                        {
                          a_flux[a_dir](face,  CPRES) = press;
                          //at axis, all non pressure fluxes are zero
                          if (a_side == Side::Hi)
                            {
                              a_flux[a_dir](face, inormMomVar) =
                                isign*dense*unorm*speed;
                            }
                        }
                      else
                        {
                          a_flux[a_dir](face,  CPRES) = press;
                          a_flux[a_dir](face, inormMomVar) =
                            isign*dense*unorm*speed;
                        }
                    }
                  else
                    {

                      a_flux[a_dir](face, inormMomVar) =
                        press + isign*dense*unorm*speed;
                    }
                }
            }
        }
    }
}
/****************************/
/****************************/
void
EBPin2PinIBC::define(const ProblemDomain&  a_domain,
                         const RealVect& a_dx)
{
  m_domain = a_domain;
  m_dx = a_dx[0];
  CH_assert(Abs(a_dx[0]-a_dx[1]) < 1.e-9);

  
  int numCons = CNUM;
  m_BL2Patch.resize(m_BLPtr->nvars(),0);
  for (int i = 0; i < numCons; i++)
    {
      if (CRHO ==  i) m_BL2Patch[m_BLPtr->densityIndex()] = i;
      if (CMOMX == i) {m_BL2Patch[m_BLPtr->momentumIndex()] = i;m_BL2Patch[m_BLPtr->momentumIndex()+1] = i+1;}
      if (CENG ==  i) m_BL2Patch[m_BLPtr->energyIndex()] = i;
    }

  m_isDefined = true;
}
/****************************/
/****************************/
void
EBPin2PinIBC::
initialize(LevelData<EBCellFAB>& a_conState,
           const EBISLayout& a_ebisl) const
{
 CH_assert(m_isDefined);
 CH_assert(m_isFortranCommonSet);

 Real init_dx = (m_BLinit) ? (m_dx) : 0e0;

 if(!m_BLinit) goto uniformFlow;

  // Iterator of all grids in this level
  for (DataIterator dit = a_conState.dataIterator(); dit.ok(); ++dit)
    {
      const EBISBox& ebisBox = a_ebisl[dit()];
      // Storage for current grid
      EBCellFAB& conFAB = a_conState[dit()];

      BaseFab<Real>& regConFAB = conFAB.getSingleValuedFAB();
      // Box of current grid
      Box uBox = regConFAB.box();
      uBox &= m_domain;
      IntVectSet ivs = ebisBox.getMultiCells(uBox); // to be used later on

     
      int nvars = m_BLPtr->nvars();
      const Vector<Vector<Vector<Real> > >& BLspline = m_BLPtr->exportCoeffs();
      const Vector<Vector<Real> >& BLcons = m_BLPtr->exportCons();
      const Vector<Real>& eta = m_BLPtr->exportEta();
      Real x0 = m_BLPtr->getX0(); // note it assumes that theta is assigned at x=0
      Real y0 = m_BLPtr->mp_boundLoc;
      Real ReRoot2 = m_BLPtr->getReLRoot2();
      
      // only 2-Dimensional intialization
      // set BL index to patch index
      if (SpaceDim == 3) conFAB.setVal(CMOMX+2,0.0);
      
      
      RealVect vectDx = m_dx*(RealVect::Unit);
      int ny = uBox.bigEnd(1);
      Vector<Real> y(ny+1); // y is overdimensioned on the bottom
      // lucacancel :: here should have the distance funct intead i set y=ypoint-ywall
      int i = uBox.smallEnd(0);
      int k;
      if (SpaceDim == 3) k = uBox.smallEnd(2);
      for (int j = 0; j <= ny; j++)
	{
	  if (j <= uBox.bigEnd(1) && j>= uBox.smallEnd(1))
	    {
	      IntVect iv(D_DECL(i,j,k));
	      Vector< VolIndex > vofs = ebisBox.getVoFs(iv);
	      if (vofs.size() <= 0)
		y[j] = (j + 0.5e0)*m_dx - y0;
	      else
	      {
		RealVect centroid = ebisBox.centroid(vofs[0])*m_dx;
		RealVect vofloc = EBArith::getVofLocation(vofs[0], vectDx, centroid);
		y[j] = vofloc[1]- y0;
	      }
	    }
	  else
	    y[j] = (j + 0.5e0)*m_dx - y0;
	}

      for (int i = 0; i < nvars; i++)
	{
	  int icomp =  m_BL2Patch[i];

	  // initialize one comp at the time
	  FORT_BLINIT(CHF_CONST_FRA1(regConFAB,icomp),
                      CHF_CONST_INT(i),
		      CHF_CONST_REAL(init_dx),
		      CHF_CONST_REAL(x0),
		      CHF_CONST_REAL(ReRoot2),
		      CHF_CONST_VR(y),
		      CHF_CONST_VR(eta),
		      CHF_CONST_VR(BLcons[i]),
		      CHF_CONST_VR(BLspline[i][0]),
		      CHF_CONST_VR(BLspline[i][1]),
		      CHF_CONST_VR(BLspline[i][2]),
		      CHF_CONST_VR(BLspline[i][3]),
		      CHF_CONST_VR(BLspline[i][4]),
		      CHF_BOX(uBox));

	  for (VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
	    {
	      const VolIndex& vof = vofit();
	      RealVect centroid = ebisBox.centroid(vof)*m_dx;
	      RealVect vofloc = EBArith::getVofLocation(vof, vectDx, centroid);
	      Real yy = vofloc[1]- y0;
	      Real ee = yy*ReRoot2/sqrt(vofloc[0] - x0);
	      Real& yeval = conFAB(vof, icomp);
	      FORT_SPLINE5(CHF_REAL(yeval),
			  CHF_CONST_REAL(ee),
			  CHF_CONST_VR(eta),
			  CHF_CONST_VR(BLcons[i]),
			  CHF_CONST_VR(BLspline[i][0]),
			  CHF_CONST_VR(BLspline[i][1]),
			  CHF_CONST_VR(BLspline[i][2]),
			  CHF_CONST_VR(BLspline[i][3]),
			  CHF_CONST_VR(BLspline[i][4]));
	      if(i == 2) yeval = yeval/(ReRoot2*sqrt(vofloc[0] - x0));
	    }
	}
      
      FORT_ADDKETOU(CHF_CONST_FRA(regConFAB), CHF_BOX(uBox));
      for (VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
	{
	  const VolIndex& vof = vofit();
	  Real Ke = 0;
	  Real Invrhoh = 1.0;
	  for (int i = 0; i < nvars; i++)
	    {
	      int icomp =  m_BL2Patch[i];
	      Real& yeval = conFAB(vof, icomp);
	      if (i == 0)
		Invrhoh = 0.5/yeval;
	      else if(i == 1)
		Ke = (yeval*yeval);
	      else if(i == 2)
		Ke += (yeval*yeval);
	      else if (i == 3)
		yeval += Ke *Invrhoh;
	    }
	}	
	
    } //end loop over boxes
  
 uniformFlow:
  // Iterator of all grids in this level
  for(DataIterator dit = a_conState.dataIterator(); dit.ok(); ++dit)
  {
    const EBISBox& ebisBox = a_ebisl[dit()];
    // Storage for current grid
    EBCellFAB& conFAB = a_conState[dit()];

    BaseFab<Real>& regConFAB = conFAB.getSingleValuedFAB();
    // Box of current grid
    Box uBox = regConFAB.box();
    uBox &= m_domain;
    // Set up initial condition in this grid
    /**/
    FORT_PIN2PININIT(CHF_CONST_FRA(regConFAB), CHF_CONST_REAL(m_dx), CHF_BOX(uBox));
    /**/

    //now for the multivalued cells.  Since it is pointwise,
    //the regular calc is correct for all single-valued cells.
    IntVectSet ivs = ebisBox.getMultiCells(uBox);
    for(VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
      {
        const VolIndex& vof = vofit();
        const IntVect& iv = vof.gridIndex();
        RealVect momentum;
        Real energy, density;
        /**/
        FORT_POINTPIN2PININIT(CHF_REAL(density), CHF_REALVECT(momentum), CHF_REAL(energy), CHF_CONST_INTVECT(iv), CHF_CONST_REAL(m_dx) );
        /**/
        conFAB(vof, CRHO) = density;
        conFAB(vof, CENG) = energy;
        for(int idir = 0; idir < SpaceDim; idir++)
          {
            conFAB(vof, CMOMX+idir) = momentum[idir];
          }
      }//end loop over multivalued cells
  } //end loop over boxes
}
/****************************/
/****************************/
EBPin2PinIBC::
~EBPin2PinIBC()
{
}
/****************************/
/****************************/
