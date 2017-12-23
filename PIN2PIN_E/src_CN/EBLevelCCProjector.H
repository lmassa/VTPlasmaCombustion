#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBLEVELCCPROJECTOR_H_
#define _EBLEVELCCPROJECTOR_H_
#include "DisjointBoxLayout.H"
#include "EBISLayout.H"
#include "Box.H"
#include "REAL.H"
#include "LevelData.H"
#include "EBFluxFAB.H"
#include "EBCellFAB.H"
#include "EBSimpleSolver.H"
#include "EBAMRPoissonOp.H"
#include "MultiGrid.H"
#include "EBQuadCFInterp.H"
#include "NamespaceHeader.H"

///
/**
   Class to project a face-centered. velocity field on a level.  \\
   u -= G (DG^-1)(D u)
   This is done as interface to the level mac projector.  Velocity is averaged to faces
   then a mac projection is executed.  The gradient is then averaged back to cells
   and subtracted off the velocity.  So it really looks more like \\
   u -=     AveFtoC(G(DG^-1)(D(AveCToF(u))))
   Ghost cells over the coarse-fine interface are filled
   by quadradic coarse-fine interpolation.
   This class does not assume that the boundary conditions at the embedded boundary are no-flow.
 */

//reusable functions made external
extern void
ccpAverageVelocityToFaces(LevelData<EBFluxFAB>&         a_macVeloc,
                          const LevelData<EBCellFAB>&   a_cellVeloc,
                          const DisjointBoxLayout&      a_grids,
                          const EBISLayout&             a_ebisl,
                          const ProblemDomain&          a_domain,
                          const RealVect&               a_dx,
                          const LayoutData<IntVectSet>& a_cfivs);

extern void
ccpAverageVelocityToFaces(EBFaceFAB&             a_faceVel,
                          const EBCellFAB&       a_cellVel,
                          const EBGraph&         a_ebGraph,
                          const Box&             a_grid,
                          const int&             a_idir,
                          const ProblemDomain&   a_domain,
                          const RealVect&        a_dx);

extern void
ccpAverageVelocityToFaces(LevelData<EBFluxFAB>&         a_macVeloc,
                          const LevelData<EBCellFAB>&   a_cellVeloc,
                          const DisjointBoxLayout&      a_grids,
                          const EBISLayout&             a_ebisl,
                          const ProblemDomain&          a_domain,
                          const RealVect&               a_dx,
                          const LayoutData<IntVectSet>& a_cfivs,
                          const int&                    a_comp);

extern void
ccpAverageVelocityToFaces(EBFaceFAB&             a_faceVel,
                          const EBCellFAB&       a_cellVel,
                          const EBGraph&         a_ebGraph,
                          const Box&             a_grid,
                          const int&             a_idir,
                          const ProblemDomain&   a_domain,
                          const RealVect&        a_dx,
                          const int&             a_comp);

extern void
ccpAverageStressToFaces(LevelData<EBFluxFAB>&         a_macVeloc,
                          const LevelData<EBCellFAB>&   a_cellVeloc,
                          const DisjointBoxLayout&      a_grids,
                          const EBISLayout&             a_ebisl,
                          const ProblemDomain&          a_domain,
                          const RealVect&               a_dx,
                          const LayoutData<IntVectSet>& a_cfivs);

extern void
ccpAverageStressToFaces(EBFaceFAB&             a_faceVel,
                        const EBCellFAB&       a_cellVel,
                        const EBGraph&         a_ebGraph,
                        const Box&             a_grid,
                        const int&             a_idir,
                        const ProblemDomain&   a_domain,
                        const RealVect&        a_dx);

extern void
ccpAverageFaceToCells(LevelData<EBCellFAB>&         a_cellData,
                      const LevelData<EBFluxFAB>&   a_macData,
                      const DisjointBoxLayout&      a_grids,
                      const EBISLayout&             a_ebisl,
                      const ProblemDomain&          a_domain,
                      const RealVect&               a_dx);

extern void
ccpAverageFaceToCells(EBCellFAB&             a_cellData,
                      const EBFluxFAB&       a_fluxData,
                      const EBGraph&         a_ebGraph,
                      const Box&             a_grid,
                      const ProblemDomain&   a_domain,
                      const RealVect&        a_dx);

extern void
ccpAverageFaceToCellsScalar(LevelData<EBCellFAB>&         a_cellData,
                            const LevelData<EBFluxFAB>&   a_macData,
                            const DisjointBoxLayout&      a_grids,
                            const EBISLayout&             a_ebisl,
                            const ProblemDomain&          a_domain,
                            const RealVect&               a_dx,
                            const int &                   a_dir);

extern void
ccpAverageFaceToCellsScalar(EBCellFAB&             a_cellData,
                            const EBFluxFAB&       a_fluxData,
                            const EBGraph&         a_ebGraph,
                            const Box&             a_grid,
                            const ProblemDomain&   a_domain,
                            const RealVect&        a_dx,
                            const int &            a_dir);

extern void
ccpCellGradientFromFaceData(LevelData<EBCellFAB>&         a_cellData,
                            const LevelData<EBFluxFAB>&   a_macData,
                            const DisjointBoxLayout&      a_grids,
                            const EBISLayout&             a_ebisl,
                            const ProblemDomain&          a_domain,
                            const RealVect&               a_dx,
                            const int &                   a_dir);

extern void
ccpCellGradientFromFaceData(EBCellFAB&             a_cellData,
                            const EBFluxFAB&       a_fluxData,
                            const EBGraph&         a_ebGraph,
                            const Box&             a_grid,
                            const ProblemDomain&   a_domain,
                            const RealVect&        a_dx,
                            const int &            a_dir);

extern void
ccpExtrapolateToDomainBoundaries(LevelData<EBFluxFAB>&    a_macData,
                                 const DisjointBoxLayout& a_grids,
                                 const EBISLayout&        a_ebisl,
                                 const ProblemDomain&     a_domain,
                                 const RealVect&          a_dx);

extern void
ccpExtrapolateToDomainBoundaries(EBFaceFAB&             a_faceData,
                                 const EBGraph&         a_ebGraph,
                                 const Box&             a_grid,
                                 const int&             a_idir,
                                 const ProblemDomain&   a_domain,
                                 const RealVect&        a_dx);

extern void
ccpLinearInterp(Real &                    a_dataOnLine,
                const Vector<RealVect>&   a_faceLoc,
                const RealVect&           a_intersectLoc,
                const Vector<Real>&       a_interpolationData,
                const int&                a_planeDir);

extern void
ccpBilinearInterp(Real &                             a_dataOnLine,
                  const Vector<RealVect>&            a_faceLoc,
                  const RealVect&                    a_intersectLoc,
                  const Vector<Real>&                a_interpolationData,
                  const Tuple<int, SpaceDim-1>&      a_planeDir);


extern void
ccpExtrapFaceToCovered(bool&                   a_dropOrder,
                           Real&                   a_extrapVal,
                           const EBFaceFAB&        a_primFace,
                           const VolIndex&         a_vof,
                           const int&              a_faceDir,
                           const Side::LoHiSide&   a_sd,
                           const RealVect&         a_normal,
                           const RealVect&         a_dx,
                           const int&              a_icomp);

///
/**
   Multi-D extrapolation to covered face
*/
extern void
ccpJohansenExtrapFaceToCovered(bool&                   a_dropOrder,
                               Real&                   a_extrapVal,
                               const EBFaceFAB&        a_primFace,
                               const EBISBox&          a_ebisBox,
                               const VolIndex&         a_vof,
                               const int&              a_faceDir,
                               const Side::LoHiSide&   a_sd,
                               const RealVect&         a_normal,
                               const RealVect&         a_dx,
                               const int&              a_icomp) ;

///
/**
   do one dimensional extrapolation to covered face
   preferring face direction
*/
extern Real
ccpOneDCoveredExtrapValue(const VolIndex&         a_vof,
                          const int&              a_idir,
                          const Side::LoHiSide    a_side,
                          const EBFaceFAB &       a_faceGrad,
                          const EBISBox &         a_ebisBox,
                          const Box &             a_grid,
                          const ProblemDomain &   a_domain,
                          const RealVect &        a_dx,
                          const int&              a_icomp);

///
/**
   do multi-D  extrapolation to covered face
   if possible.  if not, drop to 1d
   preferring face direction
*/
extern Real
ccpGetCoveredExtrapValue(const VolIndex&         a_vof,
                         const int&              a_idir,
                         const Side::LoHiSide    a_side,
                         const EBFaceFAB &       a_faceGrad,
                         const EBISBox &         a_ebisBox,
                         const Box &             a_grid,
                         const ProblemDomain &   a_domain,
                         const RealVect &        a_dx,
                         const int&              a_icomp);
/*****/
extern void
macGradient(LevelData<EBFluxFAB>&        a_gradPhi,
            const LevelData<EBCellFAB>&  a_phi,
            const DisjointBoxLayout&     a_grids,
            const EBISLayout&            a_ebisl,
            const ProblemDomain&         a_domain,
            const RealVect&              a_dx);
/****/
extern void
macGradient(EBFluxFAB&           a_gradPhi,
            const EBCellFAB&     a_phi,
            const EBISBox&       a_ebisBox,
            const Box&           a_box,
            const ProblemDomain& a_domain,
            const RealVect&      a_dx);
extern Real
extrapFaceGradToBoundary(const FaceIndex&      a_bndryFace,
			 const Side::LoHiSide& a_side,
			 const int&            a_idir,
			 const EBGraph&        a_ebGraph,
			 const EBFaceFAB&      a_faceData,
			 const int&            a_comp);
#include "NamespaceFooter.H"
#endif
