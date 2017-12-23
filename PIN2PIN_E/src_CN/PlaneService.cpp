#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//  ANAG, LBNL, DTG

#include "PlaneService.H"
#include "BoxIterator.H"
#include "VoFIterator.H"
#include "PolyGeom.H"
#include "NamespaceHeader.H"

/******************/
//[MD] This PlaneService generates a plan EB. A box with 2 cell height in the y-direction is 
//provided as input (a_coveredRegion). The program generates the plan between these 2 cell rows.
/******************/

PlaneService::PlaneService(const Box& a_coveredRegion)
  :GeometryService()
{
  m_coveredRegion = a_coveredRegion;
  PolyGeom::setVectDx(RealVect::Unit);
}
/******************/
/******************/
PlaneService::~PlaneService()
{
}
/******************/
/******************/
bool
PlaneService::isRegular(const Box& a_region,
                       const ProblemDomain& a_domain,
                       const RealVect& a_origin,
                       const Real& a_dx) const
{
  //use grown box in case the region and
  //the box are colinear
  //[MD] grow box only in directions different from y since input Box already 
  //includes the irregular cells adjacent to the plane in this direction 
	Box grownBox = a_region;
    for (int idir = 0; ((idir < SpaceDim) && (idir !=1)); idir++)
		{
			grownBox.grow(idir, 1);
		}
  //
  grownBox &= a_domain;
  Box interBox = m_coveredRegion & grownBox;
  return (interBox.isEmpty());
}
/******************/
/******************/
bool
PlaneService::isCovered(const Box& a_region,
                       const ProblemDomain& a_domain,
                       const RealVect& a_origin,
                       const Real& a_dx) const
{
  return (m_coveredRegion.contains(a_region));
}
/******************/
/******************/

GeometryService::InOut PlaneService::InsideOutside(const Box&           a_region,
                                                  const ProblemDomain& a_domain,
                                                  const RealVect&      a_origin,
                                                  const Real&          a_dx) const
{
  if (isRegular(a_region, a_domain, a_origin, a_dx)) return GeometryService::Regular;
  if (isCovered(a_region, a_domain, a_origin, a_dx)) return GeometryService::Covered;
  return GeometryService::Irregular;
}

/********************/

void
PlaneService::fillGraph(BaseFab<int>&     a_regIrregCovered,
                       Vector<IrregNode>& a_nodes,
                       const Box&         a_validRegion,
                       const Box&         a_ghostRegion,
                       const ProblemDomain&         a_domain,
                       const RealVect&    a_origin,
                       const Real&        a_dx) const
{

  Box grownCovBox = m_coveredRegion;
  //[MD] grow box only in directions different from y	
	for (int idir = 0; ((idir < SpaceDim) && (idir !=1)); idir++)
		{
			grownCovBox.grow(idir,1);
		} 
  grownCovBox &= a_ghostRegion;
  IntVectSet ivsIrreg(grownCovBox);
  //[MD] commented out:
  //ivsIrreg -= m_coveredRegion;
  ivsIrreg &= a_domain;

  //set regirregcoverred flags
  CH_assert(a_regIrregCovered.box().contains(a_ghostRegion));

  //set every cell to regular
  a_regIrregCovered.setVal(1);

  for (BoxIterator bit(a_ghostRegion); bit.ok(); ++bit)
    {
      if (m_coveredRegion.contains(bit()))
        {
          //set covered cells to -1
          //[MD] changed to zero (irregular cells)
          a_regIrregCovered(bit(), 0) = 0;
        }
      else if (ivsIrreg.contains(bit()))
        {
          //set irreg cells to 0
          a_regIrregCovered(bit(), 0) = 0;
        }
    }

  //now loop through irreg cells and make Nodes for them
  a_nodes.resize(0);
  for (IVSIterator ivsit(ivsIrreg); ivsit.ok(); ++ivsit)
    {
      const IntVect& iv = ivsit();
      if (a_validRegion.contains(iv))
        {
          IrregNode node;
          //first the obvious
          node.m_cell = iv;
          node.m_volFrac = 1.0;
          node.m_cellIndex = 0;
          node.m_volCentroid    = RealVect::Zero;
          //any time the next cell over is in the covered
          //region, there is no face.  If there is a cell
          //but it is outside the domain, the arc=-1 to signify
          //a boundary face.  Otherwise, the arc is -2 if it
          //is to a regular cell and 0 if it is to an irregular cell
		
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              for (SideIterator sit; sit.ok(); ++sit)
                {
                  int arcIndex = node.index(idir, sit());
                  Vector<int>&      arcs     = node.m_arc[arcIndex];
                  Vector<Real>&     areaFracs= node.m_areaFrac[arcIndex];
                  Vector<RealVect>& faceCents= node.m_faceCentroid[arcIndex];
                  IntVect otherIV = iv + sign(sit())*BASISV(idir);
				  
				  int otherCellIndex;
			      Real areaFrac;
    			  
                  //if (m_coveredRegion.contains(otherIV)) 
                  //[MD] The EB face is recognized by its adjacency to 2 irregular cells 
				  //in the y-direction; for such faces, the area fraction is set to zero:
				  if ((ivsIrreg.contains(otherIV)) && (idir==1))
                    {
					  //do nothing, covered face. leave the vector empty
                    }
                  else
                    {
                      //int otherCellIndex;
					  //Real areaFrac;
                      if (ivsIrreg.contains(otherIV))
                        {
                          //arc irregular cell inside the domain
                          otherCellIndex = 0;
						  areaFrac = 1.0;
                        }
                      else if (!a_domain.contains(otherIV))
                        {
                          //boundary face
                          otherCellIndex = -1;
						  areaFrac = 1.0;
                        }
                      else
                        {
                          //arc to regular cell
                          otherCellIndex = -2;
						  areaFrac = 1.0;
                        }
                      arcs.push_back(otherCellIndex);
                      //Real     areaFrac = 1.0;
                      RealVect faceCent = RealVect::Zero;

                      areaFracs.push_back(areaFrac);
                      faceCents.push_back(faceCent);
                    } //end otherIV not covered
                }
            }

          //the boundary centroid and normal
          //depend on which side is covered

          for (int idir = 0; idir < SpaceDim; idir++)
            {
              for (SideIterator sit; sit.ok(); ++sit)
                {
                  int arcIndex = node.index(idir, sit());
                  Vector<int>& arcs = node.m_arc[arcIndex];
				  //[MD] checking 0-area fraction to recognize EB faces:
                  if (arcs.size() == 0)				  
                    {
                      int isign = sign(sit());
                      node.m_bndryCentroid[idir] = Real(isign)*0.5;
                    }
                }
            }
          a_nodes.push_back(node);
        }
    }
}
/******************/
/******************/
#include "NamespaceFooter.H"
