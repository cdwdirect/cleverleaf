/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2018 Lawrence Livermore National Security, LLC
 * Description:   Weighted averaging operator for outerface complex data on
 *                a Cartesian mesh.
 *
 ************************************************************************/

#ifndef included_geom_CartesianOuterfaceComplexWeightedAverage
#define included_geom_CartesianOuterfaceComplexWeightedAverage

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/CoarsenOperator.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"

#include <string>
#include <memory>

namespace SAMRAI {
namespace geom {

/**
 * Class CartesianOuterfaceComplexWeightedAverage implements conservative
 * face-weighted averaging for outerface complex patch data defined over
 * a Cartesian mesh.  It is derived from the hier::CoarsenOperator base class.
 * The numerical operations for theaveraging use FORTRAN numerical routines.
 *
 * @see hier::CoarsenOperator
 */

class CartesianOuterfaceComplexWeightedAverage:
   public hier::CoarsenOperator
{
public:
   /**
    * Uninteresting default constructor.
    */
   CartesianOuterfaceComplexWeightedAverage();

   /**
    * Uninteresting virtual destructor.
    */
   virtual ~CartesianOuterfaceComplexWeightedAverage();

   /**
    * The priority of outerface complex weighted averaging is 0.
    * It will be performed before any user-defined coarsen operations.
    */
   int
   getOperatorPriority() const;

   /**
    * The stencil width of the weighted averaging operator is the vector of
    * zeros.  That is, its stencil does not extend outside the fine box.
    */
   hier::IntVector
   getStencilWidth(
      const tbox::Dimension& dim) const;

   /**
    * Coarsen the source component on the fine patch to the destination
    * component on the coarse patch using the outerface complex weighted
    * averaging operator.  Coarsening is performed on the intersection of
    * the destination patch and the coarse box.  It is assumed that the
    * fine patch contains sufficient data for the stencil width of the
    * coarsening operator.
    *
    * @pre (fine.getDim() == coarse.getDim()) &&
    *      (fine.getDim() == coarse_box.getDim()) &&
    *      (fine.getDim() == ratio.getDim())
    * @pre fine.getPatchData(src_component) is actually a std::shared_ptr<pdat::OuterfaceData<dcomplex> >
    * @pre coarse.getPatchData(dst_component) is actually a std::shared_ptr<pdat::OuterfaceData<dcomplex> >
    * @pre fine.getPatchData(src_component)->getDepth() == coarse.getPatchData(dst_component)->getDepth()
    * @pre (fine.getDim().getValue() == 1) ||
    *      (fine.getDim().getValue() == 2) || (fine.getDim().getValue() == 3)
    */
   void
   coarsen(
      hier::Patch& coarse,
      const hier::Patch& fine,
      const int dst_component,
      const int src_component,
      const hier::Box& coarse_box,
      const hier::IntVector& ratio) const;
};

}
}

#endif
