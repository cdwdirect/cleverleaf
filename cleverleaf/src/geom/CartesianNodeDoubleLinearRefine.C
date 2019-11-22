//////////////////////////////////////////////////////////////////////////////
// Crown Copyright 2014 AWE, Copyright 2014 David Beckingsale.
//
// This file is part of CleverLeaf.
//
// CleverLeaf is free software: you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version.
//
// CleverLeaf is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE. See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with
// CleverLeaf. If not, see http://www.gnu.org/licenses/.
//////////////////////////////////////////////////////////////////////////////
#include "CartesianNodeDoubleLinearRefine.h"

#include "SAMRAI/pdat/ForAll.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/pdat/NodeOverlap.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"

#include "pdat/detail.h"

namespace clever {
namespace geom {

using NodeData = SAMRAI::pdat::NodeData<Real>;

CartesianNodeDoubleLinearRefine::
CartesianNodeDoubleLinearRefine():
  SAMRAI::hier::RefineOperator("RAJA_LINEAR_REFINE")
{
}

CartesianNodeDoubleLinearRefine::
~CartesianNodeDoubleLinearRefine()
{
}

int CartesianNodeDoubleLinearRefine::getOperatorPriority() const
{
  return 0;
}

SAMRAI::hier::IntVector
CartesianNodeDoubleLinearRefine::getStencilWidth(
  const SAMRAI::tbox::Dimension &dim) const
{
  return SAMRAI::hier::IntVector::getZero(dim);
}

void CartesianNodeDoubleLinearRefine::refine(
  SAMRAI::hier::Patch& fine,
  const SAMRAI::hier::Patch& coarse,
  const int dst_component,
  const int src_component,
  const SAMRAI::hier::BoxOverlap& fine_overlap,
  const SAMRAI::hier::IntVector& ratio) const
{
  const SAMRAI::pdat::NodeOverlap* cell_overlap =
    PTR_CAST<const SAMRAI::pdat::NodeOverlap *>(&fine_overlap);

  const SAMRAI::hier::BoxContainer& boxes =
    cell_overlap->getDestinationBoxContainer();

  for (SAMRAI::hier::BoxContainer::const_iterator box = boxes.begin();
       box != boxes.end();
       ++box) {
    refine(fine,
           coarse,
           dst_component,
           src_component,
           *box,
           ratio);
  }
}

void CartesianNodeDoubleLinearRefine::refine(
  SAMRAI::hier::Patch& fine,
  const SAMRAI::hier::Patch& coarse,
  const int dst_component,
  const int src_component,
  const SAMRAI::hier::Box& fine_box,
  const SAMRAI::hier::IntVector& ratio) const
{
  const auto& dim = fine.getDim();

  const auto coarse_geometry =
    SAMRAI_SHARED_PTR_CAST<SAMRAI::geom::CartesianPatchGeometry>(coarse.getPatchGeometry());

  const auto fine_geometry =
    SAMRAI_SHARED_PTR_CAST<SAMRAI::geom::CartesianPatchGeometry>(fine.getPatchGeometry());

  const auto coarse_box = SAMRAI::hier::Box::coarsen(fine_box, ratio);

  const intNd r{ratio[0], ratio[1]};

  auto fine_data =
    SAMRAI_SHARED_PTR_CAST<NodeData>(fine.getPatchData(dst_component));
  auto coarse_data =
    SAMRAI_SHARED_PTR_CAST<NodeData>(coarse.getPatchData(src_component));

  for(int depth = 0; depth < fine_data->getDepth(); depth++) {
    auto fine_array = SAMRAI::pdat::get_view<Dim>(*fine_data, depth);
    auto coarse_array = SAMRAI::pdat::get_const_view<Dim>(*coarse_data, depth);

    if (dim == SAMRAI::tbox::Dimension(2)) {
      const Real realrat0 = 1.0/static_cast<Real>(ratio(0));
      const Real realrat1 = 1.0/static_cast<Real>(ratio(1));

      // TODO I'm not sure if this should be cells or nodes...
  FORCEINLINE_MACRO
      SAMRAI::pdat::parallel_for_all(fine_box, [=] SAMRAI_HOST_DEVICE (int k, int j) {
        const int ic0 = floor(static_cast<Real>(j)/r.x);
        const int ic1 = floor(static_cast<Real>(k)/r.y);

        const int ir0 = j - (ic0*r.x);
        const int ir1 = k - (ic1*r.y);

        const Real x = static_cast<Real>(ir0)*realrat0;
        const Real y = static_cast<Real>(ir1)*realrat1;

        fine_array(j,k) = (coarse_array(ic0,ic1)*(1.0-x)
                           + coarse_array(ic0+1,ic1)*x)*(1.0-y)
          + (coarse_array(ic0,ic1+1)*(1.0-x)
             + coarse_array(ic0+1,ic1+1)*x)*y;
      });
    } else {
      TBOX_ERROR("CartesianNodeDoubleLinearRefine error...\n"
                 << "dim != 2 not supported." << std::endl);
    }
  }
}

}
}
