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
#include "CartesianSideDoubleFirstOrderRefine.h"

#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/SideOverlap.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"

#include "pdat/detail.h"

namespace clever {
namespace geom {

using SideData = SAMRAI::pdat::SideData<double>;

CartesianSideDoubleFirstOrderRefine::
CartesianSideDoubleFirstOrderRefine():
  SAMRAI::hier::RefineOperator("RAJA_FIRST_ORDER_REFINE")
{
}

CartesianSideDoubleFirstOrderRefine::
~CartesianSideDoubleFirstOrderRefine()
{
}

int CartesianSideDoubleFirstOrderRefine::getOperatorPriority() const
{
  return 0;
}

SAMRAI::hier::IntVector
CartesianSideDoubleFirstOrderRefine::getStencilWidth(
  const SAMRAI::tbox::Dimension &dim) const
{
  return SAMRAI::hier::IntVector::getZero(dim);
}

void CartesianSideDoubleFirstOrderRefine::refine(
  SAMRAI::hier::Patch& fine,
  const SAMRAI::hier::Patch& coarse,
  const int dst_component,
  const int src_component,
  const SAMRAI::hier::BoxOverlap& fine_overlap,
  const SAMRAI::hier::IntVector& ratio) const
{
  const auto& dim = fine.getDim();

  const SAMRAI::pdat::SideOverlap* t_overlap =
    PTR_CAST<const SAMRAI::pdat::SideOverlap *>(&fine_overlap);

  const intNd r{ratio[0], ratio[1]};

  auto fine_data =
    SAMRAI_SHARED_PTR_CAST<SideData>(fine.getPatchData(dst_component));
  auto coarse_data =
    SAMRAI_SHARED_PTR_CAST<SideData>(coarse.getPatchData(src_component));

  for (int axis = 0; axis < dim.getValue(); axis++) {
    const SAMRAI::hier::BoxContainer& boxes =
      t_overlap->getDestinationBoxContainer(axis);

    for (SAMRAI::hier::BoxContainer::const_iterator b = boxes.begin();
         b != boxes.end(); ++b) {

      SAMRAI::hier::Box fine_box(*b);
      /* fine_box.setUpper(axis, fine_box.upper(axis)-1); */

      const auto fine_side_box = SAMRAI::pdat::SideGeometry::toSideBox(fine_box, axis);

      for (int depth = 0; depth < fine_data->getDepth(); depth++) {
        if (dim == SAMRAI::tbox::Dimension(2)) {
          const Real denominator = (axis == 0) ? r.x : r.y;

          auto fine_array = SAMRAI::pdat::get_view<Dim>(*fine_data, axis, depth);
          auto coarse_array = SAMRAI::pdat::get_const_view<Dim>(*coarse_data, axis, depth);

  FORCEINLINE_MACRO
          SAMRAI::pdat::parallel_for_all(fine_side_box, [=] SAMRAI_HOST_DEVICE (int k, int j) {
            const int ic1 = (k < 0) ? (k+1)/r.y-1 : k/r.y;
            const int ic0 = (j < 0) ? (j+1)/r.x-1 : j/r.x;
            fine_array(j,k)=coarse_array(ic0,ic1)/denominator;
          });

          /*
          if (axis == 0) {
            SideView fine_array = getView(fine_data, Side::X, depth);
            // View2d<Real> fine_array(fine_data->getPointer(0, depth),
            //   fine_data->getGhostBox().numberCells(0)+1,
            //   fine_data->getGhostBox().numberCells(1),
            //   static_cast<int>(fine_data->getGhostBox().lower()[0]),
            //   static_cast<int>(fine_data->getGhostBox().lower()[1]));

            SideView coarse_array = getView(coarse_data, Side::X, depth);
            // View2d<Real> coarse_array(coarse_data->getPointer(0, depth),
            //   coarse_data->getGhostBox().numberCells(0)+1,
            //   coarse_data->getGhostBox().numberCells(1),
            //   static_cast<int>(coarse_data->getGhostBox().lower()[0]),
            //   static_cast<int>(coarse_data->getGhostBox().lower()[1]));

            // TODO: This may need to be getIndexBox<false>(fine_data, Side::X)


            parallel_for_all2(fine_box, [=] SAMRAI_HOST_DEVICE (int k, int j) {
              // RAJA::forallN<PatchPolicy>(isets.getCells(1), isets.getNodes(0), [=] RAJA_DEVICE (int k, int j) {
              const int ic1 = (k < 0) ? (k+1)/r.y-1 : k/r.y;
              const int ic0 = (j < 0) ? (j+1)/r.x-1 : j/r.x;
              fine_array(j,k)=coarse_array(ic0,ic1)/static_cast<Real>(r.x);
            });
          }
          if (axis == 1) {
            SideView fine_array = getView(fine_data, Side::Y, depth);
            // View2d<Real> fine_array(fine_data->getPointer(1, depth),
            //   fine_data->getGhostBox().numberCells(0),
            //   fine_data->getGhostBox().numberCells(1)+1,
            //   static_cast<int>(fine_data->getGhostBox().lower()[0]),
            //   static_cast<int>(fine_data->getGhostBox().lower()[1]));

            SideView coarse_array = getView(coarse_data, Side::Y, depth);
            // View2d<Real> coarse_array(coarse_data->getPointer(1, depth),
            //   coarse_data->getGhostBox().numberCells(0),
            //   coarse_data->getGhostBox().numberCells(1)+1,
            //   static_cast<int>(coarse_data->getGhostBox().lower()[0]),
            //   static_cast<int>(coarse_data->getGhostBox().lower()[1]));

            // TODO: This may need to be getIndexBox<false>(fine_data, Side::Y)

            parallel_for_all2(fine_box, [=] SAMRAI_HOST_DEVICE (int k, int j) {
              // RAJA::forallN<PatchPolicy>(isets.getNodes(1), isets.getCells(0), [=] RAJA_DEVICE (int k, int j) {
              const int ic1 = (k < 0) ? (k+1)/r.y-1 : k/r.y;
              const int ic0 = (j < 0) ? (j+1)/r.x-1 : j/r.x;
              fine_array(j,k)=coarse_array(ic0,ic1)/static_cast<Real>(r.y);
            });
          }
          */
        } else {
          TBOX_ERROR( "CartesianSideFirstOrderRefine::refine dimension != 2 not supported"
                      << std::endl);
        }
      }
    }
  }
}

}
}
