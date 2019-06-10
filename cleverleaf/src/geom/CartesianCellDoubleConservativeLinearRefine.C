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
#include "CartesianCellDoubleConservativeLinearRefine.h"

#include "SAMRAI/pdat/ForAll.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellOverlap.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/tbox/AllocatorDatabase.h"

#include "pdat/detail.h"

namespace clever {
namespace geom {

using CellData = SAMRAI::pdat::CellData<Real>;

CartesianCellDoubleConservativeLinearRefine::
CartesianCellDoubleConservativeLinearRefine():
  SAMRAI::hier::RefineOperator("RAJA_CONSERVATIVE_LINEAR_REFINE")
{
}

CartesianCellDoubleConservativeLinearRefine::
~CartesianCellDoubleConservativeLinearRefine()
{
}

int CartesianCellDoubleConservativeLinearRefine::getOperatorPriority() const
{
  return 0;
}

SAMRAI::hier::IntVector
CartesianCellDoubleConservativeLinearRefine::getStencilWidth(
  const SAMRAI::tbox::Dimension &dim) const
{
  return SAMRAI::hier::IntVector::getOne(dim);
}

void CartesianCellDoubleConservativeLinearRefine::refine(
  SAMRAI::hier::Patch& fine,
  const SAMRAI::hier::Patch& coarse,
  const int dst_component,
  const int src_component,
  const SAMRAI::hier::BoxOverlap& fine_overlap,
  const SAMRAI::hier::IntVector& ratio) const
{
  const SAMRAI::pdat::CellOverlap* cell_overlap =
    PTR_CAST<const SAMRAI::pdat::CellOverlap *>(&fine_overlap);

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

void CartesianCellDoubleConservativeLinearRefine::refine(
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

  const auto cgbox = coarse.getBox();

  const SAMRAI::hier::IntVector tmp_ghosts(dim, 0);

  auto alloc_db = SAMRAI::tbox::AllocatorDatabase::getDatabase();

  CellData diff_data(SAMRAI::pdat::NodeGeometry::toNodeBox(cgbox), Dim, tmp_ghosts, alloc_db->getDevicePool());
  CellData slope_data(cgbox, Dim, tmp_ghosts, alloc_db->getDevicePool());

  const intNd r{ratio[0], ratio[1]};

  auto diff0 = SAMRAI::pdat::get_view<Dim>(diff_data, Coord::X);
  auto diff1 = SAMRAI::pdat::get_view<Dim>(diff_data, Coord::Y);

  auto slope0 = SAMRAI::pdat::get_view<Dim>(slope_data, Coord::X);
  auto slope1 = SAMRAI::pdat::get_view<Dim>(slope_data, Coord::Y);

  const auto coarse_box = SAMRAI::hier::Box::coarsen(fine_box, ratio);

  const double* dxfine = fine_geometry->getDx();
  const double* dxcoarse = coarse_geometry->getDx();

  const realNd dxf{dxfine[0], dxfine[1]};
  const realNd dxc{dxcoarse[0], dxcoarse[1]};

  const auto coarse_node_box = SAMRAI::pdat::NodeGeometry::toNodeBox(coarse_box);

  auto fine_data =
    SAMRAI_SHARED_PTR_CAST<CellData>(fine.getPatchData(dst_component));
  auto coarse_data =
    SAMRAI_SHARED_PTR_CAST<CellData>(coarse.getPatchData(src_component));

  for(int depth = 0; depth < fine_data->getDepth(); depth++) {
    const auto fine_array = SAMRAI::pdat::get_view<Dim>(*fine_data, depth);
    auto coarse_array = SAMRAI::pdat::get_view<Dim>(*coarse_data, depth);

    if (dim == SAMRAI::tbox::Dimension(2)) {

      SAMRAI::pdat::parallel_for_all(coarse_node_box, [=] SAMRAI_HOST_DEVICE (int k, int j) {
        diff0(j,k) = coarse_array(j,k) - coarse_array(j-1,k);
        diff1(j,k) = coarse_array(j,k) - coarse_array(j,k-1);
      });

      SAMRAI::pdat::parallel_for_all(coarse_box, [=] SAMRAI_HOST_DEVICE (int k, int j) {
        const Real coef2j = 0.5*(diff0(j+1,k)+diff0(j,k));
        const Real boundj = 2.0*MIN(abs(diff0(j+1,k)),abs(diff0(j,k)));

        if (diff0(j,k)*diff0(j+1,k) > 0.0) {
          slope0(j,k) = copysign(MIN(abs(coef2j),boundj),coef2j)/dxc.x;
        } else {
          slope0(j,k) = 0.0;
        }

        const Real coef2k = 0.5*(diff1(j,k+1)+diff1(j,k));
        const Real boundk = 2.0*MIN(abs(diff1(j,k+1)),abs(diff1(j,k)));

        if (diff1(j,k)*diff1(j,k+1) > 0.0) {
          slope1(j,k) = copysign(MIN(abs(coef2k),boundk),coef2k)/dxc.y;
        } else {
          slope1(j,k) = 0.0;
        }
      });

      SAMRAI::pdat::parallel_for_all(fine_box, [=] SAMRAI_HOST_DEVICE (int k, int j) {
        const int ic1 = (k < 0) ? (k+1)/r.y-1 : k/r.y;
        const int ic0 = (j < 0) ? (j+1)/r.x-1 : j/r.x;

        const int ir0 = j - ic0*r.x;
        const int ir1 = k - ic1*r.y;

        const Real deltax1 = (static_cast<Real>(ir1)+0.5)*dxf.y-dxc.y*0.5;
        const Real deltax0 = (static_cast<Real>(ir0)+0.5)*dxf.x-dxc.x*0.5;

        fine_array(j,k) = coarse_array(ic0,ic1) + slope0(ic0, ic1)*deltax0 + slope1(ic0,ic1)*deltax1;
      });

    } else {
      TBOX_ERROR("CartesianCellDoubleConservativeLinearRefine error...\n"
                 << "dim != 2 not supported." << std::endl);
    }
  }
}

}
}
