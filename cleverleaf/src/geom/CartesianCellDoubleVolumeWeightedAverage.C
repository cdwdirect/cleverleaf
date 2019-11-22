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
#include "CartesianCellDoubleVolumeWeightedAverage.h"

#include "SAMRAI/pdat/ForAll.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/tbox/Utilities.h"

#include "pdat/detail.h"

namespace clever {
namespace geom {

using CellVariable = SAMRAI::pdat::CellVariable<Real>;
using CellData = SAMRAI::pdat::CellData<Real>;

SAMRAI::tbox::StartupShutdownManager::Handler
CartesianCellDoubleVolumeWeightedAverage::s_initialize_handler(
  CartesianCellDoubleVolumeWeightedAverage::initializeCallback,
  0,
  0,
  CartesianCellDoubleVolumeWeightedAverage::finalizeCallback,
  SAMRAI::tbox::StartupShutdownManager::priorityTimers);

std::shared_ptr<SAMRAI::tbox::Timer>
CartesianCellDoubleVolumeWeightedAverage::t_coarsen;

CartesianCellDoubleVolumeWeightedAverage::CartesianCellDoubleVolumeWeightedAverage(
  const SAMRAI::tbox::Dimension& dim):
  SAMRAI::hier::CoarsenOperator("RAJA_VOLUME_WEIGHTED_COARSEN")
{
}

CartesianCellDoubleVolumeWeightedAverage::~CartesianCellDoubleVolumeWeightedAverage()
{
}

bool CartesianCellDoubleVolumeWeightedAverage::findCoarsenOperator(
  const std::shared_ptr<SAMRAI::hier::Variable>& var,
  const std::string& op_name) const
{
  const std::shared_ptr<CellVariable > cast_var(
    SAMRAI_SHARED_PTR_CAST<CellVariable>(var));

  return (cast_var && (op_name == getOperatorName()));
}

int CartesianCellDoubleVolumeWeightedAverage::getOperatorPriority() const
{
  return 1;
}

SAMRAI::hier::IntVector CartesianCellDoubleVolumeWeightedAverage::getStencilWidth(
  const SAMRAI::tbox::Dimension &dim) const
{
  return SAMRAI::hier::IntVector::getZero(dim);
}

void CartesianCellDoubleVolumeWeightedAverage::coarsen(
  SAMRAI::hier::Patch& coarse,
  const SAMRAI::hier::Patch& fine,
  const int dst_component,
  const int src_component,
  const SAMRAI::hier::Box& coarse_box,
  const SAMRAI::hier::IntVector& ratio) const
{
  t_coarsen->start();
  const auto& dim = fine.getDim();

  const auto fgeom =
    SAMRAI_SHARED_PTR_CAST<SAMRAI::geom::CartesianPatchGeometry>(fine.getPatchGeometry());
  const auto cgeom =
    SAMRAI_SHARED_PTR_CAST<SAMRAI::geom::CartesianPatchGeometry>(coarse.getPatchGeometry());

  const double* fdx = fgeom->getDx();
  const double* cdx = cgeom->getDx();

  auto fine_data =
    SAMRAI_SHARED_PTR_CAST<CellData>(fine.getPatchData(src_component));
  auto coarse_data =
    SAMRAI_SHARED_PTR_CAST<CellData>(coarse.getPatchData(dst_component));

  for (int depth = 0; depth < coarse_data->getDepth(); depth++) {
    auto fine_array = SAMRAI::pdat::get_const_view<Dim>(*fine_data, depth);
    auto coarse_array = SAMRAI::pdat::get_view<Dim>(*coarse_data, depth);

    const intNd r{ratio[0], ratio[1]};

    if ((dim == SAMRAI::tbox::Dimension(2))) {
      const Real Vf = fdx[0]*fdx[1];
      const Real Vc = cdx[0]*cdx[1];

  FORCEINLINE_MACRO
      SAMRAI::pdat::parallel_for_all(coarse_box, [=] SAMRAI_HOST_DEVICE (int k, int j) {
        Real spv = 0.0;
        for (int r0 = 0; r0 < r.x; r0++) {
          for (int r1 = 0; r1 < r.y; r1++) {
            const int jf = j*r.x+r0;
            const int kf = k*r.y+r1;
            spv = spv + fine_array(jf,kf)*Vf;
          }
        }
        coarse_array(j,k) = spv/Vc;
      });

    } else {
      TBOX_ERROR("CartesianCellDoubleVolumeWeightedAverage error...\n"
                 << "dim != 2 not supported." << std::endl);
    }
  }
  t_coarsen->stop();
}

void CartesianCellDoubleVolumeWeightedAverage::initializeCallback()
{
  t_coarsen = SAMRAI::tbox::TimerManager::getManager()->getTimer(
    "CartesianCellDoubleVolumeWeightedAverage::t_coarsen");
}

void CartesianCellDoubleVolumeWeightedAverage::finalizeCallback()
{
  t_coarsen.reset();
}

}
}
