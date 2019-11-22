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
#include "CartesianCellIntConstantCoarsen.h"

#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/VariableDatabase.h"

#include "pdat/detail.h"

namespace clever {
namespace geom {

using IndicatorData = SAMRAI::pdat::CellData<int>;

SAMRAI::tbox::StartupShutdownManager::Handler
CartesianCellIntConstantCoarsen::s_initialize_handler(
  CartesianCellIntConstantCoarsen::initializeCallback,
  0,
  0,
  CartesianCellIntConstantCoarsen::finalizeCallback,
  SAMRAI::tbox::StartupShutdownManager::priorityTimers);

std::shared_ptr<SAMRAI::tbox::Timer> CartesianCellIntConstantCoarsen::t_coarsen;

CartesianCellIntConstantCoarsen::CartesianCellIntConstantCoarsen(
  const SAMRAI::tbox::Dimension& dim):
  SAMRAI::hier::CoarsenOperator("RAJA_CONSTANT_INDICATOR_COARSEN")
{
}

CartesianCellIntConstantCoarsen::~CartesianCellIntConstantCoarsen()
{
}

bool CartesianCellIntConstantCoarsen::findCoarsenOperator(
  const std::shared_ptr<SAMRAI::hier::Variable>& var,
  const std::string& op_name) const
{
  const std::shared_ptr<TagVariable > cast_var(
    SAMRAI_SHARED_PTR_CAST<TagVariable>(var));

  return (cast_var && (op_name == getOperatorName()));
}

int CartesianCellIntConstantCoarsen::getOperatorPriority() const
{
  return 0;
}

SAMRAI::hier::IntVector CartesianCellIntConstantCoarsen::getStencilWidth(
  const SAMRAI::tbox::Dimension &dim) const
{
  return SAMRAI::hier::IntVector::getZero(dim);
}

void CartesianCellIntConstantCoarsen::coarsen(
  SAMRAI::hier::Patch& coarse,
  const SAMRAI::hier::Patch& fine,
  const int dst_component,
  const int src_component,
  const SAMRAI::hier::Box& coarse_box,
  const SAMRAI::hier::IntVector& ratio) const
{
  t_coarsen->start();
  const auto& dim = fine.getDim();

  const intNd r{ratio[0], ratio[1]};

  auto fine_data =
    SAMRAI_SHARED_PTR_CAST<NodeData>(fine.getPatchData(src_component));
  auto coarse_data =
    SAMRAI_SHARED_PTR_CAST<NodeData>(coarse.getPatchData(dst_component));

  for (int d = 0; d < coarse_data->getDepth(); d++) {
    if ((dim == SAMRAI::tbox::Dimension(2))) {
      auto fine_array = SAMRAI::pdat::get_const_view<Dim>(*fine_data, d);
      auto coarse_array = SAMRAI::pdat::get_view<Dim>(*coarse_data, d);

  FORCEINLINE_MACRO
      SAMRAI::pdat::parallel_for_all(coarse_box, [=] SAMRAI_HOST_DEVICE (int k, int j) {
        const int j_fine = j*r.x;
        const int k_fine = k*r.y;
        coarse_array(j,k) = fine_array(j_fine, k_fine);
      });
    } else {
      TBOX_ERROR("CartesianCellIntConstantCoarsen error...\n"
                 << "dim != 2 not supported." << std::endl);
    }
  }
  t_coarsen->stop();
}

void CartesianCellIntConstantCoarsen::initializeCallback()
{
  t_coarsen = SAMRAI::tbox::TimerManager::getManager()->getTimer(
    "CartesianCellIntConstantCoarsen::t_coarsen");
}

void CartesianCellIntConstantCoarsen::finalizeCallback()
{
  t_coarsen.reset();
}

}
}
