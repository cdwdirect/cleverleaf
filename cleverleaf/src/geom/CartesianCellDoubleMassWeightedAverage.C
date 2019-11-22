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
#include "CartesianCellDoubleMassWeightedAverage.h"

#include "SAMRAI/pdat/ForAll.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/VariableDatabase.h"

#include "SAMRAI/tbox/Utilities.h"

#include "pdat/detail.h"

namespace clever {
namespace geom {

using CellVariable = SAMRAI::pdat::CellVariable<Real>;
using CellData = SAMRAI::pdat::CellData<Real>;

SAMRAI::tbox::StartupShutdownManager::Handler
CartesianCellDoubleMassWeightedAverage::s_initialize_handler(
  CartesianCellDoubleMassWeightedAverage::initializeCallback,
  0,
  0,
  CartesianCellDoubleMassWeightedAverage::finalizeCallback,
  SAMRAI::tbox::StartupShutdownManager::priorityTimers);

std::shared_ptr<SAMRAI::tbox::Timer>
CartesianCellDoubleMassWeightedAverage::t_coarsen;

CartesianCellDoubleMassWeightedAverage::CartesianCellDoubleMassWeightedAverage(
  const SAMRAI::tbox::Dimension& dim):
  SAMRAI::hier::CoarsenOperator("RAJA_MASS_WEIGHTED_COARSEN")
{
}

CartesianCellDoubleMassWeightedAverage::~CartesianCellDoubleMassWeightedAverage()
{
}

bool CartesianCellDoubleMassWeightedAverage::findCoarsenOperator(
  const std::shared_ptr<SAMRAI::hier::Variable>& var,
  const std::string& op_name) const
{
  const std::shared_ptr<CellVariable > cast_var(
    SAMRAI_SHARED_PTR_CAST<CellVariable>(var));

  return (cast_var && (op_name == getOperatorName()));
}

int CartesianCellDoubleMassWeightedAverage::getOperatorPriority() const
{
  /*
   * This priority is lower than the priority of the volume_weighted coarsen,
   * since we need to do that first.
   */
  return 2;
}

SAMRAI::hier::IntVector
CartesianCellDoubleMassWeightedAverage::getStencilWidth(
  const SAMRAI::tbox::Dimension &dim) const {
  return SAMRAI::hier::IntVector::getZero(dim);
}

void CartesianCellDoubleMassWeightedAverage::coarsen(
  SAMRAI::hier::Patch& coarse,
  const SAMRAI::hier::Patch& fine,
  const int dst_component,
  const int src_component,
  const SAMRAI::hier::Box& coarse_box,
  const SAMRAI::hier::IntVector& ratio) const
{
  t_coarsen->start();
  const auto& dim = fine.getDim();

  const auto variable_db = SAMRAI::hier::VariableDatabase::getDatabase();

  const int density_id = variable_db->mapVariableAndContextToIndex(
    variable_db->getVariable("density"),
    variable_db->getContext("CURRENT"));

  const auto fgeom =
    SAMRAI_SHARED_PTR_CAST<SAMRAI::geom::CartesianPatchGeometry>(fine.getPatchGeometry());
  const auto cgeom =
    SAMRAI_SHARED_PTR_CAST<SAMRAI::geom::CartesianPatchGeometry>(coarse.getPatchGeometry());

  auto fine_data =
    SAMRAI_SHARED_PTR_CAST<CellData>(fine.getPatchData(src_component));
  auto coarse_data =
    SAMRAI_SHARED_PTR_CAST<CellData>(coarse.getPatchData(dst_component));

  const double* fdx = fgeom->getDx();
  const double* cdx = cgeom->getDx();

  for (int depth = 0; depth < coarse_data->getDepth(); depth++) {
    auto fine_array = SAMRAI::pdat::get_const_view<Dim>(*fine_data, depth);
    auto coarse_array = SAMRAI::pdat::get_view<Dim>(*coarse_data, depth);

    auto fine_mass = SAMRAI::pdat::get_const_view<Dim, CellData>(fine.getPatchData(density_id), depth);
    auto coarse_mass = SAMRAI::pdat::get_const_view<Dim, CellData>(coarse.getPatchData(density_id), depth);

    const intNd r{ratio[0], ratio[1]};

    if ((dim == SAMRAI::tbox::Dimension(2))) {
      const Real Vf = fdx[0]*fdx[1];
      const Real Vc = cdx[0]*cdx[1];

  FORCEINLINE_MACRO
      SAMRAI::pdat::parallel_for_all(coarse_box, [=] SAMRAI_HOST_DEVICE (int k, int j) {
        Real seM = 0.0;
        for (int r0 = 0; r0 < r.x; r0++) {
          for (int r1 = 0; r1 < r.y; r1++) {
            const int jf = j*r.x+r0;
            const int kf = k*r.y+r1;
            seM = seM + fine_array(jf,kf)*fine_mass(jf,kf)*Vf;
          }
        }
        coarse_array(j,k) = seM/(coarse_mass(j,k)*Vc);
      });
    } else {
      TBOX_ERROR("CartesianCellDoubleMassWeightedAverage error...\n"
                 << "dim != 2 not supported." << std::endl);
    }
  }
  t_coarsen->stop();
}

void CartesianCellDoubleMassWeightedAverage::initializeCallback()
{
  t_coarsen = SAMRAI::tbox::TimerManager::getManager()->getTimer(
    "CartesianCellDoubleMassWeightedAverage::t_coarsen");
}

void CartesianCellDoubleMassWeightedAverage::finalizeCallback()
{
  t_coarsen.reset();
}

}
}
