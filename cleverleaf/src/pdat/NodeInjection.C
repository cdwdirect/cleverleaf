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
#ifndef CLEVERLEAF_PDAT_NODEINJECTION_C_
#define CLEVERLEAF_PDAT_NODEINJECTION_C_

#include "NodeInjection.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/pdat/ForAll.h"

#include "pdat/detail.h"

namespace clever {
namespace pdat {

template<typename TYPE>
NodeInjection<TYPE>::NodeInjection():
  SAMRAI::hier::CoarsenOperator("RAJA_CONSTANT_COARSEN")
{
}

template<typename TYPE>
NodeInjection<TYPE>::~NodeInjection()
{
}

template<typename TYPE>
int NodeInjection<TYPE>::getOperatorPriority() const
{
  return 0;
}

template<typename TYPE>
SAMRAI::hier::IntVector NodeInjection<TYPE>::getStencilWidth(
  const SAMRAI::tbox::Dimension &dim) const
{
  return SAMRAI::hier::IntVector::getZero(dim);
}

template<typename TYPE>
void NodeInjection<TYPE>::coarsen(
  SAMRAI::hier::Patch& coarse,
  const SAMRAI::hier::Patch& fine,
  const int dst_component,
  const int src_component,
  const SAMRAI::hier::Box& coarse_box,
  const SAMRAI::hier::IntVector& ratio) const
{
  using NodeData = SAMRAI::pdat::NodeData<TYPE>;

  const auto& dim = fine.getDim();

  const intNd r{ratio[0], ratio[1]};

  auto fine_data =
    SAMRAI_SHARED_PTR_CAST<NodeData>(fine.getPatchData(src_component));
  auto coarse_data =
    SAMRAI_SHARED_PTR_CAST<NodeData>(coarse.getPatchData(dst_component));

  for(int depth = 0; depth < coarse_data->getDepth(); depth++) {
    if (fine.getDim() == SAMRAI::tbox::Dimension(2)) {
      auto coarse_array = SAMRAI::pdat::get_view<Dim, NodeData>(coarse_data, depth);
      auto fine_array = SAMRAI::pdat::get_const_view<Dim, NodeData>(fine_data, depth);

  FORCEINLINE_MACRO
      SAMRAI::pdat::parallel_for_all(SAMRAI::pdat::NodeGeometry::toNodeBox(coarse_box), [=] SAMRAI_HOST_DEVICE (int k, int j) {
        const int j_fine = j*r.x;
        const int k_fine = k*r.y;
        coarse_array(j,k) = fine_array(j_fine, k_fine);
      });
    } else {
      TBOX_ERROR("NodeInjection error...\n"
                 << "dim != 2 not supported." << std::endl);
    }
  }
}
}
}

#endif
