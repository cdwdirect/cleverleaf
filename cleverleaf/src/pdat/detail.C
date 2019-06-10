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
#include <cstdlib>
#include "detail.h"
#if defined(ENABLE_UVM)
#include "cuda_runtime.h"
#endif

namespace clever {

void* allocate(size_t size)
{
  void* ret = NULL;
#if defined(ENABLE_UVM)
  cudaMallocManaged(&ret, size);
#else
  ret = std::malloc(size);
#endif

  return ret;
}

void free(void* ptr)
{
#if defined(ENABLE_UVM)
  cudaFree(ptr);
#else
  std::free(ptr);
#endif
}

template <>
SAMRAI::hier::Box getIndexBox<true>(const CellData& data)
{
  return data.getArrayData().getBox();
}

template <>
SAMRAI::hier::Box getIndexBox<true>(const NodeData& data)
{
  return data.getArrayData().getBox();
}

template <>
SAMRAI::hier::Box getIndexBox<true>(const SideData& data, int side_normal)
{
  return data.getArrayData(side_normal).getBox();
}

template <>
SAMRAI::hier::Box getIndexBox<false>(const CellData& data)
{
  return SAMRAI::pdat::CellGeometry::toCellBox(data.getBox());
}

template <>
SAMRAI::hier::Box getIndexBox<false>(const NodeData& data)
{
  return SAMRAI::pdat::NodeGeometry::toNodeBox(data.getBox());
}

template <>
SAMRAI::hier::Box getIndexBox<false>(const SideData& data, int side_normal)
{
  // The following includes the ghosts, since ArrayData is allocated with ghosts
  return SAMRAI::pdat::SideGeometry::toSideBox(data.getBox(), side_normal);
}

} // end namespace clever
