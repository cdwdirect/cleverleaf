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
#ifndef CLEVERLEAF_PDAT_DETAIL_H_
#define CLEVERLEAF_PDAT_DETAIL_H_

#ifndef NDEBUG
#define PTR_CAST dynamic_cast
#else
#define PTR_CAST static_cast
#endif // DEBUG

#include <cmath>
#include <type_traits>

#define MAX(a, b) (((b) > (a)) ? (b) : (a))
#define MIN(a, b) (((b) < (a)) ? (b) : (a))

#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/NodeVariable.h"
#include "SAMRAI/pdat/SideVariable.h"

#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/SideData.h"

namespace clever {

typedef double Real;

struct intNd { int x; int y; };
struct realNd { Real x; Real y; };

void* allocate(size_t size);
void free(void* ptr);

struct Coord { enum { X = 0, Y = 1 }; };
struct Side { enum { X = 0, Y = 1 }; };
constexpr int Dim = 2;

typedef SAMRAI::pdat::CellVariable<Real> CellVariable;
typedef SAMRAI::pdat::NodeVariable<Real> NodeVariable;
typedef SAMRAI::pdat::SideVariable<Real> SideVariable;
typedef SAMRAI::pdat::CellVariable<int> IndicatorVariable;
typedef SAMRAI::pdat::CellVariable<int> TagVariable;

typedef SAMRAI::pdat::CellData<Real> CellData;
typedef SAMRAI::pdat::NodeData<Real> NodeData;
typedef SAMRAI::pdat::SideData<Real> SideData;
typedef SAMRAI::pdat::CellData<int> IndicatorData;

typedef CellData::View<Dim> CellView;
typedef NodeData::View<Dim> NodeView;
typedef SideData::View<Dim> SideView;
typedef IndicatorData::View<Dim> TagView;
typedef IndicatorData::View<Dim> IndicatorView;

template<typename Variable>
struct VariableTraits {};

template <>
struct VariableTraits<CellVariable> { typedef CellData DataType; };

template <>
struct VariableTraits<NodeVariable> { typedef NodeData DataType; };

template <>
struct VariableTraits<SideVariable> { typedef SideData DataType; };

template <>
struct VariableTraits<IndicatorVariable> { typedef IndicatorData DataType; };

template <typename Variable>
std::shared_ptr<typename VariableTraits<Variable>::DataType>
getPatchData(
  SAMRAI::hier::Patch& patch,
  std::shared_ptr<Variable>& var,
  std::shared_ptr<SAMRAI::hier::VariableContext>& context)
{
  return SAMRAI_SHARED_PTR_CAST<typename VariableTraits<Variable>::DataType>(
    patch.getPatchData(var, context));
}

template <bool Ghosts>
SAMRAI::hier::Box getIndexBox(const CellData& data);

template <bool Ghosts>
SAMRAI::hier::Box getIndexBox(const NodeData& data);

template <bool Ghosts>
SAMRAI::hier::Box getIndexBox(const SideData& data, int side_normal);

template <bool Ghosts>
SAMRAI::hier::Box getIndexBox(const std::shared_ptr<CellData>& data)
{ return getIndexBox<Ghosts>(*data); }

template <bool Ghosts>
SAMRAI::hier::Box getIndexBox(const std::shared_ptr<NodeData>& data)
{ return getIndexBox<Ghosts>(*data); }

template <bool Ghosts>
SAMRAI::hier::Box getIndexBox(const std::shared_ptr<SideData>& data, int side_normal)
{ return getIndexBox<Ghosts>(*data, side_normal); }

}

#include "pdat/NodeInjection.h"
#include "geom/CartesianCellDoubleVolumeWeightedAverage.h"
#include "geom/CartesianCellDoubleMassWeightedAverage.h"
#include "geom/CartesianCellIntConstantCoarsen.h"
#include "geom/CartesianSideDoubleFirstOrderRefine.h"
#include "geom/CartesianCellDoubleConservativeLinearRefine.h"
#include "geom/CartesianNodeDoubleLinearRefine.h"

namespace clever {

// TODO Replace these by SAMRAI versions where possible
typedef clever::geom::CartesianCellDoubleVolumeWeightedAverage VolumeWeightedAverage;
typedef clever::geom::CartesianCellDoubleMassWeightedAverage MassWeightedAverage;
typedef clever::geom::CartesianCellIntConstantCoarsen ConstantCoarsen;

typedef clever::pdat::NodeInjection<Real> NodeInjection;
typedef clever::geom::CartesianNodeDoubleLinearRefine NodeLinearRefine;
typedef clever::geom::CartesianSideDoubleFirstOrderRefine SideFirstOrderRefine;
typedef clever::geom::CartesianCellDoubleConservativeLinearRefine CellConservativeLinearRefine;

}

#endif
