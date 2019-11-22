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
#include "Cleverleaf.h"

#include <iostream>
#include <cmath>

#include "SAMRAI/tbox/Collectives.h"
#include "SAMRAI/pdat/ForAll.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"


namespace clever {
namespace hydro {

const int Cleverleaf::g_rectangle;
const int Cleverleaf::g_circle;
const int Cleverleaf::g_point;

const Real Cleverleaf::g_small = 1.0e-16;
const Real Cleverleaf::g_big = 1.0e+21;

SAMRAI::tbox::StartupShutdownManager::Handler
Cleverleaf::s_initialize_handler(
  Cleverleaf::initializeCallback,
  0,
  0,
  Cleverleaf::finalizeCallback,
  SAMRAI::tbox::StartupShutdownManager::priorityTimers);

std::shared_ptr<SAMRAI::tbox::Timer> Cleverleaf::t_fill_boundary;

Cleverleaf::Cleverleaf(
  std::shared_ptr<SAMRAI::tbox::Database> input_database,
  std::shared_ptr<SAMRAI::hier::PatchHierarchy> hierarchy,
  const tbox::Dimension& dim,
  std::shared_ptr<SAMRAI::geom::CartesianGridGeometry> grid_geometry):
  LagrangianEulerianPatchStrategy(dim),
  d_dim(dim),
  d_nghosts(d_dim, 2),
  input_db(input_database),
  state_prefix("state"),
  d_velocity(new NodeVariable(d_dim, "velocity", dim.getValue())),
  d_massflux(new SideVariable(d_dim, "massflux", 1)),
  d_volflux(new SideVariable(d_dim, "volflux", 1)),
  d_pressure(new CellVariable(d_dim, "pressure", 1)),
  d_viscosity(new CellVariable(d_dim, "viscosity", 1)),
  d_soundspeed(new CellVariable(d_dim, "soundspeed", 1)),
  d_density(new CellVariable(d_dim, "density", 1)),
  d_energy(new CellVariable(d_dim, "energy", 1)),
  d_volume(new CellVariable(d_dim, "volume", 1)),
  d_celldeltas(new CellVariable(d_dim, "celldelta", dim.getValue())),
  d_cellcoords(new CellVariable(d_dim, "cellcoords", dim.getValue())),
  d_vertexdeltas(new NodeVariable(d_dim, "vertexdeltas",dim.getValue())),
  d_vertexcoords(new NodeVariable(d_dim, "vertexcoords", dim.getValue())),
  d_workarray1(new NodeVariable(d_dim, "workarray 1", 1)),
  d_workarray2(new NodeVariable(d_dim, "workarray 2", 1)),
  d_workarray3(new NodeVariable(d_dim, "workarray 3", 1)),
  d_workarray4(new NodeVariable(d_dim, "workarray 4", 1)),
  d_workarray5(new NodeVariable(d_dim, "workarray 5", 1)),
  d_workarray6(new NodeVariable(d_dim, "workarray 6", 1)),
  d_workarray7(new NodeVariable(d_dim, "workarray 7", 1)),
  d_tags(new TagVariable(d_dim, "tags", 1)),
  d_level_indicator(new IndicatorVariable(d_dim, "level_indicator", 1))
{
  d_hierarchy = hierarchy;
  d_grid_geometry = grid_geometry;

  /*
   * Add our coarsen operators to the registry.
   */
  std::shared_ptr<SAMRAI::hier::CoarsenOperator> vol_weighted_avg(
    new VolumeWeightedAverage(dim));
  std::shared_ptr<SAMRAI::hier::CoarsenOperator> mass_weighted_avg(
    new MassWeightedAverage(dim));
  std::shared_ptr<SAMRAI::hier::CoarsenOperator> constant_cell_coarsen(
    new ConstantCoarsen(dim));

  std::shared_ptr<SAMRAI::hier::CoarsenOperator> ndi(new NodeInjection());
  std::shared_ptr<SAMRAI::hier::RefineOperator> cndlr(new NodeLinearRefine());
  std::shared_ptr<SAMRAI::hier::RefineOperator> cedclr(
    new SideFirstOrderRefine());
  std::shared_ptr<SAMRAI::hier::RefineOperator> ccdclr(
    new CellConservativeLinearRefine());

  d_grid_geometry->addCoarsenOperator( typeid(CellVariable).name(),
                                       vol_weighted_avg);
  d_grid_geometry->addCoarsenOperator( typeid(CellVariable).name(),
                                       mass_weighted_avg);
  d_grid_geometry->addCoarsenOperator( typeid(IndicatorVariable).name(),
                                       constant_cell_coarsen);
  d_grid_geometry->addCoarsenOperator( typeid(NodeVariable).name(), ndi);
  d_grid_geometry->addRefineOperator( typeid(NodeVariable).name(), cndlr);
  d_grid_geometry->addRefineOperator( typeid(SideVariable).name(), cedclr);
  d_grid_geometry->addRefineOperator( typeid(CellVariable).name(), ccdclr);

  d_tag_all = input_database->getBoolWithDefault("tag_all", false);
  d_tag_q = input_database->getBoolWithDefault("tag_q", true);
  d_tag_density = input_database->getBoolWithDefault("tag_density", true);
  d_tag_energy = input_database->getBoolWithDefault("tag_energy", true);
  d_tag_pressure = input_database->getBoolWithDefault("tag_pressure", true);

  d_tag_q_threshold = input_database->getDoubleWithDefault(
    "tag_q_threshold", 0.001);
  d_tag_density_gradient = input_database->getDoubleWithDefault(
    "tag_density_threshold", 0.1);
  d_tag_energy_gradient = input_database->getDoubleWithDefault(
    "tag_energy_threshold", 0.1);
  d_tag_pressure_gradient = input_database->getDoubleWithDefault(
    "tag_pressure_threshold", 0.1);

  d_pdv_weight = input_database->getIntegerWithDefault(
    "physics_weight", 1);
}

void Cleverleaf::registerModelVariables(
  LagrangianEulerianLevelIntegrator* integrator)
{
  integrator->registerVariable(
    d_tags,
    LagrangianEulerianLevelIntegrator::NORMAL,
    LagrangianEulerianLevelIntegrator::NO_EXCH,
    d_nghosts,
    d_grid_geometry);

  integrator->registerVariable(
    d_velocity,
    LagrangianEulerianLevelIntegrator::FIELD,
    LagrangianEulerianLevelIntegrator::PRIME_CELLS_EXCH |
    LagrangianEulerianLevelIntegrator::PRE_LAGRANGE_EXCH |
    LagrangianEulerianLevelIntegrator::PRE_SWEEP_1_MOM_EXCH |
    LagrangianEulerianLevelIntegrator::PRE_SWEEP_2_MOM_EXCH,
    d_nghosts,
    d_grid_geometry);

  integrator->registerVariable(
    d_massflux,
    LagrangianEulerianLevelIntegrator::NORMAL,
    LagrangianEulerianLevelIntegrator::PRE_SWEEP_1_MOM_EXCH |
    LagrangianEulerianLevelIntegrator::PRE_SWEEP_2_MOM_EXCH,
    d_nghosts,
    d_grid_geometry);

  integrator->registerVariable(
    d_volflux,
    LagrangianEulerianLevelIntegrator::NORMAL,
    LagrangianEulerianLevelIntegrator::PRE_SWEEP_1_CELL_EXCH,
    d_nghosts,
    d_grid_geometry);

  integrator->registerVariable(
    d_pressure,
    LagrangianEulerianLevelIntegrator::NORMAL,
    LagrangianEulerianLevelIntegrator::PRIME_CELLS_EXCH |
    LagrangianEulerianLevelIntegrator::PRE_LAGRANGE_EXCH |
    LagrangianEulerianLevelIntegrator::HALF_STEP_EXCH,
    d_nghosts,
    d_grid_geometry);

  integrator->registerVariable(
    d_viscosity,
    LagrangianEulerianLevelIntegrator::NORMAL,
    LagrangianEulerianLevelIntegrator::PRIME_CELLS_EXCH |
    LagrangianEulerianLevelIntegrator::POST_VISCOSITY_EXCH,
    d_nghosts,
    d_grid_geometry);

  integrator->registerVariable(
    d_soundspeed,
    LagrangianEulerianLevelIntegrator::NORMAL,
    LagrangianEulerianLevelIntegrator::NO_EXCH,
    d_nghosts,
    d_grid_geometry);

  integrator->registerVariable(
    d_density,
    LagrangianEulerianLevelIntegrator::FIELD |
    LagrangianEulerianLevelIntegrator::REVERT,
    LagrangianEulerianLevelIntegrator::PRIME_CELLS_EXCH |
    LagrangianEulerianLevelIntegrator::PRE_LAGRANGE_EXCH |
    LagrangianEulerianLevelIntegrator::PRE_SWEEP_1_CELL_EXCH |
    LagrangianEulerianLevelIntegrator::PRE_SWEEP_1_MOM_EXCH |
    LagrangianEulerianLevelIntegrator::PRE_SWEEP_2_MOM_EXCH,
    d_nghosts,
    d_grid_geometry);

  integrator->registerVariable(
    d_energy,
    LagrangianEulerianLevelIntegrator::FIELD |
    LagrangianEulerianLevelIntegrator::REVERT,
    LagrangianEulerianLevelIntegrator::PRIME_CELLS_EXCH |
    LagrangianEulerianLevelIntegrator::PRE_LAGRANGE_EXCH |
    LagrangianEulerianLevelIntegrator::PRE_SWEEP_1_CELL_EXCH |
    LagrangianEulerianLevelIntegrator::PRE_SWEEP_1_MOM_EXCH |
    LagrangianEulerianLevelIntegrator::PRE_SWEEP_2_MOM_EXCH,
    d_nghosts,
    d_grid_geometry);

  integrator->registerVariable(
    d_volume,
    LagrangianEulerianLevelIntegrator::NORMAL,
    LagrangianEulerianLevelIntegrator::NO_EXCH,
    d_nghosts,
    d_grid_geometry);

  integrator->registerVariable(
    d_celldeltas,
    LagrangianEulerianLevelIntegrator::NORMAL,
    LagrangianEulerianLevelIntegrator::NO_EXCH,
    d_nghosts,
    d_grid_geometry);
  integrator->registerVariable(
    d_cellcoords,
    LagrangianEulerianLevelIntegrator::NORMAL,
    LagrangianEulerianLevelIntegrator::NO_EXCH,
    d_nghosts,
    d_grid_geometry);

  integrator->registerVariable(
    d_vertexdeltas,
    LagrangianEulerianLevelIntegrator::NORMAL,
    LagrangianEulerianLevelIntegrator::NO_EXCH,
    d_nghosts,
    d_grid_geometry);

  integrator->registerVariable(
    d_vertexcoords,
    LagrangianEulerianLevelIntegrator::NORMAL,
    LagrangianEulerianLevelIntegrator::NO_EXCH,
    d_nghosts,
    d_grid_geometry);

  integrator->registerVariable(
    d_workarray1,
    LagrangianEulerianLevelIntegrator::NORMAL,
    LagrangianEulerianLevelIntegrator::NO_EXCH,
    d_nghosts,
    d_grid_geometry);

  integrator->registerVariable(
    d_workarray2,
    LagrangianEulerianLevelIntegrator::NORMAL,
    LagrangianEulerianLevelIntegrator::NO_EXCH,
    d_nghosts,
    d_grid_geometry);

  integrator->registerVariable(
    d_workarray3,
    LagrangianEulerianLevelIntegrator::NORMAL,
    LagrangianEulerianLevelIntegrator::NO_EXCH,
    d_nghosts,
    d_grid_geometry);

  integrator->registerVariable(
    d_workarray4,
    LagrangianEulerianLevelIntegrator::NORMAL,
    LagrangianEulerianLevelIntegrator::NO_EXCH,
    d_nghosts,
    d_grid_geometry);

  integrator->registerVariable(
    d_workarray5,
    LagrangianEulerianLevelIntegrator::NORMAL,
    LagrangianEulerianLevelIntegrator::NO_EXCH,
    d_nghosts,
    d_grid_geometry);

  integrator->registerVariable(
    d_workarray6,
    LagrangianEulerianLevelIntegrator::NORMAL,
    LagrangianEulerianLevelIntegrator::NO_EXCH,
    d_nghosts,
    d_grid_geometry);

  integrator->registerVariable(
    d_workarray7,
    LagrangianEulerianLevelIntegrator::NORMAL,
    LagrangianEulerianLevelIntegrator::NO_EXCH,
    d_nghosts,
    d_grid_geometry);

  integrator->registerVariable(
    d_level_indicator,
    LagrangianEulerianLevelIntegrator::INDICATOR,
    LagrangianEulerianLevelIntegrator::NO_EXCH,
    d_nghosts,
    d_grid_geometry);

  hier::VariableDatabase* vardb = hier::VariableDatabase::getDatabase();

  d_plot_context = integrator->getPlotContext();

  if (d_visit_writer) {
    d_visit_writer->registerPlotQuantity(
      "Pressure",
      "SCALAR",
      vardb->mapVariableAndContextToIndex(
        d_pressure, d_plot_context));

    d_visit_writer->registerPlotQuantity(
      "Viscosity",
      "SCALAR",
      vardb->mapVariableAndContextToIndex(
        d_viscosity, d_plot_context));

    d_visit_writer->registerPlotQuantity(
      "Soundspeed",
      "SCALAR",
      vardb->mapVariableAndContextToIndex(
        d_soundspeed, d_plot_context));

    d_visit_writer->registerPlotQuantity(
      "Density",
      "SCALAR",
      vardb->mapVariableAndContextToIndex(
        d_density, d_plot_context));

    d_visit_writer->registerPlotQuantity(
      "Energy",
      "SCALAR",
      vardb->mapVariableAndContextToIndex(
        d_energy, d_plot_context));

    d_visit_writer->registerPlotQuantity(
      "Velocity",
      "VECTOR",
      vardb->mapVariableAndContextToIndex(
        d_velocity, d_plot_context));
  }
}

void Cleverleaf::registerVisItDataWriter(
  std::shared_ptr<SAMRAI::appu::VisItDataWriter> writer)
{
  d_visit_writer = writer;
}

void Cleverleaf::initializeDataOnPatch(
  SAMRAI::hier::Patch& patch,
  double init_data_time,
  bool initial_time)
{
  auto context = getCurrentDataContext();

  auto volume = SAMRAI::pdat::get_view<Dim, CellData>(getPatchData(patch, d_volume, context));

  std::shared_ptr<CellData> coord_data = getPatchData(patch, d_cellcoords, context);
  auto cellx = SAMRAI::pdat::get_view<Dim, CellData>(coord_data, Coord::X);
  auto celly = SAMRAI::pdat::get_view<Dim, CellData>(coord_data, Coord::Y);

  std::shared_ptr<CellData> delta_data = getPatchData(patch, d_celldeltas, context);
  auto celldx = SAMRAI::pdat::get_view<Dim, CellData>(delta_data, Coord::X);
  auto celldy = SAMRAI::pdat::get_view<Dim, CellData>(delta_data, Coord::Y);

  std::shared_ptr<NodeData> vertex_data = getPatchData(patch, d_vertexcoords, context);
  auto vertexx = SAMRAI::pdat::get_view<Dim, NodeData>(vertex_data, Coord::X);
  auto vertexy = SAMRAI::pdat::get_view<Dim, NodeData>(vertex_data, Coord::Y);

  std::shared_ptr<NodeData> vdelta_data = getPatchData(patch, d_vertexdeltas, context);
  auto vertexdx = SAMRAI::pdat::get_view<Dim, NodeData>(vdelta_data, Coord::X);
  auto vertexdy = SAMRAI::pdat::get_view<Dim, NodeData>(vdelta_data, Coord::Y);

  // std::cout << "patch.getBox().numberCells(0) + 4" << patch.getBox().numberCells(0) + 4 << std::endl;
  // std::cout << "patch.getBox().numberCells(1) + 4" << patch.getBox().numberCells(1) + 4 << std::endl;
  // std::cout << "static_cast<int>(patch.getBox().lower()[0]) - 2" << static_cast<int>(patch.getBox().lower()[0]) - 2 << std::endl;
  // std::cout << "static_cast<int>(patch.getBox().lower()[1]) - 2" << static_cast<int>(patch.getBox().lower()[1]) - 2 << std::endl;

  /*
   * Get the patch geometry - this stores the coordinates of the patch corner
   * in our Cartesian geometry, and the as well as the cell sizes.
   */
  const std::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> pgeom =
    SAMRAI_SHARED_PTR_CAST<SAMRAI::geom::CartesianPatchGeometry>(patch.getPatchGeometry());

  /* Array containing dx and dy */
  const double* dxs = pgeom->getDx();

  const Real dx = dxs[0];
  const Real dy = dxs[1];

  /* Lower left coordinate of the patch */
  const double* coords = pgeom->getXLower();

  const Real physical_xmin = coords[0];
  const Real physical_ymin = coords[1];

  // hier::Box node_box = patch.getBox();
  // node_box.grow(d_nghosts);
  // node_box.growUpper(hier::IntVector::getOne(d_dim));
  const SAMRAI::hier::Box& node_ghost_box = getIndexBox<true>(vertex_data);
  const SAMRAI::hier::Box& cell_ghost_box = getIndexBox<true>(coord_data);

  const int j_min = patch.getBox().lower(0);
  const int k_min = patch.getBox().lower(1);

  // SAMRAI::hier::Box tmp_box = node_ghost_box;
  // tmp_box.shorten(0, 2);
  // tmp_box.shorten(1, 2);

  FORCEINLINE_MACRO
  SAMRAI::pdat::parallel_for_all(node_ghost_box, [=] SAMRAI_HOST_DEVICE (int k, int j) {
    vertexx(j,k) = physical_xmin+(dx*(j-j_min));
    vertexdx(j,k) = dx;
    vertexy(j,k) = physical_ymin+(dy*(k-k_min));
    vertexdy(j,k) = dy;
  });

  // hier::Box cell_box = patch.getBox();
  // cell_box.grow(d_nghosts);

  FORCEINLINE_MACRO
  SAMRAI::pdat::parallel_for_all(cell_ghost_box, [=] SAMRAI_HOST_DEVICE (int k, int j) {
    cellx(j,k) = 0.5*(vertexx(j,k)+vertexx(j+1,k));
    celldx(j,k) = dx;
    celly(j,k) = 0.5*(vertexy(j,k)+vertexy(j,k+1));
    celldy(j,k) = dy;
    volume(j,k) = dx*dy;
  });

  if (initial_time) {
    std::shared_ptr<NodeData> velocity_data = getPatchData(patch, d_velocity, context);
    auto xvel0 = SAMRAI::pdat::get_view<Dim, NodeData>(velocity_data, Coord::X);
    auto yvel0 = SAMRAI::pdat::get_view<Dim, NodeData>(velocity_data, Coord::Y);

    auto pressure = SAMRAI::pdat::get_view<Dim, CellData>(getPatchData(patch, d_pressure, context));
    auto viscosity = SAMRAI::pdat::get_view<Dim, CellData>(getPatchData(patch, d_viscosity, context));
    auto soundspeed = SAMRAI::pdat::get_view<Dim, CellData>(getPatchData(patch, d_soundspeed, context));
    auto density = SAMRAI::pdat::get_view<Dim, CellData>(getPatchData(patch, d_density, context));
    auto energy = SAMRAI::pdat::get_view<Dim, CellData>(getPatchData(patch, d_energy, context));

    std::shared_ptr<SideData> mass_flux_data = getPatchData(patch, d_massflux, context);
    auto mass_flux_x = SAMRAI::pdat::get_view<Dim, SideData>(mass_flux_data, Side::X);
    auto mass_flux_y = SAMRAI::pdat::get_view<Dim, SideData>(mass_flux_data, Side::Y);

    std::shared_ptr<SideData> vol_flux_data = getPatchData(patch, d_volflux, context);
    auto vol_flux_x = SAMRAI::pdat::get_view<Dim, SideData>(getPatchData(patch, d_volflux, context), Side::X);
    auto vol_flux_y = SAMRAI::pdat::get_view<Dim, SideData>(getPatchData(patch, d_volflux, context), Side::Y);

  FORCEINLINE_MACRO
    SAMRAI::pdat::parallel_for_all(node_ghost_box, [=] SAMRAI_HOST_DEVICE (int k, int j) {
      xvel0(j,k) = 0.0;
      yvel0(j,k) = 0.0;
    });

  FORCEINLINE_MACRO
    SAMRAI::pdat::parallel_for_all(cell_ghost_box, [=] SAMRAI_HOST_DEVICE (int k, int j) {
      viscosity(j,k) = 0.0;
      soundspeed(j,k) = 0.0;
      pressure(j,k) = 0.0;
      energy(j,k) = 0.0;
      density(j,k) = 0.0;
    });

    // hier::Box side_x_box = patch.getBox();
    // side_x_box.grow(d_nghosts);
    // side_x_box.growUpper(0,1);

  FORCEINLINE_MACRO
    SAMRAI::pdat::parallel_for_all(getIndexBox<true>(mass_flux_data, Side::X), [=] SAMRAI_HOST_DEVICE (int k, int j) {
      //RAJA::forallN<PatchPolicy>(isets->getCells(1,true), isets->getNodes(0,true), [=] SAMRAI_HOST_DEVICE (int k, int j) {
      mass_flux_x(j,k) = 0.0;
      vol_flux_x(j,k) = 0.0;
    });

    // hier::Box side_y_box = patch.getBox();
    // side_y_box.grow(d_nghosts);
    // side_y_box.growUpper(1,1);

  FORCEINLINE_MACRO
    SAMRAI::pdat::parallel_for_all(getIndexBox<true>(mass_flux_data, Side::Y), [=] SAMRAI_HOST_DEVICE (int k, int j) {
      //RAJA::forallN<PatchPolicy>(isets->getNodes(1,true), isets->getCells(0,true), [=] SAMRAI_HOST_DEVICE (int k, int j) {
      mass_flux_y(j,k) = 0.0;
      vol_flux_y(j,k) = 0.0;
    });

    SAMRAI::tbox::parallel_synchronize();

    std::shared_ptr<tbox::Database> states_db = input_db->getDatabase("states");

    const int number_of_states = states_db->getInteger("num_states");

    Real* state_density  = (Real *) clever::allocate(sizeof(Real)*number_of_states);
    Real* state_energy   = (Real *) clever::allocate(sizeof(Real)*number_of_states);
    Real* state_xvel     = (Real *) clever::allocate(sizeof(Real)*number_of_states);
    Real* state_yvel     = (Real *) clever::allocate(sizeof(Real)*number_of_states);
    Real* state_xmin     = (Real *) clever::allocate(sizeof(Real)*number_of_states);
    Real* state_ymin     = (Real *) clever::allocate(sizeof(Real)*number_of_states);
    Real* state_xmax     = (Real *) clever::allocate(sizeof(Real)*number_of_states);
    Real* state_ymax     = (Real *) clever::allocate(sizeof(Real)*number_of_states);
    Real* state_radius   = (Real *) clever::allocate(sizeof(Real)*number_of_states);
    int* state_geometry  = (int *)  clever::allocate(sizeof(int) *number_of_states);

    for(int state = 0; state < number_of_states; state++) {
      std::ostringstream state_stream;
      state_stream << state_prefix << state;

      std::shared_ptr<tbox::Database> current_state = states_db->getDatabase(
        state_stream.str());

      const std::string state_geometry_string = current_state->getStringWithDefault("geometry", "RECTANGLE");

      if (state_geometry_string.compare("RECTANGLE") == 0) {
        state_geometry[state] = g_rectangle;
      } else if (state_geometry_string.compare("CIRCLE") == 0) {
        state_geometry[state] = g_circle;
      } else if (state_geometry_string.compare("POINT") == 0) {
        state_geometry[state] = g_point;
      }
      else {
        SAMRAI::tbox::perr << "[ERROR] An error occurred.... " << std::endl;
        exit(-1);
      }

      if(state_geometry[state] == g_circle ||
         state_geometry[state] == g_point) {
        double center[2];

        current_state->getDoubleArray("center", center, 2);

        state_xmin[state] = center[0];
        state_ymin[state] = center[1];

        state_xmax[state] = -1;
        state_ymax[state] = -1;
      } else {
        if(state == 0) {
          state_xmin[state] = -1;
          state_ymin[state] = -1;

          state_xmax[state] = -1;
          state_ymax[state] = -1;
        } else {
          double state_min[2];
          double state_max[2];

          current_state->getDoubleArray("min", state_min, 2);
          current_state->getDoubleArray("max", state_max, 2);

          state_xmin[state] = state_min[0];
          state_ymin[state] = state_min[1];

          state_xmax[state] = state_max[0];
          state_ymax[state] = state_max[1];
        }
      }

      state_density[state] = current_state->getDouble("density");
      state_energy[state] = current_state->getDouble("energy");
      state_xvel[state] = current_state->getDoubleWithDefault("xvel", 0.0);
      state_yvel[state] = current_state->getDoubleWithDefault("yvel", 0.0);
      state_radius[state] = current_state->getDoubleWithDefault("radius", -1.0);
    }

  FORCEINLINE_MACRO
    SAMRAI::pdat::parallel_for_all(cell_ghost_box,  [=] SAMRAI_HOST_DEVICE (int k, int j) {
      //RAJA::forallN<PatchPolicy>(isets->getCells(1,true), isets->getCells(0,true), [=] SAMRAI_HOST_DEVICE (int k, int j) {
      energy(j,k) = state_energy[0];
      density(j,k) = state_density[0];
    });

  FORCEINLINE_MACRO
    SAMRAI::pdat::parallel_for_all(node_ghost_box, [=] SAMRAI_HOST_DEVICE (int k, int j) {
      //RAJA::forallN<PatchPolicy>(isets->getNodes(1,true), isets->getNodes(0,true), [=] SAMRAI_HOST_DEVICE (int k, int j) {
      xvel0(j,k) = state_xvel[0];
      yvel0(j,k) = state_yvel[0];
    });

    for (int state = 1; state < number_of_states; ++state) {
      SAMRAI::tbox::parallel_synchronize();

      switch (state_geometry[state]) {
      case g_rectangle:
  FORCEINLINE_MACRO
        SAMRAI::pdat::parallel_for_all(cell_ghost_box,  [=] SAMRAI_HOST_DEVICE (int k, int j) {
          if (vertexx(j+1,k) >= state_xmin[state] && vertexx(j,k) < state_xmax[state]) {
            if (vertexy(j,k+1) >= state_ymin[state] && vertexy(j,k) < state_ymax[state]) {
              energy(j,k)=state_energy[state];
              density(j,k)=state_density[state];
              xvel0(j,k)=state_xvel[state];
              xvel0(j+1,k)=state_xvel[state];
              xvel0(j,k+1)=state_xvel[state];
              xvel0(j+1,k+1)=state_xvel[state];
              yvel0(j,k)=state_yvel[state];
              yvel0(j+1,k)=state_yvel[state];
              yvel0(j,k+1)=state_yvel[state];
              yvel0(j+1,k+1)=state_yvel[state];
            }
          }
        }); break;
      case g_circle: {
        const Real x_cent=state_xmin[state];
        const Real y_cent=state_ymin[state];

  FORCEINLINE_MACRO
        SAMRAI::pdat::parallel_for_all(cell_ghost_box,  [=] SAMRAI_HOST_DEVICE (int k, int j) {
          const Real radius = sqrt(
            (cellx(j,k)-x_cent)*(cellx(j,k)-x_cent)
            +(celly(j,k)-y_cent)*(celly(j,k)-y_cent));

          if (radius <= state_radius[state]) {
            energy(j,k)=state_energy[state];
            density(j,k)=state_density[state];
            xvel0(j,k)=state_xvel[state];
            xvel0(j+1,k)=state_xvel[state];
            xvel0(j,k+1)=state_xvel[state];
            xvel0(j+1,k+1)=state_xvel[state];
            yvel0(j,k)=state_yvel[state];
            yvel0(j+1,k)=state_yvel[state];
            yvel0(j,k+1)=state_yvel[state];
            yvel0(j+1,k+1)=state_yvel[state];
          }
        }); } break;
      case g_point: {
        const Real x_cent=state_xmin[state];
        const Real y_cent=state_ymin[state];

  FORCEINLINE_MACRO
        SAMRAI::pdat::parallel_for_all(cell_ghost_box,  [=] SAMRAI_HOST_DEVICE (int k, int j) {
          if (vertexx(j,k) == x_cent && vertexy(j,k) == y_cent) {
            energy(j,k)=state_energy[state];
            density(j,k)=state_density[state];
            xvel0(j,k)=state_xvel[state];
            xvel0(j+1,k)=state_xvel[state];
            xvel0(j,k+1)=state_xvel[state];
            xvel0(j+1,k+1)=state_xvel[state];
            yvel0(j,k)=state_yvel[state];
            yvel0(j+1,k)=state_yvel[state];
            yvel0(j,k+1)=state_yvel[state];
            yvel0(j+1,k+1)=state_yvel[state];
          }
        }); } break;
      }
    }

    clever::free(state_density);
    clever::free(state_energy);
    clever::free(state_xvel);
    clever::free(state_yvel);
    clever::free(state_xmin);
    clever::free(state_ymin);
    clever::free(state_xmax);
    clever::free(state_ymax);
    clever::free(state_radius);
    clever::free(state_geometry);
  }

  ideal_gas_knl(patch, false);
}

void Cleverleaf::accelerate(SAMRAI::hier::Patch& patch, double dt)
{
  auto current_context = getCurrentDataContext();
  auto new_context = getNewDataContext();

  auto density0 = SAMRAI::pdat::get_view<Dim, CellData>(getPatchData(patch, d_density, current_context));
  auto volume = SAMRAI::pdat::get_view<Dim, CellData>(getPatchData(patch, d_volume, current_context));
  auto pressure = SAMRAI::pdat::get_view<Dim, CellData>(getPatchData(patch, d_pressure, current_context));
  auto viscosity = SAMRAI::pdat::get_view<Dim, CellData>(getPatchData(patch, d_viscosity, current_context));

  std::shared_ptr<CellData> delta_data = getPatchData(patch, d_celldeltas, current_context);
  auto yarea = SAMRAI::pdat::get_view<Dim, CellData>(delta_data, Coord::X);
  auto xarea = SAMRAI::pdat::get_view<Dim, CellData>(delta_data, Coord::Y);

  std::shared_ptr<NodeData> current_velocity_data = getPatchData(patch, d_velocity, current_context);
  std::shared_ptr<NodeData> new_velocity_data = getPatchData(patch, d_velocity, new_context);
  auto xvel0 = SAMRAI::pdat::get_view<Dim, NodeData>(current_velocity_data, Coord::X);
  auto xvel1 = SAMRAI::pdat::get_view<Dim, NodeData>(new_velocity_data, Coord::X);

  auto yvel0 = SAMRAI::pdat::get_view<Dim, NodeData>(current_velocity_data, Coord::Y);
  auto yvel1 = SAMRAI::pdat::get_view<Dim, NodeData>(new_velocity_data, Coord::Y);

  auto stepbymass = SAMRAI::pdat::get_view<Dim, NodeData>(getPatchData(patch, d_workarray1, current_context));

  const SAMRAI::hier::Box &node_box = getIndexBox<false>(current_velocity_data);
  // hier::Box node_box = patch.getBox();
  // node_box.growUpper(hier::IntVector::getOne(d_dim));

  FORCEINLINE_MACRO
  SAMRAI::pdat::parallel_for_all(node_box,  [=] SAMRAI_HOST_DEVICE (int k, int j) {
    const Real nodal_mass = (density0(j-1,k-1)*volume(j-1,k-1)
                             + density0(j,k-1)*volume(j,k-1)
                             + density0(j,k)*volume(j,k)
                             + density0(j-1,k)*volume(j-1,k))*0.25;

    stepbymass(j,k)=(0.5*dt)/nodal_mass;
  });

  FORCEINLINE_MACRO
  SAMRAI::pdat::parallel_for_all(node_box,  [=] SAMRAI_HOST_DEVICE (int k, int j) {
    xvel1(j,k)=xvel0(j,k)-stepbymass(j,k)*(xarea(j  ,k  )*(pressure(j  ,k  )-pressure(j-1,k  ))
                                           +xarea(j  ,k-1)*(pressure(j  ,k-1)-pressure(j-1,k-1)));
    yvel1(j,k)=yvel0(j,k)-stepbymass(j,k)*(yarea(j  ,k  )*(pressure(j  ,k  )-pressure(j  ,k-1))
                                           +yarea(j-1,k  )*(pressure(j-1,k  )-pressure(j-1,k-1)));
  });

  FORCEINLINE_MACRO
  SAMRAI::pdat::parallel_for_all(node_box,  [=] SAMRAI_HOST_DEVICE (int k, int j) {
    xvel1(j,k)=xvel1(j,k)-stepbymass(j,k)*(xarea(j  ,k  )*(viscosity(j  ,k  )-viscosity(j-1,k  ))
                                           +xarea(j  ,k-1)*(viscosity(j  ,k-1)-viscosity(j-1,k-1)));
    yvel1(j,k)=yvel1(j,k)-stepbymass(j,k)*(yarea(j  ,k  )*(viscosity(j  ,k  )-viscosity(j  ,k-1))
                                           +yarea(j-1,k  )*(viscosity(j-1,k  )-viscosity(j-1,k-1)));
  });
}

void Cleverleaf::ideal_gas_knl(hier::Patch& patch, const bool predict)
{
  auto current_context = getCurrentDataContext();
  auto new_context = getNewDataContext();

  std::shared_ptr<CellData> density_data =
    predict ? getPatchData(patch, d_density, new_context) : getPatchData(patch, d_density, current_context);
  auto density = SAMRAI::pdat::get_view<Dim, CellData>(density_data);

  std::shared_ptr<CellData> energy_data =
    predict ? getPatchData(patch, d_energy, new_context) : getPatchData(patch, d_energy, current_context);
  auto energy = SAMRAI::pdat::get_view<Dim, CellData>(energy_data);

  auto pressure = SAMRAI::pdat::get_view<Dim, CellData>(getPatchData(patch, d_pressure, current_context));
  auto soundspeed = SAMRAI::pdat::get_view<Dim, CellData>(getPatchData(patch, d_soundspeed, current_context));

  FORCEINLINE_MACRO
  SAMRAI::pdat::parallel_for_all(patch.getBox(),  [=] SAMRAI_HOST_DEVICE (int k, int j) {
    const Real v=1.0/density(j,k);
    pressure(j,k)=(1.4-1.0)*density(j,k)*energy(j,k);
    const Real pressurebyenergy=(1.4-1.0)*density(j,k);
    const Real pressurebyvolume=-density(j,k)*pressure(j,k);
    const Real sound_speed_squared=v*v*(pressure(j,k)*pressurebyenergy-pressurebyvolume);
    soundspeed(j,k)=sqrt(sound_speed_squared);
  });
}

void Cleverleaf::viscosity_knl(SAMRAI::hier::Patch& patch)
{
  auto context = getCurrentDataContext();
  auto density0 = SAMRAI::pdat::get_view<Dim, CellData>(getPatchData(patch, d_density, context));
  auto pressure = SAMRAI::pdat::get_view<Dim, CellData>(getPatchData(patch, d_pressure, context));
  auto viscosity = SAMRAI::pdat::get_view<Dim, CellData>(getPatchData(patch, d_viscosity, context));

  std::shared_ptr<NodeData> velocity_data = getPatchData(patch, d_velocity, context);
  auto xvel0 = SAMRAI::pdat::get_view<Dim, NodeData>(velocity_data, Coord::X);
  auto yvel0 = SAMRAI::pdat::get_view<Dim, NodeData>(velocity_data, Coord::Y);

  std::shared_ptr<CellData> delta_data = getPatchData(patch, d_celldeltas, context);
  auto celldx = SAMRAI::pdat::get_view<Dim, CellData>(delta_data, Coord::X);
  auto celldy = SAMRAI::pdat::get_view<Dim, CellData>(delta_data, Coord::Y);

  FORCEINLINE_MACRO
  SAMRAI::pdat::parallel_for_all(patch.getBox(),  [=] SAMRAI_HOST_DEVICE (int k, int j) {
    const Real ugrad=(xvel0(j+1,k  )+xvel0(j+1,k+1))-(xvel0(j  ,k  )+xvel0(j  ,k+1));
    const Real vgrad=(yvel0(j  ,k+1)+yvel0(j+1,k+1))-(yvel0(j  ,k  )+yvel0(j+1,k  ));

    const Real div = (celldx(j,k)*(ugrad)+  celldy(j,k)*(vgrad));

    const Real strain2 = 0.5*(xvel0(j,  k+1) + xvel0(j+1,k+1)-xvel0(j  ,k  )-xvel0(j+1,k  ))/celldy(j,k)
      + 0.5*(yvel0(j+1,k  ) + yvel0(j+1,k+1)-yvel0(j  ,k  )-yvel0(j  ,k+1))/celldx(j,k);

    Real pgradx=(pressure(j+1,k)-pressure(j-1,k))/(celldx(j,k)+celldx(j+1,k));
    Real pgrady=(pressure(j,k+1)-pressure(j,k-1))/(celldy(j,k)+celldy(j,k+1));

    const Real pgradx2 = pgradx*pgradx;
    const Real pgrady2 = pgrady*pgrady;

    const Real limiter = ((0.5*(ugrad)/celldx(j,k))*pgradx2+(0.5*(vgrad)/celldy(j,k))*pgrady2+strain2*pgradx*pgrady)
      /MAX(pgradx2+pgrady2,1.0e-16);
    if (limiter > 0.0 || div >= 0.0) {
      viscosity(j,k) = 0.0;
    } else {
      const Real dirx = pgradx < 0.0 ? 1.0 : -1.0;
      pgradx = dirx*MAX(1.0e-16,fabs(pgradx));
      const Real diry= pgrady < 0.0 ? 1.0 : -1.0;
      pgrady = diry*MAX(1.0e-16,fabs(pgrady));

      const Real pgrad = sqrt(pgradx*pgradx+pgrady*pgrady);
      const Real xgrad = fabs(celldx(j,k)*pgrad/pgradx);
      const Real ygrad = fabs(celldy(j,k)*pgrad/pgrady);
      const Real grad  = MIN(xgrad,ygrad);
      const Real grad2 = grad*grad;
      viscosity(j,k)=2.0*density0(j,k)*grad2*limiter*limiter;
    }
  });
}

double Cleverleaf::calc_dt_knl(SAMRAI::hier::Patch& patch)
{
  auto context = getCurrentDataContext();

  auto density0 = SAMRAI::pdat::get_view<Dim, CellData>(getPatchData(patch, d_density, context));
  auto energy0 = SAMRAI::pdat::get_view<Dim, CellData>(getPatchData(patch, d_energy, context));
  auto pressure = SAMRAI::pdat::get_view<Dim, CellData>(getPatchData(patch, d_pressure, context));
  auto soundspeed = SAMRAI::pdat::get_view<Dim, CellData>(getPatchData(patch, d_soundspeed, context));
  auto viscosity = SAMRAI::pdat::get_view<Dim, CellData>(getPatchData(patch, d_viscosity, context));

  std::shared_ptr<NodeData> velocity_data = getPatchData(patch, d_velocity, context);
  auto xvel0 = SAMRAI::pdat::get_view<Dim, NodeData>(velocity_data, Coord::X);
  auto yvel0 = SAMRAI::pdat::get_view<Dim, NodeData>(velocity_data, Coord::Y);

  std::shared_ptr<CellData> delta_data = getPatchData(patch, d_celldeltas, context);
  auto celldx = SAMRAI::pdat::get_view<Dim, CellData>(delta_data, Coord::X);
  auto celldy = SAMRAI::pdat::get_view<Dim, CellData>(delta_data, Coord::Y);

  auto volume = SAMRAI::pdat::get_view<Dim, CellData>(getPatchData(patch, d_volume, context));

  auto xarea = SAMRAI::pdat::get_view<Dim, CellData>(delta_data, Coord::Y);
  auto yarea = SAMRAI::pdat::get_view<Dim, CellData>(delta_data, Coord::X);

  Real dt_min_val = g_big;
  const Real dtc_safe = 0.7;
  const Real dtu_safe = 0.5;
  const Real dtv_safe = 0.5;
  const Real dtdiv_safe = 0.7;
  const Real dt_min = 0.0000001;

  SAMRAI::tbox::parallel_reduction_variable_t<tbox::Reduction::Min, Real> min_dt(1.0e+20);

  FORCEINLINE_MACRO
  SAMRAI::pdat::parallel_for_all(patch.getBox(),  [=] SAMRAI_HOST_DEVICE (int k, int j) {
    const Real dsx=celldx(j,k);
    const Real dsy=celldy(j,k);

    Real cc=soundspeed(j,k)*soundspeed(j,k);
    cc=cc+2.0*viscosity(j,k)/density0(j,k);
    cc=MAX(sqrt(cc),g_small);

    const Real dtct=dtc_safe*MIN(dsx,dsy)/cc;

    Real div=0.0;

    Real dv1=(xvel0(j  ,k)+xvel0(j  ,k+1))*xarea(j  ,k);
    Real dv2=(xvel0(j+1,k)+xvel0(j+1,k+1))*xarea(j+1,k);

    div=div+dv2-dv1;

    const Real dtut=dtu_safe*2.0*volume(j,k)/MAX(fabs(dv1),MAX(fabs(dv2),g_small*volume(j,k)));

    dv1=(yvel0(j,k  )+yvel0(j+1,k  ))*yarea(j,k  );
    dv2=(yvel0(j,k+1)+yvel0(j+1,k+1))*yarea(j,k+1);

    div=div+dv2-dv1;

    const Real dtvt=dtv_safe*2.0*volume(j,k)/MAX(fabs(dv1),MAX(fabs(dv2),g_small*volume(j,k)));

    div=div/(2.0*volume(j,k));

    const Real dtdivt = (div < -g_small) ? dtdiv_safe*(-1.0/div) : g_big;

    min_dt.min(MIN(dtct,MIN(dtut,MIN(dtvt,dtdivt))));
  });

  dt_min_val = Real(min_dt);

  if (dt_min_val < dt_min) {
    SAMRAI::tbox::perr << "SMALL TIMESTEP: " << dt_min_val << " we should die here..." << std::endl;
  }

  return dt_min_val;
}

void Cleverleaf::pdv_knl(SAMRAI::hier::Patch& patch, const double dt, const bool predict)
{
  auto current_context = getCurrentDataContext();
  auto new_context = getNewDataContext();

  auto volume = SAMRAI::pdat::get_view<Dim, CellData>(getPatchData(patch, d_volume, current_context));
  auto density0 = SAMRAI::pdat::get_view<Dim, CellData>(getPatchData(patch, d_density, current_context));
  auto density1 = SAMRAI::pdat::get_view<Dim, CellData>(getPatchData(patch, d_density, new_context));

  auto energy0 = SAMRAI::pdat::get_view<Dim, CellData>(getPatchData(patch, d_energy, current_context));
  auto energy1 = SAMRAI::pdat::get_view<Dim, CellData>(getPatchData(patch, d_energy, new_context));

  auto pressure = SAMRAI::pdat::get_view<Dim, CellData>(getPatchData(patch, d_pressure, current_context));
  auto viscosity = SAMRAI::pdat::get_view<Dim, CellData>(getPatchData(patch, d_viscosity, current_context));

  std::shared_ptr<CellData> delta_data = getPatchData(patch, d_celldeltas, current_context);
  auto xarea = SAMRAI::pdat::get_view<Dim, CellData>(delta_data, Coord::Y);
  auto yarea = SAMRAI::pdat::get_view<Dim, CellData>(delta_data, Coord::X);

  std::shared_ptr<NodeData> current_velocity_data = getPatchData(patch, d_velocity, current_context);
  std::shared_ptr<NodeData> new_velocity_data = getPatchData(patch, d_velocity, new_context);
  auto xvel0 = SAMRAI::pdat::get_view<Dim, NodeData>(current_velocity_data, Coord::X);
  auto yvel0 = SAMRAI::pdat::get_view<Dim, NodeData>(current_velocity_data, Coord::Y);

  auto xvel1 = SAMRAI::pdat::get_view<Dim, NodeData>(new_velocity_data, Coord::X);
  auto yvel1 = SAMRAI::pdat::get_view<Dim, NodeData>(new_velocity_data, Coord::Y);

  auto volume_change = SAMRAI::pdat::get_view<Dim, NodeData>(getPatchData(patch, d_workarray1, current_context));

  for (int i = 0; i < d_pdv_weight; i++) {
    if (predict) {
  FORCEINLINE_MACRO
            SAMRAI::pdat::parallel_for_all(patch.getBox(),  [=] SAMRAI_HOST_DEVICE (int k, int j) {
        const Real left_flux=  (xarea(j  ,k  )*(xvel0(j  ,k  )+xvel0(j  ,k+1)
                                                +xvel0(j  ,k  )+xvel0(j  ,k+1)))*0.25*dt*0.5;
        const Real right_flux= (xarea(j+1,k  )*(xvel0(j+1,k  )+xvel0(j+1,k+1)
                                                +xvel0(j+1,k  )+xvel0(j+1,k+1)))*0.25*dt*0.5;
        const Real bottom_flux=(yarea(j  ,k  )*(yvel0(j  ,k  )+yvel0(j+1,k  )
                                                +yvel0(j  ,k  )+yvel0(j+1,k  )))*0.25*dt*0.5;
        const Real top_flux=   (yarea(j  ,k+1)*(yvel0(j  ,k+1)+yvel0(j+1,k+1)
                                                +yvel0(j  ,k+1)+yvel0(j+1,k+1)))*0.25*dt*0.5;
        const Real total_flux=right_flux-left_flux+top_flux-bottom_flux;

        volume_change(j,k)=volume(j,k)/(volume(j,k)+total_flux);

        const Real min_cell_volume=MIN(
          volume(j,k)+right_flux-left_flux+top_flux-bottom_flux,
          MIN(volume(j,k)+right_flux-left_flux,volume(j,k)+top_flux-bottom_flux));

        const Real recip_volume=1.0/volume(j,k);
        const Real energy_change=(pressure(j,k)/density0(j,k)+viscosity(j,k)/density0(j,k))*total_flux*recip_volume;

        energy1(j,k)=energy0(j,k)-energy_change;
        density1(j,k)=density0(j,k)*volume_change(j,k);
      });
    } else {
  FORCEINLINE_MACRO
            SAMRAI::pdat::parallel_for_all(patch.getBox(),  [=] SAMRAI_HOST_DEVICE (int k, int j) {
        const Real left_flux=  (xarea(j  ,k  )*(xvel0(j  ,k  )+xvel0(j  ,k+1)
                                                +xvel1(j  ,k  )+xvel1(j  ,k+1)))*0.25*dt;
        const Real right_flux= (xarea(j+1,k  )*(xvel0(j+1,k  )+xvel0(j+1,k+1)
                                                +xvel1(j+1,k  )+xvel1(j+1,k+1)))*0.25*dt;
        const Real bottom_flux=(yarea(j  ,k  )*(yvel0(j  ,k  )+yvel0(j+1,k  )
                                                +yvel1(j  ,k  )+yvel1(j+1,k  )))*0.25*dt;
        const Real top_flux=   (yarea(j  ,k+1)*(yvel0(j  ,k+1)+yvel0(j+1,k+1)
                                                +yvel1(j  ,k+1)+yvel1(j+1,k+1)))*0.25*dt;
        const Real total_flux=right_flux-left_flux+top_flux-bottom_flux;

        volume_change(j,k)=volume(j,k)/(volume(j,k)+total_flux);

        const Real min_cell_volume=MIN(
          volume(j,k)+right_flux-left_flux+top_flux-bottom_flux,
          MIN(volume(j,k)+right_flux-left_flux,volume(j,k)+top_flux-bottom_flux));

        const Real recip_volume=1.0/volume(j,k) ;

        const Real energy_change=(pressure(j,k)/density0(j,k)+viscosity(j,k)/density0(j,k))*total_flux*recip_volume;

        energy1(j,k)=energy0(j,k)-energy_change;
        density1(j,k)=density0(j,k)*volume_change(j,k);
      });
    }
  }
}

void Cleverleaf::flux_calc_knl(SAMRAI::hier::Patch& patch, double dt)
{
  auto current_context = getCurrentDataContext();
  auto new_context = getNewDataContext();

  std::shared_ptr<NodeData> current_velocity_data = getPatchData(patch, d_velocity, current_context);
  std::shared_ptr<NodeData> new_velocity_data = getPatchData(patch, d_velocity, new_context);
  auto xvel0 = SAMRAI::pdat::get_view<Dim, NodeData>(current_velocity_data, Coord::X);
  auto yvel0 = SAMRAI::pdat::get_view<Dim, NodeData>(current_velocity_data, Coord::Y);
  auto xvel1 = SAMRAI::pdat::get_view<Dim, NodeData>(new_velocity_data, Coord::X);
  auto yvel1 = SAMRAI::pdat::get_view<Dim, NodeData>(new_velocity_data, Coord::Y);

  std::shared_ptr<SideData> volflux_data = getPatchData(patch, d_volflux, current_context);
  auto vol_flux_x = SAMRAI::pdat::get_view<Dim, SideData>(volflux_data, Side::X);
  auto vol_flux_y = SAMRAI::pdat::get_view<Dim, SideData>(volflux_data, Side::Y);

  std::shared_ptr<CellData> delta_data = getPatchData(patch, d_celldeltas, current_context);
  auto xarea = SAMRAI::pdat::get_view<Dim, CellData>(delta_data, Coord::Y);
  auto yarea = SAMRAI::pdat::get_view<Dim, CellData>(delta_data, Coord::X);

  // SAMRAI::hier::Box side_x_box = patch.getBox();
  // side_x_box.growUpper(0,1);

  FORCEINLINE_MACRO
  SAMRAI::pdat::parallel_for_all(getIndexBox<false>(volflux_data, Side::X),  [=] SAMRAI_HOST_DEVICE (int k, int j) {
    vol_flux_x(j,k) = 0.25*dt*xarea(j,k)
      *(xvel0(j,k)+xvel0(j,k+1)+xvel1(j,k)+xvel1(j,k+1));
  });

  // SAMRAI::hier::Box side_y_box = patch.getBox();
  // side_x_box.growUpper(1,1);

  FORCEINLINE_MACRO
  SAMRAI::pdat::parallel_for_all(getIndexBox<false>(volflux_data, Side::Y),  [=] SAMRAI_HOST_DEVICE (int k, int j) {
    vol_flux_y(j,k)=0.25*dt*yarea(j,k)*(yvel0(j,k)+yvel0(j+1,k)+yvel1(j,k)+yvel1(j+1,k));
  });
}

void Cleverleaf::advec_cell(SAMRAI::hier::Patch& patch, int sweep_number, const int dir)
{
  auto current_context = getCurrentDataContext();
  auto new_context = getNewDataContext();

  std::shared_ptr<CellData> volume_data = getPatchData(patch, d_volume, current_context);
  auto volume = SAMRAI::pdat::get_view<Dim, CellData>(volume_data);
  auto density1 = SAMRAI::pdat::get_view<Dim, CellData>(getPatchData(patch, d_density, new_context));
  auto energy1 = SAMRAI::pdat::get_view<Dim, CellData>(getPatchData(patch, d_energy, new_context));

  std::shared_ptr<SideData> vol_flux_data = getPatchData(patch, d_volflux, current_context);
  auto vol_flux_x = SAMRAI::pdat::get_view<Dim, SideData>(vol_flux_data, Side::X);
  auto vol_flux_y = SAMRAI::pdat::get_view<Dim, SideData>(vol_flux_data, Side::Y);

  std::shared_ptr<SideData> mass_flux_data = getPatchData(patch, d_massflux, current_context);
  auto mass_flux_x = SAMRAI::pdat::get_view<Dim, SideData>(mass_flux_data, Side::X);
  auto mass_flux_y = SAMRAI::pdat::get_view<Dim, SideData>(mass_flux_data, Side::Y);

  auto pre_vol = SAMRAI::pdat::get_view<Dim, NodeData>(getPatchData(patch, d_workarray1, current_context));
  auto post_vol = SAMRAI::pdat::get_view<Dim, NodeData>(getPatchData(patch, d_workarray2, current_context));
  auto pre_mass = SAMRAI::pdat::get_view<Dim, NodeData>(getPatchData(patch, d_workarray3, current_context));
  auto post_mass = SAMRAI::pdat::get_view<Dim, NodeData>(getPatchData(patch, d_workarray4, current_context));
  auto advec_vol = SAMRAI::pdat::get_view<Dim, NodeData>(getPatchData(patch, d_workarray5, current_context));
  auto post_ener = SAMRAI::pdat::get_view<Dim, NodeData>(getPatchData(patch, d_workarray6, current_context));
  auto ener_flux = SAMRAI::pdat::get_view<Dim, NodeData>(getPatchData(patch, d_workarray7, current_context));

  std::shared_ptr<NodeData> delta_data = getPatchData(patch, d_vertexdeltas, current_context);
  auto vertexdx = SAMRAI::pdat::get_view<Dim, NodeData>(delta_data, Coord::X);
  auto vertexdy = SAMRAI::pdat::get_view<Dim, NodeData>(delta_data, Coord::Y);

  const Real one_by_six = 1.0/6.0;

  const int x_max = patch.getBox().upper(0);
  const int y_max = patch.getBox().upper(1);

  const SAMRAI::hier::Box& cell_ghost_box = getIndexBox<true>(volume_data);
  const SAMRAI::hier::Box& cell_box = getIndexBox<false>(volume_data);

  if (dir == Coord::X) {
    if (sweep_number == 1) {
  FORCEINLINE_MACRO
      SAMRAI::pdat::parallel_for_all(cell_ghost_box,  [=] SAMRAI_HOST_DEVICE (int k, int j) {
        pre_vol(j,k)=volume(j,k)+(vol_flux_x(j+1,k  )-vol_flux_x(j,k)+vol_flux_y(j  ,k+1)-vol_flux_y(j,k));
        post_vol(j,k)=pre_vol(j,k)-(vol_flux_x(j+1,k  )-vol_flux_x(j,k));
      });
    } else {
  FORCEINLINE_MACRO
      SAMRAI::pdat::parallel_for_all(cell_ghost_box,  [=] SAMRAI_HOST_DEVICE (int k, int j) {
        pre_vol(j,k)=volume(j,k)+vol_flux_x(j+1,k)-vol_flux_x(j,k);
        post_vol(j,k)=volume(j,k);
      });
    }

    SAMRAI::hier::Box advec_cell_box = patch.getBox();
    advec_cell_box.growUpper(Coord::X, d_nghosts[Coord::X]);

  FORCEINLINE_MACRO
    SAMRAI::pdat::parallel_for_all(advec_cell_box, [=] SAMRAI_HOST_DEVICE (int k, int j) {
      //RAJA::forallN<PatchPolicy>(isets->getCells(1), isets->getAdvecCells(0), [=] SAMRAI_HOST_DEVICE (int k, int j) {
      const int upwind = vol_flux_x(j,k) > 0.0 ? j-2 : MIN(j+1, x_max+2);
      const int donor = vol_flux_x(j,k) > 0.0 ? j-1 : j;
      const int downwind = vol_flux_x(j,k) > 0.0 ? j : j-1;
      const int diff = vol_flux_x(j,k) > 0.0 ? donor : upwind;

      const Real sigmat=fabs(vol_flux_x(j,k))/pre_vol(donor,k);
      const Real sigma3=(1.0+sigmat)*(vertexdx(j,k)/vertexdx(diff,k));
      const Real sigma4=2.0-sigmat;

      const Real sigma=sigmat;
      const Real sigmav=sigmat;

      Real diffuw=density1(donor,k)-density1(upwind,k);
      Real diffdw=density1(downwind,k)-density1(donor,k);
      Real wind = diffdw <= 0.0 ? -1.0 : 1.0;
      Real limiter = (diffuw*diffdw) > 0.0 ?
        (1.0-sigmav)*wind*MIN(fabs(diffuw),MIN(fabs(diffdw),one_by_six*(sigma3*fabs(diffuw)+sigma4*fabs(diffdw))))
        : 0.0;

      mass_flux_x(j,k)=vol_flux_x(j,k)*(density1(donor,k)+limiter);

      const Real sigmam=fabs(mass_flux_x(j,k))/(density1(donor,k)*pre_vol(donor,k));
      diffuw=energy1(donor,k)-energy1(upwind,k);
      diffdw=energy1(downwind,k)-energy1(donor,k);

      wind = diffdw <= 0.0 ? -1.0 : 1.0;

      limiter = diffuw*diffdw > 0.0 ? (1.0-sigmam)*wind*MIN(fabs(diffuw),MIN(fabs(diffdw),one_by_six*(sigma3*fabs(diffuw)+sigma4*fabs(diffdw)))) : 0.0;

      ener_flux(j,k)=mass_flux_x(j,k)*(energy1(donor,k)+limiter);
    });

  FORCEINLINE_MACRO
    SAMRAI::pdat::parallel_for_all(cell_box,  [=] SAMRAI_HOST_DEVICE (int k, int j) {
      pre_mass(j,k)=density1(j,k)*pre_vol(j,k);
      post_mass(j,k)=pre_mass(j,k)+mass_flux_x(j,k)-mass_flux_x(j+1,k);
      post_ener(j,k)=(energy1(j,k)*pre_mass(j,k)+ener_flux(j,k)-ener_flux(j+1,k))/post_mass(j,k);
      advec_vol(j,k)=pre_vol(j,k)+vol_flux_x(j,k)-vol_flux_x(j+1,k);
      density1(j,k)=post_mass(j,k)/advec_vol(j,k);
      energy1(j,k)=post_ener(j,k);
    });
  } else if (dir == Coord::Y) {
    if (sweep_number == 1) {
  FORCEINLINE_MACRO
            SAMRAI::pdat::parallel_for_all(cell_ghost_box,  [=] SAMRAI_HOST_DEVICE (int k, int j) {
        //RAJA::forallN<PatchPolicy>(isets->getCells(1,true), isets->getCells(0,true), [=] SAMRAI_HOST_DEVICE (int k, int j) {
        pre_vol(j,k)=volume(j,k)+(vol_flux_y(j  ,k+1)-vol_flux_y(j,k)+vol_flux_x(j+1,k  )-vol_flux_x(j,k));
        post_vol(j,k)=pre_vol(j,k)-(vol_flux_y(j  ,k+1)-vol_flux_y(j,k));
      });
    } else {
  FORCEINLINE_MACRO
            SAMRAI::pdat::parallel_for_all(cell_ghost_box,  [=] SAMRAI_HOST_DEVICE (int k, int j) {
        //RAJA::forallN<PatchPolicy>(isets->getCells(1,true), isets->getCells(0,true), [=] SAMRAI_HOST_DEVICE (int k, int j) {
        pre_vol(j,k)=volume(j,k)+vol_flux_y(j  ,k+1)-vol_flux_y(j,k);
        post_vol(j,k)=volume(j,k);
      });
    }

    SAMRAI::hier::Box advec_cell_box = patch.getBox();
    advec_cell_box.growUpper(Coord::Y, d_nghosts[Coord::Y]);

  FORCEINLINE_MACRO
    SAMRAI::pdat::parallel_for_all(advec_cell_box,  [=] SAMRAI_HOST_DEVICE (int k, int j) {
      //RAJA::forallN<PatchPolicy>(isets->getAdvecCells(1), isets->getCells(0), [=] SAMRAI_HOST_DEVICE (int k, int j) {
      const int upwind = vol_flux_y(j,k) > 0.0 ? k-2 : MIN(k+1, y_max+2);
      const int donor = vol_flux_y(j,k) > 0.0 ? k-1 : k;
      const int downwind = vol_flux_y(j,k) > 0.0 ? k : k-1;
      const int diff = vol_flux_y(j,k) > 0.0 ? donor : upwind;

      const Real sigmat=fabs(vol_flux_y(j,k))/pre_vol(j,donor);
      const Real sigma3=(1.0+sigmat)*(vertexdy(j,k)/vertexdy(j,diff));
      const Real sigma4=2.0-sigmat;

      const Real sigma=sigmat;
      const Real sigmav=sigmat;

      Real diffuw=density1(j,donor)-density1(j,upwind);
      Real diffdw=density1(j,downwind)-density1(j,donor);

      Real wind = diffdw <= 0.0 ? -1.0 : 1.0;

      Real limiter = diffuw*diffdw > 0.0 ? (1.0-sigmav)*wind*MIN(fabs(diffuw),MIN(fabs(diffdw),one_by_six*(sigma3*fabs(diffuw)+sigma4*fabs(diffdw)))) : 0.0;

      mass_flux_y(j,k)=vol_flux_y(j,k)*(density1(j,donor)+limiter);

      const Real sigmam=fabs(mass_flux_y(j,k))/(density1(j,donor)*pre_vol(j,donor));
      diffuw=energy1(j,donor)-energy1(j,upwind);
      diffdw=energy1(j,downwind)-energy1(j,donor);

      wind = diffdw <= 0.0 ? -1.0 : 1.0;
      limiter = diffuw*diffdw > 0.0 ? (1.0-sigmam)*wind*MIN(fabs(diffuw),MIN(fabs(diffdw),one_by_six*(sigma3*fabs(diffuw)+sigma4*fabs(diffdw)))) : 0.0;

      ener_flux(j,k)=mass_flux_y(j,k)*(energy1(j,donor)+limiter);
    });

  FORCEINLINE_MACRO
    SAMRAI::pdat::parallel_for_all(cell_box,  [=] SAMRAI_HOST_DEVICE (int k, int j) {
      pre_mass(j,k)=density1(j,k)*pre_vol(j,k);
      post_mass(j,k)=pre_mass(j,k)+mass_flux_y(j,k)-mass_flux_y(j,k+1);
      post_ener(j,k)=(energy1(j,k)*pre_mass(j,k)+ener_flux(j,k)-ener_flux(j,k+1))/post_mass(j,k);
      advec_vol(j,k)=pre_vol(j,k)+vol_flux_y(j,k)-vol_flux_y(j,k+1);
      density1(j,k)=post_mass(j,k)/advec_vol(j,k);
      energy1(j,k)=post_ener(j,k);
    });
  }
}

void Cleverleaf::advec_mom(
  SAMRAI::hier::Patch& patch,
  const int sweep_number,
  const int direction,
  const int which_vel)
{
  auto current_context = getCurrentDataContext();
  auto new_context = getNewDataContext();

  std::shared_ptr<CellData> volume_data = getPatchData(patch, d_volume, current_context);
  auto volume = SAMRAI::pdat::get_view<Dim, CellData>(volume_data);

  auto density1 = SAMRAI::pdat::get_view<Dim, CellData>(getPatchData(patch, d_density, new_context));

  std::shared_ptr<SideData> vol_flux_data = getPatchData(patch, d_volflux, current_context);
  auto vol_flux_x = SAMRAI::pdat::get_view<Dim, SideData>(vol_flux_data, Side::X);
  auto vol_flux_y = SAMRAI::pdat::get_view<Dim, SideData>(vol_flux_data, Side::Y);

  std::shared_ptr<SideData> mass_flux_data = getPatchData(patch, d_massflux, current_context);
  auto mass_flux_x = SAMRAI::pdat::get_view<Dim, SideData>(mass_flux_data, Side::X);
  auto mass_flux_y = SAMRAI::pdat::get_view<Dim, SideData>(mass_flux_data, Side::Y);

  std::shared_ptr<NodeData> node_flux_data = getPatchData(patch, d_workarray1, current_context);
  auto node_flux = SAMRAI::pdat::get_view<Dim, NodeData>(node_flux_data);
  auto node_mass_post = SAMRAI::pdat::get_view<Dim, NodeData>(getPatchData(patch, d_workarray2, current_context));
  auto node_mass_pre = SAMRAI::pdat::get_view<Dim, NodeData>(getPatchData(patch, d_workarray3, current_context));

  auto advec_vel = SAMRAI::pdat::get_view<Dim, NodeData>(getPatchData(patch, d_workarray4, current_context));
  auto mom_flux = SAMRAI::pdat::get_view<Dim, NodeData>(getPatchData(patch, d_workarray5, current_context));

  auto pre_vol = SAMRAI::pdat::get_view<Dim, NodeData>(getPatchData(patch, d_workarray6, current_context));
  auto post_vol = SAMRAI::pdat::get_view<Dim, NodeData>(getPatchData(patch, d_workarray7, current_context));

  std::shared_ptr<CellData> delta_data = getPatchData(patch, d_celldeltas, current_context);
  auto celldx = SAMRAI::pdat::get_view<Dim, CellData>(delta_data, Coord::X);
  auto celldy = SAMRAI::pdat::get_view<Dim, CellData>(delta_data, Coord::Y);

  std::shared_ptr<NodeData> new_velocity_data = getPatchData(patch, d_velocity, new_context);
  auto xvel1 = SAMRAI::pdat::get_view<Dim, NodeData>(new_velocity_data, Coord::X);
  auto yvel1 = SAMRAI::pdat::get_view<Dim, NodeData>(new_velocity_data, Coord::Y);

  const int mom_sweep = (direction+1) + 2*(sweep_number-1);

  const SAMRAI::hier::Box& cell_ghost_box = getIndexBox<true>(volume_data);
  const SAMRAI::hier::Box& node_ghost_box = getIndexBox<true>(node_flux_data);

  const SAMRAI::hier::Box& node_box = getIndexBox<false>(node_flux_data);

  if (which_vel == Coord::X) {
          auto vel1 = SAMRAI::pdat::get_view<Dim, NodeData>(getPatchData(patch, d_velocity, new_context), Coord::X);

    if (mom_sweep == 1) {
  FORCEINLINE_MACRO
            SAMRAI::pdat::parallel_for_all(cell_ghost_box,  [=] SAMRAI_HOST_DEVICE (int k, int j) {
        post_vol(j,k)= volume(j,k)+vol_flux_y(j  ,k+1)-vol_flux_y(j,k);
        pre_vol(j,k)=post_vol(j,k)+vol_flux_x(j+1,k  )-vol_flux_x(j,k);
      });
    } else if (mom_sweep == 2) {
  FORCEINLINE_MACRO
            SAMRAI::pdat::parallel_for_all(cell_ghost_box,  [=] SAMRAI_HOST_DEVICE (int k, int j) {
        post_vol(j,k)= volume(j,k)+vol_flux_x(j+1,k  )-vol_flux_x(j,k);
        pre_vol(j,k)=post_vol(j,k)+vol_flux_y(j  ,k+1)-vol_flux_y(j,k);
      });
    } else if (mom_sweep == 3) {
  FORCEINLINE_MACRO
            SAMRAI::pdat::parallel_for_all(cell_ghost_box,  [=] SAMRAI_HOST_DEVICE (int k, int j) {
        post_vol(j,k)=volume(j,k);
        pre_vol(j,k)=post_vol(j,k)+vol_flux_y(j  ,k+1)-vol_flux_y(j,k);
      });
    } else if (mom_sweep == 4) {
  FORCEINLINE_MACRO
            SAMRAI::pdat::parallel_for_all(cell_ghost_box,  [=] SAMRAI_HOST_DEVICE (int k, int j) {
        post_vol(j,k)=volume(j,k);
        pre_vol(j,k)=post_vol(j,k)+vol_flux_x(j+1,k  )-vol_flux_x(j,k);
      });
    }

    if (direction == Coord::X) {
      if (which_vel == Coord::X) {
        SAMRAI::hier::Box cell_node_box = patch.getBox();
        cell_node_box.grow(Coord::X, d_nghosts[Coord::X]);
        cell_node_box.growUpper(Coord::Y, 1);
  FORCEINLINE_MACRO
        SAMRAI::pdat::parallel_for_all(cell_node_box, [=] SAMRAI_HOST_DEVICE (int k, int j) {
          // RAJA::forallN<PatchPolicy>(isets->getNodes(1), isets->getCells(0,true), [=] SAMRAI_HOST_DEVICE (int k, int j) {
          node_flux(j,k)=0.25*(mass_flux_x(j,k-1  )+mass_flux_x(j  ,k)
                               +mass_flux_x(j+1,k-1)+mass_flux_x(j+1,k));
        });

        SAMRAI::hier::Box staggered_cell_node_box = patch.getBox();
        staggered_cell_node_box.growLower(Coord::X, 1);
        staggered_cell_node_box.growUpper(Coord::X, d_nghosts[Coord::X]);
        staggered_cell_node_box.growUpper(Coord::Y, 1);

  FORCEINLINE_MACRO
        SAMRAI::pdat::parallel_for_all(staggered_cell_node_box, [=] SAMRAI_HOST_DEVICE (int k, int j) {
          // RAJA::forallN<PatchPolicy>(isets->getNodes(1), isets->getStaggeredCells(0), [=] SAMRAI_HOST_DEVICE (int k, int j) {
          node_mass_post(j,k)=0.25*(density1(j  ,k-1)*post_vol(j  ,k-1)
                                    +density1(j  ,k  )*post_vol(j  ,k  )
                                    +density1(j-1,k-1)*post_vol(j-1,k-1)
                                    +density1(j-1,k  )*post_vol(j-1,k  ));
        });

  FORCEINLINE_MACRO
        SAMRAI::pdat::parallel_for_all(staggered_cell_node_box, [=] SAMRAI_HOST_DEVICE (int k, int j) {
          // RAJA::forallN<PatchPolicy>(isets->getNodes(1), isets->getStaggeredCells(0), [=] SAMRAI_HOST_DEVICE (int k, int j) {
          node_mass_pre(j,k)=node_mass_post(j,k)-node_flux(j-1,k)+node_flux(j,k);
        });
      }

      SAMRAI::hier::Box oneghost_node_box = patch.getBox();
      oneghost_node_box.grow(Coord::X, 1);
      oneghost_node_box.growUpper(Coord::Y, 1);
  FORCEINLINE_MACRO
      SAMRAI::pdat::parallel_for_all(oneghost_node_box, [=] SAMRAI_HOST_DEVICE (int k, int j) {
        // RAJA::forallN<PatchPolicy>(isets->getNodes(1), isets->getOneGhostCells(0), [=] SAMRAI_HOST_DEVICE (int k, int j) {
        const int upwind = (node_flux(j,k) < 0.0) ? j+2 : j-1;
        const int donor = (node_flux(j,k) < 0.0) ? j+1 : j;
        const int downwind = (node_flux(j,k) < 0.0) ? j : j+1;
        const int dif = (node_flux(j,k) < 0.0) ? donor : upwind;

        const Real sigma=fabs(node_flux(j,k))/(node_mass_pre(donor,k));
        const Real width=celldx(j,k);
        const Real vdiffuw=vel1(donor,k)-vel1(upwind,k);
        const Real vdiffdw=vel1(downwind,k)-vel1(donor,k);

        Real limiter=0.0;

        if (vdiffuw*vdiffdw > 0.0) {
          const Real auw=fabs(vdiffuw);
          const Real adw=fabs(vdiffdw);
          const Real wind = (vdiffdw <= 0.0) ? -1.0 : 1.0;
          limiter=wind*MIN(width*((2.0-sigma)*adw/width+(1.0+sigma)*auw/celldx(dif,k))/6.0,MIN(auw,adw));
        }
        advec_vel(j,k)=vel1(donor,k)+(1.0-sigma)*limiter;
        mom_flux(j,k)=advec_vel(j,k)*node_flux(j,k);
      });

  FORCEINLINE_MACRO
      SAMRAI::pdat::parallel_for_all(node_box, [=] SAMRAI_HOST_DEVICE (int k, int j) {
        // RAJA::forallN<PatchPolicy>(isets->getNodes(1), isets->getNodes(0), [=] SAMRAI_HOST_DEVICE (int k, int j) {
        vel1 (j,k)=(vel1 (j,k)*node_mass_pre(j,k)+mom_flux(j-1,k)-mom_flux(j,k))/node_mass_post(j,k);
      });

    } else if (direction == Coord::Y) {
      if (which_vel == Coord::X) {
        SAMRAI::hier::Box node_cell_ghost_box = patch.getBox();
        node_cell_ghost_box.growUpper(Coord::X, 1);
        node_cell_ghost_box.grow(Coord::Y, d_nghosts[Coord::Y]);
  FORCEINLINE_MACRO
        SAMRAI::pdat::parallel_for_all(node_cell_ghost_box, [=] SAMRAI_HOST_DEVICE (int k, int j) {
          // RAJA::forallN<PatchPolicy>(isets->getCells(1,true), isets->getNodes(0), [=] SAMRAI_HOST_DEVICE (int k, int j) {
          node_flux(j,k)=0.25*(mass_flux_y(j-1,k  )+mass_flux_y(j  ,k  )
                               +mass_flux_y(j-1,k+1)+mass_flux_y(j  ,k+1));
        });

        SAMRAI::hier::Box node_staggered_cell_box = patch.getBox();
        node_staggered_cell_box.growUpper(Coord::X, 1);
        node_staggered_cell_box.growLower(Coord::Y, 1);
        node_staggered_cell_box.growUpper(Coord::Y, d_nghosts[Coord::Y]);

  FORCEINLINE_MACRO
        SAMRAI::pdat::parallel_for_all(node_staggered_cell_box, [=] SAMRAI_HOST_DEVICE (int k, int j) {
          // RAJA::forallN<PatchPolicy>(isets->getStaggeredCells(1), isets->getNodes(0), [=] SAMRAI_HOST_DEVICE (int k, int j) {
          node_mass_post(j,k)=0.25*(density1(j  ,k-1)*post_vol(j  ,k-1)
                                    +density1(j  ,k  )*post_vol(j  ,k  )
                                    +density1(j-1,k-1)*post_vol(j-1,k-1)
                                    +density1(j-1,k  )*post_vol(j-1,k  ));
        });

  FORCEINLINE_MACRO
        SAMRAI::pdat::parallel_for_all(node_staggered_cell_box, [=] SAMRAI_HOST_DEVICE (int k, int j) {
          // RAJA::forallN<PatchPolicy>(isets->getStaggeredCells(1), isets->getNodes(0), [=] SAMRAI_HOST_DEVICE (int k, int j) {
          node_mass_pre(j,k)=node_mass_post(j,k)-node_flux(j,k-1)+node_flux(j,k);
        });
      }

      SAMRAI::hier::Box node_oneghost_box = patch.getBox();
      node_oneghost_box.growUpper(Coord::X, 1);
      node_oneghost_box.grow(Coord::Y, 1);
  FORCEINLINE_MACRO
      SAMRAI::pdat::parallel_for_all(node_oneghost_box, [=] SAMRAI_HOST_DEVICE (int k, int j) {
        // RAJA::forallN<PatchPolicy>(isets->getOneGhostCells(1), isets->getNodes(0), [=] SAMRAI_HOST_DEVICE (int k, int j) {
        const int upwind = (node_flux(j,k) < 0.0) ?  k+2 : k-1;
        const int donor = (node_flux(j,k) < 0.0) ?  k+1 : k;
        const int downwind = (node_flux(j,k) < 0.0) ?  k : k+1;
        const int dif = (node_flux(j,k) < 0.0) ?  donor : upwind;

        const Real sigma=fabs(node_flux(j,k))/(node_mass_pre(j,donor));
        const Real width=celldy(j,k);
        const Real vdiffuw=vel1(j,donor)-vel1(j,upwind);
        const Real vdiffdw=vel1(j,downwind)-vel1(j,donor);

        Real limiter=0.0;

        if (vdiffuw*vdiffdw > 0.0) {
          const Real auw=fabs(vdiffuw);
          const Real adw=fabs(vdiffdw);
          const Real wind = (vdiffdw <= 0.0) ? -1.0 : 1.0;
          limiter=wind*MIN(width*((2.0-sigma)*adw/width+(1.0+sigma)*auw/celldy(j,dif))/6.0,MIN(auw,adw));
        }

        advec_vel(j,k)=vel1(j,donor)+(1.0-sigma)*limiter;
        mom_flux(j,k)=advec_vel(j,k)*node_flux(j,k);
      });

  FORCEINLINE_MACRO
      SAMRAI::pdat::parallel_for_all(node_box, [=] SAMRAI_HOST_DEVICE (int k, int j) {
        // RAJA::forallN<PatchPolicy>(isets->getNodes(1), isets->getNodes(0), [=] SAMRAI_HOST_DEVICE (int k, int j) {
        vel1 (j,k)=(vel1(j,k)*node_mass_pre(j,k)+mom_flux(j,k-1)-mom_flux(j,k))/node_mass_post(j,k);
      });
    }
  } else {
          auto vel1 = SAMRAI::pdat::get_view<Dim, NodeData>(getPatchData(patch, d_velocity, new_context), Coord::Y);

    if (mom_sweep == 1) {
  FORCEINLINE_MACRO
            SAMRAI::pdat::parallel_for_all(cell_ghost_box, [=] SAMRAI_HOST_DEVICE (int k, int j) {
        // RAJA::forallN<PatchPolicy>(isets->getCells(1,true), isets->getCells(0,true), [=] SAMRAI_HOST_DEVICE (int k, int j) {
        post_vol(j,k)= volume(j,k)+vol_flux_y(j  ,k+1)-vol_flux_y(j,k);
        pre_vol(j,k)=post_vol(j,k)+vol_flux_x(j+1,k  )-vol_flux_x(j,k);
      });
    } else if (mom_sweep == 2) {
  FORCEINLINE_MACRO
            SAMRAI::pdat::parallel_for_all(cell_ghost_box, [=] SAMRAI_HOST_DEVICE (int k, int j) {
        // RAJA::forallN<PatchPolicy>(isets->getCells(1,true), isets->getCells(0,true), [=] SAMRAI_HOST_DEVICE (int k, int j) {
        post_vol(j,k)= volume(j,k)+vol_flux_x(j+1,k  )-vol_flux_x(j,k);
        pre_vol(j,k)=post_vol(j,k)+vol_flux_y(j  ,k+1)-vol_flux_y(j,k);
      });
    } else if (mom_sweep == 3) {
  FORCEINLINE_MACRO
            SAMRAI::pdat::parallel_for_all(cell_ghost_box, [=] SAMRAI_HOST_DEVICE (int k, int j) {
        // RAJA::forallN<PatchPolicy>(isets->getCells(1,true), isets->getCells(0,true), [=] SAMRAI_HOST_DEVICE (int k, int j) {
        post_vol(j,k)=volume(j,k);
        pre_vol(j,k)=post_vol(j,k)+vol_flux_y(j  ,k+1)-vol_flux_y(j,k);
      });
    } else if (mom_sweep == 4) {
  FORCEINLINE_MACRO
            SAMRAI::pdat::parallel_for_all(cell_ghost_box, [=] SAMRAI_HOST_DEVICE (int k, int j) {
        // RAJA::forallN<PatchPolicy>(isets->getCells(1,true), isets->getCells(0,true), [=] SAMRAI_HOST_DEVICE (int k, int j) {
        post_vol(j,k)=volume(j,k);
        pre_vol(j,k)=post_vol(j,k)+vol_flux_x(j+1,k  )-vol_flux_x(j,k);
      });
    }

    if (direction == Coord::X) {
      if (which_vel == Coord::X) {
        SAMRAI::hier::Box cell_ghost_node_box = patch.getBox();
        cell_ghost_node_box.grow(Coord::X, d_nghosts[Coord::X]);
        cell_ghost_node_box.growUpper(Coord::Y, 1);
  FORCEINLINE_MACRO
        SAMRAI::pdat::parallel_for_all(cell_ghost_node_box, [=] SAMRAI_HOST_DEVICE (int k, int j) {
          // RAJA::forallN<PatchPolicy>(isets->getNodes(1), isets->getCells(0,true), [=] SAMRAI_HOST_DEVICE (int k, int j) {
          node_flux(j,k)=0.25*(mass_flux_x(j,k-1  )+mass_flux_x(j  ,k)
                               +mass_flux_x(j+1,k-1)+mass_flux_x(j+1,k));
        });

        SAMRAI::hier::Box staggered_cell_node_box = patch.getBox();
        staggered_cell_node_box.growLower(Coord::X, 1);
        staggered_cell_node_box.growUpper(Coord::X, d_nghosts[Coord::X]);
        staggered_cell_node_box.growUpper(Coord::Y, 1);

  FORCEINLINE_MACRO
        SAMRAI::pdat::parallel_for_all(staggered_cell_node_box, [=] SAMRAI_HOST_DEVICE (int k, int j) {
          // RAJA::forallN<PatchPolicy>(isets->getNodes(1), isets->getStaggeredCells(0), [=] SAMRAI_HOST_DEVICE (int k, int j) {
          node_mass_post(j,k)=0.25*(density1(j  ,k-1)*post_vol(j  ,k-1)
                                    +density1(j  ,k  )*post_vol(j  ,k  )
                                    +density1(j-1,k-1)*post_vol(j-1,k-1)
                                    +density1(j-1,k  )*post_vol(j-1,k  ));
        });

  FORCEINLINE_MACRO
        SAMRAI::pdat::parallel_for_all(staggered_cell_node_box, [=] SAMRAI_HOST_DEVICE (int k, int j) {
          // RAJA::forallN<PatchPolicy>(isets->getNodes(1), isets->getStaggeredCells(0), [=] SAMRAI_HOST_DEVICE (int k, int j) {
          node_mass_pre(j,k)=node_mass_post(j,k)-node_flux(j-1,k)+node_flux(j,k);
        });
      }

      SAMRAI::hier::Box oneghost_node_box = patch.getBox();
      oneghost_node_box.grow(Coord::X, 1);
      oneghost_node_box.growUpper(Coord::Y, 1);

  FORCEINLINE_MACRO
      SAMRAI::pdat::parallel_for_all(oneghost_node_box, [=] SAMRAI_HOST_DEVICE (int k, int j) {
        // RAJA::forallN<PatchPolicy>(isets->getNodes(1), isets->getOneGhostCells(0), [=] SAMRAI_HOST_DEVICE (int k, int j) {
        const int upwind = (node_flux(j,k) < 0.0) ? j+2 : j-1;
        const int donor = (node_flux(j,k) < 0.0) ? j+1 : j;
        const int downwind = (node_flux(j,k) < 0.0) ? j : j+1;
        const int dif = (node_flux(j,k) < 0.0) ? donor : upwind;

        const Real sigma=fabs(node_flux(j,k))/(node_mass_pre(donor,k));
        const Real width=celldx(j,k);
        const Real vdiffuw=vel1(donor,k)-vel1(upwind,k);
        const Real vdiffdw=vel1(downwind,k)-vel1(donor,k);

        Real limiter=0.0;

        if (vdiffuw*vdiffdw > 0.0) {
          const Real auw=fabs(vdiffuw);
          const Real adw=fabs(vdiffdw);
          const Real wind = (vdiffdw <= 0.0) ? -1.0 : 1.0;
          limiter=wind*MIN(width*((2.0-sigma)*adw/width+(1.0+sigma)*auw/celldx(dif,k))/6.0,MIN(auw,adw));
        }
        advec_vel(j,k)=vel1(donor,k)+(1.0-sigma)*limiter;
        mom_flux(j,k)=advec_vel(j,k)*node_flux(j,k);
      });

  FORCEINLINE_MACRO
      SAMRAI::pdat::parallel_for_all(node_box, [=] SAMRAI_HOST_DEVICE (int k, int j) {
        // RAJA::forallN<PatchPolicy>(isets->getNodes(1), isets->getNodes(0), [=] SAMRAI_HOST_DEVICE (int k, int j) {
        vel1 (j,k)=(vel1 (j,k)*node_mass_pre(j,k)+mom_flux(j-1,k)-mom_flux(j,k))/node_mass_post(j,k);
      });

    } else if (direction == Coord::Y) {
      if (which_vel == Coord::X) {
        SAMRAI::hier::Box node_cell_ghost_box = patch.getBox();
        node_cell_ghost_box.growUpper(Coord::X, 1);
        node_cell_ghost_box.grow(Coord::Y, d_nghosts[Coord::Y]);
  FORCEINLINE_MACRO
        SAMRAI::pdat::parallel_for_all(node_cell_ghost_box, [=] SAMRAI_HOST_DEVICE (int k, int j) {
          // RAJA::forallN<PatchPolicy>(isets->getCells(1,true), isets->getNodes(0), [=] SAMRAI_HOST_DEVICE (int k, int j) {
          node_flux(j,k)=0.25*(mass_flux_y(j-1,k  )+mass_flux_y(j  ,k  )
                               +mass_flux_y(j-1,k+1)+mass_flux_y(j  ,k+1));
        });

        SAMRAI::hier::Box node_staggered_cell_box = patch.getBox();
        node_staggered_cell_box.growUpper(Coord::X, 1);
        node_staggered_cell_box.growLower(Coord::Y, 1);
        node_staggered_cell_box.growUpper(Coord::Y, d_nghosts[Coord::Y]);

  FORCEINLINE_MACRO
        SAMRAI::pdat::parallel_for_all(node_staggered_cell_box, [=] SAMRAI_HOST_DEVICE (int k, int j) {
          // RAJA::forallN<PatchPolicy>(isets->getStaggeredCells(1), isets->getNodes(0), [=] SAMRAI_HOST_DEVICE (int k, int j) {
          node_mass_post(j,k)=0.25*(density1(j  ,k-1)*post_vol(j  ,k-1)
                                    +density1(j  ,k  )*post_vol(j  ,k  )
                                    +density1(j-1,k-1)*post_vol(j-1,k-1)
                                    +density1(j-1,k  )*post_vol(j-1,k  ));
        });

  FORCEINLINE_MACRO
        SAMRAI::pdat::parallel_for_all(node_staggered_cell_box, [=] SAMRAI_HOST_DEVICE (int k, int j) {
          // RAJA::forallN<PatchPolicy>(isets->getStaggeredCells(1), isets->getNodes(0), [=] SAMRAI_HOST_DEVICE (int k, int j) {
          node_mass_pre(j,k)=node_mass_post(j,k)-node_flux(j,k-1)+node_flux(j,k);
        });
      }

      SAMRAI::hier::Box node_oneghost_box = patch.getBox();
      node_oneghost_box.growUpper(Coord::X, 1);
      node_oneghost_box.grow(Coord::Y, 1);

  FORCEINLINE_MACRO
      SAMRAI::pdat::parallel_for_all(node_oneghost_box, [=] SAMRAI_HOST_DEVICE (int k, int j) {
        // RAJA::forallN<PatchPolicy>(isets->getOneGhostCells(1), isets->getNodes(0), [=] SAMRAI_HOST_DEVICE (int k, int j) {
        const int upwind = (node_flux(j,k) < 0.0) ?  k+2 : k-1;
        const int donor = (node_flux(j,k) < 0.0) ?  k+1 : k;
        const int downwind = (node_flux(j,k) < 0.0) ?  k : k+1;
        const int dif = (node_flux(j,k) < 0.0) ?  donor : upwind;

        const Real sigma=fabs(node_flux(j,k))/(node_mass_pre(j,donor));
        const Real width=celldy(j,k);
        const Real vdiffuw=vel1(j,donor)-vel1(j,upwind);
        const Real vdiffdw=vel1(j,downwind)-vel1(j,donor);

        Real limiter=0.0;

        if (vdiffuw*vdiffdw > 0.0) {
          const Real auw=fabs(vdiffuw);
          const Real adw=fabs(vdiffdw);
          const Real wind = (vdiffdw <= 0.0) ? -1.0 : 1.0;
          limiter=wind*MIN(width*((2.0-sigma)*adw/width+(1.0+sigma)*auw/celldy(j,dif))/6.0,MIN(auw,adw));
        }

        advec_vel(j,k)=vel1(j,donor)+(1.0-sigma)*limiter;
        mom_flux(j,k)=advec_vel(j,k)*node_flux(j,k);
      });

  FORCEINLINE_MACRO
      SAMRAI::pdat::parallel_for_all(node_box, [=] SAMRAI_HOST_DEVICE (int k, int j) {
        // RAJA::forallN<PatchPolicy>(isets->getNodes(1), isets->getNodes(0), [=] SAMRAI_HOST_DEVICE (int k, int j) {
        vel1 (j,k)=(vel1(j,k)*node_mass_pre(j,k)+mom_flux(j,k-1)-mom_flux(j,k))/node_mass_post(j,k);
      });
    }
  }
}

void Cleverleaf::setPhysicalBoundaryConditions(
  SAMRAI::hier::Patch& patch,
  const double fill_time,
  const SAMRAI::hier::IntVector& ghost_width_to_fill)
{
  auto current_context = getCurrentDataContext();
  auto new_context = getNewDataContext();

  t_fill_boundary->start();

  // const int depth = ghost_width_to_fill[0];

  const std::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> pgeom(
    SAMRAI_SHARED_PTR_CAST<SAMRAI::geom::CartesianPatchGeometry>(patch.getPatchGeometry()));

  const std::vector<hier::BoundaryBox>& edge_bdry(
    pgeom->getCodimensionBoundaries(Bdry::EDGE2D));

  for (const hier::BoundaryBox& ebb : edge_bdry) { // small loop (4 cases)
    const int edge = ebb.getLocationIndex();
    const int edge_normal_dim = edge / Dim; // See SAMRAI/hier/BoundaryBox.h getLocationIndex() comment

    SAMRAI::hier::Box box = patch.getBox();

    SAMRAI::hier::IntVector gcw_to_fill = hier::IntVector::min(d_nghosts, ghost_width_to_fill);
    box.grow(gcw_to_fill);
    const int ghosts = gcw_to_fill[edge_normal_dim];

    SAMRAI::hier::Box cell_bbox  = SAMRAI::pdat::CellGeometry::toCellBox(box);
    SAMRAI::hier::Box node_bbox  = SAMRAI::pdat::NodeGeometry::toNodeBox(box);
    SAMRAI::hier::Box side_bbox0 = SAMRAI::pdat::SideGeometry::toSideBox(box, Side::X);
    SAMRAI::hier::Box side_bbox1 = SAMRAI::pdat::SideGeometry::toSideBox(box, Side::Y);

    switch(d_which_exchange) {
    case LagrangianEulerianLevelIntegrator::FIELD_EXCH:
      setBCs(getPatchData(patch, d_density, current_context), cell_bbox, edge, ghosts);
      setBCs(getPatchData(patch, d_energy, current_context), cell_bbox, edge, ghosts);
      setBCs(getPatchData(patch, d_velocity, current_context), node_bbox, edge, ghosts);
      break;
    case LagrangianEulerianLevelIntegrator::PRIME_CELLS_EXCH:
      setBCs(getPatchData(patch, d_density, current_context), cell_bbox, edge, ghosts);
      setBCs(getPatchData(patch, d_density, new_context), cell_bbox, edge, ghosts);

      setBCs(getPatchData(patch, d_energy, current_context), cell_bbox, edge, ghosts);
      setBCs(getPatchData(patch, d_energy, new_context), cell_bbox, edge, ghosts);

      setBCs(getPatchData(patch, d_pressure, current_context), cell_bbox, edge, ghosts);
      setBCs(getPatchData(patch, d_viscosity, current_context), cell_bbox, edge, ghosts);

      setBCs(getPatchData(patch, d_velocity, current_context), node_bbox, edge, ghosts);
      setBCs(getPatchData(patch, d_velocity, new_context), node_bbox, edge, ghosts);
      break;
    case LagrangianEulerianLevelIntegrator::PRE_LAGRANGE_EXCH:
      setBCs(getPatchData(patch, d_density, current_context), cell_bbox, edge, ghosts);
      setBCs(getPatchData(patch, d_energy, current_context), cell_bbox, edge, ghosts);
      setBCs(getPatchData(patch, d_pressure, current_context), cell_bbox, edge, ghosts);
      setBCs(getPatchData(patch, d_velocity, current_context), node_bbox, edge, ghosts);
      break;
    case LagrangianEulerianLevelIntegrator::POST_VISCOSITY_EXCH:
      setBCs(getPatchData(patch, d_viscosity, current_context), cell_bbox, edge, ghosts);
      break;
    case LagrangianEulerianLevelIntegrator::HALF_STEP_EXCH:
      setBCs(getPatchData(patch, d_pressure, current_context), cell_bbox, edge, ghosts);
      break;
    case LagrangianEulerianLevelIntegrator::PRE_SWEEP_1_CELL_EXCH:
      setBCs(getPatchData(patch, d_density, new_context), cell_bbox, edge, ghosts);
      setBCs(getPatchData(patch, d_energy, new_context), cell_bbox, edge, ghosts);
      setBCs(getPatchData(patch, d_volflux, current_context), side_bbox0, side_bbox1, edge, ghosts);
      break;
    case LagrangianEulerianLevelIntegrator::PRE_SWEEP_1_MOM_EXCH:
      setBCs(getPatchData(patch, d_density, new_context), cell_bbox, edge, ghosts);
      setBCs(getPatchData(patch, d_energy, new_context), cell_bbox, edge, ghosts);
      setBCs(getPatchData(patch, d_massflux, current_context), side_bbox0, side_bbox1, edge, ghosts);
      setBCs(getPatchData(patch, d_velocity, new_context), node_bbox, edge, ghosts);
      break;
    case LagrangianEulerianLevelIntegrator::PRE_SWEEP_2_MOM_EXCH:
      setBCs(getPatchData(patch, d_density, new_context), cell_bbox, edge, ghosts);
      setBCs(getPatchData(patch, d_energy, new_context), cell_bbox, edge, ghosts);
      setBCs(getPatchData(patch, d_massflux, current_context), side_bbox0, side_bbox1, edge, ghosts);
      setBCs(getPatchData(patch, d_velocity, new_context), node_bbox, edge, ghosts);
      break;
    default:
      SAMRAI::tbox::perr << "[ERROR] Unknown exchange id in setPhysicalBoundaryConditions... "
                         << std::endl;
      exit(-1);
    }
  }

  t_fill_boundary->stop();
}

void Cleverleaf::field_summary(
  hier::Patch& patch,
  double* total_volume,
  double* total_mass,
  double* total_pressure,
  double* total_internal_energy,
  double* total_kinetic_energy,
  int* total_effective_cells)
{
  auto context = getCurrentDataContext();

  auto volume = SAMRAI::pdat::get_view<Dim, CellData>(getPatchData(patch, d_volume, context));
  auto density0 = SAMRAI::pdat::get_view<Dim, CellData>(getPatchData(patch, d_density, context));
  auto energy0 = SAMRAI::pdat::get_view<Dim, CellData>(getPatchData(patch, d_energy, context));
  auto pressure = SAMRAI::pdat::get_view<Dim, CellData>(getPatchData(patch, d_pressure, context));
  auto xvel0 = SAMRAI::pdat::get_view<Dim, NodeData>(getPatchData(patch, d_velocity, context), Coord::X);
  auto yvel0 = SAMRAI::pdat::get_view<Dim, NodeData>(getPatchData(patch, d_velocity, context), Coord::Y);
  auto level_indicator = SAMRAI::pdat::get_view<Dim, IndicatorData>(getPatchData(patch, d_level_indicator, context));

  const int level_number = patch.getPatchLevelNumber();

  SAMRAI::tbox::parallel_reduction_variable_t<tbox::Reduction::Sum, Real> vol(0.0);
  SAMRAI::tbox::parallel_reduction_variable_t<tbox::Reduction::Sum, Real> mass(0.0);
  SAMRAI::tbox::parallel_reduction_variable_t<tbox::Reduction::Sum, Real> ie(0.0);
  SAMRAI::tbox::parallel_reduction_variable_t<tbox::Reduction::Sum, Real> ke(0.0);
  SAMRAI::tbox::parallel_reduction_variable_t<tbox::Reduction::Sum, Real> press(0.0);
  SAMRAI::tbox::parallel_reduction_variable_t<tbox::Reduction::Sum, int> cells(0);

  FORCEINLINE_MACRO
  SAMRAI::pdat::parallel_for_all(patch.getBox(), [=] SAMRAI_HOST_DEVICE (int k, int j) {
    // RAJA::forallN<PatchPolicy>(isets->getCells(1), isets->getCells(0), [=] SAMRAI_HOST_DEVICE (int k, int j) {
    if (level_indicator(j,k) == level_number) {
      double vsqrd = 0.25*(xvel0(j,k)*xvel0(j,k)+yvel0(j,k)*yvel0(j,k));
      vsqrd += 0.25*(xvel0(j+1,k)*xvel0(j+1,k)+yvel0(j+1,k)*yvel0(j+1,k));
      vsqrd += 0.25*(xvel0(j,k+1)*xvel0(j,k+1)+yvel0(j,k+1)*yvel0(j,k+1));
      vsqrd += 0.25*(xvel0(j+1,k+1)*xvel0(j+1,k+1)+yvel0(j+1,k+1)*yvel0(j+1,k+1));

      const double cell_vol=volume(j,k);
      const double cell_mass=cell_vol*density0(j,k);

      vol += cell_vol;
      mass += cell_mass;
      ie += cell_mass*energy0(j,k);
      ke += cell_mass*0.5*vsqrd;
      press += cell_vol*pressure(j,k);
      cells += 1;
    }
  });

  *total_volume = double(vol);
  *total_mass = double(mass);
  *total_pressure = double(press);
  *total_internal_energy = double(ie);
  *total_kinetic_energy = double(ke);
  *total_effective_cells = int(cells);
}

void Cleverleaf::tagGradientDetectorCells(
  SAMRAI::hier::Patch& patch,
  const double regrid_time,
  const bool initial_error,
  const int tag_index)
{
  auto context = getCurrentDataContext();

  std::shared_ptr<IndicatorData> samrai_tags_data(
    SAMRAI_SHARED_PTR_CAST<IndicatorData>(patch.getPatchData(tag_index)));
  samrai_tags_data->fillAll(0);

  auto tags = SAMRAI::pdat::get_view<Dim, IndicatorData>(getPatchData(patch, d_tags, context));
  auto density = SAMRAI::pdat::get_view<Dim, CellData>(getPatchData(patch, d_density, context));
  auto energy = SAMRAI::pdat::get_view<Dim, CellData>(getPatchData(patch, d_energy, context));
  auto pressure = SAMRAI::pdat::get_view<Dim, CellData>(getPatchData(patch, d_pressure, context));
  auto viscosity = SAMRAI::pdat::get_view<Dim, CellData>(getPatchData(patch, d_viscosity, context));

  FORCEINLINE_MACRO
  SAMRAI::pdat::parallel_for_all(patch.getBox(), [=] SAMRAI_HOST_DEVICE (int k, int j) {
    // RAJA::forallN<PatchPolicy>(isets->getCells(1), isets->getCells(0), [=] SAMRAI_HOST_DEVICE (int k, int j) {
    tags(j,k) = 0;
  });


  const Real tag_q = d_tag_q_threshold;
  const Real tag_density = d_tag_density_gradient;
  const Real tag_energy = d_tag_energy_gradient;
  const Real tag_pressure = d_tag_pressure_gradient;

  if (d_tag_q) {
  FORCEINLINE_MACRO
          SAMRAI::pdat::parallel_for_all(patch.getBox(), [=] SAMRAI_HOST_DEVICE (int k, int j) {
      // RAJA::forallN<PatchPolicy>(isets->getCells(1), isets->getCells(0), [=] SAMRAI_HOST_DEVICE (int k, int j) {
      if (viscosity(j,k) > tag_q) {
        tags(j,k) = 1;
      }
    });
  }

  if (d_tag_density) {
  FORCEINLINE_MACRO
          SAMRAI::pdat::parallel_for_all(patch.getBox(), [=] SAMRAI_HOST_DEVICE (int k, int j) {
      // RAJA::forallN<PatchPolicy>(isets->getCells(1), isets->getCells(0), [=] SAMRAI_HOST_DEVICE (int k, int j) {
      const Real d2x = fabs(density(j+1,k) - 2.0*density(j,k) + density(j-1,k));
      const Real d2y = fabs(density(j,k+1) - 2.0*density(j,k) + density(j,k-1));

      const Real dxy = fabs(density(j+1,k+1) - 2.0*density(j,k) + density(j-1,k-1));
      const Real dyx = fabs(density(j-1,k+1) - 2.0*density(j,k) + density(j+1,k-1));

      const Real dd = MAX(d2x,MAX(d2y,MAX(dxy,dyx)));

      if (dd > tag_density) { tags(j,k) = 1; }
    });
  }

  if (d_tag_energy) {
  FORCEINLINE_MACRO
          SAMRAI::pdat::parallel_for_all(patch.getBox(), [=] SAMRAI_HOST_DEVICE (int k, int j) {
      // RAJA::forallN<PatchPolicy>(isets->getCells(1), isets->getCells(0), [=] SAMRAI_HOST_DEVICE (int k, int j) {
      const Real d2x = fabs(energy(j+1,k) - 2.0*energy(j,k) + energy(j-1,k));
      const Real d2y = fabs(energy(j,k+1) - 2.0*energy(j,k) + energy(j,k-1));

      const Real dxy = fabs(energy(j+1,k+1) - 2.0*energy(j,k) + energy(j-1,k-1));
      const Real dyx = fabs(energy(j-1,k+1) - 2.0*energy(j,k) + energy(j+1,k-1));

      const Real dd = MAX(d2x,MAX(d2y,MAX(dxy,dyx)));

      if (dd > tag_energy) { tags(j,k) = 1; }
    });
  }

  if (d_tag_pressure) {
  FORCEINLINE_MACRO
          SAMRAI::pdat::parallel_for_all(patch.getBox(), [=] SAMRAI_HOST_DEVICE (int k, int j) {
      // RAJA::forallN<PatchPolicy>(isets->getCells(1), isets->getCells(0), [=] SAMRAI_HOST_DEVICE (int k, int j) {
      const Real d2x = fabs(pressure(j+1,k) - 2.0*pressure(j,k) + pressure(j-1,k));
      const Real d2y = fabs(pressure(j,k+1) - 2.0*pressure(j,k) + pressure(j,k-1));

      const Real dxy = fabs(pressure(j+1,k+1) - 2.0*pressure(j,k) + pressure(j-1,k-1));
      const Real dyx = fabs(pressure(j-1,k+1) - 2.0*pressure(j,k) + pressure(j+1,k-1));

      const Real dd = MAX(d2x,MAX(d2y,MAX(dxy,dyx)));

      if (dd > tag_pressure) { tags(j,k) = 1; }
    });
  }

  if(d_tag_all) {
  FORCEINLINE_MACRO
          SAMRAI::pdat::parallel_for_all(patch.getBox(), [=] SAMRAI_HOST_DEVICE (int k, int j) {
      // RAJA::forallN<PatchPolicy>(isets->getCells(1), isets->getCells(0), [=] SAMRAI_HOST_DEVICE (int k, int j) {
      tags(j,k) = 1;
    });
  }

  auto samrai_tags = SAMRAI::pdat::get_view<Dim, IndicatorData>(samrai_tags_data);

  FORCEINLINE_MACRO
  SAMRAI::pdat::parallel_for_all(patch.getBox(), [=] SAMRAI_HOST_DEVICE (int k, int j) {
    // RAJA::forallN<PatchPolicy>(
    //     isets->getCells(1,false),
    //     isets->getCells(0,false), [=] SAMRAI_HOST_DEVICE (int k, int j) {
    // int tag_index = SAMRAI_TAGS->getGhostBox().offset(hier::Index(j,k));
    //SAMRAI_TAGS->getPointer()[tag_index] = tags(j,k);
    samrai_tags(j,k) = tags(j,k);
  });
}

void Cleverleaf::fillLevelIndicator(
  SAMRAI::hier::Patch& patch,
  const int level_number)
{
  auto context = getCurrentDataContext();
  auto indicator = SAMRAI::pdat::get_view<Dim, IndicatorData>(getPatchData(patch, d_level_indicator, context));
  FORCEINLINE_MACRO
  SAMRAI::pdat::parallel_for_all(patch.getBox(), [=] SAMRAI_HOST_DEVICE (int k, int j) {
    // RAJA::forallN<PatchPolicy>(isets->getCells(1), isets->getCells(0), [=] SAMRAI_HOST_DEVICE (int k, int j) {
    indicator(j,k) = level_number;
  });
}

void Cleverleaf::initializeCallback()
{
  t_fill_boundary = tbox::TimerManager::getManager()->getTimer(
    "Cleverleaf::setPhysicalBoundaryConditions()");
}

void Cleverleaf::finalizeCallback()
{
  t_fill_boundary.reset();
}

void Cleverleaf::setBCs(
  std::shared_ptr<CellData> data,
  const SAMRAI::hier::Box& bbox,
  const int edge,
  const int ghosts)
{
  const int depth = data->getDepth();
  if (depth != 1) {
    SAMRAI::tbox::perr << "[ERROR] Depth != 1... " << std::endl;
    exit(-1);
  }

  auto field = SAMRAI::pdat::get_view<Dim, CellData>(data);

  const SAMRAI::hier::Box& box = data->getBox();

  switch(edge) {
  case BdryLoc::YLO:
  {
    const int y_min = box.lower(Coord::Y);
    // const int y_min = isets->getMin(1);
  FORCEINLINE_MACRO
    SAMRAI::pdat::parallel_for_all(bbox, Coord::X, [=] SAMRAI_HOST_DEVICE (int j) {
      // RAJA::forall<BoundaryPolicy>(isets->getBoundaryCells(0, depth-1), [=] SAMRAI_HOST_DEVICE (int j) {
      field(j,y_min-1) = field(j,y_min);
      if (ghosts == 2) {
        field(j,y_min-2) = field(j,y_min+1);
      }
    });
  }
  break;
  case BdryLoc::YHI:
  {
    const int y_max = box.upper(Coord::Y);
    // const int y_max = isets->getMax(1);
  FORCEINLINE_MACRO
    SAMRAI::pdat::parallel_for_all(bbox, Coord::X, [=] SAMRAI_HOST_DEVICE (int j) {
      // RAJA::forall<BoundaryPolicy>(isets->getBoundaryCells(0, depth-1), [=] SAMRAI_HOST_DEVICE (int j) {
      field(j,y_max+1) = field(j,y_max);
      if (ghosts == 2) {
        field(j,y_max+2) = field(j,y_max-1);
      }
    });
  }
  break;
  case BdryLoc::XLO:
  {
    const int x_min = box.lower(Coord::X);
    // const int x_min = isets->getMin(0);
  FORCEINLINE_MACRO
    SAMRAI::pdat::parallel_for_all(bbox, Coord::Y, [=] SAMRAI_HOST_DEVICE (int k) {
      // RAJA::forall<BoundaryPolicy>(isets->getBoundaryCells(1, depth-1), [=] SAMRAI_HOST_DEVICE (int k) {
      field(x_min-1,k) = field(x_min,k);
      if (ghosts == 2) {
        field(x_min-2,k) = field(x_min+1,k);
      }
    });
  }
  break;
  case BdryLoc::XHI:
  {
    const int x_max = box.upper(Coord::X);
    // const int x_max = isets->getMax(0);
  FORCEINLINE_MACRO
    SAMRAI::pdat::parallel_for_all(bbox, Coord::Y, [=] SAMRAI_HOST_DEVICE (int k) {
      // RAJA::forall<BoundaryPolicy>(isets->getBoundaryCells(1, depth-1), [=] SAMRAI_HOST_DEVICE (int k) {
      field(x_max+1,k) = field(x_max,k);
      if (ghosts == 2) {
        field(x_max+2,k) = field(x_max-1,k);
      }
    });
  }
  break;
  }
}

void Cleverleaf::setBCs(
  std::shared_ptr<NodeData> data,
  const SAMRAI::hier::Box& bbox,
  const int edge,
  const int ghosts)
{

        auto field0 = SAMRAI::pdat::get_view<Dim, NodeData>(data, Coord::X);
        auto field1 = SAMRAI::pdat::get_view<Dim, NodeData>(data, Coord::Y);

  const SAMRAI::hier::Box& box = data->getBox();

  switch(edge) {
  case BdryLoc::YLO:
  {
    const int y_min = box.lower(Coord::Y);
  FORCEINLINE_MACRO
    SAMRAI::pdat::parallel_for_all(bbox, Coord::X, [=] SAMRAI_HOST_DEVICE (int j) {
      // RAJA::forall<BoundaryPolicy>(isets->getBoundaryNodes(0, depth-1), [=] SAMRAI_HOST_DEVICE (int j) {
      field0(j,y_min-1) =  1.0*field0(j,y_min+1);
      field1(j,y_min-1) = -1.0*field1(j,y_min+1);
      if (ghosts == 2) {
        field0(j,y_min-2) =  1.0*field0(j,y_min+2);
        field1(j,y_min-2) = -1.0*field1(j,y_min+2);
      }
    });
  }
  break;
  case BdryLoc::YHI:
  {
    const int y_max = box.upper(Coord::Y);
  FORCEINLINE_MACRO
    SAMRAI::pdat::parallel_for_all(bbox, Coord::X, [=] SAMRAI_HOST_DEVICE (int j) {
      // RAJA::forall<BoundaryPolicy>(isets->getBoundaryNodes(0, depth-1), [=] SAMRAI_HOST_DEVICE (int j) {
      field0(j,y_max+2) =  1.0*field0(j,y_max);
      field1(j,y_max+2) = -1.0*field1(j,y_max);
      if (ghosts == 2) {
        field0(j,y_max+3) =  1.0*field0(j,y_max-1);
        field1(j,y_max+3) = -1.0*field1(j,y_max-1);
      }
    });
  }
  break;
  case BdryLoc::XLO:
  {
    const int x_min = box.lower(Coord::X);
  FORCEINLINE_MACRO
    SAMRAI::pdat::parallel_for_all(bbox, Coord::Y, [=] SAMRAI_HOST_DEVICE (int k) {
      // RAJA::forall<BoundaryPolicy>(isets->getBoundaryNodes(1, depth-1), [=] SAMRAI_HOST_DEVICE (int k) {
      field0(x_min-1,k) = -1.0*field0(x_min+1,k);
      field1(x_min-1,k) =  1.0*field1(x_min+1,k);
      if (ghosts == 2) {
        field0(x_min-2,k) = -1.0*field0(x_min+2,k);
        field1(x_min-2,k) =  1.0*field1(x_min+2,k);
      }
    });
  }
  break;
  case BdryLoc::XHI:
  {
    const int x_max = box.upper(Coord::X);
  FORCEINLINE_MACRO
    SAMRAI::pdat::parallel_for_all(bbox, Coord::Y, [=] SAMRAI_HOST_DEVICE (int k) {
      // RAJA::forall<BoundaryPolicy>(isets->getBoundaryNodes(1, depth-1), [=] SAMRAI_HOST_DEVICE (int k) {
      field0(x_max+2,k) = -1.0*field0(x_max,k);
      field1(x_max+2,k) =  1.0*field1(x_max,k);
      if (ghosts == 2) {
        field0(x_max+3,k) = -1.0*field0(x_max-1,k);
        field1(x_max+3,k) =  1.0*field1(x_max-1,k);
      }
    });
  }
  break;
  }
}

void Cleverleaf::setBCs(
  std::shared_ptr<SideData> data,
  const SAMRAI::hier::Box& bbox0,
  const SAMRAI::hier::Box& bbox1,
  const int edge,
  const int ghosts)
{
        auto field0 = SAMRAI::pdat::get_view<Dim, SideData>(data, Side::X);
        auto field1 = SAMRAI::pdat::get_view<Dim, SideData>(data, Side::Y);

  const SAMRAI::hier::Box& box = data->getBox();

  switch(edge) {
  case BdryLoc::YLO:
  {
    const int y_min = box.lower(Coord::Y);
  FORCEINLINE_MACRO
    SAMRAI::pdat::parallel_for_all(bbox0, Coord::X, [=] SAMRAI_HOST_DEVICE (int j) {
      // RAJA::forall<BoundaryPolicy>(isets->getBoundaryNodes(0, depth-1), [=] SAMRAI_HOST_DEVICE (int j) {
      field0(j,y_min-1) = field0(j,y_min);
      if (ghosts == 2) {
        field0(j,y_min-2) = field0(j,y_min+1);
      }
    });
  FORCEINLINE_MACRO
    SAMRAI::pdat::parallel_for_all(bbox1, Coord::X, [=] SAMRAI_HOST_DEVICE (int j) {
      // RAJA::forall<BoundaryPolicy>(isets->getBoundaryCells(0, depth-1), [=] SAMRAI_HOST_DEVICE (int j) {
      field1(j,y_min-1) = -1.0*field1(j,y_min+1);
      if (ghosts == 2) {
        field1(j,y_min-2) = -1.0*field1(j,y_min+2);
      }
    });
  }
  break;
  case BdryLoc::YHI:
  {
    const int y_max = box.upper(Coord::Y);
  FORCEINLINE_MACRO
    SAMRAI::pdat::parallel_for_all(bbox0, Coord::X, [=] SAMRAI_HOST_DEVICE (int j) {
      // RAJA::forall<BoundaryPolicy>(isets->getBoundaryNodes(0, depth-1), [=] SAMRAI_HOST_DEVICE (int j) {
      field0(j,y_max+1) = field0(j,y_max);
      if (ghosts == 2) {
        field0(j,y_max+2) = field0(j,y_max-1);
      }
    });
  FORCEINLINE_MACRO
    SAMRAI::pdat::parallel_for_all(bbox1, Coord::X, [=] SAMRAI_HOST_DEVICE (int j) {
      // RAJA::forall<BoundaryPolicy>(isets->getBoundaryCells(0, depth-1), [=] SAMRAI_HOST_DEVICE (int j) {
      field1(j,y_max+2) = -1.0*field1(j,y_max);
      if (ghosts == 2) {
        field1(j,y_max+3) = -1.0*field1(j,y_max-1);
      }
    });
  }
  break;
  case BdryLoc::XLO:
  {
    const int x_min = box.lower(Coord::X);
  FORCEINLINE_MACRO
    SAMRAI::pdat::parallel_for_all(bbox0, Coord::Y, [=] SAMRAI_HOST_DEVICE (int k) {
      // RAJA::forall<BoundaryPolicy>(isets->getBoundaryCells(1, depth-1), [=] SAMRAI_HOST_DEVICE (int k) {
      field0(x_min-1,k) = -1.0*field0(x_min+1,k);
      if (ghosts == 2) {
        field0(x_min-2,k) = -1.0*field0(x_min+2,k);
      }
    });
  FORCEINLINE_MACRO
    SAMRAI::pdat::parallel_for_all(bbox1, Coord::Y, [=] SAMRAI_HOST_DEVICE (int k) {
      // RAJA::forall<BoundaryPolicy>(isets->getBoundaryNodes(1, depth-1), [=] SAMRAI_HOST_DEVICE (int k) {
      field1(x_min-1,k) = field1(x_min,k);
      if (ghosts == 2) {
        field1(x_min-2,k) = field1(x_min+1,k);
      }
    });
  }
  break;
  case BdryLoc::XHI:
  {
    const int x_max = box.upper(Coord::X);
  FORCEINLINE_MACRO
    SAMRAI::pdat::parallel_for_all(bbox0, Coord::Y, [=] SAMRAI_HOST_DEVICE (int k) {
      // RAJA::forall<BoundaryPolicy>(isets->getBoundaryCells(1, depth-1), [=] SAMRAI_HOST_DEVICE (int k) {
      field0(x_max+2,k) = -1.0*field0(x_max,k);
      if (ghosts == 2) {
        field0(x_max+3,k) = -1.0*field0(x_max-1,k);
      }
    });
  FORCEINLINE_MACRO
    SAMRAI::pdat::parallel_for_all(bbox1, Coord::Y, [=] SAMRAI_HOST_DEVICE (int k) {
      // RAJA::forall<BoundaryPolicy>(isets->getBoundaryNodes(1, depth-1), [=] SAMRAI_HOST_DEVICE (int k) {
      field1(x_max+1,k) = field1(x_max,k);
      if (ghosts == 2) {
        field1(x_max+2,k) = field1(x_max-1,k);
      }
    });
  }
  break;
  }
}

void Cleverleaf::debug_knl(hier::Patch& patch)
{
  auto current_context = getCurrentDataContext();
  auto new_context = getNewDataContext();

  auto density0 = getPatchData(patch, d_density, current_context);
  auto density1 = getPatchData(patch, d_density, new_context);

  auto energy0 = getPatchData(patch, d_energy, current_context);
  auto energy1 = getPatchData(patch, d_energy, new_context);

  auto pressure = getPatchData(patch, d_pressure, current_context);
  auto soundspeed = getPatchData(patch, d_soundspeed, current_context);
  auto viscosity = getPatchData(patch, d_viscosity, current_context);
  auto volume = getPatchData(patch, d_volume, current_context);

  auto velocity0 = getPatchData(patch, d_velocity, current_context);
  auto velocity1 = getPatchData(patch, d_velocity, new_context);

  auto massflux = getPatchData(patch, d_massflux, current_context);
  auto volflux = getPatchData(patch, d_volflux, current_context);

  std::cout << "density0" << std::endl << std::flush;
  debug(patch, density0);
  SAMRAI::tbox::parallel_synchronize();
  std::cout << "density1" << std::endl << std::flush;
  debug(patch, density1);
  SAMRAI::tbox::parallel_synchronize();
  std::cout << "energy0" << std::endl << std::flush;
  debug(patch, energy0);
  SAMRAI::tbox::parallel_synchronize();
  std::cout << "energy1" << std::endl << std::flush;
  debug(patch, energy1);
  SAMRAI::tbox::parallel_synchronize();
  std::cout << "pressure" << std::endl << std::flush;
  debug(patch, pressure);
  SAMRAI::tbox::parallel_synchronize();
  std::cout << "soundspeed" << std::endl << std::flush;
  debug(patch, soundspeed);
  SAMRAI::tbox::parallel_synchronize();
  std::cout << "viscosity" << std::endl << std::flush;
  debug(patch, viscosity);
  SAMRAI::tbox::parallel_synchronize();
  std::cout << "volume" << std::endl << std::flush;
  debug(patch, volume);
  SAMRAI::tbox::parallel_synchronize();
  std::cout << "velocity0" << std::endl << std::flush;
  debug(patch, velocity0);
  SAMRAI::tbox::parallel_synchronize();
  std::cout << "velocity1" << std::endl << std::flush;
  debug(patch, velocity1);
  SAMRAI::tbox::parallel_synchronize();
  std::cout << "massflux" << std::endl << std::flush;
  debug(patch, massflux);
  SAMRAI::tbox::parallel_synchronize();
  std::cout << "volflux" << std::endl << std::flush;
  debug(patch, volflux);
  SAMRAI::tbox::parallel_synchronize();
}

void Cleverleaf::debug(
  SAMRAI::hier::Patch& patch,
  std::shared_ptr<CellData> data)
{
  //for (int i = 0; i < data->getDepth(); i++) {
  //  auto field = SAMRAI::pdat::get_view<Dim, CellData>(data, i);
  //  SAMRAI::pdat::for_all<tbox::policy::sequential>(getIndexBox<true>(data), [=] SAMRAI_HOST_DEVICE (int k, int j) {
  //    // RAJA::forallN<PatchPolicy>(isets->getCells(1,true), isets->getCells(0,true), [=] SAMRAI_HOST_DEVICE (int k, int j) {
  //    printf("(%d,%d) = %f\n", j,k, field(j,k) );
  //    if (field(j,k) != field(j,k)) {
  //      printf("here\n");
  //    }
  //  });
  //  std::cout << "---------------" << std::endl;
  //}
}

void Cleverleaf::debug(
  SAMRAI::hier::Patch& patch,
  std::shared_ptr<NodeData> data)
{
  //for (int i = 0; i < data->getDepth(); i++) {
  //  auto field = SAMRAI::pdat::get_view<Dim, NodeData>(data, i);
  //  SAMRAI::pdat::for_all<tbox::policy::sequential>(getIndexBox<true>(data), [=] SAMRAI_HOST_DEVICE (int k, int j) {
  //    // RAJA::forallN<PatchPolicy>(isets->getCells(1,true), isets->getCells(0,true), [=] SAMRAI_HOST_DEVICE (int k, int j) {
  //    printf("(%d,%d) = %f\n", j,k, field(j,k) );
  //    if (field(j,k) != field(j,k)) {
  //      printf("here\n");
  //    }
  //  });
  //  std::cout << "---------------" << std::endl;
  //}
}

void Cleverleaf::debug(
  SAMRAI::hier::Patch& patch,
  std::shared_ptr<SideData> data)
{
  //for (int i = 0; i < data->getDepth(); i++) {
  //  for (int j = Side::X; j <= Side::Y; j++) {
  //    auto field = SAMRAI::pdat::get_view<Dim, SideData>(data, j, i);
  //    SAMRAI::pdat::for_all<tbox::policy::sequential>(getIndexBox<true>(data, j), [=] SAMRAI_HOST_DEVICE (int k, int j) {
  //      // RAJA::forallN<PatchPolicy>(isets->getCells(1,true), isets->getCells(0,true), [=] SAMRAI_HOST_DEVICE (int k, int j) {
  //      printf("(%d,%d) = %f\n", j,k, field(j,k) );
  //      if (field(j,k) != field(j,k)) {
  //        printf("here\n");
  //      }
  //    });
  //    std::cout << "---------------" << std::endl;
  //  }
  //}
}

}
}
