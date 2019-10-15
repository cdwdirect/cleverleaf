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
#include <iostream>
#include <cstdlib>
#include <string>
#include <fstream>
#include <unistd.h>
#include <sys/time.h>

#include "hydro/Cleverleaf.h"
#include "hydro/LagrangianEulerianLevelIntegrator.h"
#include "hydro/LagrangianEulerianIntegrator.h"

#include "SAMRAI/mesh/BergerRigoutsos.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/mesh/ChopAndPackLoadBalancer.h"
#include "SAMRAI/mesh/TreeLoadBalancer.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/mesh/StandardTagAndInitialize.h"
#include "SAMRAI/appu/VisItDataWriter.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/tbox/TimerManager.h"

#ifdef ENABLE_APOLLO
#include "apollo/Apollo.h"
#ifndef VERSION
#define VERSION   "Apollo"
#endif
#define APOLLO_TIME(__APOLLO_dbl_var)                                \
    {                                                                \
        struct timeval t;                                            \
        gettimeofday(&t, NULL);                                      \
        __APOLLO_dbl_var = (double)(t.tv_sec + (t.tv_usec * 1e-6));  \
    }
#endif

#ifndef VERSION
#define VERSION   "Normal"
#endif
#ifndef HOST_NAME
#define HOST_NAME "localhost"
#endif

using namespace SAMRAI;
using namespace clever;

static std::string getFilenameFromPath(const std::string& path);
static std::string createUniqueFilename(const std::string& name, const std::string& extension);
static bool fileExists(const std::string& path);

int main(int argc, char* argv[]) {
  tbox::SAMRAI_MPI::init(&argc, &argv);
  tbox::SAMRAIManager::initialize();
  tbox::SAMRAIManager::startup();

  const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());

  {
    if (argc == 1) {
      tbox::pout << "USAGE: " << argv[0]
                 << " [-l max_level] [-r refinement_ratio] [-f regrid_frequency]"
                 << " <input filename>" << std::endl;
      return 1;
    }

    std::string input_path;

    int max_level_number = -1;
    int refinement_ratio = -1;
    int regrid_interval = -1;

    int opt;

    while ((opt = getopt(argc,argv,"l:r:f:i:")) != EOF ) {
      switch (opt) {
      case 'l':
        max_level_number = atoi(optarg);
        tbox::plog << "max_level_number overide using command-line: "
                   << max_level_number << std::endl;
        break;
      case 'r':
        refinement_ratio = atoi(optarg);
        tbox::plog << "refinement_ratio overide using command-line: "
                   << refinement_ratio << std::endl;
        break;
      case 'f':
        regrid_interval = atoi(optarg);
        tbox::plog << "regrid_interval overide using command-line: "
                   << regrid_interval << std::endl;
        break;
      case 'i':
        input_path = optarg;
        break;
      case '?':
        if (optopt == 'i')
          tbox::perr << "Option -i requires input file argument" << std::endl;
        else
          tbox::perr << "Unknown option " << optopt << std::endl;

        tbox::pout << "USAGE: " << argv[0]
                   << " [-l max_level] [-r refinement_ratio] [-f regrid_frequency]"
                   << " <input filename>" << std::endl;

        return 1;
      }
    }

    if (input_path.empty() && optind < argc) {
      input_path = argv[optind];
    } else {
      tbox::perr << "Input file required" << std::endl;
      tbox::pout << "USAGE: " << argv[0]
                 << " [-l max_level] [-r refinement_ratio] [-f regrid_frequency]"
                 << " <input filename>" << std::endl;
      return 1;
    }

    //tbox::pout << "CleverLeaf version " << VERSION
    //           << " compiled on " << HOST_NAME << std::endl;
    //tbox::pout << "Running with " << mpi.getSize() << " tasks" << std::endl;
//#if defined(_OPENMP)
//#pragma omp parallel
//    {
//#pragma omp master
//      { tbox::pout << " and " << omp_get_num_threads() << " of ";
//          tbox::pout << omp_get_max_threads() << " threads";
//      }
//    }
//#endif
//    tbox::pout << std::endl;

    tbox::plog << "Reading input from: " << input_path << std::endl;

    std::shared_ptr<tbox::InputDatabase> input_db(
      new tbox::InputDatabase("input_db"));
    tbox::InputManager::getManager()->parseInputFile(input_path, input_db);

    std::shared_ptr<tbox::Database> main_db = input_db->getDatabase(
      "Cleverleaf");

    const std::string input_file = getFilenameFromPath(input_path);
    const std::string basename = main_db->getStringWithDefault("basename", input_file);

    const tbox::Dimension dim(
      static_cast<unsigned short>(main_db->getInteger("dim")));

    int vis_dump_interval = main_db->getIntegerWithDefault(
      "vis_dump_interval", 1);
    int field_summary_interval = main_db->getIntegerWithDefault(
      "field_summary_interval", 10);
    int restart_interval = main_db->getIntegerWithDefault(
      "restart_interval", 1);

    const std::string log_basename = main_db->getStringWithDefault(
      "log_filename", basename);
    const std::string log_filename = createUniqueFilename(log_basename, "log");

    bool log_all_nodes = false;
    log_all_nodes = main_db->getBoolWithDefault("log_all_nodes", log_all_nodes);

    if (log_all_nodes) {
      tbox::PIO::logAllNodes(log_filename);
    } else {
      tbox::PIO::logOnlyNodeZero(log_filename);
    }

    if(input_db->isDatabase("TimerManager")) {
      tbox::TimerManager::createManager(input_db->getDatabase("TimerManager"));
    }

    std::shared_ptr<SAMRAI::geom::CartesianGridGeometry> grid_geometry(
      new SAMRAI::geom::CartesianGridGeometry(
        dim,
        "CartesianGeometry",
        input_db->getDatabase("CartesianGeometry")));

    /*
     * Overwrite max_levels and ratio_to_coarser if passed as command line
     * arguments.
     */
    if(max_level_number != -1) {
      input_db->getDatabase("PatchHierarchy")
        ->putInteger("max_levels", max_level_number);
    }

    if(refinement_ratio != -1) {
      std::vector<int> ratio_vector(dim.getValue());

      for(int i = 0; i < dim.getValue(); i++) {
        ratio_vector[i] = refinement_ratio;
      }

      input_db->getDatabase("PatchHierarchy")
        ->getDatabase("ratio_to_coarser")
        ->putIntegerVector("level_1", ratio_vector);
    }

    std::shared_ptr<hier::PatchHierarchy> patch_hierarchy(
      new hier::PatchHierarchy(
        "PatchHierarchy",
        grid_geometry,
        input_db->getDatabase("PatchHierarchy")));

    hydro::Cleverleaf* cleverleaf = new hydro::Cleverleaf(
      main_db,
      patch_hierarchy,
      dim,
      grid_geometry);

    std::shared_ptr<hydro::LagrangianEulerianLevelIntegrator>
      lagrangian_eulerian_level_integrator(
        new hydro::LagrangianEulerianLevelIntegrator(
          input_db->getDatabase("LagrangianEulerianLevelIntegrator"),
          cleverleaf));

    std::shared_ptr<mesh::StandardTagAndInitialize> error_detector(
      new mesh::StandardTagAndInitialize(
        "StandardTagAndInitialize",
        lagrangian_eulerian_level_integrator.get(),
        input_db->getDatabase("StandardTagAndInitialize")));

    std::shared_ptr<mesh::BergerRigoutsos> box_generator(
      new mesh::BergerRigoutsos(
        dim,
        input_db->getDatabaseWithDefault(
          "BergerRigoutsos",
          std::shared_ptr<SAMRAI::tbox::Database>())));

    bool DEV_use_chop_and_pack = false;
    DEV_use_chop_and_pack = main_db->getBoolWithDefault(
      "DEV_use_chop_and_pack",
      DEV_use_chop_and_pack);

    std::shared_ptr<mesh::LoadBalanceStrategy> load_balancer;

    if (DEV_use_chop_and_pack) {
      load_balancer.reset(new mesh::ChopAndPackLoadBalancer(
                            dim,
                            "LoadBalancer",
                            input_db->getDatabase("LoadBalancer")));
    } else {
      mesh::TreeLoadBalancer* tree_load_balancer = new mesh::TreeLoadBalancer(
        dim,
        "LoadBalancer",
        input_db->getDatabase("LoadBalancer"));
      tree_load_balancer->setSAMRAI_MPI(SAMRAI::tbox::SAMRAI_MPI::getSAMRAIWorld());

      load_balancer.reset(tree_load_balancer);
    }

    std::shared_ptr<mesh::GriddingAlgorithm> gridding_algorithm(
      new mesh::GriddingAlgorithm(
        patch_hierarchy,
        "GriddingAlgorithm",
        input_db->getDatabase("GriddingAlgorithm"),
        error_detector,
        box_generator,
        load_balancer));

    /*
     * Overwrite regrid_interval option if necessary.
     */
    if(regrid_interval != -1) {
      input_db->getDatabase("LagrangianEulerianIntegrator")
        ->putInteger("regrid_interval", regrid_interval);
    }

    std::shared_ptr<hydro::LagrangianEulerianIntegrator>
      lagrangian_eulerian_integrator(
        new hydro::LagrangianEulerianIntegrator(
          input_db->getDatabase("LagrangianEulerianIntegrator"),
          patch_hierarchy,
          lagrangian_eulerian_level_integrator,
          gridding_algorithm));

    const int visit_number_procs_per_file = 1;
    const std::string visit_dirname = createUniqueFilename(basename, "visit");
    const std::string restart_write_dirname = createUniqueFilename(basename, "restart");

    std::shared_ptr<appu::VisItDataWriter> visit_data_writer(
      new appu::VisItDataWriter(
        dim,
        "Cleverleaf VisIt Writer",
        visit_dirname,
        visit_number_procs_per_file));

    cleverleaf->registerVisItDataWriter(visit_data_writer);

    Real dt_now = lagrangian_eulerian_integrator->initializeHierarchy();

    if(vis_dump_interval > 0) {
      visit_data_writer->writePlotData(patch_hierarchy, 0, 0.0);
    }

    tbox::plog << "\nCheck input data and variables before simulation:" << std::endl;
    tbox::plog << "Input database..." << std::endl;
    input_db->printClassData(tbox::plog);
    tbox::plog << "\nVariable database..." << std::endl;
    hier::VariableDatabase::getDatabase()->printClassData(tbox::plog);

    Real loop_time = lagrangian_eulerian_integrator->getIntegratorTime();
    const Real loop_time_end = lagrangian_eulerian_integrator->getEndTime();

    MPI_Barrier(MPI_COMM_WORLD);
    double start = MPI_Wtime();

#ifdef ENABLE_APOLLO
    double APOLLO_time_before_flush = 0.0;
    double APOLLO_time_after_flush  = 0.0;
    double APOLLO_time_this_step    = 0.0;
    double APOLLO_time_cumulative   = 0.0;
    Apollo *apollo = Apollo::instance();
#endif

    Real loop_time_start = loop_time;
    Real loop_time_stop  = loop_time;

    double step_exec_start = 0.0;
    double step_exec_total = 0.0;

#ifdef ENABLE_APOLLO
    tbox::pout << "CSV,build,policy_index,schedule,step,step_exec_time,sim_time_start,sim_time_done,current_dt,apollo_xmit_time" << std::endl;
#else
    tbox::pout << "CSV,build,omp_threads,schedule,step,step_exec_time,sim_time_start,sim_time_done,current_dt,apollo_xmit_time" << std::endl;
#endif
    while ((loop_time < loop_time_end) &&
           lagrangian_eulerian_integrator->stepsRemaining()) {

      step_exec_start = MPI_Wtime();
      loop_time_start = loop_time;

      int iteration_num = lagrangian_eulerian_integrator->getIntegratorStep()+1;

      Real dt_new = lagrangian_eulerian_integrator->advanceHierarchy(dt_now);

      loop_time += dt_now;
      dt_now = dt_new;

      loop_time_stop = loop_time;

      if ((field_summary_interval > 0)
          && (iteration_num % field_summary_interval) == 0) {
        lagrangian_eulerian_integrator->printFieldSummary();
      }

      step_exec_total = MPI_Wtime() - step_exec_start;

#ifdef ENABLE_APOLLO
      APOLLO_TIME(APOLLO_time_before_flush);

      // TODO: Explore the round-robin'ing of who is sending data in...
      //       (various ways of determining how often we send, to not send every timestep.)
      apollo->flushAllRegionMeasurements(iteration_num);

      APOLLO_TIME(APOLLO_time_after_flush);
      APOLLO_time_this_step = (APOLLO_time_after_flush - APOLLO_time_before_flush);
      APOLLO_time_cumulative += APOLLO_time_this_step;
      tbox::pout << "CSV,APOLLO" \
          << "," << getenv("APOLLO_POLICY_INDEX") \
          << ",apollo" \
          << "," << (iteration_num - 1) \
          << "," << step_exec_total \
          << "," << loop_time_start \
          << "," << loop_time_stop \
          << "," << dt_now \
          << "," << APOLLO_time_this_step \
          << std::endl;
#else
      tbox::pout << "CSV, NORMAL" \
          << "," << getenv("OMP_NUM_THREADS") \
          << "," << getenv("OMP_SCHEDULE") \
          << "," << (iteration_num - 1) \
          << "," << step_exec_total \
          << "," << loop_time_start \
          << "," << loop_time_stop \
          << "," << dt_now \
          << ",0.0" << std::endl;
#endif

      if ((vis_dump_interval > 0) && (iteration_num % vis_dump_interval) == 0) {
        tbox::pout << "== CLEVERLEAF: Dumping vis data..." << std::endl;
        visit_data_writer->writePlotData(
          patch_hierarchy,
          iteration_num,
          loop_time);
      }

      if ((restart_interval > 0) && (iteration_num % restart_interval) == 0) {
        tbox::pout << "== CLEVERLEAF: Dumping restart file..." << std::endl;
        tbox::RestartManager::getManager()->
          writeRestartFile(restart_write_dirname,
                           iteration_num);
      }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double end = MPI_Wtime();
    //SAMRAI::tbox::pout << "== CLEVERLEAF: Total time elapsed = " << (end - start) << std::endl;
#ifdef ENABLE_APOLLO
    //SAMRAI::tbox::pout << "== CLEVERLEAF: Total time spent flushing data to Apollo = " << APOLLO_time_cumulative << std::endl;
#endif

    /*
     * Write final visualization and restart dumps.
     */
    if ((vis_dump_interval > 0) &&
        (lagrangian_eulerian_integrator->getIntegratorStep()
         % vis_dump_interval) != 0)
    {
      visit_data_writer->writePlotData(
        patch_hierarchy,
        lagrangian_eulerian_integrator->getIntegratorStep(),
        loop_time);
    }

    if ((restart_interval) &&
        (lagrangian_eulerian_integrator->getIntegratorStep()
         % restart_interval) != 0)
    {
      tbox::RestartManager::getManager()->
        writeRestartFile(restart_write_dirname,
                         lagrangian_eulerian_integrator->getIntegratorStep());
    }

    /*
     * Write final field summary.
     */
    if ((field_summary_interval > 0)
        && (lagrangian_eulerian_integrator->getIntegratorStep()
            % field_summary_interval) != 0)
    {
      lagrangian_eulerian_integrator->printFieldSummary();
    }

    tbox::TimerManager::getManager()->print(tbox::plog);

    patch_hierarchy.reset();
    grid_geometry.reset();
    box_generator.reset();
    load_balancer.reset();
    error_detector.reset();
    gridding_algorithm.reset();

    delete cleverleaf;
  }

  tbox::SAMRAIManager::shutdown();
  tbox::SAMRAIManager::finalize();
  tbox::SAMRAI_MPI::finalize();

  return 0;
}

static std::string getFilenameFromPath(const std::string& path)
{
  std::size_t start = path.find_last_of("/");
  start += 1;

  std::size_t end = path.find_last_of(".");

  if(start == std::string::npos)
    start = 0;

  if(end == std::string::npos)
    end = path.length();

  std::size_t count = end - start;

  return path.substr(start, count);
}

static std::string createUniqueFilename(const std::string& name, const std::string& extension)
{
  int unique_id = -1;
  std::stringstream ss;
  std::string filename;

  do {
    ss.str("");
    ss.clear();
    unique_id++;
    ss << name << "." << unique_id << "." << extension;
    filename = ss.str();
  } while (fileExists(filename));

  return filename;
}

static inline bool fileExists(const std::string& path)
{
  std::ifstream ifile(path.c_str());
  return ifile.good();
}
