/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2018 Lawrence Livermore National Security, LLC
 * Description:   Blueprint utilities.
 *
 ************************************************************************/
#ifndef included_hier_BlueprintUtils
#define included_hier_BlueprintUtils

#include "SAMRAI/SAMRAI_config.h"

#ifdef HAVE_CONDUIT
#include "SAMRAI/hier/BlueprintUtilsStrategy.h"
#include "SAMRAI/tbox/Database.h"

#include <memory>

namespace SAMRAI {
namespace hier {

/*!
 * @brief Utilities for performing common tasks for describing the mesh
 * in the Conduit blueprint format and writing to files.
 *
 * See https://llnl-conduit.readthedocs.io for documentation of the Conduit library.
 *
 * This class can use a BlueprintUtilsStrategy to allow a call-back to
 * application code to fill a blueprint with problem-specific coordinate
 * information.
 *
 * Combined with a BlueprintUtilsStrategy, this loops over a hierarchy,
 * fills Blueprint topology and coordinate entries, calls back to user code
 * for specific coord choices:  uniform, rectilinear, explicit
 */

class PatchHierarchy;

class BlueprintUtils
{

public:

   /*!
    * @brief Constructor
    *
    * @param strategy    Strategy pointer for callback to application code.
    */
   BlueprintUtils(BlueprintUtilsStrategy* strategy);

   /*!
    * @brief Destructor
    */
   virtual ~BlueprintUtils();

   /*!
    * @brief Put topology and coordinates to the database
    *
    * Combined with a BlueprintUtilsStrategy, this loops over a hierarchy,
    * fills blueprint topology and coordinate entries, calls back to user code
    * for specific types of coordinates recognized by the blueprint:
    * uniform, rectilinear, or explicit
    *
    * @param blueprint_db  Top-level blueprint database holding all local
    *                      domain information
    * @param hierarchy     The hierarchy being described
    * @param topology_name Name of the topology
    */
   void putTopologyAndCoordinatesToDatabase(
      const std::shared_ptr<tbox::Database>& blueprint_db,
      const PatchHierarchy& hierarchy,
      const std::string& topology_name) const;

   /*!
    * @brief Write blueprint to files.
    *
    * Given a conduit node that holds a blueprint, write it in parallel files.
    * There will be a root file holding the blueprint index and a data files
    * for each rank holding all of the parallel domain data from each rank's
    * local patches.
    *
    * The output files will be {rootfile_name}.root and
    * {data_name}NNNNNN.{io_protocol}, with the Ns representing a six-digit
    * rank number.
    *
    * @param blueprint  Node holding the local blueprint on all ranks
    * @param samrai_mpi SAMRAI_MPI object containing all ranks in use
    * @param num_global_domains  Global number of domains (patches) in the
    *                            hierarchy described by the blueprint
    * @param mesh_name  A name for the mesh
    * @param data_name  File name stub for the files holding the domain data
    * @param rootfile_name  File name sub for the root file.
    * @param io_protocol    I/O protocol string identifier--must be a valid
    *                       protocol as described in Conduit documentation
    */
   void writeBlueprintMesh(
      const conduit::Node& blueprint,
      const tbox::SAMRAI_MPI& samrai_mpi,
      const int num_global_domains,
      const std::string& mesh_name,
      const std::string& data_name,
      const std::string& rootfile_name,
      const std::string& io_protocol) const;

private:

   BlueprintUtilsStrategy* d_strategy;

};

}
}

#endif // HAVE_CONDUIT

#endif  // included_hier_BlueprintUtils
