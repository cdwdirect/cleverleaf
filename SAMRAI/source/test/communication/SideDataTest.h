/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2018 Lawrence Livermore National Security, LLC
 * Description:   AMR communication tests for side-centered patch data
 *
 ************************************************************************/

#ifndef included_pdat_SideDataTest
#define included_pdat_SideDataTest

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "PatchDataTestStrategy.h"
#ifndef included_String
#include <string>
using namespace std;
#define included_String
#endif
#include "SAMRAI/hier/Variable.h"

#include <memory>

namespace SAMRAI {

class CommTester;

/**
 * Class SideDataTest provides routines to test communication operations
 * for side-centered patch data on an AMR patch hierarchy.
 *
 * Required input keys and data types for test:
 *
 *   NONE...
 *
 * See PatchDataTestStrategy header file comments for variable and
 * refinement input data description.  Additionally, there are two
 * optional input parameters for each side variable.  These are:
 *
 *
 *
 *
 *
 *    - \b  test_direction   side directions to test
 *                             (default = -1 ie, all directions)
 *    - \b  use_fine_value_at_interface   which values to use at coarse-
 *                                          fine interface (default = TRUE)
 *
 *
 *
 *
 *
 */

class SideDataTest:public PatchDataTestStrategy
{
public:
   /**
    * The constructor initializes variable data arrays to zero length.
    */
   SideDataTest(
      const std::string& object_name,
      const tbox::Dimension& dim,
      std::shared_ptr<tbox::Database> main_input_db,
      bool do_refine,
      bool do_coarsen,
      const std::string& refine_option);

   /**
    * Virtual destructor for SideDataTest.
    */
   ~SideDataTest();

   /**
    * User-supplied boundary conditions.  Note that we do not implement
    * user-defined coarsen and refine operations.
    */
   virtual void
   setPhysicalBoundaryConditions(
      const hier::Patch& patch,
      const double time,
      const hier::IntVector& gcw) const;

   /**
    * This function is called from the CommTester constructor.  Its
    * purpose is to register variables used in the patch data test
    * and appropriate communication parameters (ghost cell widths,
    * coarsen/refine operations) with the CommTester object, which
    * manages the variable storage.
    */
   void
   registerVariables(
      CommTester* commtest);

   /**
    * Function for setting data on new patch in hierarchy.
    *
    * @param src_or_dst Flag set to 's' for source or 'd' for destination
    *        to indicate variables to set data for.
    */
   virtual void
   initializeDataOnPatch(
      const hier::Patch& patch,
      const std::shared_ptr<hier::PatchHierarchy> hierarchy,
      int level_number,
      char src_or_dst);

   /**
    * Function for checking results of communication operations.
    */
   bool
   verifyResults(
      const hier::Patch& patch,
      const std::shared_ptr<hier::PatchHierarchy> hierarchy,
      int level_number);

private:
   /*
    * Function for reading test data from input file.
    */
   void
   readTestInput(
      std::shared_ptr<tbox::Database> db);

   void
   setLinearData(
      std::shared_ptr<pdat::SideData<double> > data,
      const hier::Box& box,
      const hier::Patch& patch) const;

   /**
    * Set periodic linear function data for testing interpolation in
    * periodic domains.
    */
   void
   setPeriodicData(
      std::shared_ptr<pdat::SideData<double> > data,
      const hier::Box& box,
      const hier::Patch& patch) const;

   /*
    * Set conservative linear function data for testing coarsening
    */
   void
   setConservativeData(
      std::shared_ptr<pdat::SideData<double> > data,
      const hier::Box& box,
      const hier::Patch& patch,
      const std::shared_ptr<hier::PatchHierarchy> hierarchy,
      int level_number) const;

   void
   checkPatchInteriorData(
      const std::shared_ptr<pdat::SideData<double> >& data,
      const hier::Box& interior,
      const hier::Patch& patch) const;

   const tbox::Dimension d_dim;

   /*
    * Object std::string identifier for error reporting
    */
   std::string d_object_name;

   /*
    * Data members specific to this side data test.
    */
   std::shared_ptr<geom::CartesianGridGeometry> d_cart_grid_geometry;

   /*
    * Data members specific to this side data test.
    */
   std::vector<hier::IntVector> d_test_direction;
   std::vector<bool> d_use_fine_value_at_interface;

   double d_Acoef;
   double d_Bcoef;
   double d_Ccoef;
   double d_Dcoef;

   bool d_do_refine;
   bool d_do_coarsen;
   std::string d_refine_option;
   int d_finest_level_number;

   std::vector<std::shared_ptr<hier::Variable> > d_variables;

};

}
#endif
