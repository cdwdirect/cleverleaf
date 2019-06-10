/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2018 Lawrence Livermore National Security, LLC
 * Description:   Default fill pattern class
 *
 ************************************************************************/

#ifndef included_xfer_BoxGeometryVariableFillPattern
#define included_xfer_BoxGeometryVariableFillPattern

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/xfer/VariableFillPattern.h"

namespace SAMRAI {
namespace xfer {

/*!
 * @brief Class BoxGeometryVariableFillPattern is a default implementation of
 * the abstract base class VariableFillPattern.
 *
 * It is used to calculate overlaps that consist of the full intersection
 * between source and destination patches, including all ghost regions.  If
 * no VariableFillPattern object is provided when a refine operation is
 * registered with a RefineAlgorithm, this class is used by default.
 *
 * @see RefineAlgorithm
 */

class BoxGeometryVariableFillPattern:
   public VariableFillPattern
{
public:
   /*!
    * @brief Default constructor
    */
   BoxGeometryVariableFillPattern();

   /*!
    * @brief Destructor
    */
   virtual ~BoxGeometryVariableFillPattern();

   /*!
    * @brief Calculate overlap between the destination and source geometries
    * using the geometries' own overlap calculation methods.
    *
    * The intersection between the given dst_geometry and src_geometry
    * will be calculated according to the properties of those geometries.
    *
    * @param[in] dst_geometry    geometry object for destination box
    * @param[in] src_geometry    geometry object for source box
    * @param[in] dst_patch_box   box for the destination patch
    * @param[in] src_mask        the source mask, the box resulting from
    *                            transforming the source box
    * @param[in] fill_box        the box to be filled
    * @param[in] overwrite_interior  controls whether or not to include the
    *                                destination box interior in the overlap.
    * @param[in] transformation  the transformation from source to
    *                            destination index space.
    *
    * @return                    std::shared_ptr to the calculated overlap
    *                            object
    *
    * @pre dst_patch_box.getDim() == src_mask.getDim()
    */
   std::shared_ptr<hier::BoxOverlap>
   calculateOverlap(
      const hier::BoxGeometry& dst_geometry,
      const hier::BoxGeometry& src_geometry,
      const hier::Box& dst_patch_box,
      const hier::Box& src_mask,
      const hier::Box& fill_box,
      const bool overwrite_interior,
      const hier::Transformation& transformation) const
   {
#ifndef DEBUG_CHECK_DIM_ASSERTIONS
      NULL_USE(dst_patch_box);
#endif
      TBOX_ASSERT_OBJDIM_EQUALITY2(dst_patch_box, src_mask);
      return dst_geometry.calculateOverlap(src_geometry, src_mask, fill_box,
         overwrite_interior, transformation);
   }

   /*!
    * Computes a BoxOverlap object which defines the space to be filled by
    * a refinement operation.  For this implementation, that space is the
    * intersection between fill_boxes (computed by the RefineSchedule) and
    * data_box, which specifies the extent of the destination data.  The
    * patch data factory is used to compute the overlap with the appropriate
    * data centering, consistent with the centering of the data to be filled.
    *
    * @param[in] fill_boxes  list representing the all of the space on a patch
    *                        or its ghost region that may be filled by a
    *                        refine operator (cell-centered representation)
    * @param[in] node_fill_boxes node-centered representation of fill_boxes
    * @param[in] patch_box   box representing the patch where a refine operator
    *                        will fill data.  (cell-centered representation)
    * @param[in] data_box    box representing the full extent of the region
    *                        covered by a patch data object, including all
    *                        ghosts (cell-centered representation)
    * @param[in] pdf         patch data factory for the data that is to be
    *                        filled
    */
   std::shared_ptr<hier::BoxOverlap>
   computeFillBoxesOverlap(
      const hier::BoxContainer& fill_boxes,
      const hier::BoxContainer& node_fill_boxes,
      const hier::Box& patch_box,
      const hier::Box& data_box,
      const hier::PatchDataFactory& pdf) const;

   /*!
    * @brief Implementation of interface to get stencil width of a
    * VariableFillPattern.
    *
    * For this class BoxGeometryVariableFillPattern, this method should
    * never be called, since overlaps are computed based on BoxGeometry
    * objects and not on any stencil.  An error will result if this method
    * is invoked.
    */
   const hier::IntVector&
   getStencilWidth();

   /*!
    * @brief Returns a string name identifier "BOX_GEOMETRY_FILL_PATTERN".
    */
   const std::string&
   getPatternName() const
   {
      return s_name_id;
   }

private:
   BoxGeometryVariableFillPattern(
      const BoxGeometryVariableFillPattern&);    // not implemented
   BoxGeometryVariableFillPattern&
   operator = (
      const BoxGeometryVariableFillPattern&);    // not implemented

   /*!
    * Static string containing string name identifier for this class
    */
   static const std::string s_name_id;
};

}
}

#endif
