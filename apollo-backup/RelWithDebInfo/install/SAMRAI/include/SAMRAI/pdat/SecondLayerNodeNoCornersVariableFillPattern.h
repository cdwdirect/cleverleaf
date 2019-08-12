/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2018 Lawrence Livermore National Security, LLC
 * Description:   Second layer node fill pattern class
 *
 ************************************************************************/

#ifndef included_pdat_SecondLayerNodeNoCornersVariableFillPattern
#define included_pdat_SecondLayerNodeNoCornersVariableFillPattern

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/BoxGeometry.h"
#include "SAMRAI/hier/BoxOverlap.h"
#include "SAMRAI/xfer/VariableFillPattern.h"

#include <memory>

namespace SAMRAI {
namespace pdat {

/*!
 * @brief Class SecondLayerNodeNoCornersVariableFillPattern is a concrete
 * implementation of the abstract base class VariableFillPattern.
 *
 * It is used to calculate overlaps according to a pattern which limits the
 * overlaps to the second layer of boundary nodes, that is the nodes lying
 * one cell width away from the patch boundary, but excluding nodes along
 * the corners and edges of the patch.
 */

class SecondLayerNodeNoCornersVariableFillPattern:
   public xfer::VariableFillPattern
{
public:
   /*!
    * @brief Constructor
    *
    * @param[in] dim     Dimension
    */
   explicit SecondLayerNodeNoCornersVariableFillPattern(
      const tbox::Dimension& dim);

   /*!
    * @brief Destructor
    */
   virtual ~SecondLayerNodeNoCornersVariableFillPattern();

   /*!
    * @brief Calculate overlaps between the destination and source geometries
    * according to the fill pattern.
    *
    * This will return the portion of the intersection of the node geometries
    * of the intersection of the node geometries that lies in the second
    * layer of boundary nodes, that is the nodes lying one cell width away
    * from the patch boundary, but excluding nodes along the corners and
    * edges of the patch.  The patch is identified by the argument
    * dst_patch_box.
    *
    * @param[in] dst_geometry    geometry object for destination box
    * @param[in] src_geometry    geometry object for source box
    * @param[in] dst_patch_box   box for the destination patch
    * @param[in] src_mask        the source mask, the box resulting from
    *                            transforming the source box
    * @param[in] fill_box        the box to be filled
    * @param[in] overwrite_interior  controls whether or not to include the
    *                                destination box interior in the overlap.
    *                                This has no effect in this implementation,
    *                                since all calculated overlaps will be
    *                                on the destination box exterior.
    * @param[in] transformation  the transformation from source to
    *                            destination index space.
    *
    * @return                std::shared_ptr to the calculated overlap object
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
      const hier::Transformation& transformation) const;

   /*!
    * @brief Computes a BoxOverlap object which defines the space to be filled
    * by a refinement operation.
    *
    * For this implementation, that space is the node-centered intersection
    * between fill_boxes, (computed by the RefineSchedule), data_box, which
    * specifies the extent of the destination data, and the boundary nodes
    * that define the stencil of this fill pattern.
    *
    * @param[in] fill_boxes  list representing the all of the space on a patch
    *                        or its ghost region that may be filled by a
    *                        refine operator (cell-centered represtentation)
    * @param[in] node_fill_boxes  node-centered represenation of fill_boxes
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
    * @brief Returns the stencil width of 1 in all directions.
    */
   const hier::IntVector&
   getStencilWidth();

   /*!
    * @brief Returns a string name identifier
    * "SECOND_LAYER_NODE_NO_CORNERS_FILL_PATTERN".
    */
   const std::string&
   getPatternName() const;

private:
   SecondLayerNodeNoCornersVariableFillPattern(
      const SecondLayerNodeNoCornersVariableFillPattern&);  // not implemented
   SecondLayerNodeNoCornersVariableFillPattern&
   operator = (
      const SecondLayerNodeNoCornersVariableFillPattern&);  // not implemented

   /*!
    * @brief Computes a BoxContainer defining the stencil around a given box
    * for this fill pattern
    *
    * The computed BoxContainer will consist of node-centered boxes representing
    * the nodes lying one cell width away from the box boundary, not including
    * those along the corners or edges of the box.
    *
    * @param[out] stencil_boxes   The computed stencil BoxContainer.
    * @param[in]  dst_box         Input box around which stencil is computed.
    *
    * @pre stencil_boxes.size() == 0
    */
   void
   computeStencilBoxes(
      hier::BoxContainer& stencil_boxes,
      const hier::Box& dst_box) const;

   /*!
    * @brief The dimension of this object.
    */
   const tbox::Dimension d_dim;

   /*!
    * @brief Static string holding name identifier for this class.
    */
   static const std::string s_name_id;

};

}
}
#endif
