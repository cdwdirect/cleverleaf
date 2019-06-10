/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2018 Lawrence Livermore National Security, LLC
 * Description:   Factory class for creating outerface data objects
 *
 ************************************************************************/

#ifndef included_pdat_OuterfaceDataFactory
#define included_pdat_OuterfaceDataFactory

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxGeometry.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/PatchDataFactory.h"
#include "SAMRAI/tbox/Complex.h"

#include <memory>

namespace SAMRAI {
namespace pdat {

/**
 * Class OuterfaceDataFactory is a factory class used to allocate new
 * instances of OuterfaceData objects.  It is a subclass of the patch
 * data factory class and outerface data is a subclass of patch data.  Both
 * the factory and data classes are templated on the type of the contained
 * object (e.g., double or int).
 *
 * @see OuterfaceData
 * @see PatchDataFactory
 */

template<class TYPE>
class OuterfaceDataFactory:public hier::PatchDataFactory
{
public:
   /**
    * The default constructor for the outerface data factory class.
    * The depth (number of components) gives the default for all of
    * the outerface data objects created with this factory.
    *
    * @pre depth > 0
    */
   OuterfaceDataFactory(
      const tbox::Dimension& dim,
      int depth);

   /**
    * Virtual destructor for the outerface data factory class.
    */
   virtual ~OuterfaceDataFactory<TYPE>();

   /**
    * @brief Abstract virtual function to clone a patch data factory.
    *
    * This will return a new instantiation of the abstract factory
    * with the same properties.  The properties of the cloned factory
    * can then be changed without modifying the original.
    *
    * @param ghosts default ghost cell width for concrete classes created from
    * the factory.
    *
    * @pre getDim() == ghosts.getDim()
    */
   virtual std::shared_ptr<hier::PatchDataFactory>
   cloneFactory(
      const hier::IntVector& ghosts);

   /**
    * Virtual factory function to allocate a concrete outerface data object.
    * The default information about the object (e.g., depth) is taken from
    * the factory.
    *
    * @pre getDim() == patch.getDim()
    */
   virtual std::shared_ptr<hier::PatchData>
   allocate(
      const hier::Patch& patch) const;

   /**
    * Allocate the box geometry object associated with the patch data.
    * This information will be used in the computation of intersections
    * and data dependencies between objects.
    *
    * @pre getDim() == box.getDim()
    */
   virtual std::shared_ptr<hier::BoxGeometry>
   getBoxGeometry(
      const hier::Box& box) const;

   /**
    * Get the depth (number of components).  This is the depth that
    * will be used in the instantiation of outerface data objects.
    */
   int
   getDepth() const;

   /**
    * Calculate the amount of memory needed to store the outerface data
    * object, including object data and dynamically allocated data.
    *
    * @pre getDim() == box.getDim()
    */
   virtual size_t
   getSizeOfMemory(
      const hier::Box& box) const;

   /**
    * Return a boolean true value indicating that fine data for the outerface
    * quantity will take precedence on coarse-fine interfaces.  See the
    * OuterfaceVariable<TYPE> class header file for more information.
    */
   bool
   fineBoundaryRepresentsVariable() const;

   /**
    * Return true since the outerface data index space extends beyond the
    * interior of patches.  That is, outerface data lives on patch borders.
    */
   bool
   dataLivesOnPatchBorder() const;

   /**
    * Return whether it is valid to copy this OuterfaceDataFactory to the
    * supplied destination patch data factory.  It will return true if
    * dst_pdf is FaceDataFactory or OuterfaceDataFactory, false otherwise.
    *
    * @pre getDim() == dst_pdf->getDim()
    */
   bool
   validCopyTo(
      const std::shared_ptr<hier::PatchDataFactory>& dst_pdf) const;

private:
   int d_depth;

   hier::IntVector d_no_ghosts;
};

}
}

#include "SAMRAI/pdat/OuterfaceDataFactory.C"

#endif
