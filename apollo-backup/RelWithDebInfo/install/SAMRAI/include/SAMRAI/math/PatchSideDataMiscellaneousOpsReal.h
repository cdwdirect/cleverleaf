/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2018 Lawrence Livermore National Security, LLC
 * Description:   Templated miscellaneous operations for real side-centered data.
 *
 ************************************************************************/

#ifndef included_math_PatchSideDataMiscellaneousOpsReal
#define included_math_PatchSideDataMiscellaneousOpsReal

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/math/ArrayDataMiscellaneousOpsReal.h"
#include "SAMRAI/hier/Box.h"

#include <memory>


namespace SAMRAI {
namespace math {

/**
 * Class PatchSideDataMiscellaneousOpsReal provides access to a
 * collection of operations that may be applied to numerical side-centered
 * patch data of type double and float.  The primary intent of this class is
 * to provide the interface to these operations for the class
 * PatchSideDataOpsReal which provides access to a more complete
 * set of operations that may be used to manipulate side-centered
 * patch data.  Each member function accepts a box argument
 * indicating the region of index space on which the operation should be
 * performed.  The operation will be performed on the intersection of this
 * box and those boxes corresponding to the patch data objects.  Also, each
 * operation allows an additional side-centered patch data object to be used
 * to represent a control volume that weights the contribution of each data
 * entry in the given norm calculation.  Note that the control volume patch
 * data must be of type double and have side-centered geometry (i.e., the
 * same as the data itself).  The use of control volumes is important when
 * these operations are used in vector kernels where the data resides over
 * multiple levels of spatial resolution in an AMR hierarchy.  If the control
 * volume is not given in the function call, it will be ignored in the
 * calculation.  Also, note that the depth of the control volume patch data
 * object must be either 1 or be equal to the depth of the other data objects.
 *
 * Since these operations are used only by the vector kernels for the KINSOL
 * and CVODE solver packages at this time, they are intended to be instantiated
 * for the standard built-in types double and float (since those solvers only
 * treat double and float data).  To extend this class to other data types or
 * to include other operations, the member functions must be specialized or the
 * new operations must be added.
 *
 * @see ArrayDataMiscellaneousOpsReal
 */

template<class TYPE>
class PatchSideDataMiscellaneousOpsReal
{
public:
   /**
    * Empty constructor and destructor.
    */
   PatchSideDataMiscellaneousOpsReal();

   virtual ~PatchSideDataMiscellaneousOpsReal<TYPE>();

   /**
    * Return 1 if \f$\|data2_i\| > 0\f$ and \f$data1_i * data2_i \leq 0\f$, for
    * any \f$i\f$ in the index region, where \f$cvol_i > 0\f$.  Otherwise return 0.
    * If the control volume is NULL, all values in the index set are used.
    *
    * @pre data1 && data2
    * @pre data1->getDirectionVector() == data2->getDirectionVector()
    * @pre !cvol || (data1->getDirectionVector() == hier::IntVector::min(data1->getDirectionVector(), cvol->getDirectionVector()))
    */
   int
   computeConstrProdPos(
      const std::shared_ptr<pdat::SideData<TYPE> >& data1,
      const std::shared_ptr<pdat::SideData<TYPE> >& data2,
      const hier::Box& box,
      const std::shared_ptr<pdat::SideData<double> >& cvol =
         std::shared_ptr<pdat::SideData<double> >()) const;

   /**
    * Wherever \f$cvol_i > 0\f$ in the index region, set \f$dst_i = 1\f$
    * if \f$\|src_i\| > \alpha\f$, and \f$dst_i = 0\f$ otherwise.  If the control
    * volume is NULL, all values in the index set are considered.
    *
    * @pre dst && src
    * @pre dst->getDirectionVector() == src->getDirectionVector()
    * @pre !cvol || (dst->getDirectionVector() == hier::IntVector::min(dst->getDirectionVector(), cvol->getDirectionVector()))
    */
   void
   compareToScalar(
      const std::shared_ptr<pdat::SideData<TYPE> >& dst,
      const std::shared_ptr<pdat::SideData<TYPE> >& src,
      const TYPE& alpha,
      const hier::Box& box,
      const std::shared_ptr<pdat::SideData<double> >& cvol =
         std::shared_ptr<pdat::SideData<double> >()) const;

   /**
    * Wherever \f$cvol_i > 0\f$ in the index region, set \f$dst_i = 1/src_i\f$ if
    * \f$src_i \neq 0\f$, and \f$dst_i = 0\f$ otherwise.  If \f$dst_i = 0\f$ anywhere,
    * 0 is the return value.  Otherwise 1 is returned.  If the control volume
    * all values in the index set are considered.
    *
    * @pre dst && src
    * @pre dst->getDirectionVector() == src->getDirectionVector()
    * @pre !cvol || (dst->getDirectionVector() == hier::IntVector::min(dst->getDirectionVector(), cvol->getDirectionVector()))
    */
   int
   testReciprocal(
      const std::shared_ptr<pdat::SideData<TYPE> >& dst,
      const std::shared_ptr<pdat::SideData<TYPE> >& src,
      const hier::Box& box,
      const std::shared_ptr<pdat::SideData<double> >& cvol =
         std::shared_ptr<pdat::SideData<double> >()) const;

   /*!
    * @brief Compute max of "conditional" quotients of two arrays.
    *
    * Return the maximum of pointwise "conditional" quotients of the numerator
    * and denominator.
    *
    * The "conditional" quotient is defined as |numerator/denominator|
    * where the denominator is nonzero.  Otherwise, it is defined as
    * |numerator|.
    *
    * @b Note: This method is currently intended to support the
    * PETSc-2.1.6 vector wrapper only.  Please do not use it!
    *
    * @pre numer && denom
    */
   TYPE
   maxPointwiseDivide(
      const std::shared_ptr<pdat::SideData<TYPE> >& numer,
      const std::shared_ptr<pdat::SideData<TYPE> >& denom,
      const hier::Box& box) const;

   /*!
    * @brief Compute min of quotients of two arrays.
    *
    * Return the minimum of pointwise quotients of the numerator
    * and denominator.
    *
    * The quotient is defined as (numerator/denominator)
    * where the denominator is nonzero.  When the denominator is zero, the
    * entry is skipped.  If the denominator is always zero, the value of
    * tbox::IEEE::getFLT_MAX() is returned (see @ref SAMRAI::tbox::IEEE).
    *
    * @b Note: This method is currently intended to support the
    * SUNDIALS vector wrapper only.  Please do not use it!
    *
    * @pre numer && denom
    */
   TYPE
   minPointwiseDivide(
      const std::shared_ptr<pdat::SideData<TYPE> >& numer,
      const std::shared_ptr<pdat::SideData<TYPE> >& denom,
      const hier::Box& box) const;

private:
   // The following are not implemented:
   PatchSideDataMiscellaneousOpsReal(
      const PatchSideDataMiscellaneousOpsReal&);
   PatchSideDataMiscellaneousOpsReal&
   operator = (
      const PatchSideDataMiscellaneousOpsReal&);

   ArrayDataMiscellaneousOpsReal<TYPE> d_array_ops;
};

}
}

#include "SAMRAI/math/PatchSideDataMiscellaneousOpsReal.C"

#endif
