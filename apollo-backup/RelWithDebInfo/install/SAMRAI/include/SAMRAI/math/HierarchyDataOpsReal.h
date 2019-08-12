/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2018 Lawrence Livermore National Security, LLC
 * Description:   Interface to templated operations for real data on hierarchy.
 *
 ************************************************************************/

#ifndef included_math_HierarchyDataOpsReal
#define included_math_HierarchyDataOpsReal

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/PatchHierarchy.h"

#include <iostream>
#include <memory>

namespace SAMRAI {
namespace math {

/**
 * Class HierarchyDataOpsReal defines the interface to a collection
 * of operations that may be used to manipulate real numerical (float
 * or double) patch data components over multiple levels in an AMR hierarchy.
 * It serves as a base class for subclasses which implement the operations
 * for cell-centered, face-centered, or node-centered data types.
 * The patch hierarchy and set of levels within that hierarcy over which the
 * operations will be performed are set in the constructor of the subclass.
 * However, these data members may be changed at any time via the virtual
 * access functions setPatchHierarchy() and resetLevels() below.  The
 * operations include basic arithmetic, norms and ordering, and assorted
 * miscellaneous operations.
 *
 * Note that, when it makes sense, an operation accept a boolean argument
 * which indicates whether the operation should be performed on all of the
 * data or just those data elements corresponding to the patch interiors.
 * If no boolean argument is provided, the default behavior is to treat only
 * the patch interiors.  Also, interfaces for similar sets of operations
 * for complex and integer hierarchy data are defined in the classes
 * HierarchyDataOpsComplex and HierarchyDataOpsInteger,
 * respectively.
 */

template<class TYPE>
class HierarchyDataOpsReal
{
public:
   /**
    * The constructor for the HierarchyDataOpsReal class.
    */
   HierarchyDataOpsReal();

   /**
    * Virtual destructor for the HierarchyDataOpsReal class.
    */
   virtual ~HierarchyDataOpsReal<TYPE>();

   /**
    * Reset patch hierarchy over which operations occur.
    */
   virtual
   void
   setPatchHierarchy(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy) = 0;

   /**
    * Reset range of patch levels over which operations occur.
    */
   virtual void
   resetLevels(
      const int coarsest_level,
      const int finest_level) = 0;

   /**
    * Return const pointer to patch hierarchy associated with operations.
    */
   virtual
   const std::shared_ptr<hier::PatchHierarchy>
   getPatchHierarchy() const = 0;

   /**
    * Copy source data to destination data.
    */
   virtual void
   copyData(
      const int dst_id,
      const int src_id,
      const bool interior_only = true) const = 0;

   /**
    * Swap data pointers (i.e., storage) between two data components.
    */
   virtual void
   swapData(
      const int data1_id,
      const int data2_id) const = 0;

   /**
    * Print data over multiple levels to specified output stream.
    */
   virtual void
   printData(
      const int data_id,
      std::ostream& s,
      const bool interior_only = true) const = 0;

   /**
    * Set data component to given scalar.
    */
   virtual void
   setToScalar(
      const int data_id,
      const TYPE& alpha,
      const bool interior_only = true) const = 0;

   /**
    * Set destination to source multiplied by given scalar, pointwise.
    */
   virtual void
   scale(
      const int dst_id,
      const TYPE& alpha,
      const int src_id,
      const bool interior_only = true) const = 0;

   /**
    * Add scalar to each entry in source data and set destination to result.
    */
   virtual void
   addScalar(
      const int dst_id,
      const int src_id,
      const TYPE& alpha,
      const bool interior_only = true) const = 0;

   /**
    * Set destination to sum of two source components, pointwise.
    */
   virtual void
   add(
      const int dst_id,
      const int src1_id,
      const int src2_id,
      const bool interior_only = true) const = 0;

   /**
    * Subtract second source component from first source component pointwise
    * and set destination data component to result.
    */
   virtual void
   subtract(
      const int dst_id,
      const int src1_id,
      const int src2_id,
      const bool interior_only = true) const = 0;

   /**
    * Set destination component to product of two source components, pointwise.
    */
   virtual void
   multiply(
      const int dst_id,
      const int src1_id,
      const int src2_id,
      const bool interior_only = true) const = 0;

   /**
    * Divide first data component by second source component pointwise
    * and set destination data component to result.
    */
   virtual void
   divide(
      const int dst_id,
      const int src1_id,
      const int src2_id,
      const bool interior_only = true) const = 0;

   /**
    * Set each entry of destination component to reciprocal of corresponding
    * source data component entry.
    */
   virtual void
   reciprocal(
      const int dst_id,
      const int src_id,
      const bool interior_only = true) const = 0;

   /**
    * Set \f$d = \alpha s_1 + \beta s_2\f$, where \f$d\f$ is the destination patch
    * data component and \f$s_1, s_2\f$ are the first and second source components,
    * respectively.  Here \f$\alpha, \beta\f$ are scalar values.
    */
   virtual void
   linearSum(
      const int dst_id,
      const TYPE& alpha,
      const int src1_id,
      const TYPE& beta,
      const int src2_id,
      const bool interior_only = true) const = 0;

   /**
    * Set \f$d = \alpha s_1 + s_2\f$, where \f$d\f$ is the destination patch data
    * component and \f$s_1, s_2\f$ are the first and second source components,
    * respectively.  Here \f$\alpha\f$ is a scalar.
    */
   virtual void
   axpy(
      const int dst_id,
      const TYPE& alpha,
      const int src1_id,
      const int src2_id,
      const bool interior_only = true) const = 0;

   /**
    * Set \f$d = \alpha s_1 - s_2\f$, where \f$d\f$ is the destination patch data
    * component and \f$s_1, s_2\f$ are the first and second source components,
    * respectively.  Here \f$\alpha\f$ is a scalar.
    */
   virtual void
   axmy(
      const int dst_id,
      const TYPE& alpha,
      const int src1_id,
      const int src2_id,
      const bool interior_only = true) const = 0;

   /**
    * Set destination data to absolute value of source data, pointwise.
    */
   virtual void
   abs(
      const int dst_id,
      const int src_id,
      const bool interior_only = true) const = 0;

   /**
    * Return minimum data value over all patches in the collection of levels.
    */
   virtual TYPE
   min(
      const int data_id,
      const bool interior_only = true) const = 0;

   /**
    * Return maximum data value over all patches in the collection of levels.
    */
   virtual TYPE
   max(
      const int data_id,
      const bool interior_only = true) const = 0;

   /**
    * Set data entries to random values.  See the operations in the
    * array data operation classes for each data type for details on the
    * generation of the random values.
    */
   virtual void
   setRandomValues(
      const int data_id,
      const TYPE& width,
      const TYPE& low,
      const bool interior_only = true) const = 0;

   /**
    * Return the total number of data values for the component on the set
    * of hierarchy levels.  If the boolean argument is true, the number of
    * elements will be summed over patch interiors.  If the boolean argument
    * is false, all elements will be counted (including ghost values)
    * over all patches.
    */
   virtual size_t
   numberOfEntries(
      const int data_id,
      const bool interior_only = true) const = 0;

   /**
    * Return sum of the control volumes associated with the data component.
    * Note that if the control volumes are set properly, this is equivalent to
    * integrating a data component containing all ones over the collection of
    * hierarchy levels.
    */
   virtual double
   sumControlVolumes(
      const int data_id,
      const int vol_id) const = 0;

   /**
    * Return discrete \f$L_1\f$-norm of the data using the control volume to
    * weight the contribution of each data entry to the sum.  That is, the
    * return value is the sum \f$\sum_i ( \| data_i \| cvol_i )\f$.  If the
    * control volume is not defined (vol_id < 0), the return value is
    * \f$\sum_i ( \| data_i \| )\f$.
    * If local_only is true, the global reduction is not
    * performed (thus each process will get only local results).
    */
   virtual double
   L1Norm(
      const int data_id,
      const int vol_id = -1,
      bool local_only = false) const = 0;

   /**
    * Return discrete \f$L_2\f$-norm of the data using the control volume to
    * weight the contribution of each data entry to the sum.  That is, the
    * return value is the sum \f$\sqrt{ \sum_i ( (data_i)^2 cvol_i ) }\f$.
    * If the control volume is  not defined (vol_id < 0), the return value
    * is \f$\sqrt{ \sum_i ( (data_i)^2 cvol_i ) }\f$.
    * If local_only is true, the global reduction is not
    * performed (thus each process will get only local results).
    */
   virtual double
   L2Norm(
      const int data_id,
      const int vol_id = -1,
      bool local_only = false) const = 0;

   /**
    * Return discrete weighted \f$L_2\f$-norm of the data using the control
    * volume to weight the contribution of the data and weight entries to
    * the sum.  That is, the return value is the sum \f$\sqrt{ \sum_i (
    * (data_i * weight_i)^2 cvol_i ) }\f$.  If the control volume is not defined
    * (vol_id < 0), the return value is
    * \f$\sqrt{ \sum_i ( (data_i * weight_i)^2 ) }\f$.
    */
   virtual double
   weightedL2Norm(
      const int data_id,
      const int weight_id,
      const int vol_id = -1) const = 0;

   /**
    * Return discrete root mean squared norm of the data.  If the control
    * volume is defined, the return value is the \f$L_2\f$-norm divided by
    * the square root of the sum of the control volumes.  Otherwise, the
    * return value is the \f$L_2\f$-norm divided by the square root of the
    * number of data entries.
    */
   virtual double
   RMSNorm(
      const int data_id,
      const int vol_id = -1) const = 0;

   /**
    * Return discrete weighted root mean squared norm of the data.  If the
    * control volume is define, the return value is the weighted \f$L_2\f$-norm
    * divided by the square root of the sum of the control volumes.  Otherwise,
    * the return value is the weighted \f$L_2\f$-norm divided by the square root
    * of the number of data entries.
    */
   virtual double
   weightedRMSNorm(
      const int data_id,
      const int weight_id,
      const int vol_id = -1) const = 0;

   /**
    * Return the \f$\max\f$-norm of the data using the control volume to weight
    * the contribution of each data entry to the maximum.  That is, the return
    * value is \f$\max_i ( \| data_i \| )\f$, where the max is over the data
    * elements where \f$cvol_i > 0\f$.  If the control volume is undefined
    * (vol_id < 0), it is ignored during the computation of the maximum.
    * If local_only is true, the global reduction is not
    * performed (thus each process will get only local results).
    */
   virtual double
   maxNorm(
      const int data_id,
      const int vol_id = -1,
      bool local_only = false) const = 0;

   /**
    * Return the dot product of the two data arrays using the control volume
    * to weight the contribution of each product to the sum.  That is, the
    * return value is the sum \f$\sum_i ( data1_i * data2_i * cvol_i )\f$.
    * If the control volume is undefined (vol_id < 0), it is ignored during
    * the summation.  If local_only is true, the global reduction is not
    * performed (thus each process will get only local results).
    */
   virtual TYPE
   dot(
      const int data1_id,
      const int data2_id,
      const int vol_id = -1,
      bool local_only = false) const = 0;

   /**
    * Return 1 if \f$\|data2_i\| > 0\f$ and \f$data1_i * data2_i \leq 0\f$, for
    * any \f$i\f$ in the set of patch data indices, where \f$cvol_i > 0\f$.  Otherwise,
    * return 0.  If the control volume is undefined (vol_id < 0), all values
    * on the patch interiors are considered.
    */
   virtual int
   computeConstrProdPos(
      const int data1_id,
      const int data2_id,
      const int vol_id = -1) const = 0;

   /**
    * Wherever \f$cvol_i > 0\f$ in the set of patch data indices, set \f$dst_i = 1\f$
    * if \f$\|src_i\| > \alpha\f$, and \f$dst_i = 0\f$ otherwise.  If the control
    * volume is undefined (vol_id < 0), all values on the patch interiors
    * are considered.
    */
   virtual void
   compareToScalar(
      const int dst_id,
      const int src_id,
      const TYPE& alpha,
      const int vol_id = -1) const = 0;

   /**
    * Wherever \f$cvol_i > 0\f$ in the set of patch data indices, set
    * \f$dst_i = 1/src_i\f$ if \f$src_i \neq 0\f$, and \f$dst_i = 0\f$ otherwise.
    * If \f$dst_i = 0\f$ anywhere, 0 is the return value.  Otherwise 1 is returned.
    * If the control volume is undefined (vol_id < 0), all values on the
    * patch interiors are considered.
    */
   virtual int
   testReciprocal(
      const int dst_id,
      const int src_id,
      const int vol_id = -1) const = 0;

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
    */
   virtual TYPE
   maxPointwiseDivide(
      const int numer_id,
      const int denom_id,
      bool local_only = false) const = 0;

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
    */
   virtual TYPE
   minPointwiseDivide(
      const int numer_id,
      const int denom_id,
      bool local_only = false) const = 0;

private:
   // The following are not implemented
   HierarchyDataOpsReal(
      const HierarchyDataOpsReal&);
   HierarchyDataOpsReal&
   operator = (
      const HierarchyDataOpsReal&);

};

}
}

#include "SAMRAI/math/HierarchyDataOpsReal.C"

#endif
