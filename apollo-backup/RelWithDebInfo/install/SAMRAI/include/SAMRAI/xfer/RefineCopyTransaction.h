/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2018 Lawrence Livermore National Security, LLC
 * Description:   Communication transaction for data copies during data
 *                refining
 *
 ************************************************************************/

#ifndef included_xfer_RefineCopyTransaction
#define included_xfer_RefineCopyTransaction

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/Transaction.h"
#include "SAMRAI/hier/BaseGridGeometry.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/xfer/RefineClasses.h"

#include <iostream>

namespace SAMRAI {
namespace xfer {

/*!
 * @brief Class RefineCopyTransaction represents a single copy communication
 * transaction between two processors or a local data copy for refine
 * schedules.  Note that to there is an implicit hand-shaking between objects
 * of this class and the RefineSchedule object that constructs them.  Following
 * the refine schedule implementation, the source patch data index for a copy
 * transaction always refers to the source data and the destination patch data
 * index for a copy transaction is always the scratch data, all as defined in
 * the RefineClasses class.
 *
 * @see RefineSchedule
 * @see RefineClasses
 * @see tbox::Schedule
 * @see tbox::Transaction
 */

class RefineCopyTransaction:public tbox::Transaction
{
public:
   /*!
    * Construct a transaction with the specified source and destination
    * levels, patches, and patch data components found in the refine class
    * item with the given id owned by the calling refine schedule.  In general,
    * this constructor is called by a RefineSchedule object for each data
    * transaction (not involving time interpolation) that must occur.  This
    * transaction will be responsible for one of the following: (1) a local
    * data copy, (2) packing a message stream with source patch data, or (3)
    * unpacking destination patch data from a message stream.
    *
    * @param dst_level        std::shared_ptr to destination patch level.
    * @param src_level        std::shared_ptr to source patch level.
    * @param overlap          std::shared_ptr to overlap region between
    *                         patches.
    * @param dst_box          Destination Box in destination patch level.
    * @param src_box          Source Box in source patch level.
    * @param refine_data      Pointer to array of refine data items
    * @param item_id          Integer id of refine data item owned by refine
    *                         schedule.
    *
    * @pre dst_level
    * @pre src_level
    * @pre overlap
    * @pre dst_box.getLocalId() >= 0
    * @pre src_box.getLocalId() >= 0
    * @pre refine_data != 0
    * @pre refine_item_id >= 0
    * @pre (dst_level->getDim() == src_level->getDim()) &&
    *      (dst_level->getDim() == dst_box.getDim()) &&
    *      (dst_level->getDim() == src_box.getDim())
    */
   RefineCopyTransaction(
      const std::shared_ptr<hier::PatchLevel>& dst_level,
      const std::shared_ptr<hier::PatchLevel>& src_level,
      const std::shared_ptr<hier::BoxOverlap>& overlap,
      const hier::Box& dst_box,
      const hier::Box& src_box,
      const RefineClasses::Data ** refine_data,
      int item_id);

   /*!
    * The virtual destructor for the copy transaction releases all
    * memory associated with the transaction.
    */
   virtual ~RefineCopyTransaction();

   /*!
    * Return a boolean indicating whether this transaction can estimate
    * the size of an incoming message.  If this is false, then a different
    * communication protocol kicks in and the message size is transmitted
    * between boxes.
    */
   virtual bool
   canEstimateIncomingMessageSize();

   /*!
    * Return the integer buffer space (in bytes) needed for the incoming
    * message.  This routine is only called if the transaction can estimate the
    * size of the incoming message.  See canEstimateIncomingMessageSize().
    */
   virtual size_t
   computeIncomingMessageSize();

   /*!
    * Return the integer buffer space (in bytes) needed for the outgoing
    * message.
    */
   virtual size_t
   computeOutgoingMessageSize();

   /*!
    * Return the sending processor number for the communications transaction.
    */
   virtual int
   getSourceProcessor();

   /*!
    * Return the receiving processor number for the communications transaction.
    */
   virtual int
   getDestinationProcessor();

   /*!
    * Pack the transaction data into the message stream.
    */
   virtual void
   packStream(
      tbox::MessageStream& stream);

   /*!
    * Unpack the transaction data from the message stream.
    */
   virtual void
   unpackStream(
      tbox::MessageStream& stream);

   /*!
    * Perform the local data copy for the transaction.
    */
   virtual void
   copyLocalData();

   /*!
    * Print out transaction information.
    */
   virtual void
   printClassData(
      std::ostream& stream) const;

private:
   RefineCopyTransaction(
      const RefineCopyTransaction&);                    // not implemented
   RefineCopyTransaction&
   operator = (
      const RefineCopyTransaction&);                    // not implemented

   std::shared_ptr<hier::Patch> d_dst_patch;
   int d_dst_patch_rank;
   std::shared_ptr<hier::Patch> d_src_patch;
   int d_src_patch_rank;
   std::shared_ptr<hier::BoxOverlap> d_overlap;
   const RefineClasses::Data** d_refine_data;
   int d_item_id;
   size_t d_incoming_bytes;
   size_t d_outgoing_bytes;

};

}
}

#endif
