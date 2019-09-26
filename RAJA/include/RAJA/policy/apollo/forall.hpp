/*!
 ******************************************************************************
 *
 * \file
 *
 * \brief   Header file containing RAJA index set and segment iteration
 *          template methods for Apollo-guided execution.
 *
 *          These methods should work on any platform.
 *
 ******************************************************************************
 */

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Copyright (c) 2016-18, Lawrence Livermore National Security, LLC.
//
// Produced at the Lawrence Livermore National Laboratory
//
// LLNL-CODE-689114
//
// All rights reserved.
//
// This file is part of RAJA.
//
// For details about use and distribution, please read RAJA/LICENSE.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

#ifndef RAJA_forall_apollo_HPP
#define RAJA_forall_apollo_HPP

#include <string>
#include <sstream>
#include <functional>
#include <unordered_set>

#include <omp.h>

#include "RAJA/config.hpp"
#include "RAJA/util/types.hpp"
#include "RAJA/policy/apollo/policy.hpp"
#include "RAJA/internal/fault_tolerance.hpp"

#include "apollo/Apollo.h"
#include "apollo/Region.h"

#include "CallpathRuntime.h"

namespace RAJA
{
namespace policy
{
namespace apollo
{

//template <typename T>
//int getTotalSize(T iterable) { return 0; }
//template <>
//int getTotalSize(RAJA::IndexSet idxSet)
//{
//    return idxSet.getLength();
//}
//
//template <typename T>
//int getSegmentCount(T iterable) { return 0; }
//template <>
//int getSegmentCount(RAJA::IndexSet idxSet)
//{
//    return idxSet.getNumSegments();
//}





//
//////////////////////////////////////////////////////////////////////
//
// The following function template switches between various RAJA
// execution policies based on feedback from the Apollo system.
//
//////////////////////////////////////////////////////////////////////
//

#ifndef RAJA_ENABLE_OPENMP
#error *** RAJA_ENABLE_OPENMP isn't defined! \
    This build of RAJA requires OpenMP to be enabled! ***
#endif

using apolloPolicySeq      = RAJA::seq_exec;
using apolloPolicySIMD     = RAJA::simd_exec;
using apolloPolicyLoopExec = RAJA::loop_exec;
using apolloPolicyOpenMP   = RAJA::omp_parallel_for_exec;

#if(0)
#define APOLLO_OMP_EXEC(__threads, __sched, __chunksize, __body) \
{                                                           \
    if (new_policy != old_policy) {                         \
        omp_set_num_threads(__threads);                     \
        omp_set_schedule(__sched, __chunksize);             \
        apollo->setFeature("num_threads", (double) __threads); \
    }                                                       \
    __body(apolloPolicyOpenMP{});                           \
};
#else
#define APOLLO_OMP_EXEC(__threads, __sched, __chunksize, __body) \
{                                                           \
    __body(apolloPolicyOpenMP{});                           \
};
#endif

template <typename BODY>
RAJA_INLINE void apolloPolicySwitcher(int new_policy, int old_policy, int tc[], BODY body) {
    static Apollo *apollo = Apollo::instance();
    switch(new_policy) {
        case   0: // The 0th policy is always a "safe" choice in Apollo as a
                  // default, or fail-safe when models are broken or partial..
                  // In the case of this OpenMP exploration template, the
                  // 0'th policy uses whatever was already set by the previous
                  // Apollo::Region's model, or the system defaults, if it is
                  // the first loop to get executed.
                  body(apolloPolicyOpenMP{});
                  break;
        case   1: // The 1st policy is a Sequential option, which will come into
                  // play for iterations when the number of elements a loop is
                  // operating over is low enough that the overhead of distrubuting
                  // the tasks to OpenMP is not worth paying. Learning will disrupt
                  // the performance of the application more, when this option is
                  // available, but the learned model will be able to make
                  // more significant performance improvements for applications
                  // with ocassional sparse inputs to loops.
                  body(apolloPolicySeq{});
                  break;
        case   2: APOLLO_OMP_EXEC( tc[0],     omp_sched_auto, -1, body); break;
        case   3: APOLLO_OMP_EXEC( tc[1],     omp_sched_auto, -1, body); break;
        case   4: APOLLO_OMP_EXEC( tc[2],     omp_sched_auto, -1, body); break;
        case   5: APOLLO_OMP_EXEC( tc[3],     omp_sched_auto, -1, body); break;
        case   6: APOLLO_OMP_EXEC( tc[4],     omp_sched_auto, -1, body); break;
        case   7: APOLLO_OMP_EXEC( tc[5],     omp_sched_auto, -1, body); break;
        case   8: APOLLO_OMP_EXEC( tc[0],   omp_sched_static, -1, body); break;
        case   9: APOLLO_OMP_EXEC( tc[1],   omp_sched_static, -1, body); break;
        case  10: APOLLO_OMP_EXEC( tc[2],   omp_sched_static, -1, body); break;
        case  11: APOLLO_OMP_EXEC( tc[3],   omp_sched_static, -1, body); break;
        case  12: APOLLO_OMP_EXEC( tc[4],   omp_sched_static, -1, body); break;
        case  13: APOLLO_OMP_EXEC( tc[5],   omp_sched_static, -1, body); break;
        case  14: APOLLO_OMP_EXEC( tc[0],  omp_sched_dynamic, -1, body); break;
        case  15: APOLLO_OMP_EXEC( tc[1],  omp_sched_dynamic, -1, body); break;
        case  16: APOLLO_OMP_EXEC( tc[2],  omp_sched_dynamic, -1, body); break;
        case  17: APOLLO_OMP_EXEC( tc[3],  omp_sched_dynamic, -1, body); break;
        case  18: APOLLO_OMP_EXEC( tc[4],  omp_sched_dynamic, -1, body); break;
        case  19: APOLLO_OMP_EXEC( tc[5],  omp_sched_dynamic, -1, body); break;
        case  20: APOLLO_OMP_EXEC( tc[0],   omp_sched_guided, -1, body); break;
        case  21: APOLLO_OMP_EXEC( tc[1],   omp_sched_guided, -1, body); break;
        case  22: APOLLO_OMP_EXEC( tc[2],   omp_sched_guided, -1, body); break;
        case  23: APOLLO_OMP_EXEC( tc[3],   omp_sched_guided, -1, body); break;
        case  24: APOLLO_OMP_EXEC( tc[4],   omp_sched_guided, -1, body); break;
        case  25: APOLLO_OMP_EXEC( tc[5],   omp_sched_guided, -1, body); break;
    }
    return;
}

const int POLICY_COUNT = 26;

//
// ///////////////////// Classic template ///////////////////////
//
// using apolloPolicySeq      = RAJA::seq_exec;
// using apolloPolicySIMD     = RAJA::simd_exec;
// using apolloPolicyLoopExec = RAJA::loop_exec;
// #if defined(RAJA_ENABLE_OPENMP)
//   using apolloPolicyOpenMP   = RAJA::omp_parallel_for_exec;
//   #if defined(RAJA_ENABLE_CUDA)
//     using apolloPolicyCUDA     = RAJA::cuda_exec<CUDA_BLOCK_SIZE>;
//     const int POLICY_COUNT = 5;
//   #else
//     const int POLICY_COUNT = 4;
//   #endif
// #else
//   #if defined(RAJA_ENABLE_CUDA)
//     using apolloPolicyCUDA     = RAJA::cuda_exec<CUDA_BLOCK_SIZE>;
//     const int POLICY_COUNT = 4;
//   #else
//     const int POLICY_COUNT = 3;
//   #endif
// #endif
//
// template <typename BODY>
// RAJA_INLINE void apolloPolicySwitcher(int choice, BODY body) {
//     switch (choice) {
//     case 1: body(apolloPolicySeq{});      break;
//     case 2: body(apolloPolicySIMD{});     break;
//     case 3: body(apolloPolicyLoopExec{}); break;
// #if defined(RAJA_ENABLE_OPENMP)
//     case 4: body(apolloPolicyOpenMP{});   break;
//     #if defined(RAJA_ENABLE_CUDA)
//     case 5: body(apolloPolicyCUDA{});     break;
//     #endif
// #else
//     #if defined(RAJA_ENABLE_CUDA)
//     case 4: body(apolloPolicyCUDA{});     break;
//     #endif
// #endif
//     case 0:
//     default: body(apolloPolicySeq{});     break;
//     }
// }

template <typename Iterable, typename Func>
RAJA_INLINE void forall_impl(const apollo_exec &, Iterable &&iter, Func &&body)
{
    static Apollo         *apollo             = Apollo::instance();
    static Apollo::Region *apolloRegion       = nullptr;
    static int             apollo_exec_count  = 0;
    static int             new_policy         = 0;
    static int             old_policy         = 0;
    static int             th_count_opts[6]   = {2, 2, 2, 2, 2, 2};
    if (apolloRegion == nullptr) {
        // Set up this Apollo::Region for the first time:       (Runs only once)
        std::stringstream ss_location;
        ss_location << apollo->callpath.doStackwalk().get(1);

        apolloRegion = new Apollo::Region(
            Apollo::instance(),
            ss_location.str().c_str(),
            RAJA::policy::apollo::POLICY_COUNT);
        // Set the range of thread counts we want to make available for
        // bootstrapping and use by this Apollo::Region.
        th_count_opts[0] = 2;
        th_count_opts[1] = std::min(4,  std::max(2, (int)(apollo->numThreadsPerProcCap * 0.25)));
        th_count_opts[2] = std::min(8,  std::max(2, (int)(apollo->numThreadsPerProcCap * 0.50)));
        th_count_opts[3] = std::min(16, std::max(2, (int)(apollo->numThreadsPerProcCap * 0.75)));
        th_count_opts[4] = std::min(32, std::max(2, apollo->numThreadsPerProcCap));
        th_count_opts[5] = std::max(2, apollo->numThreadsPerProcCap);
    }

    // Count the number of elements.
    double num_elements = 0.0;
    num_elements = (double) std::distance(std::begin(iter), std::end(iter));

    //
    // TODO: Also allow (elsewhere, just noting) the flushAllRegions event
    //       to be set up in different modes, where it step % rank sends, or
    //       only sends every N steps, etc.  I should have a JSON wrapper
    //       around the current JSON package so I can process it as either a
    //       Directives or ModelPackage
    //
    // TODO: Allow messages to be sent to specific ranks only
    // TODO: New magic word: __ANY_REGION__ --> __NEW_REGION__ that will only
    //       assign a policy to regions that have not been directly targeted
    //       before.
    //
    // TODO: Keep models around, when a loop gets encountered for the first time,
    //       make sure it has a chance to check out settings for itself from
    //       the prior package. (Applies when an application starts up using
    //       a prior learned model. This may already be working, just verify.)
    //
    //
    // If there isn't much to do this time...
    //if (num_elements < apolloRegion->minimum_elements_to_evaluate_model) {
    //    // Don't bother with any timing/reporting about this loop invocation.
    //    // Run the loop body with whatever OpenMP settings are in place from
    //    // this or some other previous Apollo model's evaluation:
    //    apolloPolicySwitcher(1, th_count_opts, [=] (auto pol) mutable {
    //        forall_impl(pol, iter, body); });
    //} else {
        // ...otherwise, we DO have enough work to bother checking out our
        // model and selecting an appropriate policy. Start timing, note
        // features, and pick an optimal solution:
        apolloRegion->begin(apollo_exec_count++);
        apollo->setFeature("num_elements", num_elements);

        old_policy = new_policy;
        new_policy = apolloRegion->getPolicyIndex();

        apolloPolicySwitcher(new_policy, old_policy, th_count_opts, [=] (auto pol) mutable {
            forall_impl(pol, iter, body); });

        apolloRegion->end();
    //}

}

//////////
}  // closing brace for apollo namespace
}  // closing brace for policy namespace
}  // closing brace for RAJA namespace

#endif  // closing endif for header file include guard
