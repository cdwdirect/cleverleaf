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

#define APOLLO_OMP_EXEC(__threads, __sched, __chunksize, __body) \
{                                                           \
    omp_set_num_threads(__threads);                         \
    omp_set_schedule(__sched, __chunksize);                 \
    __body(apolloPolicyOpenMP{});                           \
};


template <typename BODY>
RAJA_INLINE void apolloPolicySwitcher(int choice, BODY body) {
    switch(choice) {
        case   0: APOLLO_OMP_EXEC(   1,     omp_sched_auto, -1, body); break;
        case   1: APOLLO_OMP_EXEC(   2,     omp_sched_auto, -1, body); break;
        case   2: APOLLO_OMP_EXEC(   4,     omp_sched_auto, -1, body); break;
        case   3: APOLLO_OMP_EXEC(   6,     omp_sched_auto, -1, body); break;
        case   4: APOLLO_OMP_EXEC(   8,     omp_sched_auto, -1, body); break;
        case   5: APOLLO_OMP_EXEC(  12,     omp_sched_auto, -1, body); break;
        case   6: APOLLO_OMP_EXEC(  16,     omp_sched_auto, -1, body); break;
        case   7: APOLLO_OMP_EXEC(  24,     omp_sched_auto, -1, body); break;
        case   8: APOLLO_OMP_EXEC(  48,     omp_sched_auto, -1, body); break;
        case   9: APOLLO_OMP_EXEC(  64,     omp_sched_auto, -1, body); break;
        case  10: APOLLO_OMP_EXEC(   1,   omp_sched_static, -1, body); break;
        case  11: APOLLO_OMP_EXEC(   2,   omp_sched_static, -1, body); break;
        case  12: APOLLO_OMP_EXEC(   4,   omp_sched_static, -1, body); break;
        case  13: APOLLO_OMP_EXEC(   6,   omp_sched_static, -1, body); break;
        case  14: APOLLO_OMP_EXEC(   8,   omp_sched_static, -1, body); break;
        case  15: APOLLO_OMP_EXEC(  12,   omp_sched_static, -1, body); break;
        case  16: APOLLO_OMP_EXEC(  16,   omp_sched_static, -1, body); break;
        case  17: APOLLO_OMP_EXEC(  24,   omp_sched_static, -1, body); break;
        case  18: APOLLO_OMP_EXEC(  48,   omp_sched_static, -1, body); break;
        case  19: APOLLO_OMP_EXEC(  64,   omp_sched_static, -1, body); break;
        case  20: APOLLO_OMP_EXEC(   1,  omp_sched_dynamic, -1, body); break;
        case  21: APOLLO_OMP_EXEC(   2,  omp_sched_dynamic, -1, body); break;
        case  22: APOLLO_OMP_EXEC(   4,  omp_sched_dynamic, -1, body); break;
        case  23: APOLLO_OMP_EXEC(   6,  omp_sched_dynamic, -1, body); break;
        case  24: APOLLO_OMP_EXEC(   8,  omp_sched_dynamic, -1, body); break;
        case  25: APOLLO_OMP_EXEC(  12,  omp_sched_dynamic, -1, body); break;
        case  26: APOLLO_OMP_EXEC(  16,  omp_sched_dynamic, -1, body); break;
        case  27: APOLLO_OMP_EXEC(  24,  omp_sched_dynamic, -1, body); break;
        case  28: APOLLO_OMP_EXEC(  48,  omp_sched_dynamic, -1, body); break;
        case  29: APOLLO_OMP_EXEC(  64,  omp_sched_dynamic, -1, body); break;
        case  30: APOLLO_OMP_EXEC(   1,   omp_sched_guided, -1, body); break;
        case  31: APOLLO_OMP_EXEC(   2,   omp_sched_guided, -1, body); break;
        case  32: APOLLO_OMP_EXEC(   4,   omp_sched_guided, -1, body); break;
        case  33: APOLLO_OMP_EXEC(   6,   omp_sched_guided, -1, body); break;
        case  34: APOLLO_OMP_EXEC(   8,   omp_sched_guided, -1, body); break;
        case  35: APOLLO_OMP_EXEC(  12,   omp_sched_guided, -1, body); break;
        case  36: APOLLO_OMP_EXEC(  16,   omp_sched_guided, -1, body); break;
        case  37: APOLLO_OMP_EXEC(  24,   omp_sched_guided, -1, body); break;
        case  38: APOLLO_OMP_EXEC(  48,   omp_sched_guided, -1, body); break;
        case  39: APOLLO_OMP_EXEC(  64,   omp_sched_guided, -1, body); break;
        case  40: body(apolloPolicySeq{}); break;
    }
    return;
}

const int POLICY_COUNT = 40;

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
    static Apollo         *apollo            = Apollo::instance();
    static Apollo::Region *apolloRegion      = nullptr;
    static int             apolloExecCount   = 0;
    static int             policyIndex       = 0;
    if (apolloRegion == nullptr) {
        // ----------
        // NOTE: This section runs *once* the first time the
        //       region is encountered
        std::stringstream ss_location;
        ss_location << (const void *) &body;
        apolloRegion = new Apollo::Region(
            Apollo::instance(),
            ss_location.str().c_str(),
            RAJA::policy::apollo::POLICY_COUNT);
        // ----------
    }

    double num_elements = 0.0;
    num_elements = (double) std::distance(std::begin(iter), std::end(iter));

    if (num_elements < 10.0) {
         apolloPolicySwitcher(40 , [=] (auto pol) mutable {
             forall_impl(pol, iter, body); });
    } else {
        apolloRegion->begin(apolloExecCount++);
        apollo->setFeature("num_elements", num_elements);
        policyIndex = apolloRegion->getPolicyIndex();

        apolloPolicySwitcher(policyIndex , [=] (auto pol) mutable {
            forall_impl(pol, iter, body); });

        apolloRegion->end();
    }

}

//////////
}  // closing brace for apollo namespace
}  // closing brace for policy namespace
}  // closing brace for RAJA namespace

#endif  // closing endif for header file include guard
