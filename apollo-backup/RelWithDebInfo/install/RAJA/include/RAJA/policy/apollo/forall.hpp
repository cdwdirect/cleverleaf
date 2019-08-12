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


using apolloPolicySeq      = RAJA::seq_exec;
using apolloPolicySIMD     = RAJA::simd_exec;
using apolloPolicyLoopExec = RAJA::loop_exec;
#if defined(RAJA_ENABLE_OPENMP)
  using apolloPolicyOpenMP   = RAJA::omp_parallel_for_exec;
  #if defined(RAJA_ENABLE_CUDA)
    using apolloPolicyCUDA     = RAJA::cuda_exec<CUDA_BLOCK_SIZE>;
    const int POLICY_COUNT = 5;
  #else
    const int POLICY_COUNT = 4;
  #endif
#else
  #if defined(RAJA_ENABLE_CUDA)
    using apolloPolicyCUDA     = RAJA::cuda_exec<CUDA_BLOCK_SIZE>;
    const int POLICY_COUNT = 4;
  #else
    const int POLICY_COUNT = 3;
  #endif
#endif



//
//////////////////////////////////////////////////////////////////////
//
// The following function template switches between various RAJA
// execution policies based on feedback from the Apollo system.
//
//////////////////////////////////////////////////////////////////////
//

template <typename BODY>
RAJA_INLINE void apolloPolicySwitcher(int choice, BODY body) {
    switch (choice) {
    case 1: body(apolloPolicySeq{});      break;
    case 2: body(apolloPolicySIMD{});     break;
    case 3: body(apolloPolicyLoopExec{}); break;
#if defined(RAJA_ENABLE_OPENMP)
    case 4: body(apolloPolicyOpenMP{});   break;
    #if defined(RAJA_ENABLE_CUDA)
    case 5: body(apolloPolicyCUDA{});     break;
    #endif
#else
    #if defined(RAJA_ENABLE_CUDA)
    case 4: body(apolloPolicyCUDA{});     break;
    #endif
#endif
    case 0: 
    default: body(apolloPolicySeq{});     break;
    }
}

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


    apolloRegion->begin(apolloExecCount++);
    apollo->setFeature("num_elements", (double) std::distance(std::begin(iter), std::end(iter)));
    policyIndex = apolloRegion->getPolicyIndex();

    apolloPolicySwitcher(policyIndex , [=] (auto pol) mutable {
        forall_impl(pol, iter, body); });
    
    apolloRegion->end();
    
}

//////////
}  // closing brace for apollo namespace
}  // closing brace for policy namespace
}  // closing brace for RAJA namespace

#endif  // closing endif for header file include guard
