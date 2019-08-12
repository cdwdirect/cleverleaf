//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2018-2019, Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory
//
// Created by David Beckingsale, david@llnl.gov
// LLNL-CODE-747640
//
// All rights reserved.
//
// This file is part of Umpire.
//
// For details, see https://github.com/LLNL/Umpire
// Please also see the LICENSE file for MIT license.
//////////////////////////////////////////////////////////////////////////////
#ifndef UMPIRE_config_HPP
#define UMPIRE_config_HPP

/* #undef UMPIRE_ENABLE_CUDA */
/* #undef UMPIRE_ENABLE_NUMA */
/* #undef UMPIRE_ENABLE_SLIC */
/* #undef UMPIRE_ENABLE_LOGGING */
#define UMPIRE_ENABLE_ASSERTS
/* #undef UMPIRE_ENABLE_STATISTICS */
/* #undef UMPIRE_ENABLE_HCC */

/* #undef UMPIRE_ENABLE_DEVICE */
/* #undef UMPIRE_ENABLE_PINNED */
/* #undef UMPIRE_ENABLE_UM */

constexpr int UMPIRE_VERSION_MAJOR = 0;
constexpr int UMPIRE_VERSION_MINOR = 3;
constexpr int UMPIRE_VERSION_PATCH = 3;

#endif
