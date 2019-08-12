/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2018 Lawrence Livermore National Security, LLC
 * Description:   Class to record statistics during program execution.
 *
 ************************************************************************/

#ifndef included_tbox_NVTXUtilities
#define included_tbox_NVTXUtilities

#if defined(ENABLE_NVTX_REGIONS) && defined(ENABLE_CUDA)

#include "cuda_runtime.h"
#include "nvToolsExt.h"
#include "cuda_profiler_api.h"

const uint32_t colors[] = { 0x0000ff00, 0x000000ff, 0x00ffff00, 0x00ff00ff, 0x0000ffff, 0x00ff0000, 0x00ffffff };
const int num_colors = sizeof(colors)/sizeof(uint32_t);

#define RANGE_PUSH(name,cid) { \
   int color_id = cid; \
   color_id = color_id%num_colors;\
   nvtxEventAttributes_t eventAttrib = {0}; \
   eventAttrib.version = NVTX_VERSION; \
   eventAttrib.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE; \
   eventAttrib.colorType = NVTX_COLOR_ARGB; \
   eventAttrib.color = colors[color_id]; \
   eventAttrib.messageType = NVTX_MESSAGE_TYPE_ASCII; \
   eventAttrib.message.ascii = name; \
   cudaDeviceSynchronize();\
   nvtxRangePushEx(&eventAttrib); \
}

#define RANGE_POP { \
   cudaDeviceSynchronize(); \
   nvtxRangePop(); \
}

#else

#define RANGE_PUSH(name,cid)

#define RANGE_POP

#endif // ENABLE_NVTX_REGIONS

#endif // included_tbox_NVTXUtilities
