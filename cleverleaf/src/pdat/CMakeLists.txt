set(pdat_sources
  detail.C
  NodeInjection.C)

string(REGEX REPLACE "\\.C(;|$)" ".h\\1" pdat_headers "${pdat_sources}")

set(pdat_depends ${SAMRAI_LIBRARIES})

if (ENABLE_CUDA)
  set(pdat_depends ${pdat_depends} cuda)
endif ()

blt_add_library(
  NAME cleverleaf_pdat
  HEADERS ${pdat_headers}
  SOURCES ${pdat_sources}
  DEPENDS_ON ${pdat_depends})
