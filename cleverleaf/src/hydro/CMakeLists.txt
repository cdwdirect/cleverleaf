set(hydro_sources
  LagrangianEulerianIntegrator.C
  LagrangianEulerianLevelIntegrator.C
  LagrangianEulerianPatchStrategy.C
  Cleverleaf.C)

string(REGEX REPLACE "\\.C(;|$)" ".h\\1" hydro_headers "${hydro_sources}")

set(hydro_depends ${SAMRAI_LIBRARIES} cleverleaf_pdat cleverleaf_geom)

if (ENABLE_CUDA)
  set(hydro_depends ${hydro_depends} cuda)
endif ()

# TODO Only compile Cleverleaf.C with NVCC if needed, otherwise use host compiler to generate object files
blt_add_library(
  NAME cleverleaf_hydro
  HEADERS ${hydro_headers}
  SOURCES ${hydro_sources}
  DEPENDS_ON ${hydro_depends})
