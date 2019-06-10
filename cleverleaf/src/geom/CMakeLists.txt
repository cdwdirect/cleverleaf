set(geom_sources
  CartesianCellDoubleConservativeLinearRefine.C
  CartesianCellDoubleMassWeightedAverage.C
  CartesianCellDoubleVolumeWeightedAverage.C
  CartesianCellIntConstantCoarsen.C
  CartesianNodeDoubleLinearRefine.C
  CartesianSideDoubleFirstOrderRefine.C)

string(REGEX REPLACE "\\.C(;|$)" ".h\\1" geom_headers "${geom_sources}")

set(geom_depends ${SAMRAI_LIBRARIES} cleverleaf_pdat)

if (ENABLE_CUDA)
  set(geom_depends ${geom_depends} cuda)
endif ()

blt_add_library(
  NAME cleverleaf_geom
  HEADERS ${geom_headers}
  SOURCES ${geom_sources}
  DEPENDS_ON ${geom_depends})
