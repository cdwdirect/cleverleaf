include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(pdat)
add_subdirectory(geom)
add_subdirectory(hydro)

set(cleverleaf_depends ${SAMRAI_LIBRARIES} cleverleaf_pdat cleverleaf_geom cleverleaf_hydro)

if (ENABLE_APOLLO)
  set(cleverleaf_depends ${cleverleaf_depends} apollo)
  set(cleverleaf_depends ${cleverleaf_depends} callpath)
endif ()

if (ENABLE_OPENMP)
    set(cleverleaf_depends ${cleverleaf_depends} gomp)
endif()

if (ENABLE_CUDA)
  set(cleverleaf_depends ${cleverleaf_depends} cuda)
endif ()

blt_add_executable(
  NAME cleverleaf
  SOURCES main.C
  DEPENDS_ON ${cleverleaf_depends})

# TODO I think this is set by blt...
# set_target_properties(cleverleaf PROPERTIES LINKER_LANGUAGE CXX)

install(
  TARGETS cleverleaf
  RUNTIME DESTINATION bin)
