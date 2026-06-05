set(CMAKE_MACOSX_RPATH 1)

execute_process(
  COMMAND julia -e "using CxxWrap; print(CxxWrap.prefix_path())"
  OUTPUT_VARIABLE CxxWrap_prefix
)

set(JlCxx_DIR "${CxxWrap_prefix}/lib/cmake/JlCxx")
message(">> JlCxx_DIR: ${JlCxx_DIR}")

find_package(JlCxx)
get_target_property(JlCxx_location JlCxx::cxxwrap_julia LOCATION)
get_filename_component(JlCxx_location ${JlCxx_location} DIRECTORY)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${JlCxx_location}")

add_library(darkpropjl SHARED darkpropjl.cpp)

target_include_directories(
  darkpropjl PRIVATE
  ${PROJECT_SOURCE_DIR}
  ${GENIE_INCLUDES}
)

target_link_libraries(
  darkpropjl PUBLIC
  JlCxx::cxxwrap_julia
  JlCxx::cxxwrap_julia_stl
  GSL::gsl GSL::gslcblas
  spdlog::spdlog
  ${GENIE_LIBS}
)

install(TARGETS darkpropjl LIBRARY)
