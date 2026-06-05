find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
find_package(pybind11 CONFIG REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS C)

set(CORE_TARGET "_darkpropy")

FILE(GLOB TARGET_SRCS src/*.[hc]pp)

python_add_library(${CORE_TARGET} MODULE ${TARGET_SRCS} WITH_SOABI)

target_include_directories(
  ${CORE_TARGET} PRIVATE
  ${PROJECT_SOURCE_DIR}
  ${GENIE_INCLUDES}
)

target_link_libraries(
  ${CORE_TARGET} PRIVATE
  pybind11::headers
  GSL::gsl GSL::gslcblas
  hdf5::hdf5
  spdlog::spdlog
  ${GENIE_LIBS}
)
target_compile_definitions(${CORE_TARGET} PRIVATE VERSION_INFO=${PROJECT_VERSION})

if(SKBUILD)
  install(TARGETS ${CORE_TARGET} LIBRARY DESTINATION darkprop)
endif()
