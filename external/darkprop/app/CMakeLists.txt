find_package(MPI REQUIRED COMPONENTS CXX)
find_package(HDF5 REQUIRED COMPONENTS C)
set(HDF5_PREFER_PARALLEL TRUE)

add_executable(darkprop-bin main.cpp)
set_target_properties(darkprop-bin PROPERTIES OUTPUT_NAME darkprop)

target_include_directories(
    darkprop-bin PRIVATE
    ${PROJECT_SOURCE_DIR}
    ${GENIE_INCLUDES}
)
target_link_libraries(
    darkprop-bin PUBLIC
    MPI::MPI_CXX
    GSL::gsl GSL::gslcblas
    hdf5::hdf5
    spdlog::spdlog
    ${GENIE_LIBS}
)

install(TARGETS darkprop-bin RUNTIME)
