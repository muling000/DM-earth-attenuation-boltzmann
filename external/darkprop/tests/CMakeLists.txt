find_package(GSL REQUIRED)

file(GLOB test_files *.cpp)
add_executable(tests ${test_files})
target_include_directories(
    tests PRIVATE
    ${PROJECT_SOURCE_DIR}
)
target_link_libraries(
    tests PUBLIC
    spdlog::spdlog
    GSL::gsl GSL::gslcblas
)

add_test(NAME doctest COMMAND tests)
