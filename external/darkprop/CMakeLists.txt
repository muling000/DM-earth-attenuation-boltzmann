cmake_minimum_required(VERSION 3.20)

project(
  DarkProp
  VERSION 0.3.0
  LANGUAGES C CXX
  DESCRIPTION "Monte Carlo simulation code for the propagation of dark matter particles in matter."
)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG")
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
message(">> CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" "${CMAKE_MODULE_PATH}")

set(DARKPROP_HEADER_DIR ${PROJECT_SOURCE_DIR}/darkprop)

if(NOT DEFINED INSTALL_HEADERS)
  set(INSTALL_HEADERS ON)
endif()

if(NOT DEFINED BUILD_APP)
  set(BUILD_APP ON)
endif()

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

find_package(GSL REQUIRED)
find_package(spdlog REQUIRED)

if(GENIE)
  message(">> GENIE support enabled for ${PROJECT_NAME}")
  find_package(OpenSSL COMPONENTS Crypto)  # for cluster compatible
  find_package(ROOT REQUIRED COMPONENTS Geom EGPythia6 MathMore)
  find_package(LibXml2 REQUIRED)
  find_package(Log4cpp REQUIRED)
  find_package(LHAPDF REQUIRED)
  find_package(GENIE REQUIRED)
  add_definitions(-D__WITH_GENIE__)
  set(GENIE_INCLUDES
      ${ROOT_INCLUDE_DIRS}
      ${Log4cpp_INCLUDE_DIRS}
      ${GENIE_INCLUDE_DIRS}
  )
  set(GENIE_LIBS
      LibXml2::LibXml2
      LHAPDF::LHAPDF
      OpenSSL::Crypto
      ${ROOT_LIBRARIES}
      ${Log4cpp_LIBRARIES}
      ${GENIE_LIBRARIES}
  )
endif()

if(BUILD_APP)
  add_subdirectory(app)
endif()

if(Julia)
  message(">> Julia binding enabled for ${PROJECT_NAME}")
  add_subdirectory(julia)
endif()

if(Python OR SKBUILD)
  message(">> Python binding enabled for ${PROJECT_NAME}")
  add_subdirectory(python)
endif()

if(BUILD_TESTING)
  enable_testing()
  add_subdirectory(tests)
  message(">> Test enabled for ${PROJECT_NAME}")
else()
  message(">> Test disabled for ${PROJECT_NAME}")
endif()

if(BUILD_DOCS)
  message(">> Build documentation for ${PROJECT_NAME}")
  add_subdirectory(docs)
endif()

# install the darkprop header only library
if(INSTALL_HEADERS)
  install(DIRECTORY ${DARKPROP_HEADER_DIR} TYPE INCLUDE)
endif()

if(NOT TARGET uninstall)
  configure_file(
    "${PROJECT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY
  )

  add_custom_target(
    uninstall COMMAND
    ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake
  )
endif()
