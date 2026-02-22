cmake_minimum_required(VERSION 3.13)
project(capdDDEs LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(capdDDEs_BUILD_TESTS "Build capdDDEs tests" ON)
option(capdDDEs_BUILD_EXAMPLES "Build capdDDEs examples" ON)
option(capdDDEs_BUILD_UTILS "Build capdDDEs utility programs" ON)
option(
  capdDDEs_USE_SYSTEM_CAPD
  "Build with the bundled version of the CAPD library" ON
)
option(capdDDEs_INSTALL "Enable installation of capdDDEs library" ON)

list(
  APPEND CMAKE_MODULE_PATH
  "${CMAKE_CURRENT_SOURCE_DIR}/external/cmake/BoostTestDiscoverTests/modules"
)

include(CTest)

if(capdDDEs_USE_SYSTEM_CAPD)
  find_package(capd)
  if(NOT capd_FOUND)
    message(FATAL_ERROR
"System CAPD library not found!
To use the bundled version of the CAPD library instead, set \
capdDDEs_USE_SYSTEM_CAPD=OFF"
    )
  endif()
  message("-- Using system CAPD.")
else()
  set(CAPD_ENABLE_MULTIPRECISION ON)
  add_subdirectory(external/capd)
  add_library(capd::capd ALIAS capd)
  message("-- Using bundled CAPD.")
endif()

add_library(capdDDEs)
add_library(capdDDEs::capdDDEs ALIAS capdDDEs)

target_sources(capdDDEs
  PRIVATE
    "src/capd/ddes/DDECommon.cpp"
    "src/capd/ddeshelper/DDEHelperCommon.cpp"
    "src/capd/ddeshelper/DDEHelperDrawing.cpp"
    "src/capd/ddeshelper/DDESystemOSHelper.cpp"
    "src/capd/ddeshelper/ddeshelperlib.cpp"
)

target_include_directories(capdDDEs
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
)

target_link_libraries(capdDDEs PUBLIC capd::capd)

if(capdDDEs_BUILD_TESTS AND BUILD_TESTING)
  add_subdirectory(tests)
endif()

if(capd_DDEs_BUILD_EXAMPLES)
  add_subdirectory(programs/examples)
endif()

if(capdDDEs_BUILD_UTILS)
  add_subdirectory(programs/utils)
endif()

if(capdDDEs_INSTALL)
  include(CMakePackageConfigHelpers)
  include(GNUInstallDirs)

  install(
    TARGETS capdDDEs
    EXPORT capdddes-targets
  )

  install(
    DIRECTORY "include/"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  )

  install(
    EXPORT capdddes-targets
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
    NAMESPACE capdDDEs::
  )

  configure_package_config_file(
    "cmake/capdddes-config.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake/capdddes-config.cmake"
    INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
  )

  install(
    FILES "${CMAKE_CURRENT_BINARY_DIR}/cmake/capdddes-config.cmake"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
  )
endif()
