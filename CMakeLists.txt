cmake_minimum_required(VERSION 3.26...3.28)
project(${SKBUILD_PROJECT_NAME} LANGUAGES C CXX)

find_package(
  Python
  COMPONENTS Interpreter Development.SABIModule
  # NumPy
  REQUIRED)

# webpack
set(WEBPACK_DIR ${CMAKE_CURRENT_SOURCE_DIR}/reacnetgenerator/static/webpack)
# find node
find_program(
  NODE_EXECUTABLE
  NAMES node
  DOC "Node.js executable" REQUIRED)
message(STATUS "Found Node.js: ${NODE_EXECUTABLE}")
# use python to get yarn path from .yarnrc.yml
execute_process(
  COMMAND ${Python_EXECUTABLE} -c
          "import yaml; print(yaml.safe_load(open('.yarnrc.yml'))['yarnPath'])"
  WORKING_DIRECTORY ${WEBPACK_DIR}
  OUTPUT_VARIABLE YARN_EXECUTABLE
  OUTPUT_STRIP_TRAILING_WHITESPACE)
if(NOT YARN_EXECUTABLE)
  message(FATAL_ERROR "Yarn executable not found")
endif()
message(STATUS "Found Yarn: ${YARN_EXECUTABLE}")

add_custom_command(
  OUTPUT bundle.html
  COMMAND ${NODE_EXECUTABLE} ${YARN_EXECUTABLE}
  COMMAND ${NODE_EXECUTABLE} ${YARN_EXECUTABLE} start
  COMMENT "Building bundle.html using webpack"
  WORKING_DIRECTORY ${WEBPACK_DIR}
  VERBATIM)
add_custom_target(run ALL DEPENDS bundle.html)
# install generated bundle.html
install(
  FILES ${CMAKE_CURRENT_SOURCE_DIR}/reacnetgenerator/static/webpack/bundle.html
  DESTINATION reacnetgenerator/static/webpack/)

# dps
add_custom_command(
  OUTPUT dps.cpp
  COMMENT
    "Making ${CMAKE_CURRENT_BINARY_DIR}/dps.cpp from ${CMAKE_CURRENT_SOURCE_DIR}/reacnetgenerator/dps.pyx"
  COMMAND
    Python::Interpreter -m cython
    "${CMAKE_CURRENT_SOURCE_DIR}/reacnetgenerator/dps.pyx" --output-file dps.cpp
  DEPENDS reacnetgenerator/dps.pyx
  VERBATIM)

python_add_library(
  dps
  MODULE
  USE_SABI
  3.7
  ${CMAKE_CURRENT_BINARY_DIR}/dps.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/reacnetgenerator/c_stack.cpp
  WITH_SOABI)
target_include_directories(dps
                           PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/reacnetgenerator)
target_compile_definitions(dps PRIVATE CYTHON_LIMITED_API=1)
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  target_compile_definitions(dps PRIVATE CYTHON_TRACE=1 CYTHON_TRACE_NOGIL=1)
endif()

# utils_np temporary disable until limited api is supported
if(FALSE)
  add_custom_command(
    OUTPUT utils_np.c
    COMMENT
      "Making ${CMAKE_CURRENT_BINARY_DIR}/utils_np.cpp from ${CMAKE_CURRENT_SOURCE_DIR}/reacnetgenerator/utils_np.pyx"
    COMMAND
      Python::Interpreter -m cython
      "${CMAKE_CURRENT_SOURCE_DIR}/reacnetgenerator/utils_np.pyx" --output-file
      utils_np.c
    DEPENDS reacnetgenerator/utils_np.pyx
    VERBATIM)

  python_add_library(utils_np MODULE USE_SABI 3.11
                     ${CMAKE_CURRENT_BINARY_DIR}/utils_np.c WITH_SOABI)
  target_link_libraries(utils_np PRIVATE Python::NumPy)
  target_compile_definitions(utils_np PRIVATE CYTHON_LIMITED_API=1)
  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_definitions(utils_np PRIVATE CYTHON_TRACE=1
                                                CYTHON_TRACE_NOGIL=1)
  endif()
  install(TARGETS utils_np DESTINATION reacnetgenerator/)
endif()
install(TARGETS dps DESTINATION reacnetgenerator/)
