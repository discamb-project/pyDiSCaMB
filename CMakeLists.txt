cmake_minimum_required(VERSION 3.15)

project(
  ${SKBUILD_PROJECT_NAME}
  VERSION ${SKBUILD_PROJECT_VERSION}
)
set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)

# pybind11 
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

# discamb
FIND_LIBRARY(DISCAMB_LIB discamb PATH ${PROJECT_BINARY_DIR}/build/lib)
if(NOT DISCAMB_LIB)
  message("DiSCaMB not found, building...")
  option(BUILD_EXAMPLES OFF)
  option(BUILD_SHARED_LIBS ON)
  add_subdirectory(lib/discamb)
else()
  message("DiSCaMB found: ${DISCAMB_LIB}")
  add_library(discamb SHARED IMPORTED)
  set_target_properties(discamb PROPERTIES IMPORTED_LOCATION ${DISCAMB_LIB})
  configure_file ("${PROJECT_SOURCE_DIR}/lib/discamb/config.h.in" "${PROJECT_BINARY_DIR}/lib/discamb/build/include/discamb/config.h")
endif()
find_package(OpenMP)

# python module
python_add_library(_taam_sf MODULE src/TAAM_SF.cpp src/discamb_wrapper.cpp)
target_link_libraries(_taam_sf PRIVATE pybind11::headers)
target_link_libraries(_taam_sf PUBLIC discamb PRIVATE OpenMP::OpenMP_CXX)
target_include_directories(_taam_sf PRIVATE lib/discamb/include)
target_include_directories(_taam_sf PRIVATE lib/discamb/external)
target_include_directories(_taam_sf PRIVATE ${PROJECT_BINARY_DIR}/lib/discamb/build/include) # For the config.h
target_include_directories(_taam_sf PUBLIC include)
install(TARGETS _taam_sf DESTINATION taam_sf)
