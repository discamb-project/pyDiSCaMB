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
option(BUILD_EXAMPLES OFF)
option(BUILD_SHARED_LIBS ON)
add_subdirectory(lib/discamb)
find_package(OpenMP)
target_link_libraries(discamb PUBLIC OpenMP::OpenMP_CXX)

# discamb wrapper
add_library(discamb_wrapper src/discamb_wrapper.cpp)
target_link_libraries(discamb_wrapper PUBLIC discamb)
target_include_directories(discamb_wrapper PUBLIC include)
target_include_directories(discamb_wrapper PUBLIC lib/discamb/include)
target_include_directories(discamb_wrapper PUBLIC lib/discamb/external)
target_include_directories(discamb_wrapper PUBLIC lib/discamb/build/build/include)

# python module
python_add_library(_taam_sf MODULE src/TAAM_SF.cpp)
target_link_libraries(_taam_sf PRIVATE pybind11::headers)
target_link_libraries(_taam_sf PRIVATE discamb_wrapper)
install(TARGETS _taam_sf DESTINATION taam_sf)
