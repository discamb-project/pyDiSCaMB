cmake_minimum_required(VERSION 3.15)

project(
  ${SKBUILD_PROJECT_NAME}
  VERSION ${SKBUILD_PROJECT_VERSION}
)

# pybind11 
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

# discamb
add_library(discamb_wrapper src/discamb_wrapper.cpp)
target_include_directories(discamb_wrapper PUBLIC include)
target_include_directories(discamb_wrapper PRIVATE lib/discamb/include)
link_directories(lib/discamb/build/build/lib)
target_link_libraries(discamb_wrapper discamb)

# python module
python_add_library(_taam_sf MODULE src/TAAM_SF.cpp)
target_link_libraries(_taam_sf PRIVATE pybind11::headers)
target_link_libraries(_taam_sf PRIVATE discamb_wrapper)
install(TARGETS _taam_sf DESTINATION taam_sf)
