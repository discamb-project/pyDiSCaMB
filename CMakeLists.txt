cmake_minimum_required(VERSION 3.15)

project(
  ${SKBUILD_PROJECT_NAME}
  VERSION ${SKBUILD_PROJECT_VERSION}
  )
set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)
  
# Windows setup
if(MSVC)
  add_definitions(/MP /MT)
endif()
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# pybind11 
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

# discamb
find_package(OpenMP)
option(USE_PRECOMPILED_DISCAMB "Whether to look for a compiled DiSCaMB library. Will compile as normal if not found" OFF)

if (USE_PRECOMPILED_DISCAMB)
  FIND_LIBRARY(DISCAMB_LIB discamb PATHS 
    ${PROJECT_BINARY_DIR}/build/lib 
    ${PROJECT_SOURCE_DIR}/lib/discamb/build/build/lib 
    ${PROJECT_SOURCE_DIR}/lib/discamb/build/build/lib/Release
  )
  if(NOT DISCAMB_LIB)
    message("DiSCaMB not found, building...")
    option(DISCAMB_BUILD_EXAMPLES OFF)
    option(BUILD_SHARED_LIBS ON)
    add_subdirectory(lib/discamb)
  else()
    message("DiSCaMB found: ${DISCAMB_LIB}")
    add_library(discamb SHARED IMPORTED)
    set_target_properties(discamb PROPERTIES IMPORTED_LOCATION ${DISCAMB_LIB})
    configure_file ("${PROJECT_SOURCE_DIR}/lib/discamb/config.h.in" "${PROJECT_BINARY_DIR}/lib/discamb/build/include/discamb/config.h")
  endif()
else()
  option(DISCAMB_BUILD_EXAMPLES OFF)
  option(BUILD_SHARED_LIBS ON)
  add_subdirectory(lib/discamb)
endif()

# python module
python_add_library(_wrapper MODULE 
  src/python_module.cpp 
  src/DiscambWrapper.cpp
  src/DiscambStructureFactorCalculator.cpp 
  src/scattering_table.cpp
  src/atom_assignment.cpp
  src/read_structure.cpp
  src/tests.cpp
)
target_link_libraries(_wrapper PRIVATE pybind11::headers)
target_link_libraries(_wrapper PRIVATE OpenMP::OpenMP_CXX)
target_link_libraries(_wrapper PUBLIC discamb)
# For some reason, the include directories are not propagated when linking
get_target_property(DISCAMB_INCLUDE_DIRECTORIES discamb INCLUDE_DIRECTORIES)
target_include_directories(_wrapper PUBLIC ${DISCAMB_INCLUDE_DIRECTORIES})
target_include_directories(_wrapper PUBLIC include)
install(TARGETS _wrapper DESTINATION pydiscamb)

# Generate stub
if (NOT MSVC) # Does not work properly on windows. Might be due to windows adding "Release" to the folder structure.
add_custom_command(TARGET _wrapper POST_BUILD
COMMAND stubgen -m _wrapper --include-docstrings -o ${CMAKE_BINARY_DIR}/stubs
BYPRODUCTS ${CMAKE_BINARY_DIR}/stubs
COMMENT "Use stubgen to create .pyi for statement completion"
)
install(DIRECTORY ${CMAKE_BINARY_DIR}/stubs/ DESTINATION pydiscamb)
endif()

# TAAM databank
file(GLOB_RECURSE TAAM_DATA_FILES "${PROJECT_SOURCE_DIR}/data/*databank.txt") 
foreach(TAAM_BANK ${TAAM_DATA_FILES})
  install(FILES ${TAAM_BANK} DESTINATION pydiscamb/data)
endforeach(TAAM_BANK ${TAAM_DATA_FILES})
