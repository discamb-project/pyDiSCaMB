cmake_minimum_required(VERSION 3.15)

project(TAAM_SF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(lib/discamb)

add_executable(discamb_wrapper src/discamb_wrapper.cpp)
target_include_directories(discamb_wrapper PUBLIC include)
target_include_directories(discamb_wrapper PRIVATE lib/discamb/include)
target_link_libraries(discamb_wrapper discamb)
