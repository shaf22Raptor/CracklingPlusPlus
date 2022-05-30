# CMakeList.txt : CMake project for CracklingPlusPlus
cmake_minimum_required (VERSION 3.8)

# Project name, version and description
project ("CracklingPlusPlus" VERSION 1.0.0 DESCRIPTION "C++ version of Crackling")

# Add include directory 
include_directories (${CMAKE_SOURCE_DIR}/include)

# Add library subdirectory
add_subdirectory (${CMAKE_SOURCE_DIR}/lib)

# Link library subdirectory
link_directories (${CMAKE_SOURCE_DIR}/lib)

# Link HelperLib
link_libraries (helperLib)

# Add source to this project's executable.
add_executable (${PROJECT_NAME} ${CMAKE_SOURCE_DIR}/src/CracklingPlusPlus.cpp)

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET CracklingPlusPlus PROPERTY CXX_STANDARD 20)
endif()

# Add tests subdirectory
add_subdirectory (${CMAKE_SOURCE_DIR}/tests)