# CMakeList.txt : CMake project for CracklingPlusPlus

cmake_minimum_required (VERSION 3.8)

# Project name, version and description
project ("CracklingPlusPlus" VERSION 1.0.0 DESCRIPTION "C++ version of Crackling")

# Added Helperlib
add_library( HelperLib STATIC ${CMAKE_SOURCE_DIR}/lib/Helpers.cpp )

# Add source to this project's executable.
add_executable (${CMAKE_PROJECT_NAME} "src/CracklingPlusPlus.cpp" "include/Helpers.h" "lib/Helpers.cpp")

# Include Helper.h
include_directories(${CMAKE_PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/include)

# Link HelperLib
target_link_libraries(${CMAKE_PROJECT_NAME} HelperLib )


if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET CracklingPlusPlus PROPERTY CXX_STANDARD 20)
endif()

# TODO: Add tests and install targets if needed.
