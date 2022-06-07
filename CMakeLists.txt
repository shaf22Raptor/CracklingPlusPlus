# CMakeList.txt : CMake project for CracklingPlusPlus
cmake_minimum_required (VERSION 3.8)

set (CMAKE_CXX_STANDARD 17)

# Project name, version and description
project ("CracklingPlusPlus" VERSION 1.0.0 DESCRIPTION "C++ version of Crackling")

# Copy file
configure_file(${CMAKE_SOURCE_DIR}/src/test_config.ini ${CMAKE_BINARY_DIR} COPYONLY)

# Add include directory 
include_directories (${CMAKE_SOURCE_DIR}/include)

# Add library subdirectory
add_subdirectory (${CMAKE_SOURCE_DIR}/lib)

# Link library subdirectory
link_directories (${CMAKE_SOURCE_DIR}/lib)

# Link HelperLib
link_libraries (helperLib)

# Link ConfigManager
link_libraries (ConfigManager)

# Link Logger Class
link_libraries (Logger)

# Link CHOPCHOP Class
link_libraries (CHOPCHOP)

# Add source to this project's executable.
add_executable (${PROJECT_NAME} ${CMAKE_SOURCE_DIR}/src/CracklingPlusPlus.cpp)

# Add tests subdirectory
add_subdirectory (${CMAKE_SOURCE_DIR}/tests)