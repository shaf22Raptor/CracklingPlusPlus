# CMakeList.txt : CMake project for CracklingPlusPlus
cmake_minimum_required (VERSION 3.8)

# Add include directory 
include_directories (${CMAKE_SOURCE_DIR}/include)

# Enable testing
enable_testing()

# Seperate file with doctest main() entry point
add_library(testMain OBJECT testMain.cpp)

# Adding Helpers tests excutable
add_executable(helpersTest Helpers.cpp $<TARGET_OBJECTS:testMain>)

# Link library subdirectory
target_link_directories (helpersTest PUBLIC ${CMAKE_SOURCE_DIR}/lib)

# Link HelperLib
target_link_libraries (helpersTest PUBLIC helperLib)

# Added Helpers test
add_test(NAME helpersTest COMMAND helpersTest)
