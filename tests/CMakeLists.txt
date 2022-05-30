# CMakeList.txt : CMake project for CracklingPlusPlus
cmake_minimum_required (VERSION 3.8)

# Enable testing
enable_testing()

# Seperate file with doctest main() entry point
add_library(testMain OBJECT testMain.cpp)

# Adding Helpers tests excutable
add_executable(testHelpers Helpers.cpp $<TARGET_OBJECTS:testMain>)

# Added Helpers test
add_test(NAME testHelpers COMMAND testHelpers)
