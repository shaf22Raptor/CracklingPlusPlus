# CMakeList.txt : CMake project for CracklingPlusPlus
cmake_minimum_required (VERSION 3.8)

enable_testing()

add_library(testMain OBJECT testMain.cpp)

add_executable(helpersTest Helpers.cpp $<TARGET_OBJECTS:testMain>)
# Link library subdirectory
target_link_directories (helpersTest PUBLIC ${CMAKE_SOURCE_DIR}/lib)

# Link HelperLib
target_link_libraries (helpersTest PUBLIC helperLib)
add_test(NAME helpersTest COMMAND helpersTest)
