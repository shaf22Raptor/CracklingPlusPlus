# CMakeList.txt : CMake project for CracklingPlusPlus
cmake_minimum_required (VERSION 3.8)

set (CMAKE_CXX_STANDARD 17)

# Project name, version and description
project ("CracklingPlusPlus" VERSION 1.0.0 DESCRIPTION "C++ version of Crackling")

# Add source to this project's executable.
add_executable (${PROJECT_NAME} ${CMAKE_SOURCE_DIR}/src/CracklingPlusPlus.cpp)

# Add include directory 
include_directories (${CMAKE_SOURCE_DIR}/include)

# Add library subdirectory
add_subdirectory (${CMAKE_SOURCE_DIR}/lib)

# Link library subdirectory
link_directories (${CMAKE_SOURCE_DIR}/lib)

# Link Libraries 
target_link_libraries (${PROJECT_NAME} helperLib ConfigManager Logger cas9InputProcessor CHOPCHOP mm10db sgrnascorer2 bowtie2 offTargetScoring)

# Add tests subdirectory
add_subdirectory (${CMAKE_SOURCE_DIR}/tests)

# Copy testing data
file(COPY ${CMAKE_SOURCE_DIR}/tests/data/ DESTINATION data/)

set_target_properties(${PROJECT_NAME} PROPERTIES 
    VS_DEBUGGER_COMMAND ${PROJECT_NAME}
    VS_DEBUGGER_COMMAND_ARGUMENTS "data/test_config.ini")