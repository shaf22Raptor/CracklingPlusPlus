# CMakeList.txt : Main CMake file for CracklingPlusPlus
cmake_minimum_required (VERSION 3.8)

# Set CXX Standard for project
set(CMAKE_CXX_STANDARD 17 CACHE STRING "The C++ standard to use")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# TODO: Add flags for debug and release
# Set Compiler Flags for project
# add_compile_options(/-O2)

# Project name, version and description
project ("CracklingPlusPlus" LANGUAGES CXX C VERSION 1.0.0 DESCRIPTION "C++ version of Crackling" )

# Add source to this project's executable.
add_executable (${PROJECT_NAME} ${CMAKE_SOURCE_DIR}/src/CracklingPlusPlus.cpp)

# Add library subdirectory
add_subdirectory (${CMAKE_SOURCE_DIR}/lib)

# Link Libraries to main executable
target_link_libraries (
	${PROJECT_NAME} 
	PRIVATE Helpers 
	PRIVATE ConfigManager 
	PRIVATE Logger 
	PRIVATE cas9InputProcessor 
	PRIVATE CHOPCHOP 
	PRIVATE mm10db 
	PRIVATE sgrnascorer2 
	PRIVATE bowtie2 
	PRIVATE ISSLOffTargetScoring
)

# Added ISSLCreateIndex subdirectory
add_subdirectory (${CMAKE_SOURCE_DIR}/ISSLCreateIndex)

# Add tests subdirectory
add_subdirectory (${CMAKE_SOURCE_DIR}/tests)

# TODO: Sort out testing data in single location
# Copy testing data
file(COPY ${CMAKE_SOURCE_DIR}/tests/data/ DESTINATION data/)
