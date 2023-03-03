# CMakeList.txt : Main CMake file for CracklingPlusPlus
cmake_minimum_required (VERSION 3.8)

# Set CXX Standard for project
set(CMAKE_CXX_STANDARD 17 CACHE STRING "The C++ standard to use")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Update to use target compile flags. temp solution
#if(NOT CMAKE_BUILD_TYPE)
	#set(CMAKE_BUILD_TYPE Release)
#endif()

#set(CMAKE_CXX_FLAGS "")
#set(CMAKE_CXX_FLAGS_DEBUG "-G")
#set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# Project name, version and description
project ("CracklingPlusPlus" LANGUAGES CXX C VERSION 1.0.0 DESCRIPTION "C++ version of Crackling" )

# Add source to this project's executable.
add_executable (${PROJECT_NAME} ${CMAKE_SOURCE_DIR}/src/CracklingPlusPlus.cpp)

# Add library subdirectory
add_subdirectory (${CMAKE_SOURCE_DIR}/lib)

find_package(Boost 1.79 REQUIRED COMPONENTS regex)

# Link Libraries to main executable
target_link_libraries (
	${PROJECT_NAME} 
	#PRIVATE Helpers 
	#PRIVATE ConfigManager 
	PRIVATE Logger 
	#PRIVATE cas9InputProcessor 
	#PRIVATE CHOPCHOP 
	#PRIVATE mm10db 
	#PRIVATE sgrnascorer2 
	#PRIVATE bowtie2 
	#PRIVATE ISSLOffTargetScoring
	#PRIVATE ISSLClustering
	#PRIVATE ISSL2Stage
	PRIVATE Boost::boost
	PRIVATE Boost::regex
	PRIVATE configParserModule
	PRIVATE cas9InputModule
	PRIVATE chopchopModule
	PRIVATE mm10dbModule
	PRIVATE sgrnascorer2Module
	PRIVATE bowtie2Module
	PRIVATE ISSLScoringModule
	PRIVATE ISSLScoringModuleMMF
	PRIVATE outputModule
)

# Added ISSLCreateIndex subdirectory
add_subdirectory (${CMAKE_SOURCE_DIR}/ISSLCreateIndex)

# Add tests subdirectory
add_subdirectory (${CMAKE_SOURCE_DIR}/tests)

# TODO: Sort out testing data in single location
# Copy testing data
file(COPY ${CMAKE_SOURCE_DIR}/tests/data/ DESTINATION data/)
