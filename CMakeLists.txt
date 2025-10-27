# CMakeList.txt : Main CMake file for CracklingPlusPlus
cmake_minimum_required (VERSION 3.8)

# Project name, version and description
project ("CracklingPlusPlus" LANGUAGES CXX C VERSION 1.0.0 DESCRIPTION "C++ version of Crackling" )

# Set CXX Standard for project
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Optimisation flags
# Optimisation and debug info for profiling (cross-platform)
if (MSVC)
  # Compiler
  set(CMAKE_CXX_FLAGS_RELEASE
      "${CMAKE_CXX_FLAGS_RELEASE} /O2 /Zi /Zo /Oy- /Ob1 /FS /DNDEBUG")

  # Linker (EXE)
  set(CMAKE_EXE_LINKER_FLAGS_RELEASE
      "${CMAKE_EXE_LINKER_FLAGS_RELEASE} /DEBUG:FULL /INCREMENTAL:NO /OPT:NOICF /OPT:NOREF")

  # Linker (DLL)
  set(CMAKE_SHARED_LINKER_FLAGS_RELEASE
      "${CMAKE_SHARED_LINKER_FLAGS_RELEASE} /DEBUG:FULL /INCREMENTAL:NO /OPT:NOICF /OPT:NOREF")

  # Ensure no LTCG/IPO at CMake level
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE OFF)

  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE OFF)

  # Define NOMINMAX for all targets in this tree (propagates to subdirs)
  add_compile_definitions(NOMINMAX)

endif()

# Add library subdirectory 
add_subdirectory (${CMAKE_SOURCE_DIR}/lib)

# Boost
find_package(Boost REQUIRED COMPONENTS regex iostreams)

# Add source to this project's executable.
add_executable (${PROJECT_NAME} ${CMAKE_SOURCE_DIR}/src/CracklingPlusPlus.cpp)

# Link Libraries to main executable
target_link_libraries (
	${PROJECT_NAME}
	PRIVATE Boost::dynamic_linking
	PRIVATE Boost::boost
	PRIVATE Boost::regex
	PRIVATE Boost::iostreams
	PRIVATE Logger 
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

# Added ISSLScoreOfftargets subdirectory
add_subdirectory (${CMAKE_SOURCE_DIR}/ISSLScoreOfftargets)

target_include_directories(ISSLScoreOfftargets PRIVATE
    ${CMAKE_SOURCE_DIR}/third_party/cvmarkers
)

# Added ISSLScoreOfftargetsMMF subdirectory
add_subdirectory (${CMAKE_SOURCE_DIR}/ISSLScoreOfftargetsMMF)

# Added ExtractOffTargets subdirectory
add_subdirectory (${CMAKE_SOURCE_DIR}/ExtractOfftargets)

# Add tests subdirectory
# add_subdirectory (${CMAKE_SOURCE_DIR}/tests)

target_compile_definitions(ISSLScoreOfftargets PRIVATE MANPROF_ENABLE)

# TODO: Sort out testing data in single location
# Copy testing data
file(COPY ${CMAKE_SOURCE_DIR}/tests/data/ DESTINATION data/)
