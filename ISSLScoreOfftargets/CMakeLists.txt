# Boost
find_package(Boost REQUIRED COMPONENTS regex iostreams)

# OpenMP
find_package(OpenMP REQUIRED)

# Add ISSLScoreOfftargets executable
add_executable(ISSLScoreOfftargets
    ISSLScoreOfftargets.cpp
    ISSLScoreOfftargets.hpp
)

# Always compile the scoring stub (has no OpenCL deps)
target_sources(ISSLScoreOfftargets PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/cl_scoring.cpp
)
target_include_directories(ISSLScoreOfftargets PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(ISSLScoreOfftargets PRIVATE
    Boost::boost
    Boost::dynamic_linking
    Boost::regex
    Boost::iostreams
    OpenMP::OpenMP_CXX
    utils
    otScorePenalties
)

# ---- Toggle GPU seq->sig (OpenCL) ----
option(SEQ2SIG_OPENCL "Enable OpenCL path for sequence->signature" ON)

if(SEQ2SIG_OPENCL)
  find_package(OpenCL CONFIG QUIET)
  if(NOT OpenCL_FOUND)
    find_package(OpenCL REQUIRED)
  endif()

  target_compile_definitions(ISSLScoreOfftargets PRIVATE USE_OPENCL_SEQ2SIG=1)

  target_sources(ISSLScoreOfftargets PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}/cl_seqsig.cpp
  )

  if(TARGET OpenCL::OpenCL)
    target_link_libraries(ISSLScoreOfftargets PRIVATE OpenCL::OpenCL)
  else()
    target_include_directories(ISSLScoreOfftargets PRIVATE ${OpenCL_INCLUDE_DIR})
    target_link_libraries(ISSLScoreOfftargets PRIVATE ${OpenCL_LIBRARY})
  endif()
else()
  target_compile_definitions(ISSLScoreOfftargets PRIVATE USE_OPENCL_SEQ2SIG=0)
endif()
# --------------------------------------

# ---- Toggle GPU scoring (OpenCL) ----
option(SCORING_OPENCL "Enable OpenCL path for scoring" OFF)

if(SCORING_OPENCL)
  # Only find/link OpenCL if not already added by seq2sig
  if(NOT TARGET OpenCL::OpenCL AND NOT OpenCL_FOUND)
    find_package(OpenCL CONFIG QUIET)
    if(NOT OpenCL_FOUND)
      find_package(OpenCL REQUIRED)
    endif()
  endif()

  target_compile_definitions(ISSLScoreOfftargets PRIVATE USE_OPENCL_SCORING=1)

  if(TARGET OpenCL::OpenCL)
    target_link_libraries(ISSLScoreOfftargets PRIVATE OpenCL::OpenCL)
  else()
    target_include_directories(ISSLScoreOfftargets PRIVATE ${OpenCL_INCLUDE_DIR})
    target_link_libraries(ISSLScoreOfftargets PRIVATE ${OpenCL_LIBRARY})
  endif()
else()
  target_compile_definitions(ISSLScoreOfftargets PRIVATE USE_OPENCL_SCORING=0)
endif()
# ------------------------------------
