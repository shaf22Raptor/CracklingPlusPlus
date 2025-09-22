# Boost
find_package(Boost REQUIRED COMPONENTS regex iostreams)

# OpenMP
find_package(OpenMP REQUIRED)

# Add ISSLScoreOfftargets executable
add_executable(ISSLScoreOfftargets
    ISSLScoreOfftargets.cpp
    ISSLScoreOfftargets.hpp
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
  # Try vcpkg config mode first, fallback to module mode
  find_package(OpenCL CONFIG QUIET)
  if(NOT OpenCL_FOUND)
    find_package(OpenCL REQUIRED)
  endif()

  target_compile_definitions(ISSLScoreOfftargets PRIVATE USE_OPENCL_SEQ2SIG=1)

  # Add the OpenCL host wrapper
  target_sources(ISSLScoreOfftargets PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}/cl_seqsig.cpp
  )
  target_include_directories(ISSLScoreOfftargets PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}
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
