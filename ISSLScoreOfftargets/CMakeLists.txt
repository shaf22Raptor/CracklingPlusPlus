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
option(SCORING_OPENCL "Enable OpenCL path for scoring" ON)

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

# -------- Kernel resource handling (runtime-loaded .cl) --------
# Keep the .cl next to the built executable so read_text_file("scoring_kernels.cl") works.
set(KERNEL_SRC ${CMAKE_CURRENT_SOURCE_DIR}/scoring_kernels.cl)

add_custom_command(
  TARGET ISSLScoreOfftargets POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
          ${KERNEL_SRC}
          $<TARGET_FILE_DIR:ISSLScoreOfftargets>/scoring_kernels.cl
  COMMENT "Copy scoring_kernels.cl next to the binary"
)

# When debugging in Visual Studio, run the exe from its output dir
# so relative paths to scoring_kernels.cl resolve.
set_target_properties(ISSLScoreOfftargets PROPERTIES
  VS_DEBUGGER_WORKING_DIRECTORY $<TARGET_FILE_DIR:ISSLScoreOfftargets>
)

# Optional: install the kernel alongside the binary on `cmake --install`
install(TARGETS ISSLScoreOfftargets RUNTIME DESTINATION bin)
install(FILES ${KERNEL_SRC} DESTINATION bin)
