# Add ISSLCreateIndex executable
add_executable (ISSLCreateIndex ISSLCreateIndex.cpp ISSLCreateIndex.hpp)

# Find OpenMP
find_package(OpenMP)
# Link OpenMP if available
if(OpenMP_CXX_FOUND)
    target_link_libraries(ISSLCreateIndex PRIVATE OpenMP::OpenMP_CXX)
endif()