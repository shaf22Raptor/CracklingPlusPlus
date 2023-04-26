# Add ExtractOfftarget executable
add_executable (ExtractOfftargets ExtractOfftargets.cpp ExtractOfftargets.hpp)

#########
# Boost #
#########
find_package(Boost REQUIRED COMPONENTS regex)
target_link_libraries(ExtractOfftargets PRIVATE Boost::dynamic_linking PRIVATE Boost::boost PRIVATE Boost::regex PRIVATE Boost::iostreams PRIVATE utils)

# Find OpenMP
find_package(OpenMP)
# Link OpenMP if available
if(OpenMP_CXX_FOUND)
    target_link_libraries(ExtractOfftargets PRIVATE OpenMP::OpenMP_CXX)
endif()