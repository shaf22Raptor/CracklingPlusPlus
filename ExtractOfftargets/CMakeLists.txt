# Boost
find_package(Boost REQUIRED COMPONENTS regex)

# OpenMP
find_package(OpenMP REQUIRED)

# Add ExtractOfftarget executable
add_executable (ExtractOfftargets ExtractOfftargets.cpp ExtractOfftargets.hpp)
target_link_libraries(ExtractOfftargets PRIVATE Boost::boost Boost::regex OpenMP::OpenMP_CXX utils)