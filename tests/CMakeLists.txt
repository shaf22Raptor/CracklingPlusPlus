# CMakeList.txt : CMake file to build testing suite
# Copy testing data
file (COPY ${CMAKE_SOURCE_DIR}/tests/data/ DESTINATION data/)

# Enable testing
enable_testing()

# Seperate file with doctest main() entry point
add_library (testMain OBJECT testMain.cpp)

# Add Constants test executable
add_executable (testConstants Constants.cpp $<TARGET_OBJECTS:testMain>)

# Add Constants test
add_test (NAME testConstants COMMAND testConstants)

# Add Helpers tests executable
add_executable (testHelpers Helpers.cpp $<TARGET_OBJECTS:testMain>)

# Link Libraries 
target_link_libraries (testHelpers Helpers)

# Add Helpers test
add_test (NAME testHelpers COMMAND testHelpers)

# Add CHOPCHOP tests executable
add_executable (testCHOPCHOP CHOPCHOP.cpp $<TARGET_OBJECTS:testMain>)

# Link Libraries 
target_link_libraries (testCHOPCHOP CHOPCHOP)

# Add CHOPCHOP test
add_test (NAME testCHOPCHOP COMMAND testCHOPCHOP)

# Add mm10db tests executable
add_executable (testmm10db mm10db.cpp $<TARGET_OBJECTS:testMain>)

# Link Libraries 
target_link_libraries (testmm10db mm10db)

# Add mm10db test
add_test (NAME testmm10db COMMAND testmm10db)

# Add sgrnascorer2 tests executable
add_executable (testsgrnascorer2 sgrnascorer2.cpp $<TARGET_OBJECTS:testMain>)

# Link Libraries 
target_link_libraries (testsgrnascorer2 sgrnascorer2)

# Add sgrnascorer2 test
add_test (NAME testsgrnascorer2 COMMAND testsgrnascorer2)

# Add bowtie2 tests executable
add_executable (testbowtie2 bowtie2.cpp $<TARGET_OBJECTS:testMain>)

# Link Libraries 
target_link_libraries (testbowtie2 bowtie2)

# Add sgrnascorer2 test
add_test (NAME testbowtie2 COMMAND testbowtie2)
