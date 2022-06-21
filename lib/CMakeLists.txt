# Add Helperlib
add_library(helperLib STATIC ${CMAKE_SOURCE_DIR}/lib/Helpers.cpp ${CMAKE_SOURCE_DIR}/include/Helpers.hpp)
# Link Libraries 
link_libraries (helperLib)

# Add ConfigManager
add_library(ConfigManager STATIC ${CMAKE_SOURCE_DIR}/lib/ConfigManager.cpp ${CMAKE_SOURCE_DIR}/include/ConfigManager.hpp)
# Link Libraries 
link_libraries (ConfigManager)

# Add Logger
add_library(Logger STATIC ${CMAKE_SOURCE_DIR}/lib/Logger.cpp ${CMAKE_SOURCE_DIR}/include/Logger.hpp)

# Add inputProcessor
add_library(inputProcessor STATIC ${CMAKE_SOURCE_DIR}/lib/inputProcessor.cpp ${CMAKE_SOURCE_DIR}/include/inputProcessor.hpp)

# Add cas9InputProcessor
add_library(cas9InputProcessor STATIC ${CMAKE_SOURCE_DIR}/lib/cas9InputProcessor.cpp ${CMAKE_SOURCE_DIR}/include/cas9InputProcessor.hpp)
# Link Libraries 
target_link_libraries (cas9InputProcessor inputProcessor)

# Add CHOPCHOP
add_library(CHOPCHOP STATIC ${CMAKE_SOURCE_DIR}/lib/CHOPCHOP.cpp ${CMAKE_SOURCE_DIR}/include/CHOPCHOP.hpp)

# Add mm10db
add_library(mm10db STATIC ${CMAKE_SOURCE_DIR}/lib/mm10db.cpp ${CMAKE_SOURCE_DIR}/include/mm10db.hpp)

# Add bowtie2
add_library(bowtie2 STATIC ${CMAKE_SOURCE_DIR}/lib/bowtie2.cpp ${CMAKE_SOURCE_DIR}/include/bowtie2.hpp)

# Add offTargetScoring
add_library(offTargetScoring STATIC ${CMAKE_SOURCE_DIR}/lib/offTargetScoring.cpp ${CMAKE_SOURCE_DIR}/include/offTargetScoring.hpp)