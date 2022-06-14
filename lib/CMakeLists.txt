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
# Link Libraries 
link_libraries (Logger)

# Add inputProcessor
add_library(inputProcessor STATIC ${CMAKE_SOURCE_DIR}/lib/inputProcessor.cpp ${CMAKE_SOURCE_DIR}/include/inputProcessor.hpp)
# Link Libraries 
link_libraries (inputProcessor)

# Add CHOPCHOP
add_library(CHOPCHOP STATIC ${CMAKE_SOURCE_DIR}/lib/CHOPCHOP.cpp ${CMAKE_SOURCE_DIR}/include/CHOPCHOP.hpp)
# Link Libraries 
link_libraries (CHOPCHOP)

# Add mm10db
add_library(mm10db STATIC ${CMAKE_SOURCE_DIR}/lib/mm10db.cpp ${CMAKE_SOURCE_DIR}/include/mm10db.hpp)
# Link Libraries 
link_libraries (mm10db)