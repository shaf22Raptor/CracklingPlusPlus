# Added Helperlib
add_library(helperLib STATIC ${CMAKE_SOURCE_DIR}/lib/Helpers.cpp ${CMAKE_SOURCE_DIR}/include/Helpers.hpp)
# Link Libraries 
link_libraries (helperLib)

# Added ConfigManager
add_library(ConfigManager STATIC ${CMAKE_SOURCE_DIR}/lib/ConfigManager.cpp ${CMAKE_SOURCE_DIR}/include/ConfigManager.hpp)
# Link Libraries 
link_libraries (ConfigManager)

# Added Logger
add_library(Logger STATIC ${CMAKE_SOURCE_DIR}/lib/Logger.cpp ${CMAKE_SOURCE_DIR}/include/Logger.hpp)
# Link Libraries 
link_libraries (Logger)

# Added CHOPCHOP
add_library(CHOPCHOP STATIC ${CMAKE_SOURCE_DIR}/lib/CHOPCHOP.cpp ${CMAKE_SOURCE_DIR}/include/CHOPCHOP.hpp)
# Link Libraries 
link_libraries (CHOPCHOP)