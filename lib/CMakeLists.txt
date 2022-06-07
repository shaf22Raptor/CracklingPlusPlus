# Added Helperlib
add_library( helperLib STATIC ${CMAKE_SOURCE_DIR}/lib/Helpers.cpp ${CMAKE_SOURCE_DIR}/include/Helpers.hpp)

# Added ConfigManager
add_library( ConfigManager STATIC ${CMAKE_SOURCE_DIR}/lib/ConfigManager.cpp ${CMAKE_SOURCE_DIR}/include/ConfigManager.hpp)

# Added Logger
add_library( Logger STATIC ${CMAKE_SOURCE_DIR}/lib/Logger.cpp ${CMAKE_SOURCE_DIR}/include/Logger.hpp)

# Added CHOPCHOP
add_library( CHOPCHOP STATIC ${CMAKE_SOURCE_DIR}/lib/CHOPCHOP.cpp ${CMAKE_SOURCE_DIR}/include/CHOPCHOP.hpp)