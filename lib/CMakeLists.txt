# CMakeList.txt : CMake file to build required libraries
###########
# HELPERS #
###########
# Add Library
add_library(Helpers STATIC ${CMAKE_SOURCE_DIR}/lib/Helpers.cpp ${CMAKE_SOURCE_DIR}/include/Helpers.hpp)


##########
# LIBSVM #
##########
# Add Library
add_library(libsvm STATIC ${CMAKE_SOURCE_DIR}/lib/svm.cpp ${CMAKE_SOURCE_DIR}/include/svm.h)


##########
# LOGGER #
##########
# Add Library
add_library(Logger STATIC ${CMAKE_SOURCE_DIR}/lib/Logger.cpp ${CMAKE_SOURCE_DIR}/include/Logger.hpp)


##################
# INPUTPROCESSOR #
##################
# Add Library
add_library(inputProcessor STATIC ${CMAKE_SOURCE_DIR}/lib/inputProcessor.cpp ${CMAKE_SOURCE_DIR}/include/inputProcessor.hpp)


#################
# CONFIGMANAGER #
#################
# Add Library
add_library(ConfigManager STATIC ${CMAKE_SOURCE_DIR}/lib/ConfigManager.cpp ${CMAKE_SOURCE_DIR}/include/ConfigManager.hpp)
# Link required librarires 
target_link_libraries (ConfigManager PRIVATE Helpers)


######################
# CAS9INPUTPROCESSOR #
######################
# Add Library
add_library(cas9InputProcessor STATIC ${CMAKE_SOURCE_DIR}/lib/cas9InputProcessor.cpp ${CMAKE_SOURCE_DIR}/include/cas9InputProcessor.hpp)
# Link required librarires 
target_link_libraries (cas9InputProcessor PRIVATE Helpers PRIVATE ConfigManager PRIVATE inputProcessor)


############
# CHOPCHOP #
############
# Add Library
add_library(CHOPCHOP STATIC ${CMAKE_SOURCE_DIR}/lib/CHOPCHOP.cpp ${CMAKE_SOURCE_DIR}/include/CHOPCHOP.hpp)
# Link required librarires 
target_link_libraries (CHOPCHOP PRIVATE Helpers PRIVATE ConfigManager)


##########
# MM10DB #
##########
# Add Library
add_library(mm10db STATIC ${CMAKE_SOURCE_DIR}/lib/mm10db.cpp ${CMAKE_SOURCE_DIR}/include/mm10db.hpp)
# Link required librarires 
target_link_libraries (mm10db PRIVATE Helpers PRIVATE ConfigManager)



################
# SGRNASCORER2 #
################
# Add Library
add_library(sgrnascorer2 STATIC ${CMAKE_SOURCE_DIR}/lib/sgrnascorer2.cpp ${CMAKE_SOURCE_DIR}/include/sgrnascorer2.hpp)
# Link required librarires 
target_link_libraries (sgrnascorer2 PRIVATE Helpers PRIVATE ConfigManager PRIVATE libsvm)


###########
# BOWTIE2 #
###########
# Add Library
add_library(bowtie2 STATIC ${CMAKE_SOURCE_DIR}/lib/bowtie2.cpp ${CMAKE_SOURCE_DIR}/include/bowtie2.hpp)
# Link required librarires 
target_link_libraries (bowtie2 PRIVATE Helpers PRIVATE ConfigManager)


########################
# ISSLOFFTARGETSCORING #
########################
# Add Library
add_library(ISSLOffTargetScoring STATIC ${CMAKE_SOURCE_DIR}/lib/ISSLOffTargetScoring.cpp ${CMAKE_SOURCE_DIR}/include/ISSLOffTargetScoring.hpp)
# Link required librarires 
target_link_libraries (ISSLOffTargetScoring PRIVATE Helpers PRIVATE ConfigManager)
# Find OpenMP
find_package(OpenMP)
# Link OpenMP if available
if(OpenMP_CXX_FOUND)
    target_link_libraries(ISSLOffTargetScoring PRIVATE OpenMP::OpenMP_CXX)
endif()