# CMakeList.txt : CMake file to build required libraries
#########
# Boost #
#########
find_package(Boost 1.79 REQUIRED COMPONENTS regex)


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

########################
# Config Parser Module #
########################
add_library(configParserModule STATIC ${CMAKE_SOURCE_DIR}/lib/configParserModule.cpp ${CMAKE_SOURCE_DIR}/include/configParserModule.hpp ${CMAKE_SOURCE_DIR}/include/pipelineModule.hpp)
target_link_libraries (configParserModule PRIVATE Boost::boost PRIVATE Boost::regex)

################
# Input Module #
################
add_library(inputModule STATIC ${CMAKE_SOURCE_DIR}/lib/inputModule.cpp ${CMAKE_SOURCE_DIR}/include/inputModule.hpp)
target_link_libraries (inputModule PRIVATE Boost::boost PRIVATE Boost::regex)

#####################
# Cas9 Input Module #
#####################
add_library(cas9InputModule STATIC ${CMAKE_SOURCE_DIR}/lib/cas9InputModule.cpp ${CMAKE_SOURCE_DIR}/include/cas9InputModule.hpp)
target_link_libraries (cas9InputModule PRIVATE Boost::regex)

# ####################
# # Consensus Module #
# ####################
# add_library(consensusModule STATIC ${CMAKE_SOURCE_DIR}/lib/consensusModule.cpp ${CMAKE_SOURCE_DIR}/include/consensusModule.hpp ${CMAKE_SOURCE_DIR}/include/pipelineModule.hpp)


# ###################
# # CHOPCHOP Module #
# ###################
# add_library(chopchopModule STATIC ${CMAKE_SOURCE_DIR}/lib/chopchopModule.cpp ${CMAKE_SOURCE_DIR}/include/chopchopModule.hpp)
# target_link_libraries (chopchopModule PRIVATE consensusModule)


# #################
# # MM10DB Module #
# #################
# add_library(mm10dbModule STATIC ${CMAKE_SOURCE_DIR}/lib/mm10dbModule.cpp ${CMAKE_SOURCE_DIR}/include/mm10dbModule.hpp)
# target_link_libraries (mm10dbModule PRIVATE consensusModule PRIVATE Boost::boost PRIVATE Boost::regex)


# #######################
# # sgRNAScorer2 Module #
# #######################
# add_library(sgrnascorer2Module STATIC ${CMAKE_SOURCE_DIR}/lib/sgrnascorer2Module.cpp ${CMAKE_SOURCE_DIR}/include/sgrnascorer2Module.hpp)
# target_link_libraries (sgrnascorer2Module PRIVATE consensusModule)


# ######################
# # Specificity Module #
# ######################
# add_library(specificityModule STATIC ${CMAKE_SOURCE_DIR}/lib/specificityModule.cpp ${CMAKE_SOURCE_DIR}/include/specificityModule.hpp ${CMAKE_SOURCE_DIR}/include/pipelineModule.hpp)

##################
# Bowtie2 Module #
##################
add_library(bowtie2Module STATIC ${CMAKE_SOURCE_DIR}/lib/bowtie2Module.cpp ${CMAKE_SOURCE_DIR}/include/bowtie2Module.hpp)
target_link_libraries (bowtie2Module PRIVATE Boost::boost)


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


##################
# ISSLClustering #
##################
# Add Library
add_library(ISSLClustering STATIC ${CMAKE_SOURCE_DIR}/lib/ISSLClustering.cpp ${CMAKE_SOURCE_DIR}/include/ISSLClustering.hpp)
# Link required librarires 
target_link_libraries (ISSLClustering PRIVATE Helpers PRIVATE ConfigManager)
# Find OpenMP
find_package(OpenMP)
# Link OpenMP if available
if(OpenMP_CXX_FOUND)
    target_link_libraries(ISSLClustering PRIVATE OpenMP::OpenMP_CXX)
endif()

##############
# ISSL2Stage #
##############
# Add Library
add_library(ISSL2Stage STATIC ${CMAKE_SOURCE_DIR}/lib/ISSL2Stage.cpp ${CMAKE_SOURCE_DIR}/include/ISSL2Stage.hpp)
# Link required librarires 
target_link_libraries (ISSL2Stage PRIVATE Helpers PRIVATE ConfigManager)
# Find OpenMP
find_package(OpenMP)
# Link OpenMP if available
if(OpenMP_CXX_FOUND)
    target_link_libraries(ISSL2Stage PRIVATE OpenMP::OpenMP_CXX)
endif()


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
    target_link_libraries(ISSLOffTargetScoring PRIVATE OpenMP::OpenMP_CXX PRIVATE Boost::boost)
endif()

########################
# ISSLOFFTARGETSCORING #
########################
# Add Library
add_library(ISSLScoringModule STATIC ${CMAKE_SOURCE_DIR}/lib/ISSLScoringModule.cpp ${CMAKE_SOURCE_DIR}/include/ISSLScoringModule.hpp)
# Link required librarires 
target_link_libraries (ISSLScoringModule PRIVATE Helpers PRIVATE ConfigManager)
# Find OpenMP
find_package(OpenMP)
# Link OpenMP if available
if(OpenMP_CXX_FOUND)
    target_link_libraries(ISSLScoringModule PRIVATE OpenMP::OpenMP_CXX PRIVATE Boost::boost)
endif()