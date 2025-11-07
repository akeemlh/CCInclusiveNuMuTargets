cmake_minimum_required( VERSION 2.8.12 )

project( MAT_IncPions )

#Add extra "modules" to CMake
include(CTest)
include(AddFileDependencies)

#Compiler flags
set( GCC_Flags_For_CXX "-std=c++1y -Wall" )
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_Flags_For_CXX}" )
set( CMAKE_CXX_FLAGS_DEBUG "-ggdb" )
set( CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG" )

#Add a PROFiling CMAKE_BUILD_TYPE
set( CMAKE_CXX_FLAGS_PROF "-ggdb -pg -DNDEBUG" )
set( CMAKE_C_FLAGS_PROF "-ggdb -pg -DNDEBUG" )
set( CMAKE_EXE_LINKER_FLAGS_PROF "-ggdb -pg -DNDEBUG" )
set( CMAKE_SHARED_LINKER_FLAGS_PROF "-ggdb -pg -DNDEBUG" )

#Tell this package where it is installed and version control status
add_definitions(-DINSTALL_DIR="${CMAKE_INSTALL_PREFIX}/")

#Let directories in this package see each other
include_directories( "${PROJECT_SOURCE_DIR}" )

#Find dependencies
list( APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net Tree Core Geom EG GenVector Minuit2 Minuit OPTIONAL_COMPONENTS Cintex)

if(${ROOT_VERSION} VERSION_LESS 6 AND NOT ${ROOT_Cintex_FOUND})
  MESSAGE(FATAL_ERROR "Cintex is optional except when it's not.  ROOT 6 has Reflex "
                      "support built in, so it doesn't have a separate component "
                      "for Cintex.  Cintex was an experimental feature in ROOT 5, "
                      "so I have to require it as a component there.  You appear to "
                      "be using ${ROOT_VERSION}, so I can't let you get away with "
                      "skipping out on Cintex support!  I need Cintex to load metadata "
                      "about MnvH1D and MnvH2D at runtime.")
endif()

if(${ROOT_VERSION} VERSION_LESS 6)
  MESSAGE("I need Cintex on this platform, so making sure to compile against it...")
else()
  MESSAGE("I don't need Cintex on this platform, so skipping it...")
  add_definitions(-DNCINTEX)
endif()

include(${ROOT_USE_FILE})

find_package(MAT REQUIRED)
include_directories(${MAT_INCLUDE_DIR})

find_package(MAT-MINERvA REQUIRED)
include_directories(${MAT-MINERvA_INCLUDE_DIR})

find_package(UnfoldUtils REQUIRED)
include_directories(${UnfoldUtils_INCLUDE_DIR})

find_package(GENIEXSecExtract REQUIRED)
include_directories(${GENIEXSecExtract_INCLUDE_DIR})

#Tell CMake to compile subdirectories before compiling main application
add_subdirectory(playlists)
add_subdirectory(util)
#TODO: Split these directories' headers into .cpp files
#add_subdirectory(event)
#add_subdirectory(cuts)
#add_subdirectory(studies)
#add_subdirectory(systematics)

#Build main executables
add_executable(runEventLoop runEventLoop.cpp)
target_link_libraries(runEventLoop ${ROOT_LIBRARIES} util MAT MAT-MINERvA UnfoldUtils) #event cuts studies systematics)
install(TARGETS runEventLoop DESTINATION bin)

#add_executable(runEventLoopTargets runEventLoopTargets.cpp)
#target_link_libraries(runEventLoopTargets ${ROOT_LIBRARIES} util MAT MAT-MINERvA UnfoldUtils) #event cuts studies systematics)
#install(TARGETS runEventLoopTargets DESTINATION bin)

add_executable(runEventLoopTargetsNew runEventLoopTargetsNew.cpp)
target_link_libraries(runEventLoopTargetsNew ${ROOT_LIBRARIES} util MAT MAT-MINERvA UnfoldUtils) #event cuts studies systematics)
install(TARGETS runEventLoopTargetsNew DESTINATION bin)

add_executable(runEventLoopTracker runEventLoopTracker.cpp)
target_link_libraries(runEventLoopTracker ${ROOT_LIBRARIES} util MAT MAT-MINERvA UnfoldUtils) #event cuts studies systematics)
install(TARGETS runEventLoopTracker DESTINATION bin)

#add_executable(runEventLoopTrackerDSWater runEventLoopTrackerDSWater.cpp)
#target_link_libraries(runEventLoopTrackerDSWater ${ROOT_LIBRARIES} util MAT MAT-MINERvA UnfoldUtils) #event cuts studies systematics)
#install(TARGETS runEventLoopTrackerDSWater DESTINATION bin)

#add_executable(runEventLoopTrackerUSWater runEventLoopTrackerUSWater.cpp)
#target_link_libraries(runEventLoopTrackerUSWater ${ROOT_LIBRARIES} util MAT MAT-MINERvA UnfoldUtils) #event cuts studies systematics)
#install(TARGETS runEventLoopTrackerUSWater DESTINATION bin)

#add_executable(runEventLoopValidations runEventLoopValidations.cpp)
#target_link_libraries(runEventLoopValidations ${ROOT_LIBRARIES} util MAT MAT-MINERvA UnfoldUtils) #event cuts studies systematics)
#install(TARGETS runEventLoopValidations DESTINATION bin)

add_executable(runEventLoopValidationsNew runEventLoopValidationsNew.cpp)
target_link_libraries(runEventLoopValidationsNew ${ROOT_LIBRARIES} util MAT MAT-MINERvA UnfoldUtils) #event cuts studies systematics)
install(TARGETS runEventLoopValidationsNew DESTINATION bin)

#add_executable(Extract1DCrossSectionTargets Extract1DCrossSectionTargets.cpp)
#target_link_libraries(Extract1DCrossSectionTargets ${ROOT_LIBRARIES} util MAT UnfoldUtils)
#install(TARGETS Extract1DCrossSectionTargets DESTINATION bin)

#add_executable(Extract1DCrossSectionTargets_ByTarget Extract1DCrossSectionTargets_ByTarget.cpp)
#target_link_libraries(Extract1DCrossSectionTargets_ByTarget ${ROOT_LIBRARIES} util MAT UnfoldUtils)
#install(TARGETS Extract1DCrossSectionTargets_ByTarget DESTINATION bin)

add_executable(Extract1DCrossSectionTargets_ByTargetNew Extract1DCrossSectionTargets_ByTargetNew.cpp)
target_link_libraries(Extract1DCrossSectionTargets_ByTargetNew ${ROOT_LIBRARIES} util MAT UnfoldUtils)
install(TARGETS Extract1DCrossSectionTargets_ByTargetNew DESTINATION bin)

#add_executable(Extract1DCrossSectionTargets_ByMaterial Extract1DCrossSectionTargets_ByMaterial.cpp)
#target_link_libraries(Extract1DCrossSectionTargets_ByMaterial ${ROOT_LIBRARIES} util MAT UnfoldUtils)
#install(TARGETS Extract1DCrossSectionTargets_ByMaterial DESTINATION bin)

#add_executable(Extract1DCrossSectionTargets_PseudoTargetsSummed Extract1DCrossSectionTargets_PseudoTargetsSummed.cpp)
#target_link_libraries(Extract1DCrossSectionTargets_PseudoTargetsSummed ${ROOT_LIBRARIES} util MAT UnfoldUtils)
#install(TARGETS Extract1DCrossSectionTargets_PseudoTargetsSummed DESTINATION bin)

add_executable(Extract1DCrossSectionTracker Extract1DCrossSectionTracker.cpp)
target_link_libraries(Extract1DCrossSectionTracker ${ROOT_LIBRARIES} util MAT UnfoldUtils)
install(TARGETS Extract1DCrossSectionTracker DESTINATION bin)

add_executable(Extract1DCrossSectionTrackerNew Extract1DCrossSectionTrackerNew.cpp)
target_link_libraries(Extract1DCrossSectionTrackerNew ${ROOT_LIBRARIES} util MAT UnfoldUtils)
install(TARGETS Extract1DCrossSectionTrackerNew DESTINATION bin)

add_executable(Extract1DCrossSectionTrackerDaisy Extract1DCrossSectionTrackerDaisy.cpp)
target_link_libraries(Extract1DCrossSectionTrackerDaisy ${ROOT_LIBRARIES} util MAT UnfoldUtils)
install(TARGETS Extract1DCrossSectionTrackerDaisy DESTINATION bin)

add_executable(Extract2DCrossSectionTrackerDaisy Extract2DCrossSectionTrackerDaisy.cpp)
target_link_libraries(Extract2DCrossSectionTrackerDaisy ${ROOT_LIBRARIES} util MAT UnfoldUtils)
install(TARGETS Extract2DCrossSectionTrackerDaisy DESTINATION bin)

add_executable(Extract2DCrossSectionTargets Extract2DCrossSectionTargets.cpp)
target_link_libraries(Extract2DCrossSectionTargets ${ROOT_LIBRARIES} util MAT UnfoldUtils)
install(TARGETS Extract2DCrossSectionTargets DESTINATION bin)

add_executable(Extract2DCrossSectionTracker Extract2DCrossSectionTracker.cpp)
target_link_libraries(Extract2DCrossSectionTracker ${ROOT_LIBRARIES} util MAT UnfoldUtils)
install(TARGETS Extract2DCrossSectionTracker DESTINATION bin)


add_executable(runXSecLooper runXSecLooper.cpp)
target_link_libraries(runXSecLooper ${ROOT_LIBRARIES} MAT GENIEXSecExtract)
install(TARGETS runXSecLooper DESTINATION bin)

add_executable(runXSecLooperTracker runXSecLooperTracker.cpp)
target_link_libraries(runXSecLooperTracker ${ROOT_LIBRARIES} MAT GENIEXSecExtract)
install(TARGETS runXSecLooperTracker DESTINATION bin)

add_executable(runXSecLooperTargets_ByTarget runXSecLooperTargets_ByTarget.cpp)
target_link_libraries(runXSecLooperTargets_ByTarget ${ROOT_LIBRARIES} MAT GENIEXSecExtract)
install(TARGETS runXSecLooperTargets_ByTarget DESTINATION bin)

add_executable(runXSecLooper1DTracker runXSecLooper1DTracker.cpp)
target_link_libraries(runXSecLooper1DTracker ${ROOT_LIBRARIES} MAT GENIEXSecExtract)
install(TARGETS runXSecLooper1DTracker DESTINATION bin)

install(FILES backgroundStack.py PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE DESTINATION bin)

install(FILES SidebandTuning.py PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE DESTINATION bin)

install(FILES runTransWarp.sh PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE DESTINATION bin)

configure_file(setup.sh.in setup_${PROJECT_NAME}.sh @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/setup_${PROJECT_NAME}.sh DESTINATION bin)
