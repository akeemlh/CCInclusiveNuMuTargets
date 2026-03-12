add_library(util SafeROOTName.cpp GetFluxIntegral.cpp GetPlaylist.cpp)
target_link_libraries(util ${ROOT_LIBRARIES} MAT-MINERvA)
install(TARGETS util DESTINATION lib)
