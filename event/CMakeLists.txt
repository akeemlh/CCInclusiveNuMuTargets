add_library(event CVUniverse.cpp)
target_link_libraries(event ${ROOT_LIBRARIES} MAT-MINERvA)
install(TARGETS event DESTINATION lib)
