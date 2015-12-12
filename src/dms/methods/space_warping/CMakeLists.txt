include_directories(${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_LIST_DIR}/../..)

file(GLOB SWM_SOURCES "*.h" "*.cpp")
add_library(swm_objlib OBJECT ${SWM_SOURCES})
