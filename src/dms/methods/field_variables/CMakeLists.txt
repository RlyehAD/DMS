include_directories(${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_LIST_DIR}/../..)

file(GLOB FVM_SOURCES "*.h" "*.cpp")
add_library(fvm_objlib OBJECT ${FVM_SOURCES})
