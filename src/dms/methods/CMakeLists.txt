include_directories(${CMAKE_CURRENT_LIST_DIR})

file(GLOB METH_SOURCES "*.h" "*.cpp")
add_library(meth_objlib OBJECT ${METH_SOURCES})
