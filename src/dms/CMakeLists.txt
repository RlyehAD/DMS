include_directories(${CMAKE_CURRENT_LIST_DIR} ${PETSC_INCLUDES})
file(GLOB DMS_SOURCES "*.h" "*.cpp")
add_library(dms_objlib OBJECT ${DMS_SOURCES})

add_subdirectory("methods")
add_subdirectory("methods/space_warping")
add_subdirectory("integrators")
#add_subdirectory("methods/field_variables")
