file(GLOB DMS_INT_SOURCES "*.h" "*.cpp")
add_library(dms_int_objlib OBJECT ${DMS_INT_SOURCES})

add_subdirectory("factorization/euler")
add_subdirectory("factorization/pade")
