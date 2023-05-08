cmake_minimum_required(VERSION 3.22)
project(EVS_codec C)

# Detect all source files
aux_source_directory(lib_enc SRC_LIB_ENC)
aux_source_directory(lib_dec SRC_LIB_DEC)
aux_source_directory(lib_com SRC_LIB_COM)

set(ENCODER_SOURCE ${SRC_LIB_ENC} ${SRC_LIB_COM})
set(DECODER_SOURCE ${SRC_LIB_DEC} ${SRC_LIB_COM})

# targets
add_library(libEVS_cod STATIC ${ENCODER_SOURCE})
add_library(libEVS_dec STATIC ${DECODER_SOURCE})

add_executable(EVS_cod encoder.c)
add_executable(EVS_dec decoder.c)

# include
include_directories(lib_enc lib_dec lib_com)

# Config for static libraries
set_target_properties(
        libEVS_cod
        PROPERTIES
        OUTPUT_NAME EVS_cod
)
target_link_libraries(libEVS_cod PRIVATE m)  # math.h

set_target_properties(
        libEVS_dec
        PROPERTIES
        OUTPUT_NAME EVS_dec
)
target_link_libraries(libEVS_dec PRIVATE m)

# link library to executable
target_link_libraries(EVS_cod PRIVATE libEVS_cod)

target_link_libraries(EVS_dec PRIVATE libEVS_dec)
