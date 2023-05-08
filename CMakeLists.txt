cmake_minimum_required(VERSION 3.22)
project(EVS_codec C)

# Detect all source files
aux_source_directory(lib_enc SRC_LIB_ENC)
aux_source_directory(lib_dec SRC_LIB_DEC)
aux_source_directory(lib_com SRC_LIB_COM)

set(ENCODER_SOURCE ${SRC_LIB_ENC} ${SRC_LIB_COM})
set(DECODER_SOURCE ${SRC_LIB_DEC} ${SRC_LIB_COM})

# Add codec executable
add_executable(EVS_cod ${ENCODER_SOURCE})
add_executable(EVS_dec ${DECODER_SOURCE})

include_directories(lib_enc lib_dec lib_com)

target_link_libraries(EVS_cod PRIVATE m)
target_link_libraries(EVS_dec PRIVATE m)
