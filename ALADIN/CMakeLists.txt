cmake_minimum_required(VERSION 3.10)
project(ALADIN_Example C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
set(CMAKE_POSITION_INDEPENDENT_CODE ON)  # -fPIC

cmake_policy(SET CMP0074 NEW)
if (APPLE AND EXISTS "/opt/homebrew/opt/libomp")
    set(OpenMP_ROOT "/opt/homebrew/opt/libomp")
endif()
find_package(OpenMP REQUIRED)

# Source files common to all executables
set(COMMON_SOURCES
    RRLB_MPC.c
)

# set(EXAMPLES toyExample robotArm helicopter)

# the list of examples is defined in the top-level CMakeLists.txt
foreach(example IN LISTS EXAMPLES)
    add_executable(ALADIN_${example} main.c ${COMMON_SOURCES})
    target_include_directories(ALADIN_${example} PRIVATE . ${example})
    target_compile_definitions(ALADIN_${example} PRIVATE "-DMPC_PROBLEM_DIR="${example})
    target_link_libraries(ALADIN_${example} PRIVATE OpenMP::OpenMP_C m)
endforeach()



# ALADIN_toyExample executable
# add_executable(ALADIN_toyExample main.c ${COMMON_SOURCES})
# target_include_directories(ALADIN_toyExample PRIVATE . toyExample)
# target_compile_definitions(ALADIN_toyExample PRIVATE "-DMPC_PROBLEM_DIR=toyExample")
# 
# # ALADIN_robotArm executable
# add_executable(ALADIN_robotArm main.c ${COMMON_SOURCES})
# target_include_directories(ALADIN_robotArm PRIVATE . robotArm)
# target_compile_definitions(ALADIN_robotArm PRIVATE "-DMPC_PROBLEM_DIR=robotArm")
# 
# add_executable(ALADIN_helicopter main.c ${COMMON_SOURCES})
# target_include_directories(ALADIN_helicopter PRIVATE . helicopter)
# target_compile_definitions(ALADIN_helicopter PRIVATE "-DMPC_PROBLEM_DIR=helicopter")

# Link libraries (if needed)
# target_link_libraries(ALADIN_toyExample PRIVATE m)
# target_link_libraries(ALADIN_robotArm PRIVATE m)
