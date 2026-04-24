
cmake_minimum_required(VERSION 3.10)

project(HPIPM_Example C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0 -g")

if(APPLE)
    # add_compile_options(-fsanitize=address -fno-omit-frame-pointer)
    # add_link_options(-fsanitize=address)
endif()

set(HPIPM_INCLUDE_DIR "/opt/hpipm/include")
set(BLASFEO_INCLUDE_DIR "/opt/blasfeo/include")
link_directories(
    "/opt/hpipm/lib" 
    "/opt/blasfeo/lib")


# Add executables for different examples
# set(EXAMPLES toyExample robotArm helicopter springMass)

foreach(example IN LISTS EXAMPLES)
    add_executable(hpipm_${example} main.c)
    target_include_directories(hpipm_${example} PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/${example}
        ${HPIPM_INCLUDE_DIR}
        ${BLASFEO_INCLUDE_DIR}
    )
    target_link_libraries(hpipm_${example} hpipm blasfeo m)
endforeach()