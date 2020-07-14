cmake_minimum_required(VERSION 2.8)
project(Sniffles)

option(STATIC "Build static binary" OFF)

 set( SNIF_VERSION_MAJOR 1 )
 set( SNIF_VERSION_MINOR 0 )
IF(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "Building in debug mode!")
    set( SNIF_VERSION_BUILD 12-debug )
else()
	set( SNIF_VERSION_BUILD 12 )
ENDIF()


set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin/sniffles-core-${SNIF_VERSION_MAJOR}.${SNIF_VERSION_MINOR}.${SNIF_VERSION_BUILD}/)
file(MAKE_DIRECTORY ${EXECUTABLE_OUTPUT_PATH})


# Set a default build type for single-configuration
# CMake generators if no build type is set.
IF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
	message(STATUS "No build type specified. Using 'release'")
	SET(CMAKE_BUILD_TYPE Release)
ENDIF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)

FIND_PACKAGE(OpenMP REQUIRED)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


add_subdirectory(lib/zlib-1.2.7)
add_subdirectory(lib/bamtools-2.3.0)

add_subdirectory(src)
