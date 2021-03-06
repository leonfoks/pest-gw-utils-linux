cmake_minimum_required(VERSION 3.8)
project(gwlib LANGUAGES Fortran)

# Optional parameters to use
# -DCMAKE_Fortran_COMPILER=/path/to/fortran/compiler
# -DCMAKE_BUILD_TYPE=[DEBUG, RELEASE]  - Build with debug flags or not.
# -DBUILD_SHARED_LIBS=[ON, OFF]  - Build and shared or static library, or link a program using the appropriate flags


# Add our local modules to the module path
if(NOT CMAKE_MODULE_PATH)
  set(CMAKE_MODULE_PATH)
endif()
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../compiler_cmakefiles")

# GNUInstallDirs is used to install libraries into correct locations
# on all platforms.
include(GNUInstallDirs)

# Setup the Fortran cmake environment
# Sets up the flags for compiling shared libraries or static libraries, and returns the libType parameter.
include("${CMAKE_MODULE_PATH}/FortranEnvironment.cmake")
# include("${CMAKE_MODULE_PATH}/blas.cmake")
# include("${CMAKE_MODULE_PATH}/lapack.cmake")
include("${CMAKE_MODULE_PATH}/mkl.cmake")

# Create the library
# ${PROJECT_NAME} is the name given for the most recent project() call
add_library(${PROJECT_NAME} "")

# Add files to the library as sources
# Each CMakeLists.txt simply uses target_sources to add the source files
# Lib
include(src_dep/CMakeLists.txt)
include(src_lib/CMakeLists.txt)

target_include_directories(${PROJECT_NAME} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}> # for generated files in build mode
                                                  $<INSTALL_INTERFACE:include> # for client in install mode
                                                  $<INSTALL_INTERFACE:lib> # for generated files in install mdoe
                                           PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}) # for source files


# ================================
# Set up the make install functionality
# ================================
# 'make install' to the correct locations (provided by GNUInstallDirs).
install(TARGETS ${PROJECT_NAME} EXPORT ${PROJECT_NAME}-config
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    INCLUDES DESTINATION include/${PROJECT_NAME})

install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION include/${PROJECT_NAME}
    FILES_MATCHING PATTERN "*.mod"
    PATTERN "CMakeFiles" EXCLUDE)

install(EXPORT ${PROJECT_NAME}-config DESTINATION lib/cmake)
install(FILES ../compiler_cmakefiles/FortranEnvironment.cmake DESTINATION lib/cmake)
install(FILES ../compiler_cmakefiles/gfortran_flags.cmake DESTINATION lib/cmake)
install(FILES ../compiler_cmakefiles/intel_flags.cmake DESTINATION lib/cmake)

export(TARGETS ${PROJECT_NAME} FILE ${PROJECT_NAME}-config.cmake)
#    RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR})  # This is for Windows

# ++++++++++++++++++++++++++++++++

# ================================
# Compile the test and scaling codes
# ================================
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/src_utils)

# ++++++++++++++++++++++++++++++++

# Thats it! Simple!
