######################################
# ---------------------------------- #
# -------- COMPILER VERSION -------- #
# ---------------------------------- #
######################################

cmake_minimum_required(VERSION 3.0)
project(Bluntifier VERSION 0.0.0)

# Compiler options.
#set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++17")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#set(CMAKE_CXX_FLAGS "-fexceptions -fno-new-ttp-matching -fsanitize=address")
set(CMAKE_CXX_FLAGS "-fexceptions -fno-new-ttp-matching")

#add_definitions(-ggdb3 -O0 -Wall)       # Debugging + No optimization
add_definitions(-O3 -Wall)              # Much optimization

#########################################
# ------------------------------------- #
# -------- SOURCES AND HEADERS -------- #
# ------------------------------------- #
#########################################

# Include header files
include_directories(
        "inc"
)

# Define our shared library sources. NOT test/executables.
set(SOURCES
        src/a.cpp
#        src/gfa_to_handle.cpp
        )


project(Bluntifier)
add_library(Bluntifier SHARED ${SOURCES})

# To make sure the library is named Name.so,
# get rid of the "lib" prefix.
set_target_properties(Bluntifier PROPERTIES PREFIX "")

# Eliminate an extraneous -D during compilation.
set_target_properties(Bluntifier PROPERTIES DEFINE_SYMBOL "")


##############################################
# ------------------------------------------ #
# --------- LINKING NATIVE LIBRARY --------- #
# ------------------ OMP ------------------- #
# ------------------------------------------ #
##############################################

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

    # assumes clang build
    # we can't reliably detect when we're using clang, so for the time being we assume
    # TODO: can't we though?

    # adapted from https://stackoverflow.com/questions/46414660/macos-cmake-and-openmp
    # find_package(OpenMP) does not work reliably on macOS, so we do its work ourselves
    set (OpenMP_C "${CMAKE_C_COMPILER}")
    set (OpenMP_C_FLAGS " -Xpreprocessor -fopenmp -I/opt/local/include/libomp -I/usr/local/include -L/opt/local/lib/libomp -L/usr/local/lib")
    set (OpenMP_C_LIB_NAMES "libomp" "libgomp" "libiomp5")
    set (OpenMP_CXX "${CMAKE_CXX_COMPILER}")
    set (OpenMP_CXX_FLAGS " -Xpreprocessor -fopenmp -I/opt/local/include/libomp -I/usr/local/include -L/opt/local/lib/libomp -L/usr/local/lib")
    set (OpenMP_CXX_LIB_NAMES "libomp" "libgomp" "libiomp5")
    set (OpenMP_libomp_LIBRARY "omp")
    set (OpenMP_libgomp_LIBRARY "gomp")
    set (OpenMP_libiomp5_LIBRARY "iomp5")

    # and now add the OpenMP parameters to the compile flags
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS} -lomp")

    # Mac needs libdl and libomp when linking the library
    set(PLATFORM_EXTRA_LIB_FLAGS -ldl -lomp)

elseif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")

    find_package(OpenMP REQUIRED)

#    # add the flags it detects to the compile flags
#    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -fopenmp")
#    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -fopenmp")
#    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#
#    # Linux only needs libdl when linking the library
#    set(PLATFORM_EXTRA_LIB_FLAGS -ldl)

endif()


##############################################
# ------------------------------------------ #
# -------- LINKING EXTERNAL LIBRARY -------- #
# ---------------- gfakluge ---------------- #
# ------------------------------------------ #
##############################################

# Need to explicitly enable ExternalProject functionality
include(ExternalProject)

# Download or update library as an external project
ExternalProject_Add(project_gfak
        GIT_REPOSITORY https://github.com/edawson/gfakluge.git
        PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external/
        CONFIGURE_COMMAND ""
        BUILD_IN_SOURCE True
        INSTALL_DIR ${CMAKE_SOURCE_DIR}/external/gfak/
        INSTALL_COMMAND make PREFIX=${CMAKE_SOURCE_DIR}/external/gfak/ install
        )

# Define INSTALL_DIR as the install directory for external library
ExternalProject_Get_Property(project_gfak INSTALL_DIR)

# Create new library for external project (so it can be linked with main library)
add_library(gfak STATIC IMPORTED)
set_property(TARGET gfak
        PROPERTY IMPORTED_LOCATION ${INSTALL_DIR}/bin/gfak)

# Define library as dependent on the downloaded project
add_dependencies(gfak project_gfak)

# Define main library as dependent on the downloaded project (transitively)
add_dependencies(Bluntifier gfak)

# Ensure that main library has access to primary dependencies' and secondary dependencies' headers
include_directories(external/gfak/include/)

# Add runtime path for main library so it can access external library
#set_property(TARGET Bluntifier PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/bin/gfak")

message(STATUS "INSTALL_DIR: ${INSTALL_DIR}")


##############################################
# ------------------------------------------ #
# -------- LINKING EXTERNAL LIBRARY -------- #
# ------------------ bdsg ------------------ #
# ------------------------------------------ #
##############################################

# Need to explicitly enable ExternalProject functionality
include(ExternalProject)

# Download or update library as an external project
ExternalProject_Add(project_bdsg
        GIT_REPOSITORY https://github.com/vgteam/libbdsg-easy.git
        PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external/
        CONFIGURE_COMMAND ""
#        UPDATE_COMMAND ""
        BUILD_IN_SOURCE True
        INSTALL_DIR ${CMAKE_SOURCE_DIR}/external/bdsg/
        INSTALL_COMMAND make INSTALL_PREFIX=${CMAKE_SOURCE_DIR}/external/bdsg/ install
        )

# Define INSTALL_DIR as the install directory for external library
ExternalProject_Get_Property(project_bdsg INSTALL_DIR)

# Create new library for external project (so it can be linked with main library)
add_library(bdsg STATIC IMPORTED)
set_property(TARGET bdsg
        PROPERTY IMPORTED_LOCATION ${INSTALL_DIR}/lib/libbdsg.a)

add_library(divsufsort STATIC IMPORTED)
set_property(TARGET divsufsort
        PROPERTY IMPORTED_LOCATION ${INSTALL_DIR}/lib/libdivsufsort64.a)

add_library(libhandlegraph STATIC IMPORTED)
set_property(TARGET libhandlegraph
        PROPERTY IMPORTED_LOCATION ${INSTALL_DIR}/lib/libhandlegraph.a)

add_library(libsdsl STATIC IMPORTED)
set_property(TARGET libsdsl
        PROPERTY IMPORTED_LOCATION ${INSTALL_DIR}/lib/libsdsl.a)

# Define library as dependent on the downloaded project
add_dependencies(bdsg project_bdsg)

# Define main library as dependent on the downloaded project (transitively)
add_dependencies(Bluntifier bdsg)

# Ensure that main library has access to primary dependencies' and secondary dependencies' headers
include_directories(external/bdsg/include/)

# Add runtime path for main library so it can access external library
#set_property(TARGET Bluntifier PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/bin/bdsg")

message(STATUS "INSTALL_DIR: ${INSTALL_DIR}")


############################################
# ---------------------------------------- #
# -------- Generating executables -------- #
# ---------------------------------------- #
############################################

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

# stupid "experimental" filesystem library for gcc <8.0
link_libraries(stdc++fs)


# -------- TESTS --------

set(FILENAME_PREFIX test_gfak)
add_executable(${FILENAME_PREFIX} src/test/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN")
target_link_libraries(${FILENAME_PREFIX}
        Bluntifier
        Threads::Threads
        stdc++fs
        gfak
        bdsg
        divsufsort
        libhandlegraph
        libsdsl
        OpenMP::OpenMP_CXX)

set(FILENAME_PREFIX test_bdsg)
add_executable(${FILENAME_PREFIX} src/test/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN")
target_link_libraries(${FILENAME_PREFIX}
        Bluntifier
        Threads::Threads
        stdc++fs
        gfak
        bdsg
        divsufsort
        libhandlegraph
        libsdsl
        OpenMP::OpenMP_CXX)

set(FILENAME_PREFIX test_handlegraph)
add_executable(${FILENAME_PREFIX} src/test/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN")
target_link_libraries(${FILENAME_PREFIX}
        Bluntifier
        Threads::Threads
        stdc++fs
        gfak
        bdsg
        divsufsort
        libhandlegraph
        libsdsl
        OpenMP::OpenMP_CXX)

# -------- EXECUTABLES --------

#set(FILENAME_PREFIX a)
#add_executable(${FILENAME_PREFIX} src/executables/${FILENAME_PREFIX}.cpp)
#set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN")
#target_link_libraries(${FILENAME_PREFIX} Bluntifier Threads::Threads ${Boost_LIBRARIES} stdc++fs)

# -------- final steps --------

# Where to install
set(BINARY_INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bin)

# The shared library goes to the bin directory.
install(TARGETS
        Bluntifier
        DESTINATION ${BINARY_INSTALL_DIR})

# Some potentially helpful messages
message(STATUS "CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}")
message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
message(STATUS "CMAKE_CURRENT_BINARY_DIR: ${CMAKE_CURRENT_BINARY_DIR}")
message(STATUS "CMAKE_CURRENT_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")

