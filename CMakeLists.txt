######################################
# ---------------------------------- #
# -------- COMPILER VERSION -------- #
# ---------------------------------- #
######################################

cmake_minimum_required(VERSION 3.0)
project(GetBlunted VERSION 0.0.0)

message(STATUS "CMAKE_CXX_COMPILER_ID: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "CMAKE_CXX_COMPILER_VERSION: ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "CMAKE_SYSTEM: ${CMAKE_SYSTEM}")

# Compiler options.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ------- Memory debug ------- #

# No ASan for static builds :'(
#set(CMAKE_CXX_FLAGS "-fexceptions -fsanitize=address")
#set(ASAN_OPTIONS=check_initialization_order=1)
#set(ASAN_OPTIONS=detect_leaks=1)

# --- Standard compilation --- #

set(CMAKE_CXX_FLAGS "-fexceptions")

# ---------------------------- #

if(dev)
    message(STATUS "--- Using dev build options ---")
    add_definitions(-ggdb3 -O0 -Wall)       # Debugging + No optimization
else()
    add_definitions(-O3 -Wall)              # Much optimization
endif()

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
        src/AdjacencyComponent.cpp
        src/apply_bulk_modifications.cpp
	    src/BicliqueCover.cpp
	    src/Biclique.cpp
	    src/BipartiteGraph.cpp
        src/Bluntifier.cpp
	    src/copy_graph.cpp
        src/Cigar.cpp
        src/Duplicator.cpp
        src/duplicate_terminus.cpp
        src/find_tips.cpp
	    src/GaloisLattice.cpp
        src/gfa_to_handle.cpp
        src/handle_to_gfa.cpp
        src/IncrementalIdMap.cpp
        src/is_single_stranded.cpp
        src/NodeInfo.cpp
        src/OverlapMap.cpp
        src/OverlappingOverlap.cpp
        src/OverlappingOverlapSplicer.cpp
        src/BluntifierAlign.cpp
	    src/ReducedDualGraph.cpp
        src/SubtractiveHandleGraph.cpp
        src/Subgraph.cpp
        src/topological_sort.cpp
        src/traverse.cpp
        src/unchop.cpp
        src/utility.cpp
        src/VertexColoring.cpp
        )


project(GetBlunted)
add_library(GetBlunted STATIC ${SOURCES}) 

# To make sure the library is named Name.so,
# get rid of the "lib" prefix.
set_target_properties(GetBlunted PROPERTIES PREFIX "")

# Eliminate an extraneous -D during compilation.
set_target_properties(GetBlunted PROPERTIES DEFINE_SYMBOL "")


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
    if(EXISTS /usr/local/Cellar/libomp/)
        set(OMPDIR /usr/local/Cellar/libomp/)
    elseif(EXISTS /opt/local/include/libomp)
        set(OMPDIR /opt/local/include/libomp)
    endif()

    set (OpenMP_C "${CMAKE_C_COMPILER}")
    set (OpenMP_C_FLAGS " -Xpreprocessor -fopenmp -I${OMPDIR} -I/usr/local/include -L${OMPDIR} -L/usr/local/lib")
    set (OpenMP_C_LIB_NAMES "libomp" "libgomp" "libiomp5")
    set (OpenMP_CXX "${CMAKE_CXX_COMPILER}")
    set (OpenMP_CXX_FLAGS " -Xpreprocessor -fopenmp -I${OMPDIR} -I/usr/local/include -L${OMPDIR} -L/usr/local/lib")
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

    # To link statically, it is necessary to call lower level commands, as opposed to just
    # target_link_library(target OpenMP::OpenMP_CXX), or else make tries to -lgomp which is apparently
    # bad because it brings in dynamic libraries dlsym, dlopen, dlwtf
    if (OPENMP_FOUND)
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif()

endif()


##############################################
# ------------------------------------------ #
# -------- LINKING EXTERNAL LIBRARY -------- #
# ---------------- gfakluge ---------------- #
# ------------------------------------------ #
##############################################

# Need to explicitly enable ExternalProject functionality
include(ExternalProject)

if (dev)
    # Download or update library as an external project
    ExternalProject_Add(project_gfak
            GIT_REPOSITORY https://github.com/edawson/gfakluge.git
            PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external/
            CONFIGURE_COMMAND ""
            DOWNLOAD_COMMAND ""
            UPDATE_COMMAND ""
            BUILD_IN_SOURCE True
            INSTALL_DIR ${CMAKE_SOURCE_DIR}/external/gfak/
            INSTALL_COMMAND make PREFIX=${CMAKE_SOURCE_DIR}/external/gfak/ install
            )
else()
    # Download or update library as an external project
    ExternalProject_Add(project_gfak
            GIT_REPOSITORY https://github.com/edawson/gfakluge.git
            PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external/
            CONFIGURE_COMMAND ""
            BUILD_IN_SOURCE True
            INSTALL_DIR ${CMAKE_SOURCE_DIR}/external/gfak/
            INSTALL_COMMAND make PREFIX=${CMAKE_SOURCE_DIR}/external/gfak/ install
            )
endif()

# Define INSTALL_DIR as the install directory for external library
ExternalProject_Get_Property(project_gfak INSTALL_DIR)


# THIS SECTION QUARANTINED FOR MAC compatibility :'(
# ---
# Create new library for external project (so it can be linked with main library)
#add_library(gfak STATIC IMPORTED)
#set_property(TARGET gfak
#        PROPERTY IMPORTED_LOCATION ${INSTALL_DIR}/bin/gfak)
#
## Define library as dependent on the downloaded project
#add_dependencies(gfak project_gfak)
# ---


# Define main library as dependent on the downloaded project (transitively)
add_dependencies(GetBlunted project_gfak)

# Ensure that main library has access to primary dependencies' and secondary dependencies' headers
include_directories(external/gfak/include/)

# Add runtime path for main library so it can access external library
#set_property(TARGET GetBlunted PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/bin/gfak")

message(STATUS "INSTALL_DIR: ${INSTALL_DIR}")

##############################################
# ------------------------------------------ #
# -------- LINKING EXTERNAL LIBRARY -------- #
# ------------------ abPOA ----------------- #
# ------------------------------------------ #
##############################################

# Need to explicitly enable ExternalProject functionality
include(ExternalProject)

if(dev)
    # Download or update library as an external project
    ExternalProject_Add(project_abPOA
            GIT_REPOSITORY https://github.com/yangao07/abPOA.git
            GIT_TAG bfe4ac0a4945ed3eadf68282776fc816b299947e
            DOWNLOAD_COMMAND ""
            UPDATE_COMMAND ""
            PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external/
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_SOURCE_DIR}/external/abPOA/ -DCMAKE_INSTALL_LIBDIR=${CMAKE_SOURCE_DIR}/external/abPOA/lib
            BUILD_IN_SOURCE True
            INSTALL_DIR ${CMAKE_SOURCE_DIR}/external/abPOA/
            INSTALL_COMMAND make install
            )
else()
    # Download or update library as an external project
    ExternalProject_Add(project_abPOA
            GIT_REPOSITORY https://github.com/yangao07/abPOA.git
            GIT_TAG bfe4ac0a4945ed3eadf68282776fc816b299947e
            PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external/
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_SOURCE_DIR}/external/abPOA/ -DCMAKE_INSTALL_LIBDIR=${CMAKE_SOURCE_DIR}/external/abPOA/lib
            BUILD_IN_SOURCE True
            INSTALL_DIR ${CMAKE_SOURCE_DIR}/external/abPOA/
            INSTALL_COMMAND make install
            )

endif()

# Define INSTALL_DIR as the install directory for external library
ExternalProject_Get_Property(project_abPOA INSTALL_DIR)

# Create new library for external project (so it can be linked with main library)
add_library(abPOA STATIC IMPORTED)
set_property(TARGET abPOA
        PROPERTY IMPORTED_LOCATION ${INSTALL_DIR}/lib/libabpoa.a)

# Define library as dependent on the downloaded project
add_dependencies(abPOA
        project_abPOA
        )

# Define main library as dependent on the downloaded project (transitively)
add_dependencies(GetBlunted abPOA)

# Add runtime path for main library so it can access external library
#set_property(TARGET GetBlunted PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/bin/abPOA")

include_directories(external/abPOA/include/)

message(STATUS "INSTALL_DIR: ${INSTALL_DIR}")


##############################################
# ------------------------------------------ #
# -------- LINKING EXTERNAL LIBRARY -------- #
# ----------------- kalign ----------------- #
# ------------------------------------------ #
##############################################

# Need to explicitly enable ExternalProject functionality
include(ExternalProject)

if(dev)
#    # Download or update library as an external project
# Download or update library as an external project
    ExternalProject_Add(project_kalign
            GIT_REPOSITORY https://github.com/TimoLassmann/kalign.git
            GIT_TAG 58ca06a51b53d76d3fb96ef335fbc7110c36cd46
            DOWNLOAD_COMMAND ""
            UPDATE_COMMAND ""
            PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external/
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_SOURCE_DIR}/external/kalign/ -DCMAKE_INSTALL_LIBDIR=${CMAKE_SOURCE_DIR}/external/kalign/lib
            BUILD_IN_SOURCE True
            INSTALL_DIR ${CMAKE_SOURCE_DIR}/external/kalign/
            INSTALL_COMMAND make install
            )
else()
    # Download or update library as an external project
    ExternalProject_Add(project_kalign
            GIT_REPOSITORY https://github.com/TimoLassmann/kalign.git
            GIT_TAG 58ca06a51b53d76d3fb96ef335fbc7110c36cd46
            PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external/
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_SOURCE_DIR}/external/kalign/ -DCMAKE_INSTALL_LIBDIR=${CMAKE_SOURCE_DIR}/external/kalign/lib
            BUILD_IN_SOURCE True
            INSTALL_DIR ${CMAKE_SOURCE_DIR}/external/kalign/
            INSTALL_COMMAND make install
            )
endif()

# Define INSTALL_DIR as the install directory for external library
ExternalProject_Get_Property(project_kalign INSTALL_DIR)

# Create new library for external project (so it can be linked with main library)
add_library(kalign SHARED IMPORTED)
set_property(TARGET kalign
        PROPERTY IMPORTED_LOCATION ${INSTALL_DIR}/lib/libkalign.so)

# Define library as dependent on the downloaded project
add_dependencies(kalign
        project_kalign
        )

# Define main library as dependent on the downloaded project (transitively)
add_dependencies(GetBlunted kalign)

# Add runtime path for main library so it can access external library
#set_property(TARGET GetBlunted PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/bin/abPOA")

include_directories(external/kalign/include/)

message(STATUS "INSTALL_DIR: ${INSTALL_DIR}")


##############################################
# ------------------------------------------ #
# -------- LINKING EXTERNAL LIBRARY -------- #
# ------------------ bdsg ------------------ #
# ------------------------------------------ #
##############################################

# Need to explicitly enable ExternalProject functionality
include(ExternalProject)

if(dev)
    # Download or update library as an external project
    ExternalProject_Add(project_bdsg
            GIT_REPOSITORY https://github.com/vgteam/libbdsg-easy.git
            GIT_TAG abe12120bbbbc0127c69de7cc01e925001a0c922
            PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external/
            CONFIGURE_COMMAND ""
            DOWNLOAD_COMMAND ""
            UPDATE_COMMAND ""
            BUILD_IN_SOURCE True
            INSTALL_DIR ${CMAKE_SOURCE_DIR}/external/bdsg/
            INSTALL_COMMAND make INSTALL_PREFIX=${CMAKE_SOURCE_DIR}/external/bdsg/ install
        )
else()
    # Download or update library as an external project
    ExternalProject_Add(project_bdsg
            GIT_REPOSITORY https://github.com/vgteam/libbdsg-easy.git
            GIT_TAG abe12120bbbbc0127c69de7cc01e925001a0c922
            PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external/
            CONFIGURE_COMMAND ""
            BUILD_IN_SOURCE True
            INSTALL_DIR ${CMAKE_SOURCE_DIR}/external/bdsg/
            INSTALL_COMMAND make INSTALL_PREFIX=${CMAKE_SOURCE_DIR}/external/bdsg/ install
    )
endif()

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
add_dependencies(bdsg
        project_bdsg
        libsdsl
        libhandlegraph
        divsufsort)

# Define main library as dependent on the downloaded project (transitively)
add_dependencies(GetBlunted bdsg)

# Ensure that main library has access to primary dependencies' and secondary dependencies' headers
include_directories(external/bdsg/include/
        external/bdsg/include/handlegraph/
        external/bdsg/include/bdsg/)

# Add runtime path for main library so it can access external library
#set_property(TARGET GetBlunted PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/bin/bdsg")

message(STATUS "INSTALL_DIR: ${INSTALL_DIR}")


##############################################
# ------------------------------------------ #
# -------- LINKING EXTERNAL LIBRARY -------- #
# ------------------ spoa ------------------ #
# ------------------------------------------ #
##############################################

# Need to explicitly enable ExternalProject functionality
include(ExternalProject)

if(dev)
    # Download or update library as an external project
    ExternalProject_Add(project_spoa
            GIT_REPOSITORY https://github.com/rvaser/spoa.git
            GIT_TAG 1ab9ab076171e5e4f5fcd4d2a2369f26c7f8f48f
            DOWNLOAD_COMMAND ""
            UPDATE_COMMAND ""
            PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external/
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_SOURCE_DIR}/external/spoa/ -DCMAKE_INSTALL_LIBDIR=${CMAKE_SOURCE_DIR}/external/spoa/lib -DCMAKE_BUILD_TYPE=Release
            BUILD_IN_SOURCE True
            INSTALL_DIR ${CMAKE_SOURCE_DIR}/external/spoa/
            INSTALL_COMMAND make install
    )

else()
    # Download or update library as an external project
    ExternalProject_Add(project_spoa
            GIT_REPOSITORY https://github.com/rvaser/spoa.git
            GIT_TAG 1ab9ab076171e5e4f5fcd4d2a2369f26c7f8f48f
            PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external/
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_SOURCE_DIR}/external/spoa/ -DCMAKE_INSTALL_LIBDIR=${CMAKE_SOURCE_DIR}/external/spoa/lib -DCMAKE_BUILD_TYPE=Release
            BUILD_IN_SOURCE True
            INSTALL_DIR ${CMAKE_SOURCE_DIR}/external/spoa/
            INSTALL_COMMAND make install
    )
endif()

# Define INSTALL_DIR as the install directory for external library
ExternalProject_Get_Property(project_spoa INSTALL_DIR)

# Create new library for external project (so it can be linked with main library)
add_library(spoa STATIC IMPORTED)
set_property(TARGET spoa
        PROPERTY IMPORTED_LOCATION ${INSTALL_DIR}/lib/libspoa.a)

# Define library as dependent on the downloaded project
add_dependencies(spoa
        project_spoa
        )

# Define main library as dependent on the downloaded project (transitively)
add_dependencies(GetBlunted spoa)

# Ensure that main library has access to primary dependencies' and secondary dependencies' headers
include_directories(external/spoa/include/
        ${CMAKE_CURRENT_BINARY_DIR}/external/src/project_spoa/vendor/cereal/include/
        ${CMAKE_CURRENT_BINARY_DIR}/external/src/project_spoa/vendor/cpu_features/include/
        ${CMAKE_CURRENT_BINARY_DIR}/external/src/project_spoa/vendor/simde/simde/
        ${CMAKE_CURRENT_BINARY_DIR}/external/src/project_spoa/vendor/simde/arm/
        ${CMAKE_CURRENT_BINARY_DIR}/external/src/project_spoa/vendor/simde/arm/neon/
        )

# Add runtime path for main library so it can access external library
#set_property(TARGET GetBlunted PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/bin/spoa")

message(STATUS "INSTALL_DIR: ${INSTALL_DIR}")

############################################
# ---------------------------------------- #
# -------- Generating executables -------- #
# ---------------------------------------- #
############################################

# Commented out for static linking!

# --- Recipe for run time linking taken from CMake docs ---#
# https://gitlab.kitware.com/cmake/community/-/wikis/doc/cmake/RPATH-handling#always-full-rpath
# use, i.e. don't skip the full RPATH for the build tree
#set(CMAKE_SKIP_BUILD_RPATH FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
#set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

#set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
#set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
#list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
#if("${isSystemDir}" STREQUAL "-1")
#    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
#endif("${isSystemDir}" STREQUAL "-1")

# ---------------------------------------------------------#




set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

# we can't do static builds in macOS
#if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
#    set(GETBLUNTED_STATIC_LIB_FLAG "")
#else()
#    SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
#    set(GETBLUNTED_STATIC_LIB_FLAG "-static")
#endif()

# kalign breaks static build
set(GETBLUNTED_STATIC_LIB_FLAG "")

find_package(ZLIB REQUIRED)

# -------- TESTS --------

set(TESTS
        test_AdjacencyComponent
        test_bdsg
	    test_BicliqueCover
        test_cigar
        test_divide_handle
        test_duplicate_terminus
        test_gfak
        test_handlegraph
        test_handle_to_gfa
        test_IncrementalIdMap
        test_kalign
        test_map_range_methods
        test_overlaps
        test_spoa
        test_abPOA
        test_utility
        )

foreach(FILENAME_PREFIX ${TESTS})
    add_executable(${FILENAME_PREFIX} src/test/${FILENAME_PREFIX}.cpp)
    target_link_libraries(${FILENAME_PREFIX}
            GetBlunted
            Threads::Threads
            bdsg
            spoa
            divsufsort
            libhandlegraph
            libsdsl
            abPOA
            ZLIB::ZLIB
            kalign
            ${GETBLUNTED_STATIC_LIB_FLAG})

#    if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
#        target_link_libraries(${FILENAME_PREFIX}
#                OpenMP::OpenMP_CXX)
#    endif()
endforeach()


# -------- EXECUTABLES --------

set(EXECUTABLES
        get_blunted
        extract_subgraph
        )

foreach(FILENAME_PREFIX ${EXECUTABLES})
    add_executable(${FILENAME_PREFIX} src/executable/${FILENAME_PREFIX}.cpp)
    target_link_libraries(${FILENAME_PREFIX}
            GetBlunted
            Threads::Threads
            bdsg
            spoa
            divsufsort
            libhandlegraph
            libsdsl
            abPOA
            ZLIB::ZLIB
            kalign
            ${GETBLUNTED_STATIC_LIB_FLAG})

#    if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
#        target_link_libraries(${FILENAME_PREFIX}
#                OpenMP::OpenMP_CXX)
#    endif()
endforeach()

#set_target_properties(get_blunted PROPERTIES LINK_FLAGS "-static" )
#SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -static-libstdc++ -static-libgcc")


# -------- final steps --------

# Where to install
set(BINARY_INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bin)

# The shared library goes to the bin directory.
install(TARGETS
        GetBlunted
        get_blunted
        DESTINATION ${BINARY_INSTALL_DIR})

# Some potentially helpful messages
message(STATUS "CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}")
message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
message(STATUS "CMAKE_CURRENT_BINARY_DIR: ${CMAKE_CURRENT_BINARY_DIR}")
message(STATUS "CMAKE_CURRENT_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")

