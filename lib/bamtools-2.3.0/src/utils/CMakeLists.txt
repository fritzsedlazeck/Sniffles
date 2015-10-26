# ==========================
# BamTools CMakeLists.txt
# (c) 2010 Derek Barnett
#
# src/utils/
# ==========================

# list include paths
include_directories( ${BamTools_SOURCE_DIR}/src/api )

# add compiler definitions 
add_definitions( -DBAMTOOLS_UTILS_LIBRARY ) # (for proper exporting of library symbols)
add_definitions( -fPIC ) # (attempt to force PIC compiling on CentOS, not being set on shared libs by CMake)

# create BamTools utils library
add_library( BamTools-utils STATIC
             bamtools_fasta.cpp
             bamtools_options.cpp
             bamtools_pileup_engine.cpp
             bamtools_utilities.cpp
           )

# link BamTools-utils library with BamTools automatically
target_link_libraries( BamTools-utils BamTools )

# set BamTools library properties
set_target_properties( BamTools-utils PROPERTIES
                       OUTPUT_NAME bamtools-utils
                       PREFIX "lib"
                     )
