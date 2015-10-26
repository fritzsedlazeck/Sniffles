# ==========================
# BamTools CMakeLists.txt
# (c) 2010 Derek Barnett
#
# src/toolkit
# ==========================

# set include path
include_directories( ${BamTools_SOURCE_DIR}/src/api
                     ${BamTools_SOURCE_DIR}/src/utils
                     ${BamTools_SOURCE_DIR}/src/third_party
                   )

# compile main bamtools application
add_executable( bamtools_cmd
                bamtools_convert.cpp
                bamtools_count.cpp
                bamtools_coverage.cpp
                bamtools_filter.cpp
                bamtools_header.cpp
                bamtools_index.cpp
                bamtools_merge.cpp
                bamtools_random.cpp
                bamtools_resolve.cpp
                bamtools_revert.cpp
                bamtools_sort.cpp
                bamtools_split.cpp
                bamtools_stats.cpp
                bamtools.cpp
              )

# set BamTools application properties
set_target_properties( bamtools_cmd PROPERTIES
                       VERSION  2.3.0
                       OUTPUT_NAME "bamtools"
                     )
# make version info available in application
configure_file( bamtools_version.h.in ${BamTools_SOURCE_DIR}/src/toolkit/bamtools_version.h )

# define libraries to link
target_link_libraries( bamtools_cmd BamTools BamTools-utils jsoncpp )

# set application install destinations
install( TARGETS bamtools_cmd DESTINATION "bin")
