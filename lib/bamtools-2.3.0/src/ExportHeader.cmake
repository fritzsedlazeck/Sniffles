#
# ExportHeader
#

function( ExportHeader MODULE FILE DEST )

    # if haven't defined our custom 'build target'
    # not exactly a build target, but lets this command get
    # checked any time build step happens
    if( NOT TARGET ${MODULE} )
        add_custom_target( ${MODULE} ALL COMMENT "Exporting ${MODULE}" )
    endif( NOT TARGET ${MODULE} )

    # get the filename (without path)
    get_filename_component( FILENAME "${FILE}" NAME )

    # copy header to destination
    add_custom_command( TARGET ${MODULE} COMMAND
        ${CMAKE_COMMAND} -E copy_if_different
        "${CMAKE_CURRENT_SOURCE_DIR}/${FILE}"
        "${CMAKE_SOURCE_DIR}/include/${DEST}/${FILENAME}" )

    # make sure files are properly 'installed'
    install( FILES "${FILE}" DESTINATION "include/bamtools/${DEST}" )

endfunction( ExportHeader )

