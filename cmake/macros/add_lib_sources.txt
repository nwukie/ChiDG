macro (add_lib_sources)

    #
    # Get relative file path between the base directory and the current directory
    #
    file (RELATIVE_PATH _relPath "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")


    #
    # Append incoming sources to CMake variable list of sources
    #
    foreach (_src ${ARGN})
        if (_relPath)
            list (APPEND LIB_SRCS "${_relPath}/${_src}")
        else()
            list (APPEND LIB_SRCS "${_src}")
        endif()
    endforeach()


    #
    # If the path is relative, propagate variable containing sources to the parent scope.
    #
    if (_relPath)
        # propagate LIB_SRCS to parent directory
        set (LIB_SRCS ${LIB_SRCS} PARENT_SCOPE)
    endif()

endmacro()
