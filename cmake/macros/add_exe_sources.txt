macro (add_exe_sources)

    #
    # Get relative file path between the base directory and the current directory
    #
    file (RELATIVE_PATH _relPath "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")


    #
    # Append incoming sources to CMake variable list of sources
    #
    foreach (_src ${ARGN})
        if (_relPath)
            list (APPEND EXE_SRCS "${_relPath}/${_src}")
        else()
            list (APPEND EXE_SRCS "${_src}")
        endif()
    endforeach()


    #
    # If the path is relative, propagate variable containing sources to the parent scope.
    #
    if (_relPath)
        # propagate EXE_SRCS to parent directory
        set (EXE_SRCS ${EXE_SRCS} PARENT_SCOPE)
    endif()


endmacro()
