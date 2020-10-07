if ("$ENV{LIBSCI_ROOT}" STREQUAL "")
    message(FATAL_ERROR "Could not find LIBSCI install directory.  Please create an environment variable LIBSCI_ROOT to the install directory.")
else()
    set(LIBSCI_ROOT "$ENV{LIBSCI_ROOT}" CACHE INTERNAL "Get the LIBSCI install directory")
    message("Using LIBSCI installed at ${LIBSCI_ROOT}")

    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -I${LIBSCI_ROOT}/include")

    if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${LIBSCI_ROOT}/lib -lsci_gnu")
    elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${LIBSCI_ROOT}/lib -lsci_intel")
    endif()
        
endif()
