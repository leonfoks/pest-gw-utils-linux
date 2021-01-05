if ("$ENV{MKL_ROOT}" STREQUAL "")
    message(FATAL_ERROR "Could not find MKL install directory.  Please create an environment variable MKL_ROOT to the install directory.")
else()
    set(SCI_ROOT "$ENV{MKL_ROOT}" CACHE INTERNAL "Get the SCI install directory")
    message("Using MKL installed at ${MKL_ROOT}")

    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -I${MKL_ROOT}/include")

    if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${MKL_ROOT}/lib/intel64 -lsci_gnu")
    elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${MKL_ROOT}/lib/intel64 -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl")
    endif()
        
endif()
