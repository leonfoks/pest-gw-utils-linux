cmake_minimum_required(VERSION 3.8)


file(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/*.f90")
foreach(file ${files})
    get_filename_component(prog ${file} NAME_WE)

    # Create the testing suite
    project(${prog} LANGUAGES Fortran)

    add_executable(${PROJECT_NAME} ${file})

    
    # Link coretran to the test suite
    target_link_libraries(${PROJECT_NAME} gwlib)

endforeach()


file(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/*.f")
foreach(file ${files})
    get_filename_component(prog ${file} NAME_WE)

    # Create the testing suite
    project(${prog} LANGUAGES Fortran)

    add_executable(${PROJECT_NAME} ${file})


    # Link coretran to the test suite
    target_link_libraries(${PROJECT_NAME} gwlib)

endforeach()
