add_executable (nucleusParametersTuner.e
    main.cpp 
    Nucleus.cpp 
    ParameterReader.cpp 
    MakeDensity.cpp
    gauss_quadrature.cpp
    )

target_link_libraries (nucleusParametersTuner.e ${LIBS})

install (TARGETS nucleusParametersTuner.e DESTINATION ${CMAKE_HOME_DIRECTORY})
