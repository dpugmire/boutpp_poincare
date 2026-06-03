cmake_minimum_required(VERSION 3.16)
project(codex_cxx2 LANGUAGES CXX)

set(CMAKE_CTEST_ARGUMENTS "--output-on-failure")
include(CTest)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(CODEX_USE_MPI "Enable MPI parallel tracing over particles" OFF)
option(CODEX_USE_VISKORES "Enable Viskores runtime integration" OFF)
option(CODEX_USE_ADIOS "Enable ADIOS2 output" OFF)
set(CODEX_VISKORES_FLOAT_PRECISION "viskores" CACHE STRING "Viskores trace scalar precision: viskores, float32, or float64")
set_property(CACHE CODEX_VISKORES_FLOAT_PRECISION PROPERTY STRINGS viskores float32 float64)

set(NETCDF_ROOT "" CACHE PATH "NetCDF C installation prefix")
if (CODEX_USE_MPI)
    find_package(MPI REQUIRED COMPONENTS CXX)
endif()
if (CODEX_USE_VISKORES)
    find_package(Viskores REQUIRED CONFIG)
    find_package(OpenMP REQUIRED COMPONENTS CXX)
    find_package(Kokkos CONFIG QUIET)
    string(TOLOWER "${CODEX_VISKORES_FLOAT_PRECISION}" CODEX_VISKORES_FLOAT_PRECISION_LOWER)
    if (NOT CODEX_VISKORES_FLOAT_PRECISION_LOWER MATCHES "^(viskores|float32|float64)$")
        message(FATAL_ERROR "CODEX_VISKORES_FLOAT_PRECISION must be one of: viskores, float32, float64")
    endif()
    if (CODEX_VISKORES_FLOAT_PRECISION_LOWER STREQUAL "float32")
        set(CODEX_VISKORES_FLOAT_PRECISION_DEFINE CODEX_VISKORES_FLOAT_PRECISION_FLOAT32=1)
    elseif (CODEX_VISKORES_FLOAT_PRECISION_LOWER STREQUAL "float64")
        set(CODEX_VISKORES_FLOAT_PRECISION_DEFINE CODEX_VISKORES_FLOAT_PRECISION_FLOAT64=1)
    else()
        set(CODEX_VISKORES_FLOAT_PRECISION_DEFINE CODEX_VISKORES_FLOAT_PRECISION_VISKORES=1)
    endif()
    message(STATUS "CodeX Viskores trace precision: ${CODEX_VISKORES_FLOAT_PRECISION_LOWER}")
    if (APPLE AND CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND TARGET viskores::openmp)
        # Some Viskores installs export -fopenmp=libomp, while this local
        # AppleClang build uses the flags detected by FindOpenMP instead.
        # Also avoid linking a second OpenMP runtime beside Viskores' runtime.
        set_target_properties(viskores::openmp PROPERTIES INTERFACE_COMPILE_OPTIONS "")
        set_target_properties(viskores::openmp PROPERTIES INTERFACE_LINK_LIBRARIES "")
    endif()
endif()
if (CODEX_USE_ADIOS)
    find_package(adios2 CONFIG QUIET)
    if (NOT adios2_FOUND)
        find_package(ADIOS2 REQUIRED CONFIG)
    endif()
endif()

find_path(NETCDF_C_INCLUDE_DIR
    NAMES netcdf.h
    HINTS
        ${NETCDF_ROOT}
        $ENV{NETCDF_ROOT}
        /opt/homebrew/include
        /usr/local/include
    PATH_SUFFIXES include)

find_library(NETCDF_C_LIBRARY
    NAMES netcdf
    HINTS
        ${NETCDF_ROOT}
        $ENV{NETCDF_ROOT}
        /opt/homebrew/lib
        /usr/local/lib
    PATH_SUFFIXES lib lib64)

if (NOT NETCDF_C_INCLUDE_DIR)
    message(FATAL_ERROR "Could not find NetCDF C headers (netcdf.h)")
endif()

if (NOT NETCDF_C_LIBRARY)
    message(FATAL_ERROR "Could not find NetCDF C library (libnetcdf)")
endif()

message(STATUS "NetCDF include dir: ${NETCDF_C_INCLUDE_DIR}")
message(STATUS "NetCDF library    : ${NETCDF_C_LIBRARY}")

add_executable(trace_fieldlines
    src/main.cpp
    src/Interpolator.cpp
    src/AparData.cpp
    src/AparFieldModel.cpp
    src/FieldLineIntegrator.cpp
    src/TracePostProcessor.cpp
    src/ViskoresFieldLineTracer.cpp
    src/PunctureDetector.cpp
    src/PoincareOutput.cpp
    src/MatlabComparator.cpp
    src/ValidationSuite.cpp)

target_include_directories(trace_fieldlines
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${NETCDF_C_INCLUDE_DIR}
)

target_link_libraries(trace_fieldlines
    PRIVATE
        ${NETCDF_C_LIBRARY}
)

if (CODEX_USE_MPI)
    target_compile_definitions(trace_fieldlines PRIVATE CODEX_USE_MPI=1)
    target_link_libraries(trace_fieldlines PRIVATE MPI::MPI_CXX)
endif()

if (CODEX_USE_VISKORES)
    target_compile_definitions(trace_fieldlines PRIVATE
        CODEX_USE_VISKORES=1
        ${CODEX_VISKORES_FLOAT_PRECISION_DEFINE})
    target_compile_options(trace_fieldlines PRIVATE $<TARGET_PROPERTY:OpenMP::OpenMP_CXX,INTERFACE_COMPILE_OPTIONS>)
    target_compile_definitions(trace_fieldlines PRIVATE $<TARGET_PROPERTY:OpenMP::OpenMP_CXX,INTERFACE_COMPILE_DEFINITIONS>)
    target_include_directories(trace_fieldlines PRIVATE $<TARGET_PROPERTY:OpenMP::OpenMP_CXX,INTERFACE_INCLUDE_DIRECTORIES>)
    # Avoid loading two different OpenMP runtimes in one process.
    # Viskores already links against the OpenMP runtime from its own build.
    target_link_libraries(trace_fieldlines PRIVATE viskores::cont viskores::worklet)
    if (TARGET Kokkos::kokkos)
        target_link_libraries(trace_fieldlines PRIVATE Kokkos::kokkos)
    endif()
    viskores_add_target_information(trace_fieldlines
        DROP_UNUSED_SYMBOLS
        MODIFY_CUDA_FLAGS
        DEVICE_SOURCES
            src/ViskoresFieldLineTracer.cpp)
endif()

if (CODEX_USE_ADIOS)
    target_sources(trace_fieldlines PRIVATE src/AdiosPoincareOutput.cpp)
    target_compile_definitions(trace_fieldlines PRIVATE CODEX_USE_ADIOS=1)
    if (CODEX_USE_MPI)
        if (TARGET adios2::cxx_mpi)
            set(CODEX_ADIOS_TARGET adios2::cxx_mpi)
        elseif (TARGET adios2::adios2_mpi)
            set(CODEX_ADIOS_TARGET adios2::adios2_mpi)
        elseif (TARGET ADIOS2::adios2_mpi)
            set(CODEX_ADIOS_TARGET ADIOS2::adios2_mpi)
        else()
            message(FATAL_ERROR "CODEX_USE_ADIOS=ON with CODEX_USE_MPI=ON requires an MPI-enabled ADIOS2 C++ target")
        endif()
    else()
        if (TARGET adios2::cxx)
            set(CODEX_ADIOS_TARGET adios2::cxx)
        elseif (TARGET adios2::adios2)
            set(CODEX_ADIOS_TARGET adios2::adios2)
        elseif (TARGET ADIOS2::adios2)
            set(CODEX_ADIOS_TARGET ADIOS2::adios2)
        else()
            message(FATAL_ERROR "CODEX_USE_ADIOS=ON requires an ADIOS2 C++ target")
        endif()
    endif()
    message(STATUS "ADIOS2 output target: ${CODEX_ADIOS_TARGET}")
    target_link_libraries(trace_fieldlines PRIVATE ${CODEX_ADIOS_TARGET})
endif()

if (BUILD_TESTING AND CODEX_USE_ADIOS)
    add_executable(test_adios_poincare_reader
        tests/adios_validate_poincare.cpp)
    target_link_libraries(test_adios_poincare_reader
        PRIVATE
            ${CODEX_ADIOS_TARGET})
    if (CODEX_USE_MPI)
        target_compile_definitions(test_adios_poincare_reader
            PRIVATE
                CODEX_USE_MPI=1)
        target_link_libraries(test_adios_poincare_reader
            PRIVATE
                MPI::MPI_CXX)
    endif()
endif()

if (BUILD_TESTING)
    add_executable(test_periodic_z_spline
        tests/test_periodic_z_spline.cpp
        src/AparFieldModel.cpp
        src/AparData.cpp
        src/Interpolator.cpp)

    target_include_directories(test_periodic_z_spline
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/include
            ${NETCDF_C_INCLUDE_DIR}
    )

    target_link_libraries(test_periodic_z_spline
        PRIVATE
            ${NETCDF_C_LIBRARY}
    )

    add_test(NAME periodic_z_spline
        COMMAND $<TARGET_FILE:test_periodic_z_spline>)
endif()

if (BUILD_TESTING AND CODEX_USE_VISKORES)
    find_package(Python3 REQUIRED COMPONENTS Interpreter)

    set(CODEX_TEST_DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/data" CACHE PATH "Directory containing apar NetCDF files and MATLAB baselines")
    set(CODEX_TEST_APAR_SINGLE "${CODEX_TEST_DATA_DIR}/apar.single.nc" CACHE FILEPATH "Single-null apar NetCDF file for tests")
    set(CODEX_TEST_APAR_CIRC "${CODEX_TEST_DATA_DIR}/apar.circ.nc" CACHE FILEPATH "Circular apar NetCDF file for tests")
    set(CODEX_TEST_REFERENCE_DIR "${CODEX_TEST_DATA_DIR}/zperiod_5" CACHE PATH "MATLAB reference directory for tests")
    set(CODEX_TEST_VISKORES_DEVICE "serial" CACHE STRING "Viskores device used by CTest trace tests")
    set(CODEX_TEST_NP_MAX "700" CACHE STRING "Maximum punctures per field line used by CTest trace tests")
    set(CODEX_TEST_TOLERANCE "1e-8" CACHE STRING "Maximum absolute XYZ tolerance for MATLAB comparison tests")
    set(CODEX_TEST_MATLAB_PARTIAL_TOLERANCE "1e-2" CACHE STRING "Maximum per-row XYZ tolerance for partial MATLAB comparison tests")
    set(CODEX_TEST_MPI_PROCS "2" CACHE STRING "Number of MPI ranks used by MPI CTest trace tests")

    set(CODEX_TEST_SINGLE_LINES "11,21,26,31,41,51,61,81,101,121,126,131,141,147,148,149,150,151,152,153,161,171,176,181,191,201")
    set(CODEX_TEST_CIRC_LINES "51,101,151,201,301,401,501")
    set(CODEX_TEST_SINGLE_MATLAB_PREFIX_LINES "26,31")
    set(CODEX_TEST_CIRC_MATLAB_PREFIX_LINES "51,101,151")
    set(CODEX_TEST_ADIOS_SINGLE_LINES "26" CACHE STRING "Single-null line used by ADIOS schema smoke tests")
    set(CODEX_TEST_ADIOS_CIRC_LINES "51" CACHE STRING "Circular line used by ADIOS schema smoke tests")
    set(CODEX_TEST_ADIOS_SINGLE_PSI_AXIS "-0.2754453" CACHE STRING "Single-null psi_axis value used by ADIOS schema smoke tests")
    set(CODEX_TEST_ADIOS_SINGLE_PSI_BNDRY "-0.1304577" CACHE STRING "Single-null psi_bndry value used by ADIOS schema smoke tests")
    set(CODEX_TEST_ADIOS_CIRC_PSI_AXIS "-1.5170289278030396" CACHE STRING "Circular psi_axis value used by ADIOS schema smoke tests")
    set(CODEX_TEST_ADIOS_CIRC_PSI_BNDRY "-0.5056424736976624" CACHE STRING "Circular psi_bndry value used by ADIOS schema smoke tests")

    function(codex_add_vp_test test_name test_mode divertor lines mpi_procs)
        set(test_args
            --mode ${test_mode}
            --exe $<TARGET_FILE:trace_fieldlines>
            --divertor ${divertor}
            --lines ${lines}
            --single-apar ${CODEX_TEST_APAR_SINGLE}
            --circ-apar ${CODEX_TEST_APAR_CIRC}
            --reference-dir ${CODEX_TEST_REFERENCE_DIR}
            --output-dir ${CMAKE_CURRENT_BINARY_DIR}/test-output/${test_name}
            --viskores-device ${CODEX_TEST_VISKORES_DEVICE}
            --np-max ${CODEX_TEST_NP_MAX}
            --tolerance ${CODEX_TEST_TOLERANCE})

        list(APPEND test_args ${ARGN})

        if (${mpi_procs} GREATER 1)
            list(APPEND test_args
                --mpi-exec ${MPIEXEC_EXECUTABLE}
                --mpi-numproc-flag=${MPIEXEC_NUMPROC_FLAG}
                --mpi-procs ${mpi_procs})
        endif()

        add_test(NAME ${test_name}
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_trace_test.py
                ${test_args})
        set_tests_properties(${test_name} PROPERTIES
            SKIP_RETURN_CODE 77
            TIMEOUT 180)
    endfunction()

    codex_add_vp_test(viskores_generate_single_text generate single ${CODEX_TEST_SINGLE_LINES} 1)
    codex_add_vp_test(viskores_generate_circ_text generate circ ${CODEX_TEST_CIRC_LINES} 1)
    codex_add_vp_test(viskores_compare_single_matlab_prefix compare single ${CODEX_TEST_SINGLE_MATLAB_PREFIX_LINES} 1
        --compare-policy leading-prefix
        --row-tolerance ${CODEX_TEST_MATLAB_PARTIAL_TOLERANCE}
        --min-leading-rows 12)
    codex_add_vp_test(viskores_compare_circ_matlab_prefix compare circ ${CODEX_TEST_CIRC_MATLAB_PREFIX_LINES} 1
        --compare-policy leading-prefix
        --row-tolerance ${CODEX_TEST_MATLAB_PARTIAL_TOLERANCE}
        --min-leading-rows 20)
    codex_add_vp_test(viskores_compare_single_matlab_exact_all compare single ${CODEX_TEST_SINGLE_LINES} 1
        --compare-policy exact)
    codex_add_vp_test(viskores_compare_circ_matlab_exact_all compare circ ${CODEX_TEST_CIRC_LINES} 1
        --compare-policy exact)

    if (CODEX_USE_ADIOS)
        codex_add_vp_test(viskores_generate_single_adios generate single ${CODEX_TEST_SINGLE_LINES} 1
            --output-format both
            --output ${CMAKE_CURRENT_BINARY_DIR}/test-output/viskores_generate_single_adios/poincare.bp)
        codex_add_vp_test(viskores_generate_circ_adios generate circ ${CODEX_TEST_CIRC_LINES} 1
            --output-format both
            --output ${CMAKE_CURRENT_BINARY_DIR}/test-output/viskores_generate_circ_adios/poincare.bp)
        codex_add_vp_test(viskores_adios_smoke_single adios-smoke single ${CODEX_TEST_ADIOS_SINGLE_LINES} 1
            --output-format adios
            --output ${CMAKE_CURRENT_BINARY_DIR}/test-output/viskores_adios_smoke_single/poincare.bp
            --psi-axis ${CODEX_TEST_ADIOS_SINGLE_PSI_AXIS}
            --psi-bndry ${CODEX_TEST_ADIOS_SINGLE_PSI_BNDRY}
            --adios-validator $<TARGET_FILE:test_adios_poincare_reader>)
        codex_add_vp_test(viskores_adios_smoke_circ adios-smoke circ ${CODEX_TEST_ADIOS_CIRC_LINES} 1
            --output-format adios
            --output ${CMAKE_CURRENT_BINARY_DIR}/test-output/viskores_adios_smoke_circ/poincare.bp
            --psi-axis ${CODEX_TEST_ADIOS_CIRC_PSI_AXIS}
            --psi-bndry ${CODEX_TEST_ADIOS_CIRC_PSI_BNDRY}
            --adios-validator $<TARGET_FILE:test_adios_poincare_reader>)
    endif()

    if (CODEX_USE_MPI AND MPIEXEC_EXECUTABLE)
        codex_add_vp_test(viskores_mpi_generate_single_text generate single ${CODEX_TEST_SINGLE_LINES} ${CODEX_TEST_MPI_PROCS})
        codex_add_vp_test(viskores_mpi_generate_circ_text generate circ ${CODEX_TEST_CIRC_LINES} ${CODEX_TEST_MPI_PROCS})
        codex_add_vp_test(viskores_mpi_compare_single_serial compare-serial single ${CODEX_TEST_SINGLE_LINES} ${CODEX_TEST_MPI_PROCS})
        codex_add_vp_test(viskores_mpi_compare_circ_serial compare-serial circ ${CODEX_TEST_CIRC_LINES} ${CODEX_TEST_MPI_PROCS})
        codex_add_vp_test(viskores_mpi_compare_single_matlab_prefix compare single ${CODEX_TEST_SINGLE_MATLAB_PREFIX_LINES} ${CODEX_TEST_MPI_PROCS}
            --compare-policy leading-prefix
            --row-tolerance ${CODEX_TEST_MATLAB_PARTIAL_TOLERANCE}
            --min-leading-rows 12)
        codex_add_vp_test(viskores_mpi_compare_circ_matlab_prefix compare circ ${CODEX_TEST_CIRC_MATLAB_PREFIX_LINES} ${CODEX_TEST_MPI_PROCS}
            --compare-policy leading-prefix
            --row-tolerance ${CODEX_TEST_MATLAB_PARTIAL_TOLERANCE}
            --min-leading-rows 20)
        if (CODEX_USE_ADIOS)
            codex_add_vp_test(viskores_mpi_generate_single_adios generate single ${CODEX_TEST_SINGLE_LINES} ${CODEX_TEST_MPI_PROCS}
                --output-format both
                --output ${CMAKE_CURRENT_BINARY_DIR}/test-output/viskores_mpi_generate_single_adios/poincare.bp)
            codex_add_vp_test(viskores_mpi_generate_circ_adios generate circ ${CODEX_TEST_CIRC_LINES} ${CODEX_TEST_MPI_PROCS}
                --output-format both
                --output ${CMAKE_CURRENT_BINARY_DIR}/test-output/viskores_mpi_generate_circ_adios/poincare.bp)
        endif()
    endif()
endif()
