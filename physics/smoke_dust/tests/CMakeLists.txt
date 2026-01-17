cmake_minimum_required(VERSION 3.10)
project(fengsha_test LANGUAGES Fortran)

# Find OpenMP
find_package(OpenMP)

# Paths to source files
set(HOOKS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../hooks)
set(SMOKE_DUST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/..)

# Source files needed for the test
set(SOURCES
    ${HOOKS_DIR}/machine.F
    ${SMOKE_DUST_DIR}/rrfs_smoke_config.F90
    ${SMOKE_DUST_DIR}/dust_data_mod.F90
    ${SMOKE_DUST_DIR}/dust_fengsha_mod.F90
    test_dust_fengsha.F90
)

add_executable(test_dust_fengsha ${SOURCES})

# Compile flags
set_target_properties(test_dust_fengsha PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules
)

if(OpenMP_Fortran_FOUND)
    target_link_libraries(test_dust_fengsha PUBLIC OpenMP::OpenMP_Fortran)
endif()

enable_testing()
add_test(NAME run_fengsha_test COMMAND test_dust_fengsha)
