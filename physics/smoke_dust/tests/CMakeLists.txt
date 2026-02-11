cmake_minimum_required(VERSION 3.10)
project(smoke_dust_tests LANGUAGES Fortran)

# Find OpenMP
find_package(OpenMP)

# Paths to source files
set(HOOKS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../hooks)
set(SMOKE_DUST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/..)

# FENGSHA Test
set(FENGSHA_SOURCES
    ${HOOKS_DIR}/machine.F
    ${SMOKE_DUST_DIR}/rrfs_smoke_config.F90
    ${SMOKE_DUST_DIR}/dust_data_mod.F90
    ${SMOKE_DUST_DIR}/dust_fengsha_mod.F90
    test_dust_fengsha.F90
)
add_executable(test_dust_fengsha ${FENGSHA_SOURCES})

# GOCART Sea Salt Test
set(SEAS_SOURCES
    ${HOOKS_DIR}/machine.F
    ${SMOKE_DUST_DIR}/seas_data_mod.F90
    ${SMOKE_DUST_DIR}/seas_ngac_mod.F90
    ${SMOKE_DUST_DIR}/seas_mod.F90
    test_seas_gocart.F90
)
add_executable(test_seas_gocart ${SEAS_SOURCES})

# Large-scale Wet Deposition Test
set(WETDEP_SOURCES
    ${HOOKS_DIR}/machine.F
    ${SMOKE_DUST_DIR}/rrfs_smoke_config.F90
    ${SMOKE_DUST_DIR}/module_wetdep_ls.F90
    test_wetdep_ls.F90
)
add_executable(test_wetdep_ls ${WETDEP_SOURCES})

# Compile flags and OpenMP linking
foreach(target test_dust_fengsha test_seas_gocart test_wetdep_ls)
    set_target_properties(${target} PROPERTIES
        Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules
    )
    if(OpenMP_Fortran_FOUND)
        target_link_libraries(${target} PUBLIC OpenMP::OpenMP_Fortran)
    endif()
endforeach()

enable_testing()
add_test(NAME run_fengsha_test COMMAND test_dust_fengsha)
add_test(NAME run_seas_test COMMAND test_seas_gocart)
add_test(NAME run_wetdep_test COMMAND test_wetdep_ls)
