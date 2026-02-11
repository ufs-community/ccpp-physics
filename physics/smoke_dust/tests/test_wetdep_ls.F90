!> \file test_wetdep_ls.F90
!! Comprehensive unit test for the large-scale wet deposition scheme.

program test_wetdep_ls
  use machine, only : kind_phys
  use module_wetdep_ls, only : wetdep_ls
  use rrfs_smoke_config, only : p_smoke, p_dust_1, p_coarse_pm, p_qc
  implicit none

  ! Parameters
  integer, parameter :: im = 1, jm = 1, km = 5
  integer, parameter :: num_chem = 20
  integer, parameter :: num_moist = 5
  integer, parameter :: ndvel = 1

  ! Global data for tests
  real(kind_phys) :: dt
  real(kind_phys), dimension(im, km, jm, num_moist) :: moist
  real(kind_phys), dimension(im, km, jm) :: rho, dz8w, vvel
  real(kind_phys), dimension(im, km, jm, num_chem) :: var
  real(kind_phys), dimension(im, jm) :: rain, wetdpr_smoke, wetdpr_dust, wetdpr_coarsepm

  integer :: total_passed, total_failed

  total_passed = 0
  total_failed = 0

  dt = 600.0_kind_phys

  print *, "Starting Large-scale Wet Deposition comprehensive unit test suite..."

  ! Case 1: Standard rain (should have wet deposition)
  call reset_inputs()
  rain = 1.0e-3_kind_phys ! Some rain
  moist(:, :, :, p_qc) = 1.0e-4_kind_phys ! some cloud water
  call run_case("Standard rain wet deposition", .true.)

  ! Case 2: No rain (should have no wet deposition)
  call reset_inputs()
  rain = 0.0_kind_phys
  moist(:, :, :, p_qc) = 1.0e-4_kind_phys
  call run_case("No rain", .false.)

  ! Case 3: No cloud water (should have no wet deposition)
  call reset_inputs()
  rain = 1.0e-3_kind_phys
  moist(:, :, :, p_qc) = 0.0_kind_phys
  call run_case("No cloud water", .false.)

  print *, "---------------------------------------"
  print *, "Test Summary: ", total_passed, " PASSED, ", total_failed, " FAILED"
  print *, "---------------------------------------"

  if (total_failed > 0) stop 1

contains

  subroutine reset_inputs()
    var = 1.0e-6_kind_phys
    moist = 0.0_kind_phys
    rho = 1.0_kind_phys
    dz8w = 100.0_kind_phys
    vvel = 1.0_kind_phys ! vertical velocity
    rain = 0.0_kind_phys
    wetdpr_smoke = 0.0_kind_phys
    wetdpr_dust = 0.0_kind_phys
    wetdpr_coarsepm = 0.0_kind_phys
  end subroutine reset_inputs

  subroutine run_case(name, expect_non_zero)
    character(len=*), intent(in) :: name
    logical, intent(in) :: expect_non_zero
    logical :: is_non_zero
    real(kind_phys) :: total_dep

    call wetdep_ls(dt, var, rain, moist,                       &
                   rho, num_chem, num_moist, ndvel, dz8w, vvel,&
                   wetdpr_smoke, wetdpr_dust, wetdpr_coarsepm, &
                   1, im, 1, jm, 1, km,                        &
                   1, im, 1, jm, 1, km,                        &
                   1, im, 1, jm, 1, km                         )

    total_dep = sum(wetdpr_smoke) + sum(wetdpr_dust) + sum(wetdpr_coarsepm)
    is_non_zero = total_dep > 1.0e-20_kind_phys

    if (is_non_zero .eqv. expect_non_zero) then
       print *, "PASS: ", name, " (Dep: ", total_dep, ")"
       total_passed = total_passed + 1
    else
       print *, "FAIL: ", name, " (Dep: ", total_dep, ", Expected non-zero: ", expect_non_zero, ")"
       total_failed = total_failed + 1
    endif
  end subroutine run_case

end program test_wetdep_ls
