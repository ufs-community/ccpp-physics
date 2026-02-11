!> \file test_dust_fengsha.F90
!! Comprehensive unit test for the FENGSHA dust emission scheme.

program test_dust_fengsha
  use machine, only : kind_phys
  use dust_fengsha_mod, only : gocart_dust_fengsha_driver
  implicit none

  ! Parameters
  integer, parameter :: im = 1, jm = 1, km = 1
  integer, parameter :: num_chem = 20
  integer, parameter :: num_emis_dust = 5
  integer, parameter :: num_soil_layers = 2
  integer, parameter :: p_dust_1 = 10
  integer, parameter :: p_dust_2 = 11
  integer, parameter :: p_dust_3 = 12
  integer, parameter :: p_dust_4 = 13
  integer, parameter :: p_dust_5 = 14

  ! Global data for tests
  real(kind_phys) :: dt, g
  real(kind_phys), dimension(im, km, jm, num_chem) :: chem
  real(kind_phys), dimension(im, km, jm) :: rho_phy
  real(kind_phys), dimension(im, num_soil_layers, jm) :: smois, stemp
  real(kind_phys), dimension(im, km+1, jm) :: p8w
  real(kind_phys), dimension(im, jm) :: ssm, snowh, xland, area, ust, znt, clay, sand, rdrag, uthr
  integer, dimension(im, jm) :: isltyp
  real(kind_phys), dimension(im, 1, jm, num_emis_dust) :: emis_dust

  character(len=128) :: errmsg
  integer :: errflg
  integer :: total_passed, total_failed

  total_passed = 0
  total_failed = 0

  dt = 600.0_kind_phys
  g = 9.80665_kind_phys

  print *, "Starting FENGSHA comprehensive unit test suite..."

  ! Case 1: Standard land emission (should be non-zero)
  call reset_inputs()
  call run_case("Standard land emission", .true.)

  ! Case 2: Water (xland=2.0) (should be zero)
  call reset_inputs()
  xland = 2.0_kind_phys
  call run_case("Water (xland=2.0)", .false.)

  ! Case 3: Frozen soil (stemp < 268) (should be zero)
  call reset_inputs()
  stemp = 260.0_kind_phys
  call run_case("Frozen soil", .false.)

  ! Case 4: Snow cover (snowh > 0) (should be zero)
  call reset_inputs()
  snowh = 0.1_kind_phys
  call run_case("Snow cover", .false.)

  ! Case 5: Roughness too high (znt > 0.2) (should be zero)
  call reset_inputs()
  znt = 0.3_kind_phys
  call run_case("High roughness length", .false.)

  ! Case 6: Invalid soil type (isltyp=15) (should be zero)
  call reset_inputs()
  isltyp = 15
  call run_case("Invalid soil type (15)", .false.)

  ! Case 7: Wind below threshold (should be zero)
  call reset_inputs()
  ust = 0.01_kind_phys
  call run_case("Wind below threshold", .false.)

  ! Case 8: Very dry soil (should have higher emission than moist)
  call reset_inputs()
  smois = 0.01_kind_phys
  call run_case("Very dry soil", .true.)

  print *, "---------------------------------------"
  print *, "Test Summary: ", total_passed, " PASSED, ", total_failed, " FAILED"
  print *, "---------------------------------------"

  if (total_failed > 0) stop 1

contains

  subroutine reset_inputs()
    chem = 1.0e-12_kind_phys
    rho_phy = 1.2_kind_phys
    smois = 0.1_kind_phys
    stemp = 300.0_kind_phys
    p8w(:,1,:) = 101325.0_kind_phys
    p8w(:,2,:) = 100000.0_kind_phys
    ssm = 0.5_kind_phys
    isltyp = 1
    snowh = 0.0_kind_phys
    xland = 1.0_kind_phys ! land
    area = 1000000.0_kind_phys
    ust = 1.0_kind_phys
    znt = 0.01_kind_phys
    clay = 0.2_kind_phys
    sand = 0.4_kind_phys
    rdrag = 0.8_kind_phys
    uthr = 0.3_kind_phys
    emis_dust = 0.0_kind_phys
    errmsg = ''
    errflg = 0
  end subroutine reset_inputs

  subroutine run_case(name, expect_non_zero)
    character(len=*), intent(in) :: name
    logical, intent(in) :: expect_non_zero
    logical :: is_non_zero
    real(kind_phys) :: total_emis

    call gocart_dust_fengsha_driver(dt,              &
         chem,rho_phy,smois,stemp,p8w,ssm,                 &
         isltyp,snowh,xland,area,g,emis_dust,              &
         ust,znt,clay,sand,rdrag,uthr,                     &
         num_emis_dust,num_chem,num_soil_layers,           &
         1, im, 1, jm, 1, km,                              &
         1, im, 1, jm, 1, km,                              &
         1, im, 1, jm, 1, km,                              &
         errmsg, errflg)

    if (errflg /= 0) then
       print *, "FAIL: ", name, " - Driver returned error: ", trim(errmsg)
       total_failed = total_failed + 1
       return
    endif

    total_emis = sum(emis_dust(1,1,1,:))
    is_non_zero = total_emis > 1.0e-20_kind_phys

    if (is_non_zero .eqv. expect_non_zero) then
       print *, "PASS: ", name, " (Emis: ", total_emis, ")"
       total_passed = total_passed + 1
    else
       print *, "FAIL: ", name, " (Emis: ", total_emis, ", Expected non-zero: ", expect_non_zero, ")"
       total_failed = total_failed + 1
    endif
  end subroutine run_case

end program test_dust_fengsha
