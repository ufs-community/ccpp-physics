!> \file test_seas_gocart.F90
!! Comprehensive unit test for the GOCART sea salt emission scheme.

program test_seas_gocart
  use machine, only : kind_phys
  use seas_mod, only : gocart_seasalt_driver
  implicit none

  ! Parameters
  integer, parameter :: im = 1, jm = 1, km = 1
  integer, parameter :: num_chem = 30
  integer, parameter :: num_emis_seas = 5

  ! Global data for tests
  real(kind_phys) :: dt, g, pi
  real(kind_phys), dimension(im, km, jm, num_chem) :: chem
  real(kind_phys), dimension(im, km, jm) :: rho_phy, alt, t_phy, dz8w, u_phy, v_phy
  real(kind_phys), dimension(im, km+1, jm) :: p8w
  real(kind_phys), dimension(im, jm) :: u10, v10, ustar, tsk, xland, xlat, xlong, area, seashelp
  real(kind_phys), dimension(im, 1, jm, num_emis_seas) :: emis_seas

  integer :: total_passed, total_failed

  total_passed = 0
  total_failed = 0

  dt = 600.0_kind_phys
  g = 9.80665_kind_phys
  pi = acos(-1.0_kind_phys)

  print *, "Starting GOCART sea salt comprehensive unit test suite..."

  ! Case 1: Standard ocean emission (xland=0.0) (should be non-zero)
  call reset_inputs()
  xland = 0.0_kind_phys
  call run_case_opt("Standard ocean emission (GOCART)", .true., 1)

  ! Case 2: Land (xland=1.0) (should be zero)
  call reset_inputs()
  xland = 1.0_kind_phys
  call run_case_opt("Land (xland=1.0)", .false., 1)

  ! Case 3: High wind speed
  call reset_inputs()
  xland = 0.0_kind_phys
  u10 = 20.0_kind_phys
  call run_case_opt("High wind speed", .true., 1)

  ! Case 4: seas_opt = 2 (NGAC scheme)
  call reset_inputs()
  xland = 0.0_kind_phys
  call run_case_opt("NGAC sea salt scheme", .true., 2)

  print *, "---------------------------------------"
  print *, "Test Summary: ", total_passed, " PASSED, ", total_failed, " FAILED"
  print *, "---------------------------------------"

  if (total_failed > 0) stop 1

contains

  subroutine reset_inputs()
    chem = 1.0e-12_kind_phys
    rho_phy = 1.2_kind_phys
    alt = 1.0_kind_phys / rho_phy
    t_phy = 300.0_kind_phys
    u_phy = 10.0_kind_phys
    v_phy = 0.0_kind_phys
    dz8w = 50.0_kind_phys
    p8w(:,1,:) = 101325.0_kind_phys
    p8w(:,2,:) = 100000.0_kind_phys
    u10 = 10.0_kind_phys
    v10 = 0.0_kind_phys
    ustar = 0.5_kind_phys
    tsk = 290.0_kind_phys
    xland = 0.0_kind_phys ! ocean
    xlat = 45.0_kind_phys
    xlong = -30.0_kind_phys
    area = 1000000.0_kind_phys
    emis_seas = 0.0_kind_phys
  end subroutine reset_inputs

  subroutine run_case_opt(name, expect_non_zero, opt)
    character(len=*), intent(in) :: name
    logical, intent(in) :: expect_non_zero
    integer, intent(in) :: opt
    logical :: is_non_zero
    real(kind_phys) :: total_emis

    call gocart_seasalt_driver(dt,alt,t_phy,u_phy,             &
         v_phy,chem,rho_phy,dz8w,u10,v10,ustar,p8w,tsk,            &
         xland,xlat,xlong,area,g,emis_seas,pi,                     &
         seashelp,num_emis_seas,num_chem,opt,                      &
         1, im, 1, jm, 1, km,                                      &
         1, im, 1, jm, 1, km,                                      &
         1, im, 1, jm, 1, km                                       )

    total_emis = sum(emis_seas(1,1,1,:))
    is_non_zero = total_emis > 1.0e-20_kind_phys

    if (is_non_zero .eqv. expect_non_zero) then
       print *, "PASS: ", name, " (Emis: ", total_emis, ")"
       total_passed = total_passed + 1
    else
       print *, "FAIL: ", name, " (Emis: ", total_emis, ", Expected non-zero: ", expect_non_zero, ")"
       total_failed = total_failed + 1
    endif
  end subroutine run_case_opt

end program test_seas_gocart
