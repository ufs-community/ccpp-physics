!> \file test_dust_fengsha.F90
!! Unit test for the FENGSHA dust emission scheme.

program test_dust_fengsha
  use machine, only : kind_phys
  use dust_fengsha_mod, only : gocart_dust_fengsha_driver
  implicit none

  ! Parameters
  integer, parameter :: im = 2, jm = 1, km = 1
  integer, parameter :: num_chem = 20
  integer, parameter :: num_emis_dust = 5
  integer, parameter :: num_soil_layers = 2
  integer, parameter :: p_dust_1 = 10
  integer, parameter :: p_dust_2 = 11
  integer, parameter :: p_dust_3 = 12
  integer, parameter :: p_dust_4 = 13
  integer, parameter :: p_dust_5 = 14

  ! Arguments
  real(kind_phys) :: dt, g
  real(kind_phys), dimension(im, km, jm, num_chem) :: chem
  real(kind_phys), dimension(im, jm) :: rho_phy_2d ! for lowest level
  real(kind_phys), dimension(im, km, jm) :: rho_phy
  real(kind_phys), dimension(im, num_soil_layers, jm) :: smois, stemp
  real(kind_phys), dimension(im, km+1, jm) :: p8w
  real(kind_phys), dimension(im, jm) :: ssm, snowh, xland, area, ust, znt, clay, sand, rdrag, uthr
  integer, dimension(im, jm) :: isltyp
  real(kind_phys), dimension(im, 1, jm, num_emis_dust) :: emis_dust

  character(len=128) :: errmsg
  integer :: errflg

  ! Initialize
  dt = 600.0_kind_phys
  g = 9.80665_kind_phys

  chem = 1.0e-12_kind_phys
  rho_phy = 1.2_kind_phys
  smois = 0.1_kind_phys
  stemp = 300.0_kind_phys
  p8w(:,1,:) = 100000.0_kind_phys
  p8w(:,2,:) = 99000.0_kind_phys
  ssm = 0.5_kind_phys
  isltyp = 1
  snowh = 0.0_kind_phys
  xland = 1.0_kind_phys ! land
  area = 1000000.0_kind_phys
  ust = 0.5_kind_phys
  znt = 0.01_kind_phys
  clay = 0.2_kind_phys
  sand = 0.4_kind_phys
  rdrag = 0.8_kind_phys
  uthr = 0.3_kind_phys
  emis_dust = 0.0_kind_phys

  print *, "Starting FENGSHA unit test..."

  call gocart_dust_fengsha_driver(dt,              &
       chem,rho_phy,smois,stemp,p8w,ssm,                 &
       isltyp,snowh,xland,area,g,emis_dust,              &
       ust,znt,clay,sand,rdrag,uthr,                     &
       num_emis_dust,num_chem,num_soil_layers,           &
       1, im, 1, jm, 1, km,                              & ! ids, ide, jds, jde, kds, kde
       1, im, 1, jm, 1, km,                              & ! ims, ime, jms, jme, kms, kme
       1, im, 1, jm, 1, km,                              & ! its, ite, jts, jte, kts, kte
       errmsg, errflg)

  if (errflg /= 0) then
     print *, "Test FAILED with error: ", trim(errmsg)
     stop 1
  else
     print *, "Test PASSED."
     print *, "Bin 1 dust concentration: ", chem(1,1,1,p_dust_1)
     print *, "Bin 1 dust emission: ", emis_dust(1,1,1,1)
  endif

end program test_dust_fengsha
