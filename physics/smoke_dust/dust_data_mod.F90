!>\file  dust_data_mod.F90
!! This file contains the data for the dust flux schemes.

module dust_data_mod

  use machine ,        only : kind_phys
  use rrfs_smoke_config, only : p_dust_1, p_dust_2, p_dust_3, p_dust_4, p_dust_5, &
                              p_edust1, p_edust2, p_edust3, p_edust4, p_edust5


  implicit none

  integer, parameter :: ndust = 5
  integer, parameter :: ndcls = 3
  integer, parameter :: ndsrc = 1
  integer, parameter :: maxstypes = 100
  integer, parameter :: nsalt = 9

  real(kind_phys),    parameter :: dyn_visc = 1.5E-5

  ! -- dust parameters
  ! never used: integer,         dimension(ndust), parameter :: ipoint    = (/ 3, 2, 2, 2, 2 /)
  real(kind_phys), dimension(ndust), parameter :: den_dust  = (/   2500.,  2650.,  2650.,  2650.,  2650. /)
  real(kind_phys), dimension(ndust), parameter :: reff_dust = (/ 0.73D-6, 1.4D-6, 2.4D-6, 4.5D-6, 8.0D-6 /)
  real(kind_phys), dimension(ndust), parameter :: frac_s    = (/     0.1,   0.25,   0.25,   0.25,   0.25 /)
  real(kind_phys), dimension(ndust), parameter :: lo_dust   = (/  0.1D-6, 1.0D-6, 1.8D-6, 3.0D-6, 6.0D-6 /)
  real(kind_phys), dimension(ndust), parameter :: up_dust   = (/  1.0D-6, 1.8D-6, 3.0D-6, 6.0D-6,10.0D-6 /)
  ! never used: real(kind_phys), dimension(ndust, 12)        :: ch_dust   = 0.8e-09_kind_phys

  ! -- default dust parameters
  ! -- AFWA & GOCART
  ! -----------+----------+-----------+
  ! Parameter  | FIM-Chem | HRRR-Chem |
  ! -----------+----------+-----------+
  ! alpha      |      1.0 |       0.5 |
  ! gamma      |      1.6 |       1.0 |
  ! -----------+----------+-----------+
  ! Never used:
  ! real(kind_phys), parameter :: afwa_alpha    = 0.2
  ! real(kind_phys), parameter :: afwa_gamma    = 1.3
  ! real(kind_phys), parameter :: gocart_alpha  = 0.3
  ! real(kind_phys), parameter :: gocart_gamma  = 1.3
  ! -- FENGSHA
  ! Never used:
  ! real(kind_phys), parameter :: fengsha_alpha = 0.3
  ! real(kind_phys), parameter :: fengsha_gamma = 1.3

  ! -- FENGSHA threshold velocities based on Dale A. Gillette's data
  integer, parameter :: fengsha_maxstypes = 13

  real(kind_phys), dimension(fengsha_maxstypes), parameter :: dust_uthres = &
    (/ 0.065,   & ! Sand            - 1
       0.18,    & ! Loamy Sand      - 2
       0.27,    & ! Sandy Loam      - 3
       0.30,    & ! Silt Loam       - 4
       0.35,    & ! Silt            - 5
       0.38,    & ! Loam            - 6
       0.35,    & ! Sandy Clay Loam - 7
       0.41,    & ! Silty Clay Loam - 8
       0.41,    & ! Clay Loam       - 9
       0.45,    & ! Sandy Clay      - 10
       0.50,    & ! Silty Clay      - 11
       0.45,    & ! Clay            - 12
       9999.0 /)   ! Other           - 13

  ! -- unified densities
  real(kind_phys), parameter :: rho_soil = 2650.0_kind_phys
  real(kind_phys), parameter :: rho_water = 1000.0_kind_phys

  ! -- FENGSHA parameters
  type :: fengsha_params_type
     real(kind_phys) :: mmd_dust = 3.4e-6_kind_phys     !< median mass diameter (m)
     real(kind_phys) :: gsd_dust = 3.0_kind_phys        !< geom. std deviation
     real(kind_phys) :: lambda = 12.0e-6_kind_phys      !< crack propagation length (m)
     real(kind_phys) :: cv = 12.62e-6_kind_phys         !< normalization constant
     real(kind_phys) :: z0s = 1.0e-4_kind_phys          !< Surface roughness for ideal bare surface (m)
     real(kind_phys) :: clay_thresh = 0.2_kind_phys     !< clay fraction threshold
     real(kind_phys) :: cmb = 1.0_kind_phys             !< constant of proportionality
     real(kind_phys) :: kvhmax = 2.0e-4_kind_phys       !< max vertical to horizontal flux ratio
     real(kind_phys) :: frozen_soil_thresh = 268.0_kind_phys !< frozen soil threshold (K)
     real(kind_phys) :: znt_limit = 0.2_kind_phys       !< roughness length limit (m)
     real(kind_phys) :: snow_limit = 0.0_kind_phys      !< snow depth limit (m)
  end type fengsha_params_type

  ! -- FENGSHA uses precalculated drag partition
  integer, parameter :: dust_calcdrag = 1
  ! -- FENGSHA dust moisture parameterization 1:fecan  -  2:shao 
  integer :: dust_moist_opt = 1
  
  real(kind_phys) :: dust_alpha = 1.0
  real(kind_phys) :: dust_gamma = 1.0
  real(kind_phys) :: dust_moist_correction = 1.0
  real(kind_phys) :: dust_drylimit_factor = 1.0

  ! -- sea salt parameters
  integer,            dimension(nsalt), parameter :: spoint    = (/ 1, 2, 2, 2, 2, 2, 3, 3, 3 /)  ! 1 Clay, 2 Silt, 3 Sand
  real(kind_phys), dimension(nsalt), parameter :: reff_salt = &
    (/ 0.71D-6, 1.37D-6, 2.63D-6, 5.00D-6, 9.50D-6, 18.1D-6, 34.5D-6, 65.5D-6, 125.D-6 /)
  real(kind_phys), dimension(nsalt), parameter :: den_salt  = &
    (/   2500.,   2650.,   2650.,   2650.,   2650.,   2650.,   2650.,   2650.,   2650. /)
  real(kind_phys), dimension(nsalt), parameter :: frac_salt = &
    (/      1.,     0.2,     0.2,     0.2,     0.2,     0.2,   0.333,   0.333,   0.333 /)


  ! -- soil vegatation parameters
  integer, parameter :: max_soiltyp = 30
  real(kind_phys), dimension(max_soiltyp), parameter :: &
    maxsmc = (/ 0.421, 0.464, 0.468, 0.434, 0.406, 0.465, &
                0.404, 0.439, 0.421, 0.000, 0.000, 0.000, &
                0.000, 0.000, 0.000, 0.000, 0.000, 0.000, &
                0.000, 0.000, 0.000, 0.000, 0.000, 0.000, &
                0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /)

  ! -- other soil parameters
  ! never used: real(kind_phys), dimension(maxstypes) :: porosity

  public

end module dust_data_mod
