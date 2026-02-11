!>\file  dust_fengsha_mod.F90
!! This file contains the FENGSHA dust scheme.

module dust_fengsha_mod
!
!  This module developed by Barry Baker (NOAA ARL)
!  For serious questions contact barry.baker@noaa.gov
!
!  07/16/2019 - Adapted for NUOPC/GOCART, R. Montuoro
!  02/01/2020 - Adapted for FV3/CCPP, Haiqin Li
!  Refactored for Fortran 2008 compliance and best practices.

  use machine ,        only : kind_phys
  use dust_data_mod,   only : ndust, reff_dust, lo_dust, up_dust, &
                              dust_calcdrag, dust_moist_opt, dust_alpha, dust_gamma, &
                              dust_moist_correction, dust_drylimit_factor, &
                              p_dust_1, p_dust_2, p_dust_3, p_dust_4, p_dust_5, &
                              p_edust1, p_edust2, p_edust3, p_edust4, p_edust5

  implicit none

  private

  public :: gocart_dust_fengsha_driver

  ! -- unified densities for FENGSHA
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

  !> Parameter set for FENGSHA scheme
  type(fengsha_params_type), parameter :: fparams = fengsha_params_type()

contains

  !> \brief Driver for the GOCART FENGSHA dust emission scheme
  !! \section arg_table_gocart_dust_fengsha_driver Arguments
  !! \htmlinclude gocart_dust_fengsha_driver.html
  !! \param[in] dt physics time step (s)
  !! \param[inout] chem constituent mixing ratio (kg kg-1)
  !! \param[in] rho_phy air density (kg m-3)
  !! \param[in] smois volumetric soil moisture (m3 m-3)
  !! \param[in] stemp soil temperature (K)
  !! \param[in] p8w air pressure at interfaces (Pa)
  !! \param[in] ssm sediment supply map (none)
  !! \param[in] isltyp dominant soil type (index)
  !! \param[in] snowh snow depth (m)
  !! \param[in] xland land mask (1 for land, 2 for water)
  !! \param[in] area grid cell area (m2)
  !! \param[in] g gravitational acceleration (m s-2)
  !! \param[inout] emis_dust optional dust emission flux (kg m-2 s-1)
  !! \param[in] ust friction velocity (m s-1)
  !! \param[in] znt surface roughness length (m)
  !! \param[in] clay clay fraction (none)
  !! \param[in] sand sand fraction (none)
  !! \param[in] rdrag drag partition correction (none)
  !! \param[in] uthr dry threshold velocity (m s-1)
  !! \param[in] num_emis_dust number of dust bins
  !! \param[in] num_chem number of chemistry tracers
  !! \param[in] num_soil_layers number of soil layers
  !! \param[in] ids horizontal dimension start index
  !! \param[in] ide horizontal dimension end index
  !! \param[in] jds horizontal dimension 2 start index
  !! \param[in] jde horizontal dimension 2 end index
  !! \param[in] kds vertical dimension start index
  !! \param[in] kde vertical dimension end index
  !! \param[in] ims horizontal dimension ims
  !! \param[in] ime horizontal dimension ime
  !! \param[in] jms horizontal dimension 2 jms
  !! \param[in] jme horizontal dimension 2 jme
  !! \param[in] kms vertical dimension kms
  !! \param[in] kme vertical dimension kme
  !! \param[in] its horizontal dimension its
  !! \param[in] ite horizontal dimension ite
  !! \param[in] jts horizontal dimension 2 jts
  !! \param[in] jte horizontal dimension 2 jte
  !! \param[in] kts vertical dimension kts
  !! \param[in] kte vertical dimension kte
  !! \param[out] errmsg ccpp error message
  !! \param[out] errflg ccpp error code
  subroutine gocart_dust_fengsha_driver(dt,              &
       chem,rho_phy,smois,stemp,p8w,ssm,                 &
       isltyp,snowh,xland,area,g,emis_dust,              &
       ust,znt,clay,sand,rdrag,uthr,                     &
       num_emis_dust,num_chem,num_soil_layers,           &
       ids,ide, jds,jde, kds,kde,                        &
       ims,ime, jms,jme, kms,kme,                        &
       its,ite, jts,jte, kts,kte,                        &
       errmsg, errflg)

    integer,      intent(in) ::                       &
         ids,ide, jds,jde, kds,kde,                      &
         ims,ime, jms,jme, kms,kme,                      &
         its,ite, jts,jte, kts,kte,                      &
         num_emis_dust,num_chem,num_soil_layers

    ! 2d input variables
    real(kind_phys), dimension( ims:ime , jms:jme ), intent(in) :: ssm     ! Sediment supply map
    real(kind_phys), dimension( ims:ime , jms:jme ), intent(in) :: snowh   ! snow height (m)
    real(kind_phys), dimension( ims:ime , jms:jme ), intent(in) :: xland   ! dominant land use type
    real(kind_phys), dimension( ims:ime , jms:jme ), intent(in) :: area    ! area of grid cell
    real(kind_phys), dimension( ims:ime , jms:jme ), intent(in) :: ust     ! friction velocity
    real(kind_phys), dimension( ims:ime , jms:jme ), intent(in) :: znt     ! Surface Roughness length (m)
    real(kind_phys), dimension( ims:ime , jms:jme ), intent(in) :: clay    ! Clay Fraction (-)
    real(kind_phys), dimension( ims:ime , jms:jme ), intent(in) :: sand    ! Sand Fraction (-)
    real(kind_phys), dimension( ims:ime , jms:jme ), intent(in) :: rdrag   ! Drag Partition (-)
    real(kind_phys), dimension( ims:ime , jms:jme ), intent(in) :: uthr    ! Dry Threshold Velocity (m/s)

    integer,         dimension( ims:ime , jms:jme ), intent(in) :: isltyp  ! soil type

    ! 3d input variables
    real(kind_phys), dimension( ims:ime , kms:kme , jms:jme ), intent(in) :: p8w
    real(kind_phys), dimension( ims:ime , kms:kme , jms:jme ), intent(in) :: rho_phy
    real(kind_phys), dimension( ims:ime, kms:kme, jms:jme, num_chem ), intent(inout) :: chem
    real(kind_phys), dimension( ims:ime, 1, jms:jme,num_emis_dust),optional, intent(inout) :: emis_dust
    real(kind_phys), dimension( ims:ime, num_soil_layers, jms:jme ), intent(in) :: smois, stemp

    !0d input variables 
    real(kind_phys), intent(in) :: dt ! time step
    real(kind_phys), intent(in) :: g  ! gravity (m/s**2)

    ! CCPP error handling
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    ! Local variables
    integer :: i, j
    real(kind_phys), parameter :: conver = 1.0e-9_kind_phys
    real(kind_phys), parameter :: converi = 1.0e9_kind_phys

    ! Initialize error handling
    errmsg = ''
    errflg = 0

    ! OpenMP parallelization over the grid
    !$omp parallel do collapse(2) default(shared) private(i, j)
    do j = jts, jte
       do i = its, ite
          block
             integer :: ilwi
             integer :: nmx
             real(kind_phys) :: airden ! air density
             real(kind_phys) :: airmas ! dry air mass
             real(kind_phys) :: dxy
             real(kind_phys) :: R ! local drag partition
             real(kind_phys) :: ustar
             real(kind_phys) :: tc(num_emis_dust)
             real(kind_phys) :: bems(num_emis_dust)
             real(kind_phys) :: massfrac(3)
             real(kind_phys) :: erodtot
             real(kind_phys) :: moist_volumetric

             nmx = ndust

             ! Don't do dust over water!!!
             ilwi = 0
             if (xland(i,j) < 1.5_kind_phys) then
                ilwi = 1

                ! Total concentration at lowest model level.
                tc(1) = chem(i,kts,j,p_dust_1) * conver
                tc(2) = chem(i,kts,j,p_dust_2) * conver
                tc(3) = chem(i,kts,j,p_dust_3) * conver
                tc(4) = chem(i,kts,j,p_dust_4) * conver
                tc(5) = chem(i,kts,j,p_dust_5) * conver

                ! Air mass and density at lowest model level.
                airmas = -(p8w(i,kts+1,j) - p8w(i,kts,j)) * area(i,j) / g
                airden = rho_phy(i,kts,j)
                ustar = ust(i,j)
                dxy = area(i,j)

                ! Mass fractions of clay, silt, and sand.
                massfrac(1) = clay(i,j)
                massfrac(2) = 1.0_kind_phys - (clay(i,j) + sand(i,j))
                massfrac(3) = sand(i,j)

                ! Total erodibility.
                erodtot = ssm(i,j)

                ! Don't allow roughness lengths greater than limit to be lofted.
                if (znt(i,j) > fparams%znt_limit) then
                   ilwi = 0
                endif

                ! limit where there is snow on the ground
                if (snowh(i,j) > fparams%snow_limit) then
                   ilwi = 0
                endif

                ! Don't emit over frozen soil
                if (stemp(i,1,j) < fparams%frozen_soil_thresh) then
                   ilwi = 0
                endif

                ! Do not allow areas with bedrock, lava, or land-ice to loft
                if (isltyp(i,j) == 15 .or. isltyp(i,j) == 16 .or. &
                     isltyp(i,j) == 18 .or. isltyp(i,j) == 0) then
                   ilwi = 0
                endif

                if (ilwi /= 0) then
                   ! get drag partition
                   if (dust_calcdrag /= 1) then
                      call fengsha_drag(znt(i,j), R)
                   else
                      ! use the precalculated version
                      if (rdrag(i,j) > 0.0_kind_phys) then
                         R = rdrag(i,j)
                      else
                         ilwi = 0
                      endif
                   endif
                endif

                if (ilwi /= 0) then
                   ! soil moisture correction factor
                   moist_volumetric = dust_moist_correction * smois(i,2,j)

                   ! Call dust emission routine.
                   call source_dust(nmx, dt, tc, ustar, massfrac, &
                        erodtot, dxy, moist_volumetric, airden, airmas, bems, g, dust_alpha, dust_gamma, &
                        R, uthr(i,j))

                   ! convert back to concentration
                   chem(i,kts,j,p_dust_1) = tc(1) * converi
                   chem(i,kts,j,p_dust_2) = tc(2) * converi
                   chem(i,kts,j,p_dust_3) = tc(3) * converi
                   chem(i,kts,j,p_dust_4) = tc(4) * converi
                   chem(i,kts,j,p_dust_5) = tc(5) * converi

                   ! For output diagnostics
                   if (present(emis_dust)) then
                      emis_dust(i,1,j,p_edust1) = bems(1)
                      emis_dust(i,1,j,p_edust2) = bems(2)
                      emis_dust(i,1,j,p_edust3) = bems(3)
                      emis_dust(i,1,j,p_edust4) = bems(4)
                      emis_dust(i,1,j,p_edust5) = bems(5)
                   endif
                endif
             endif
          end block
       enddo
    enddo

  end subroutine gocart_dust_fengsha_driver


  !> \brief Evaluates the source of each dust particle size bin
  !! \param[in] nmx number of dust bins
  !! \param[in] dt1 time step (s)
  !! \param[inout] tc total concentration of dust (kg kg-1)
  !! \param[in] ustar friction velocity (m s-1)
  !! \param[in] massfrac fraction of mass in each of 3 soil classes
  !! \param[in] erod fraction of erodible grid cell
  !! \param[in] dxy grid cell area (m2)
  !! \param[in] smois volumetric soil moisture (m3 m-3)
  !! \param[in] airden density of air (kg m-3)
  !! \param[in] airmas mass of air for each grid box (kg)
  !! \param[out] bems source of each dust type (ug m-2 s-1)
  !! \param[in] g0 gravitational acceleration (m s-2)
  !! \param[in] alpha scaling factor
  !! \param[in] gamma scaling factor
  !! \param[in] R drag partition
  !! \param[in] uthres dry threshold velocity (m s-1)
  subroutine source_dust(nmx, dt1, tc, ustar, massfrac, &
                  erod, dxy, smois, airden, airmas, bems, g0, alpha, gamma, &
                  R, uthres)

    integer,            intent(in)    :: nmx
    real(kind_phys),    intent(in)    :: dt1
    real(kind_phys),    intent(in)    :: ustar
    real(kind_phys),    intent(in)    :: massfrac(3)
    real(kind_phys),    intent(in)    :: erod
    real(kind_phys),    intent(in)    :: dxy
    real(kind_phys),    intent(in)    :: smois
    real(kind_phys),    intent(in)    :: airden
    real(kind_phys),    intent(in)    :: airmas
    real(kind_phys),    intent(in)    :: g0
    real(kind_phys),    intent(in)    :: alpha
    real(kind_phys),    intent(in)    :: gamma
    real(kind_phys),    intent(in)    :: R
    real(kind_phys),    intent(in)    :: uthres

    ! Output
    real(kind_phys),    intent(inout) :: tc(nmx)
    real(kind_phys),    intent(out)   :: bems(nmx)

    ! Local Variables
    real(kind_phys) :: dvol(nmx)
    real(kind_phys) :: distr_dust(nmx)
    real(kind_phys) :: dlndp(nmx)
    real(kind_phys) :: dsrc
    real(kind_phys) :: dvol_tot
    real(kind_phys) :: emit
    integer         :: n

    ! calculate the total vertical dust flux 
    emit = 0.0_kind_phys

    call DustEmissionFENGSHA(smois, massfrac(1), massfrac(3), massfrac(2), &
                             erod, R, airden, ustar, uthres, alpha, gamma, &
                             g0, emit)

    ! Now that we have the total dust emission, distribute into dust bins using
    ! lognormal distribution (Dr. Jasper Kok)
    dvol_tot = 0.0_kind_phys
    do n = 1, nmx
       dlndp(n) = log(up_dust(n) / lo_dust(n))
       dvol(n) = (2.0_kind_phys * reff_dust(n) / fparams%cv) * &
            (1.0_kind_phys + erf(log(2.0_kind_phys * reff_dust(n) / fparams%mmd_dust) / &
            (sqrt(2.0_kind_phys) * log(fparams%gsd_dust)))) * &
            exp(-(2.0_kind_phys * reff_dust(n) / fparams%lambda)**3.0_kind_phys) * dlndp(n)
       dvol_tot = dvol_tot + dvol(n)
    end do

    do n = 1, nmx
       distr_dust(n) = dvol(n) / dvol_tot
    end do

    ! Now distribute total vertical emission into dust bins and update concentration.
    do n = 1, nmx
       ! Calculate total mass emitted
       dsrc = emit * distr_dust(n) * dxy * dt1  ! (kg)
       if (dsrc < 0.0_kind_phys) dsrc = 0.0_kind_phys
       
       ! Update dust mixing ratio at first model level.
       tc(n) = tc(n) + dsrc / airmas ! (kg/kg)
       bems(n) = 1.0e9_kind_phys * dsrc / (dxy * dt1) ! diagnostic (ug/m2/s)
    end do

    ! RRFS-SD Kludge
    tc(1) = tc(1) + 0.286_kind_phys * tc(2)
    tc(5) = 0.714_kind_phys * tc(2) + tc(3) + tc(4)

  end subroutine source_dust


  !> \brief Calculates the MacKinnon et al. 2004 Drag Partition Correction
  !! \param[in] z0 surface roughness length (m)
  !! \param[out] R drag partition correction
  subroutine fengsha_drag(z0, R)
    real(kind_phys), intent(in)  :: z0
    real(kind_phys), intent(out) :: R

    ! Drag partition correction. See MacKinnon et al. (2004),
    ! doi:10.1016/j.geomorph.2004.03.009
    R = 1.0_kind_phys - log(z0 / fparams%z0s) / log(0.7_kind_phys * (12255.0_kind_phys / fparams%z0s)**0.8_kind_phys)

  end subroutine fengsha_drag

  !> \brief Computes dust emissions using NOAA/ARL FENGSHA model
  !! \param[in] slc volumetric soil moisture fraction
  !! \param[in] clay fractional clay content
  !! \param[in] sand fractional sand content
  !! \param[in] silt fractional silt content
  !! \param[in] ssm erosion map
  !! \param[in] rdrag drag partition
  !! \param[in] airdens air density at lowest level (kg m-3)
  !! \param[in] ustar friction velocity (m s-1)
  !! \param[in] uthrs threshold velocity (m s-1)
  !! \param[in] alpha scaling factor
  !! \param[in] gamma scaling factor
  !! \param[in] grav gravitational acceleration (m s-2)
  !! \param[inout] emissions total surface emissions (kg m-2 s-1)
  subroutine DustEmissionFENGSHA(slc, clay, sand, silt,  &
                                  ssm, rdrag, airdens, ustar, uthrs, alpha, gamma, &
                                  grav, emissions)
    
    real(kind_phys), intent(in) :: slc      ! liquid water content of soil layer, volumetric fraction [1]
    real(kind_phys), intent(in) :: clay     ! fractional clay content [1]
    real(kind_phys), intent(in) :: sand     ! fractional sand content [1]
    real(kind_phys), intent(in) :: silt     ! fractional silt content [1]
    real(kind_phys), intent(in) :: ssm      ! erosion map [1]
    real(kind_phys), intent(in) :: rdrag    ! drag partition [1/m]
    real(kind_phys), intent(in) :: airdens  ! air density at lowest level [kg/m^3]
    real(kind_phys), intent(in) :: ustar    ! friction velocity [m/sec]
    real(kind_phys), intent(in) :: uthrs    ! threshold velocity [m/s]
    real(kind_phys), intent(in) :: alpha    ! scaling factor [1]
    real(kind_phys), intent(in) :: gamma    ! scaling factor [1]
    real(kind_phys), intent(in) :: grav     ! gravity [m/sec^2]
    
    real(kind_phys), intent(inout) :: emissions ! surface emissions [kg/(m^2 sec)]
    
    ! Local Variables
    real(kind_phys) :: alpha_grav
    real(kind_phys) :: h
    real(kind_phys) :: kvh
    real(kind_phys) :: q
    real(kind_phys) :: rustar
    real(kind_phys) :: u_sum, u_thresh
    
    emissions = 0.0_kind_phys
    alpha_grav = alpha / grav

    ! Compute vertical-to-horizontal mass flux ratio
    kvh = DustFluxV2HRatioMB95(clay)

    ! Compute total emissions base
    emissions = alpha_grav * (ssm ** gamma) * airdens * kvh

    ! Compute threshold wind friction velocity using drag partition
    rustar = rdrag * ustar

    if (dust_moist_opt == 1) then
       ! Fecan moisture correction
       h = moistureCorrectionFecan(slc, sand, clay)
    else
       ! shao soil moisture correction
       h = moistureCorrectionShao(slc)
    end if

    ! Adjust threshold
    u_thresh = uthrs * h
    u_sum = rustar + u_thresh
   
    ! Compute Horizontal Saltation Flux according to Eq (9) in Webb et al. (2020)
    q = max(0.0_kind_phys, rustar - u_thresh) * u_sum * u_sum
   
    ! Distribute emissions and convert to mass flux (kg m-2 s-1)
    emissions = emissions * q

  end subroutine DustEmissionFENGSHA

  !> \brief Convert soil moisture fraction from volumetric to gravimetric
  !! \param[in] vsoil volumetric soil moisture fraction
  !! \param[in] sandfrac fractional sand content
  !! \return gravimetric soil moisture fraction
  function soilMoistureConvertVol2Grav(vsoil, sandfrac) result(res)
    real(kind_phys), intent(in) :: vsoil       ! volumetric soil moisture fraction [1]
    real(kind_phys), intent(in) :: sandfrac    ! fractional sand content [1]
    real(kind_phys)             :: res
    
    real(kind_phys) :: vsat

    ! Saturated volumetric water content (sand-dependent) [m3 m-3]
    vsat = 0.489_kind_phys - 0.126_kind_phys * sandfrac

    ! Gravimetric soil content
    res = 100.0_kind_phys * (vsoil * rho_water / rho_soil / (1.0_kind_phys - vsat))

  end function soilMoistureConvertVol2Grav

  !> \brief Compute correction factor to account for Fecan soil moisture
  !! \param[in] slc liquid water content of top soil layer
  !! \param[in] sand fractional sand content
  !! \param[in] clay fractional clay content
  !! \return correction factor
  function moistureCorrectionFecan(slc, sand, clay) result(res)
    real(kind_phys), intent(in) :: slc     ! liquid water content of top soil layer
    real(kind_phys), intent(in) :: sand    ! fractional sand content
    real(kind_phys), intent(in) :: clay    ! fractional clay content
    real(kind_phys)             :: res

    real(kind_phys) :: grvsoilm
    real(kind_phys) :: drylimit

    ! Convert soil moisture from volumetric to gravimetric
    grvsoilm = soilMoistureConvertVol2Grav(slc, sand)

    ! Compute fecan dry limit
    drylimit = dust_drylimit_factor * clay * (14.0_kind_phys * clay + 17.0_kind_phys)

    ! Compute soil moisture correction
    res = sqrt(1.0_kind_phys + 1.21_kind_phys * max(0.0_kind_phys, grvsoilm - drylimit)**0.68_kind_phys)

  end function moistureCorrectionFecan

  !> \brief Compute correction factor to account for Shao soil moisture
  !! \param[in] slc liquid water content of top soil layer
  !! \return correction factor
  function moistureCorrectionShao(slc) result(res)
    real(kind_phys), intent(in) :: slc
    real(kind_phys)             :: res

    if (slc < 0.03_kind_phys) then
       res = exp(22.7_kind_phys * slc)
    else
       res = exp(95.3_kind_phys * slc - 2.029_kind_phys)
    end if

  end function moistureCorrectionShao

  !> \brief Computes the vertical-to-horizontal dust flux ratio
  !! \param[in] clay fractional clay content
  !! \return flux ratio
  function DustFluxV2HRatioMB95(clay) result(res)
    real(kind_phys), intent(in) :: clay      ! fractional clay content
    real(kind_phys)             :: res

    if (clay > fparams%clay_thresh) then
       res = fparams%kvhmax
    else
       res = 10.0_kind_phys**(13.4_kind_phys * clay - 6.0_kind_phys)
    end if

  end function DustFluxV2HRatioMB95
  
end module dust_fengsha_mod
