!>\file mp_tempo.F90
!! This file contains aerosol-aware TEMPO MP scheme.


!>\defgroup aatempo Aerosol-Aware TEMPO MP Module
!! This module contains the aerosol-aware TEMPO microphysics scheme.
module mp_tempo

      use mpi_f08
      use machine, only : kind_phys

      use module_mp_tempo_params
      use module_mp_tempo_cfgs, only : ty_tempo_cfgs
      use module_mp_tempo_driver, only : tempo_init, tempo_run, ty_tempo_driver_diags, tempo_aerosol_surface_emissions

      implicit none

      public :: mp_tempo_init, mp_tempo_run, mp_tempo_finalize

      private

   contains

!> This subroutine is a wrapper around the actual tempo_init().
!! \section arg_table_mp_tempo_init Argument Table
!! \htmlinclude mp_tempo_init.html
!!
      subroutine mp_tempo_init(ncol, nlev, &
           imp_physics, imp_physics_tempo, &
           mpicomm, mpirank, mpiroot, &
           tgrs, prsl, phil, con_g, con_rd, con_eps, &
           restart, convert_dry_rho, is_aerosol_aware, &
           is_hail_aware, do_sat_adj, semi_sedi, &
           spechum, nwfa, nifa, nwfa2d, nifa2d, &
           tempo_cfgs, is_initialized, errmsg, errflg)
         
         ! Interface variables
         integer,                   intent(in   ) :: ncol
         integer,                   intent(in   ) :: nlev
         logical,                   intent(in   ) :: restart
         logical,                   intent(inout) :: is_initialized
         integer,                   intent(in   ) :: imp_physics
         integer,                   intent(in   ) :: imp_physics_tempo
         logical,                   intent(in   ) :: do_sat_adj
         logical,                   intent(in   ) :: semi_sedi
         logical,                   intent(in   ) :: convert_dry_rho
         logical,                   intent(in   ) :: is_aerosol_aware
         logical,                   intent(in   ) :: is_hail_aware
         real(kind_phys),           intent(in   ) :: con_g, con_rd, con_eps
         ! Hydrometeors
         real(kind_phys),           intent(inout) :: spechum(:,:)
         ! Aerosols
         real(kind_phys),           intent(inout), optional :: nwfa(:,:)
         real(kind_phys),           intent(inout), optional :: nifa(:,:)
         real(kind_phys),           intent(inout), optional :: nwfa2d(:)
         real(kind_phys),           intent(inout), optional :: nifa2d(:)

         ! State variables
         real(kind_phys),           intent(in   ) :: tgrs(:,:)
         real(kind_phys),           intent(in   ) :: prsl(:,:)
         real(kind_phys),           intent(in   ) :: phil(:,:)
         ! MPI information
         type(MPI_Comm),            intent(in   ) :: mpicomm
         integer,                   intent(in   ) :: mpirank
         integer,                   intent(in   ) :: mpiroot
         ! CCPP error handling
         character(len=*),          intent(  out) :: errmsg
         integer,                   intent(  out) :: errflg
         type(ty_tempo_cfgs),       intent(inout) :: tempo_cfgs
         
         real(kind_phys) :: qv(1:ncol,1:nlev)       ! kg kg-1 (water vapor mixing ratio)
         real(kind_phys) :: hgt(1:ncol,1:nlev)      ! m
         real(kind_phys) :: rho(1:ncol,1:nlev)      ! kg m-3
         real(kind_phys) :: orho(1:ncol,1:nlev)     ! m3 kg-1
         
         real (kind=kind_phys) :: h_01, z1, niIN3, niCCN3
         integer :: i, k
         
         ! Initialize the CCPP error handling variables
         errmsg = ''
         errflg = 0

          if (do_sat_adj) then
            if ((is_aerosol_aware) .or. (is_hail_aware)) then
               write(errmsg, fmt='((a))') 'do_sat_adj should be run with is_aerosol_aware=F and is_hail_aware=F'
               errflg = 1
               return
            endif
         end if

         if (is_initialized) return
         
         ! Consistency checks
         if (imp_physics/=imp_physics_tempo) then
            write(errmsg,'(*(a))') "Logic error: namelist choice of microphysics is different from Tempo MP"
            errflg = 1
            return
         end if

         ! Call tempo init (also sets initial default values of physical constants)
         if (mpirank==mpiroot) write(*,*) 'Calling tempo_init() with ltaerosol= ', is_aerosol_aware, &
              ' lthailaware= ', is_hail_aware, ' sedi_semi= ', semi_sedi, ' do_sat_adj= ' do_sat_adj

         ! Main call to tempo_init()
         call tempo_init(aerosolaware_flag=is_aerosol_aware, hailaware_flag=is_hail_aware, &
              semi_sedi_flag=semi_sedi, cloud_condensation_flag=(.not. do_sat_adj), &
              tempo_cfgs=tempo_cfgs)

         if (errflg /= 0) return

         ! For restart runs, the init is done here
         if (restart) then
           is_initialized = .true.
           return
         end if

         where(spechum<0) spechum = 1.0e-10
         qv = spechum/(1.0_kind_phys-spechum)         
         if (convert_dry_rho) then
           if (is_aerosol_aware) then
              nwfa = nwfa/(1.0_kind_phys-spechum)
              nifa = nifa/(1.0_kind_phys-spechum)
           end if
         end if

         ! Geopotential height in m2 s-2 to height in m
         hgt = phil/con_g
         
         ! Density of moist air in kg m-3 and inverse density of air
         rho = con_eps*prsl/(con_rd*tgrs*(qv+con_eps))
         orho = 1.0/rho

         ! Check for existing aerosol data, both CCN and IN aerosols.  If missing
         ! fill in just a basic vertical profile, somewhat boundary-layer following.
         if (is_aerosol_aware) then

           ! Potential cloud condensation nuclei (CCN)
           if (MAXVAL(nwfa) .lt. eps) then
             if (mpirank==mpiroot) write(*,*) ' There are no initial CCN aerosols. A basic vertical profile will be created.'
             do i = 1, ncol
               if (hgt(i,1).le.1000.0) then
                 h_01 = 0.8
               elseif (hgt(i,1).ge.2500.0) then
                 h_01 = 0.01
               else
                 h_01 = 0.8*cos(hgt(i,1)*0.001 - 1.0)
               endif
               niCCN3 = -1.0*ALOG(naCCN1/naCCN0)/h_01
               nwfa(i,1) = naCCN1+naCCN0*exp(-((hgt(i,2)-hgt(i,1))/1000.)*niCCN3)
               z1 = hgt(i,2)-hgt(i,1)
               nwfa2d(i) = nwfa(i,1) * 0.000196 * (50./z1)
               do k = 2, nlev
                 nwfa(i,k) = naCCN1+naCCN0*exp(-((hgt(i,k)-hgt(i,1))/1000.)*niCCN3)
               enddo
             enddo
           else
             if (mpirank==mpiroot) write(*,*) ' Initial CCN aerosols are present.'
             if (MAXVAL(nwfa2d) .lt. eps) then
               !+---+-----------------------------------------------------------------+
               !..Scale the lowest level aerosol data into an emissions rate.  This is
               !.. very far from ideal, but need higher emissions where larger amount
               !.. of (climo) existing and lesser emissions where there exists fewer to
               !.. begin as a first-order simplistic approach.  Later, proper connection to
               !.. emission inventory would be better.
               !+---+-----------------------------------------------------------------+
               if (mpirank==mpiroot) write(*,*) ' There are no initial CCN aerosol surface emission rates. Rates will be created from surface values.'
               do i = 1, ncol
                  z1 = hgt(i,2)-hgt(i,1)
                  nwfa2d(i) = nwfa(i,1) * 0.000196 * (5./z1)
               enddo
             else
                if (mpirank==mpiroot) write(*,*) ' Initial CCN aerosol surface emission rates are present.'
             endif
           endif

           ! Potential ice nuclei (IN)
           if (MAXVAL(nifa) .lt. eps) then
             if (mpirank==mpiroot) write(*,*) ' There are no initial IN aerosols. A basic vertical profile will be created.'
             do i = 1, ncol
               if (hgt(i,1).le.1000.0) then
                  h_01 = 0.8
               elseif (hgt(i,1).ge.2500.0) then
                  h_01 = 0.01
               else
                  h_01 = 0.8*cos(hgt(i,1)*0.001 - 1.0)
               endif
               niIN3 = -1.0*ALOG(naIN1/naIN0)/h_01
               nifa(i,1) = naIN1+naIN0*exp(-((hgt(i,2)-hgt(i,1))/1000.)*niIN3)
               nifa2d(i) = 0.
               do k = 2, nlev
                  nifa(i,k) = naIN1+naIN0*exp(-((hgt(i,k)-hgt(i,1))/1000.)*niIN3)
               enddo
             enddo
           else
             if (mpirank==mpiroot) write(*,*) ' Initial IN aerosols are present.'
             if (MAXVAL(nifa2d) .lt. eps) then
               if (mpirank==mpiroot) write(*,*) ' There are no initial IN aerosol surface emission rates. Rates will be set to zero.'
               ! calculate IN surface flux here, right now just set to zero
               nifa2d = 0.
             else
               if (mpirank==mpiroot) write(*,*) ' Initial IN aerosol surface emission rates are present.'
             endif
           endif

           ! Ensure non-negative aerosol number concentrations.
           where(nwfa .LE. 0.0) nwfa = 1.1E6
           where(nifa .LE. 0.0) nifa = naIN1*0.01
         end if

         if (convert_dry_rho) then
           if (is_aerosol_aware) then
              nwfa = nwfa/(1.0_kind_phys+qv)
              nifa = nifa/(1.0_kind_phys+qv)
           end if
         end if

         is_initialized = .true.

      end subroutine mp_tempo_init


!> \section arg_table_mp_tempo_run Argument Table
!! \htmlinclude mp_tempo_run.html
!!
!>\ingroup aatempo
!>\section gen_tempo TEMPO MP General Algorithm
      subroutine mp_tempo_run(ncol, nlev, &
        mpicomm, mpirank, mpiroot, blkno, &
        convert_dry_rho, dtp, dt_inner, &
        spechum, qc, qr, qi, qs, qg, ni, nr, &
        nc, nwfa, nifa, nwfa2d, nifa2d, ng, volg, &
        con_g, con_rd, con_eps, &
        tgrs, prsl, phii, omega, &
        is_aerosol_aware, is_hail_aware, &
        prcp, rain, graupel, ice, snow, sr, refl_10cm, &
        do_radar_ref, &
        is_initialized, tempo_cfgs, errmsg, errflg)


         ! Interface variables
         logical,                   intent(in   ) :: is_initialized
         logical,                   intent(in   ) :: convert_dry_rho
         logical,                   intent(in   ) :: do_radar_ref
         ! Dimensions and constants
         integer,                   intent(in   ) :: ncol
         integer,                   intent(in   ) :: nlev
         real(kind_phys),           intent(in   ) :: con_g
         real(kind_phys),           intent(in   ) :: con_rd
         real(kind_phys),           intent(in   ) :: con_eps
         ! Hydrometeors
         real(kind_phys),           intent(inout) :: spechum(:,:)
         real(kind_phys),           intent(inout) :: qc(:,:)
         real(kind_phys),           intent(inout) :: qr(:,:)
         real(kind_phys),           intent(inout) :: qi(:,:)
         real(kind_phys),           intent(inout) :: qs(:,:)
         real(kind_phys),           intent(inout) :: qg(:,:)
         real(kind_phys),           intent(inout) :: ni(:,:)
         real(kind_phys),           intent(inout) :: nr(:,:)
         real(kind_phys), volatile, optional, intent(inout) :: nc(:,:)
         real(kind_phys), volatile, optional, intent(inout) :: nwfa(:,:)
         real(kind_phys), volatile, optional, intent(inout) :: nifa(:,:)
         real(kind_phys), optional, intent(in   ) :: nwfa2d(:)
         real(kind_phys), optional, intent(in   ) :: nifa2d(:)
         real(kind_phys), volatile, optional, intent(inout) :: ng(:,:)
         real(kind_phys), volatile, optional, intent(inout) :: volg(:,:)
         logical,                   intent(in)    :: is_aerosol_aware
         logical,                   intent(in)    :: is_hail_aware
         ! Precip/rain/snow/graupel fall amounts and fraction of frozen precip
         real(kind_phys),           intent(inout) :: prcp(:)
         real(kind_phys),           intent(inout) :: rain(:)
         real(kind_phys),           intent(inout) :: graupel(:)
         real(kind_phys),           intent(inout) :: ice(:)
         real(kind_phys),           intent(inout) :: snow(:)
         real(kind_phys),           intent(  out) :: sr(:)
         ! Radar reflectivity
         real(kind_phys),           intent(inout) :: refl_10cm(:,:)         
         ! State variables and timestep information
         real(kind_phys),           intent(inout) :: tgrs(:,:)
         real(kind_phys),           intent(in   ) :: prsl(:,:)
         real(kind_phys),           intent(in   ) :: phii(:,:)
         real(kind_phys),           intent(in   ) :: omega(:,:)
         real(kind_phys),           intent(in   ) :: dtp
         real,                      intent(in   ) :: dt_inner
         logical,                   intent(in   ) :: first_time_step
         ! MPI and block information
         integer,                   intent(in)    :: blkno
         type(MPI_Comm),            intent(in)    :: mpicomm
         integer,                   intent(in)    :: mpirank
         integer,                   intent(in)    :: mpiroot
         ! CCPP error handling
         character(len=*),          intent(  out) :: errmsg
         integer,                   intent(  out) :: errflg

         type(ty_tempo_cfgs), intent(in) :: tempo_cfgs
         type(ty_tempo_driver_diags) :: tempo_driver_diags
         ! Local variables

         ! Reduced time step if dt_inner
         real(kind_phys) :: dt
         ! Air density
         real(kind_phys) :: rho(1:ncol,1:nlev)              !< kg m-3
         ! Water vapor mixing ratio (instead of specific humidity)
         real(kind_phys) :: qv(1:ncol,1:nlev)               !< kg kg-1
         ! Vertical velocity and level width
         real(kind_phys) :: w(1:ncol,1:nlev)                !< m s-1
         real(kind_phys) :: dz(1:ncol,1:nlev)               !< m
         real(kind_phys) :: xnwfa(1:ncol,1:nlev,1)
         real(kind_phys) :: xnwfa2d(1:ncol,1)

         ! Dimensions
         integer :: ndt, i, k
         integer         :: ids,ide, jds,jde, kds,kde, &
                            ims,ime, jms,jme, kms,kme, &
                            its,ite, jts,jte, kts,kte
         integer :: itimestep = 1
         
         ! Initialize the CCPP error handling variables
         errmsg = ''
         errflg = 0

         if (first_time_step .and. blkno==1) then
            ! Check initialization state
            if (.not.is_initialized) then
               write(errmsg, fmt='((a))') 'mp_tempo_run called before mp_tempo_init'
               errflg = 1
               return
            end if
         endif
         
         ndt = max(nint(dtp/dt_inner), 1)
         dt = dtp/ndt
         if (dt <= dt_inner) dt = dt_inner

         !> - Convert specific humidity to water vapor mixing ratio.
         !> - Also, hydrometeor variables are mass or number mixing ratio
         !> - either kg of species per kg of dry air, or per kg of (dry + vapor).
         qv = spechum/(1.0_kind_phys-spechum)

         if (convert_dry_rho) then
           qc = qc/(1.0_kind_phys-spechum)
           qr = qr/(1.0_kind_phys-spechum)
           qi = qi/(1.0_kind_phys-spechum)
           qs = qs/(1.0_kind_phys-spechum)
           qg = qg/(1.0_kind_phys-spechum)
           ni = ni/(1.0_kind_phys-spechum)
           nr = nr/(1.0_kind_phys-spechum)
           if (is_hail_aware) then
              ng = ng/(1.0_kind_phys-spechum)
              volg = volg/(1.0_kind_phys-spechum)
           endif
           if (is_aerosol_aware) then
              nc = nc/(1.0_kind_phys-spechum)
              nwfa = nwfa/(1.0_kind_phys-spechum)
              nifa = nifa/(1.0_kind_phys-spechum)
           end if
         end if

         !> - Density of air in kg m-3
         rho = con_eps*prsl/(con_rd*tgrs*(qv+con_eps))

         !> - Convert omega in Pa s-1 to vertical velocity w in m s-1
         w = -omega/(rho*con_g)

         !> - Layer width in m from geopotential in m2 s-2
         dz = (phii(:,2:nlev+1) - phii(:,1:nlev)) / con_g

         ! Set internal dimensions
         ids = 1
         ims = 1
         its = 1
         ide = ncol
         ime = ncol
         ite = ncol
         jds = 1
         jms = 1
         jts = 1
         jde = 1
         jme = 1
         jte = 1
         kds = 1
         kms = 1
         kts = 1
         kde = nlev
         kme = nlev
         kte = nlev

         if (present(nwfa) .and. present(nwfa2d)) then
            xnwfa(:,:,1) = nwfa(:,:)
            xnwfa2d(:,1) = nwfa2d(:)
            call tempo_aerosol_surface_emissions(dt=dt, nwfa=xnwfa, nwfa2d=xnwfa2d, ims=ims, ime=ime, &
                 jms=jms, jme=jme, kms=kms, kme=kme, kts=kts)
            nwfa(:,:) = xnwfa(:,:,1)
         endif

         call tempo_run(tempo_cfgs=tempo_cfgs, &
            dt=dt, itimestep=itimestep , &
            qv=qv, qc=qc, qr=qr, qi=qi, qs=qs, qg=qg, ni=ni, nr=nr, &
            nc=nc, nwfa=nwfa, nifa=nifa, &
            ng=ng, qb=volg, &
            w=w, t=tgrs, p=prsl, dz=dz, &
            ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde , &
            ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme , &
            its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte , &
            tempo_diags=tempo_driver_diags)

         itimestep = itimestep + 1
         
         if (errflg/=0) return

         !> - Convert water vapor mixing ratio back to specific humidity
         spechum = qv/(1.0_kind_phys+qv)

         if (convert_dry_rho) then
           qc = qc/(1.0_kind_phys+qv)
           qr = qr/(1.0_kind_phys+qv)
           qi = qi/(1.0_kind_phys+qv)
           qs = qs/(1.0_kind_phys+qv)
           qg = qg/(1.0_kind_phys+qv)
           ni = ni/(1.0_kind_phys+qv)
           nr = nr/(1.0_kind_phys+qv)
           if (is_hail_aware) then
              ng = ng/(1.0_kind_phys+qv)
              volg = volg/(1.0_kind_phys+qv)
           endif
           if (is_aerosol_aware) then
              nc = nc/(1.0_kind_phys+qv)
              nwfa = nwfa/(1.0_kind_phys+qv)
              nifa = nifa/(1.0_kind_phys+qv)
           end if
         end if

         if (do_radar_ref) then
            refl_10cm = tempo_driver_diags%refl10cm(:,:,1)
         endif
         
         ice = ice + max(0.0, tempo_driver_diags%ice_liquid_equiv_precip(:,1)/1000.0_kind_phys)
         snow = snow + (max(0.0, tempo_driver_diags%ice_liquid_equiv_precip(:,1)) + &
              max(0.0, tempo_driver_diags%snow_liquid_equiv_precip(:,1)))/1000.0_kind_phys
         graupel = graupel + max(0.0, tempo_driver_diags%graupel_liquid_equiv_precip(:,1)/1000.0_kind_phys)
         rain = rain + max(0.0, tempo_driver_diags%rain_precip(:,1)/1000.0_kind_phys)                           
         prcp = prcp + (max(0.0, tempo_driver_diags%ice_liquid_equiv_precip(:,1)) + &
              max(0.0, tempo_driver_diags%snow_liquid_equiv_precip(:,1)) + &
              max(0.0, tempo_driver_diags%graupel_liquid_equiv_precip(:,1)) + &
              max(0.0, tempo_driver_diags%rain_precip(:,1)))/1000._kind_phys
         sr = tempo_driver_diags%frozen_fraction(:,1)

      end subroutine mp_tempo_run

!> \section arg_table_mp_tempo_finalize Argument Table
!! \htmlinclude mp_tempo_finalize.html
!!
      subroutine mp_tempo_finalize(is_initialized, errmsg, errflg)
        
         logical,                   intent(inout) :: is_initialized
         character(len=*),          intent(  out) :: errmsg
         integer,                   intent(  out) :: errflg

         ! Initialize the CCPP error handling variables
         errmsg = ''
         errflg = 0

         if (.not.is_initialized) return

         is_initialized = .false.

      end subroutine mp_tempo_finalize

end module mp_tempo
