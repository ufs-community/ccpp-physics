!> \file GFS_suite_interstitial_1.f90
!!  Contains code to calculate scale-aware variables used in cs_conv, gwdc, and precpd and to reset tendencies used in the
!!  process-split section of GFS-based physics suites.

    module GFS_suite_interstitial_1

    contains

!> \section arg_table_GFS_suite_interstitial_1_run Argument Table
!! \htmlinclude GFS_suite_interstitial_1_run.html
!!
    subroutine GFS_suite_interstitial_1_run (im, levs, ntrac, dtf, dtp, slmsk, area, dxmin, dxinv, pgr, &
         islmsk, work1, work2, psurf, dudt, dvdt, dtdt, dqdt, ten_u, ten_v, ten_t, ten_q, kpbl, dusfcg, &
         dvsfcg, oc, oa4, clx, theta, gamma, sigma, elvmax, dp, dzgrid, varss, ocss, oa4ss, clxss,      &
         dycore_active, dycore_fv3, errmsg, errflg)

      use machine, only: kind_phys

      implicit none

      ! interface variables
      integer,              intent(in )                   :: im, levs, ntrac
      real(kind=kind_phys), intent(in )                   :: dtf, dtp, dxmin, dxinv
      real(kind=kind_phys), intent(in ), dimension(:)     :: slmsk, area, pgr
      integer,              intent(in )                   :: dycore_active, dycore_fv3

      integer,              intent(out), dimension(:)     :: islmsk
      real(kind=kind_phys), intent(out), dimension(:)     :: work1, work2, psurf
      real(kind=kind_phys), intent(out), dimension(:,:)   :: dudt, dvdt, dtdt
      real(kind=kind_phys), intent(out), dimension(:,:,:) :: dqdt
      real(kind=kind_phys), intent(out), dimension(:,:)   :: ten_u, ten_v, ten_t
      real(kind=kind_phys), intent(out), dimension(:,:,:) :: ten_q
      integer,              intent(out), dimension(:)     :: kpbl
      real(kind=kind_phys), intent(out), dimension(:)     :: dusfcg, dvsfcg
      ! Orography statistics used by the GWD schemes (filled from hprime each step in FV3)
      real(kind=kind_phys), intent(inout) :: oc(:), oa4(:,:), clx(:,:)
      real(kind=kind_phys), intent(inout) :: theta(:), gamma(:), sigma(:), elvmax(:)
      ! Layer pressure and height thickness (formerly Interstitial%del and %deltaZ)
      real(kind=kind_phys), intent(inout) :: dp(:,:), dzgrid(:,:)
      ! Small-scale orography statistics (allocated only for gwd_opt = 2, 3, 22, 33)
      real(kind=kind_phys), intent(inout), optional :: varss(:), ocss(:), oa4ss(:,:), clxss(:,:)

      character(len=*),     intent(out)                   :: errmsg
      integer,              intent(out)                   :: errflg

      ! local variables
      real(kind=kind_phys), parameter   :: zero = 0.0_kind_phys, one = 1.0_kind_phys
      integer :: i, k, n

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      do i = 1, im
        islmsk(i)   = nint(slmsk(i))

        work1(i) = (log(area(i)) - dxmin) * dxinv
        work1(i) = max(zero, min(one, work1(i)))
        work2(i) = one - work1(i)
        psurf(i) = pgr(i)
      end do

      ! Reset accumulated tendencies
      do k=1,levs
        do i=1,im
          dudt(i,k)  = zero
          dvdt(i,k)  = zero
          dtdt(i,k)  = zero
        enddo
      enddo
      do n=1,ntrac
        do k=1,levs
          do i=1,im
            dqdt(i,k,n) = zero
          enddo
        enddo
     enddo

     ! Reset instantaeous tendencies
     do k=1,levs
        do i=1,im
          ten_u(i,k)  = zero
          ten_v(i,k)  = zero
          ten_t(i,k)  = zero
        enddo
      enddo
      do n=1,ntrac
        do k=1,levs
          do i=1,im
            ten_q(i,k,n) = zero
          enddo
        enddo
     enddo

     ! Reset instantaeous diagnostics
     kpbl   = 0
     dusfcg = zero
     dvsfcg = zero
     dp     = zero
     dzgrid = zero

     ! DJS2026:
     ! For UFS-FV3, these fields below were defined as physics interstitials
     ! and reset to the initialization value at the start of each physics timestep
     ! in GFS_GWD_generic_pre.F90.
     ! For UFS-MPAS, these fields were moved into a persistent data container and set
     ! once during initialization, ufs_mpas_gwd_to_physics.
     ! GFS_GWD_generic_pre.F90 is not called from MPAS based SDFs.
     !
     ! Both FV3 and MPAS use this interstitial, but we don't want to zero out these
     ! for MPAS, since GFS_GWD_generic_pre.F90 is not called for MPAS.
     !
     if (dycore_active == dycore_fv3) then
        oc     = zero
        oa4    = zero
        clx    = zero
        theta  = zero
        gamma  = zero
        sigma  = zero
        elvmax = zero
 
        if (present(varss)) varss = zero
        if (present(ocss))  ocss  = zero
        if (present(oa4ss)) oa4ss = zero
        if (present(clxss)) clxss = zero
     end if

    end subroutine GFS_suite_interstitial_1_run

  end module GFS_suite_interstitial_1
