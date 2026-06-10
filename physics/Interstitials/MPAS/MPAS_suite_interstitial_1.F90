!> \file MPAS_suite_interstitial_1.f90
!!  Contains code to initialize process-split state tendencies

    module MPAS_suite_interstitial_1

    contains

!> \section arg_table_MPAS_suite_interstitial_1_run Argument Table
!! \htmlinclude MPAS_suite_interstitial_1_run.html
!!
    subroutine MPAS_suite_interstitial_1_run (im, levs, ntrac, &
      dudt, dvdt, dtdt, dqdt, errmsg, errflg)

      use machine, only: kind_phys

      implicit none

      ! interface variables
      integer,              intent(in )                   :: im, levs, ntrac

      real(kind=kind_phys), intent(out), dimension(:,:)   :: dudt, dvdt, dtdt
      real(kind=kind_phys), intent(out), dimension(:,:,:) :: dqdt
      
      character(len=*),     intent(out)                   :: errmsg
      integer,              intent(out)                   :: errflg

      ! local variables
      real(kind=kind_phys), parameter   :: zero = 0.0_kind_phys
      integer :: i, k, n

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

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

    end subroutine MPAS_suite_interstitial_1_run

  end module MPAS_suite_interstitial_1