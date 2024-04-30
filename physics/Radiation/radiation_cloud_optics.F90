module module_radiation_cloud_optics
  use machine,      only: kind_phys
  implicit none
  
  real (kind_phys), parameter :: &
       cld_limit_lower = 0.001, &
       cld_limit_ovcst = 1.0 - 1.0e-8, &
       reliq_def  = 10.0 ,      & ! Default liq radius, in stratiform cloud, to 10 micron (used when effr_in=F)
       reice_def  = 50.0,       & ! Default ice radius, in stratiform cloud, to 50 micron (used when effr_in=F)
       reliqcnv_def = 10.0 ,    & ! Default liq radius, in convective cloud, to 10 micron (used when effr_in=F)
       reicecnv_def = 50.0,     & ! Default ice radius, in convective cloud, to 50 micron (used when effr_in=F)
       rerain_def = 1000.0,     & ! Default rain radius to 1000 micron (used when effr_in=F)
       resnow_def = 250.0,      & ! Default snow radius to 250 micron (used when effr_in=F)
       reice_min  = 10.0,       & ! Minimum ice size allowed by GFDL MP scheme
       reice_max  = 150.0         ! Maximum ice size allowed by GFDL MP scheme
contains
  
!> \ingroup radiation_cloud_optics
!! Compute cloud radiative properties for SAMF convective cloud scheme.
!!
!! - The total-cloud convective mixing-ratio is partitioned by phase into liquid/ice 
!!   cloud properties. LWP and IWP are computed.
!!
!! - The liquid and ice cloud effective particle sizes are assigned reference values.
!!   (*NOTE* STUB in place to expand this using "cmp_Re")
!!
!! - If cmp_XuRndl = True, the convective cloud-fraction is computed using Xu-Randall (1996).
!!   Otherwise, the cloud-fraction provided by the convection scheme is unperturbed.
!!
!! \section cloud_mp_SAMF_gen General Algorithm
  subroutine cloud_mp_SAMF(cmp_XuRndl, cmp_Re, nCol, nLev, t_lay, p_lev, p_lay, qs_lay,  &
       relhum, cnv_mixratio, con_ttp, con_g, alpha0,                                     &
       cld_cnv_lwp, cld_cnv_reliq, cld_cnv_iwp, cld_cnv_reice, cld_cnv_frac)
    implicit none

    ! Inputs
    logical, intent(in)    :: &
         cmp_XuRndl,    & ! Compute convective cloud fraction using Xu-Randall?
         cmp_Re           ! Compute liquid/ice particle sizes using ?????
    integer, intent(in)    :: &
         nCol,          & ! Number of horizontal grid points
         nLev             ! Number of vertical layers
    real(kind_phys), intent(in) :: &
         con_g,         & ! Physical constant: gravity         (m s-2)
         con_ttp,       & ! Triple point temperature of water  (K)
         alpha0           ! Parameter for Xu-Randall scheme.   (-)
    real(kind_phys), dimension(:,:),intent(in) :: &
         t_lay,         & ! Temperature at layer-centers       (K)
         p_lev,         & ! Pressure at layer-interfaces       (Pa)
         p_lay,         & ! Presure at layer-centers           (Pa)
         qs_lay,        & ! Specific-humidity at layer-centers (kg/kg)
         relhum,        & ! Relative-humidity                  (1)
         cnv_mixratio     ! Convective cloud mixing-ratio      (kg/kg)
    ! Outputs
    real(kind_phys), dimension(:,:),intent(out) :: &
         cld_cnv_lwp,   & ! Convective cloud liquid water path
         cld_cnv_reliq, & ! Convective cloud liquid effective radius
         cld_cnv_iwp,   & ! Convective cloud ice water path
         cld_cnv_reice    ! Convective cloud ice effecive radius
    real(kind_phys), dimension(:,:),intent(inout) :: &
         cld_cnv_frac     ! Convective cloud-fraction
    ! Local
    integer :: iCol, iLay
    real(kind_phys) :: tem0, tem1, deltaP, clwc

    tem0 = 1.0e5/con_g
    do iLay = 1, nLev
       do iCol = 1, nCol
          if (cnv_mixratio(iCol,iLay) > 0._kind_phys) then
             ! Partition water paths by phase.
             tem1   = min(1.0, max(0.0, (con_ttp-t_lay(iCol,iLay))*0.05))
             deltaP = abs(p_lev(iCol,iLay+1)-p_lev(iCol,iLay))*0.01
             clwc   = max(0.0, cnv_mixratio(iCol,iLay)) * tem0 * deltaP
             cld_cnv_iwp(iCol,iLay)   = clwc * tem1
             cld_cnv_lwp(iCol,iLay)   = clwc - cld_cnv_iwp(iCol,iLay)
             
             ! Assign particles size(s).
             if (cmp_Re) then
                ! do something here a bit more fancy?
             else
                ! Assume default liquid/ice effective radius (microns)
                cld_cnv_reliq(iCol,iLay) = reliqcnv_def
                cld_cnv_reice(iCol,iLay) = reicecnv_def
             endif

             ! Recompute cloud-fraction using Xu-Randall (1996)?
             if (cmp_XuRndl) then
                cld_cnv_frac(iCol,iLay) = cld_frac_XuRandall(p_lay(iCol,iLay),           &
                     qs_lay(iCol,iLay), relhum(iCol,iLay), cnv_mixratio(iCol,iLay), alpha0)
             else
                ! Otherwise, cloud-fraction from convection scheme will pass through and
                ! be used by the radiation.
                !cld_cnv_frac(iCol,iLay) = 1._kind_phys
             endif
          endif ! No juice.
       enddo    ! Columns
    enddo       ! Layers

  end subroutine cloud_mp_SAMF

!> \ingroup radiation_cloud_optics
!! This function computes the cloud-fraction following.
!! Xu-Randall(1996) A Semiempirical Cloudiness Parameterization for Use in Climate Models
!! https://doi.org/10.1175/1520-0469(1996)053<3084:ASCPFU>2.0.CO;2
!!
!! cld_frac = {1-exp[-alpha*cld_mr/((1-relhum)*qs_lay)**lambda]}*relhum**P
!!
!! \section cld_frac_XuRandall_gen General Algorithm
  function cld_frac_XuRandall(p_lay, qs_lay, relhum, cld_mr, alpha)
    implicit none
    ! Inputs
    real(kind_phys), intent(in) :: &
       p_lay,    & ! Pressure (Pa)
       qs_lay,   & ! Saturation vapor-pressure (Pa)
       relhum,   & ! Relative humidity
       cld_mr,   & ! Total cloud mixing ratio
       alpha       ! Scheme parameter (default=100)

    ! Outputs
    real(kind_phys) :: cld_frac_XuRandall

    ! Locals
    real(kind_phys) :: clwt, clwm, onemrh, tem1, tem2, tem3

    ! Parameters
    real(kind_phys) :: &
       lambda = 0.50, & !
       P      = 0.25

    clwt = 1.0e-6 * (p_lay*0.001)
    if (cld_mr > clwt) then
       onemrh = max(1.e-10, 1.0 - relhum)
       tem1   = alpha / min(max((onemrh*qs_lay)**lambda,0.0001),1.0)
       tem2   = max(min(tem1*(cld_mr - clwt), 50.0 ), 0.0 )
       tem3   = sqrt(sqrt(relhum)) ! This assumes "p" = 0.25. Identical, but cheaper than relhum**p
       !
       cld_frac_XuRandall = max( tem3*(1.0-exp(-tem2)), 0.0 )
    else
       cld_frac_XuRandall = 0.0
    endif

    return
  end function cld_frac_XuRandall
  
end module module_radiation_cloud_optics
