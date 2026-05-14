      module progomega

        use mo_conv_kind, only : conv_wp

        implicit none

        public progomega_calc

      contains

!>\file progomega_calc.f90
!! This file contains the subroutine that calculates the prognostic
!! updraft vertical velocity that is used for closure computations in 
!! saSAS and C3 deep and shallow convection. 

!>\ingroup SAMFdeep
!>\ingroup SAMF_shal
!> This subroutine computes a prognostic updraft vertical velocity
!! used in the closure computations in the samfdeepcnv.f and cu_c3_conv.f scheme
!! This subroutine computes a prognostic updraft vertical velocity
!! used in the closure computations in the samfshalcnv. and cu_c3_shal scheme
!!\section gen_progomega progomega_calc General Algorithm
       
   subroutine progomega_calc(first_time_step,flag_restart,im,km,kbcon1,ktcon,omegain,delt,del, &
        zi,cnvflg,omegaout,grav,buo,drag,wush,bb1,bb2)
     
     use machine,  only : kind_phys
     use funcphys, only : fpvs  
     implicit none

     integer, intent(in)  :: im, km
     integer, intent(in)  :: kbcon1(im),ktcon(im)
     real(kind=conv_wp), intent(in)  :: delt,grav,bb1,bb2
     real(kind=conv_wp), intent(in)  :: omegain(im,km), del(im,km),zi(im,km)
     real(kind=conv_wp), intent(in)  :: drag(im,km),buo(im,km),wush(im,km)
     real(kind=conv_wp), intent(inout) :: omegaout(im,km)
     logical, intent(in)               :: cnvflg(im),first_time_step,flag_restart
     real(kind=conv_wp) :: termA(im,km),termB(im,km),termC(im,km),omega(im,km)
     real(kind=conv_wp) :: RHS(im,km),Kd(im,km)
     real(kind=conv_wp) :: dp,dz,discr,wush_pa,lbb1,lbb2,lbb3
     integer              :: i,k

     lbb1  = 1.5_conv_wp
     lbb2  = 0.6_conv_wp
     lbb3  = 1.2_conv_wp
     
     !Initialization 2D
     do k = 1,km
        do i = 1,im
           termA(i,k)=0.0_conv_wp
           termB(i,k)=0.0_conv_wp
           termC(i,k)=0.0_conv_wp
           RHS(i,k)=0.0_conv_wp
           omega(i,k)=omegain(i,k)
        enddo
     enddo

     do k = 1,km
        do i = 1,im
           if(cnvflg(i))then
              if(omega(i,k) < 1.0E-5) then
                 omega(i,k) = 0.
              endif
           endif
        enddo
     enddo
     
     if(first_time_step .and. .not. flag_restart)then
        do k = 1,km
           do i = 1,im
              if(cnvflg(i))then
                 omega(i,k)=-1.2_conv_wp !Pa/s 
              endif
           enddo
        enddo
     endif
     
     ! Compute RHS terms
     !Lisa Bengtsson: !  compute updraft velocity omega (Pa/s)
     !> - Expand the steady state solution of updraft velocity from Han et al.'s (2017)
     !> \cite han_et_al_2017 equation 7 to include the time-derivative, and an aerodynamic
     !> drag term from Gueremy 2016.
     !> Solve using implicit time-stepping scheme, solving the quadratic equation for omega. 
     
     do k = 2, km
        do i = 1, im
           if (cnvflg(i)) then
              if (k >= kbcon1(i) .and. k < ktcon(i)) then

                 ! Scale by dp/dz to have equation in Pa/s
                 !(dp/dz > 0)
                 dp = 1000.0_conv_wp * del(i,k)
                 dz = zi(i,k+1) - zi(i,k)
                 
                 !termA	- Ensures quadratic damping (drag).
                 !termB	- Ensures linear damping from wind shear.
                 !termC - Adds buoyancy forcing 
                 
                 !Coefficients for the quadratic equation
                 termA(i,k) = delt * ((lbb1 * drag(i,k) * (dp/dz)))
                 termB(i,k) = 1.0_conv_wp - delt * lbb3 * wush(i,k) * dp/dz
                 termC(i,k) = omega(i,k) - delt * lbb2 * buo(i,k) * (dp/dz) &
                      - delt * omega(i,k) * (omega(i,k-1) - omega(i,k)) / dp
                 !Compute the discriminant
                 discr = termB(i,k)**2 - 4.0_conv_wp * termA(i,k) * termC(i,k)

                 ! Check if discriminant is non-negative
                 if (discr >= 0.0_conv_wp) then
                 ! Solve quadratic equation, take the negative root
                 omegaout(i,k) = (-termB(i,k) - sqrt(discr)) / (2.0_conv_wp * termA(i,k))
                 else
                 omegaout(i,k) = omega(i,k)
                 endif

                 omegaout(i,k) = MAX(MIN(omegaout(i,k), -1.2_conv_wp), -80.0_conv_wp)
                
              endif
           endif
        enddo
     enddo
     
    end subroutine progomega_calc
end module progomega
