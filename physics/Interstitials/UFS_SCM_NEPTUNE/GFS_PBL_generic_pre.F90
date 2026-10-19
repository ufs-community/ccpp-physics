!> \file GFS_PBL_generic_pre.F90
!!  Contains code related to PBL schemes to be called prior to PBL schemes within GFS-based physics suites.

      module GFS_PBL_generic_pre

      contains

!> \brief This scheme sets up the vertically diffused tracer array for any PBL scheme based on the microphysics scheme chosen
!! \section arg_table_GFS_PBL_generic_pre_run Argument Table
!! \htmlinclude GFS_PBL_generic_pre_run.html
!!
      subroutine GFS_PBL_generic_pre_run (im, levs, nvdiff, ntrac, rtg_ozone_index,      &
        ntqv, ntcw, ntiw, ntrw, ntsw, ntlnc, ntinc, ntrnc, ntsnc, ntgnc,                 &
        ntwa, ntia, ntgl, ntoz, ntke, ntkev, nqrimef, trans_aero, ntchs, ntchm,          &
        ntccn, nthl, nthnc, ntgv, nthv, ntrz, ntgz, nthz,                                &
        imp_physics, imp_physics_gfdl, imp_physics_thompson, imp_physics_wsm6,           &
        imp_physics_tempo, lthailaware, imp_physics_mg, imp_physics_fer_hires,           &
        imp_physics_nssl, ltaerosol, mraerosol, nssl_ccn_on, nssl_hail_on, nssl_3moment, &
        hybedmf, do_shoc, satmedmf, sa3dtke, qgrs, vdftra, def_4, def4q,                 &
        ugrs, vgrs, tgrs, errmsg, errflg)
        
      use machine,                only : kind_phys
      use GFS_PBL_generic_common, only : set_aerosol_tracer_index

      implicit none

      integer, parameter  :: kp = kind_phys
      integer, intent(out) :: rtg_ozone_index
      integer, intent(in) :: im, levs, nvdiff, ntrac
      integer, intent(in) :: ntqv, ntcw, ntiw, ntrw, ntsw, ntlnc, ntinc, ntrnc, ntsnc, ntgnc
      integer, intent(in) :: ntwa, ntia, ntgl, ntoz, ntke, ntkev, nqrimef,ntchs, ntchm
      integer, intent(in) :: ntccn, nthl, nthnc, ntgv, nthv, ntrz, ntgz, nthz
      logical, intent(in) :: trans_aero
      integer, intent(in) :: imp_physics, imp_physics_gfdl, imp_physics_thompson, imp_physics_wsm6, imp_physics_tempo
      integer, intent(in) :: imp_physics_mg, imp_physics_fer_hires
      logical, intent(in) :: ltaerosol, hybedmf, do_shoc, satmedmf, mraerosol, lthailaware
      logical, intent(in) :: sa3dtke !flag for SA-3D-TKE scheme
      integer, intent(in) :: imp_physics_nssl
      logical, intent(in) :: nssl_hail_on, nssl_ccn_on, nssl_3moment

      real(kind=kind_phys), dimension(:,:,:), intent(in) :: qgrs, def_4
      real(kind=kind_phys), dimension(:,:), intent(in) :: ugrs, vgrs, tgrs
      real(kind=kind_phys), dimension(:,:,:), intent(inout) :: vdftra, def4q

      ! CCPP error handling variables
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      real (kind=kind_phys), parameter :: zero = 0.0_kp, one=1.0_kp

      ! Local variables
      integer :: i, k, kk, k1, n

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      rtg_ozone_index=-1
!DH: dvdftra is only used if nvdiff != ntrac or (nvdiff == ntrac .and. )
      if (nvdiff == ntrac .and. (hybedmf .or. do_shoc .or. satmedmf)) then
        vdftra = qgrs
        if(sa3dtke) then
          def4q = def_4
        endif
        rtg_ozone_index = ntoz
      else
        if (imp_physics == imp_physics_wsm6) then
  ! WSM6
          do k=1,levs
            do i=1,im
              vdftra(i,k,1) = qgrs(i,k,ntqv)
              vdftra(i,k,2) = qgrs(i,k,ntcw)
              vdftra(i,k,3) = qgrs(i,k,ntiw)
              vdftra(i,k,4) = qgrs(i,k,ntoz)
              if(sa3dtke) then
                def4q(i,k,1) = def_4(i,k,ntqv)
                def4q(i,k,2) = def_4(i,k,ntcw)
                def4q(i,k,3) = def_4(i,k,ntiw)
                def4q(i,k,4) = def_4(i,k,ntoz)
              endif
            enddo
          enddo
          rtg_ozone_index = 4

  ! Ferrier-Aligo
        elseif (imp_physics == imp_physics_fer_hires) then
          do k=1,levs
            do i=1,im
              vdftra(i,k,1) = qgrs(i,k,ntqv)
              vdftra(i,k,2) = qgrs(i,k,ntcw)
              vdftra(i,k,3) = qgrs(i,k,ntiw)
              vdftra(i,k,4) = qgrs(i,k,ntrw)
              vdftra(i,k,5) = qgrs(i,k,nqrimef)
              vdftra(i,k,6) = qgrs(i,k,ntoz)
              if(sa3dtke) then
                def4q(i,k,1) = def_4(i,k,ntqv)
                def4q(i,k,2) = def_4(i,k,ntcw)
                def4q(i,k,3) = def_4(i,k,ntiw)
                def4q(i,k,4) = def_4(i,k,ntrw)
                def4q(i,k,5) = def_4(i,k,nqrimef)
                def4q(i,k,6) = def_4(i,k,ntoz)
              endif
            enddo
          enddo
          rtg_ozone_index = 6
        
        elseif (imp_physics == imp_physics_thompson) then
  ! Thompson
          if(ltaerosol) then
            do k=1,levs
              do i=1,im
                vdftra(i,k,1)  = qgrs(i,k,ntqv)
                vdftra(i,k,2)  = qgrs(i,k,ntcw)
                vdftra(i,k,3)  = qgrs(i,k,ntiw)
                vdftra(i,k,4)  = qgrs(i,k,ntrw)
                vdftra(i,k,5)  = qgrs(i,k,ntsw)
                vdftra(i,k,6)  = qgrs(i,k,ntgl)
                vdftra(i,k,7)  = qgrs(i,k,ntlnc)
                vdftra(i,k,8)  = qgrs(i,k,ntinc)
                vdftra(i,k,9)  = qgrs(i,k,ntrnc)
                vdftra(i,k,10) = qgrs(i,k,ntoz)
                vdftra(i,k,11) = qgrs(i,k,ntwa)
                vdftra(i,k,12) = qgrs(i,k,ntia)
                if(sa3dtke) then
                  def4q(i,k,1)  = def_4(i,k,ntqv)
                  def4q(i,k,2)  = def_4(i,k,ntcw)
                  def4q(i,k,3)  = def_4(i,k,ntiw)
                  def4q(i,k,4)  = def_4(i,k,ntrw)
                  def4q(i,k,5)  = def_4(i,k,ntsw)
                  def4q(i,k,6)  = def_4(i,k,ntgl)
                  def4q(i,k,7)  = def_4(i,k,ntlnc)
                  def4q(i,k,8)  = def_4(i,k,ntinc)
                  def4q(i,k,9)  = def_4(i,k,ntrnc)
                  def4q(i,k,10) = def_4(i,k,ntoz)
                  def4q(i,k,11) = def_4(i,k,ntwa)
                  def4q(i,k,12) = def_4(i,k,ntia)
                endif
              enddo
            enddo
            rtg_ozone_index = 10
          elseif(mraerosol) then
            do k=1,levs
              do i=1,im
                vdftra(i,k,1)  = qgrs(i,k,ntqv)
                vdftra(i,k,2)  = qgrs(i,k,ntcw)
                vdftra(i,k,3)  = qgrs(i,k,ntiw)
                vdftra(i,k,4)  = qgrs(i,k,ntrw)
                vdftra(i,k,5)  = qgrs(i,k,ntsw)
                vdftra(i,k,6)  = qgrs(i,k,ntgl)
                vdftra(i,k,7)  = qgrs(i,k,ntlnc)
                vdftra(i,k,8)  = qgrs(i,k,ntinc)
                vdftra(i,k,9)  = qgrs(i,k,ntrnc)
                vdftra(i,k,10) = qgrs(i,k,ntoz)
                if(sa3dtke) then
                  def4q(i,k,1)  = def_4(i,k,ntqv)
                  def4q(i,k,1)  = def_4(i,k,ntqv)
                  def4q(i,k,2)  = def_4(i,k,ntcw)
                  def4q(i,k,3)  = def_4(i,k,ntiw)
                  def4q(i,k,4)  = def_4(i,k,ntrw)
                  def4q(i,k,5)  = def_4(i,k,ntsw)
                  def4q(i,k,6)  = def_4(i,k,ntgl)
                  def4q(i,k,7)  = def_4(i,k,ntlnc)
                  def4q(i,k,8)  = def_4(i,k,ntinc)
                  def4q(i,k,9)  = def_4(i,k,ntrnc)
                  def4q(i,k,10) = def_4(i,k,ntoz)
                endif
              enddo
            enddo
            rtg_ozone_index = 10
          else
            do k=1,levs
              do i=1,im
                vdftra(i,k,1) = qgrs(i,k,ntqv)
                vdftra(i,k,2) = qgrs(i,k,ntcw)
                vdftra(i,k,3) = qgrs(i,k,ntiw)
                vdftra(i,k,4) = qgrs(i,k,ntrw)
                vdftra(i,k,5) = qgrs(i,k,ntsw)
                vdftra(i,k,6) = qgrs(i,k,ntgl)
                vdftra(i,k,7) = qgrs(i,k,ntinc)
                vdftra(i,k,8) = qgrs(i,k,ntrnc)
                vdftra(i,k,9) = qgrs(i,k,ntoz)
                if(sa3dtke) then
                  def4q(i,k,1) = def_4(i,k,ntqv)
                  def4q(i,k,2) = def_4(i,k,ntcw)
                  def4q(i,k,3) = def_4(i,k,ntiw)
                  def4q(i,k,4) = def_4(i,k,ntrw)
                  def4q(i,k,5) = def_4(i,k,ntsw)
                  def4q(i,k,6) = def_4(i,k,ntgl)
                  def4q(i,k,7) = def_4(i,k,ntinc)
                  def4q(i,k,8) = def_4(i,k,ntrnc)
                  def4q(i,k,9) = def_4(i,k,ntoz)
                endif
              enddo
            enddo
            rtg_ozone_index = 9
          endif
        elseif (imp_physics == imp_physics_tempo) then
  ! Tempo
           do k=1,levs
              do i=1,im
                 vdftra(i,k,1)  = qgrs(i,k,ntqv)
                 vdftra(i,k,2)  = qgrs(i,k,ntcw)
                 vdftra(i,k,3)  = qgrs(i,k,ntiw)
                 vdftra(i,k,4)  = qgrs(i,k,ntrw)
                 vdftra(i,k,5)  = qgrs(i,k,ntsw)
                 vdftra(i,k,6)  = qgrs(i,k,ntgl)
                 vdftra(i,k,7)  = qgrs(i,k,ntinc)
                 vdftra(i,k,8)  = qgrs(i,k,ntrnc)
                 vdftra(i,k,9)  = qgrs(i,k,ntoz)
                 if(sa3dtke) then
                   def4q(i,k,1)  = def_4(i,k,ntqv)
                   def4q(i,k,2)  = def_4(i,k,ntcw)
                   def4q(i,k,3)  = def_4(i,k,ntiw)
                   def4q(i,k,4)  = def_4(i,k,ntrw)
                   def4q(i,k,5)  = def_4(i,k,ntsw)
                   def4q(i,k,6)  = def_4(i,k,ntgl)
                   def4q(i,k,7)  = def_4(i,k,ntinc)
                   def4q(i,k,8)  = def_4(i,k,ntrnc)
                   def4q(i,k,9)  = def_4(i,k,ntoz)
                 endif
              enddo
           enddo
           rtg_ozone_index = 9
           n = 10
           
           if (ltaerosol) then
              do k=1,levs
                 do i=1,im
                    vdftra(i,k,n) = qgrs(i,k,ntlnc)
                    vdftra(i,k,n+1) = qgrs(i,k,ntwa)
                    vdftra(i,k,n+2) = qgrs(i,k,ntia)
                    if(sa3dtke) then
                      def4q(i,k,n) = def_4(i,k,ntlnc)
                      def4q(i,k,n+1) = def_4(i,k,ntwa)
                      def4q(i,k,n+2) = def_4(i,k,ntia)
                    endif
                 enddo
              enddo
              n = 13
           endif

           if (lthailaware) then
              do k=1,levs
                 do i=1,im              
                    vdftra(i,k,n) = qgrs(i,k,ntgnc)
                    vdftra(i,k,n+1) = qgrs(i,k,ntgv)
                    if(sa3dtke) then
                      def4q(i,k,n) = def_4(i,k,ntgnc)
                      def4q(i,k,n+1) = def_4(i,k,ntgv)
                    endif
                 enddo
              enddo
           endif
  ! MG
        elseif (imp_physics == imp_physics_mg) then        ! MG3/2
          if (ntgl > 0) then                               ! MG3
            do k=1,levs
              do i=1,im
                vdftra(i,k,1)  = qgrs(i,k,ntqv)
                vdftra(i,k,2)  = qgrs(i,k,ntcw)
                vdftra(i,k,3)  = qgrs(i,k,ntiw)
                vdftra(i,k,4)  = qgrs(i,k,ntrw)
                vdftra(i,k,5)  = qgrs(i,k,ntsw)
                vdftra(i,k,6)  = qgrs(i,k,ntgl)
                vdftra(i,k,7)  = qgrs(i,k,ntlnc)
                vdftra(i,k,8)  = qgrs(i,k,ntinc)
                vdftra(i,k,9)  = qgrs(i,k,ntrnc)
                vdftra(i,k,10) = qgrs(i,k,ntsnc)
                vdftra(i,k,11) = qgrs(i,k,ntgnc)
                vdftra(i,k,12) = qgrs(i,k,ntoz)
                if(sa3dtke) then
                  def4q(i,k,1)  = def_4(i,k,ntqv)
                  def4q(i,k,2)  = def_4(i,k,ntcw)
                  def4q(i,k,3)  = def_4(i,k,ntiw)
                  def4q(i,k,4)  = def_4(i,k,ntrw)
                  def4q(i,k,5)  = def_4(i,k,ntsw)
                  def4q(i,k,6)  = def_4(i,k,ntgl)
                  def4q(i,k,7)  = def_4(i,k,ntlnc)
                  def4q(i,k,8)  = def_4(i,k,ntinc)
                  def4q(i,k,9)  = def_4(i,k,ntrnc)
                  def4q(i,k,10) = def_4(i,k,ntsnc)
                  def4q(i,k,11) = def_4(i,k,ntgnc)
                  def4q(i,k,12) = def_4(i,k,ntoz)
                endif
              enddo
            enddo
            rtg_ozone_index = 12
          else                                             ! MG2
            do k=1,levs
              do i=1,im
                vdftra(i,k,1)  = qgrs(i,k,ntqv)
                vdftra(i,k,2)  = qgrs(i,k,ntcw)
                vdftra(i,k,3)  = qgrs(i,k,ntiw)
                vdftra(i,k,4)  = qgrs(i,k,ntrw)
                vdftra(i,k,5)  = qgrs(i,k,ntsw)
                vdftra(i,k,6)  = qgrs(i,k,ntlnc)
                vdftra(i,k,7)  = qgrs(i,k,ntinc)
                vdftra(i,k,8)  = qgrs(i,k,ntrnc)
                vdftra(i,k,9)  = qgrs(i,k,ntsnc)
                vdftra(i,k,10) = qgrs(i,k,ntoz)
                if(sa3dtke) then
                  def4q(i,k,1)  = def_4(i,k,ntqv)
                  def4q(i,k,2)  = def_4(i,k,ntcw)
                  def4q(i,k,3)  = def_4(i,k,ntiw)
                  def4q(i,k,4)  = def_4(i,k,ntrw)
                  def4q(i,k,5)  = def_4(i,k,ntsw)
                  def4q(i,k,6)  = def_4(i,k,ntlnc)
                  def4q(i,k,7)  = def_4(i,k,ntinc)
                  def4q(i,k,8)  = def_4(i,k,ntrnc)
                  def4q(i,k,9)  = def_4(i,k,ntsnc)
                  def4q(i,k,10) = def_4(i,k,ntoz)
                endif
              enddo
            enddo
            rtg_ozone_index = 10
          endif
        elseif (imp_physics == imp_physics_gfdl) then
  ! GFDL MP
          do k=1,levs
            do i=1,im
              vdftra(i,k,1) = qgrs(i,k,ntqv)
              vdftra(i,k,2) = qgrs(i,k,ntcw)
              vdftra(i,k,3) = qgrs(i,k,ntiw)
              vdftra(i,k,4) = qgrs(i,k,ntrw)
              vdftra(i,k,5) = qgrs(i,k,ntsw)
              vdftra(i,k,6) = qgrs(i,k,ntgl)
              vdftra(i,k,7) = qgrs(i,k,ntoz)
              if(sa3dtke) then
                def4q(i,k,1) = def_4(i,k,ntqv)
                def4q(i,k,2) = def_4(i,k,ntcw)
                def4q(i,k,3) = def_4(i,k,ntiw)
                def4q(i,k,4) = def_4(i,k,ntrw)
                def4q(i,k,5) = def_4(i,k,ntsw)
                def4q(i,k,6) = def_4(i,k,ntgl)
                def4q(i,k,7) = def_4(i,k,ntoz)
              endif
            enddo
          enddo
          rtg_ozone_index = 7
        elseif (imp_physics == imp_physics_nssl ) then
  ! nssl
            IF ( nssl_hail_on ) THEN
            do k=1,levs
              do i=1,im
                vdftra(i,k,1)  = qgrs(i,k,ntqv)
                vdftra(i,k,2)  = qgrs(i,k,ntcw)
                vdftra(i,k,3)  = qgrs(i,k,ntiw)
                vdftra(i,k,4)  = qgrs(i,k,ntrw)
                vdftra(i,k,5)  = qgrs(i,k,ntsw)
                vdftra(i,k,6)  = qgrs(i,k,ntgl)
                vdftra(i,k,7)  = qgrs(i,k,nthl)
                vdftra(i,k,8)  = qgrs(i,k,ntlnc)
                vdftra(i,k,9)  = qgrs(i,k,ntinc)
                vdftra(i,k,10) = qgrs(i,k,ntrnc)
                vdftra(i,k,11) = qgrs(i,k,ntsnc)
                vdftra(i,k,12) = qgrs(i,k,ntgnc)
                vdftra(i,k,13) = qgrs(i,k,nthnc)
                vdftra(i,k,14) = qgrs(i,k,ntgv)
                vdftra(i,k,15) = qgrs(i,k,nthv)
                vdftra(i,k,16) = qgrs(i,k,ntoz)
                if(sa3dtke) then
                  def4q(i,k,1)  = def_4(i,k,ntqv)
                  def4q(i,k,2)  = def_4(i,k,ntcw)
                  def4q(i,k,3)  = def_4(i,k,ntiw)
                  def4q(i,k,4)  = def_4(i,k,ntrw)
                  def4q(i,k,5)  = def_4(i,k,ntsw)
                  def4q(i,k,6)  = def_4(i,k,ntgl)
                  def4q(i,k,7)  = def_4(i,k,nthl)
                  def4q(i,k,8)  = def_4(i,k,ntlnc)
                  def4q(i,k,9)  = def_4(i,k,ntinc)
                  def4q(i,k,10) = def_4(i,k,ntrnc)
                  def4q(i,k,11) = def_4(i,k,ntsnc)
                  def4q(i,k,12) = def_4(i,k,ntgnc)
                  def4q(i,k,13) = def_4(i,k,nthnc)
                  def4q(i,k,14) = def_4(i,k,ntgv)
                  def4q(i,k,15) = def_4(i,k,nthv)
                  def4q(i,k,16) = def_4(i,k,ntoz)
                endif
                n = 16
                IF ( nssl_ccn_on ) THEN
                 vdftra(i,k,n+1)  = qgrs(i,k,ntccn)
                 if(sa3dtke) then
                   def4q(i,k,n+1)  = def_4(i,k,ntccn)
                 endif
                 n = n+1
                ENDIF
                IF ( nssl_3moment ) THEN
                 vdftra(i,k,n+1)  = qgrs(i,k,ntrz)
                 vdftra(i,k,n+2)  = qgrs(i,k,ntgz)
                 vdftra(i,k,n+3)  = qgrs(i,k,nthz)
                 if(sa3dtke) then
                   def4q(i,k,n+1)  = def_4(i,k,ntrz)
                   def4q(i,k,n+2)  = def_4(i,k,ntgz)
                   def4q(i,k,n+3)  = def_4(i,k,nthz)
                 endif
                 n = n+3
                ENDIF
              enddo
            enddo
            
            ELSE
            ! no hail
            do k=1,levs
              do i=1,im
                vdftra(i,k,1)  = qgrs(i,k,ntqv)
                vdftra(i,k,2)  = qgrs(i,k,ntcw)
                vdftra(i,k,3)  = qgrs(i,k,ntiw)
                vdftra(i,k,4)  = qgrs(i,k,ntrw)
                vdftra(i,k,5)  = qgrs(i,k,ntsw)
                vdftra(i,k,6)  = qgrs(i,k,ntgl)
                vdftra(i,k,7)  = qgrs(i,k,ntlnc)
                vdftra(i,k,8)  = qgrs(i,k,ntinc)
                vdftra(i,k,9)  = qgrs(i,k,ntrnc)
                vdftra(i,k,10) = qgrs(i,k,ntsnc)
                vdftra(i,k,11) = qgrs(i,k,ntgnc)
                vdftra(i,k,12) = qgrs(i,k,ntgv)
                vdftra(i,k,13) = qgrs(i,k,ntoz)
                if(sa3dtke) then
                  def4q(i,k,1)  = def_4(i,k,ntqv)
                  def4q(i,k,2)  = def_4(i,k,ntcw)
                  def4q(i,k,3)  = def_4(i,k,ntiw)
                  def4q(i,k,4)  = def_4(i,k,ntrw)
                  def4q(i,k,5)  = def_4(i,k,ntsw)
                  def4q(i,k,6)  = def_4(i,k,ntgl)
                  def4q(i,k,7)  = def_4(i,k,ntlnc)
                  def4q(i,k,8)  = def_4(i,k,ntinc)
                  def4q(i,k,9)  = def_4(i,k,ntrnc)
                  def4q(i,k,10) = def_4(i,k,ntsnc)
                  def4q(i,k,11) = def_4(i,k,ntgnc)
                  def4q(i,k,12) = def_4(i,k,ntgv)
                  def4q(i,k,13) = def_4(i,k,ntoz)
                endif
                 n = 13
                IF ( nssl_ccn_on ) THEN
                 vdftra(i,k,n+1) = qgrs(i,k,ntccn)
                 if(sa3dtke) then
                   def4q(i,k,n+1) = def_4(i,k,ntccn)
                 endif
                 n = n+1
                ENDIF
                IF ( nssl_3moment ) THEN
                 vdftra(i,k,n+1) = qgrs(i,k,ntrz)
                 vdftra(i,k,n+2) = qgrs(i,k,ntgz)
                 if(sa3dtke) then
                   def4q(i,k,n+1) = def_4(i,k,ntrz)
                   def4q(i,k,n+2) = def_4(i,k,ntgz)
                 endif
                 n = n+2
                ENDIF
              enddo
            enddo
            
            ENDIF


        endif
!
        if (trans_aero) then
          call set_aerosol_tracer_index(imp_physics, imp_physics_wsm6,          &
                                        imp_physics_thompson, ltaerosol,mraerosol, &
                                        imp_physics_tempo, lthailaware, &
                                        imp_physics_mg, ntgl, imp_physics_gfdl, &
                                        imp_physics_nssl,                       &
                                        nssl_hail_on, nssl_ccn_on, nssl_3moment, kk, &
                                        errmsg, errflg)
          if (errflg /= 0) return
          !
          k1 = kk
          do n=ntchs,ntchm+ntchs-1
            k1 = k1 + 1
            do k=1,levs
              do i=1,im
                vdftra(i,k,k1) = qgrs(i,k,n)
                if(sa3dtke) then
                  def4q(i,k,k1) = def_4(i,k,n)
                endif
              enddo
            enddo
          enddo
        endif
!
        if (ntke>0) then
          do k=1,levs
            do i=1,im
              vdftra(i,k,ntkev) = qgrs(i,k,ntke)
              if(sa3dtke) then
                def4q(i,k,ntkev) = def_4(i,k,ntke)
              endif
            enddo
          enddo
        endif
!
      endif

    end subroutine GFS_PBL_generic_pre_run

    end module GFS_PBL_generic_pre
