!---------------------------------------------------------------------
!
! Copyright (C) 2001-2002 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Generate a pseudopotential in the Virtual Crystal Approximation:
!
!   V^{(vca)} = V_{loc)^{(vca)} + V_{nl}^{(vca)}
! where
!   V_{loc)^{(vca)} = x V_{loc}^{(1)} + (1-x) V_{loc}^{(2)}
! and
!   V_{nl)^{(vca)} = \sum_{ij} |\beta^{(1)}_i> x D^{(1)}_{ij} <\beta^{(1)}_j|
!                  + \sum_{ij} |\beta^{(2)}_i> (1-x)D^{(2)}_{ij} <\beta^{{2)}_j|
! where
!   V_{loc}^{(n)}(r) is the local part of pseudopot n
!   \beta^{{n)}_i(r) are the projectors for pseudopot n
!   D^{(n))_{ij} are the (bare) components of matrix D for pseudopot n
!
!
! virtual_v2: supports reading of UPF v2 format
! Author: Jingyang Wang (jw598@cornell.edu)
!
! virtual_v3: supports reading of UPF v2 format (SL, USPP, PAW and GIPAW)
! Author: By Student
!
PROGRAM virtual_test

  !USE pseudo_types, ONLY : pseudo_upf, nullify_pseudo_upf, &
  !                         deallocate_pseudo_upf
  USE upf_module, ONLY : read_upf
  USE emend_upf_module, ONLY: make_emended_upf_copy
  USE wrappers, ONLY: f_remove
  USE write_upf_module, ONLY : write_upf
  USE radial_grids, ONLY : radial_grid_type, nullify_radial_grid
  USE environment, ONLY: environment_start, environment_end
  USE mp_global, ONLY: mp_startup, mp_global_end
  USE io_global, ONLY: ionode, stdout
  USE pseudo_types, ONLY : paw_in_upf, pseudo_upf, &
                           nullify_paw_in_upf, nullify_pseudo_upf, &
                           deallocate_paw_in_upf, deallocate_pseudo_upf

  !
  IMPLICIT NONE
  !
  INTEGER :: is, iunps, ierr, ounps
  real(8) :: x
  CHARACTER(256) :: filein(2), fileout
  !
  !  Local variables
  !
  INTEGER :: ios
  TYPE (pseudo_upf) :: upf(2), upf_vca
  TYPE (radial_grid_type) :: grid(2)
  LOGICAL :: is_xml, exst 
#if defined(__MPI)
  CALL mp_startup()
#endif
  CALL environment_start('VIRTUAL_V3.X')
  IF (ionode) THEN
     !
     PRINT '('' '')'
     PRINT '('' Generate the UPF pseudopotential for a virtual atom '')'
     PRINT '('' combining two pseudopootentials in UPF format '')'
     PRINT '('' '')'
     !

     ! Read pseudopotentials
     !
     DO is=1,2

       PRINT '(''  Input PP file # '',i2,'' in UPF format > '',$)', is
       READ (5, '(a)', end = 20, err = 20) filein(is)

       !  nullify objects as soon as they are instantiated

       CALL nullify_pseudo_upf(upf(is))
       CALL nullify_paw_in_upf(upf(is)%paw)
       CALL nullify_radial_grid(grid(is))
       INQUIRE ( FILE = TRIM(filein(is)), EXIST = exst )  
       IF (.NOT. exst ) CALL errore ( 'virtual_v3.x: ', TRIM(filein(is)) // ' not found', 5)  
       ierr = 0
       CALL read_upf(upf(is), GRID = grid(is), IERR = ierr, FILENAME = trim(filein(is)))
       IF (ierr ==-81 ) THEN
          IF (ionode) is_xml = make_emended_upf_copy( trim(filein(is)), 'tmp.upf' )
          CALL  read_upf(upf(is), GRID = grid(is), IERR = ierr, FILENAME = 'tmp.upf' )
          IF (ionode) ios = f_remove('tmp.upf' )
       ENDIF
       IF ( ALL([0,-1,-2] /= ierr)) THEN
          PRINT *, ierr
          CALL errore('virtual_test', 'reading pseudo upf', ierr)
       ENDIF
       PRINT '('' '')'
       !
       !v2 message (old version)
       !IF ( TRIM(upf(is)%typ) == 'PAW') CALL errore('virtual_v2.x: ', &
       !                                              'Use of PAW is not implemented', 1) 
       !
     ENDDO
     ! CHECK Z-valence 
     IF ( upf(1)%zp /= upf(2)%zp ) WRITE (stdout, *) "CAUTION !!! "//& 
       "You are mixing pseudos with different number of electrons in valence"
     IF (upf(1)%lmax /= upf(2)%lmax ) WRITE ( stdout, *) "CAUTION !!! " //& 
      " You are mixing pseudos that  act on different angular momenta " 
     IF (upf(1)%lmax_rho /= upf(2)%lmax_rho ) WRITE ( stdout, *) "CAUTION !!! " //& 
      " You are mixing pseudos that  act on different angular momenta " 
     IF (upf(1)%paw%lmax_aug /= upf(2)%paw%lmax_aug ) WRITE ( stdout, *) "CAUTION !!! " //& 
      " You are mixing pseudos that  act on different angular momenta " 
     IF ( upf(1)%nbeta /= upf(2)%nbeta ) WRITE ( stdout, *) "CAUTION !!! " //&
      " You are mixing pseudos with a different number of projectors " 

     ! Choose mixing parameter x
     !
     PRINT '('' New Pseudo = x '',a,'' + (1-x) '',a)', (trim(filein(is)), is=1,2)
   10 CONTINUE
     WRITE(stdout,'('' mixing parameter x [0<x<1] = '')', advance="NO")
     READ (5,*) x
     IF (x<0.d0 .or. x>1)  GOTO 10

     ! compute virtual crystal approximation
     !
     CALL compute_virtual(x, filein, upf(1), upf_vca)

     ! write VCA pseudopotential to file
     !
     fileout='NewPseudo.UPF'
     PRINT '("Output PP file in UPF format :  ",a)', fileout
     !
     CALL write_upf ( trim(fileout), upf_vca, SCHEMA='v2' )
     !
     CLOSE(ounps)
     CALL deallocate_pseudo_upf(upf(1))
     CALL deallocate_pseudo_upf(upf(2))
     !     ----------------------------------------------------------
     WRITE (stdout,"('Pseudopotential successfully written')")
     WRITE (stdout,"('Please review the content of the PP_INFO fields')")
     WRITE (stdout,"('*** Please TEST BEFORE USING !!! ***')")
     WRITE (stdout,*) "Base UPF is upf(1) = ",filein(1)
     WRITE (stdout,*) "mesh data are got from upf(1) = ",filein(1)
     WRITE (stdout,"('*** Attention !!! SL(semi local) is development version ***')") 
     WRITE (stdout,"('*** Attention !!! USPP+GIPAW is development version ***')") 
     WRITE (stdout,"('*** Attention !!! PAW-as-GIPAW, PAW+GIPAW and PAW are development version ***')") 
     IF ( upf(1)%zp /= upf(2)%zp ) WRITE (stdout, *) "CAUTION !!! "//& 
       "You are mixing pseudos with different number of electrons in valence"
     IF (upf(1)%lmax /= upf(2)%lmax ) WRITE ( stdout, *) "CAUTION !!! " //& 
      " You are mixing pseudos that  act on different angular momenta " 
     IF (upf(1)%lmax_rho /= upf(2)%lmax_rho ) WRITE ( stdout, *) "CAUTION !!! " //& 
      " You are mixing pseudos that  act on different lmax_rho " 
     IF (upf(1)%paw%lmax_aug /= upf(2)%paw%lmax_aug ) WRITE ( stdout, *) "CAUTION !!! " //& 
      " You are mixing pseudos that  act on different lmax_aug " 
     IF ( upf(1)%nbeta /= upf(2)%nbeta ) WRITE ( stdout, *) "CAUTION !!! " //&
      " You are mixing pseudos with a different number of projectors " 
     IF ( upf(1)%paw%iraug /= upf(2)%paw%iraug ) WRITE ( stdout, *) "CAUTION !!! " //&
      " You are mixing pseudos with a different number of cutoff_r_index " 
     !     ----------------------------------------------------------
     !
   ENDIF
  CALL environment_end('VIRTUAL_V3.X')
#if defined(__MPI)
  CALL mp_global_end()
#endif

   STOP
20 CALL errore ('virtual.x', 'error reading pseudo file', 1)

END PROGRAM virtual_test

!
!---------------------------------------------------------------------
SUBROUTINE compute_virtual(x, filein, upf, upf_vca)

  USE pseudo_types, ONLY : paw_in_upf, pseudo_upf
  USE splinelib
  USE funct, ONLY : set_dft_from_name, get_iexch, get_icorr, get_igcx, get_igcc

  IMPLICIT NONE

  INTEGER :: i, j, ib, ijv, ijv2, prova
  real(8), INTENT(in) :: x
  CHARACTER (len=256) :: filein(2)
  TYPE (pseudo_upf), INTENT(in) :: upf(2)
  TYPE (pseudo_upf), INTENT(inout) :: upf_vca

  CHARACTER(len=9) :: day, hour
  CHARACTER (len=5) :: xlabel

  LOGICAL, EXTERNAL :: matches

  !
  ! All variables to be written into the UPF file
  ! (UPF = unified pseudopotential format, v.2)
  !

  ! pp_info
  INTEGER :: upf_rel
  real(8) :: upf_rcloc

  ! pp_header
  INTEGER :: upf_iexch, upf_icorr, upf_igcx, upf_igcc
  INTEGER :: upf_lmax, upf_mesh, upf_nbeta, upf_ntwfc
  INTEGER :: upf_lmax_rho, upf_paw_raug, upf_paw_iraug, upf_paw_lmax_aug
  INTEGER :: upf_gipaw_ncore_orbitals, upf_gipaw_wfs_nchannels
  INTEGER :: upf1_ind, upf2_ind, upf_ind
  INTEGER :: debag_flag
  real(8) :: upf_zp, upf_ecutrho, upf_ecutwfc, upf_etotps
  real(8) :: upf_xmin, upf_rmax, upf_dx, upf_zmesh
  real(8) :: upf_qqq_eps   ! qfunc is null if its norm is .lt. qqq_eps
  LOGICAL :: upf_nlcc
  real(8), ALLOCATABLE :: upf_ocw(:)
  CHARACTER(len=2), ALLOCATABLE :: upf_elsw(:)
  CHARACTER(len=12):: upf_paw_augshape      ! shape of augmentation charge
  INTEGER, ALLOCATABLE :: upf_lchiw(:)

  ! pp_mesh
  real(8), ALLOCATABLE :: upf_r(:), upf_rab(:)
  real(8), ALLOCATABLE :: upf_rho_atc(:)
  real(8) :: capel
  real(8), ALLOCATABLE :: aux1(:,:), aux2(:,:)
  LOGICAL :: interpolate

  ! pp_local
  real(8), ALLOCATABLE :: upf_vloc0(:)

  ! pp_nonlocal
  ! pp_beta
  real(8), ALLOCATABLE :: upf_betar(:,:)
  INTEGER, ALLOCATABLE :: upf_lll(:), upf_ikk2(:)
  ! pp_dij
  real(8), ALLOCATABLE :: upf_dion(:,:)
  ! pp_qij
  INTEGER :: l, l1, l2
  INTEGER :: upf_nqf, upf_nqlc
  real(8), ALLOCATABLE :: upf_rinner(:), upf_qqq(:,:), upf_qfunc(:,:), upf_qfuncl(:,:,:)
  ! pp_qfcoef
  real(8), ALLOCATABLE :: upf_qfcoef(:,:,:,:)
  !
  ! pp_pswfc
  real(8), ALLOCATABLE :: upf_chi(:,:)
  !
  ! pp_rhoatom
  real(8), ALLOCATABLE :: upf_rho_at(:)
  !
  ! pp_spin_orb
  real(8), ALLOCATABLE :: upf_jjj(:)
  real(8), ALLOCATABLE :: upf_nn(:)
  real(8), ALLOCATABLE :: upf_jchi(:)
  !
  ! SL (semilocal)
  ! pp_vnl
  real(8), ALLOCATABLE :: upf_vnl(:,:,:) ! vnl(i,l,s) = V(r_i)_{ls}
  !
  ! PAW
  ! pp_augmentation
  real(8), ALLOCATABLE :: upf_paw_augmom(:,:,:)
  !
  ! pp augfun
  !real(8), ALLOCATABLE :: upf_paw_augfun(:,:,:,:)
  !
  ! pp_ae_rho_atc (aeccharg)
  real(8), ALLOCATABLE :: upf_paw_ae_rho_atc(:)
  !
  ! pp_aewfc
  real(8), ALLOCATABLE :: upf_aewfc(:,:)
  !
  ! pp_pswfc
  real(8), ALLOCATABLE :: upf_pswfc(:,:)
  !
  ! pp_paw_pfunc
  real(8), ALLOCATABLE :: upf_paw_pfunc(:,:,:)
  !
  ! pp_paw_pfunc_rel
  real(8), ALLOCATABLE :: upf_paw_pfunc_rel(:,:,:)
  !
  ! pp_paw_ptfunc
  real(8), ALLOCATABLE :: upf_paw_ptfunc(:,:,:)
  !
  ! pp_paw_aewfc_rel
  real(8), ALLOCATABLE :: upf_paw_aewfc_rel(:,:)
  !
  ! pp_ae_vloc
  real(8), ALLOCATABLE :: upf_paw_ae_vloc(:)
  !
  ! pp_kdiff
  !real(8), ALLOCATABLE :: upf_paw_kdiff(:,:)
  !
  ! pp_sqr
  !real(8), ALLOCATABLE :: upf_paw_sqr(:)
  !
  ! GIPAW
  ! pp_gipaw_core_orbital_el
  CHARACTER(LEN=2), ALLOCATABLE :: upf_gipaw_core_orbital_el(:)
  !
  ! pp_gipaw_core_orbital_n
  real(8), ALLOCATABLE :: upf_gipaw_core_orbital_n(:)
  !
  ! pp_gipaw_core_orbital_l
  real(8), ALLOCATABLE :: upf_gipaw_core_orbital_l(:)
  !
  ! pp_gipaw_core_orbital
  real(8), ALLOCATABLE :: upf_gipaw_core_orbital(:,:)
  !
  ! pp_gipaw_vlocal_ae
  real(8), ALLOCATABLE :: upf_gipaw_vlocal_ae(:)
  !
  ! pp_gipaw_vlocal_ps
  real(8), ALLOCATABLE :: upf_gipaw_vlocal_ps(:)
  !
  ! upf_gipaw_wfs_el
  CHARACTER(LEN=2), ALLOCATABLE :: upf_gipaw_wfs_el(:)
  !
  ! upf_gipaw_wfs_n
  real(8), ALLOCATABLE :: upf_gipaw_wfs_n(:)
  !
  ! upf_gipaw_wfs_l
  real(8), ALLOCATABLE :: upf_gipaw_wfs_l(:)
  !
  ! upf_gipaw_wfs_ll
  real(8), ALLOCATABLE :: upf_gipaw_wfs_ll(:)
  !
  ! upf_gipaw_wfs_rcut
  real(8), ALLOCATABLE :: upf_gipaw_wfs_rcut(:)
  !
  ! upf_gipaw_wfs_rcutus
  real(8), ALLOCATABLE :: upf_gipaw_wfs_rcutus(:)
  !
  ! pp_gipaw_wfs_ae
  real(8), ALLOCATABLE :: upf_gipaw_wfs_ae(:,:)
  !
  ! pp_gipaw_wfs_ps
  real(8), ALLOCATABLE :: upf_gipaw_wfs_ps(:,:)

  interpolate = .false.

  upf_vca = upf(1)

  !pp_info
  !upf_rel = -1
  !upf_rcloc = 0.d0
  !

  !pp_header
  upf_vca%generated  = 'Generated using virtual.x code'
  upf_vca%author = 'virtual_v3'

  CALL date_and_tim(day, hour)
  upf_vca%date = trim(day)

  WRITE( xlabel, '(f5.3)' ) x
  upf_vca%comment    = 'Pseudo = x '//trim(filein(1))//&
                   ' + (1-x) '//trim(filein(2))//', with x='//xlabel
  upf_vca%psd = "Xx"
  upf_vca%typ = "NC"
  IF (matches(upf(1)%typ, "SL")) THEN
    upf_vca%typ = "SL"
  ELSE IF (matches(upf(1)%typ, "USPP") .and. matches(upf(2)%typ, "USPP") .and. &
    (upf(1)%has_gipaw .eqv. upf(2)%has_gipaw)) THEN
     upf_vca%typ = "USPP"
  ELSE IF (matches(upf(1)%typ, "PAW") .and. matches(upf(2)%typ, "PAW") .and.  &
    (upf(1)%has_gipaw .eqv. upf(2)%has_gipaw) .and. &
    (upf(1)%paw%augshape == upf(2)%paw%augshape) .and. &
    (upf(1)%tpawp .eqv. upf(2)%tpawp)) THEN
     upf_vca%typ = "PAW"
  ELSE
     CALL errore('virtual_v3.x: ', 'potential types are not match !!!', 1) 
  ENDIF

  CALL set_dft_from_name(upf(1)%dft)
  upf_iexch = get_iexch()
  upf_icorr = get_icorr()
  upf_igcx  = get_igcx()
  upf_igcc  = get_igcc()
  CALL set_dft_from_name(upf(2)%dft)
  IF (get_iexch()/=upf_iexch .or. get_icorr()/=upf_icorr .or. &
      get_igcx()/=upf_igcx .or. get_igcc()/=upf_igcc) &
      CALL errore ('virtual','conflicting DFT tionals',1)

  upf_lmax = max(upf(1)%lmax, upf(2)%lmax)

  IF ( upf(1)%mesh/=upf(2)%mesh ) THEN
     WRITE (*,*) " "
     WRITE (*,*) " pseudopotentials have different mesh "
     WRITE (*,*) "         upf(1)     upf(2)"
     WRITE (*,*) "mesh   ",upf(1)%mesh, upf(2)%mesh
     WRITE (*,*) "r(1st) ",upf(1)%r(1), upf(2)%r(1)
     WRITE (*,*) "r(last)",upf(1)%r(upf(1)%mesh), upf(2)%r(upf(2)%mesh)
     interpolate = .true.
  ENDIF

  upf_mesh = upf(1)%mesh
  upf_nbeta = upf(1)%nbeta+upf(2)%nbeta
  !upf_ntwfc = upf(1)%nwfc
  upf_nlcc  = upf(1)%nlcc.or.upf(2)%nlcc
  upf_ecutrho = upf(1)%ecutrho
  upf_ecutwfc = upf(1)%ecutwfc
  !upf_etotps  = upf(1)%etotps
  !
  upf_ntwfc = upf(1)%nwfc+upf(2)%nwfc
  upf_ecutrho = x * upf(1)%ecutrho + (1.d0-x) * upf(2)%ecutrho
  upf_ecutwfc = x * upf(1)%ecutwfc + (1.d0-x) * upf(2)%ecutwfc
  upf_etotps  = x * upf(1)%etotps + (1.d0-x) * upf(2)%etotps
  !
  
  !
  IF (matches(upf(1)%typ, "PAW") .or. matches(upf(2)%typ, "PAW")) THEN
     IF (upf(1)%paw%augshape /= upf(2)%paw%augshape) &
        CALL errore('virtual_v3.x: ', 'potential types (augshape) are not match !!!', 1) 
     upf_paw_augshape = upf(1)%paw%augshape
     upf_paw_lmax_aug = max(upf(1)%paw%lmax_aug, upf(2)%paw%lmax_aug)
     upf_paw_raug     = x * upf(1)%paw%raug + (1.d0-x) * upf(1)%paw%raug
     upf_paw_iraug    = x * upf(1)%paw%iraug + (1.d0-x) * upf(1)%paw%iraug
     upf_qqq_eps      = x * upf(1)%qqq_eps + (1.d0-x) * upf(2)%qqq_eps
     !
     WRITE (*,*) ""
     WRITE (*,*) "PAW parameters"
     WRITE (*,*) "shape          = augshape, ", upf(1)%paw%augshape, upf(1)%paw%augshape
     WRITE (*,*) "cutoff_r       = raug    , ", upf(1)%paw%raug, upf(2)%paw%raug
     WRITE (*,*) "cutoff_r_index = iraug   , ", upf(1)%paw%iraug, upf(2)%paw%iraug
     WRITE (*,*) "lmax_aug                 , ", upf(1)%paw%lmax_aug, upf(2)%paw%lmax_aug
     WRITE (*,*) "augmentation_epsilon = qqq_eps ", upf(1)%qqq_eps, upf(2)%qqq_eps
  ENDIF
  IF (upf_vca%has_gipaw) THEN
     upf_lmax_rho = max(upf(1)%lmax_rho, upf(2)%lmax_rho)
     upf_xmin     = upf(1)%xmin
     upf_rmax     = upf(1)%rmax
     upf_zmesh    = upf(1)%zmesh
     upf_dx       = upf(1)%dx
     !
     WRITE (*,*) ""
     WRITE (*,*) "GIPAW parameters"
     WRITE (*,*) "lmax_rho", upf(1)%lmax_rho, upf(2)%lmax_rho
     WRITE (*,*) "xmin    ", upf(1)%xmin, upf(2)%xmin
     WRITE (*,*) "rmax    ", upf(1)%rmax, upf(2)%rmax
     WRITE (*,*) "zmesh   ", upf(1)%zmesh, upf(2)%zmesh
     WRITE (*,*) "dx      ", upf(1)%dx, upf(2)%dx
  ENDIF

  !ALLOCATE( upf_ocw(upf_ntwfc), upf_elsw(upf_ntwfc), upf_lchiw(upf_ntwfc) )
  !upf_ocw(1:upf_ntwfc)  = upf(1)%oc(1:upf_ntwfc)
  !upf_elsw(1:upf_ntwfc) = upf(1)%els(1:upf_ntwfc)
  !upf_lchiw(1:upf_ntwfc) = upf(1)%lchi(1:upf_ntwfc)
  upf_zp    =  x * upf(1)%zp + (1.d0-x) * upf(2)%zp
  !
  !pp_mesh
  capel = 0.d0
  DO i=1,upf_mesh
     IF (i<=upf(2)%mesh) THEN
        capel = capel + abs(upf(1)%r(i)-upf(2)%r(i)) + abs(upf(1)%rab(i)-upf(2)%rab(i))
     ELSE
        capel = capel + abs(upf(1)%r(i)) + abs(upf(1)%rab(i))
     ENDIF
  ENDDO
  IF (capel>1.d-6) THEN
     WRITE (*,*) " "
     WRITE (*,*) " pseudopotentials have different mesh "
     WRITE (*,*) "capel = ", capel
     interpolate = .true.
  ENDIF

  WRITE (*,*) " "
  WRITE (*,*) "INTERPOLATE = ", interpolate
  IF (interpolate) ALLOCATE ( aux1(1,upf(1)%mesh), aux2(1,upf(2)%mesh) )

  ALLOCATE( upf_r(upf_mesh), upf_rab(upf_mesh) )
  upf_r(1:upf_mesh)   = upf(1)%r(1:upf_mesh)
  upf_rab(1:upf_mesh) = upf(1)%rab(1:upf_mesh)
  !

  !pp_nlcc
  ALLOCATE( upf_rho_atc(upf_mesh) )
  IF (interpolate) THEN
     WRITE (*,*) "interpolate rho_atc"
     aux2(1,1:upf(2)%mesh) = upf(2)%rho_atc(1:upf(2)%mesh)
     CALL dosplineint( upf(2)%r(1:upf(2)%mesh), aux2, upf_r(1:upf_mesh), aux1 )
     ! upf(2)%rho_atc(1:upf_mesh) = aux1(1,1:upf_mesh)
     upf_rho_atc(1:upf_mesh) =    x     * upf(1)%rho_atc(1:upf_mesh) + &
                               (1.d0-x) * aux1(1,1:upf_mesh)
     WRITE (*,*) " done"
  ELSE
     upf_rho_atc(1:upf_mesh) =    x     * upf(1)%rho_atc(1:upf_mesh) + &
                               (1.d0-x) * upf(2)%rho_atc(1:upf_mesh)
  ENDIF
  !

  !pp_local
  ALLOCATE( upf_vloc0(upf_mesh) )
  IF (interpolate) THEN
     WRITE (*,*) "interpolate vloc0"
     aux2(1,1:upf(2)%mesh) =  upf(2)%vloc(1:upf(2)%mesh)

     CALL dosplineint( upf(2)%r(1:upf(2)%mesh), aux2, upf_r(1:upf_mesh), aux1 )

     ! upf(2)%vloc(1:upf_mesh) = aux1(1,1:upf_mesh)

     ! Jivtesh - if the mesh of the first atom extends to a larger radius
     ! than the mesh of the second atom, then, for those radii that are
     ! greater than the maximum radius of the second atom, the local potential
     ! of the second atom is calculated using the expression
     ! v_local = (-2)*Z/r instead of using the extrapolated value.
     ! This is because, typically extrapolation leads to positive potentials.
     ! This is implemented in lines 240-242

     DO i=1,upf(1)%mesh
        IF ( upf(1)%r(i) > upf(2)%r(upf(2)%mesh) ) &
           ! upf(2)%vloc(i) = -(2.0*upf(2)%zp)/upf(1)%r(i)
           aux1(1,i) = -(2.0*upf(2)%zp)/upf(1)%r(i)
     ENDDO
     upf_vloc0(1:upf_mesh) =      x     * upf(1)%vloc(1:upf_mesh) + &
                               (1.d0-x) * aux1(1,1:upf_mesh)
     WRITE (*,*) " done"
  ELSE
     upf_vloc0(1:upf_mesh) =      x     * upf(1)%vloc(1:upf_mesh) + &
                               (1.d0-x) * upf(2)%vloc(1:upf_mesh)
  ENDIF
  !

  !pp_nonlocal
  !pp_beta
  ALLOCATE( upf_betar(upf_mesh, upf_nbeta), &
            upf_lll(upf_nbeta), upf_ikk2(upf_nbeta) )
  ib = 0
  DO i=1,upf(1)%nbeta
     ib  = ib + 1
     upf_betar(1:upf_mesh,ib) = upf(1)%beta(1:upf_mesh,i)
     upf_lll(ib)              = upf(1)%lll(i)
     upf_ikk2(ib)             = upf(1)%kbeta(i)
  ENDDO

  DO i=1,upf(2)%nbeta
     ib  = ib + 1
     IF (interpolate) THEN
        IF (i==1) WRITE (*,*) "interpolate betar"
        aux2(1,1:upf(2)%mesh) = upf(2)%beta(1:upf(2)%mesh,i)
        CALL dosplineint( upf(2)%r(1:upf(2)%mesh), aux2, upf_r(1:upf_mesh), aux1 )
        ! upf(2)%beta(1:upf_mesh,i) = aux1(1,1:upf_mesh)
        upf_betar(1:upf_mesh,ib) = aux1(1,1:upf_mesh)
     ELSE
        upf_betar(1:upf_mesh,ib) = upf(2)%beta(1:upf_mesh,i)
     ENDIF

     upf_lll(ib)              = upf(2)%lll(i)

     ! SdG - when the meshes of the two pseudo are different the ikk2 limits
     ! for the beta functions of the second one must be set properly
     ! This is done in lines 273-277
     IF (interpolate) THEN
        j = 1
        DO WHILE ( upf_r(j) < upf(2)%r(upf(2)%kbeta(i)) )
           j = j + 1
        ENDDO
        upf_ikk2(ib) = j
        IF (i==upf(2)%nbeta) WRITE (*,*) " done"
     ELSE
        upf_ikk2(ib) = upf(2)%kbeta(i)
     ENDIF

  ENDDO
  !

  WRITE (*,*) " "
  WRITE (*,*) "upf(1)%lll = ", upf(1)%lll
  WRITE (*,*) "upf(2)%lll = ", upf(2)%lll
  WRITE (*,*) " "

  !pp_dij
  ALLOCATE( upf_dion(upf_nbeta, upf_nbeta) )
  upf_dion(:,:) = 0.d0

  DO i=1,upf(1)%nbeta
     DO j=1,upf(1)%nbeta
        upf_dion(i,j) = x * upf(1)%dion(i,j)
     ENDDO
  ENDDO

  DO i=1,upf(2)%nbeta
     DO j=1,upf(2)%nbeta
        upf_dion(upf(1)%nbeta+i, upf(1)%nbeta+j) = (1.d0-x) * upf(2)%dion(i,j)
     ENDDO
  ENDDO
  !

  WRITE (*,*) "pp_dij completed."
  WRITE (*,*) " "

  IF (matches(upf_vca%typ, "USPP") .or. matches(upf_vca%typ, "PAW")) THEN

    !pp_qij
    IF (upf(1)%nqf/=upf(2)%nqf) &
        CALL errore ("Virtual","different nqf are not implemented (yet)", 1)

    IF (upf(1)%nqlc/=upf(2)%nqlc) &
        CALL errore ("Virtual","different nqlc are not implemented (yet)", 1)

    upf_nqf = upf(1)%nqf    ! number of Q coefficients
    upf_nqlc = upf(1)%nqlc  ! number of angular momenta in Q

    ALLOCATE( upf_rinner(upf_nqlc) )
    DO i=1,upf_nqlc
      IF(upf(1)%rinner(i)/=upf(2)%rinner(i)) &
         CALL errore("Virtual","different rinner are not implemented (yet)",i)
    ENDDO

    upf_rinner(1:upf_nqlc) = upf(1)%rinner(1:upf_nqlc)

    ALLOCATE( upf_qqq(upf_nbeta,upf_nbeta) )
    upf_qqq(:,:) = 0.d0
    IF ( upf(1)%q_with_l .neqv. upf(2)%q_with_l) &
       CALL errore ( 'virtual.x: ', 'Augmentation charges are written using uncompatible formats',1)
    IF( upf(1)%q_with_l ) THEN
       upf_vca%q_with_l = .true.
       ALLOCATE( upf_qfuncl(upf_mesh, upf_nbeta*(upf_nbeta+1)/2, 0:2*upf(1)%lmax) )
       upf_qfuncl(:,:,:) = 0.d0
    ELSE
       upf_vca%q_with_l = .false.
       ALLOCATE( upf_qfunc(upf_mesh, upf_nbeta*(upf_nbeta+1)/2) )
       upf_qfunc(:,:) = 0.d0
    ENDIF

    WRITE (*,*) "pp_qij"
    WRITE (*,*) " "

    DO i=1,upf(1)%nbeta
       DO j=i,upf(1)%nbeta
          ijv = j * (j-1)/2 + i
          upf_qqq(i,j) = x * upf(1)%qqq(i,j)


          IF( allocated(upf_qfuncl) ) THEN
             l1=upf(1)%lll(i)
             l2=upf(1)%lll(j)
             DO l=abs(l1-l2), l1+l2
                upf_qfuncl(1:upf_mesh,ijv,l) = x * upf(1)%qfuncl(1:upf_mesh,ijv,l)
             ENDDO
          ELSE
             upf_qfunc(1:upf_mesh,ijv) = x * upf(1)%qfunc(1:upf_mesh,ijv)
          ENDIF

       ENDDO
    ENDDO


    DO i=1,upf(2)%nbeta
       DO j=i,upf(2)%nbeta
          ijv = j * (j-1)/2 + i
          !ijv  = (upf(1)%nbeta+j)*(upf(1)%nbeta+j-1)/2 + i + upf(1)%nbeta
          upf_qqq(upf(1)%nbeta+i,upf(1)%nbeta+j) = (1.d0-x) * upf(2)%qqq(i,j)

          ijv2 = (upf(1)%nbeta+j) * (upf(1)%nbeta+j-1) / 2 + (upf(1)%nbeta+i)


          IF( allocated(upf_qfuncl) ) THEN
             l1=upf(2)%lll(i)
             l2=upf(2)%lll(j)
             DO l=abs(l1-l2), l1+l2
                IF (interpolate) THEN
                   IF (i==1 .and. j==1 .and. l==abs(l1-l2)) WRITE (*,*) "interpolate qfuncl"
                   aux2(1,1:upf(2)%mesh) = upf(2)%qfuncl(1:upf(2)%mesh,ijv,l)
                   CALL dosplineint( upf(2)%r(1:upf(2)%mesh), aux2, upf_r(1:upf_mesh), aux1 )
                   upf_qfuncl(1:upf_mesh, ijv2, l) = (1.d0-x) * aux1(1,1:upf_mesh)
                   IF (i==upf(2)%nbeta .and. j==upf(2)%nbeta .and. l==(l1+l2)) WRITE (*,*) " done"
                ELSE
                   upf_qfuncl(1:upf_mesh, ijv2, l) = (1.d0-x) * upf(2)%qfuncl(1:upf_mesh,ijv,l)

                   IF ((i==1) .and. (j==1)) THEN
                   ENDIF

                ENDIF
             ENDDO
          ELSE
             IF (interpolate) THEN
                 aux2(1,1:upf(2)%mesh) = upf(2)%qfunc(1:upf(2)%mesh,ijv)
                CALL dosplineint( upf(2)%r(1:upf(2)%mesh), aux2, upf_r(1:upf_mesh), aux1 )
                upf_qfunc(1:upf_mesh, ijv2) = (1.d0-x) * aux1(1,1:upf_mesh)
             ELSE
                upf_qfunc(1:upf_mesh, ijv2) = (1.d0-x) * upf(2)%qfunc(1:upf_mesh,ijv)
             ENDIF
          ENDIF
          ! ! upf_qfunc(1:upf_mesh, upf(1)%nbeta+i, upf(1)%nbeta+j) = (1.d0-x) * qfunc(1:upf_mesh,i,j,2)
       ENDDO
    ENDDO
    !

    !pp_qfcoef
    IF (upf_nqf/=0) THEN
       ALLOCATE( upf_qfcoef(upf_nqf,upf_nqlc,upf_nbeta,upf_nbeta) )
       upf_qfcoef(:,:,:,:) = 0.d0
       DO i=1,upf(1)%nbeta
          DO j=1,upf(1)%nbeta
             upf_qfcoef(1:upf_nqf,1:upf_nqlc,i,j) = &
                 x * upf(1)%qfcoef(1:upf_nqf,1:upf_nqlc,i,j)
          ENDDO
       ENDDO
       DO i=1,upf(2)%nbeta
          DO j=1,upf(2)%nbeta
             upf_qfcoef(1:upf_nqf,1:upf_nqlc,upf(1)%nbeta+i,upf(1)%nbeta+j) = &
                 (1.d0-x) * upf(2)%qfcoef(1:upf_nqf,1:upf_nqlc,i,j)
          ENDDO
       ENDDO
    ENDIF
    !
  ENDIF

  !pp_pswfc
  ALLOCATE ( upf_chi(upf_mesh,upf_ntwfc) )
  upf_chi(1:upf_mesh,1:upf(1)%nwfc) = x * upf(1)%chi(1:upf_mesh,1:upf(1)%nwfc)
  DO i=1,upf(2)%nwfc
     IF (interpolate) THEN
        IF (i==1) WRITE (*,*) "interpolate chi"
        aux2(1,1:upf(2)%mesh) = upf(2)%chi(1:upf(2)%mesh,i)
        CALL dosplineint( upf(2)%r(1:upf(2)%mesh), aux2, upf_r(1:upf_mesh), aux1 )
        upf_chi(1:upf_mesh,upf(1)%nwfc+i) = (1.d0-x) * aux1 (1,1:upf_mesh)
        IF (i==upf(2)%nwfc) WRITE (*,*) " done"
     ELSE
        upf_chi(1:upf_mesh,upf(1)%nwfc+i) = (1.d0-x) * upf(2)%chi(1:upf_mesh,i)
     ENDIF
  ENDDO
  !IF (upf(1)%nwfc==upf(2)%nwfc) THEN
  !   DO i=1,upf(2)%nwfc
  !      IF (interpolate) THEN
  !         WRITE (*,*) " interpolate chi"
  !         aux2(1,1:upf(2)%mesh) = upf(2)%chi(1:upf(2)%mesh,i)
  !         CALL dosplineint( upf(2)%r(1:upf(2)%mesh), aux2, upf_r(1:upf_mesh), aux1 )
  !         ! chi(1:upf_mesh,i,2) = aux1(1,1:upf_mesh)
  !         upf_chi(1:upf_mesh,i) = x     * upf(1)%chi(1:upf_mesh,i) + &
  !                              (1.d0-x) * aux1 (1,1:upf_mesh)
  !      ELSE
  !         upf_chi(1:upf_mesh,i) =    x     * upf(1)%chi(1:upf_mesh,i) + &
  !                                 (1.d0-x) * upf(2)%chi(1:upf_mesh,i)
  !      ENDIF
  !      ! Jivtesh - The wavefunctions are calculated to be the average of the
  !      ! wavefunctions of the two atoms - lines 365-366
  !   ENDDO
  !ELSE
  !   WRITE (*,*) "Number of wavefunctions not the same for the two pseudopotentials"
  !ENDIF
  !!upf_chi(1:upf_mesh,1:upf_ntwfc) = chi(1:upf_mesh,1:upf_ntwfc,1)
  !
  !pp_rhoatm
  ALLOCATE ( upf_rho_at(upf_mesh) )
  IF (interpolate) THEN
     WRITE (*,*) "interpolate rho_at"
     aux2(1,1:upf(2)%mesh) = upf(2)%rho_at(1:upf(2)%mesh)
     CALL dosplineint( upf(2)%r(1:upf(2)%mesh), aux2, upf_r(1:upf_mesh), aux1 )
     ! rho_at(1:upf_mesh,2) = aux1(1,1:upf_mesh)
     upf_rho_at(1:upf_mesh) =    x     * upf(1)%rho_at(1:upf_mesh) + &
                              (1.d0-x) * aux1(1,1:upf_mesh)
     WRITE (*,*) " done"
  ELSE
     upf_rho_at(1:upf_mesh) =    x     * upf(1)%rho_at(1:upf_mesh) + &
                              (1.d0-x) * upf(2)%rho_at(1:upf_mesh)
  ENDIF

  !pp_spin_orb
  IF (upf_vca%has_so) THEN
    WRITE (*,*) " this pseudopotential has spin orbit coupling!"
    ALLOCATE (upf_jjj(upf_nbeta))
    upf_jjj(1:upf(1)%nbeta) = upf(1)%jjj
    upf_jjj(upf(1)%nbeta+1:upf_nbeta) = upf(2)%jjj
    !
    ALLOCATE (upf_nn(upf_ntwfc))
    upf_nn(1:upf(1)%nwfc) = upf(1)%nn
    upf_nn(upf(1)%nwfc+1:upf_ntwfc) = upf(2)%nn
    !
    ALLOCATE (upf_jchi(upf_ntwfc))
    upf_jchi(1:upf(1)%nwfc) = upf(1)%jchi
    upf_jchi(upf(1)%nwfc+1:upf_ntwfc) = upf(2)%jchi
  ENDIF

  IF (matches(upf_vca%typ, "PAW")) THEN
     !
     !pp_augmentation
     ALLOCATE ( upf_paw_augmom(upf_nbeta,upf_nbeta,0:2*upf_lmax) )
     upf_paw_augmom(:,:,:) = 0.d0
     upf_paw_augmom(1:upf(1)%nbeta,1:upf(1)%nbeta,0:2*upf(1)%lmax) = &
              x * upf(1)%paw%augmom(1:upf(1)%nbeta,1:upf(1)%nbeta,0:2*upf(1)%lmax)
     upf_paw_augmom(upf(1)%nbeta+1:upf_nbeta,upf(1)%nbeta+1:upf_nbeta,0:2*upf(2)%lmax) = &
       (1.d0-x) * upf(2)%paw%augmom(1:upf(2)%nbeta,1:upf(2)%nbeta,0:2*upf(2)%lmax)
     !
     !pp_paw_ae_rho_atc
     ALLOCATE ( upf_paw_ae_rho_atc(upf_mesh) )
     IF (interpolate) THEN
        WRITE (*,*) "interpolate ae_rho_atc"
        aux2(1,1:upf(2)%mesh) = upf(2)%paw%ae_rho_atc(1:upf(2)%mesh)
        CALL dosplineint( upf(2)%r(1:upf(2)%mesh), aux2, upf_r(1:upf_mesh), aux1 )
        upf_paw_ae_rho_atc(1:upf_mesh) =  x   * upf(1)%paw%ae_rho_atc(1:upf_mesh) + &
                                     (1.d0-x) * aux1(1,1:upf_mesh)
        WRITE (*,*) " done"
     ELSE
        upf_paw_ae_rho_atc(1:upf_mesh) =  x   * upf(1)%paw%ae_rho_atc(1:upf_mesh) + &
                                     (1.d0-x) * upf(2)%paw%ae_rho_atc(1:upf_mesh)
     ENDIF
     !
     !pp_aewfc
     ALLOCATE ( upf_aewfc(upf_mesh,upf_nbeta) )
     upf_aewfc(:,:) = 0.d0
     upf_aewfc(1:upf_mesh,1:upf(1)%nbeta) = x * upf(1)%aewfc(1:upf_mesh,1:upf(1)%nbeta)
     DO i=1,upf(2)%nbeta
       IF (interpolate) THEN
          IF (i==1) WRITE (*,*) "interpolate aewfc"
          aux2(1,1:upf(2)%mesh) = upf(2)%aewfc(1:upf(2)%mesh,i)
          CALL dosplineint( upf(2)%r(1:upf(2)%mesh), aux2, upf_r(1:upf_mesh), aux1 )
          upf_aewfc(1:upf_mesh, upf(1)%nbeta+i) = (1.d0-x) * aux1(1,1:upf_mesh)
          IF (i==upf(2)%nbeta) WRITE (*,*) " done"
       ELSE
          upf_aewfc(1:upf_mesh, upf(1)%nbeta+i) = (1.d0-x) * upf(2)%aewfc(1:upf_mesh,i)
       ENDIF
     ENDDO
     !
     !pp_pswfc
     ALLOCATE ( upf_pswfc(upf_mesh,upf_nbeta) )
     upf_pswfc(:,:) = 0.d0
     upf_pswfc(1:upf_mesh,1:upf(1)%nbeta) = x * upf(1)%pswfc(1:upf_mesh,1:upf(1)%nbeta)
     DO i=1,upf(2)%nbeta
       IF (interpolate) THEN
          IF (i==1) WRITE (*,*) "interpolate pswfc"
          aux2(1,1:upf(2)%mesh) = upf(2)%pswfc(1:upf(2)%mesh,i)
          CALL dosplineint( upf(2)%r(1:upf(2)%mesh), aux2, upf_r(1:upf_mesh), aux1 )
          upf_pswfc(1:upf_mesh, upf(1)%nbeta+i) = (1.d0-x) * aux1(1,1:upf_mesh)
          IF (i==upf(2)%nbeta) WRITE (*,*) " done"
       ELSE
          upf_pswfc(1:upf_mesh, upf(1)%nbeta+i) = (1.d0-x) * upf(2)%pswfc(1:upf_mesh,i)
       ENDIF
     ENDDO
     !
     !pp_pfunc
     ALLOCATE( upf_paw_pfunc(upf_mesh, upf_nbeta, upf_nbeta) )
     upf_paw_pfunc(:,:,:) = 0.d0
     upf_paw_pfunc(1:upf_mesh,1:upf(1)%nbeta,1:upf(1)%nbeta) = &
             x  * upf(1)%paw%pfunc(1:upf_mesh,1:upf(1)%nbeta,1:upf(1)%nbeta)
     IF (interpolate) THEN
        WRITE (*,*) "interpolate pfunc"
        DO i=1,upf(2)%nbeta
            DO j=1,upf(2)%nbeta
               aux2(1,1:upf(2)%mesh) = upf(2)%paw%pfunc(1:upf(2)%mesh,i,j)
               CALL dosplineint( upf(2)%r(1:upf(2)%mesh), aux2, upf_r(1:upf_mesh), aux1 )
               upf_paw_pfunc(1:upf_mesh,upf(1)%nbeta+i,upf(1)%nbeta+j) = &
                  (1.d0-x) * aux1(1,1:upf_mesh)
            ENDDO
        ENDDO
        WRITE (*,*) " done"
     ENDIF
     !
     !pp_ptfunc
     ALLOCATE( upf_paw_ptfunc(upf_mesh, upf_nbeta, upf_nbeta) )
     upf_paw_ptfunc(:,:,:) = 0.d0
     upf_paw_ptfunc(1:upf_mesh,1:upf(1)%nbeta,1:upf(1)%nbeta) = &
             x  * upf(1)%paw%ptfunc(1:upf_mesh,1:upf(1)%nbeta,1:upf(1)%nbeta)
     IF (interpolate) THEN
        WRITE (*,*) "interpolate ptfunc"
        DO i=1,upf(2)%nbeta
            DO j=1,upf(2)%nbeta
               aux2(1,1:upf(2)%mesh) = upf(2)%paw%ptfunc(1:upf(2)%mesh,i,j)
               CALL dosplineint( upf(2)%r(1:upf(2)%mesh), aux2, upf_r(1:upf_mesh), aux1 )
               upf_paw_ptfunc(1:upf_mesh,upf(1)%nbeta+i,upf(1)%nbeta+j) = &
                  (1.d0-x) * aux1(1,1:upf_mesh)
            ENDDO
        ENDDO
        WRITE (*,*) " done"
     ENDIF
     !
     !pp_pfunc_rel
     IF (upf_vca%has_so) THEN
        ALLOCATE( upf_paw_pfunc_rel(upf_mesh, upf_nbeta, upf_nbeta) )
        upf_paw_pfunc_rel(:,:,:) = 0.d0
        upf_paw_pfunc_rel(1:upf_mesh,1:upf(1)%nbeta,1:upf(1)%nbeta) = &
                x  * upf(1)%paw%pfunc_rel(1:upf_mesh,1:upf(1)%nbeta,1:upf(1)%nbeta)
        IF (interpolate) THEN
           WRITE (*,*) "interpolate pfunc_rel"
           DO i=1,upf(2)%nbeta
              DO j=1,upf(2)%nbeta
                 aux2(1,1:upf(2)%mesh) = upf(2)%paw%pfunc_rel(1:upf(2)%mesh,i,j)
                 CALL dosplineint( upf(2)%r(1:upf(2)%mesh), aux2, upf_r(1:upf_mesh), aux1 )
                 upf_paw_pfunc_rel(1:upf_mesh,upf(1)%nbeta+i,upf(1)%nbeta+j) = &
                    (1.d0-x) * aux1(1,1:upf_mesh)
              ENDDO
           ENDDO
           WRITE (*,*) " done"
        ENDIF
        !
        !pp_aewfc_rel
        ALLOCATE ( upf_paw_aewfc_rel(upf_mesh,upf_nbeta) )
        upf_paw_aewfc_rel(:,:) = 0.d0
        upf_paw_aewfc_rel(1:upf_mesh,1:upf(1)%nbeta) = x * upf(1)%paw%aewfc_rel(1:upf_mesh,1:upf(1)%nbeta)
        DO i=1,upf(2)%nbeta
           IF (interpolate) THEN
              IF (i==1) WRITE (*,*) "interpolate aewfc_rel"
              aux2(1,1:upf(2)%mesh) = upf(2)%paw%aewfc_rel(1:upf(2)%mesh,i)
              CALL dosplineint( upf(2)%r(1:upf(2)%mesh), aux2, upf_r(1:upf_mesh), aux1 )
              upf_paw_aewfc_rel(1:upf_mesh, upf(1)%nbeta+i) = (1.d0-x) * aux1(1,1:upf_mesh)
              IF (i==upf(2)%nbeta) WRITE (*,*) " done"
           ELSE
              upf_paw_aewfc_rel(1:upf_mesh, upf(1)%nbeta+i) = (1.d0-x) * upf(2)%paw%aewfc_rel(1:upf_mesh,i)
           ENDIF
        ENDDO
     ENDIF
     !
     !pp_paw_ae_vloc
     ALLOCATE ( upf_paw_ae_vloc(upf_mesh) )
     IF (interpolate) THEN
        WRITE (*,*) "interpolate ae_vloc"
        aux2(1,1:upf(2)%mesh) = upf(2)%paw%ae_vloc(1:upf(2)%mesh)
        CALL dosplineint( upf(2)%r(1:upf(2)%mesh), aux2, upf_r(1:upf_mesh), aux1 )
        upf_paw_ae_vloc(1:upf_mesh) =   x     * upf(1)%paw%ae_vloc(1:upf_mesh) + &
                                     (1.d0-x) * aux1(1,1:upf_mesh)
        WRITE (*,*) " done"
     ELSE
        upf_paw_ae_vloc(1:upf_mesh) =   x     * upf(1)%paw%ae_vloc(1:upf_mesh) + &
                                     (1.d0-x) * upf(2)%paw%ae_vloc(1:upf_mesh)
     ENDIF
     !
  ENDIF

  IF (upf_vca%has_gipaw) THEN
     !
     ! pp_gipaw_core_orbitals
     upf_gipaw_ncore_orbitals = upf(1)%gipaw_ncore_orbitals+upf(2)%gipaw_ncore_orbitals
     !
     ! pp_gipaw_core_orbitals_el
     ALLOCATE( upf_gipaw_core_orbital_el(upf_gipaw_ncore_orbitals) )
     upf_gipaw_core_orbital_el(1:upf(1)%gipaw_ncore_orbitals) = upf(1)%gipaw_core_orbital_el
     upf_gipaw_core_orbital_el(upf(1)%gipaw_ncore_orbitals+1:upf_vca%gipaw_ncore_orbitals) = &
       upf(2)%gipaw_core_orbital_el
     !
     ! pp_gipaw_core_orbitals_n
     ALLOCATE( upf_gipaw_core_orbital_n(upf_gipaw_ncore_orbitals) )
     upf_gipaw_core_orbital_n(1:upf(1)%gipaw_ncore_orbitals) = upf(1)%gipaw_core_orbital_n
     upf_gipaw_core_orbital_n(upf(1)%gipaw_ncore_orbitals+1:upf_vca%gipaw_ncore_orbitals) = &
       upf(2)%gipaw_core_orbital_n
     !
     ! pp_gipaw_core_orbitals_l
     ALLOCATE( upf_gipaw_core_orbital_l(upf_gipaw_ncore_orbitals) )
     upf_gipaw_core_orbital_l(1:upf(1)%gipaw_ncore_orbitals) = upf(1)%gipaw_core_orbital_l
     upf_gipaw_core_orbital_l(upf(1)%gipaw_ncore_orbitals+1:upf_vca%gipaw_ncore_orbitals) = &
       upf(2)%gipaw_core_orbital_l
     !
     ! pp_gipaw_core_orbital
     ALLOCATE ( upf_gipaw_core_orbital(upf_mesh,upf_gipaw_ncore_orbitals) )
     upf_gipaw_core_orbital(1:upf_mesh,1:upf(1)%gipaw_ncore_orbitals) = &
        x * upf(1)%gipaw_core_orbital(1:upf_mesh,1:upf(1)%gipaw_ncore_orbitals)
     DO i=1,upf(2)%gipaw_ncore_orbitals
       IF (interpolate) THEN
          IF (i==1) WRITE (*,*) "interpolate gipaw_wfs_ae"
          aux2(1,1:upf(2)%mesh) = upf(2)%gipaw_core_orbital(1:upf(2)%mesh,i)
          CALL dosplineint( upf(2)%r(1:upf(2)%mesh), aux2, upf_r(1:upf_mesh), aux1 )
          upf_gipaw_core_orbital(1:upf_mesh, upf(1)%gipaw_ncore_orbitals+i) = &
             (1.d0-x) * aux1(1,1:upf_mesh)
          IF (i==upf(2)%gipaw_ncore_orbitals) WRITE (*,*) " done"
       ELSE
          upf_gipaw_core_orbital(1:upf_mesh, upf(1)%gipaw_ncore_orbitals+i) = &
             (1.d0-x) * upf(2)%gipaw_core_orbital(1:upf_mesh,i)
       ENDIF
     ENDDO
     !
     ! pp_gipaw_vlocal_ae
     ALLOCATE ( upf_gipaw_vlocal_ae(upf_mesh) )
     IF (interpolate) THEN
        WRITE (*,*) "interpolate gipaw_vlocal_ae"
        aux2(1,1:upf(2)%mesh) =  upf(2)%gipaw_vlocal_ae(1:upf(2)%mesh)

        CALL dosplineint( upf(2)%r(1:upf(2)%mesh), aux2, upf_r(1:upf_mesh), aux1 )

        ! upf(2)%gipaw_vlocal_ae(1:upf_mesh) = aux1(1,1:upf_mesh)

        ! Jivtesh - if the mesh of the first atom extends to a larger radius
        ! than the mesh of the second atom, then, for those radii that are
        ! greater than the maximum radius of the second atom, the local potential
        ! of the second atom is calculated using the expression
        ! v_local = (-2)*Z/r instead of using the extrapolated value.
        ! This is because, typically extrapolation leads to positive potentials.
        ! This is implemented in lines 240-242

        DO i=1,upf(1)%mesh
           IF ( upf(1)%r(i) > upf(2)%r(upf(2)%mesh) ) &
              ! upf(2)%gipaw_vlocal_ae(i) = -(2.0*upf(2)%zp)/upf(1)%r(i)
              aux1(1,i) = -(2.0*upf(2)%zp)/upf(1)%r(i)
        ENDDO
        upf_gipaw_vlocal_ae(1:upf_mesh) =    x     * upf(1)%gipaw_vlocal_ae(1:upf_mesh) + &
                                          (1.d0-x) * aux1(1,1:upf_mesh)
        WRITE (*,*) " done"
     ELSE
        upf_gipaw_vlocal_ae(1:upf_mesh) =    x     * upf(1)%gipaw_vlocal_ae(1:upf_mesh) + &
                                          (1.d0-x) * upf(2)%gipaw_vlocal_ae(1:upf_mesh)
     ENDIF
     !
     ! pp_gipaw_vlocal_ps
     ALLOCATE ( upf_gipaw_vlocal_ps(upf_mesh) )
     IF (interpolate) THEN
        WRITE (*,*) "interpolate gipaw_vlocal_ps"
        aux2(1,1:upf(2)%mesh) =  upf(2)%gipaw_vlocal_ps(1:upf(2)%mesh)

        CALL dosplineint( upf(2)%r(1:upf(2)%mesh), aux2, upf_r(1:upf_mesh), aux1 )

        ! upf(2)%gipaw_vlocal_ps(1:upf_mesh) = aux1(1,1:upf_mesh)

        ! Jivtesh - if the mesh of the first atom extends to a larger radius
        ! than the mesh of the second atom, then, for those radii that are
        ! greater than the maximum radius of the second atom, the local potential
        ! of the second atom is calculated using the expression
        ! v_local = (-2)*Z/r instead of using the extrapolated value.
        ! This is because, typically extrapolation leads to positive potentials.
        ! This is implemented in lines 240-242

        DO i=1,upf(1)%mesh
           IF ( upf(1)%r(i) > upf(2)%r(upf(2)%mesh) ) &
              ! upf(2)%gipaw_vlocal_ps(i) = -(2.0*upf(2)%zp)/upf(1)%r(i)
              aux1(1,i) = -(2.0*upf(2)%zp)/upf(1)%r(i)
        ENDDO
        upf_gipaw_vlocal_ps(1:upf_mesh) =    x     * upf(1)%gipaw_vlocal_ps(1:upf_mesh) + &
                                          (1.d0-x) * aux1(1,1:upf_mesh)
        WRITE (*,*) " done"
     ELSE
        upf_gipaw_vlocal_ps(1:upf_mesh) =    x     * upf(1)%gipaw_vlocal_ps(1:upf_mesh) + &
                                          (1.d0-x) * upf(2)%gipaw_vlocal_ps(1:upf_mesh)
     ENDIF
     !
     !pp_gipaw_orbitals
     ! pp_gipaw_wfs_nchannels
     IF (upf(1)%paw_as_gipaw.or.upf(2)%paw_as_gipaw) THEN
        upf_gipaw_wfs_nchannels = upf_nbeta
     ELSEIF (upf(1)%tcoulombp) THEN
        upf_gipaw_wfs_nchannels = upf(1)%gipaw_wfs_nchannels + upf(2)%gipaw_wfs_nchannels
     ELSE
        upf_gipaw_wfs_nchannels = upf(1)%gipaw_wfs_nchannels + upf(2)%gipaw_wfs_nchannels
     ENDIF
     !
     !pp_gipaw_wfs_el
     ALLOCATE (upf_gipaw_wfs_el(upf_gipaw_wfs_nchannels))
     upf_gipaw_wfs_el(1:upf(1)%gipaw_wfs_nchannels) = upf(1)%gipaw_wfs_el
     upf_gipaw_wfs_el(upf(1)%gipaw_wfs_nchannels+1:upf_nbeta) = upf(2)%gipaw_wfs_el
     !
     !pp_gipaw_wfs_n
     !ALLOCATE (upf_gipaw_wfs_n(upf_gipaw_wfs_nchannels))
     !upf_gipaw_wfs_n(1:upf(1)%gipaw_wfs_nchannels) = upf(1)%gipaw_wfs_n
     !upf_gipaw_wfs_n(upf(1)%gipaw_wfs_nchannels+1:upf_nbeta) = upf(2)%gipaw_wfs_n
     !
     !pp_gipaw_wfs_l
     !ALLOCATE (upf_gipaw_wfs_l(upf_gipaw_wfs_nchannels))
     !upf_gipaw_wfs_l(1:upf(1)%gipaw_wfs_nchannels) = upf(1)%gipaw_wfs_l
     !upf_gipaw_wfs_l(upf(1)%gipaw_wfs_nchannels+1:upf_nbeta) = upf(2)%gipaw_wfs_l
     !
     !pp_gipaw_wfs_ll
     ALLOCATE (upf_gipaw_wfs_ll(upf_gipaw_wfs_nchannels))
     upf_gipaw_wfs_ll(1:upf(1)%gipaw_wfs_nchannels) = upf(1)%gipaw_wfs_ll
     upf_gipaw_wfs_ll(upf(1)%gipaw_wfs_nchannels+1:upf_nbeta) = upf(2)%gipaw_wfs_ll
     !
     !pp_gipaw_wfs_rcut
     ALLOCATE (upf_gipaw_wfs_rcut(upf_gipaw_wfs_nchannels))
     upf_gipaw_wfs_rcut(1:upf(1)%gipaw_wfs_nchannels) = upf(1)%gipaw_wfs_rcut
     upf_gipaw_wfs_rcut(upf(1)%gipaw_wfs_nchannels+1:upf_nbeta) = upf(2)%gipaw_wfs_rcut
     !
     !pp_gipaw_wfs_rcutus
     ALLOCATE (upf_gipaw_wfs_rcutus(upf_gipaw_wfs_nchannels))
     upf_gipaw_wfs_rcutus(1:upf(1)%gipaw_wfs_nchannels) = upf(1)%gipaw_wfs_rcutus
     upf_gipaw_wfs_rcutus(upf(1)%gipaw_wfs_nchannels+1:upf_nbeta) = upf(2)%gipaw_wfs_rcutus
     !
     !pp_gipaw_wfs_ae
     ALLOCATE ( upf_gipaw_wfs_ae(upf_mesh,upf_gipaw_wfs_nchannels) )
     upf_gipaw_wfs_ae(1:upf_mesh,1:upf(1)%gipaw_wfs_nchannels) = &
        x * upf(1)%gipaw_wfs_ae(1:upf_mesh,1:upf(1)%gipaw_wfs_nchannels)
     DO i=1,upf(2)%gipaw_wfs_nchannels
       IF (interpolate) THEN
          IF (i==1) WRITE (*,*) "interpolate gipaw_wfs_ae"
          aux2(1,1:upf(2)%mesh) = upf(2)%gipaw_wfs_ae(1:upf(2)%mesh,i)
          CALL dosplineint( upf(2)%r(1:upf(2)%mesh), aux2, upf_r(1:upf_mesh), aux1 )
          upf_gipaw_wfs_ae(1:upf_mesh, upf(1)%gipaw_wfs_nchannels+i) = &
             (1.d0-x) * aux1(1,1:upf_mesh)
          IF (i==upf(2)%gipaw_wfs_nchannels) WRITE (*,*) " done"
       ELSE
          upf_gipaw_wfs_ae(1:upf_mesh, upf(1)%gipaw_wfs_nchannels+i) = &
             (1.d0-x) * upf(2)%gipaw_wfs_ae(1:upf_mesh,i)
       ENDIF
     ENDDO
     !
     !pp_gipaw_wfs_ps
     ALLOCATE ( upf_gipaw_wfs_ps(upf_mesh,upf_gipaw_wfs_nchannels) )
     upf_gipaw_wfs_ps(1:upf_mesh,1:upf(1)%gipaw_wfs_nchannels) = &
        x * upf(1)%gipaw_wfs_ps(1:upf_mesh,1:upf(1)%gipaw_wfs_nchannels)
     DO i=1,upf(2)%gipaw_wfs_nchannels
       IF (interpolate) THEN
          IF (i==1) WRITE (*,*) "interpolate gipaw_wfs_ps"
          aux2(1,1:upf(2)%mesh) = upf(2)%gipaw_wfs_ps(1:upf(2)%mesh,i)
          CALL dosplineint( upf(2)%r(1:upf(2)%mesh), aux2, upf_r(1:upf_mesh), aux1 )
          upf_gipaw_wfs_ps(1:upf_mesh, upf(1)%gipaw_wfs_nchannels+i) = &
             (1.d0-x) * aux1(1,1:upf_mesh)
          IF (i==upf(2)%gipaw_wfs_nchannels) WRITE (*,*) " done"
       ELSE
          upf_gipaw_wfs_ps(1:upf_mesh, upf(1)%gipaw_wfs_nchannels+i) = &
             (1.d0-x) * upf(2)%gipaw_wfs_ps(1:upf_mesh,i)
       ENDIF
     ENDDO
     !
  ENDIF
  WRITE (*,*) " "

  ! pp_info
  ! pp_header
  upf_vca%lmax = upf_lmax
  upf_vca%nlcc = upf_nlcc
  upf_vca%zp = upf_zp
  ! pp_nlcc
  upf_vca%rho_atc = upf_rho_atc
  ! pp_local
  upf_vca%vloc = upf_vloc0

  ! pp_pswfc
  NULLIFY( upf_vca%chi )
  ALLOCATE( upf_vca%chi(upf_mesh, upf_ntwfc) )
  upf_vca%chi = upf_chi
  ! pp_rhoatom
  upf_vca%rho_at = upf_rho_at

  ! pp_nonlocal
  ! pp_beta
  upf_vca%nbeta = upf_nbeta

  NULLIFY( upf_vca%kbeta,          &
           upf_vca%lll,            &
           upf_vca%beta,           &
           upf_vca%els_beta,       &
           upf_vca%dion,           &
           upf_vca%qqq,            &
           upf_vca%rcut,           &
           upf_vca%rcutus          &
  )

  ALLOCATE( upf_vca%kbeta(upf_nbeta),          &
            upf_vca%lll(upf_nbeta),            &
            upf_vca%beta(upf_mesh, upf_nbeta), &
            upf_vca%els_beta(upf_nbeta),       &
            upf_vca%dion(upf_nbeta, upf_nbeta),&
            upf_vca%qqq(upf_nbeta, upf_nbeta), &
            upf_vca%rcut(upf_nbeta),           &
            upf_vca%rcutus(upf_nbeta)          &
  )

  upf_vca%kbeta = upf_ikk2
  upf_vca%beta = upf_betar
  upf_vca%lll = upf_lll
  ! pp_dij
  upf_vca%dion = upf_dion

  IF (matches(upf_vca%typ, "USPP") .or. matches(upf_vca%typ, "PAW")) THEN
    ! pp_qij
    upf_vca%qqq = upf_qqq
    IF( upf_vca%q_with_l ) THEN
       NULLIFY( upf_vca%qfuncl )
       ALLOCATE( upf_vca%qfuncl(upf_mesh, upf_nbeta*(upf_nbeta+1)/2, 0:2*upf_lmax) )
       upf_vca%qfuncl = upf_qfuncl
    ELSE
       NULLIFY( upf_vca%qfunc )
       ALLOCATE( upf_vca%qfunc(upf_mesh, upf_nbeta*(upf_nbeta+1)/2) )
       upf_vca%qfunc = upf_qfunc
    ENDIF
    ! pp_qfcoef
    IF ( allocated(upf_qfcoef) ) THEN
       NULLIFY( upf_vca%qfcoef )
       ALLOCATE( upf_vca%qfcoef(upf_nqf, upf_nqlc, upf_nbeta, upf_nbeta) )
       upf_vca%qfcoef = upf_qfcoef
    ENDIF
  ENDIF
  
  upf_vca%qqq_eps = min(upf(1)%qqq_eps,upf(2)%qqq_eps) 
  ! qqq_eps ! qfunc is null if its norm is .lt. qqq_eps

  upf_vca%kkbeta = maxval(upf_vca%kbeta)
  ! For PAW augmentation charge may extend a bit further:
  !IF(upf(1)%tpawp) upf_vca%kkbeta = MAX(upf(1)%kkbeta, upf(1)%paw%iraug)
  !IF(upf(2)%tpawp) upf_vca%kkbeta = MAX(upf(2)%kkbeta, upf(2)%paw%iraug)
  IF(upf(1)%tpawp.or.upf(1)%tpawp) upf_vca%kkbeta = max(upf(1)%kkbeta,upf(2)%kkbeta)

  upf_vca%els_beta(1:upf(1)%nbeta) = upf(1)%els_beta
  upf_vca%els_beta(upf(1)%nbeta+1:upf_nbeta) = upf(2)%els_beta

  upf_vca%rcut(1:upf(1)%nbeta) = upf(1)%rcut
  upf_vca%rcut(upf(1)%nbeta+1:upf_nbeta) = upf(2)%rcut

  upf_vca%rcutus(1:upf(1)%nbeta) = upf(1)%rcutus
  upf_vca%rcutus(upf(1)%nbeta+1:upf_nbeta) = upf(2)%rcutus

  IF (upf_vca%has_so) THEN
     NULLIFY( upf_vca%jjj )
     ALLOCATE( upf_vca%jjj(upf_nbeta) )
     upf_vca%jjj = upf_jjj
     !
     NULLIFY( upf_vca%nn )
     ALLOCATE( upf_vca%nn(upf_ntwfc) )
     upf_vca%nn = upf_nn
     !
     NULLIFY( upf_vca%jchi )
     ALLOCATE(upf_vca%jchi(upf_ntwfc) )
     upf_vca%jchi = upf_jchi
  ENDIF

  IF (matches(upf_vca%typ, "SL")) THEN
     upf1_ind = 1
     upf2_ind = 1
     upf_ind  = 1
     IF (upf_vca%has_so) THEN
        DO i = 1,upf(1)%nbeta
           IF ( l > 0 .AND. ABS (upf(1)%jjj(i)-l-0.5d0) < 0.001d0) upf1_ind = 2
        ENDDO
        DO i = 1,upf(2)%nbeta
           IF ( l > 0 .AND. ABS (upf(2)%jjj(i)-l-0.5d0) < 0.001d0) upf2_ind = 2
        ENDDO
        upf_ind = max(upf1_ind,upf2_ind)
     ENDIF
     ALLOCATE( upf_vnl(upf_mesh,0:upf_lmax,upf_ind) )
     upf_vnl(:,:,:) = 0.d0
     DO i = 0,upf_lmax
        DO j = 1,upf1_ind
           upf_vnl(1:upf_mesh,i,j) = x * upf(1)%vnl(1:upf_mesh,i,j)
        ENDDO
     ENDDO
     DO i = 0,upf_lmax
        DO j = 1,upf2_ind
              IF (interpolate) THEN
              WRITE (*,*) " interpolate gipaw_wfs_ps"
              aux2(1,1:upf(2)%mesh) = upf(2)%vnl(1:upf(2)%mesh,i,j)
              CALL dosplineint( upf(2)%r(1:upf(2)%mesh), aux2, upf_r(1:upf_mesh), aux1 )
              ! upf_vnl(1:upf_mesh,i,2) = aux1(1,1:upf_mesh)
              upf_vnl(1:upf_mesh,i,j)  = upf_vnl(1:upf_mesh,i,j) + &
                                  (1.d0-x) * aux1(1,1:upf_mesh)
           ELSE
              upf_vnl(1:upf_mesh,i,j) = upf_vnl(1:upf_mesh,i,j) + &
                                 (1.d0-x) * upf(2)%vnl(1:upf_mesh,i,j)
           ENDIF
        ENDDO
     ENDDO
     NULLIFY( upf_vca%vnl )
     ALLOCATE( upf_vca%vnl(upf_mesh,0:upf_lmax,upf_ind) )
     upf_vca%vnl = upf_vnl
     ! (semilocal form) only for single-channel NC PP
     ! Wavefunctions and projectors
  ENDIF

  IF (matches(upf_vca%typ, "PAW")) THEN
     NULLIFY( upf_vca%paw%augmom )
     NULLIFY( upf_vca%paw%ae_rho_atc )
     NULLIFY( upf_vca%aewfc )
     NULLIFY( upf_vca%pswfc )
     NULLIFY( upf_vca%paw%pfunc )
     NULLIFY( upf_vca%paw%ptfunc )
     NULLIFY( upf_vca%paw%ae_vloc )
     NULLIFY( upf_vca%paw%oc )
     !
     ALLOCATE( upf_vca%paw%augmom(upf_nbeta, upf_nbeta, 0:2*upf_lmax))
     ALLOCATE( upf_vca%paw%ae_rho_atc(upf_mesh) )
     ALLOCATE( upf_vca%aewfc(upf_mesh,upf_nbeta) )
     ALLOCATE( upf_vca%pswfc(upf_mesh,upf_nbeta) )
     ALLOCATE( upf_vca%paw%pfunc(upf_mesh, upf_nbeta, upf_nbeta) )
     ALLOCATE( upf_vca%paw%ptfunc(upf_mesh, upf_nbeta, upf_nbeta) )
     ALLOCATE( upf_vca%paw%ae_vloc(upf_mesh) )
     ALLOCATE( upf_vca%paw%oc(upf_nbeta) )
     !
     IF (matches(upf_vca%typ, "PAW")) THEN
        upf_vca%tvanp = .false.
     ELSE
        upf_vca%tvanp = upf(1)%tvanp.or.upf(2)%tvanp
     ENDIF
     !upf_vca%tpawp = upf(1)%tpawp.or.upf(2)%tpawp
     !upf_vca%tcoulombp = upf(1)%tcoulombp.or.upf(2)%tcoulombp
     !upf_vca%nlcc    = upf(1)%nlcc.or.upf(2)%nlcc
     !
     upf_vca%paw%augshape = upf_paw_augshape ! shape
     upf_vca%paw%raug     = upf_paw_raug     ! cutoff_r
     upf_vca%paw%iraug    = upf_paw_iraug    ! cutoff_r_index, near upf%mesh
     upf_vca%paw%lmax_aug = upf_paw_lmax_aug ! near upf%lmax_rho = 2*upf%lmax
     upf_vca%qqq_eps      = upf_qqq_eps      ! augmentation_epsilon
     !
     upf_vca%paw%core_energy = x * upf(1)%paw%core_energy + (1.d0-x) * upf(2)%paw%core_energy
     !
     upf_vca%paw%augmom = upf_paw_augmom
     upf_vca%paw%ae_rho_atc = upf_paw_ae_rho_atc
     upf_vca%aewfc = upf_aewfc
     upf_vca%pswfc = upf_pswfc
     upf_vca%paw%pfunc = upf_paw_pfunc
     upf_vca%paw%ptfunc = upf_paw_ptfunc
     upf_vca%paw%ae_vloc = upf_paw_ae_vloc
     !
     ! occ, f6,2
     upf_vca%paw%oc(1:upf(1)%nbeta) = x * upf(1)%paw%oc
     upf_vca%paw%oc(upf(1)%nbeta+1:upf_nbeta) = (1.d0-x) * upf(2)%paw%oc
     !
     IF (upf_vca%has_so) THEN
        NULLIFY( upf_vca%paw%pfunc_rel )
        ALLOCATE( upf_vca%paw%pfunc_rel(upf_mesh, upf_nbeta, upf_nbeta) )
        upf_vca%paw%pfunc_rel = upf_paw_pfunc_rel
        !
        NULLIFY( upf_vca%paw%aewfc_rel )
        ALLOCATE( upf_vca%paw%aewfc_rel(upf_mesh,upf_nbeta) )
        upf_vca%paw%aewfc_rel = upf_paw_aewfc_rel
     ELSE
        NULLIFY( upf_vca%paw%pfunc_rel )
        ALLOCATE( upf_vca%paw%pfunc_rel(upf_mesh, upf_nbeta, upf_nbeta) )
        upf_vca%paw%pfunc_rel(:,:,:) = 0.0d0
        !
        NULLIFY( upf_vca%paw%aewfc_rel )
        ALLOCATE( upf_vca%paw%aewfc_rel(upf_mesh,upf_nbeta) )
        upf_vca%paw%aewfc_rel(:,:) = 0.0d0
     ENDIF
  ENDIF

  IF (upf_vca%has_gipaw) THEN
     NULLIFY( upf_vca%els, upf_vca%nchi, upf_vca%lchi, upf_vca%oc, &
              upf_vca%epseu, upf_vca%rcut_chi, upf_vca%rcutus_chi)
     NULLIFY( upf_vca%rinner )
     !
     NULLIFY ( upf_vca%gipaw_core_orbital_el )
     NULLIFY ( upf_vca%gipaw_core_orbital_n )
     NULLIFY ( upf_vca%gipaw_core_orbital_l )
     NULLIFY ( upf_vca%gipaw_core_orbital )
     NULLIFY ( upf_vca%gipaw_vlocal_ae )
     NULLIFY ( upf_vca%gipaw_vlocal_ps )
     !
     NULLIFY ( upf_vca%gipaw_wfs_el )
     NULLIFY ( upf_vca%gipaw_wfs_ll )
     NULLIFY ( upf_vca%gipaw_wfs_rcut )
     NULLIFY ( upf_vca%gipaw_wfs_rcutus )
     NULLIFY ( upf_vca%gipaw_wfs_ae )
     NULLIFY ( upf_vca%gipaw_wfs_ps )
     !
     ALLOCATE( upf_vca%els(upf_ntwfc), &
               upf_vca%nchi(upf_ntwfc), &
               upf_vca%lchi(upf_ntwfc), &
               upf_vca%oc(upf_ntwfc), &
               upf_vca%rcut_chi(upf_ntwfc), &
               upf_vca%rcutus_chi(upf_ntwfc), &
               upf_vca%epseu(upf_ntwfc) &
             )
     ALLOCATE( upf_vca%rinner(0:2*upf_lmax) )
     !
     ALLOCATE( upf_vca%gipaw_core_orbital_el(upf_gipaw_ncore_orbitals) )
     ALLOCATE( upf_vca%gipaw_core_orbital_n(upf_gipaw_ncore_orbitals) )
     ALLOCATE( upf_vca%gipaw_core_orbital_l(upf_gipaw_ncore_orbitals) )
     ALLOCATE( upf_vca%gipaw_core_orbital(upf_mesh,upf_gipaw_ncore_orbitals) )
     ALLOCATE( upf_vca%gipaw_vlocal_ae(upf_mesh) )
     ALLOCATE( upf_vca%gipaw_vlocal_ps(upf_mesh) )
     !
     ALLOCATE( upf_vca%gipaw_wfs_el(upf_gipaw_wfs_nchannels) )
     ALLOCATE( upf_vca%gipaw_wfs_ll(upf_gipaw_wfs_nchannels) )
     ALLOCATE( upf_vca%gipaw_wfs_rcut(upf_gipaw_wfs_nchannels) )
     ALLOCATE( upf_vca%gipaw_wfs_rcutus(upf_gipaw_wfs_nchannels) )
     ALLOCATE( upf_vca%gipaw_wfs_ae(upf_mesh,upf_gipaw_wfs_nchannels) )
     ALLOCATE( upf_vca%gipaw_wfs_ps(upf_mesh,upf_gipaw_wfs_nchannels) )
     !
     IF (matches(upf_vca%typ, "PAW")) THEN
        upf_vca%tvanp = .false.
     ELSE
        upf_vca%tvanp = upf(1)%tvanp.or.upf(2)%tvanp
     ENDIF
     !upf_vca%tpawp = upf(1)%tpawp.or.upf(2)%tpawp
     !upf_vca%tcoulombp = upf(1)%tcoulombp.or.upf(2)%tcoulombp
     !upf_vca%nlcc    = upf(1)%nlcc.or.upf(2)%nlcc
     !upf_vca%nv = upf_nv ! UPF file three-digit version i.e. 2.0.0, CHARACTER
     upf_vca%lmax    = upf_lmax
     upf_vca%lmax_rho = upf_lmax_rho
     upf_vca%nwfc    = upf_ntwfc
     !
     upf_vca%rinner  = upf_rinner
     !
     upf_vca%lloc    = max(upf(1)%lloc,upf(2)%lloc)   ! l_local
     upf_vca%rcloc   = max(upf(1)%rcloc,upf(2)%rcloc)
     !
     ! pp mesh
     upf_vca%dx      = upf_dx
     upf_vca%mesh    = upf_mesh
     upf_vca%xmin    = upf_xmin
     upf_vca%rmax    = upf_rmax
     upf_vca%zmesh   = upf_zmesh
     !
     upf_vca%tpawp   = upf(1)%tpawp.or.upf(2)%tpawp
     upf_vca%has_wfc = upf(1)%has_wfc.or.upf(2)%has_wfc
     upf_vca%has_gipaw = upf(1)%has_gipaw.or.upf(2)%has_gipaw
     upf_vca%paw_as_gipaw = upf(1)%paw_as_gipaw.or.upf(2)%paw_as_gipaw
     !
     ! nl, a2
     upf_vca%els(1:upf(1)%nwfc) = upf(1)%els
     upf_vca%els(upf(1)%nwfc+1:upf_ntwfc) = upf(2)%els
     !
     ! pn, i3
     upf_vca%nchi(1:upf(1)%nwfc) = upf(1)%nchi
     upf_vca%nchi(upf(1)%nwfc+1:upf_ntwfc) = upf(2)%nchi
     !
     ! l, i3
     upf_vca%lchi(1:upf(1)%nwfc) = upf(1)%lchi
     upf_vca%lchi(upf(1)%nwfc+1:upf_ntwfc) = upf(2)%lchi
     !
     ! occ, f6,2
     upf_vca%oc(1:upf(1)%nwfc) = x * upf(1)%oc
     upf_vca%oc(upf(1)%nwfc+1:upf_ntwfc) = (1.d0-x) * upf(2)%oc
     !
     ! Rcut, f11.3
     upf_vca%rcut_chi(1:upf(1)%nwfc) = upf(1)%rcut_chi
     upf_vca%rcut_chi(upf(1)%nwfc+1:upf_ntwfc) = upf(2)%rcut_chi
     !
     ! Rcut US, f11.3
     upf_vca%rcutus_chi(1:upf(1)%nwfc) = upf(1)%rcutus_chi
     upf_vca%rcutus_chi(upf(1)%nwfc+1:upf_ntwfc) = upf(2)%rcutus_chi
     !
     ! E pseu, f11.3
     upf_vca%epseu(1:upf(1)%nwfc) = upf(1)%epseu
     upf_vca%epseu(upf(1)%nwfc+1:upf_ntwfc) = upf(2)%epseu
     !
     !upf_vca%gipaw_data_format = upf_gipaw_data_format    ! The version of the format
     !
     upf_vca%gipaw_ncore_orbitals  = upf_gipaw_ncore_orbitals
     upf_vca%gipaw_core_orbital_el = upf_gipaw_core_orbital_el
     upf_vca%gipaw_core_orbital_n  = upf_gipaw_core_orbital_n
     upf_vca%gipaw_core_orbital_l  = upf_gipaw_core_orbital_l
     upf_vca%gipaw_core_orbital    = upf_gipaw_core_orbital
     upf_vca%gipaw_vlocal_ae       = upf_gipaw_vlocal_ae
     upf_vca%gipaw_vlocal_ps       = upf_gipaw_vlocal_ps
     !
     upf_vca%gipaw_wfs_nchannels = upf_gipaw_wfs_nchannels
     upf_vca%gipaw_wfs_el     = upf_gipaw_wfs_el
     upf_vca%gipaw_wfs_ll     = upf_gipaw_wfs_ll
     upf_vca%gipaw_wfs_rcut   = upf_gipaw_wfs_rcut
     upf_vca%gipaw_wfs_rcutus = upf_gipaw_wfs_rcutus
     upf_vca%gipaw_wfs_ae     = upf_gipaw_wfs_ae
     upf_vca%gipaw_wfs_ps     = upf_gipaw_wfs_ps
  ENDIF
  
  !! DEBUG
  debag_flag = 1
  IF (debag_flag==1) THEN
     WRITE (*,*) "****************************************"
     !!! upf(1)
     WRITE (*,*) "upf(1)%nbeta = ", upf(1)%nbeta
     WRITE (*,*) "shape of upf(1)%kbeta = ", shape(upf(1)%kbeta)
     WRITE (*,*) "upf(1)%kbeta = ", upf(1)%kbeta
     WRITE (*,*) "upf(1)%kkbeta = ", upf(1)%kkbeta
     WRITE (*,*) "shape of upf(1)%lll = ", shape(upf(1)%lll)
     WRITE (*,*) "shape of upf(1)%beta = ", shape(upf(1)%beta)
     WRITE (*,*) "shape of upf(1)%els_beta = ", shape(upf(1)%els_beta)
     WRITE (*,*) "shape of upf(1)%dion = ", shape(upf(1)%dion)
     WRITE (*,*) "shape of upf(1)%qqq = ", shape(upf(1)%qqq)
     WRITE (*,*) "shape of upf(1)%qfuncl = ", shape(upf(1)%qfuncl)
     WRITE (*,*) "shape of upf(1)%qfcoef = ", shape(upf(1)%qfcoef)
     WRITE (*,*) "shape of upf(1)%aewfc = ", shape(upf(1)%aewfc)
     WRITE (*,*) "shape of upf(1)%pswfc = ", shape(upf(1)%pswfc)
     WRITE (*,*) "shape of upf(1)%rcut = ", shape(upf(1)%rcut)
     WRITE (*,*) "shape of upf(1)%rcutus = ", shape(upf(1)%rcutus)
     WRITE (*,*) " "
     WRITE (*,*) "is_ultrasoft,upf(1)%tvanp = ", upf(1)%tvanp
     WRITE (*,*) "is_paw,upf(1)%tpawp = ", upf(1)%tpawp
     WRITE (*,*) "is_coulomb,upf(1)%tcoulombp = ", upf(1)%tcoulombp
     WRITE (*,*) " "
     !
     IF (matches(upf(1)%typ, "PAW")) THEN
        WRITE (*,*) "PAW"
        WRITE (*,*) "upf(1)%nbeta = ", upf(1)%nbeta
        WRITE (*,*) "upf(1)%lmax = ", upf(1)%lmax
        WRITE (*,*) " (2*lmax+1) = ", 2*upf(1)%lmax+1
        WRITE (*,*) "upf(1)%mesh = ", upf(1)%mesh
        WRITE (*,*) "shape of upf(1)%paw%augmom = ", shape(upf(1)%paw%augmom)
       !WRITE (*,*) "shape of upf(1)%paw%augfun = ", shape(upf(1)%paw%augfun)
        WRITE (*,*) "shape of upf(1)%paw%ae_rho_atc = ", shape(upf(1)%paw%ae_rho_atc)
        WRITE (*,*) "shape of upf(1)%aewfc = ", shape(upf(1)%aewfc)
        WRITE (*,*) "shape of upf(1)%paw%aewfc_rel = ", shape(upf(1)%paw%aewfc_rel)
        WRITE (*,*) "shape of upf(1)%pswfc = ", shape(upf(1)%pswfc)
        WRITE (*,*) "shape of upf(1)%paw%pfunc = ", shape(upf(1)%paw%pfunc)
        WRITE (*,*) "shape of upf(1)%paw%pfunc_rel = ", shape(upf(1)%paw%pfunc_rel)
        WRITE (*,*) "shape of upf(1)%paw%ptfunc = ", shape(upf(1)%paw%ptfunc)
        WRITE (*,*) "shape of upf(1)%paw%ae_vloc = ", shape(upf(1)%paw%ae_vloc)
        WRITE (*,*) "shape of upf(1)%paw%oc = ", shape(upf(1)%paw%oc)
        !WRITE (*,*) "shape of upf(1)%paw%kdiff = ", shape(upf(1)%paw%kdiff)
        !WRITE (*,*) "shape of upf(1)%paw%sqr = ", shape(upf(1)%paw%sqr)
        WRITE (*,*) "upf(1)%paw%augshape = ", upf(1)%paw%augshape
        WRITE (*,*) "upf(1)%paw%raug = ", upf(1)%paw%raug
        WRITE (*,*) "upf(1)%paw%iraug = ", upf(1)%paw%iraug
        WRITE (*,*) "upf(1)%paw%lmax_aug = ", upf(1)%paw%lmax_aug
        WRITE (*,*) "upf(1)%qqq_eps  = ", upf(1)%qqq_eps 
        WRITE (*,*) "upf(1)%paw%core_energy  = ", upf(1)%paw%core_energy
        WRITE (*,*) " "
     ENDIF
     !
     IF (upf(2)%has_gipaw) THEN
        WRITE (*,*) "GIPAW"
        WRITE (*,*) "upf(1)%nwfc = ", upf(1)%nwfc
        WRITE (*,*) "shape of upf(1)%els = ", shape(upf(1)%els)
        WRITE (*,*) "shape of upf(1)%nchi = ", shape(upf(1)%nchi)
        WRITE (*,*) "shape of upf(1)%lchi = ", shape(upf(1)%lchi)
        WRITE (*,*) "shape of upf(1)%oc = ", shape(upf(1)%oc)
        WRITE (*,*) "shape of upf(1)%rcut_chi = ", shape(upf(1)%rcut_chi)
        WRITE (*,*) "shape of upf(1)%rcutus_chi = ", shape(upf(1)%rcutus_chi)
        WRITE (*,*) "shape of upf(1)%epseu = ", shape(upf(1)%epseu)
        WRITE (*,*) "upf(1)%lmax = ", upf(1)%lmax
        WRITE (*,*) " (2*lmax+1) = ", 2*upf(1)%lmax+1
        WRITE (*,*) "shape of upf(1)%rinner = ", shape(upf(1)%rinner)
        !
        WRITE (*,*) "upf(1)%gipaw_ncore_orbitals = ", upf(1)%gipaw_ncore_orbitals
        WRITE (*,*) "upf(1)%mesh = ", upf(1)%mesh
        WRITE (*,*) "shape of upf(1)%gipaw_core_orbital_el = ", shape(upf(1)%gipaw_core_orbital_el)
        WRITE (*,*) "shape of upf(1)%gipaw_core_orbital_n = ", shape(upf(1)%gipaw_core_orbital_n)
        WRITE (*,*) "shape of upf(1)%gipaw_core_orbital_l = ", shape(upf(1)%gipaw_core_orbital_l)
        WRITE (*,*) "shape of upf(1)%gipaw_core_orbital = ", shape(upf(1)%gipaw_core_orbital)
        WRITE (*,*) "shape of upf(1)%gipaw_vlocal_ae = ", shape(upf(1)%gipaw_vlocal_ae)
        WRITE (*,*) "shape of upf(1)%gipaw_vlocal_ps = ", shape(upf(1)%gipaw_vlocal_ps)
        !
        WRITE (*,*) "upf(1)%gipaw_wfs_nchannels = ", upf(1)%gipaw_wfs_nchannels
        WRITE (*,*) "upf(1)%mesh = ", upf(1)%mesh
        WRITE (*,*) "shape of upf(1)%gipaw_wfs_el = ", shape(upf(1)%gipaw_wfs_el)
        !WRITE (*,*) "shape of upf(1)%gipaw_wfs_n = ", shape(upf(1)%gipaw_wfs_n)
        !WRITE (*,*) "shape of upf(1)%gipaw_wfs_l = ", shape(upf(1)%gipaw_wfs_l)
        WRITE (*,*) "shape of upf(1)%gipaw_wfs_ll = ", shape(upf(1)%gipaw_wfs_ll)
        WRITE (*,*) "shape of upf(1)%gipaw_wfs_rcut = ", shape(upf(1)%gipaw_wfs_rcut)
        WRITE (*,*) "shape of upf(1)%gipaw_wfs_rcutus = ", shape(upf(1)%gipaw_wfs_rcutus)
        WRITE (*,*) "shape of upf(1)%gipaw_wfs_ae = ", shape(upf(1)%gipaw_wfs_ae)
        WRITE (*,*) "shape of upf(1)%gipaw_wfs_ps = ", shape(upf(1)%gipaw_wfs_ps)
     ENDIF
     !
     WRITE (*,*) "****************************************"
     !!! upf(2)
     WRITE (*,*) "upf(2)%nbeta = ", upf(2)%nbeta
     WRITE (*,*) "shape of upf(2)%kbeta = ", shape(upf(2)%kbeta)
     WRITE (*,*) "upf(2)%kbeta = ", upf(2)%kbeta
     WRITE (*,*) "upf(2)%kkbeta = ", upf(2)%kkbeta
     WRITE (*,*) "shape of upf(2)%lll = ", shape(upf(2)%lll)
     WRITE (*,*) "shape of upf(2)%beta = ", shape(upf(2)%beta)
     WRITE (*,*) "shape of upf(2)%els_beta = ", shape(upf(2)%els_beta)
     WRITE (*,*) "shape of upf(2)%dion = ", shape(upf(2)%dion)
     WRITE (*,*) "shape of upf(2)%qqq = ", shape(upf(2)%qqq)
     WRITE (*,*) "shape of upf(2)%qfuncl = ", shape(upf(2)%qfuncl)
     WRITE (*,*) "shape of upf(2)%qfcoef = ", shape(upf(2)%qfcoef)
     WRITE (*,*) "shape of upf(2)%aewfc = ", shape(upf(2)%aewfc)
     WRITE (*,*) "shape of upf(2)%pswfc = ", shape(upf(2)%pswfc)
     WRITE (*,*) "shape of upf(2)%rcut = ", shape(upf(2)%rcut)
     WRITE (*,*) "shape of upf(2)%rcutus = ", shape(upf(2)%rcutus)
     WRITE (*,*) " "
     WRITE (*,*) "is_ultrasoft,upf(2)%tvanp = ", upf(2)%tvanp
     WRITE (*,*) "is_paw,upf(2)%tpawp = ", upf(2)%tpawp
     WRITE (*,*) "is_coulomb,upf(2)%tcoulombp = ", upf(2)%tcoulombp
     WRITE (*,*) " "
     !
     IF (matches(upf(2)%typ, "PAW")) THEN
        WRITE (*,*) "PAW"
        WRITE (*,*) "upf(2)%nbeta = ", upf(2)%nbeta
        WRITE (*,*) "upf(2)%lmax = ", upf(2)%lmax
        WRITE (*,*) " (2*lmax+1) = ", 2*upf(2)%lmax+1
        WRITE (*,*) "upf(2)%mesh = ", upf(2)%mesh
        WRITE (*,*) "shape of upf(2)%paw%augmom = ", shape(upf(2)%paw%augmom)
        !WRITE (*,*) "shape of upf(2)%paw%augfun = ", shape(upf(2)%paw%augfun)
        WRITE (*,*) "shape of upf(2)%paw%ae_rho_atc = ", shape(upf(2)%paw%ae_rho_atc)
        WRITE (*,*) "shape of upf(2)%aewfc = ", shape(upf(2)%aewfc)
        WRITE (*,*) "shape of upf(2)%paw%aewfc_rel = ", shape(upf(2)%paw%aewfc_rel)
        WRITE (*,*) "shape of upf(2)%pswfc = ", shape(upf(2)%pswfc)
        WRITE (*,*) "shape of upf(2)%paw%pfunc = ", shape(upf(2)%paw%pfunc)
        WRITE (*,*) "shape of upf(2)%paw%pfunc_rel = ", shape(upf(2)%paw%pfunc_rel)
        WRITE (*,*) "shape of upf(2)%paw%ptfunc = ", shape(upf(2)%paw%ptfunc)
        WRITE (*,*) "shape of upf(2)%paw%ae_vloc = ", shape(upf(2)%paw%ae_vloc)
        WRITE (*,*) "shape of upf(2)%paw%oc = ", shape(upf(2)%paw%oc)
        !WRITE (*,*) "shape of upf(2)%kdiff = ", shape(upf(2)%kdiff)
        !WRITE (*,*) "shape of upf(2)%sqr = ", shape(upf(2)%sqr)
        WRITE (*,*) "upf(2)%paw%augshape = ", upf(2)%paw%augshape
        WRITE (*,*) "upf(2)%paw%raug = ", upf(2)%paw%raug
        WRITE (*,*) "upf(2)%paw%iraug = ", upf(2)%paw%iraug
        WRITE (*,*) "upf(2)%paw%lmax_aug = ", upf(2)%paw%lmax_aug
        WRITE (*,*) "upf(2)%qqq_eps  = ", upf(2)%qqq_eps 
        WRITE (*,*) "upf(2)%paw%core_energy  = ", upf(2)%paw%core_energy
        WRITE (*,*) " "
     ENDIF
     !
     IF (upf(2)%has_gipaw) THEN
        WRITE (*,*) "GIPAW"
        WRITE (*,*) "upf(2)%nwfc = ", upf(2)%nwfc
        WRITE (*,*) "shape of upf(2)%els = ", shape(upf(2)%els)
        WRITE (*,*) "shape of upf(2)%nchi = ", shape(upf(2)%nchi)
        WRITE (*,*) "shape of upf(2)%lchi = ", shape(upf(2)%lchi)
        WRITE (*,*) "shape of upf(2)%oc = ", shape(upf(2)%oc)
        WRITE (*,*) "shape of upf(2)%rcut_chi = ", shape(upf(2)%rcut_chi)
        WRITE (*,*) "shape of upf(2)%rcutus_chi = ", shape(upf(2)%rcutus_chi)
        WRITE (*,*) "shape of upf(2)%epseu = ", shape(upf(2)%epseu)
        WRITE (*,*) "upf(2)%lmax = ", upf(2)%lmax
        WRITE (*,*) " (2*lmax+1) = ", 2*upf(2)%lmax+1
        WRITE (*,*) "shape of upf(2)%rinner = ", shape(upf(2)%rinner)
        !
        WRITE (*,*) "upf(2)%gipaw_ncore_orbitals = ", upf(2)%gipaw_ncore_orbitals
        WRITE (*,*) "upf(2)%mesh = ", upf(2)%mesh
        WRITE (*,*) "shape of upf(2)%gipaw_core_orbital_el = ", shape(upf(2)%gipaw_core_orbital_el)
        WRITE (*,*) "shape of upf(2)%gipaw_core_orbital_n = ", shape(upf(2)%gipaw_core_orbital_n)
        WRITE (*,*) "shape of upf(2)%gipaw_core_orbital_l = ", shape(upf(2)%gipaw_core_orbital_l)
        WRITE (*,*) "shape of upf(2)%gipaw_core_orbital = ", shape(upf(2)%gipaw_core_orbital)
        WRITE (*,*) "shape of upf(2)%gipaw_vlocal_ae = ", shape(upf(2)%gipaw_vlocal_ae)
        WRITE (*,*) "shape of upf(2)%gipaw_vlocal_ps = ", shape(upf(2)%gipaw_vlocal_ps)
        !
        WRITE (*,*) "upf(2)%gipaw_wfs_nchannels = ", upf(2)%gipaw_wfs_nchannels
        WRITE (*,*) "upf(2)%mesh = ", upf(2)%mesh
        WRITE (*,*) "shape of upf(2)%gipaw_wfs_el = ", shape(upf(2)%gipaw_wfs_el)
        !WRITE (*,*) "shape of upf(2)%gipaw_wfs_n = ", shape(upf(2)%gipaw_wfs_n)
        !WRITE (*,*) "shape of upf(2)%gipaw_wfs_l = ", shape(upf(2)%gipaw_wfs_l)
        WRITE (*,*) "shape of upf(2)%gipaw_wfs_ll = ", shape(upf(2)%gipaw_wfs_ll)
        WRITE (*,*) "shape of upf(2)%gipaw_wfs_rcut = ", shape(upf(2)%gipaw_wfs_rcut)
        WRITE (*,*) "shape of upf(2)%gipaw_wfs_rcutus = ", shape(upf(2)%gipaw_wfs_rcutus)
        WRITE (*,*) "shape of upf(2)%gipaw_wfs_ae = ", shape(upf(2)%gipaw_wfs_ae)
        WRITE (*,*) "shape of upf(2)%gipaw_wfs_ps = ", shape(upf(2)%gipaw_wfs_ps)
     ENDIF
     !
     WRITE (*,*) "****************************************"
     !!! upf_vca
     WRITE (*,*) "upf_vca%nbeta = ", upf_vca%nbeta
     WRITE (*,*) "shape of upf_vca%kbeta = ", shape(upf_vca%kbeta)
     WRITE (*,*) "upf_vca%kbeta = ", upf_vca%kbeta
     WRITE (*,*) "upf_vca%kkbeta = ", upf_vca%kkbeta
     WRITE (*,*) "shape of upf_vca%lll = ", shape(upf_vca%lll)
     WRITE (*,*) "shape of upf_vca%beta = ", shape(upf_vca%beta)
     WRITE (*,*) "shape of upf_vca%els_beta = ", shape(upf_vca%els_beta)
     WRITE (*,*) "shape of upf_vca%dion = ", shape(upf_vca%dion)
     WRITE (*,*) "shape of upf_vca%qqq = ", shape(upf_vca%qqq)
     WRITE (*,*) "shape of upf_vca%qfuncl = ", shape(upf_vca%qfuncl)
     WRITE (*,*) "shape of upf_vca%qfcoef = ", shape(upf_vca%qfcoef)
     WRITE (*,*) "shape of upf_vca%rcut = ", shape(upf_vca%rcut)
     WRITE (*,*) "shape of upf_vca%rcutus = ", shape(upf_vca%rcutus)
     WRITE (*,*) " "
     WRITE (*,*) "is_ultrasoft,upf_vca%tvanp = ", upf_vca%tvanp
     WRITE (*,*) "is_paw,upf_vca%tpawp = ", upf_vca%tpawp
     WRITE (*,*) "is_coulomb,upf_vca%tcoulombp = ", upf_vca%tcoulombp
     WRITE (*,*) " "
     !
     IF (matches(upf_vca%typ, "PAW")) THEN
        WRITE (*,*) "PAW"
        WRITE (*,*) "upf_vca%nbeta = ", upf_vca%nbeta, " = ", upf(1)%nbeta, " + ", upf(2)%nbeta, &
           " = ", upf(1)%nbeta+upf(2)%nbeta
        WRITE (*,*) "upf_vca%lmax = ", upf_vca%lmax
        WRITE (*,*) "  (2*lmax+1) = ", 2*upf_vca%lmax+1
        WRITE (*,*) "upf_vca%mesh = ", upf_vca%mesh
        WRITE (*,*) "shape of upf_vca%paw%augmom = ", shape(upf_vca%paw%augmom)
        !WRITE (*,*) "shape of upf_vca%paw%augfun = ", shape(upf_vca%paw%augfun)
        WRITE (*,*) "shape of upf_vca%paw%ae_rho_atc = ", shape(upf_vca%paw%ae_rho_atc)
        WRITE (*,*) "shape of upf_vca%aewfc = ", shape(upf_vca%aewfc)
        WRITE (*,*) "shape of upf_vca%pswfc = ", shape(upf_vca%pswfc)
        WRITE (*,*) "shape of upf_vca%paw%pfunc = ", shape(upf_vca%paw%pfunc)
        WRITE (*,*) "shape of upf_vca%paw%pfunc_rel = ", shape(upf_vca%paw%pfunc_rel)
        WRITE (*,*) "shape of upf_vca%paw%ptfunc = ", shape(upf_vca%paw%ptfunc)
        WRITE (*,*) "shape of upf_vca%paw%aewfc_rel = ", shape(upf_vca%paw%aewfc_rel)
        WRITE (*,*) "shape of upf_vca%paw%ae_vloc = ", shape(upf_vca%paw%ae_vloc)
        WRITE (*,*) "shape of upf_vca%paw%oc = ", shape(upf_vca%paw%oc)
        !WRITE (*,*) "shape of upf_vca%kdiff = ", shape(upf_vca%kdiff)
        !WRITE (*,*) "shape of upf_vca%sqr = ", shape(upf_vca%sqr)
        WRITE (*,*) "upf_vca%paw%augshape = ", upf_vca%paw%augshape
        WRITE (*,*) "upf_vca%paw%raug = ", upf_vca%paw%raug
        WRITE (*,*) "upf_vca%paw%iraug = ", upf_vca%paw%iraug
        WRITE (*,*) "upf_vca%paw%lmax_aug = ", upf_vca%paw%lmax_aug
        WRITE (*,*) "upf_vca%qqq_eps  = ", upf_vca%qqq_eps 
        WRITE (*,*) "upf_vca%paw%core_energy  = ", upf_vca%paw%core_energy
        WRITE (*,*) " "
     ENDIF
     !
     IF (upf_vca%has_gipaw) THEN
        WRITE (*,*) "GIPAW"
        WRITE (*,*) "upf_vca%nwfc = ", upf_vca%nwfc
        WRITE (*,*) "shape of upf_vca%els = ", shape(upf_vca%els)
        WRITE (*,*) "shape of upf_vca%nchi = ", shape(upf_vca%nchi)
        WRITE (*,*) "shape of upf_vca%lchi = ", shape(upf_vca%lchi)
        WRITE (*,*) "shape of upf_vca%oc = ", shape(upf_vca%oc)
        WRITE (*,*) "shape of upf_vca%rcut_chi = ", shape(upf_vca%rcut_chi)
        WRITE (*,*) "shape of upf_vca%rcutus_chi = ", shape(upf_vca%rcutus_chi)
        WRITE (*,*) "shape of upf_vca%epseu = ", shape(upf_vca%epseu)
        WRITE (*,*) "upf_vca%lmax = ", upf_vca%lmax
        WRITE (*,*) " (2*lmax+1) = ", 2*upf_vca%lmax+1
        WRITE (*,*) "shape of upf_vca%rinner = ", shape(upf_vca%rinner)
        !
        WRITE (*,*) "upf_vca%gipaw_ncore_orbitals = ", upf_vca%gipaw_ncore_orbitals
        WRITE (*,*) "upf_vca%mesh = ", upf_vca%mesh
        WRITE (*,*) "shape of upf_vca%gipaw_core_orbital_el = ", shape(upf_vca%gipaw_core_orbital_el)
        WRITE (*,*) "shape of upf_vca%gipaw_core_orbital_n = ", shape(upf_vca%gipaw_core_orbital_n)
        WRITE (*,*) "shape of upf_vca%gipaw_core_orbital_l = ", shape(upf_vca%gipaw_core_orbital_l)
        WRITE (*,*) "shape of upf_vca%gipaw_core_orbital = ", shape(upf_vca%gipaw_core_orbital)
        WRITE (*,*) "shape of upf_vca%gipaw_vlocal_ae = ", shape(upf_vca%gipaw_vlocal_ae)
        WRITE (*,*) "shape of upf_vca%gipaw_vlocal_ps = ", shape(upf_vca%gipaw_vlocal_ps)
        !
        WRITE (*,*) "upf_vca%gipaw_wfs_nchannels = ", upf_vca%gipaw_wfs_nchannels
        WRITE (*,*) "upf_vca%mesh = ", upf_vca%mesh
        WRITE (*,*) "shape of upf_vca%gipaw_wfs_el = ", shape(upf_vca%gipaw_wfs_el)
        !WRITE (*,*) "shape of upf_vca%gipaw_wfs_n = ", shape(upf_vca%gipaw_wfs_n)
        !WRITE (*,*) "shape of upf_vca%gipaw_wfs_l = ", shape(upf_vca%gipaw_wfs_l)
        WRITE (*,*) "shape of upf_vca%gipaw_wfs_ll = ", shape(upf_vca%gipaw_wfs_ll)
        WRITE (*,*) "shape of upf_vca%gipaw_wfs_rcut = ", shape(upf_vca%gipaw_wfs_rcut)
        WRITE (*,*) "shape of upf_vca%gipaw_wfs_rcutus = ", shape(upf_vca%gipaw_wfs_rcutus)
        WRITE (*,*) "shape of upf_vca%gipaw_wfs_ae = ", shape(upf_vca%gipaw_wfs_ae)
        WRITE (*,*) "shape of upf_vca%gipaw_wfs_ps = ", shape(upf_vca%gipaw_wfs_ps)
        WRITE (*,*) " "
     ENDIF
  ENDIF
  !!! DEBUG

END SUBROUTINE compute_virtual
