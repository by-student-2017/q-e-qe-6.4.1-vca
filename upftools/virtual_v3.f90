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
! virtual_v3: supports reading of UPF v2 format (USPP, PAW and GIPAW)
! Author: By Student
!
PROGRAM virtual_test

  !USE pseudo_types, ONLY : pseudo_upf, nullify_pseudo_upf, &
  !                         deallocate_pseudo_upf
  !USE upf_module, ONLY : read_upf
  !USE emend_upf_module, ONLY: make_emended_upf_copy
  !USE wrappers, ONLY: f_remove
  !USE write_upf_module, ONLY : write_upf
  !USE radial_grids, ONLY : radial_grid_type, nullify_radial_grid
  !USE environment, ONLY: environment_start, environment_end
  !USE mp_global, ONLY: mp_startup, mp_global_end
  !USE io_global, ONLY: ionode, stdout
  USE pseudo_types
  USE emend_upf_module
  USE wrappers
  USE upf_module
  USE write_upf_module
  USE radial_grids
  USE environment
  USE mp_global
  USE io_global

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
       CALL nullify_radial_grid(grid(is))
       INQUIRE ( FILE = TRIM(filein(is)), EXIST = exst )  
       IF (.NOT. exst ) CALL errore ( 'virtual_v2.x: ', TRIM(filein(is)) // ' not found', 5)  
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
     CALL write_upf ( trim(fileout), upf_vca, SCHEMA='v2')
     !
     CLOSE(ounps)
     CALL deallocate_pseudo_upf(upf(1))
     CALL deallocate_pseudo_upf(upf(2))
     !     ----------------------------------------------------------
     WRITE (stdout,"('Pseudopotential successfully written')")
     WRITE (stdout,"('Please review the content of the PP_INFO fields')")
     WRITE (stdout,"('*** Please TEST BEFORE USING !!! ***')")
     !     ----------------------------------------------------------
     !
   ENDIF
  CALL environment_end('VIRTUAL_V2.X')
#if defined(__MPI)
  CALL mp_global_end()
#endif

   STOP
20 CALL errore ('virtual.x', 'error reading pseudo file', 1)

END PROGRAM virtual_test

!
!---------------------------------------------------------------------
SUBROUTINE compute_virtual(x, filein, upf, upf_vca)

  !USE pseudo_types, ONLY : pseudo_upf
  USE pseudo_types
  USE splinelib
  !USE funct, ONLY : set_dft_from_name, get_iexch, get_icorr, get_igcx, get_igcc
  USE funct

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
  INTEGER :: upf_paw_data_format      ! The version of the format
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
  !
  ! PAW
  ! pp_augmentation
  real(8), ALLOCATABLE :: upf_paw_augmom(:,:,:)
  !
  ! pp augfun
  real(8), ALLOCATABLE :: upf_paw_augfun(:,:,:,:)
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
  real(8), ALLOCATABLE :: upf_paw_kdiff(:,:)
  !
  ! pp_sqr
  real(8), ALLOCATABLE :: upf_paw_sqr(:)
  !
  ! GIPAW
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
  upf_rel = -1
  upf_rcloc = 0.d0
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
  IF (matches(upf(1)%typ, "USPP") .and. matches(upf(2)%typ, "USPP") .and. &
    (upf(1)%has_gipaw .eqv. upf(2)%has_gipaw)) THEN
     upf_vca%typ = "USPP"
  ELSE
     CALL errore('virtual_v3.x: ', 'potential types are not match !!!', 1) 
  ENDIF
  IF (matches(upf(1)%typ, "PAW") .and. matches(upf(2)%typ, "PAW") .and.  &
    (upf(1)%has_gipaw .eqv. upf(2)%has_gipaw) .and. (upf(1)%paw%augshape == upf(2)%paw%augshape)) THEN
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
     WRITE (*,*) " pseudopotentials have different mesh "
     WRITE (*,*) upf(1)%mesh, upf(2)%mesh
     WRITE (*,*) upf(1)%r(1), upf(2)%r(1)
     WRITE (*,*) upf(1)%r(upf(1)%mesh), upf(2)%r(upf(2)%mesh)
     interpolate = .true.
  ENDIF

  upf_mesh = upf(1)%mesh
  upf_nbeta = upf(1)%nbeta+upf(2)%nbeta
  upf_ntwfc = upf(1)%nwfc
  upf_nlcc  = upf(1)%nlcc.or.upf(2)%nlcc
  upf_ecutrho = upf(1)%ecutrho
  upf_ecutwfc = upf(1)%ecutwfc
  upf_etotps  = upf(1)%etotps
  !
  
  !
  IF (matches(upf(1)%typ, "PAW") .or. matches(upf(2)%typ, "PAW")) THEN
     upf_lmax_rho = max(upf(1)%lmax_rho, upf(2)%lmax_rho)
     upf_paw_lmax_aug = max(upf(1)%paw%lmax_aug, upf(2)%paw%lmax_aug)
     upf_paw_raug     = max(upf(1)%paw%raug, upf(2)%paw%raug) ! r_match_augfun
     upf_paw_iraug    = max(upf(1)%paw%iraug, upf(2)%paw%iraug)
     upf_xmin     = min(upf(1)%xmin, upf(2)%xmin)
     upf_rmax     = max(upf(1)%rmax, upf(2)%rmax)
     upf_zmesh    = max(upf(1)%zmesh, upf(2)%zmesh)
     upf_dx       = max(upf(1)%dx, upf(2)%dx)
  ENDIF

  ALLOCATE( upf_ocw(upf_ntwfc), upf_elsw(upf_ntwfc), upf_lchiw(upf_ntwfc) )
  upf_ocw(1:upf_ntwfc)  = upf(1)%oc(1:upf_ntwfc)
  upf_elsw(1:upf_ntwfc) = upf(1)%els(1:upf_ntwfc)
  upf_lchiw(1:upf_ntwfc) = upf(1)%lchi(1:upf_ntwfc)
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
     WRITE (*,*) " pseudopotentials have different mesh "
     WRITE (*,*) "capel = ", capel
     interpolate = .true.
  ENDIF

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
     WRITE (*,*) " interpolate vloc0"
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
        WRITE (*,*) " interpolate betar"
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
     ELSE
        upf_ikk2(ib) = upf(2)%kbeta(i)
     ENDIF

  ENDDO
  !

  WRITE (*,*) "upf(1)%lll = ", upf(1)%lll
  WRITE (*,*) "upf(2)%lll = ", upf(2)%lll






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


  IF (matches(upf_vca%typ, "USPP") .or. matches(upf_vca%typ, "PAW")) THEN

    !pp_qij
    IF (upf(1)%nqf/=upf(2)%nqf) &
        CALL errore ("Virtual","different nqf are not implemented (yet)", 1)

    IF (upf(1)%nqlc/=upf(2)%nqlc) &
        CALL errore ("Virtual","different nqlc are not implemented (yet)", 1)

    upf_nqf = upf(1)%nqf
    upf_nqlc = upf(1)%nqlc

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
                   WRITE (*,*) " interpolate qfunc"
                   aux2(1,1:upf(2)%mesh) = upf(2)%qfuncl(1:upf(2)%mesh,ijv,l)
                   CALL dosplineint( upf(2)%r(1:upf(2)%mesh), aux2, upf_r(1:upf_mesh), aux1 )
                   upf_qfuncl(1:upf_mesh, ijv2, l) = (1.d0-x) * aux1(1,1:upf_mesh)
                   WRITE (*,*) " done"
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
  IF (upf(1)%nwfc==upf(2)%nwfc) THEN
     DO i=1,upf(2)%nwfc
        IF (interpolate) THEN
           WRITE (*,*) " interpolate chi"
           aux2(1,1:upf(2)%mesh) = upf(2)%chi(1:upf(2)%mesh,i)
           CALL dosplineint( upf(2)%r(1:upf(2)%mesh), aux2, upf_r(1:upf_mesh), aux1 )
           ! chi(1:upf_mesh,i,2) = aux1(1,1:upf_mesh)
           upf_chi(1:upf_mesh,i) =    x     * upf(1)%chi(1:upf_mesh,i) + &
              (1.d0-x) * aux1 (1,1:upf_mesh)
        ELSE
        upf_chi(1:upf_mesh,i) =    x     * upf(1)%chi(1:upf_mesh,i) + &
                                (1.d0-x) * upf(2)%chi(1:upf_mesh,i)
        ENDIF
        ! Jivtesh - The wavefunctions are calculated to be the average of the
        ! wavefunctions of the two atoms - lines 365-366
     ENDDO
  ELSE
     WRITE (*,*) "Number of wavefunctions not the same for the two pseudopotentials"
  ENDIF
  !upf_chi(1:upf_mesh,1:upf_ntwfc) = chi(1:upf_mesh,1:upf_ntwfc,1)
  !
  !pp_rhoatm
  ALLOCATE ( upf_rho_at(upf_mesh) )
  IF (interpolate) THEN
     WRITE (*,*) " interpolate rho_at"
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
  ENDIF

  IF (matches(upf_vca%typ, "PAW")) THEN
     !
     !pp_augmentation
     ALLOCATE ( upf_paw_augmom(upf_nbeta,upf_nbeta,upf_paw_lmax_aug) )
     upf_paw_augmom(1:upf_nbeta,1:upf_nbeta,0:upf_paw_lmax_aug) = 0.0
     upf_paw_augmom(1:upf(1)%nbeta,1:upf(1)%nbeta,0:upf(1)%paw%lmax_aug) = &
       upf(1)%paw%augmom(1:upf(1)%nbeta,1:upf(1)%nbeta,0:upf(1)%paw%lmax_aug)
     upf_paw_augmom(upf(1)%nbeta+1:upf_nbeta,upf(1)%nbeta+1:upf_nbeta,0:upf(2)%paw%lmax_aug) = &
       upf(2)%paw%augmom(1:upf(2)%nbeta,1:upf(2)%nbeta,0:upf(2)%paw%lmax_aug)
     !
     !ALLOCATE ( upf_paw_augfun(1:upf_mesh,upf_nbeta,upf_nbeta,upf_paw_lmax_aug) )
     !upf_paw_augfun(1:upf_mesh,1:upf_nbeta,1:upf_nbeta,0:upf_paw_lmax_aug) = 0.0
     !upf_paw_augfun(1:upf_mesh,1:upf(1)%nbeta,1:upf(1)%nbeta,0:upf(1)%paw%lmax_aug) = &
     !  upf(1)%paw%augfun(1:upf_mesh,1:upf(1)%nbeta,1:upf(1)%nbeta,0:upf(1)%paw%lmax_aug)
     !upf_paw_augfun(1:upf_mesh,upf(1)%nbeta+1:upf_nbeta,upf(1)%nbeta+1:upf_nbeta,0:upf(2)%paw%lmax_aug) = &
     !  upf(2)%paw%augfun(1:upf_mesh,1:upf(2)%nbeta,1:upf(2)%nbeta,0:upf(2)%paw%lmax_aug)
     !
     !pp_paw_ae_rho_atc
     ALLOCATE ( upf_paw_ae_rho_atc(upf_mesh) )
     IF (interpolate) THEN
        WRITE (*,*) " interpolate ae_rho_atc"
        aux2(1,1:upf(2)%mesh) = upf(2)%paw%ae_rho_atc(1:upf(2)%mesh)
        CALL dosplineint( upf(2)%r(1:upf(2)%mesh), aux2, upf_r(1:upf_mesh), aux1 )
        ! upf_paw_ae_rho_atc(1:upf_mesh,2) = aux1(1,1:upf_mesh)
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
     DO j=1,upf(1)%nbeta
        upf_aewfc(1:upf_mesh,j) = x * upf(1)%aewfc(1:upf_mesh,j)
     ENDDO
     DO j=1,upf(2)%nbeta
       IF (interpolate) THEN
          WRITE (*,*) " interpolate aewfc"
          aux2(1,1:upf(2)%mesh) = upf(2)%aewfc(1:upf(2)%mesh,j)
          CALL dosplineint( upf(2)%r(1:upf(2)%mesh), aux2, upf_r(1:upf_mesh), aux1 )
          ! upf_aewfc(1:upf_mesh,i,2) = aux1(1,1:upf_mesh)
          upf_aewfc(1:upf_mesh, upf(1)%nbeta+j) = (1.d0-x) * aux1(1,1:upf_mesh)
       ELSE
          upf_aewfc(1:upf_mesh, upf(1)%nbeta+j) = (1.d0-x) * upf(2)%aewfc(1:upf_mesh,j)
       ENDIF
     ENDDO
     !
     !pp_pswfc
     ALLOCATE ( upf_pswfc(upf_mesh,upf_nbeta) )
     DO j=1,upf(1)%nbeta
        upf_pswfc(1:upf_mesh,j) = x * upf(1)%pswfc(1:upf_mesh,j)
     ENDDO
     DO j=1,upf(2)%nbeta
       IF (interpolate) THEN
          WRITE (*,*) " interpolate pswfc"
          aux2(1,1:upf(2)%mesh) = upf(2)%pswfc(1:upf(2)%mesh,j)
          CALL dosplineint( upf(2)%r(1:upf(2)%mesh), aux2, upf_r(1:upf_mesh), aux1 )
          ! upf_pswfc(1:upf_mesh,i,2) = aux1(1,1:upf_mesh)
          upf_pswfc(1:upf_mesh, upf(1)%nbeta+j) = (1.d0-x) * aux1(1,1:upf_mesh)
       ELSE
          upf_pswfc(1:upf_mesh, upf(1)%nbeta+j) = (1.d0-x) * upf(2)%pswfc(1:upf_mesh,j)
       ENDIF
     ENDDO
     !
     !pp_paw_ae_vloc
     ALLOCATE ( upf_paw_ae_vloc(upf_mesh) )
     IF (interpolate) THEN
        WRITE (*,*) " interpolate ae_vloc"
        aux2(1,1:upf(2)%mesh) = upf(2)%paw%ae_vloc(1:upf(2)%mesh)
        CALL dosplineint( upf(2)%r(1:upf(2)%mesh), aux2, upf_r(1:upf_mesh), aux1 )
        ! upf_paw_ae_vloc(1:upf_mesh,2) = aux1(1,1:upf_mesh)
        upf_paw_ae_vloc(1:upf_mesh) =   x     * upf(1)%paw%ae_vloc(1:upf_mesh) + &
                                     (1.d0-x) * aux1(1,1:upf_mesh)
        WRITE (*,*) " done"
     ELSE
        upf_paw_ae_vloc(1:upf_mesh) =   x     * upf(1)%paw%ae_vloc(1:upf_mesh) + &
                                     (1.d0-x) * upf(2)%paw%ae_vloc(1:upf_mesh)
     ENDIF
     !
     !pp_kdiff
     !ALLOCATE ( upf_kdiff(upf_nbeta,upf_nbeta) )
     !upf_kdiff(1:upf_nbeta,1:upf_nbeta) = 0.0
     !upf_kdiff(1:upf(1)%nbeta,1:upf(1)%nbeta) = &
     !  upf(1)%kdiff(1:upf(1)%nbeta,1:upf(1)%nbeta)
     !upf_kdiff(upf(1)%nbeta+1:upf_nbeta,upf(1)%nbeta+1:upf_nbeta) = &
     !  upf(2)%kdiff(1:upf(2)%nbeta,1:upf(2)%nbeta)
     !
     !pp_grid_recon
     !pp_sqrt
     !ALLOCATE ( upf_sqr(upf_mesh) )
     !IF (interpolate) THEN
     !   WRITE (*,*) " interpolate sar"
     !   aux2(1,1:upf(2)%mesh) = upf(2)%sqr(1:upf(2)%mesh)
     !   CALL dosplineint( upf(2)%r(1:upf(2)%mesh), aux2, upf_r(1:upf_mesh), aux1 )
     !   ! upf_sqr(1:upf_mesh,2) = aux1(1,1:upf_mesh)
     !   upf_sqr(1:upf_mesh) =   x     * upf(1)%sqr(1:upf_mesh) + &
     !                        (1.d0-x) * aux1(1,1:upf_mesh)
     !   WRITE (*,*) " done"
     !ELSE
     !   upf_sar(1:upf_mesh) =   x     * upf(1)%sqr(1:upf_mesh) + &
     !                        (1.d0-x) * upf(2)%sqr(1:upf_mesh)
     !ENDIF
     !
  ENDIF

  IF (upf_vca%has_gipaw) THEN
     !
     !pp_gipaw_core_orbital
     ALLOCATE ( upf_gipaw_core_orbital(upf_mesh,upf_nbeta) )
     DO j=1,upf(1)%nbeta
        upf_gipaw_core_orbital(1:upf_mesh,j) = x * upf(1)%gipaw_core_orbital(1:upf_mesh,j)
     ENDDO
     DO j=1,upf(2)%nbeta
       IF (interpolate) THEN
          WRITE (*,*) " interpolate gipaw_core_orbital"
          aux2(1,1:upf(2)%mesh) = upf(2)%gipaw_core_orbital(1:upf(2)%mesh,j)
          CALL dosplineint( upf(2)%r(1:upf(2)%mesh), aux2, upf_r(1:upf_mesh), aux1 )
          ! upf_gipaw_core_orbital(1:upf_mesh,i,2) = aux1(1,1:upf_mesh)
          upf_gipaw_core_orbital(1:upf_mesh, upf(1)%nbeta+j) = (1.d0-x) * aux1(1,1:upf_mesh)
       ELSE
          upf_gipaw_core_orbital(1:upf_mesh, upf(1)%nbeta+j) = (1.d0-x) * upf(2)%gipaw_core_orbital(1:upf_mesh,j)
       ENDIF
     ENDDO
     !
     !pp_gipaw_local_data
     !pp_gipaw_vlocal_ae
     ALLOCATE( upf_gipaw_vlocal_ae(upf_mesh) )
     IF (interpolate) THEN
        WRITE (*,*) " interpolate gipaw_vlocal_ae"
        aux2(1,1:upf(2)%mesh) =  upf(2)%gipaw_vlocal_ae(1:upf(2)%mesh)

        CALL dosplineint( upf(2)%r(1:upf(2)%mesh), aux2, upf_r(1:upf_mesh), aux1 )

        ! upf(2)%vloc_ae(1:upf_mesh) = aux1(1,1:upf_mesh)

        ! Jivtesh - if the mesh of the first atom extends to a larger radius
        ! than the mesh of the second atom, then, for those radii that are
        ! greater than the maximum radius of the second atom, the local potential
        ! of the second atom is calculated using the expression
        ! v_local = (-2)*Z/r instead of using the extrapolated value.
        ! This is because, typically extrapolation leads to positive potentials.
        ! This is implemented in lines 240-242

        DO i=1,upf(1)%mesh
           IF ( upf(1)%r(i) > upf(2)%r(upf(2)%mesh) ) &
              ! upf(2)%vloc_ae(i) = -(2.0*upf(2)%zp)/upf(1)%r(i)
              aux1(1,i) = -(2.0*upf(2)%zp)/upf(1)%r(i)
        ENDDO
        upf_gipaw_vlocal_ae(1:upf_mesh) =    x     * upf(1)%gipaw_vlocal_ae(1:upf_mesh) + &
                                          (1.d0-x) * aux1(1,1:upf_mesh)
     ELSE
        upf_gipaw_vlocal_ae(1:upf_mesh) =    x     * upf(1)%gipaw_vlocal_ae(1:upf_mesh) + &
                                          (1.d0-x) * upf(2)%gipaw_vlocal_ae(1:upf_mesh)
     ENDIF
     !pp_gipaw_vlocal_ps
     ALLOCATE( upf_gipaw_vlocal_ps(upf_mesh) )
     IF (interpolate) THEN
        WRITE (*,*) " interpolate gipaw_vlocal_ps"
        aux2(1,1:upf(2)%mesh) =  upf(2)%gipaw_vlocal_ps(1:upf(2)%mesh)

        CALL dosplineint( upf(2)%r(1:upf(2)%mesh), aux2, upf_r(1:upf_mesh), aux1 )

        ! upf(2)%vloc_ps(1:upf_mesh) = aux1(1,1:upf_mesh)

        ! Jivtesh - if the mesh of the first atom extends to a larger radius
        ! than the mesh of the second atom, then, for those radii that are
        ! greater than the maximum radius of the second atom, the local potential
        ! of the second atom is calculated using the expression
        ! v_local = (-2)*Z/r instead of using the extrapolated value.
        ! This is because, typically extrapolation leads to positive potentials.
        ! This is implemented in lines 240-242

        DO i=1,upf(1)%mesh
           IF ( upf(1)%r(i) > upf(2)%r(upf(2)%mesh) ) &
              ! upf(2)%vloc_ps(i) = -(2.0*upf(2)%zp)/upf(1)%r(i)
              aux1(1,i) = -(2.0*upf(2)%zp)/upf(1)%r(i)
        ENDDO
        upf_gipaw_vlocal_ps(1:upf_mesh) =    x     * upf(1)%gipaw_vlocal_ps(1:upf_mesh) + &
                                          (1.d0-x) * aux1(1,1:upf_mesh)
     ELSE
        upf_gipaw_vlocal_ps(1:upf_mesh) =    x     * upf(1)%gipaw_vlocal_ps(1:upf_mesh) + &
                                          (1.d0-x) * upf(2)%gipaw_vlocal_ps(1:upf_mesh)
     ENDIF
     !
     !pp_gipaw_orbitals
     !pp_gipaw_wfs_el
     ALLOCATE (upf_gipaw_wfs_el(upf_nbeta))
     upf_gipaw_wfs_el(1:upf(1)%nbeta) = upf(1)%gipaw_wfs_el
     upf_gipaw_wfs_el(upf(1)%nbeta+1:upf_nbeta) = upf(2)%gipaw_wfs_el
     !
     !pp_gipaw_wfs_n
     !ALLOCATE (upf_gipaw_wfs_n(upf_nbeta))
     !upf_gipaw_wfs_n(1:upf(1)%nbeta) = upf(1)%gipaw_wfs_n
     !upf_gipaw_wfs_n(upf(1)%nbeta+1:upf_nbeta) = upf(2)%gipaw_wfs_n
     !
     !pp_gipaw_wfs_l
     !ALLOCATE (upf_gipaw_wfs_l(upf_nbeta))
     !upf_gipaw_wfs_l(1:upf(1)%nbeta) = upf(1)%gipaw_wfs_l
     !upf_gipaw_wfs_l(upf(1)%nbeta+1:upf_nbeta) = upf(2)%gipaw_wfs_l
     !
     !pp_gipaw_wfs_ll
     ALLOCATE (upf_gipaw_wfs_ll(upf_nbeta))
     upf_gipaw_wfs_ll(1:upf(1)%nbeta) = upf(1)%gipaw_wfs_ll
     upf_gipaw_wfs_ll(upf(1)%nbeta+1:upf_nbeta) = upf(2)%gipaw_wfs_ll
     !
     !pp_gipaw_wfs_rcut
     ALLOCATE (upf_gipaw_wfs_rcut(upf_nbeta))
     upf_gipaw_wfs_rcut(1:upf(1)%nbeta) = upf(1)%gipaw_wfs_rcut
     upf_gipaw_wfs_rcut(upf(1)%nbeta+1:upf_nbeta) = upf(2)%gipaw_wfs_rcut
     !
     !pp_gipaw_wfs_rcutus
     ALLOCATE (upf_gipaw_wfs_rcutus(upf_nbeta))
     upf_gipaw_wfs_rcutus(1:upf(1)%nbeta) = upf(1)%gipaw_wfs_rcutus
     upf_gipaw_wfs_rcutus(upf(1)%nbeta+1:upf_nbeta) = upf(2)%gipaw_wfs_rcutus
     !
     !pp_gipaw_wfs_ae
     ALLOCATE ( upf_gipaw_wfs_ae(upf_mesh,upf_nbeta) )
     DO j=1,upf(1)%nbeta
        upf_gipaw_wfs_ae(1:upf_mesh,j) = x * upf(1)%gipaw_wfs_ae(1:upf_mesh,j)
     ENDDO
     DO j=1,upf(2)%nbeta
       IF (interpolate) THEN
          WRITE (*,*) " interpolate gipaw_wfs_ae"
          aux2(1,1:upf(2)%mesh) = upf(2)%gipaw_wfs_ae(1:upf(2)%mesh,j)
          CALL dosplineint( upf(2)%r(1:upf(2)%mesh), aux2, upf_r(1:upf_mesh), aux1 )
          ! upf_gipaw_wfs_ae(1:upf_mesh,i,2) = aux1(1,1:upf_mesh)
          upf_gipaw_wfs_ae(1:upf_mesh, upf(1)%nbeta+j) = (1.d0-x) * aux1(1,1:upf_mesh)
       ELSE
          upf_gipaw_wfs_ae(1:upf_mesh, upf(1)%nbeta+j) = (1.d0-x) * upf(2)%gipaw_wfs_ae(1:upf_mesh,j)
       ENDIF
     ENDDO
     !
     !pp_gipaw_wfs_ps
     ALLOCATE ( upf_gipaw_wfs_ps(upf_mesh,upf_nbeta) )
     DO j=1,upf(1)%nbeta
        upf_gipaw_wfs_ps(1:upf_mesh,j) = x * upf(1)%gipaw_wfs_ps(1:upf_mesh,j)
     ENDDO
     DO j=1,upf(2)%nbeta
       IF (interpolate) THEN
          WRITE (*,*) " interpolate gipaw_wfs_ps"
          aux2(1,1:upf(2)%mesh) = upf(2)%gipaw_wfs_ps(1:upf(2)%mesh,j)
          CALL dosplineint( upf(2)%r(1:upf(2)%mesh), aux2, upf_r(1:upf_mesh), aux1 )
          ! upf_gipaw_wfs_ps(1:upf_mesh,i,2) = aux1(1,1:upf_mesh)
          upf_gipaw_wfs_ps(1:upf_mesh, upf(1)%nbeta+j) = (1.d0-x) * aux1(1,1:upf_mesh)
       ELSE
          upf_gipaw_wfs_ps(1:upf_mesh, upf(1)%nbeta+j) = (1.d0-x) * upf(2)%gipaw_wfs_ps(1:upf_mesh,j)
       ENDIF
     ENDDO
     !
  ENDIF

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

  upf_vca%kkbeta = maxval(upf_vca%kbeta)

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
  ENDIF

  IF (matches(upf_vca%typ, "PAW")) THEN
     NULLIFY( upf_vca%paw%augmom )
     NULLIFY( upf_vca%paw%ae_rho_atc )
     NULLIFY( upf_vca%paw%aewfc_rel )
     NULLIFY( upf_vca%paw%pfunc )
     NULLIFY( upf_vca%paw%pfunc_rel )
     NULLIFY( upf_vca%paw%ptfunc )
     NULLIFY( upf_vca%paw%ae_vloc )
     NULLIFY( upf_vca%els, &
              upf_vca%lchi, &
              upf_vca%nchi, &
              upf_vca%paw%oc, &
              upf_vca%rcut_chi, &
              upf_vca%rcutus_chi, &
              upf_vca%epseu &
            )
     !
     ALLOCATE( upf_vca%paw%augmom(upf_nbeta,upf_nbeta, 0:2*upf_lmax))
     ALLOCATE( upf_vca%paw%ae_rho_atc(upf_mesh) )
     ALLOCATE( upf_vca%paw%pfunc(upf_mesh, upf_nbeta, upf_nbeta) )
     ALLOCATE( upf_vca%paw%pfunc_rel(upf_mesh, upf_nbeta, upf_nbeta) )
     ALLOCATE( upf_vca%paw%ptfunc(upf_mesh, upf_nbeta, upf_nbeta) )
     ALLOCATE( upf_vca%paw%ae_vloc(upf_mesh) )
     ALLOCATE( upf_vca%els(upf_ntwfc), &
               upf_vca%lchi(upf_ntwfc), &
               upf_vca%nchi(upf_ntwfc), &
               upf_vca%paw%oc(upf_ntwfc), &
               upf_vca%rcut_chi(upf_ntwfc), &
               upf_vca%rcutus_chi(upf_ntwfc), &
               upf_vca%epseu(upf_ntwfc) &
             )
     !
     upf_vca%paw%augshape = upf_paw_augshape
     upf_vca%paw%raug = upf_paw_raug
     !upf_vca%core_enery = upf_core_energy 
     upf_vca%paw%iraug = upf_paw_iraug
     !
     upf_vca%lmax_rho = upf_lmax_rho
     upf_vca%paw%lmax_aug = upf_paw_lmax_aug
     !
     upf_vca%paw%augmom = upf_paw_augmom
     upf_vca%paw%ae_rho_atc = upf_paw_ae_rho_atc
     upf_vca%aewfc = upf_aewfc
     !
     ALLOCATE( upf_vca%paw%ptfunc(upf_mesh, upf_nbeta, upf_nbeta) )
     upf_vca%paw%ptfunc = upf_paw_ptfunc
     !
     upf_vca%paw%raug    = upf_paw_raug
     upf_vca%paw%iraug   = upf_paw_iraug
     upf_vca%dx       = upf_dx
     upf_vca%xmin     = upf_xmin
     upf_vca%zmesh    = upf_zmesh
     !
     ! nl, a2
     upf_vca%els(1:upf(1)%nbeta) = upf(1)%els
     upf_vca%els(upf(1)%nbeta+1:upf_nbeta) = upf(2)%els
     !
     ! pn, i3
     upf_vca%nchi(1:upf(1)%nbeta) = x * upf(1)%nchi
     upf_vca%nchi(upf(1)%nbeta+1:upf_nbeta) = (1.d0-x) * upf(2)%nchi
     !
     ! l, i3
     upf_vca%lchi(1:upf(1)%nbeta) = x * upf(1)%lchi
     upf_vca%lchi(upf(1)%nbeta+1:upf_nbeta) = (1.d0-x) * upf(2)%lchi
     !
     ! occ, f6,2
     upf_vca%oc(1:upf(1)%nbeta) = x * upf(1)%oc
     upf_vca%oc(upf(1)%nbeta+1:upf_nbeta) = (1.d0-x) * upf(2)%oc
     !
     ! Rcut, f11.3
     upf_vca%rcut_chi(1:upf(1)%nbeta) = upf(1)%rcut_chi
     upf_vca%rcut_chi(upf(1)%nbeta+1:upf_nbeta) = upf(2)%rcut_chi
     !
     ! Rcut US, f11.3
     upf_vca%rcutus_chi(1:upf(1)%nbeta) = upf(1)%rcutus_chi
     upf_vca%rcutus_chi(upf(1)%nbeta+1:upf_nbeta) = upf(2)%rcutus_chi
     !
     ! E pseu, f11.3
     upf_vca%epseu(1:upf(1)%nbeta) = upf(1)%epseu
     upf_vca%epseu(upf(1)%nbeta+1:upf_nbeta) = upf(2)%epseu
     !
     IF (upf_vca%has_so) THEN
        upf_vca%paw%pfunc_rel = upf_paw_pfunc_rel
        !
        upf_vca%paw%aewfc_rel = upf_paw_aewfc_rel
     ENDIF
  ENDIF

  IF (upf_vca%has_gipaw) THEN
     NULLIFY( upf_vca%grid ) 
     NULLIFY( upf_vca%els, upf_vca%lchi, upf_vca%nchi, upf_vca%jchi, upf_vca%oc )
     NULLIFY( upf_vca%r, upf_vca%rab )
     NULLIFY( upf_vca%rho_atc, upf_vca%vloc )
     NULLIFY( upf_vca%nn)
     NULLIFY( upf_vca%els_beta)
     NULLIFY( upf_vca%rcut_chi, upf_vca%rcutus_chi )
     NULLIFY( upf_vca%epseu)
     NULLIFY( upf_vca%vnl)
     NULLIFY( upf_vca%aewfc, upf_vca%pswfc )
     NULLIFY( upf_vca%rinner )
     NULLIFY( upf_vca%chi )
     NULLIFY( upf_vca%rho_at )
     NULLIFY ( upf_vca%gipaw_core_orbital_n )
     NULLIFY ( upf_vca%gipaw_core_orbital_l )
     NULLIFY ( upf_vca%gipaw_core_orbital_el )
     NULLIFY ( upf_vca%gipaw_core_orbital )
     NULLIFY ( upf_vca%gipaw_vlocal_ae )
     NULLIFY ( upf_vca%gipaw_vlocal_ps )
     NULLIFY ( upf_vca%gipaw_wfs_el )
     NULLIFY ( upf_vca%gipaw_wfs_ll )
     NULLIFY ( upf_vca%gipaw_wfs_ae )
     NULLIFY ( upf_vca%gipaw_wfs_rcut )
     NULLIFY ( upf_vca%gipaw_wfs_rcutus )
     NULLIFY ( upf_vca%gipaw_wfs_ps )
     ALLOCATE( upf_vca%els(upf_ntwfc), &
               upf_vca%oc(upf_ntwfc), &
               upf_vca%lchi(upf_ntwfc), &
               upf_vca%nchi(upf_ntwfc), &
               upf_vca%rcut_chi(upf_ntwfc), &
               upf_vca%rcutus_chi(upf_ntwfc), &
               upf_vca%epseu(upf_ntwfc) &
             )
     !
     upf_vca%tvanp = upf(1)%tvanp.or.upf(2)%tvanp
     upf_vca%tcoulombp = upf(1)%tcoulombp.or.upf(2)%tcoulombp
     upf_vca%nlcc    = upf_nlcc
     !upf_vca%dft     = upf_dft
     upf_vca%zp      = upf_zp
     upf_vca%etotps  = upf_etotps
     upf_vca%ecutwfc = upf_ecutwfc
     upf_vca%ecutrho = upf_ecutrho
     !upf_vca%nv      = upf_nv
     upf_vca%lmax    = upf_lmax
     upf_vca%lmax_rho = upf_lmax_rho
     upf_vca%nwfc    = upf_ntwfc
     upf_vca%nbeta   = upf_nbeta
     !upf_vca%kkbeta  = upf_kkbeta
     upf_vca%mesh    = upf_mesh
     upf_vca%xmin    = upf_xmin
     upf_vca%rmax    = upf_rmax
     !upf_vca%zmesh   = upf_zmesh
     upf_vca%dx      = upf_dx
     !upf_vca%lloc    = upf_lloc
     upf_vca%rcloc   = upf_rcloc
     !upf_vca%q_with_l = upf_q_with_l
     upf_vca%nqf     = upf_nqf
     upf_vca%nqlc    = upf_nqlc
     upf_vca%qqq_eps = upf_qqq_eps
     upf_vca%has_wfc = upf(1)%has_wfc.or.upf(2)%has_wfc
     upf_vca%paw_data_format = upf_paw_data_format
     upf_vca%tpawp   = upf(1)%tpawp.or.upf(2)%tpawp
     upf_vca%has_gipaw = upf(1)%has_gipaw.or.upf(2)%has_gipaw
     upf_vca%paw_as_gipaw = upf(1)%paw_as_gipaw.or.upf(2)%paw_as_gipaw
     !upf_vca%gipaw_data_format = upf_gipaw_data_format
     !upf_vca%gipaw_ncore_orbitals = upf_gipaw_ncore_orbitals
     !upf_vca%gipaw_wfs_nchannels = upf_gipaw_wfs_nchannels
     !
     upf_vca%els(1:upf(1)%nbeta) = upf(1)%els
     upf_vca%els(upf(1)%nbeta+1:upf_nbeta) = upf(2)%els
     !
     upf_vca%nchi(1:upf(1)%nbeta) = x * upf(1)%nchi
     upf_vca%nchi(upf(1)%nbeta+1:upf_nbeta) = (1.d0-x) * upf(2)%nchi
     !
     upf_vca%lchi(1:upf(1)%nbeta) = x * upf(1)%lchi
     upf_vca%lchi(upf(1)%nbeta+1:upf_nbeta) = (1.d0-x) * upf(2)%lchi
     !
     upf_vca%oc(1:upf(1)%nbeta) = x * upf(1)%oc
     upf_vca%oc(upf(1)%nbeta+1:upf_nbeta) = (1.d0-x) * upf(2)%oc
     !
     upf_vca%rcut_chi(1:upf(1)%nbeta) = upf(1)%rcut_chi
     upf_vca%rcut_chi(upf(1)%nbeta+1:upf_nbeta) = upf(2)%rcut_chi
     !
     upf_vca%rcutus_chi(1:upf(1)%nbeta) = upf(1)%rcutus_chi
     upf_vca%rcutus_chi(upf(1)%nbeta+1:upf_nbeta) = upf(2)%rcutus_chi
  ENDIF

  ! !! DEBUG
  ! !! upf(1)
  ! WRITE (*,*) "upf(1)%nbeta = ", upf(1)%nbeta
  ! WRITE (*,*) "shape of upf(1)%kbeta = ", shape(upf(1)%kbeta)
  ! WRITE (*,*) "upf(1)%kbeta = ", upf(1)%kbeta
  ! WRITE (*,*) "upf(1)%kkbeta = ", upf(1)%kkbeta
  ! WRITE (*,*) "shape of upf(1)%lll = ", shape(upf(1)%lll)
  ! WRITE (*,*) "shape of upf(1)%beta = ", shape(upf(1)%beta)
  ! WRITE (*,*) "shape of upf(1)%els_beta = ", shape(upf(1)%els_beta)
  ! WRITE (*,*) "shape of upf(1)%dion = ", shape(upf(1)%dion)
  ! WRITE (*,*) "shape of upf(1)%qqq = ", shape(upf(1)%qqq)
  ! WRITE (*,*) "shape of upf(1)%qfuncl = ", shape(upf(1)%qfuncl)
  ! WRITE (*,*) "shape of upf(1)%qfcoef = ", shape(upf(1)%qfcoef)
  ! WRITE (*,*) "upf(1)%aewfc = ", upf(1)%aewfc
  ! WRITE (*,*) "upf(1)%pswfc = ", upf(1)%pswfc
  ! WRITE (*,*) "shape of upf(1)%rcut = ", shape(upf(1)%rcut)
  ! WRITE (*,*) "shape of upf(1)%rcutus = ", shape(upf(1)%rcutus)
  ! WRITE (*,*) "PAW"
  ! WRITE (*,*) "shape of upf(1)%augmom = ", shape(upf(1)%augmom)
  ! WRITE (*,*) "shape of upf(1)%augfun = ", shape(upf(1)%augfun)
  ! WRITE (*,*) "shape of upf(1)%upf_ae_rho_atc = ", shape(upf(1)%upf_ae_rho_atc)
  ! WRITE (*,*) "shape of upf(1)%upf_aewfc = ", shape(upf(1)%upf_aewfc)
  ! WRITE (*,*) "shape of upf(1)%upf_pswfc = ", shape(upf(1)%upf_pswfc)
  ! WRITE (*,*) "shape of upf(1)%ae_vloc = ", shape(upf(1)%ae_vloc)
  ! WRITE (*,*) "shape of upf(1)%kdiff = ", shape(upf(1)%kdiff)
  ! WRITE (*,*) "shape of upf(1)%sqr = ", shape(upf(1)%sqr)
  ! WRITE (*,*) "GIPAW"
  ! WRITE (*,*) "shape of upf(1)%gipaw_core_orbital = ", shape(upf(1)%gipaw_core_orbital)
  ! WRITE (*,*) "shape of upf(1)%vlocal_ae = ", shape(upf(1)%vlocal_ae)
  ! WRITE (*,*) "shape of upf(1)%vlocal_ps = ", shape(upf(1)%vlocal_ps)
  ! WRITE (*,*) ""
  ! !!! upf(2)
  ! WRITE (*,*) "upf(2)%nbeta = ", upf(2)%nbeta
  ! WRITE (*,*) "shape of upf(2)%kbeta = ", shape(upf(2)%kbeta)
  ! WRITE (*,*) "upf(2)%kbeta = ", upf(2)%kbeta
  ! WRITE (*,*) "upf(2)%kkbeta = ", upf(2)%kkbeta
  ! WRITE (*,*) "shape of upf(2)%lll = ", shape(upf(2)%lll)
  ! WRITE (*,*) "shape of upf(2)%beta = ", shape(upf(2)%beta)
  ! WRITE (*,*) "shape of upf(2)%els_beta = ", shape(upf(2)%els_beta)
  ! WRITE (*,*) "shape of upf(2)%dion = ", shape(upf(2)%dion)
  ! WRITE (*,*) "shape of upf(2)%qqq = ", shape(upf(2)%qqq)
  ! WRITE (*,*) "shape of upf(2)%qfuncl = ", shape(upf(2)%qfuncl)
  ! WRITE (*,*) "shape of upf(2)%qfcoef = ", shape(upf(2)%qfcoef)
  ! WRITE (*,*) "upf(2)%aewfc = ", upf(2)%aewfc
  ! WRITE (*,*) "upf(2)%pswfc = ", upf(2)%pswfc
  ! WRITE (*,*) "shape of upf(2)%rcut = ", shape(upf(2)%rcut)
  ! WRITE (*,*) "shape of upf(2)%rcutus = ", shape(upf(2)%rcutus)
  ! WRITE (*,*) "PAW"
  ! WRITE (*,*) "shape of upf(2)%augmom = ", shape(upf(2)%augmom)
  ! WRITE (*,*) "shape of upf(2)%augfun = ", shape(upf(2)%augfun)
  ! WRITE (*,*) "shape of upf(2)%upf_ae_rho_atc = ", shape(upf(2)%upf_ae_rho_atc)
  ! WRITE (*,*) "shape of upf(2)%upf_aewfc = ", shape(upf(2)%upf_aewfc)
  ! WRITE (*,*) "shape of upf(2)%upf_pswfc = ", shape(upf(2)%upf_pswfc)
  ! WRITE (*,*) "shape of upf(2)%ae_vloc = ", shape(upf(2)%ae_vloc)
  ! WRITE (*,*) "shape of upf(2)%kdiff = ", shape(upf(2)%kdiff)
  ! WRITE (*,*) "shape of upf(2)%sqr = ", shape(upf(2)%sqr)
  ! WRITE (*,*) "GIPAW"
  ! WRITE (*,*) "shape of upf(2)%gipaw_core_orbital = ", shape(upf(2)%gipaw_core_orbital)
  ! WRITE (*,*) "shape of upf(2)%vlocal_ae = ", shape(upf(2)%vlocal_ae)
  ! WRITE (*,*) "shape of upf(2)%vlocal_ps = ", shape(upf(2)%vlocal_ps)
  ! WRITE (*,*) ""
  ! !!! upf_vca
  ! WRITE (*,*) "upf_vca%nbeta = ", upf_vca%nbeta
  ! WRITE (*,*) "shape of upf_vca%kbeta = ", shape(upf_vca%kbeta)
  ! WRITE (*,*) "upf_vca%kbeta = ", upf_vca%kbeta
  ! WRITE (*,*) "upf_vca%kkbeta = ", upf_vca%kkbeta
  ! WRITE (*,*) "shape of upf_vca%lll = ", shape(upf_vca%lll)
  ! WRITE (*,*) "shape of upf_vca%beta = ", shape(upf_vca%beta)
  ! WRITE (*,*) "shape of upf_vca%els_beta = ", shape(upf_vca%els_beta)
  ! WRITE (*,*) "shape of upf_vca%dion = ", shape(upf_vca%dion)
  ! WRITE (*,*) "shape of upf_vca%qqq = ", shape(upf_vca%qqq)
  ! WRITE (*,*) "shape of upf_vca%qfuncl = ", shape(upf_vca%qfuncl)
  ! WRITE (*,*) "shape of upf_vca%qfcoef = ", shape(upf_vca%qfcoef)
  ! WRITE (*,*) "shape of upf_vca%rcut = ", shape(upf_vca%rcut)
  ! WRITE (*,*) "shape of upf_vca%rcutus = ", shape(upf_vca%rcutus)
  ! WRITE (*,*) "PAW"
  ! WRITE (*,*) "shape of upf_vca%augmom = ", shape(upf_vca%augmom)
  ! WRITE (*,*) "shape of upf_vca%augfun = ", shape(upf_vca%augfun)
  ! WRITE (*,*) "shape of upf_vca%upf_ae_rho_atc = ", shape(upf_vca%upf_ae_rho_atc)
  ! WRITE (*,*) "shape of upf_vca%upf_aewfc = ", shape(upf_vca%upf_aewfc)
  ! WRITE (*,*) "shape of upf_vca%upf_pswfc = ", shape(upf_vca%upf_pswfc)
  ! WRITE (*,*) "shape of upf_vca%ae_vloc = ", shape(upf_vca%ae_vloc)
  ! WRITE (*,*) "shape of upf_vca%kdiff = ", shape(upf_vca%kdiff)
  ! WRITE (*,*) "shape of upf_vca%sqr = ", shape(upf_vca%sqr)
  ! WRITE (*,*) "GIPAW"
  ! WRITE (*,*) "shape of upf_vca%gipaw_core_orbital = ", shape(upf_vca%gipaw_core_orbital)
  ! WRITE (*,*) "shape of upf_vca%vlocal_ae = ", shape(upf_vca%vlocal_ae)
  ! WRITE (*,*) "shape of upf_vca%vlocal_ps = ", shape(upf_vca%vlocal_ps)
  ! !! DEBUG


END SUBROUTINE compute_virtual
