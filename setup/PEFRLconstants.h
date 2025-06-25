
!!-----------------------------------------------------------
!!
!! symplectic time scheme - Position-Extended Forest-Ruth-Like (PEFRL) scheme (4th order)
!!
!!-----------------------------------------------------------

! Position-Extended Forest-Ruth-Like (PEFRL) scheme
!
! reference:
! Omelyan I P, Mryglod I M, Folk R. Optimized Forest–Ruth-and Suzuki-like algorithms for integration of motion in many-body systems. Computer Physics Communications, 2002, 146(2): 188-202.
!

! number of stages
  integer, parameter :: NSTAGE_PEFRL = 4

  double precision, parameter :: PEFRL_xi = 0.1786178958448091d0
  double precision, parameter :: PEFRL_lambda = -0.2123418310626054d0
  double precision, parameter :: PEFRL_chi = -0.06626458266981849d0

  real(kind=CUSTOM_REAL), dimension(NSTAGE_PEFRL), parameter :: ALPHA_PEFRL = &
    (/ PEFRL_xi, PEFRL_chi, 1.0_CUSTOM_REAL - 2.0_CUSTOM_REAL*(PEFRL_chi + PEFRL_xi), PEFRL_chi /)

  real(kind=CUSTOM_REAL), dimension(NSTAGE_PEFRL), parameter :: BETA_PEFRL = &
    (/ 0.5_CUSTOM_REAL - PEFRL_lambda, PEFRL_lambda, PEFRL_lambda, 0.5_CUSTOM_REAL - PEFRL_lambda /)

  real(kind=CUSTOM_REAL), dimension(NSTAGE_PEFRL), parameter :: C_PEFRL = &
    (/ 0.0_CUSTOM_REAL, 0.28765816893_CUSTOM_REAL, 0.5_CUSTOM_REAL, 0.71234183106_CUSTOM_REAL/)
  




