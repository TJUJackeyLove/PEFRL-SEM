
subroutine update_displ_PEFRL()

! low-memory Runge-Kutta time scheme

    use specfem_par
    use specfem_par_acoustic
    use specfem_par_elastic
    use specfem_par_poroelastic
    use pml_par

  implicit none

  if (ACOUSTIC_SIMULATION) then
    potential_dot_dot_acoustic = 0._CUSTOM_REAL
    if (FIX_UNDERFLOW_PROBLEM) potential_dot_dot_acoustic = VERYSMALLVAL
  endif

  if (ELASTIC_SIMULATION) then
    accel = 0._CUSTOM_REAL
    call update_displacement_elastic_PEFRL()
    if (FIX_UNDERFLOW_PROBLEM) accel = VERYSMALLVAL
  endif

  end subroutine update_displ_PEFRL

!
!-------------------------------------------------------------------------------------------------

  subroutine update_displacement_elastic_PEFRL()


    use specfem_par
    use specfem_par_elastic
    use pml_par
  
    implicit none
    real(kind=CUSTOM_REAL) :: alpha,beta

    ! current Runge-Kutta coefficients
    alpha = ALPHA_PEFRL(istage)
    beta = BETA_PEFRL(istage)

    call update_displelastic_PEFRL(NGLOB_AB,displ,veloc,deltat,alpha,istage)

    
  end subroutine update_displacement_elastic_PEFRL

!-------------------------------------------------------------------------------------------------
  subroutine update_displelastic_PEFRL(NGLOB,displ,veloc,deltat,alpha,istage)

  use constants, only: CUSTOM_REAL,NDIM

  implicit none
  integer,intent(in) :: istage
  integer,intent(in) :: NGLOB


  ! wavefields
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB),intent(inout) :: displ
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB),intent(in) :: veloc


  real(kind=CUSTOM_REAL),intent(in) :: deltat
  ! Runge-Kutta coefficients
  real(kind=CUSTOM_REAL),intent(in) :: alpha

  ! local parameters
  integer :: i

  do i = 1,NGLOB
    ! PEFRL: intermediate storage wavefields
    displ(:,i)= displ(:,i) + alpha * deltat * veloc(:,i)
    

    
  enddo

  end subroutine update_displelastic_PEFRL
  

!---------------------

  subroutine update_displ_PEFRL_backward()

! low-memory Runge-Kutta time scheme

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic

  end subroutine update_displ_PEFRL_backward

!
!-------------------------------------------------------------------------------------------------
!
! acoustic domains
!
!-------------------------------------------------------------------------------------------------
!

  subroutine update_potential_dot_acoustic_PEFRL()

! updates acceleration, velocity and displacement in acoustic region (outer core)

  use specfem_par
  use specfem_par_acoustic

  implicit none
  ! local parameters
  real(kind=CUSTOM_REAL) :: alpha,beta

  ! current Runge-Kutta coefficients
  alpha = ALPHA_PEFRL(istage)
  beta = BETA_PEFRL(istage)

  ! forward wavefields
  call update_acoustic_PEFRL(NGLOB_AB,NGLOB_AB_PEFRL, &
                             potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                             potential_acoustic_PEFRL,potential_dot_acoustic_PEFRL, &
                             deltat,alpha,beta)

  end subroutine update_potential_dot_acoustic_PEFRL

!
!-------------------------------------------------------------------------------------------------
!

  subroutine update_acoustic_PEFRL(NGLOB,NGLOB_PEFRL, &
                                   potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                                   potential_acoustic_PEFRL,potential_dot_acoustic_PEFRL, &
                                   deltat,alpha,beta)

  use constants, only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: NGLOB,NGLOB_PEFRL

  ! wavefields
  real(kind=CUSTOM_REAL), dimension(NGLOB),intent(inout) :: &
            potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic
  real(kind=CUSTOM_REAL), dimension(NGLOB_PEFRL),intent(inout) :: &
            potential_acoustic_PEFRL,potential_dot_acoustic_PEFRL

  real(kind=CUSTOM_REAL),intent(in) :: deltat
  ! Runge-Kutta coefficients
  real(kind=CUSTOM_REAL),intent(in) :: alpha,beta

  ! local parameters
  integer :: i

  if (NGLOB /= NGLOB_PEFRL) stop 'error in definition of NGLOB_PEFRL'

  ! low-memory Runge-Kutta scheme

  do i = 1,NGLOB
    ! low-memory Runge-Kutta: intermediate storage wavefields
    potential_dot_acoustic_PEFRL(i) = alpha * potential_dot_acoustic_PEFRL(i) + deltat * potential_dot_dot_acoustic(i)
    potential_acoustic_PEFRL(i) = alpha * potential_acoustic_PEFRL(i) + deltat * potential_dot_acoustic(i)
    ! updates wavefields
    potential_dot_acoustic(i) = potential_dot_acoustic(i) + beta * potential_dot_acoustic_PEFRL(i)
    potential_acoustic(i) = potential_acoustic(i) + beta * potential_acoustic_PEFRL(i)
  enddo

  end subroutine update_acoustic_PEFRL

!
!-------------------------------------------------------------------------------------------------
!
! elastic domains
!
!-------------------------------------------------------------------------------------------------
!

  subroutine update_veloc_elastic_PEFRL()

! updates acceleration,velocity and displacement in elastic regions (crust/mantle,inner core)

  use specfem_par
  use specfem_par_elastic

  implicit none
  ! local parameters
  real(kind=CUSTOM_REAL) :: alpha,beta

  ! current Runge-Kutta coefficients
  alpha = ALPHA_PEFRL(istage)
  beta = BETA_PEFRL(istage)

  ! forward wavefields
  call update_elastic_PEFRL(NGLOB_AB,NGLOB_AB_PEFRL,displ,veloc,accel, &
                            displ_PEFRL,veloc_PEFRL,deltat,alpha,beta,istage)

  end subroutine update_veloc_elastic_PEFRL

!
!-------------------------------------------------------------------------------------------------
!

! not used yet...
!  subroutine update_veloc_elastic_PEFRL_backward()
!
!! updates acceleration,velocity and displacement in elastic regions
!
!  implicit none
!
!  end subroutine update_veloc_elastic_PEFRL_backward

!
!-------------------------------------------------------------------------------------------------
!

  subroutine update_elastic_PEFRL(NGLOB,NGLOB_PEFRL,displ,veloc,accel, &
                                  displ_PEFRL,veloc_PEFRL,deltat,alpha,beta,istage)

  use constants, only: CUSTOM_REAL,NDIM

  implicit none

  integer,intent(in) :: NGLOB,NGLOB_PEFRL,istage

  ! wavefields
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB),intent(inout) :: displ,veloc,accel
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_PEFRL),intent(inout) :: displ_PEFRL,veloc_PEFRL

  real(kind=CUSTOM_REAL),intent(in) :: deltat
  ! Runge-Kutta coefficients
  real(kind=CUSTOM_REAL),intent(in) :: alpha,beta
  ! local parameters
  integer :: i

  if (NGLOB /= NGLOB_PEFRL) stop 'error in definition of NGLOB_PEFRL'
  
  ! PEFRL scheme

  do i = 1,NGLOB
    ! PEFRL: intermediate storage wavefields
    displ_PEFRL(:,i) = displ(:,i)
    veloc_PEFRL(:,i) = veloc_PEFRL(:,i) + beta * deltat * accel(:,i)
    
    ! updates wavefields 

    veloc(:,i) = veloc_PEFRL(:,i)
    
  enddo

   if ( istage  ==  4 ) then
    do i = 1,NGLOB
      displ(:,i) = displ_PEFRL(:,i) + 0.1786178958448091d0 * deltat * veloc(:,i)
    enddo
  end if
  
  end subroutine update_elastic_PEFRL


