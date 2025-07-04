

  subroutine print_stf_file()

  use constants, only: MAX_STRING_LEN,CUSTOM_REAL,IMAIN,IO_STF, &
    C_PEFRL,myrank

  use specfem_par, only: NSTEP,SIMULATION_TYPE,OUTPUT_FILES, &
    NSOURCES,islice_selected_source,ispec_selected_source,tshift_src,t0, &
    DT,USE_PEFRL,UNDO_ATTENUATION_AND_OR_PML

  use specfem_par_acoustic, only: ispec_is_acoustic
  use specfem_par_elastic, only: ispec_is_elastic
  use specfem_par_poroelastic, only: ispec_is_poroelastic

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL) :: stf_used,time_source
  real(kind=CUSTOM_REAL),dimension(NSTEP) :: source_time_function

  double precision :: stf,time_source_dble
  double precision,external :: get_stf_acoustic,get_stf_viscoelastic,get_stf_poroelastic

  integer :: isource,ispec,ier,it,istage,it_tmp_ext
  character(len=MAX_STRING_LEN) :: plot_file

  ! check
  if (SIMULATION_TYPE /= 1 .and. SIMULATION_TYPE /= 2 .and. SIMULATION_TYPE /= 3) &
    stop 'unrecognized SIMULATION_TYPE value in printing stf file'

  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'printing the source-time function'
    call flush_IMAIN()
  endif

  ! note: the source time function will be output for each source separately,
  !       instead of summing up all stf from single source contributions

  ! source contributions
  do isource = 1,NSOURCES

    ! initializes
    source_time_function(:) = 0._CUSTOM_REAL

    ! compute the source contribution (only if this proc carries the source)
    if (myrank == islice_selected_source(isource)) then

      ! time loop
      do it = 1,NSTEP
        ! current source time
        if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 2) then
          if (USE_PEFRL) then
            ! PEFRL
            ! note: the PEFRL scheme updates displacement after the stiffness computations and
            !       after adding boundary/coupling/source terms.
            !       thus, at each time loop step it, displ(:) is still at (n) and not (n+1) like for the Newmark scheme
            !       when entering this routine. we therefore at an additional -DT to have the corresponding timing for the source.
            istage = 1  ! only 1. stage output (corresponds to u(n))
            time_source_dble = dble(it)*DT  - t0 - tshift_src(isource)
          else
            time_source_dble = dble(it)*DT - t0 - tshift_src(isource)
          endif
          ! for external stf
          it_tmp_ext = it
        else
          ! backward simulation (SIMULATION_TYPE == 3)
          if (USE_PEFRL) then
            ! PEFRL
            ! note: the PEFRL scheme updates displacement after the stiffness computations and
            !       after adding boundary/coupling/source terms.
            !       thus, at each time loop step it, displ(:) is still at (n) and not (n+1) like for the Newmark scheme
            !       when entering this routine. we therefore at an additional -DT to have the corresponding timing for the source.
            istage = 1  ! only 1. stage output (corresponds to u(n))
            if (UNDO_ATTENUATION_AND_OR_PML) then
              ! stepping moves forward from snapshot position
              time_source_dble = dble(NSTEP-it)*DT- t0 - tshift_src(isource)
            else
              time_source_dble = dble(NSTEP-it)*DT - t0 - tshift_src(isource)
            endif
          else
            time_source_dble = dble(NSTEP-it)*DT - t0 - tshift_src(isource)
          endif
          ! for external stf
          it_tmp_ext = NSTEP - it + 1
        endif

        ispec = ispec_selected_source(isource)

        ! determines source time function value
        if (ispec_is_acoustic(ispec)) then
          stf = get_stf_acoustic(time_source_dble,isource,it_tmp_ext)
        else if (ispec_is_elastic(ispec)) then
          stf = get_stf_viscoelastic(time_source_dble,isource,it_tmp_ext)
        else if (ispec_is_poroelastic(ispec)) then
          stf = get_stf_poroelastic(time_source_dble,isource,it_tmp_ext)
        else
          call exit_MPI(myrank,'Invalid source element type, please check your mesh...')
        endif

        ! distinguishes between single and double precision for reals
        stf_used = real(stf,kind=CUSTOM_REAL)

        ! for file output
        source_time_function(it) = stf_used
      enddo
    endif

    ! main collects stf (if it does not already have it, i.e. if this source is not on the main)
    if (islice_selected_source(isource) /= 0) then
      if (myrank == 0) then
        ! main collects
        call recvv_cr(source_time_function,NSTEP,islice_selected_source(isource),0)
      else if (myrank == islice_selected_source(isource)) then
        ! secondary sends to main
        call sendv_cr(source_time_function,NSTEP,0,0)
      endif
    endif

    ! main prints out to file
    if (myrank == 0) then
      ! opens source time function file
      if (NSOURCES == 1) then
        plot_file = '/plot_source_time_function.txt'
      else if (isource < 10) then
        write(plot_file,"('/plot_source_time_function',i1,'.txt')") isource
      else if (isource < 100) then
        write(plot_file,"('/plot_source_time_function',i2,'.txt')") isource
      else
        write(plot_file,"('/plot_source_time_functionA',i7.7,'.txt')") isource
      endif
      open(unit=IO_STF,file=trim(OUTPUT_FILES)//trim(plot_file),status='unknown',iostat=ier)
      if (ier /= 0) call exit_mpi(myrank,'Error opening plot_source_time_function file')

      do it = 1,NSTEP
        ! overall time, note that tshift_src will start at zero for simulation
        if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 2) then
          if (USE_PEFRL) then
            ! PEFRL
            ! note: the PEFRL scheme updates displacement after the stiffness computations and
            !       after adding boundary/coupling/source terms.
            !       thus, at each time loop step it, displ(:) is still at (n) and not (n+1) like for the Newmark scheme
            !       when entering this routine. we therefore at an additional -DT to have the corresponding timing for the source.
            istage = 1  ! only 1. stage output (corresponds to u(n))
            time_source = dble(it-1)*DT - t0
          else
            time_source = dble(it-1)*DT - t0
          endif
        else
          ! backward simulation (SIMULATION_TYPE == 3)
          if (USE_PEFRL) then
            ! PEFRL
            ! note: the PEFRL scheme updates displacement after the stiffness computations and
            !       after adding boundary/coupling/source terms.
            !       thus, at each time loop step it, displ(:) is still at (n) and not (n+1) like for the Newmark scheme
            !       when entering this routine. we therefore at an additional -DT to have the corresponding timing for the source.
            istage = 1  ! only 1. stage output (corresponds to u(n))
            time_source = dble(NSTEP-it)*DT - t0
          else
            time_source = dble(NSTEP-it)*DT - t0
          endif
        endif

        ! file output
        write(IO_STF,*) time_source,source_time_function(it)
      enddo

      close(IO_STF)

    endif ! of if (myrank == 0) then

  enddo ! NSOURCES

  end subroutine print_stf_file
